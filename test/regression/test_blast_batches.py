"""Query batching preserves exon evidence and searches one complete target DB."""
import argparse
from contextlib import closing
import gzip
import io
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

from test_exon_transcript_v3 import ROOT, rows
from test_identical_paralogs import ReferenceFixture
import align_exon_blastdb_v2 as aligner
import build_exon_blastdb_v2 as builder


class QueryBatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.query = self.folder / 'exons.fa.gz'
        with gzip.open(self.query, 'wt') as handle:
            handle.write('>E1.1 description\nAc\ngT\n>E2\nGG\n'
                         '>long_exon\nACGTacgtACGT\n>last\nTT\n')

    def test_boundaries_preserve_every_record_and_masking(self):
        batches = []
        for path, count, bases in aligner.iter_exon_query_batches(
                str(self.query), str(self.folder), 24):
            records = list(aligner.iter_query_records(path))
            batches.append(records)
            self.assertEqual(count, len(records))
            self.assertEqual(bases, sum(len(seq) for _, seq in records))
            self.assertLessEqual(Path(path).stat().st_size, 24)
        self.assertEqual(batches, [
            [('>E1.1', 'AcgT'), ('>E2', 'GG')],
            [('>long_exon', 'ACGTacgtACGT')], [('>last', 'TT')]])
        self.assertEqual(list(self.folder.iterdir()), [self.query])

    def test_early_close_cleans_scratch_and_disabled_mode_keeps_input(self):
        with closing(aligner.iter_exon_query_batches(
                str(self.query), str(self.folder), 24)) as batches:
            path, _, _ = next(batches)
            self.assertTrue(Path(path).exists())
        self.assertFalse(Path(path).exists())
        self.assertEqual(list(aligner.iter_exon_query_batches(
            str(self.query), str(self.folder), 0)), [(str(self.query), 0, 0)])
        self.assertTrue(self.query.exists())

    def test_invalid_input_fails_and_cleans_scratch(self):
        for text in ('', '>empty\n', 'ACGT\n', '>\nACGT\n'):
            path = self.folder / 'invalid.fa'
            path.write_text(text)
            with self.assertRaises(ValueError):
                list(aligner.iter_exon_query_batches(str(path), str(self.folder), 24))
            self.assertFalse((self.folder / 'exon_query_batch.fa').exists())

    def test_oversized_exon_fails_instead_of_splitting_or_exceeding_limit(self):
        path = self.folder / 'oversized.fa'
        path.write_text('>exon\nACGT\n')
        with self.assertRaisesRegex(ValueError, 'exon needs 11 FASTA bytes'):
            list(aligner.iter_exon_query_batches(str(path), str(self.folder), 10))
        self.assertFalse((self.folder / 'exon_query_batch.fa').exists())

    def args(self):
        return argparse.Namespace(
            blast_tabular=None, exon_fasta=str(self.query), db='reference',
            tmp_dir=str(self.folder), query='assembly.fa', output='alignments.tsv',
            makeblastdb='makeblastdb', blastn='blastn', threads=32, word_size=19,
            evalue='1e-30', blast_perc_identity=95, qcov_hsp_perc=None,
            max_target_seqs=100, blast_query_batch_bytes=24)

    def test_one_database_all_queries_and_requested_threads(self):
        commands, processes, identifiers = [], [], []

        def launch(command, **kwargs):
            if processes:
                processes[-1].wait.assert_called()
            commands.append(command)
            records = list(aligner.iter_query_records(command[command.index('-query') + 1]))
            ids = [header[1:].split()[0] for header, _ in records]
            identifiers.extend(ids)
            proc = Mock(stdout=io.StringIO(''.join(name + '\n' for name in ids)))
            proc.wait.return_value = proc.poll.return_value = 0
            processes.append(proc)
            return proc

        with patch.object(aligner.subprocess, 'run') as build, \
                patch.object(aligner.subprocess, 'Popen', side_effect=launch):
            output = list(aligner.iter_exon_query_lines(self.args()))
        build.assert_called_once()
        self.assertEqual(identifiers, ['E1.1', 'E2', 'long_exon', 'last'])
        self.assertEqual(output, [name + '\n' for name in identifiers])
        self.assertEqual(len(commands), 3)
        dbs = {cmd[cmd.index('-db') + 1] for cmd in commands}
        self.assertEqual(len(dbs), 1)
        for cmd in commands:
            self.assertEqual(cmd[cmd.index('-num_threads') + 1], '32')
        self.assertEqual(list(self.folder.iterdir()), [self.query])

    def test_internal_query_chunks_reach_blast_and_preserve_environment(self):
        for batch_bytes in (0, 24):
            for override in (None, '2000000'):
                with self.subTest(batch_bytes=batch_bytes, override=override):
                    args = self.args()
                    args.blast_query_batch_bytes = batch_bytes
                    environment = {'PATH': '/custom/blast/bin'}
                    if override is not None:
                        environment['BLAST_MT_QUERY_BATCH_SIZE'] = override
                    expected = dict(environment, BLAST_MT_QUERY_BATCH_SIZE=override or '1000000')

                    def launch(command, **kwargs):
                        self.assertEqual(kwargs['env'], expected)
                        proc = Mock(stdout=io.StringIO(''))
                        proc.wait.return_value = proc.poll.return_value = 0
                        return proc

                    with patch.dict(aligner.os.environ, environment, clear=True), \
                            patch.object(aligner.subprocess, 'run'), \
                            patch.object(aligner.subprocess, 'Popen', side_effect=launch) as blast, \
                            patch('sys.stderr', new_callable=io.StringIO) as log:
                        list(aligner.iter_exon_query_lines(args))
                        self.assertEqual(dict(aligner.os.environ), environment)
                    self.assertEqual(blast.call_count, 3 if batch_bytes else 1)
                    self.assertIn(
                        f"BLAST ThreadByQuery chunk size: {expected['BLAST_MT_QUERY_BATCH_SIZE']} bases",
                        log.getvalue(),
                    )

    def test_failed_batch_stops_before_later_queries(self):
        processes = []

        def launch(command, **kwargs):
            proc = Mock(stdout=io.StringIO('one row\n'))
            proc.wait.return_value = proc.poll.return_value = 2 if processes else 0
            processes.append(proc)
            return proc

        with patch.object(aligner.subprocess, 'run'), \
                patch.object(aligner.subprocess, 'Popen', side_effect=launch):
            with self.assertRaises(subprocess.CalledProcessError):
                list(aligner.iter_exon_query_lines(self.args()))
        self.assertEqual(len(processes), 2)
        self.assertEqual(list(self.folder.iterdir()), [self.query])


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb') and shutil.which('minimap2'),
                     'BLAST+ is required for the live equivalence check')
class LiveBlastBatchTests(unittest.TestCase):
    def test_batched_and_unbatched_alignments_and_both_call_tables_match(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            rng = random.Random(761)
            sequences = [''.join(rng.choices('ACGT', k=200)) for _ in range(5)]
            fixture = ReferenceFixture(folder)
            fixture.transcript('GA', 'Alpha', 'TA', sequences[:2])
            fixture.transcript('GB', 'Beta', 'TB', sequences[2:4])
            fixture.transcript('GC', 'Gamma', 'TC', [sequences[0], sequences[4]])
            fixture.write()
            full_gene = fixture.contigs[0][1]
            fragment = fixture.contigs[1][1][:208]
            assembly = folder / 'assembly.fa'
            assembly.write_text('>chrA\n' + full_gene + 'N' * 100 + full_gene +
                                '\n>chrB\n' + builder.revcomp(full_gene) +
                                '\n>fragment\n' + fragment + '\n')
            queries = folder / 'queries.txt'
            queries.write_text(f'sample {assembly}\n')
            outputs = []
            for batch_bytes, threads in ((0, 1), (250, 2)):
                output = folder / f'out_{batch_bytes}'
                command = [sys.executable, str(ROOT / 'annotate_assemblies.py'),
                           '--reference', str(fixture.genome), '--gff3', str(fixture.gff),
                           '--query-list', str(queries), '--exon-database-dir', str(folder / 'db'),
                           '--output', str(output), '--anchor-target-length', '14', '--min-unmasked', '0',
                           '--blast-threads', str(threads), '--caller-threads', '1',
                           '--blast-query-batch-bytes', str(batch_bytes)]
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                if batch_bytes:
                    self.assertIn('Local BLAST query batch 2:', result.stderr)
                    self.assertIn('-num_threads 1', result.stderr)
                    self.assertIn('up to 2 concurrent genes', result.stderr)
                call_dir = output / 'sample'
                calls = call_dir / 'sample.transcript_calls.tsv'
                fragments = next(call_dir.glob('*pseudo*.tsv'))
                self.assertTrue(rows(calls))
                self.assertTrue(rows(fragments))
                alignments = output / 'temp/sample/sample.exon_alignments.tsv'
                outputs.append((sorted(alignments.read_text().splitlines()),
                                calls.read_bytes(), fragments.read_bytes()))
            self.assertEqual(outputs[0], outputs[1])


if __name__ == '__main__':
    unittest.main()
