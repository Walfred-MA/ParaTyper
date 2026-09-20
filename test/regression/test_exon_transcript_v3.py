"""Regression tests for MANE-aware ambiguity reporting; no BLAST installation needed."""
import csv
from dataclasses import asdict, replace
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2] / 'scripts'
sys.path.insert(0, str(ROOT))
import call_genes_from_exon_alignments_v3 as caller
import annotate_assemblies as runner
import calls_to_bed as bed


def make_call(tid, gene='G1', intervals=((100, 200), (300, 400)), mane=True,
              identical=100, query='sample', strand='+', expected=None, transcript_length=None):
    info = caller.TranscriptInfo(tid, tid + '.1', gene, gene + '.1', gene,
                                 'protein_coding', mane, 100 * len(intervals) if transcript_length is None else transcript_length)
    hits = []
    for i, (start, end) in enumerate(intervals, 1):
        aln = caller.ExonAlignment(query, identical, start, end, strand,
                                   gene + 'E' + str(i), 100, 0, 100, 100,
                                   100.0, identical, 100.0)
        hits.append(caller.TranscriptHit(i, aln, 'chr1', i * 1000, i * 1000 + 100))
    return caller.TranscriptCall(tid, info, hits, expected or list(range(1, len(hits) + 1)),
                                 10.0, 2.0, 1.0)


def resolve(calls, **kwargs):
    return caller.resolve_call_overlaps(calls, **kwargs)


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def cli(*args, check=True):
    result = subprocess.run([sys.executable, str(ROOT / 'call_genes_from_exon_alignments_v3.py'),
                             *map(str, args)], text=True, capture_output=True)
    if check and result.returncode:
        raise AssertionError(result.stderr)
    return result


def write_fixture(folder):
    """Two MANE isoforms of one gene, a tied MANE paralog, and a non-MANE isoform."""
    gff = folder / 'fixture.gff3'
    text = ['##gff-version 3']
    for tid, gene, mane, exons in [('TA', 'GA', True, ('E1', 'E2')),
                                  ('TB', 'GA', True, ('E1', 'E2')),
                                  ('TC', 'GC', True, ('E3', 'E4')),
                                  ('TN', 'GA', False, ('E1', 'E2'))]:
        attrs = f'gene_id={gene}.1;transcript_id={tid}.1;gene_name={gene};transcript_type=protein_coding'
        if mane:
            attrs += ';tag=MANE_Select'
        text.append(f'chr1\ttest\ttranscript\t1001\t2100\t.\t+\t.\tID={tid}.1;{attrs}')
        for i, eid in enumerate(exons, 1):
            # These are competing paralogs at distinct reference loci, not
            # overlapping annotations of the same reference exon sequence.
            start = i * 1000 + (10000 if gene == 'GC' else 0)
            text.append(f'chr1\ttest\texon\t{start+1}\t{start+100}\t.\t+\t.\tParent={tid}.1;{attrs};exon_number={i};exon_id={eid}.1')
    gff.write_text('\n'.join(text) + '\n')
    alignments = folder / 'alignments.tsv'
    fields = list(caller.ExonAlignment.__dataclass_fields__)
    align_rows = []
    for query in ('sampleA', 'sampleB'):
        for gene, exons in [('GA', ('E1', 'E2')), ('GC', ('E3', 'E4'))]:
            call = make_call('dummy', gene=gene, query=query)
            for hit, eid in zip(call.hits, exons):
                align_rows.append(asdict(replace(hit.alignment, exon_id=eid)))
    with alignments.open('w') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(align_rows)
    aliases = folder / 'reference_exons.exon_aliases.tsv'
    aliases.write_text('exon_id\nE1\nE2\nE3\nE4\n')
    return gff, alignments


class RankingTests(unittest.TestCase):
    def test_isoforms_and_genes_tie_after_mane_priority(self):
        calls = [make_call('TA'), make_call('TB'), make_call('TC', gene='G2'),
                 make_call('TN', mane=False)]
        out = resolve(calls)
        self.assertEqual({c.transcript_id for c in out}, {'TA', 'TC'})
        self.assertEqual({c.tie_count for c in out}, {2})
        self.assertEqual(len({c.tie_group_id for c in out}), 1)
        self.assertTrue(out[0].tie_group_id)

    def test_mane_priority_and_explicit_opt_out(self):
        calls = [make_call('M', identical=96), make_call('N', mane=False)]
        self.assertEqual([c.transcript_id for c in resolve(calls)], ['M'])
        self.assertEqual([c.transcript_id for c in resolve(calls, prefer_mane=False)], ['N'])

    def test_equal_length_same_gene_isoforms_pick_stable_id(self):
        out = resolve([make_call('A', mane=False), make_call('B', mane=False)])
        self.assertEqual([c.transcript_id for c in out], ['A'])
        self.assertEqual(out[0].tie_count, 1)

    def test_unequal_scores_do_not_tie(self):
        out = resolve([make_call('A'), make_call('Z', identical=96)])
        self.assertEqual([c.transcript_id for c in out], ['A'])
        self.assertEqual(out[0].tie_count, 1)
        self.assertEqual(out[0].tie_group_id, '')

    def test_gene_block_maxima_cannot_discard_the_best_isoform(self):
        # Gene A has perfect short alternatives in both overlap blocks, but
        # its full MANE isoform is worse than gene B's. Pooling those maxima
        # and assigning a gene first used to incorrectly discard gene B.
        a = make_call('A_main', gene='GA', identical=98)
        b = make_call('B_main', gene='GB', identical=99)
        alternatives = []
        for i, start in enumerate((100, 300)):
            alt = make_call(f'A_alt{i}', gene='GA', mane=False,
                            intervals=((start, start + 20),))
            hit = alt.hits[0]
            aln = replace(hit.alignment, exon_length=20, exon_end=20,
                          aligned_exon_bases=20, identical_bases=20)
            alternatives.append(replace(alt, hits=[replace(hit, alignment=aln)]))
        out = resolve([a, b, *alternatives])
        self.assertEqual([c.transcript_id for c in out], ['B_main'])

    def test_mane_priority_precedes_cross_gene_elimination(self):
        out = resolve([make_call('M', gene='GA', identical=96),
                       make_call('N', gene='GB', mane=False)])
        self.assertEqual([c.transcript_id for c in out], ['M'])

    def test_separate_loci_and_strands_do_not_tie(self):
        out = resolve([make_call('A'), make_call('B', intervals=((1000, 1100), (1200, 1300))),
                       make_call('C', strand='-'), make_call('D', query='other')])
        self.assertEqual(len(out), 4)
        self.assertTrue(all(c.tie_count == 1 and c.tie_group_id == '' for c in out))

    def test_intron_nested_call_does_not_tie(self):
        out = resolve([make_call('A', intervals=((100, 200), (600, 700))),
                       make_call('B', intervals=((300, 400), (450, 550)))])
        self.assertEqual([c.tie_count for c in out], [1, 1])

    def test_tied_interval_union_blocks_lower_rank_nested_hits(self):
        calls = [make_call('A', intervals=((100, 300),)),
                 make_call('B', gene='G2', intervals=((150, 200),)),
                 make_call('C', intervals=((250, 290),), identical=96),
                 make_call('D', intervals=((400, 500),), identical=96)]
        out = resolve(calls)
        self.assertEqual({c.transcript_id for c in out}, {'A', 'B', 'D'})
        self.assertEqual({c.tie_count for c in out if c.transcript_id != 'D'}, {2})

    def test_residual_ties_after_stronger_call(self):
        calls = [make_call('top', intervals=((100, 200), (500, 600), (800, 900))),
                 make_call('A'), make_call('B', gene='G2')]
        out = {c.transcript_id: c for c in resolve(calls)}
        self.assertEqual(set(out), {'top', 'A', 'B'})
        self.assertEqual(out['A'].found_exon_numbers, [2])
        self.assertEqual(out['A'].tie_count, 2)
        self.assertEqual(out['A'].tie_group_id, out['B'].tie_group_id)

    def test_tie_ids_stable_across_input_order(self):
        calls = [make_call('A'), make_call('B', gene='G2'), make_call('C', gene='G3')]
        observed = set()
        for order in itertools.permutations(calls):
            observed.add(tuple(tuple(caller.call_to_output_row(c)) for c in resolve(order)))
        self.assertEqual(len(observed), 1)

    def test_transitive_ambiguity_component_and_adjacent_intervals(self):
        out = resolve([make_call('A', intervals=((100, 200),)),
                       make_call('B', gene='G2', intervals=((150, 250),)),
                       make_call('C', gene='G3', intervals=((200, 300),)),
                       make_call('D', intervals=((300, 400),))])
        by_id = {c.transcript_id: c for c in out}
        self.assertEqual([by_id[x].tie_count for x in 'ABC'], [3, 3, 3])
        self.assertEqual(by_id['D'].tie_count, 1)

    def test_tolerates_floating_roundoff(self):
        a, b = make_call('A'), make_call('B', gene='G2')
        b = replace(b, protein_bonus=10.0 + 1e-13)
        self.assertEqual([c.tie_count for c in resolve([a, b])], [2, 2])

    def test_legacy_does_not_silently_lose_ambiguity(self):
        out = resolve([make_call('A'), make_call('B', gene='G2')])
        with self.assertRaisesRegex(SystemExit, 'legacy output cannot represent tied'):
            caller.call_to_legacy_row(out[0], False)
        self.assertEqual(len(caller.call_to_legacy_row(make_call('A'), False)), 16)

    def test_bed_keeps_tie_metadata(self):
        call = resolve([make_call('A'), make_call('B', gene='G2')])[0]
        row = dict(zip(caller.OUTPUT_HEADER, caller.call_to_output_row(call)))
        record = bed.convert(row, 'h1')
        extras = dict(zip(bed.EXTRA_COLUMNS, record.fields[12:]))
        self.assertEqual(len(record.fields), 12 + len(bed.EXTRA_COLUMNS))
        self.assertEqual(extras['tie_group_id'], call.tie_group_id)
        self.assertEqual(extras['tie_count'], '2')
        self.assertEqual(extras['assignment_status'], 'tied')
        self.assertEqual(extras['pipeline_version'], '3.9.0')


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.gff, self.alignments = write_fixture(self.folder)
        self.blast_scratch = []

    def run_caller(self, name, threads=1, extra=()):
        out = self.folder / name
        cli('-i', self.alignments, '-g', self.gff, '-o', out, '--threads', threads,
            '--query-coordinate-mode', 'local', '--shard-storage', 'disk', *extra)
        return out

    def test_serial_parallel_and_opt_out(self):
        serial = self.run_caller('serial.tsv')
        parallel = self.run_caller('parallel.tsv', 2)
        self.assertEqual(serial.read_bytes(), parallel.read_bytes())
        data = rows(serial)
        self.assertEqual(len(data), 4)
        self.assertEqual({r['transcript_id'] for r in data}, {'GA;GC', 'TA;TC'})
        self.assertEqual({r['model_type'] for r in data}, {'full_gene;full_gene', 'transcript;transcript'})
        self.assertEqual(len({r['GENE_index'] for r in data}), 2)
        self.assertEqual({r['tie_count'] for r in data}, {'2'})
        self.assertEqual(len({r['tie_group_id'] for r in data}), 4)
        no_mane = rows(self.run_caller('no_mane.tsv', extra=('--no-prefer-mane',)))
        self.assertEqual(len(no_mane), 4)
        self.assertEqual({r['tie_count'] for r in no_mane}, {'2'})

    def test_smn_real_fixture_serial_parallel(self):
        outputs = []
        for threads in (1, 2):
            out = self.folder / f'smn{threads}.tsv'
            cli('-i', ROOT.parent / 'test/regression/data/SMN_align_CHM13.tsv',
                '-g', ROOT.parent / 'test/regression/data/SMN.gff3', '-o', out,
                '--threads', threads, '--shard-storage', 'disk',
                '--query-coordinate-mode', 'local')
            outputs.append(out)
        self.assertEqual(outputs[0].read_bytes(), outputs[1].read_bytes())
        data = rows(outputs[0])
        self.assertTrue(data)
        self.assertTrue(all(set(r['gene_name'].split(';')) <= {'SMN1', 'SMN2'} for r in data))
        self.assertTrue(all(r['pipeline_version'] == '3.9.0' for r in data))

    def test_chm13_longest_exon_scores_separate_main_mane_calls(self):
        data_dir = ROOT.parent / 'test/regression/data'
        outputs = []
        for threads in (1, 2):
            out = self.folder / f'smn_gene_assignment_{threads}.tsv'
            cli('-i', data_dir / 'SMN_CHM13_gene_filter.alignments.tsv',
                '-g', data_dir / 'SMN_CHM13_gene_filter.gff3', '-o', out,
                '--no-full-gene-transcripts',
                '--eligible-exon-info', data_dir / 'SMN_CHM13_gene_filter.exon_ids.tsv',
                '--threads', threads, '--shard-storage', 'disk',
                '--query-coordinate-mode', 'local')
            outputs.append(out)
        self.assertEqual(outputs[0].read_bytes(), outputs[1].read_bytes())
        mane = [r for r in rows(outputs[0]) if r['ifmane_transcript'] == '1']
        self.assertEqual(
            {(r['gene_name'], r['transcript_id'], r['query_start'], r['strand']) for r in mane},
            {('SMN1', 'ENST00000380707', '71381874', '+'),
             ('SMN2', 'ENST00000380743', '70809743', '-')},
        )
        self.assertTrue(all(r['found_exons'] == '9' and r['call_status'] == 'complete'
                            for r in mane))
        positive = [r for r in mane if r['strand'] == '+']
        negative = [r for r in mane if r['strand'] == '-']
        self.assertEqual({r['tie_count'] for r in positive}, {'1'})
        self.assertEqual(len({r['tie_group_id'] for r in positive}), 1)
        self.assertEqual(len({r['weighted_score'] for r in positive}), 1)
        self.assertEqual(negative[0]['assignment_status'], 'unique')
        self.assertAlmostEqual(float(positive[0]['weighted_score']), 15960.199005, places=5)
        self.assertAlmostEqual(float(negative[0]['weighted_score']), 15841.095017, places=5)
        self.assertTrue(all(len(r['merged_interval_coordinates'].split(',')) == 8
                            for r in mane))

    def runner_args(self):
        query = self.folder / 'query.fa'
        query.write_text('>sampleA\nACGT\n>sampleB\nACGT\n>empty\nACGT\n')
        argv = ['annotate_assemblies.py', '-r', str(query), '-g', str(self.gff),
                '--query-fasta', str(query), '-d', str(self.folder),
                '-o', str(self.folder / 'out'), '--caller-threads', '2']
        with patch.object(sys, 'argv', argv):
            return runner.parse_args()

    def create_fake_blast_database(self, command):
        scratch = Path(command[command.index('--tmp-dir') + 1])
        self.assertTrue(scratch.is_dir())
        database = scratch / 'exon_blast_assembly_test' / 'assembly.nsq'
        database.parent.mkdir()
        database.write_bytes(b'temporary assembly database')
        self.blast_scratch.append(scratch)

    def fake_alignment_real_caller(self, command, label):
        if Path(command[1]).name == 'align_exon_blastdb_v2.py':
            self.create_fake_blast_database(command)
            Path(command[command.index('--output') + 1]).write_bytes(self.alignments.read_bytes())
        else:
            self.assertTrue(self.blast_scratch)
            self.assertTrue(all(not path.exists() for path in self.blast_scratch))
            result = subprocess.run(command, text=True, capture_output=True)
            if result.returncode:
                raise AssertionError(result.stderr)

    def test_combined_runner_splits_all_ties_and_empty_record(self):
        args = self.runner_args()
        scripts = runner.load_scripts(ROOT)
        with patch.object(runner, 'run_command', side_effect=self.fake_alignment_real_caller) as run:
            paths = runner.run_combined_fasta(args.query_fasta, ['sampleA', 'sampleB', 'empty'],
                                             args, scripts, self.folder / 'reference_exons',
                                             args.output, args.output / 'temp')
            self.assertEqual(run.call_count, 2)
            self.assertEqual([len(rows(p)) for p in paths], [2, 2, 0])
            self.assertTrue(all(runner.valid_call_table(p) for p in paths))
            self.assertEqual({r['tie_count'] for r in rows(paths[0])}, {'2'})
            self.assertEqual((args.output / 'temp/combined_fasta/combined.exon_alignments.tsv').read_bytes(),
                             self.alignments.read_bytes())
            self.assertTrue(args.query_fasta.is_file())
            runner.run_combined_fasta(args.query_fasta, ['sampleA', 'sampleB', 'empty'],
                                     args, scripts, self.folder / 'reference_exons',
                                     args.output, args.output / 'temp')
            self.assertEqual(run.call_count, 2)

    def test_named_assembly_runner_and_mane_option(self):
        args = self.runner_args()
        scripts = runner.load_scripts(ROOT)
        args.prefer_mane = False
        with patch.object(runner, 'run_command', side_effect=self.fake_alignment_real_caller):
            output = runner.run_one_sample(runner.AssemblyQuery('assembly', args.query_fasta),
                                           args, scripts, self.folder / 'reference_exons',
                                           args.output, args.output / 'temp')
        self.assertEqual(len(rows(output)), 4)
        self.assertEqual({r['tie_count'] for r in rows(output)}, {'2'})
        self.assertEqual((args.output / 'temp/assembly/assembly.exon_alignments.tsv').read_bytes(),
                         self.alignments.read_bytes())
        self.assertTrue(args.query_fasta.is_file())

    def test_alignment_failure_cleans_only_its_own_blast_scratch(self):
        args = self.runner_args()
        scripts = runner.load_scripts(ROOT)
        database_prefix = self.folder / 'reference_exons'
        shared_database = database_prefix.with_suffix('.nsq')
        shared_database.write_bytes(b'shared reference database')

        def fail_alignment(command, label):
            self.create_fake_blast_database(command)
            raise subprocess.CalledProcessError(1, command)

        for combined in (False, True):
            with self.subTest(combined=combined):
                sample_temp = args.output / 'temp' / ('combined_fasta' if combined else 'assembly')
                unrelated = sample_temp / 'blast_tmp_other_run' / 'assembly.nsq'
                unrelated.parent.mkdir(parents=True)
                unrelated.write_bytes(b'another run owns this database')
                with patch.object(runner, 'run_command', side_effect=fail_alignment) as run:
                    with self.assertRaises(subprocess.CalledProcessError):
                        if combined:
                            runner.run_combined_fasta(args.query_fasta, ['sampleA', 'sampleB', 'empty'],
                                                      args, scripts, database_prefix,
                                                      args.output, args.output / 'temp')
                        else:
                            runner.run_one_sample(runner.AssemblyQuery('assembly', args.query_fasta),
                                                  args, scripts, database_prefix,
                                                  args.output, args.output / 'temp')
                self.assertEqual(run.call_count, 1)
                self.assertTrue(all(not path.exists() for path in self.blast_scratch))
                self.assertEqual(unrelated.read_bytes(), b'another run owns this database')
                self.assertEqual(shared_database.read_bytes(), b'shared reference database')
                self.assertTrue(args.query_fasta.is_file())

    def test_old_output_is_not_reused(self):
        old = self.folder / 'old.tsv'
        old.write_text('query_contig\ttranscript_id\tcall_status\nsample\tTA\tcomplete\n')
        self.assertFalse(runner.valid_call_table(old))
        self.assertTrue(runner.valid_call_table(self.run_caller('new.tsv')))

    def test_previous_gene_filtered_v3_output_is_not_reused(self):
        output = self.run_caller('previous_version.tsv')
        data = rows(output)
        for row in data:
            row['pipeline_version'] = '3.0.0'
        with output.open('w') as handle:
            writer = csv.DictWriter(handle, caller.OUTPUT_HEADER, delimiter='\t')
            writer.writeheader()
            writer.writerows(data)
        self.assertFalse(runner.valid_call_table(output))

    def test_empty_tables_with_removed_columns_are_not_reused(self):
        output = self.folder / 'old_empty.tsv'
        for removed in ('transcript_id_full', 'gene_id_full', 'alternatives_json'):
            output.write_text('\t'.join(caller.OUTPUT_HEADER + [removed]) + '\n')
            self.assertFalse(runner.valid_call_table(output))
        output.write_text('\t'.join(caller.OUTPUT_HEADER) + '\n')
        self.assertTrue(runner.valid_call_table(output))

    def test_no_full_gene_liftover_entrypoints_or_dependencies(self):
        for script in ('annotate_assemblies.py', 'install.py'):
            result = subprocess.run([sys.executable, str(ROOT / script), '--help'],
                                    text=True, capture_output=True, check=True)
            for removed in ('--lift-genome', '--with-full-gene',
                            '--samtools', '--stretcher'):
                self.assertNotIn(removed, result.stdout)
        self.assertFalse((ROOT / 'GeneLiftover.py').exists())
        self.assertFalse((ROOT / 'GeneGlobalAlign.py').exists())


if __name__ == '__main__':
    unittest.main()
