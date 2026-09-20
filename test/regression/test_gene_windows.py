"""Dynamic anchors and coarse-to-local search behavior, including real BLAST."""
import argparse
from dataclasses import replace
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT, rows
from test_overlap_query_merge import extracted
import align_exon_blastdb_v2 as a
import build_exon_blastdb_v2 as b
import blast_gene_windows as w


class WindowTests(unittest.TestCase):
    def test_sorted_merge_and_padding_on_both_sides(self):
        gene = w.Gene(100, 1101)
        self.assertEqual(gene.padding, 1502)
        hits = [(100, 200), (1200, 1300), (9000, 9100)]
        padded = [(max(0, start-gene.padding), min(10000, end+gene.padding)) for start, end in hits]
        self.assertEqual(list(w.merged_intervals(sorted(padded))), [(0, 2802), (7498, 10000)])

    def test_each_window_independent_all_exons_positive_equal(self):
        def counts(x, y):
            return sorted([(exon, '+', i*100, i*100+80)
                           for exon, n in [('A', x), ('B', y)] for i in range(n)])
        self.assertTrue(w.balanced_hits({'A', 'B'}, counts(1, 1)))
        self.assertTrue(w.balanced_hits({'A', 'B'}, counts(2, 2)))
        self.assertFalse(w.balanced_hits({'A', 'B'}, counts(1, 0)))
        self.assertFalse(w.balanced_hits({'A', 'B'}, counts(1, 2)))
        self.assertFalse(w.balanced_hits({'A', 'B'}, []))
        self.assertTrue(w.balanced_hits({'A', 'B'}, [
            ('A', '+', 0, 80), ('A', '+', 0, 80), ('A', '+', 10, 90), ('B', '-', 50, 130)]))

    def test_gene_name_fallback_and_alt_locus_separation(self):
        alias = a.ExonAlias('E1.1', 'E1', 'chr1', 10, 110, '+', 100,
                            gene_name='G', gene_id='ID1', gene_start0=1, gene_end0=1001)
        aliases = [alias, replace(alias, exon_id_full='E2.1', gene_id='ID2', start0=200, end0=300)]
        mapping = {'Q.1': aliases, 'Q': aliases,
                   'ALT': [replace(alias, chrom='alt')],
                   'UNNAMED': [replace(alias, gene_name='.')],
                   'OTHER': [replace(alias, gene_name='', gene_id='ID2')]}
        genes, query_genes = w.collect_genes(mapping)
        self.assertEqual(len(genes), 4)
        self.assertEqual(genes[0].exons, {'E1.1', 'E2.1'})
        self.assertEqual(list(genes[0].queries), ['Q.1'])
        self.assertEqual(genes[0].padding, 1500)

    def test_dynamic_anchor_lengths_edges_strand_and_union_offsets(self):
        rng = random.Random(906)
        genome = ''.join(rng.choices('ACGT', k=5000))
        for strand in ('+', '-'):
            with self.subTest(strand=strand), tempfile.TemporaryDirectory() as tmp:
                folder = Path(tmp)
                fasta = folder/'ref.fa'
                fasta.write_text('>chr1\n'+genome+'\n')
                specs = [('ONE', 100, 101), ('149', 500, 649), ('150', 900, 1050),
                         ('151', 1300, 1451), ('EDGE', 0, 50), ('LAST', 4960, 5000),
                         ('LONG', 3000, 3300), ('SHORT', 2980, 3060)]
                records = [extracted(name, lo, hi, genome, strand=strand).record for name, lo, hi in specs]
                prefix = str(folder/'db')
                b.write_exon_fasta(str(fasta), records, prefix+'.exons.fa', prefix+'.seq',
                                  prefix+'.exon_info.tsv', prefix+'.exon_aliases.tsv', 150, 0)
                aliases = {row['exon_id']: row for row in rows(Path(prefix+'.exon_aliases.tsv'))}
                expected = {'ONE': (75,75), '149': (1,1), '150': (0,0), '151': (0,0),
                            'EDGE': (0,50), 'LAST': (55,0), 'LONG': (0,0), 'SHORT': (35,35)}
                queries = {h[1:].split()[0]: seq for h, seq in a.iter_query_records(prefix+'.exons.fa')}
                for name, lo, hi in specs:
                    row = aliases[name]
                    anchors = expected[name] if strand == '+' else expected[name][::-1]
                    self.assertEqual((int(row['left_anchor_length']), int(row['right_anchor_length'])), anchors)
                    qseq = queries[row['blast_exon_id_full']]
                    core = qseq[int(row['query_core_start0']):int(row['query_core_end0'])]
                    self.assertEqual(core, genome[lo:hi] if strand == '+' else b.revcomp(genome[lo:hi]))
                    self.assertEqual(int(row['gene_start0']), 0)
                    self.assertEqual(int(row['gene_end0']), 5000)

    def test_extraction_failure_is_not_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = argparse.Namespace(blastdbcmd='blastdbcmd')
            with patch.object(w.subprocess, 'run', side_effect=subprocess.CalledProcessError(1, 'blastdbcmd')):
                with self.assertRaises(subprocess.CalledProcessError):
                    w.extract_windows(args, a, 'db', tmp, [('chr1', 0, 100)])


@unittest.skipUnless(all(shutil.which(cmd) for cmd in ('blastn','makeblastdb','blastdbcmd','minimap2')), 'BLAST+ required')
class LiveWindowTests(unittest.TestCase):
    def test_local_rescue_reverse_remapping_balanced_skip_and_unseeded_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            rng = random.Random(3129)
            genome = ''.join(rng.choices('ACGT', k=5000))
            fasta = folder/'ref.fa'; fasta.write_text('>chr1\n'+genome+'\n')
            specs = [('A1',100,400,'A'), ('A2',700,1000,'A'),
                     ('B1',2000,2300,'B'), ('B2',2600,2680,'B'),
                     ('C1',3500,3580,'C')]
            records = [extracted(name, lo, hi, genome, gene=gene).record for name,lo,hi,gene in specs]
            prefix = str(folder/'db')
            b.write_exon_fasta(str(fasta), records, prefix+'.exons.fa', prefix+'.seq',
                              prefix+'.exon_info.tsv', prefix+'.exon_aliases.tsv', 150, 0)
            a_copy = genome[:1100]
            mutant = list(a_copy)
            for pos in range(720, 1000, 35):
                mutant[pos] = next(x for x in 'ACGT' if x != mutant[pos])
            mutant = ''.join(mutant)
            b_copy = genome[1900:2800]
            target = folder/'assembly.fa'
            target.write_text('>NC_060925.1\n' + a_copy + 'N'*10000 + mutant +
                              '\n>reverse\n' + b.revcomp(mutant) +
                              '\n>short\n' + b_copy +
                              '\n>unseeded\n' + genome[3450:3630] + '\n')
            results = {}
            for mode in ('coarse', 'two_stage', 'minimap_serial', 'minimap_parallel'):
                out = folder/(mode+'.tsv')
                command = [sys.executable, str(ROOT/'align_exon_blastdb_v2.py'),
                           '--query', str(target), '--db', prefix, '--output', str(out),
                           '--exons-as-query', '--threads', '1' if mode == 'minimap_serial' else '2',
                           '--blast-query-batch-bytes','500']
                if not mode.startswith('minimap'):
                    command += ['--candidate-aligner', 'blast']
                if mode == 'coarse':
                    command += ['--no-local-realignment']
                result = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
                results[mode] = rows(out)
                if mode == 'two_stage':
                    self.assertIn('-word_size 50', result.stderr)
                    self.assertIn('-evalue 1e-100', result.stderr)
                    self.assertIn('-word_size 19', result.stderr)
                    self.assertNotIn("'num_threads' is currently ignored", result.stderr)
                    self.assertNotIn('-subject ', result.stderr)
                    self.assertIn('-evalue 1e-30', result.stderr)
                    self.assertIn('1 balanced (skip), 3 require local realignment', result.stderr)
                    self.assertIn('1 gene loci have no seed', result.stderr)
                    self.assertIn('Local realignment completed 2/2 genes', result.stderr)
                if mode.startswith('minimap'):
                    self.assertIn('First pass (minimap2)', result.stderr)
                    self.assertIn('balanced-window skipping disabled', result.stderr)
                    self.assertIn('Local realignment completed 3/3 genes', result.stderr)
                    self.assertIn('1 BLAST thread per gene', result.stderr)
                    self.assertNotIn('-evalue 1e-100', result.stderr)
                    self.assertIn('-num_threads 1', result.stderr)
                    self.assertNotIn('-num_threads 2', result.stderr)
                    self.assertIn('up to 1 concurrent genes' if mode == 'minimap_serial' else 'up to 2 concurrent genes', result.stderr)
            normalize = lambda data: sorted(tuple(sorted(row.items())) for row in data)
            self.assertEqual(normalize(results['minimap_serial']), normalize(results['minimap_parallel']))
            coarse = results['coarse']; final = results['two_stage']
            def coords(data, exon, contig):
                return sorted((int(r['query_start']), int(r['query_end']), r['strand'])
                              for r in data if r['exon_id']==exon and r['query_id']==contig)
            self.assertEqual(coords(coarse,'A2','NC_060925.1'), [(700,1000,'+')])
            self.assertEqual(coords(final,'A2','NC_060925.1'), [(700,1000,'+'), (11800,12100,'+')])
            self.assertEqual(coords(final,'A2','reverse'), [(100,400,'-')])
            self.assertFalse(coords(coarse,'B2','short'))
            self.assertEqual(coords(final,'B2','short'), [(700,780,'+')])
            self.assertFalse(any(r['exon_id']=='C1' for r in final))
            self.assertEqual(coords(results['minimap_parallel'], 'C1', 'unseeded'), [(50,130,'+')])
            self.assertEqual(coords(results['minimap_parallel'], 'B2', 'short'), [(700,780,'+')])
            self.assertEqual(coords(results['minimap_parallel'], 'A2', 'reverse'), [(100,400,'-')])
            # Refinement replaces coarse evidence rather than duplicating it.
            self.assertEqual(coords(final,'A1','NC_060925.1'), [(100,400,'+'), (11200,11500,'+')])
            self.assertEqual(len({tuple(sorted(row.items())) for row in final}), len(final))


if __name__ == '__main__':
    unittest.main()
