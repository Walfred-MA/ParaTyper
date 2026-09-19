"""Length resolves within-gene isoform ties, never between-gene ambiguity."""
from dataclasses import replace
import itertools
import tempfile
from pathlib import Path
import unittest

from test_exon_transcript_v3 import caller, make_call, resolve


class IsoformTieTests(unittest.TestCase):
    def test_longest_isoform_then_cross_gene_tie(self):
        calls = [make_call('short', transcript_length=200),
                 make_call('long', transcript_length=600),
                 make_call('paralog', gene='G2', transcript_length=1000)]
        for permutation in itertools.permutations(calls):
            out = resolve(permutation)
            self.assertEqual({c.transcript_id for c in out}, {'long', 'paralog'})
            self.assertEqual({c.tie_count for c in out}, {2})

    def test_mane_and_score_take_precedence_over_isoform_length(self):
        mane = make_call('mane', transcript_length=200)
        longer = make_call('longer', mane=False, transcript_length=1000)
        self.assertEqual([c.transcript_id for c in resolve([mane, longer])], ['mane'])
        self.assertEqual([c.transcript_id for c in resolve([mane, longer], prefer_mane=False)], ['longer'])
        lower_score = make_call('lower', identical=99, transcript_length=2000)
        self.assertEqual([c.transcript_id for c in resolve([mane, lower_score])], ['mane'])

    def test_disjoint_copies_of_same_gene_are_not_collapsed(self):
        left = make_call('short', transcript_length=200)
        right = make_call('long', intervals=((1000, 1100), (1300, 1400)), transcript_length=800)
        self.assertEqual({c.transcript_id for c in resolve([left, right])}, {'short', 'long'})

    def test_losing_isoform_can_retain_unoccupied_residual_evidence(self):
        longer = make_call('long', intervals=((100, 200), (300, 400)), transcript_length=800)
        shorter = make_call('short', intervals=((100, 200), (500, 600)), transcript_length=600)
        out = {c.transcript_id: c for c in resolve([shorter, longer])}
        self.assertEqual(set(out), {'long', 'short'})
        self.assertEqual(out['short'].found_exon_numbers, [2])
        self.assertEqual(out['short'].tie_count, 1)
        self.assertLess(out['short'].weighted_score, out['long'].weighted_score)

    def test_equal_lengths_choose_id_reproducibly_without_isoform_ambiguity(self):
        a, z = make_call('A'), make_call('Z')
        for order in ([a, z], [z, a]):
            out = resolve(order)
            self.assertEqual([c.transcript_id for c in out], ['A'])
            self.assertEqual(out[0].tie_count, 1)
            self.assertFalse(out[0].tie_group_id)
            self.assertEqual(len(caller.call_to_legacy_row(out[0], False)), 16)

    def test_spliced_annotation_length_does_not_include_introns_or_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'length.gff3'
            lines = ['##gff-version 3']
            for tid, coordinates in [('A', [(1,100), (10001,10100), (1,100)]),
                                     ('B', [(500,750)])]:
                for index, (start,end) in enumerate(coordinates, 1):
                    num = 1 if index == 3 else index
                    lines.append(f'chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\t'
                                 f'transcript_id={tid};gene_id=G;gene_name=G;transcript_type=protein_coding;'
                                 f'exon_number={num};exon_id={tid}E{num}')
            path.write_text('\n'.join(lines) + '\n')
            ann = caller.parse_gencode_gff3(str(path))
            self.assertEqual(ann.transcripts['A'].exon_length, 200)
            self.assertEqual(ann.transcripts['B'].exon_length, 251)
            a, b = make_call('A'), make_call('B')
            a = replace(a, transcript_info=ann.transcripts['A'])
            b = replace(b, transcript_info=ann.transcripts['B'])
            self.assertEqual([c.transcript_id for c in resolve([a, b])], ['B'])
