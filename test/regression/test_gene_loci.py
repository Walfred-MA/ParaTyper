"""Gene chains cannot collect complementary exon fragments across distant loci."""
from dataclasses import replace
import unittest

from test_exon_transcript_v3 import caller as c, make_call
from test_copy_chaining import hit


class GeneLocusTests(unittest.TestCase):
    def test_disconnected_windows_cannot_form_one_chain_on_either_strand(self):
        for strand in ('+', '-'):
            hits = [hit(1, 16_173_421, 16_173_521, strand=strand),
                    hit(2, 119_875_251, 119_875_393, strand=strand),
                    hit(3, 144_542_404, 144_542_504, strand=strand)]
            remaining = hits[:]
            calls = []
            while remaining:
                chain = c.best_gene_chain_for_hits(
                    remaining, strand, all_hits=hits, max_gap=407_176)
                calls.append(chain)
                used = {h.competition_interval for h in chain}
                remaining = [h for h in remaining if h.competition_interval not in used]
            self.assertEqual([len(chain) for chain in calls], [1, 1, 1])

    def test_touching_padded_windows_join_but_disconnected_windows_do_not(self):
        for strand in ('+', '-'):
            for gap, expected in ((300, 2), (301, 1)):
                hits = [hit(1, 100, 120, strand=strand),
                        hit(2, 120 + gap, 140 + gap, strand=strand)]
                chain = c.best_gene_chain_for_hits(hits, strand, max_gap=300)
                self.assertEqual(len(chain), expected)

    def test_insertions_still_work_inside_a_locus(self):
        for strand in ('+', '-'):
            hits = [hit(n, i * 100, i * 100 + 20, strand=strand)
                    for i, n in enumerate((1, 2, 3, 3, 4, 5, 6))]
            hits.append(hit(7, 100_000_000, 100_000_020, strand=strand))
            chain = c.best_gene_chain_for_hits(hits, strand, max_gap=1000)
            self.assertEqual([h.exon_number for h in chain], [1, 2, 3, 3, 4, 5, 6])
            self.assertEqual(sum(h.is_insertion for h in chain), 1)

    @staticmethod
    def fixture(reference_starts, target_intervals, strand):
        template = make_call('T', intervals=target_intervals, strand=strand)
        exons = []
        associations = {}
        for number, start in enumerate(reference_starts, 1):
            eid = f'G1E{number}'
            exon = c.GffExon(eid, eid + '.1', 'T', 'T.1', number,
                             'G1', 'G1.1', 'G1', 'protein_coding', 'chr1',
                             start, start + 100, '+', False)
            exons.append(exon)
            associations[eid] = [('T', number, 'chr1', start, start + 100)]
        ann = c.Annotation({'T': template.transcript_info}, {'T': exons}, associations)
        c.add_full_gene_transcripts(ann)
        alignments = [h.alignment for h in template.hits]
        if strand == '-':
            alignments = [replace(a, query_start=1_000_000_000 - a.query_end,
                                  query_end=1_000_000_000 - a.query_start) for a in alignments]
        return ann, alignments

    def test_hierarchy_keeps_distant_fragments_in_separate_parents(self):
        for strand in ('+', '-'):
            ann, alignments = self.fixture(
                [0, 135_625], [(16_173_421, 16_173_521), (144_542_404, 144_542_504)], strand)
            calls, _ = c.build_hierarchical_calls(
                alignments, ann, c.truncate_expected_exons(ann, None), 10, 2, 1, 10, True)
            parents = [call for call in calls if call.transcript_info.model_type == 'full_gene']
            self.assertEqual(len(parents), 2)
            self.assertTrue(all(call.found_expected_count == 1 for call in parents))
            self.assertEqual(len({call.gene_index for call in parents}), 2)
            self.assertTrue(all(c.call_query_bounds(call)[1] - c.call_query_bounds(call)[0] == 100
                                for call in calls))

    def test_gene_span_includes_unaligned_exons_and_introns(self):
        # The third exon is absent from the eligible database, but the full
        # 2 Mb reference span must still allow the first two hits to connect.
        for strand in ('+', '-'):
            ann, alignments = self.fixture(
                [0, 1000, 2_000_000], [(100, 200), (3_000_000, 3_000_100)], strand)
            expected = c.truncate_expected_exons(ann, None, eligible_exon_ids={'G1E1', 'G1E2'})
            calls, _ = c.build_hierarchical_calls(alignments, ann, expected, 10, 2, 1, 10, True)
            parents = [call for call in calls if call.transcript_info.model_type == 'full_gene']
            self.assertEqual(len(parents), 1)
            self.assertEqual(parents[0].found_expected_count, 2)


if __name__ == '__main__':
    unittest.main()
