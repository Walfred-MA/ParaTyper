"""Behavioral tests for merged interval -> per-gene score -> isoform DP."""
from dataclasses import replace
import itertools
import unittest

from test_exon_transcript_v3 import caller, make_call


def with_interval(call, start, end, score):
    evidence = caller.MergedIntervalGeneScore(start, end, call.transcript_info.gene_id, score)
    return replace(call, hits=[replace(h, interval_score=evidence) for h in call.hits])


class MergedIntervalTests(unittest.TestCase):
    def test_equal_length_best_score_is_shared_across_isoforms_but_not_genes(self):
        a1 = make_call('A1', gene='GA', identical=98)
        a2 = make_call('A2', gene='GA', identical=100,
                       intervals=((150, 250), (350, 450)))
        b = make_call('B', gene='GB', identical=99)
        calls = [a1, a2, b]
        grouped = {(c.query_id, c.strand, c.transcript_id): c.hits for c in calls}
        ann = caller.Annotation({c.transcript_id: c.transcript_info for c in calls}, {}, {})
        caller.assign_merged_interval_gene_scores(grouped, ann)
        self.assertIs(a1.hits[0].interval_score, a2.hits[0].interval_score)
        self.assertEqual(a1.hits[0].competition_interval, (100, 250))
        self.assertEqual([h.score for h in a1.hits], [100, 100])
        self.assertEqual([h.score for h in b.hits], [96, 96])
        self.assertIsNot(a1.hits[0].interval_score, b.hits[0].interval_score)
        # No gene candidate was removed by constructing the shared intervals.
        self.assertEqual(len(grouped), 3)

    def test_longest_exon_supplies_score_even_when_shorter_exon_is_perfect(self):
        long = make_call('long', identical=98)
        short = make_call('short', intervals=((100, 150), (300, 350)))
        short = replace(short, hits=[replace(h, alignment=replace(
            h.alignment, exon_id='SHORT', exon_length=50, identical_bases=50,
            aligned_exon_bases=50, exon_end=50)) for h in short.hits])
        calls = [short, long]
        grouped = {(c.query_id, c.strand, c.transcript_id): c.hits for c in calls}
        ann = caller.Annotation({c.transcript_id: c.transcript_info for c in calls}, {}, {})
        caller.assign_merged_interval_gene_scores(grouped, ann)
        self.assertEqual([h.score for h in short.hits], [92, 92])
        self.assertEqual(short.hits[0].interval_score.source_exon_length, 100)
        self.assertNotEqual(short.hits[0].interval_score.source_exon_id, 'SHORT')

    def test_interval_counted_once_for_multiple_exons(self):
        call = with_interval(make_call('T', intervals=((100, 150), (200, 250))),
                             100, 250, 90)
        chain = caller.best_chain_for_hits(call.hits, '+')
        self.assertEqual([h.exon_number for h in chain], [1, 2])
        self.assertEqual(call.raw_score, 90)
        self.assertEqual(call.weighted_score, 1800)

    def test_same_gene_score_tie_prefers_longer_isoform(self):
        two = with_interval(make_call('two', intervals=((100, 150), (200, 250))),
                            100, 250, 90)
        one = with_interval(make_call('one', intervals=((100, 250),)), 100, 250, 90)
        out = caller.resolve_call_overlaps([two, one])
        self.assertEqual([c.transcript_id for c in out], ['two'])
        self.assertEqual(out[0].tie_count, 1)

    def test_selection_consumes_merged_interval_not_just_exon_span(self):
        a = with_interval(make_call('A', gene='GA', intervals=((100, 120),)),
                          100, 250, 100)
        b = with_interval(make_call('B', gene='GB', intervals=((200, 220),)),
                          100, 250, 95)
        out = caller.resolve_call_overlaps([a, b])
        self.assertEqual([c.transcript_id for c in out], ['A'])
        b = with_interval(b, 100, 250, 100)
        out = caller.resolve_call_overlaps([a, b])
        self.assertEqual({c.transcript_id for c in out}, {'A', 'B'})
        self.assertEqual({c.tie_count for c in out}, {2})

    def test_interval_scores_not_multiplied_by_gene_length(self):
        a = with_interval(make_call('short', intervals=((100, 150),)), 100, 500, 99)
        b = with_interval(make_call('long', gene='GB', intervals=((100, 500),)),
                          100, 500, 98)
        self.assertEqual([c.transcript_id for c in caller.resolve_call_overlaps([a, b])],
                         ['short'])

    def test_dp_matches_exhaustive_valid_chains_on_both_strands(self):
        specifications = [(1,100,130,0), (3,130,160,0), (2,140,170,0),
                          (3,300,330,1), (4,330,350,1), (2,350,380,1),
                          (5,500,530,2), (4,540,560,2)]
        blocks = [(100,200,50), (300,400,60), (500,600,70)]
        for strand in ('+', '-'):
            hits = []
            for exon, start, end, bi in specifications:
                low, high, score = blocks[bi]
                if strand == '-':
                    start, end = 1000-end, 1000-start
                    low, high = 1000-high, 1000-low
                hit = make_call('T', intervals=((start,end),), strand=strand).hits[0]
                hits.append(replace(hit, exon_number=exon,
                                    interval_score=caller.MergedIntervalGeneScore(low,high,'G1',score)))
            def oriented(h):
                a = h.alignment
                return (-a.query_end,-a.query_start) if strand == '-' else (a.query_start,a.query_end)
            ordered = sorted(hits, key=lambda h: (*oriented(h),h.exon_number))
            def objective(chain):
                distinct = {h.competition_interval:h.score for h in chain}
                return sum(distinct.values()), len(chain)
            best = (0,0)
            for size in range(1,len(hits)+1):
                for chain in itertools.combinations(ordered,size):
                    if all(a.exon_number < b.exon_number and oriented(a)[1] <= oriented(b)[0]
                           for a,b in zip(chain,chain[1:])):
                        best = max(best,objective(chain))
            actual = caller.best_chain_for_hits(hits,strand)
            self.assertEqual(objective(actual),best,strand)


if __name__ == '__main__':
    unittest.main()
