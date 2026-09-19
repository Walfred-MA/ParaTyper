"""Prevent one transcript chain from stitching together tandem gene copies."""
from dataclasses import replace
import itertools
import math
import random
import tempfile
from pathlib import Path
import unittest

from test_exon_transcript_v3 import ROOT, caller, cli, make_call, rows


def hit(number, start, end, score=100, strand='+', block=None):
    low, high = block or (start, end)
    if strand == '-':
        start, end = 1_000_000_000-end, 1_000_000_000-start
        low, high = 1_000_000_000-high, 1_000_000_000-low
    template = make_call('T', strand=strand, intervals=((start, end),)).hits[0]
    return replace(template, exon_number=number,
                   interval_score=caller.MergedIntervalGeneScore(low, high, 'G1', score))


class CopyChainingTests(unittest.TestCase):
    def test_scores_cannot_stitch_the_best_half_of_each_copy(self):
        for strand in ('+', '-'):
            first = [hit(i, 100*i, 100*i+20, s, strand) for i, s in enumerate((100,100,90),1)]
            second = [hit(i, 1000+100*i, 1020+100*i, s, strand) for i, s in enumerate((90,90,100),1)]
            hits = first+second
            windows = caller.exon_copy_windows(hits, strand)
            best = caller.best_chain_for_hits(hits, strand, windows)
            self.assertEqual([id(h) for h in best], [id(h) for h in first])
            remaining = [h for h in hits if id(h) not in {id(s) for s in best}]
            self.assertEqual(caller.best_chain_for_hits(remaining, strand, windows), second)

    def test_extracted_copy_remains_a_barrier_for_residual_fragments(self):
        for strand in ('+', '-'):
            hits = [hit(1,100,120,strand=strand), hit(1,300,320,strand=strand),
                    hit(2,400,420,strand=strand), hit(2,600,620,strand=strand)]
            windows = caller.exon_copy_windows(hits, strand)
            best = caller.best_chain_for_hits(hits, strand, windows)
            self.assertEqual(best, hits[1:3])
            self.assertEqual(len(caller.best_chain_for_hits([hits[0],hits[3]],strand,windows)),1)

    def test_long_intron_is_allowed_without_intervening_endpoint_copies(self):
        for strand in ('+', '-'):
            hits = [hit(1,100,120,strand=strand), hit(2,20_000_000,20_000_020,strand=strand),
                    hit(1,30_000_000,30_000_020,strand=strand)]
            self.assertEqual(caller.best_chain_for_hits(hits,strand), hits[:2])

    def test_overlapping_alternatives_are_not_copy_barriers(self):
        for strand in ('+', '-'):
            hits = [hit(1,100,140,100,strand,block=(100,170)),
                    hit(1,130,170,100,strand,block=(100,170)), hit(2,180,200,strand=strand)]
            chain = caller.best_chain_for_hits(hits,strand)
            self.assertEqual([id(h) for h in chain],[id(hits[0]),id(hits[2])])

    def test_dp_matches_exhaustive_chains_with_copy_barriers(self):
        rng = random.Random(314159)
        for _case in range(60):
            specifications = [(rng.randrange(1,5), rng.randrange(4), rng.choice((0,30,60)),
                               rng.choice((15,45))) for _ in range(8)]
            for strand in ('+', '-'):
                hits = [hit(n, b*200+s, b*200+s+length, 60+10*b, strand,
                            block=(b*200,b*200+120)) for n,b,s,length in specifications]
                def oriented(h):
                    a = h.alignment
                    return (-a.query_end,-a.query_start) if strand=='-' else (a.query_start,a.query_end)
                ordered = sorted(hits,key=lambda h:(*oriented(h),h.exon_number))
                def valid_link(a,b):
                    left_end,right_start=oriented(a)[1],oriented(b)[0]
                    return (a.exon_number < b.exon_number and left_end <= right_start
                            and not any(h.exon_number in (a.exon_number,b.exon_number)
                                        and left_end <= oriented(h)[0] and oriented(h)[1] <= right_start
                                        for h in hits))
                def objective(chain):
                    return sum({h.competition_interval:h.score for h in chain}.values()),len(chain)
                best=(-math.inf,0)
                for size in range(1,len(hits)+1):
                    for chain in itertools.combinations(ordered,size):
                        if all(valid_link(a,b) for a,b in zip(chain,chain[1:])):
                            best=max(best,objective(chain))
                self.assertEqual(objective(caller.best_chain_for_hits(hits,strand)),best)

    def test_real_c4_copies_are_complete_separate_and_serial_parallel_equal(self):
        data=ROOT.parent/'test/regression/data'
        with tempfile.TemporaryDirectory() as temporary:
            outputs=[]
            for threads in (1,2):
                out=Path(temporary)/f'calls{threads}.tsv'
                cli('-i',data/'C4_CHM13_tandem.alignments.tsv',
                    '-g',data/'C4_CHM13_tandem.gff3',
                    '--no-full-gene-transcripts',
                    '--eligible-exon-info',data/'C4_CHM13_tandem.exon_ids.tsv',
                    '--query-coordinate-mode','local','--threads',threads,'--shard-storage','disk','-o',out)
                outputs.append(out)
            self.assertEqual(outputs[0].read_bytes(),outputs[1].read_bytes())
            calls=rows(outputs[0])
            self.assertEqual(len(calls),2)
            actual={r['gene_name']:(int(r['query_start']),int(r['query_end']),r['found_exons'],r['call_status'])
                    for r in calls}
            self.assertEqual(actual,{'C4A':(31835262,31855887,'41','complete'),
                                     'C4B':(31868000,31888625,'41','complete')})
            self.assertTrue(all(r['found_exon_numbers']==','.join(map(str,range(1,42))) for r in calls))
            self.assertLess(actual['C4A'][1],actual['C4B'][0])


if __name__=='__main__':
    unittest.main()
