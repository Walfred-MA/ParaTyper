"""Hierarchical output and weighted gene edit alignment regressions."""
from dataclasses import replace
from pathlib import Path
import itertools
import random
import tempfile
import unittest

from test_exon_transcript_v3 import caller as c, make_call, write_fixture, cli, rows, bed
from test_copy_chaining import hit


def gene_call(hits, expected=None):
    base = make_call('G1_full_gene', mane=False)
    return replace(base, hits=hits, expected_exon_numbers=expected or list(range(1,9)),
                   transcript_info=replace(base.transcript_info, model_type='full_gene'))


class GeneEditTests(unittest.TestCase):
    def test_insertions_cost_50_instead_of_receiving_a_match_reward(self):
        for strand in ('+', '-'):
            for numbers in ([1,2,3,4,4,5,6,7,8], [1,2,3,4,2,5,6,7,8]):
                hits = [hit(n,i*100,i*100+20,strand=strand) for i,n in enumerate(numbers)]
                chain = c.best_gene_chain_for_hits(hits,strand)
                self.assertEqual([h.exon_number for h in chain],numbers)
                self.assertEqual([h.exon_number for h in chain if h.is_insertion],[numbers[4]])
                call = gene_call(chain)
                self.assertEqual(call.raw_score,750)
                self.assertEqual(call.insertion_penalty,50)
                self.assertEqual(call.found_expected_count,8)
                self.assertEqual(call.weighted_score,15000)
                row = dict(zip(c.OUTPUT_HEADER,c.call_to_output_row(call,False)))
                self.assertEqual(row['inserted_exons'],'1')
                self.assertEqual(row['insertion_penalty'],'50.000000')
                self.assertTrue(row['inserted_exon_query_coordinates'])

    def test_insertion_cost_is_applied_before_score_multipliers(self):
        numbers = [1,2,3,3,4,5,6,7]
        hits = [hit(n,i*100,i*100+20) for i,n in enumerate(numbers)]
        call = gene_call(c.best_gene_chain_for_hits(hits,'+'),list(range(1,8)))
        self.assertEqual(call.raw_score,650)
        self.assertEqual(call.insertion_penalty,50)
        self.assertEqual(call.found_expected_count,7)
        self.assertEqual(call.raw_score,sum(h.score for h in hits[3:])+150)

    def test_distinct_inserted_exons_are_charged_but_reference_gaps_are_free(self):
        numbers = [*range(1,7),2,3,*range(7,13)]
        hits = [hit(n,i*100,i*100+20) for i,n in enumerate(numbers)]
        call = gene_call(c.best_gene_chain_for_hits(hits,'+'),list(range(1,13)))
        self.assertEqual(call.raw_score,1100)
        self.assertEqual(call.insertion_penalty,100)
        self.assertEqual([h.exon_number for h in call.hits if h.is_insertion],[2,3])
        missing = [hit(1,0,20),hit(3,100,120),hit(5,200,220)]
        self.assertEqual(gene_call(c.best_gene_chain_for_hits(missing,'+')).raw_score,300)

    def test_tandem_copies_and_previously_extracted_blocks_remain_separate(self):
        for strand in ('+', '-'):
            hits = [hit(n,i*100,i*100+20,strand=strand)
                    for i,n in enumerate([1,2,3,1,2,3])]
            first = c.best_gene_chain_for_hits(hits,strand)
            self.assertEqual(len(first),3)
            self.assertFalse(any(h.is_insertion for h in first))
            used = {h.competition_interval for h in first}
            remaining = [h for h in hits if h.competition_interval not in used]
            self.assertEqual(len(c.best_gene_chain_for_hits(remaining,strand,all_hits=hits)),3)
            barrier = [hit(1,0,20,strand=strand),hit(2,100,120,strand=strand),hit(3,200,220,strand=strand)]
            self.assertEqual(len(c.best_gene_chain_for_hits([barrier[0],barrier[2]],strand,all_hits=barrier)),1)

    def test_only_insertions_do_not_count_as_matched_reference_exons(self):
        call = gene_call([hit(3,0,20),replace(hit(2,100,120),is_insertion=True),hit(4,200,220)], [2,3,4])
        self.assertEqual(call.found_expected_count,2)
        self.assertEqual(call.missing_count,1)
        self.assertEqual(call.raw_score,150)

    def test_removing_an_insertion_reranks_residual_gene_and_drops_insertion_only_calls(self):
        hits=[hit(1,0,20),hit(2,100,120),replace(hit(2,200,220),is_insertion=True),
              hit(3,300,320),hit(4,400,420)]
        gene=replace(gene_call(hits,[1,2,3,4]),protein_bonus=1,complete_bonus=1)
        blocker=replace(make_call('blocker',gene='GB',intervals=((200,220),)),complete_bonus=1)
        rival=replace(make_call('rival',gene='GC',intervals=((0,20),)),protein_bonus=3.75,complete_bonus=1)
        selected=c.resolve_call_overlaps([gene,blocker,rival],False)
        self.assertEqual({x.transcript_id for x in selected},{'G1_full_gene','blocker'})
        remaining=next(x for x in selected if x.transcript_id=='G1_full_gene')
        self.assertEqual(remaining.raw_score,400)
        self.assertEqual(remaining.insertion_penalty,0)
        all_matches=replace(make_call('matches',gene='GB'),hits=[h for h in hits if not h.is_insertion],
                            protein_bonus=10,complete_bonus=1)
        selected=c.resolve_call_overlaps([gene,all_matches],False)
        self.assertEqual([x.transcript_id for x in selected],['matches'])

    def test_repeated_exon_counts_once_per_run_and_again_in_later_runs(self):
        for strand in ('+','-'):
            for numbers, unique_counts in [([*range(1,7),2,2,3,2,*range(7,13)],[2]),
                                          ([1,2,3,4,2,2,5,6,7,8,2,2,9,10,11,12],[1,1])]:
                hits=[hit(n,i*100,i*100+20,strand=strand) for i,n in enumerate(numbers)]
                call=gene_call(c.best_gene_chain_for_hits(hits,strand),list(range(1,13)))
                self.assertEqual([len(run) for run in call.insertion_run_exons],unique_counts)
                self.assertEqual(call.inserted_exon_count,2)
                self.assertEqual(len(call.insertion_intervals),4)
                self.assertEqual(call.insertion_penalty,100)
                self.assertEqual(call.raw_score,1100)
                row=dict(zip(c.OUTPUT_HEADER,c.call_to_output_row(call,False)))
                self.assertEqual(row['inserted_exons'],'2')
                self.assertEqual(row['inserted_query_blocks'],'4')
                self.assertEqual(row['insertion_run_unique_exons'],','.join(map(str,unique_counts)))

    def test_cutoff_is_strictly_more_than_20_unique_exons_per_run(self):
        for strand in ('+','-'):
            for inserted in (list(range(20,0,-1)), list(range(21,0,-1)), [1]*45):
                numbers=[*range(101,126),*inserted,*range(126,151)]
                hits=[hit(n,i*100,i*100+20,strand=strand) for i,n in enumerate(numbers)]
                chain=c.best_gene_chain_for_hits(hits,strand)
                call=gene_call(chain,list(range(101,151)))
                self.assertTrue(call.valid_insertion_runs)
                if len(set(inserted))<=20:
                    self.assertEqual(len(chain),len(hits))
                    self.assertEqual(call.inserted_exon_count,len(set(inserted)))
                    self.assertEqual(call.raw_score,5000-50*len(set(inserted)))
                else:
                    self.assertFalse(chain[0].competition_interval==hits[0].competition_interval
                                     and chain[-1].competition_interval==hits[-1].competition_interval)

    def test_long_runs_are_pruned_independently_not_by_total_insertions(self):
        numbers=[*range(101,131),*range(20,0,-1),*range(131,161),*range(20,0,-1),*range(161,191)]
        hits=[hit(n,i*100,i*100+20) for i,n in enumerate(numbers)]
        call=gene_call(c.best_gene_chain_for_hits(hits,'+'),list(range(101,191)))
        self.assertEqual([len(run) for run in call.insertion_run_exons],[20,20])
        self.assertEqual(call.inserted_exon_count,40)
        self.assertTrue(call.valid_insertion_runs)
        self.assertEqual(call.raw_score,7000)

    def test_trimming_recounts_and_rejects_merged_insertion_runs_over_limit(self):
        # The separator is a match, so the repeated exon is charged twice.
        call=gene_call([hit(100,0,20),replace(hit(1,100,120),is_insertion=True),
                       hit(101,200,220),replace(hit(1,300,320),is_insertion=True),hit(102,400,420)])
        trimmed=call.recalculated_with_hits([h for h in call.hits if h.exon_number!=101])
        self.assertEqual(call.inserted_exon_count,2)
        self.assertEqual(trimmed.inserted_exon_count,1)
        # Removing a match can merge two individually valid runs into >20.
        numbers=[100,*range(1,12),101,*range(12,23),102]
        hits=[replace(hit(n,i*100,i*100+20),is_insertion=n<100) for i,n in enumerate(numbers)]
        candidate=replace(gene_call(hits,[100,101,102]),protein_bonus=1,complete_bonus=1)
        blocker=replace(make_call('blocker',gene='GB',intervals=((1200,1220),)),complete_bonus=1)
        selected=c.resolve_call_overlaps([candidate,blocker],False)
        self.assertEqual([x.transcript_id for x in selected],['blocker'])

    def test_weighted_edit_dp_matches_exhaustive_local_alignments(self):
        rng=random.Random(23019)
        for case in range(70):
            specs=[(rng.randrange(1,6),rng.randrange(5),rng.choice((0,40)),rng.choice((-250,60,90,100)))
                   for _ in range(8)]
            for strand in ('+','-'):
                hits=[hit(n,b*200+x,b*200+x+20,score,strand,block=(b*200,b*200+100)) for n,b,x,score in specs]
                # All aliases of one gene/merged interval share one score.
                scores={h.competition_interval:h.score for h in hits}
                hits=[replace(h,interval_score=replace(h.interval_score,score=scores[h.competition_interval])) for h in hits]
                ordered=sorted(hits,key=lambda h:(*c.oriented_hit_interval(h,strand),h.exon_number))
                intervals=sorted(scores,reverse=strand=='-')
                ranks={v:i for i,v in enumerate(intervals)}
                labels={interval:min((h for h in hits if h.competition_interval==interval),
                    key=lambda h:(c.interval_score_source_key(h.alignment),h.exon_number)).exon_number
                    for interval in intervals}
                best=(-float('inf'),0)
                for size in range(1,len(ordered)+1):
                    for chain in itertools.combinations(ordered,size):
                        if not all(a.exon_number<b.exon_number and c.oriented_hit_interval(a,strand)[1]<=c.oriented_hit_interval(b,strand)[0]
                                   for a,b in zip(chain,chain[1:])):
                            continue
                        blocks={h.competition_interval for h in chain}
                        positions=sorted({ranks[b] for b in blocks})
                        runs=[{labels[intervals[i]] for i in range(a+1,b)} for a,b in zip(positions,positions[1:])]
                        if any(len(run)>20 for run in runs):
                            continue
                        objective=sum(scores[b] for b in blocks)-50*sum(map(len,runs))
                        best=max(best,(objective,len(chain)))
                found=c.best_gene_chain_for_hits(hits,strand)
                actual=gene_call(found)
                self.assertEqual((actual.raw_score,len([h for h in found if not h.is_insertion])),best,(case,strand))


class HierarchyTests(unittest.TestCase):
    def test_gene_tie_pools_all_real_isoforms_but_never_ties_within_gene(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            gff,alignments=write_fixture(folder)
            ann=c.parse_gencode_gff3(str(gff))
            # Gene-stage tie remains; only GA has MANE at transcript stage.
            ann.transcripts['TC']=replace(ann.transcripts['TC'],is_mane=False)
            c.add_full_gene_transcripts(ann)
            calls,_=c.build_hierarchical_calls(c.read_alignments(str(alignments)),ann,
                c.truncate_expected_exons(ann,None),10,2,1,10,True)
            out=folder/'calls.tsv'
            c.write_transcript_calls(calls,str(out),'extended',False,False,False)
            data=rows(out)
            self.assertEqual(len(data),4)
            for parent,child in zip(data[::2],data[1::2]):
                self.assertEqual(parent['transcript_id'],'GA;GC')
                self.assertEqual(parent['gene_id'],'GA;GC')
                self.assertEqual(child['transcript_id'],'TA')
                self.assertEqual(child['gene_id'],'GA')
                self.assertEqual(child['tie_count'],'1')
                self.assertEqual(parent['GENE_index'],child['GENE_index'])
                metadata=dict(zip(bed.EXTRA_COLUMNS,bed.convert(child,'sample').fields[12:]))
                self.assertEqual(metadata['GENE_index'],parent['GENE_index'])

    def test_legacy_mode_rejects_hierarchy_before_writing(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            gff,alignments=write_fixture(folder)
            out=folder/'legacy.tsv'
            result=cli('-i',alignments,'-g',gff,'-o',out,'--output-format','legacy',check=False)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('two-stage',result.stderr)
            self.assertFalse(out.exists())


if __name__=='__main__':
    unittest.main()
