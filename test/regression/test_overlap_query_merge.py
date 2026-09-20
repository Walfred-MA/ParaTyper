"""Same-gene unions shrink BLAST input while retaining per-exon evidence."""
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

from test_exon_transcript_v3 import ROOT, rows
import build_exon_blastdb_v2 as b
import align_exon_blastdb_v2 as a


def extracted(identifier, start, end, genome, gene='G', chrom='chr1', strand='+', anchor=60):
    rec = b.ExonRecord(chrom, start, end, strand, identifier + '.1', identifier,
                      [identifier + 'T.1'], [identifier + 'T'], ['1'],
                      [gene + '.1'], [gene], [gene], ['protein_coding'], False)
    lo, hi = max(0, start - anchor), min(len(genome), end + anchor)
    sequence = genome[lo:hi]
    left, right = start - lo, hi - end
    if strand == '-':
        sequence = b.revcomp(sequence)
        left, right = right, left
    return b.ExtractedExon(rec, sequence, b.unmasked_acgt_count(sequence), left, right)


class OverlapQueryTests(unittest.TestCase):
    def setUp(self):
        rng = random.Random(372)
        self.genome = ''.join(rng.choices('ACGT', k=13000))

    def test_user_eight_exons_become_three_queries_with_all_aliases(self):
        spans = [(11120,11211), (11124,11211), (11409,11671), (11410,11671),
                 (11769,11844), (11818,11844), (11824,11844), (11827,11844)]
        exons = [extracted(f'E{i}', start, end, self.genome, gene='DDX11L1&DDX11L16')
                 for i, (start, end) in enumerate(spans)]
        queries = b.merge_overlapping_exon_queries(exons)
        self.assertEqual([(q.target.record.start0, q.target.record.end0) for q in queries],
                         [(11120,11211), (11409,11671), (11769,11844)])
        self.assertEqual(sum(len(x.sequence) for x in exons), 1799)
        self.assertEqual(sum(len(q.target.sequence) for q in queries), 788)
        self.assertEqual([len(q.aliases) for q in queries], [2, 2, 4])
        self.assertEqual({x.exon.record.exon_id for q in queries for x in q.aliases},
                         {f'E{i}' for i in range(8)})
        for query in queries:
            for alias in query.aliases:
                self.assertEqual(query.target.sequence[alias.anchor_start0:alias.anchor_end0],
                                 alias.exon.sequence)

    def test_partial_overlap_unions_both_ends_and_reverse_strand_offsets(self):
        for strand in ('+', '-'):
            exons = [extracted('E1',100,200,self.genome,strand=strand),
                     extracted('E2',180,260,self.genome,strand=strand),
                     extracted('E3',250,320,self.genome,strand=strand)]
            query, = b.merge_overlapping_exon_queries(exons)
            self.assertEqual((query.target.record.start0,query.target.record.end0),(100,320))
            expected = self.genome[40:380]
            self.assertEqual(query.target.sequence, b.revcomp(expected) if strand=='-' else expected)
            for alias in query.aliases:
                expected = self.genome[alias.exon.record.start0:alias.exon.record.end0]
                self.assertEqual(query.target.sequence[alias.core_start0:alias.core_end0],
                                 b.revcomp(expected) if strand=='-' else expected)

    def test_different_genes_contigs_strands_and_flank_only_overlap_stay_separate(self):
        exons = [extracted('E1',100,200,self.genome),
                 extracted('E2',120,190,self.genome,gene='OTHER'),
                 extracted('E3',120,190,self.genome,chrom='chr2'),
                 extracted('E4',120,190,self.genome,strand='-'),
                 extracted('E5',200,250,self.genome),
                 extracted('E6',260,300,self.genome)]
        self.assertEqual(len(b.merge_overlapping_exon_queries(exons)),6)

    def test_contig_edge_clips_anchors_without_losing_original_metadata(self):
        exons = [extracted('E1',5,100,self.genome), extracted('E2',10,110,self.genome)]
        query, = b.merge_overlapping_exon_queries(exons)
        self.assertEqual(query.target.left_anchor_length,5)
        self.assertEqual([x.exon.record.transcript_ids for x in query.aliases],[['E1T'],['E2T']])
        for alias in query.aliases:
            self.assertEqual(alias.anchor_start0,0)
            self.assertEqual(query.target.sequence[alias.anchor_start0:alias.anchor_end0],alias.exon.sequence)

    def test_same_name_merges_different_gene_ids_and_preserves_alias_metadata(self):
        exons = [extracted('E1',100,200,self.genome,gene='G1'),
                 extracted('E2',180,260,self.genome,gene='G2')]
        for exon in exons:
            exon.record.gene_names = ['SYMBOL']
        query, = b.merge_overlapping_exon_queries(exons)
        self.assertEqual((query.target.record.start0,query.target.record.end0),(100,260))
        self.assertEqual(set(query.target.record.gene_ids),{'G1','G2'})
        self.assertEqual([x.exon.record.gene_ids for x in query.aliases],[['G1'],['G2']])
        self.assertEqual([x.exon.record.transcript_ids for x in query.aliases],[['E1T'],['E2T']])
        for alias in query.aliases:
            self.assertEqual(query.target.sequence[alias.anchor_start0:alias.anchor_end0],
                             alias.exon.sequence)
        # The query identity follows the common name, independent of ID order.
        query_id = query.target.record.exon_id
        exons[0].record.gene_ids, exons[1].record.gene_ids = ['G2'], ['G1']
        reordered, = b.merge_overlapping_exon_queries(list(reversed(exons)))
        self.assertEqual(reordered.target.record.exon_id,query_id)

    def test_gene_name_takes_precedence_over_id(self):
        exons = [extracted('E1',100,200,self.genome),
                 extracted('E2',180,260,self.genome)]
        exons[0].record.gene_names = ['FIRST']
        exons[1].record.gene_names = ['SECOND']
        self.assertEqual(len(b.merge_overlapping_exon_queries(exons)),2)

    def test_missing_names_fall_back_to_ids_without_merging_unrelated_records(self):
        exons = [extracted(f'E{i}',100+i,200+i,self.genome,gene=gene)
                 for i,gene in enumerate(['G1','G1','G2','G3','G4','G5'])]
        for exon in exons:
            exon.record.gene_names = []
        exons[1].record.gene_names = ['', '.', ' ']
        # A name matching another record's ID must not collide with its fallback.
        exons[3].record.gene_names = ['G2']
        exons[4].record.gene_ids = []
        exons[5].record.gene_ids = []
        groups = {frozenset(x.exon.record.exon_id for x in q.aliases)
                  for q in b.merge_overlapping_exon_queries(exons)}
        self.assertEqual(groups,{frozenset({'E0','E1'}),frozenset({'E2'}),
                                 frozenset({'E3'}),frozenset({'E4'}),frozenset({'E5'})})

    def test_gene_name_sets_merge_independent_of_metadata_order(self):
        exons = [extracted('E1',100,200,self.genome,gene='G1'),
                 extracted('E2',180,260,self.genome,gene='G2')]
        exons[0].record.gene_names = ['B','A','B']
        exons[1].record.gene_names = ['A','B']
        query, = b.merge_overlapping_exon_queries(exons)
        self.assertEqual(len(query.aliases),2)


def alias(identifier, core_start, core_end, anchor_start, anchor_end):
    return a.ExonAlias(identifier+'.1',identifier,'chr1',core_start,core_end,'+',
                       core_end-core_start,'same_gene_overlap',core_start-anchor_start,
                       anchor_end-core_end,anchor_end-anchor_start,
                       core_start,core_end,anchor_start,anchor_end)


def hsp(qseq,sseq,qstart,qend,sstart,send,qlen,score=100):
    return '\t'.join(map(str,['UNION','target',100,len(qseq),0,0,qstart,qend,
                              sstart,send,'1e-50',100,score,qlen,10000,qseq,sseq]))


class AliasProjectionTests(unittest.TestCase):
    def test_projected_score_handles_ambiguous_nucleotides(self):
        projector=a.HspProjector('AANRYYC','AAARRAC',1,7,1,7)
        self.assertEqual(projector.project(0,7)[10],-3)

    def test_gap_projection_updates_each_exons_target_coordinates_and_score(self):
        qseq,sseq='AACCGG-TTAACC','AACCGGATTAACT'
        aliases=[alias('SHORT',4,10,2,12)]
        for sstart,send,expected,strand in [(101,113,(104,111),'+'),(113,101,(102,109),'-')]:
            row,=a.tabular_line_to_alias_alignments(hsp(qseq,sseq,1,12,sstart,send,12),{}, {}, {'UNION':aliases})
            self.assertEqual((row.query_start,row.query_end),expected)
            self.assertEqual((row.exon_start,row.exon_end,row.exon_length),(0,6,6))
            self.assertEqual((row.cigar,row.identical_bases,row.NM,row.AS),('2=1I4=',6,1,3))
            self.assertEqual(row.strand,strand)
            self.assertEqual(row.selection_AS,4)
        row,=a.tabular_line_to_alias_alignments(
            hsp(b.revcomp(qseq),b.revcomp(sseq),12,1,201,213,12),{}, {}, {'UNION':aliases})
        self.assertEqual((row.query_start,row.query_end,row.strand),(202,209,'-'))
        self.assertEqual(row.cigar,'4=1I2=')

    def test_short_complete_exon_survives_a_partial_union_hit(self):
        aliases=[alias('LONG',60,240,40,260),alias('SHORT',100,140,80,160)]
        output=list(a.tabular_line_to_alias_alignments(
            hsp('ACGT'*20,'ACGT'*20,81,160,1001,1080,300,80),{}, {}, {'UNION':aliases}))
        by_id={row.exon_id:row for row in output}
        self.assertLess(by_id['LONG'].selection_coverage,90)
        self.assertEqual(by_id['SHORT'].selection_coverage,100)
        self.assertEqual((by_id['SHORT'].query_start,by_id['SHORT'].query_end),(1020,1060))
        self.assertEqual(by_id['SHORT'].selection_AS,80)

    def test_mismatches_outside_alias_do_not_lower_its_identity(self):
        qseq='A'*300
        sseq='C'*80+'A'*120+'C'*100
        row,=a.tabular_line_to_alias_alignments(
            hsp(qseq,sseq,1,300,1,300,300),{}, {}, {'UNION':[alias('SHORT',100,180,80,200)]})
        self.assertEqual((row.selection_percent_identity,row.selection_coverage,row.selection_AS),(100,100,120))


@unittest.skipUnless(shutil.which('blastn') and shutil.which('makeblastdb'), 'BLAST+ required')
class LiveOverlapMergeTests(unittest.TestCase):
    def test_union_search_matches_original_exons_with_indel_reverse_copy_and_fragment(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            rng=random.Random(809)
            genome=''.join(rng.choices('ACGT',k=900))
            fasta=folder/'ref.fa'; fasta.write_text('>chr1\n'+genome+'\n')
            specs=[('LONG',120,420),('SHIFT',126,420),('SHORT',245,360),('OVERLAP',400,540)]
            lines=['##gff-version 3']
            for name,start,end in specs:
                attrs=f'gene_id=G.1;gene_name=G;transcript_id=T{name}.1;transcript_type=protein_coding'
                lines.append(f'chr1\ttest\ttranscript\t{start+1}\t{end}\t.\t+\t.\tID=T{name}.1;{attrs}')
                lines.append(f'chr1\ttest\texon\t{start+1}\t{end}\t.\t+\t.\t{attrs};exon_id=E{name}.1;exon_number=1')
            gff=folder/'ref.gff3'; gff.write_text('\n'.join(lines)+'\n')
            records=b.collapse_exons(b.parse_gencode_gff3(str(gff)))
            mutated=genome[:310]+'AAA'+genome[310:]
            assembly=folder/'assembly.fa'
            assembly.write_text('>copy\n'+mutated+'\n>reverse\n'+b.revcomp(mutated)+
                                '\n>fragment\n'+genome[185:420]+'\n')
            results=[]
            for merged in (False,True):
                prefix=folder/f'db{merged}'
                b.write_exon_fasta(str(fasta),records,str(prefix)+'.exons.fa',str(prefix)+'.seq',
                    str(prefix)+'.exon_info.tsv',str(prefix)+'.exon_aliases.tsv',150,0,merged)
                count=sum(line.startswith('>') for line in Path(str(prefix)+'.exons.fa').read_text().splitlines())
                self.assertEqual(count,1 if merged else 4)
                aln=folder/f'align{merged}.tsv'
                command=[sys.executable,str(ROOT/'align_exon_blastdb_v2.py'),'-q',str(assembly),
                         '-d',str(prefix),'-o',str(aln),'--exons-as-query','--threads','2',
                         '--word-size','19','--evalue','1e-30','--candidate-aligner','blast','--no-local-realignment']
                result=subprocess.run(command,text=True,capture_output=True)
                self.assertEqual(result.returncode,0,result.stderr)
                self.assertIn('-perc_identity 95.0',result.stderr)
                alignments=rows(aln)
                self.assertTrue(any(r['query_id']=='fragment' and r['exon_id']=='ESHORT' for r in alignments))
                self.assertTrue(any('I' in r['cigar'] for r in alignments))
                self.assertEqual({r['strand'] for r in alignments},{'+','-'})
                calls=folder/f'calls{merged}.transcript_calls.tsv'
                result=subprocess.run([sys.executable,str(ROOT/'call_genes_from_exon_alignments_v3.py'),
                    '-i',str(aln),'-g',str(gff),'-o',str(calls),'--query-coordinate-mode','local',
                    '--eligible-exon-info',str(prefix)+'.exon_aliases.tsv'],text=True,capture_output=True)
                self.assertEqual(result.returncode,0,result.stderr)
                fragments=next(folder.glob(f'calls{merged}.*pseudo*.tsv'))
                normalized=sorted(json.dumps({k:v for k,v in row.items() if k!='evalue'},sort_keys=True)
                                  for row in alignments)
                results.append((normalized,calls.read_bytes(),fragments.read_bytes()))
            self.assertEqual(results[0],results[1])


if __name__=='__main__':
    unittest.main()
