"""Gene identity must not be biased by transcript protein-coding bonuses."""
from pathlib import Path
import tempfile
import unittest

from test_exon_transcript_v3 import ROOT, caller as c, cli, rows


class AmyPseudogeneTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.data=ROOT.parent/'test/regression/data'

    def run_fixture(self,threads=1):
        output=self.folder/f'AMY_{threads}.tsv'
        cli('-i',self.data/'AMY_CHM13.alignments.tsv','-g',self.data/'AMY_CHM13.gff3',
            '--eligible-exon-info',self.data/'AMY_CHM13.exon_ids.tsv',
            '--identical-paralogs',self.data/'AMY_CHM13.identical_paralogs.tsv',
            '-o',output,'--threads',threads,'--shard-storage','disk','--query-coordinate-mode','local')
        return output

    def test_complete_amyp1_copies_win_and_do_not_form_the_long_amy2a_parent(self):
        serial=self.run_fixture(1);parallel=self.run_fixture(4)
        self.assertEqual(serial.read_bytes(),parallel.read_bytes())
        calls=rows(serial)
        parents=[r for r in calls if r['model_type']=='full_gene']
        pseudogenes=[r for r in parents if r['gene_name']=='AMYP1']
        self.assertEqual(len(pseudogenes),3)
        for row in pseudogenes:
            self.assertEqual((row['found_exons'],row['expected_exons'],row['call_status']),('7','7','complete'))
            self.assertEqual(row['raw_exon_score'],'600.000000')
            self.assertEqual(row['weighted_score'],'1200.000000')
            if row['strand']=='-':
                self.assertEqual(row['mean_identity'],'100.000000')
        transcripts=[r for r in calls if r['gene_name']=='AMYP1' and r['model_type']=='transcript']
        self.assertEqual({(int(r['query_start']),int(r['query_end']),r['strand']) for r in transcripts},
            {(103575109,103581258,'-'),(103656732,103662881,'+'),(103750862,103757011,'+')})
        self.assertEqual({r['GENE_index'] for r in pseudogenes},{r['GENE_index'] for r in transcripts})
        self.assertFalse(any(r['gene_name']=='AMY2A' and r['strand']=='-'
                             and int(r['query_start'])<103581258 and int(r['query_end'])>103575109
                             for r in parents))
        # The genuine coding copies must also remain identifiable.
        self.assertEqual(len([r for r in parents if r['gene_name']=='AMY1Amerged' and r['call_status']=='complete']),7)
        coding=[r for r in parents if r['gene_name']=='AMY2A' and r['call_status']=='complete']
        self.assertEqual(len(coding),1)
        self.assertEqual((coding[0]['query_start'],coding[0]['query_end']),('103466340','103474647'))
        self.assertEqual(coding[0]['weighted_score'],'1800.000000')

    def test_protein_bonus_changes_transcript_scores_but_not_gene_assignments(self):
        ann=c.parse_gencode_gff3(str(self.data/'AMY_CHM13.gff3'))
        c.apply_identical_paralogs(ann,self.data/'AMY_CHM13.identical_paralogs.tsv')
        c.add_full_gene_transcripts(ann)
        expected=c.truncate_expected_exons(ann,None,c.load_eligible_exon_ids(str(self.data/'AMY_CHM13.exon_ids.tsv')))
        alignments=c.dedup_same_exon_overlaps(c.read_alignments(str(self.data/'AMY_CHM13.alignments.tsv')))
        outcomes=[]
        for bonus in (1.,10.,100.):
            batches=[]
            calls,_=c.build_hierarchical_calls(alignments,ann,expected,bonus,2.,1.,10,True,diagnostic_batches=batches)
            genes=[v for v in calls if v.transcript_info.model_type=='full_gene']
            self.assertTrue(all(v.protein_bonus==1 for v in batches[0]['calls']))
            outcomes.append([(v.transcript_id,c.call_query_bounds(v),v.weighted_score,v.tie_group_id) for v in genes])
            coding=next(v for v in calls if v.transcript_info.model_type=='transcript'
                        and v.transcript_info.gene_name=='AMY2A' and v.missing_count==0)
            self.assertEqual(coding.protein_bonus,bonus)
            self.assertEqual(coding.weighted_score,coding.raw_score*bonus*2)
        self.assertEqual(outcomes[0],outcomes[1])
        self.assertEqual(outcomes[1],outcomes[2])


if __name__=='__main__':
    unittest.main()
