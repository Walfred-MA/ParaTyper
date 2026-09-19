"""Synthetic exon unions, full-gene/MANE priority, and real-locus regressions."""
from dataclasses import asdict, replace
from pathlib import Path
import tempfile
import unittest

from test_exon_transcript_v3 import ROOT, caller as c, cli, make_call, rows, bed


class FullGeneTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def annotation(self, strand='+', extra_locus=False):
        lines = ['##gff-version 3']
        specs = [('TA','GA','chr1',[(100,200,'E1'),(400,450,'E3')],True),
                 ('TB','GA','chr1',[(150,250,'E2'),(400,450,'E3'),(500,530,'E4')],False),
                 ('TC','GC','chr1',[(100,200,'E1')],False)]
        if extra_locus:
            specs.append(('TD','GA','chr2',[(100,200,'E5')],False))
        for tid,gene,chrom,exons,mane in specs:
            attrs=f'gene_id={gene}.1;gene_name={gene};transcript_id={tid}.1;transcript_type=protein_coding'
            if mane:
                attrs+=';tag=MANE_Select'
            lines.append(f'{chrom}\ttest\ttranscript\t101\t530\t.\t{strand}\t.\tID={tid}.1;{attrs}')
            for num,(start,end,eid) in enumerate(sorted(exons,reverse=strand=='-'),1):
                lines.append(f'{chrom}\ttest\texon\t{start+1}\t{end}\t.\t{strand}\t.\t'
                             f'{attrs};exon_id={eid}.1;exon_number={num}')
        path=self.folder/'input.gff3'
        path.write_text('\n'.join(lines)+'\n')
        return c.parse_gencode_gff3(str(path))

    def test_union_contains_all_isoforms_without_counting_overlap_twice(self):
        ann=self.annotation()
        before={tid:asdict(info) for tid,info in ann.transcripts.items()}
        report=c.add_full_gene_transcripts(ann)
        model=[r for r in report if r['transcript_id']=='GA_full_gene']
        self.assertEqual([(r['start'],r['end'],r['exon_number']) for r in model],
                         [(100,250,1),(400,450,2),(500,530,3)])
        self.assertEqual(model[0]['source_exon_ids'],'E1,E2')
        self.assertEqual(model[0]['source_transcript_ids'],'TA,TB')
        self.assertEqual(ann.transcripts['GA_full_gene'].exon_length,230)
        self.assertFalse(ann.transcripts['GA_full_gene'].is_mane)
        self.assertEqual(ann.transcripts['GA_full_gene'].model_type,'full_gene')
        self.assertEqual(before,{tid:asdict(ann.transcripts[tid]) for tid in before})
        self.assertEqual({e.exon_number for e in ann.exons_by_transcript['GC_full_gene']},{1})
        self.assertEqual({item[0] for item in ann.exon_to_transcripts['E1']},
                         {'TA','TC','GA_full_gene','GC_full_gene'})

    def test_reverse_strand_and_distinct_reference_loci(self):
        ann=self.annotation('-',extra_locus=True)
        report=c.add_full_gene_transcripts(ann)
        model=[r for r in report if r['transcript_id']=='GA_full_gene_1']
        self.assertEqual([(r['start'],r['exon_number']) for r in model],[(500,1),(400,2),(100,3)])
        self.assertEqual({r['reference_contig'] for r in model},{'chr1'})
        other=[r for r in report if r['transcript_id']=='GA_full_gene_2']
        self.assertEqual(len(other),1)
        self.assertEqual(other[0]['reference_contig'],'chr2')

    def test_eligibility_and_truncation_use_original_exon_aliases(self):
        ann=self.annotation()
        c.add_full_gene_transcripts(ann)
        self.assertEqual(c.truncate_expected_exons(ann,None,{'E2'})['GA_full_gene'],[1])
        self.assertEqual(c.truncate_expected_exons(ann,c.Region('chr1',425,535))['GA_full_gene'],[2,3])

    def test_gene_stage_uses_score_and_transcript_stage_preserves_mane(self):
        full=make_call('GA_full_gene',identical=98,mane=False)
        full=replace(full,transcript_info=replace(full.transcript_info,model_type='full_gene'))
        other=replace(full,transcript_id='GB_full_gene',hits=make_call('other',identical=100).hits)
        self.assertEqual(c.transcript_call_rank_key(full,False)[0],0)
        self.assertEqual([v.transcript_id for v in c.resolve_call_overlaps([full,other],False)],['GB_full_gene'])
        mane=make_call('MANE',identical=95)
        ordinary=make_call('ordinary',mane=False,identical=100)
        self.assertEqual([v.transcript_id for v in c.resolve_call_overlaps([mane,ordinary],True)],['MANE'])
        self.assertEqual([v.transcript_id for v in c.resolve_call_overlaps([mane,ordinary],False)],['ordinary'])

    def test_equal_score_prefers_longer_model_within_gene_and_bed_labels_it(self):
        mane=make_call('MANE',transcript_length=200)
        full=replace(mane,transcript_id='GA_full_gene',transcript_info=replace(
            mane.transcript_info,transcript_id='GA_full_gene',transcript_id_full='GA_full_gene',
            is_mane=False,model_type='full_gene',exon_length=230))
        self.assertEqual(c.resolve_call_overlaps([mane,full],False)[0].transcript_id,'GA_full_gene')
        shorter_full=replace(full,transcript_id='A_short_full',transcript_info=replace(full.transcript_info,exon_length=200))
        self.assertEqual(c.resolve_call_overlaps([shorter_full,full])[0].transcript_id,'GA_full_gene')
        row=dict(zip(c.OUTPUT_HEADER,c.call_to_output_row(full)))
        self.assertEqual(row['ifmane_transcript'],'0')
        self.assertEqual(row['model_type'],'full_gene')
        record=bed.convert(row,'sample')
        self.assertEqual(dict(zip(bed.EXTRA_COLUMNS,record.fields[12:]))['model_type'],'full_gene')

    def run_fixture(self,name,threads=1,extra=()):
        data=ROOT.parent/'test/regression/data'
        out=self.folder/f'{name}_{threads}_{len(extra)}.tsv'
        cli('-i',data/f'{name}.alignments.tsv','-g',data/f'{name}.gff3',
            '--eligible-exon-info',data/f'{name}.exon_ids.tsv','-o',out,
            '--threads',threads,'--shard-storage','disk','--query-coordinate-mode','local',*extra)
        return out

    def test_c4orf50_is_one_full_gene_and_opt_out_preserves_old_transcript_calls(self):
        report=self.folder/'models.tsv'
        serial=self.run_fixture('C4orf50_CHM13',1,('--full-gene-models-output',report))
        parallel=self.run_fixture('C4orf50_CHM13',2)
        self.assertEqual(serial.read_bytes(),parallel.read_bytes())
        calls=rows(serial)
        self.assertEqual(len(calls),4)
        self.assertEqual(len({r['GENE_index'] for r in calls}),1)
        call=calls[0]
        self.assertEqual(call['transcript_id'],'ENSG00000181215')
        self.assertEqual(call['gene_id'],call['transcript_id'])
        self.assertEqual((call['found_exons'],call['expected_exons'],call['call_status']),('35','35','complete'))
        self.assertEqual((call['query_start'],call['query_end']),('5869135','6174552'))
        self.assertEqual(call['ifmane_transcript'],'0')
        self.assertEqual(len(rows(report)),35)
        old=rows(self.run_fixture('C4orf50_CHM13',1,('--no-full-gene-transcripts',)))
        self.assertEqual({r['transcript_id'] for r in old},
                         {'ENST00000711657','ENST00000531445','ENST00000639345'})
        self.assertTrue(all(r['model_type']=='transcript' for r in old))

    def test_c4_full_gene_models_stay_in_separate_copies(self):
        all_calls=rows(self.run_fixture('C4_CHM13_tandem',2))
        calls=[r for r in all_calls if r['model_type']=='full_gene']
        self.assertEqual(len(calls),2)
        self.assertEqual(len(all_calls),4)
        for child in [r for r in all_calls if r['model_type']=='transcript']:
            parent=next(r for r in calls if r['GENE_index']==child['GENE_index'])
            self.assertEqual(parent['gene_id'],child['gene_id'])
            self.assertLessEqual(int(parent['query_start']),int(child['query_start']))
            self.assertGreaterEqual(int(parent['query_end']),int(child['query_end']))
        by_gene={r['gene_name']:r for r in calls}
        self.assertEqual(set(by_gene),{'C4A','C4B'})
        self.assertTrue(all(r['model_type']=='full_gene' and r['call_status']=='complete' for r in calls))
        self.assertLess(int(by_gene['C4A']['query_end']),int(by_gene['C4B']['query_start']))

    def test_smn_main_gene_assignments_are_preserved_with_full_gene_candidates(self):
        calls=rows(self.run_fixture('SMN_CHM13_gene_filter',2))
        full=[r for r in calls if r['model_type']=='full_gene' and r['call_status']=='complete']
        self.assertEqual({(r['gene_name'],r['strand']) for r in full},{('SMN1','+'),('SMN2','-')})
        self.assertTrue(all(r['assignment_status']=='unique' for r in full))


if __name__=='__main__':
    unittest.main()
