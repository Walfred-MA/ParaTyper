"""Candidate gene-like fragment partition, source names and two-table runners."""
import csv
from dataclasses import replace
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT, caller as c, runner, cli, make_call, rows, write_fixture


def parent(name, found=1, expected=3, biotype='protein_coding'):
    base=make_call(name+'_full_gene',gene=name,mane=False,
                   intervals=tuple((100+i*200,200+i*200) for i in range(found)),
                   expected=list(range(1,expected+1)))
    return replace(base,transcript_info=replace(base.transcript_info,model_type='full_gene',
                   transcript_type=biotype),protein_bonus=1.0,gene_index=name,
                   gene_sort_score=1000,gene_start=100,gene_end=found*200)


def child(p, mane=False, complete=False):
    tid=p.transcript_info.gene_id+'_T'
    return replace(p,transcript_id=tid,transcript_info=replace(p.transcript_info,
        transcript_id=tid,transcript_id_full=tid,model_type='transcript',is_mane=mane),
        expected_exon_numbers=list(range(1,len(p.hits)+(1 if complete else 2))))


class PartitionTests(unittest.TestCase):
    def test_all_three_criteria_and_reference_biotype_independence(self):
        for biotype in ('protein_coding','lncRNA','processed_pseudogene'):
            p=parent('newcopy',biotype=biotype)
            self.assertTrue(c.is_pseudofragment_assignment([[p],[child(p)]]))
            self.assertTrue(c.is_pseudofragment_assignment([[p],[child(p,mane=True)]]))
            self.assertFalse(c.is_pseudofragment_assignment([[p],[child(p,mane=True,complete=True)]]))
            self.assertTrue(c.is_pseudofragment_assignment([[p],[child(p,complete=True)]]))
        # Known complete single-exon pseudogenes remain ordinary gene calls.
        p=parent('known',found=1,expected=1,biotype='processed_pseudogene')
        self.assertFalse(c.is_pseudofragment_assignment([[p],[child(p,complete=True)]]))
        self.assertFalse(c.is_pseudofragment_assignment([[parent('multi',found=2)]]))
        self.assertFalse(c.is_pseudofragment_assignment([[child(parent('nogene'))]]))

    def test_tied_genes_move_only_when_all_alternatives_qualify(self):
        a,b=parent('A'),parent('B')
        self.assertTrue(c.is_pseudofragment_assignment([[a,b],[child(a),child(b)]]))
        self.assertFalse(c.is_pseudofragment_assignment([[a,parent('B',expected=1)]]))
        self.assertFalse(c.is_pseudofragment_assignment([[a,parent('B',found=2)]]))
        self.assertFalse(c.is_pseudofragment_assignment([[a,b],[child(b,mane=True,complete=True)]]))

    def test_routes_whole_parent_with_global_indices_and_preserves_ties(self):
        a=replace(parent('A'),gene_index='AB',tie_group_id='gene_tie',tie_count=2)
        b=replace(parent('B'),gene_index='AB',tie_group_id='gene_tie',tie_count=2)
        a_child=replace(child(a),tie_group_id='',tie_count=1)
        known=parent('known',expected=1,biotype='processed_pseudogene')
        protected=parent('mane')
        calls=[a,b,a_child,known,child(known,complete=True),protected,child(protected,mane=True,complete=True)]
        calls.sort(key=c.call_output_sort_key,reverse=True)
        with tempfile.TemporaryDirectory() as tmp:
            main=Path(tmp)/'calls.tsv'
            processed=Path(tmp)/'pseudofragments.tsv'
            written=c.write_transcript_calls(iter(calls),str(main),'extended',False,False,False,
                                             pseudofragments_output=str(processed))
            ordinary,pseudo=rows(main),rows(processed)
            self.assertEqual(written,len(ordinary)+len(pseudo))
            self.assertEqual(len(pseudo),2)
            self.assertEqual(pseudo[0]['gene_id'],'A;B')
            self.assertEqual(pseudo[0]['tie_count'],'2')
            self.assertEqual(pseudo[0]['found_exon_numbers'],'1;1')
            self.assertEqual(pseudo[0]['GENE_index'],pseudo[1]['GENE_index'])
            self.assertFalse({r['GENE_index'] for r in ordinary}&{r['GENE_index'] for r in pseudo})
            self.assertEqual({r['gene_id'] for r in ordinary},{'known','mane'})
            self.assertEqual(main.read_text().splitlines()[0],processed.read_text().splitlines()[0])

    def test_transcript_only_stays_main_and_empty_processed_has_header(self):
        call=replace(child(parent('A')),gene_index='')
        with tempfile.TemporaryDirectory() as tmp:
            main,processed=Path(tmp)/'calls.tsv',Path(tmp)/'pseudofragments.tsv'
            c.write_transcript_calls([call],str(main),'extended',False,False,False,
                                    pseudofragments_output=str(processed))
            self.assertEqual(len(rows(main)),1)
            self.assertEqual(rows(processed),[])
            self.assertEqual(processed.read_text().strip(),'\t'.join(c.OUTPUT_HEADER))

    def test_duplicate_output_path_is_rejected_before_overwriting(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'calls.tsv'
            out.write_text('existing results')
            with self.assertRaisesRegex(ValueError,'paths must differ'):
                c.write_transcript_calls([],str(out),'extended',False,False,False,
                                        pseudofragments_output=str(out))
            self.assertEqual(out.read_text(),'existing results')


class TwoTablePipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)
        self.gff,self.alignments=write_fixture(self.folder)
        data=rows(self.alignments)
        data=[r for r in data if not (r['query_id']=='sampleA' and r['exon_id'] in {'E2','E4'})]
        with self.alignments.open('w') as handle:
            writer=csv.DictWriter(handle,list(data[0]),delimiter='\t',lineterminator='\n')
            writer.writeheader()
            writer.writerows(data)

    def test_serial_and_parallel_tables_match_for_complete_and_partial_copies(self):
        pairs=[]
        for threads in (1,2):
            out=self.folder/f'calls{threads}.transcript_calls.tsv'
            cli('-i',self.alignments,'-g',self.gff,'-o',out,'--threads',threads,
                '--shard-storage','disk','--query-coordinate-mode','local')
            pairs.append((out,Path(c.default_pseudofragments_path(str(out)))))
        self.assertEqual(pairs[0][0].read_bytes(),pairs[1][0].read_bytes())
        self.assertEqual(pairs[0][1].read_bytes(),pairs[1][1].read_bytes())
        self.assertEqual({r['query_id'] for r in rows(pairs[0][0])},{'sampleB'})
        self.assertEqual({r['query_id'] for r in rows(pairs[0][1])},{'sampleA'})
        self.assertEqual({r['gene_id'] for r in rows(pairs[0][1])},{'GA;GC'})

    def runner_args(self):
        query=self.folder/'query.fa'
        query.write_text('>sampleA\nACGT\n>sampleB\nACGT\n>empty\nACGT\n')
        with patch.object(sys,'argv',['annotate_assemblies.py','-r',str(query),'-g',str(self.gff),
            '--query-fasta',str(query),'-d',str(self.folder),'-o',str(self.folder/'out'),
            '--caller-threads','2']):
            return runner.parse_args()

    def fake_alignment_real_caller(self,command,label):
        if Path(command[1]).name=='align_exon_blastdb_v2.py':
            Path(command[command.index('--output')+1]).write_bytes(self.alignments.read_bytes())
        else:
            result=subprocess.run(command,text=True,capture_output=True)
            if result.returncode:
                raise AssertionError(result.stderr)

    def test_named_runner_requires_both_outputs_to_resume(self):
        args=self.runner_args()
        def run():
            return runner.run_one_sample(runner.AssemblyQuery('assembly',args.query_fasta),args,
                runner.load_scripts(ROOT),self.folder/'reference_exons',args.output,args.output/'temp')
        with patch.object(runner,'run_command',side_effect=self.fake_alignment_real_caller) as command:
            main=run()
            processed=main.with_name('assembly.pseudofragments.tsv')
            self.assertEqual({p.name for p in main.parent.iterdir()},
                             {'assembly.transcript_calls.tsv', 'assembly.pseudofragments.tsv'})
            self.assertEqual({r['query_id'] for r in rows(main)},{'sampleB'})
            self.assertEqual({r['query_id'] for r in rows(processed)},{'sampleA'})
            run()
            self.assertEqual(command.call_count,2)
            processed.unlink()
            run()
            self.assertEqual(command.call_count,4)
            self.assertTrue(runner.valid_call_table(processed))

    def test_combined_runner_splits_both_tables_including_empty_samples(self):
        args=self.runner_args()
        def run():
            return runner.run_combined_fasta(args.query_fasta,['sampleA','sampleB','empty'],args,
                runner.load_scripts(ROOT),self.folder/'reference_exons',args.output,args.output/'temp')
        with patch.object(runner,'run_command',side_effect=self.fake_alignment_real_caller) as command:
            paths=run()
            processed=[p.with_name(p.parent.name+'.pseudofragments.tsv') for p in paths]
            for path in paths:
                sample=path.parent.name
                self.assertEqual({p.name for p in path.parent.iterdir()},
                                 {sample+'.transcript_calls.tsv', sample+'.pseudofragments.tsv'})
            self.assertEqual([len(rows(p)) for p in paths],[0,2,0])
            self.assertEqual([len(rows(p)) for p in processed],[2,0,0])
            self.assertTrue(all(runner.valid_call_table(p) for p in paths+processed))
            run()
            self.assertEqual(command.call_count,2)
            processed[0].unlink()
            run()
            self.assertEqual(command.call_count,4)


if __name__=='__main__':
    unittest.main()
