"""Reference-overlap gene units, builder propagation and two-stage selection."""
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import sys
import tempfile
import unittest

from test_exon_transcript_v3 import ROOT, caller as c, runner, cli, rows
import shared_exon_genes as shared


def exon(gene, start, end, chrom='chr1', strand='+', tid=None, name=None):
    tid = tid or 'T' + gene
    return c.GffExon(f'{tid}E{start}', f'{tid}E{start}.1', tid, tid + '.1', 1,
                     gene, gene + '.1', name or gene, 'protein_coding',
                     chrom, start, end, strand, False)


class SharedExonTests(unittest.TestCase):
    def test_strict_total_threshold_and_no_isoform_double_counting(self):
        for size, merged in ((50, False), (51, True), (60, True)):
            with self.subTest(shared_bases=2*size):
                evidence = [exon(g, start, start+size) for g in ('A','B') for start in (0,200)]
                evidence += [replace(e, transcript_id=e.transcript_id+'alt') for e in evidence]
                result, groups = shared.merge_shared_exon_rows(evidence)
                self.assertEqual(bool(groups), merged)
                self.assertEqual(len(result), len(evidence))
                if merged:
                    self.assertEqual(groups[0]['max_pair_shared_exonic_bp'], 2*size)
                    self.assertEqual({e.gene_id for e in result}, {'A&B'})
        # Overlapping alternative exon boundaries still cover only 100 bases.
        evidence = [exon(g,0,70) for g in ('A','B')] + [exon(g,30,100) for g in ('A','B')]
        self.assertFalse(shared.merge_shared_exon_rows(evidence)[1])
        self.assertTrue(shared.merge_shared_exon_rows([exon('A',0,200),exon('B',99,400)])[1])

    def test_components_are_transitive_deterministic_and_locus_specific(self):
        evidence = [exon('A',0,150,name='Zulu'), exon('B',0,150,name='Alpha'),
                    exon('B',300,450,name='Alpha'), exon('C',300,450,name='Middle'),
                    exon('A',0,150,chrom='chr2',name='Zulu'),
                    exon('D',0,150,strand='-'), exon('E',600,900)]
        result, groups = shared.merge_shared_exon_rows(evidence)
        self.assertEqual(groups, shared.merge_shared_exon_rows(list(reversed(evidence)))[1])
        self.assertEqual(len(groups),1)
        self.assertEqual(groups[0]['merged_gene_name'],'Alpha&Middle&Zulu')
        self.assertEqual(groups[0]['merged_gene_id'],'B&C&A')
        self.assertEqual({e.gene_id for e in result if e.chrom=='chr2'}, {'A'})
        self.assertEqual({e.gene_id for e in result if e.strand=='-'}, {'D'})
        self.assertEqual(result[-1].gene_id,'E')
        self.assertFalse(shared.merge_shared_exon_rows(result)[1])

    def test_intron_overlap_and_sequence_similarity_are_not_shared_reference_bases(self):
        evidence = [exon('A',0,150),exon('A',600,750),exon('B',300,450),
                    exon('C',0,150,chrom='chr2'),exon('D',0,150,strand='-')]
        self.assertFalse(shared.merge_shared_exon_rows(evidence)[1])

    def test_real_muc20_annotations_form_one_48_block_unit(self):
        ann=c.parse_gencode_gff3(str(ROOT.parent/'test/regression/data/MUC20_shared_exons.gff3'))
        original=set(ann.transcripts)
        models=c.add_full_gene_transcripts(ann)
        self.assertEqual(len(ann.shared_exon_gene_groups),1)
        group=ann.shared_exon_gene_groups[0]
        self.assertEqual(group['merged_gene_name'],'MUC20&MUC20-OT1')
        self.assertEqual(group['merged_gene_id'],'ENSG00000176945&ENSG00000242086')
        self.assertEqual(group['max_pair_shared_exonic_bp'],1851)
        self.assertEqual(len(models),48)
        self.assertEqual(len(original),258)
        self.assertTrue(original <= set(ann.transcripts))
        self.assertEqual({ann.transcripts[t].gene_id for t in original},{group['merged_gene_id']})
        self.assertTrue(ann.transcripts['ENST00000447234'].is_mane)
        self.assertEqual(c.apply_shared_exon_genes(ann),[group])


class SharedGenePipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name)

    def fixture(self, mane=True, equal_isoforms=False):
        # Shared sequence is 60 + 60 = 120 bp, below the threshold per exon.
        specs=[('TA','GA','A',[(100,180),(200,260),(400,460)],mane,'protein_coding'),
               ('TB','GB','B',[(200,260),(400,460),(800,920)],False,'lncRNA')]
        if equal_isoforms:
            specs=[('TA','GA','A',[(200,260),(400,460)],False,'protein_coding'),
                   ('TB','GB','B',[(190,270),(390,470)],False,'protein_coding')]
        gff=self.folder/'ref.gff3'
        text=['##gff-version 3']
        alignments=[]
        for tid,gene,name,spans,is_mane,kind in specs:
            attrs=f'gene_id={gene}.1;gene_name={name};transcript_id={tid}.1;transcript_type={kind}'
            if is_mane:
                attrs+=';tag=MANE_Select'
            text.append(f'chr1\ttest\ttranscript\t{spans[0][0]+1}\t{spans[-1][1]}\t.\t+\t.\tID={tid}.1;{attrs}')
            for number,(start,end) in enumerate(spans,1):
                eid=f'{tid}E{number}'
                text.append(f'chr1\ttest\texon\t{start+1}\t{end}\t.\t+\t.\t{attrs};exon_number={number};exon_id={eid}.1')
                for query in ('copy1','copy2'):
                    alignments.append(c.ExonAlignment(query,100,start+10000,end+10000,'+',eid,
                        end-start,0,end-start,end-start,100,end-start,end-start))
        gff.write_text('\n'.join(text)+'\n')
        aln=self.folder/'alignments.tsv'
        with aln.open('w') as h:
            writer=csv.DictWriter(h,list(c.ExonAlignment.__dataclass_fields__),delimiter='\t')
            writer.writeheader()
            writer.writerows(map(asdict,alignments))
        return gff,aln

    def test_merged_parent_keeps_all_member_isoforms_and_mane_recovers_full_transcript(self):
        gff,aln=self.fixture()
        outputs=[]
        for threads in (1,2):
            output=self.folder/f'calls{threads}.tsv'
            report=self.folder/f'groups{threads}.tsv'
            cli('-i',aln,'-g',gff,'-o',output,'--threads',threads,'--shard-storage','disk',
                '--shared-exon-genes-output',report,'--query-coordinate-mode','local')
            outputs.append(output)
            self.assertEqual(shared.read_report(report)[0]['member_gene_ids'],'GA;GB')
        self.assertEqual(outputs[0].read_bytes(),outputs[1].read_bytes())
        data=rows(outputs[0])
        parents=[r for r in data if r['model_type']=='full_gene']
        children=[r for r in data if r['model_type']=='transcript']
        self.assertEqual({r['gene_id'] for r in parents},{'GA&GB'})
        self.assertEqual({r['gene_name'] for r in parents},{'A&B'})
        self.assertEqual({(r['transcript_id'],r['gene_id'],r['gene_name']) for r in children},
                         {('TA','GA','A'),('TB','GB','B')})
        self.assertEqual({r['tie_count'] for r in data},{'1'})
        for query in ('copy1','copy2'):
            group=[r for r in data if r['query_id']==query]
            self.assertEqual(len({r['GENE_index'] for r in group}),1)
            parent=next(r for r in group if r['model_type']=='full_gene')
            self.assertEqual(parent['transcript_id'],'GA&GB')
            self.assertEqual((parent['found_exons'],parent['expected_exons']),('4','4'))
            mane=next(r for r in group if r['transcript_id']=='TA')
            self.assertEqual((mane['call_status'],mane['found_exons']),('complete','3'))
            self.assertTrue(any(r['transcript_id']=='TB' for r in group))

    def test_equal_member_isoforms_are_one_gene_tiebreak_not_cross_gene_ties(self):
        gff,aln=self.fixture(equal_isoforms=True)
        for extra in ((),('--no-full-gene-transcripts',)):
            out=self.folder/f'tie{len(extra)}.tsv'
            cli('-i',aln,'-g',gff,'-o',out,'--query-coordinate-mode','local',*extra)
            transcripts=[r for r in rows(out) if r['model_type']=='transcript']
            self.assertEqual({r['transcript_id'] for r in transcripts},{'TB'})
            self.assertEqual({r['gene_id'] for r in transcripts},{'GB'})
            self.assertEqual({r['gene_name'] for r in transcripts},{'B'})
            self.assertEqual({r['tie_count'] for r in transcripts},{'1'})

    def test_builder_keeps_both_genes_and_cache_tracks_shared_gene_report(self):
        gff,_aln=self.fixture()
        fasta=self.folder/'ref.fa'
        fasta.write_text('>chr1\n'+'ACGT'*300+'\n')
        prefix=self.folder/'db/exons'
        runner.ensure_database(prefix,runner.load_scripts(ROOT),sys.executable,
                               fasta,gff,3,0,99.,False)
        report=prefix.parent/shared.REPORT_NAME
        self.assertTrue(runner.validate_database(prefix)[0])
        self.assertEqual(shared.read_report(report)[0]['merged_gene_name'],'A&B')
        aliases=rows(str(prefix)+'.exon_aliases.tsv')
        self.assertEqual({r['gene_id'] for r in aliases},{'GA&GB'})
        self.assertEqual({tid for r in aliases for tid in r['transcript_id'].split(',')},{'TA','TB'})
        manifest_path=Path(str(prefix)+'.manifest.json')
        manifest=json.loads(manifest_path.read_text())
        manifest['shared_exon_gene_policy']='old-policy'
        manifest_path.write_text(json.dumps(manifest))
        self.assertFalse(runner.validate_database(prefix)[0])
        manifest['shared_exon_gene_policy']=shared.MERGE_POLICY
        manifest_path.write_text(json.dumps(manifest))
        report.write_text(report.read_text().replace('A&B','Changed'))
        self.assertFalse(runner.validate_database(prefix)[0])

    def test_identical_mane_shortcut_cannot_discard_shared_gene_members(self):
        gff,aln=self.fixture(equal_isoforms=True)
        # The two isoforms now have identical coordinates, DNA and MANE tags.
        text=gff.read_text().replace('191\t270','201\t260').replace('391\t470','401\t460')
        text=text.replace('191\t470','201\t460')
        gff.write_text('\n'.join(line+';tag=MANE_Select' if not line.startswith('#') else line
                                 for line in text.splitlines())+'\n')
        fasta=self.folder/'ref.fa'
        fasta.write_text('>chr1\n'+'ACGT'*300+'\n')
        prefix=self.folder/'db/exons'
        runner.ensure_database(prefix,runner.load_scripts(ROOT),sys.executable,
                               fasta,gff,3,0,99.,False)
        self.assertEqual(rows(prefix.parent/'indenticalparalogs.tsv'),[])
        self.assertEqual(shared.read_report(prefix.parent/shared.REPORT_NAME)[0]['merged_gene_id'],'GA&GB')
        aliases=rows(str(prefix)+'.exon_aliases.tsv')
        self.assertEqual({tid for r in aliases for tid in r['transcript_id'].split(',')},{'TA','TB'})
        # A saved report from an older builder must not delete either member.
        from identical_paralogs import write_report
        legacy=self.folder/'legacy.tsv'
        write_report(legacy,[dict(representative_gene_id='GA',representative_gene_name='A',
            merged_gene_name='Amerged',gene_ids='GA;GB',gene_names='A;B',
            mane_transcript_ids='TA;TB',mane_sequence_lengths='120;120',
            mane_sequence_sha256='same;same')])
        ann=c.parse_gencode_gff3(str(gff))
        c.apply_identical_paralogs(ann,legacy)
        self.assertEqual(set(ann.transcripts),{'TA','TB'})
        self.assertEqual({info.gene_name for info in ann.transcripts.values()},{'A&B'})

    def test_saved_real_muc20_alignments_recover_mane_under_merged_parent(self):
        data=ROOT.parent/'test/regression/data'
        out=self.folder/'muc.tsv'
        cli('-i',data/'MUC20_CHM13_target.alignments.tsv',
            '-g',data/'MUC20_shared_exons.gff3','-o',out,
            '--eligible-exon-info',data/'MUC20_CHM13_target.exon_ids.tsv')
        calls=rows(out)
        parents=[r for r in calls if r['model_type']=='full_gene']
        self.assertEqual(len(parents),1)
        parent=parents[0]
        self.assertEqual(parent['gene_name'],'MUC20&MUC20-OT1')
        self.assertEqual((parent['found_exons'],parent['expected_exons']),('46','48'))
        self.assertEqual((parent['query_start'],parent['query_end']),('198625730','198709445'))
        mane=next(r for r in calls if r['transcript_id']=='ENST00000447234')
        self.assertEqual((mane['found_exons'],mane['expected_exons']),('3','4'))
        self.assertEqual(mane['GENE_index'],parent['GENE_index'])
        self.assertEqual(parent['gene_id'],'ENSG00000176945&ENSG00000242086')
        self.assertEqual((mane['gene_id'],mane['gene_name']),('ENSG00000176945','MUC20'))
        self.assertEqual({r['gene_name'] for r in calls if r['model_type']=='transcript'},
                         {'MUC20','MUC20-OT1'})


if __name__=='__main__':
    unittest.main()
