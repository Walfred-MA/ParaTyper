"""End-to-end MANE sequence merging, caller propagation and database cache checks."""
import csv
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT, caller, runner, rows
import build_exon_blastdb_v2 as builder
from identical_paralogs import REPORT_NAME, read_report


class ReferenceFixture:
    def __init__(self, folder):
        self.folder = folder
        self.genome = folder / 'ref.fa'
        self.gff = folder / 'annotation.gff3'
        self.lines = ['##gff-version 3']
        self.contigs = []
        self.locations = {}

    def transcript(self, gene, name, tid, exons, mane=True, strand='+', shared_ids=None):
        chrom = 'chr_' + tid
        pieces = ['CCCCC']
        positions = {}
        order = list(enumerate(exons, 1))
        if strand == '-':
            order.reverse()
        for number, sequence in order:
            start = sum(map(len, pieces))
            pieces.append(builder.revcomp(sequence) if strand == '-' else sequence)
            positions[number] = (start, start + len(sequence))
            pieces.append('TGTGTGT')
        sequence = ''.join(pieces)
        self.contigs.append((chrom, sequence))
        attrs = f'gene_id={gene}.1;gene_name={name};transcript_id={tid}.1;transcript_type=protein_coding'
        if mane:
            attrs += ';tag=MANE_Select'
        self.lines.append(f'{chrom}\ttest\ttranscript\t1\t{len(sequence)}\t.\t{strand}\t.\tID={tid}.1;{attrs}')
        for number in sorted(positions):
            start, end = positions[number]
            eid = shared_ids[number-1] if shared_ids else f'{tid}E{number}'
            self.lines.append(f'{chrom}\ttest\texon\t{start+1}\t{end}\t.\t{strand}\t.\t'
                              f'{attrs};exon_number={number};exon_id={eid}.1')
            self.locations[(tid, number)] = (eid, end-start)

    def write(self):
        self.genome.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in self.contigs))
        self.gff.write_text('\n'.join(self.lines) + '\n')

    def merge(self):
        self.write()
        mane, order = {}, []
        exons = builder.parse_gencode_gff3(str(self.gff), mane, order)
        return builder.merge_identical_mane_genes(exons, str(self.genome), mane, order)

    def build(self, prefix):
        self.write()
        subprocess.run([sys.executable, str(ROOT / 'build_exon_blastdb_v2.py'),
                        '--genome', str(self.genome), '--gff3', str(self.gff), '--out', str(prefix),
                        '--anchor-target-length', '14', '--min-unmasked', '0', '--no-makeblastdb'],
                       check=True, capture_output=True, text=True)


class IdenticalParalogTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.fixture = ReferenceFixture(self.folder)
        self.sequence = ['ACGTTGCA', 'GATTACAA']

    def test_first_gene_wins_across_strands_case_and_exon_boundaries(self):
        f = self.fixture
        f.transcript('GZ', 'Zeta', 'TZ', self.sequence)
        f.transcript('GA', 'Alpha', 'TA', [''.join(self.sequence).lower()], strand='-')
        f.transcript('GZ', 'Zeta', 'TZalt', ['CCCCAAAA'], mane=False)
        retained, report = f.merge()
        self.assertEqual(len(report), 1)
        self.assertEqual(report[0]['gene_ids'], 'GZ;GA')
        self.assertEqual(report[0]['merged_gene_name'], 'Zetamerged')
        self.assertEqual(report[0]['mane_sequence_lengths'], '16;16')
        self.assertEqual({r.gene_id for r in retained}, {'GZ'})
        self.assertEqual({r.gene_name for r in retained}, {'Zetamerged'})
        self.assertEqual({r.transcript_id for r in retained}, {'TZ', 'TZalt'})

    def test_utr_difference_no_mane_and_ambiguous_bases_do_not_merge(self):
        f = self.fixture
        f.transcript('G1', 'One', 'T1', self.sequence)
        f.transcript('G2', 'Two', 'T2', [self.sequence[0], 'GATTACAT'])
        f.transcript('G3', 'Three', 'T3', self.sequence, mane=False)
        f.transcript('N1', 'Ambiguous1', 'TN1', ['NNNNAAAA'])
        f.transcript('N2', 'Ambiguous2', 'TN2', ['NNNNAAAA'])
        retained, report = f.merge()
        self.assertFalse(report)
        self.assertEqual(len({r.gene_id for r in retained}), 5)

    def test_all_mane_sequences_must_match_without_transitive_partial_merge(self):
        f = self.fixture
        f.transcript('G1', 'One', 'T1', ['AAAA'])
        f.transcript('G2', 'Two', 'T2a', ['AAAA'])
        f.transcript('G2', 'Two', 'T2b', ['CCCC'])
        f.transcript('G3', 'Three', 'T3', ['CCCC'])
        f.transcript('G4', 'Four', 'T4b', ['CCCC'])
        f.transcript('G4', 'Four', 'T4a', ['AAAA'])
        retained, report = f.merge()
        self.assertEqual([r['gene_ids'] for r in report], ['G2;G4'])
        self.assertEqual({r.gene_id for r in retained}, {'G1', 'G2', 'G3'})

    def test_missing_mane_exon_prevents_merge_and_duplicate_rows_do_not_inflate(self):
        f = self.fixture
        f.transcript('G1', 'One', 'T1', self.sequence)
        f.lines.append(f.lines[-1])
        f.transcript('G2', 'Two', 'T2', self.sequence)
        f.lines.append('chr_T2\ttest\ttranscript\t1\t20\t.\t+\t.\t'
                       'transcript_id=Tmissing;gene_id=G2;gene_name=Two;tag=MANE_Select')
        _, report = f.merge()
        self.assertFalse(report)
        f.lines.pop()
        _, report = f.merge()
        self.assertEqual(report[0]['mane_sequence_lengths'], '16;16')

    def test_builder_caller_report_propagation_and_shared_exon_ids(self):
        f = self.fixture
        f.transcript('GZ', 'Zeta', 'TZ', self.sequence, shared_ids=['E1','E2'])
        f.transcript('GA', 'Alpha', 'TA', self.sequence, strand='-', shared_ids=['E1','E2'])
        f.transcript('GZ', 'Zeta', 'TZalt', ['AACCGGTT'], mane=False)
        prefix = self.folder / 'db' / 'exons'
        f.build(prefix)
        report = read_report(prefix.parent / REPORT_NAME)
        self.assertEqual(report[0]['gene_ids'], 'GZ;GA')
        aliases = rows(str(prefix) + '.exon_aliases.tsv')
        self.assertEqual({r['gene_name'] for r in aliases}, {'Zetamerged'})
        self.assertEqual({r['gene_id'] for r in aliases}, {'GZ'})
        self.assertTrue(runner.validate_database(prefix)[0])
        alignments = self.folder / 'alignments.tsv'
        with alignments.open('w') as handle:
            writer = csv.DictWriter(handle, list(caller.ExonAlignment.__dataclass_fields__), delimiter='\t')
            writer.writeheader()
            for query in ('copy1', 'copy2'):
                for index, sequence in enumerate(self.sequence, 1):
                    length = len(sequence)
                    writer.writerow(asdict(caller.ExonAlignment(query,100,100*index,100*index+length,'+',
                        f'E{index}',length,0,length,length,100,length,length)))
        outputs = []
        for threads in (1, 2):
            output = self.folder / f'calls{threads}.tsv'
            subprocess.run([sys.executable,str(ROOT/'call_genes_from_exon_alignments_v3.py'),
                '-i',str(alignments),'-g',str(f.gff),'-o',str(output),
                '--eligible-exon-info',str(prefix)+'.exon_aliases.tsv','--threads',str(threads),
                '--query-coordinate-mode','local','--shard-storage','disk'],check=True,capture_output=True)
            outputs.append(output)
        self.assertEqual(outputs[0].read_bytes(),outputs[1].read_bytes())
        calls = rows(outputs[0])
        self.assertEqual(len(calls),4)
        self.assertEqual({r['gene_name'] for r in calls},{'Zetamerged'})
        self.assertEqual({r['transcript_id'] for r in calls if r['model_type']=='full_gene'},{'GZ'})
        self.assertEqual({r['model_type'] for r in calls},{'full_gene','transcript'})
        self.assertEqual({r['gene_id'] for r in calls},{'GZ'})
        self.assertEqual({r['tie_count'] for r in calls},{'1'})

    def test_cache_requires_report_policy_and_current_inputs(self):
        f = self.fixture
        f.transcript('G1','One','T1',self.sequence)
        f.write()
        prefix = self.folder/'db'/'exons'
        def ensure():
            return runner.ensure_database(prefix,runner.load_scripts(ROOT),sys.executable,f.genome,f.gff,3,0,99.,False)
        with patch.object(runner,'run_command',wraps=runner.run_command) as command:
            self.assertTrue(ensure())
            self.assertEqual(command.call_count,1)
            self.assertFalse(ensure())
            self.assertEqual(command.call_count,1)
            self.assertEqual(read_report(prefix.parent/REPORT_NAME),[])
            (prefix.parent/REPORT_NAME).unlink()
            self.assertFalse(runner.validate_database(prefix)[0])
            self.assertTrue(ensure())
            self.assertEqual(command.call_count,2)
            f.gff.write_text(f.gff.read_text()+'# changed annotation\n')
            self.assertTrue(ensure())
            self.assertEqual(command.call_count,3)
        manifest_path=Path(str(prefix)+'.manifest.json')
        manifest=json.loads(manifest_path.read_text())
        manifest['identical_mane_policy']='old-policy'
        manifest_path.write_text(json.dumps(manifest))
        self.assertFalse(runner.validate_database(prefix)[0])
