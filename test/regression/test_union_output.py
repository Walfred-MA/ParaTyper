"""Compact union output must preserve paired IDs without listing full alternatives."""
import unittest
from dataclasses import replace
from test_exon_transcript_v3 import caller, bed, make_call, resolve


class UnionTests(unittest.TestCase):
    def test_only_unversioned_transcript_and_gene_id_columns(self):
        row = dict(zip(caller.OUTPUT_HEADER, caller.call_to_output_row(make_call('T', gene='G'))))
        self.assertEqual(caller.OUTPUT_HEADER[:3], ['GENE_index', 'transcript_id', 'gene_id'])
        self.assertEqual(row['transcript_id'], 'T')
        self.assertEqual(row['gene_id'], 'G')
        self.assertNotIn('transcript_id_full', row)
        self.assertNotIn('gene_id_full', row)
        self.assertNotIn('alternatives_json', row)

    def test_semicolon_ids_preserve_gene_transcript_pairing(self):
        calls = resolve([make_call('Alpha', gene='G2'), make_call('Zulu', gene='G1')])
        row = dict(zip(caller.OUTPUT_HEADER, caller.assignment_to_output_row(calls)))
        self.assertEqual(row['transcript_id'], 'Zulu;Alpha')
        self.assertEqual(row['gene_id'], 'G1;G2')
        self.assertEqual(row['gene_name'], 'G1;G2')
        self.assertEqual(row['weighted_score'], '4000.000000')
        self.assertEqual(row['exon_query_coordinates'], '100-200,300-400;100-200,300-400')
        self.assertEqual(row['found_exon_numbers'], '1,2;1,2')
        self.assertEqual(row['exon_ids'], 'G1E1,G1E2;G2E1,G2E2')
        self.assertEqual(row['total_aligned_bases'], '200;200')
        self.assertEqual(row['mean_identity'], '100.000000;100.000000')

    def test_tied_exon_indices_keep_per_gene_lists_and_order(self):
        calls = []
        for tid, gene, indices in [('Alpha', 'G2', [3, 8]), ('Zulu', 'G1', [2, 7])]:
            call = make_call(tid, gene=gene)
            calls.append(replace(call, expected_exon_numbers=indices,
                hits=[replace(h, exon_number=i) for h, i in zip(call.hits, indices)]))
        row = dict(zip(caller.OUTPUT_HEADER, caller.assignment_to_output_row(resolve(calls))))
        self.assertEqual(row['gene_id'], 'G1;G2')
        self.assertEqual(row['transcript_id'], 'Zulu;Alpha')
        self.assertEqual(row['found_exon_numbers'], '2,7;3,8')
        record = bed.convert(row, 'test')
        metadata = dict(zip(bed.EXTRA_COLUMNS, record.fields[12:]))
        self.assertEqual(metadata['matched_exon_numbers'], '2,7;3,8')

    def test_transitive_component_is_one_compact_union(self):
        calls = resolve([make_call('A', intervals=((100, 200),)),
                         make_call('B', gene='G2', intervals=((150, 250),)),
                         make_call('C', gene='G3', intervals=((200, 300),))])
        groups = list(caller.iter_assignment_groups(calls))
        self.assertEqual(len(groups), 1)
        row = dict(zip(caller.OUTPUT_HEADER, caller.assignment_to_output_row(groups[0])))
        self.assertEqual(row['exon_query_coordinates'], '100-200;150-250;200-300')
        self.assertEqual(row['gene_name'], 'G1;G2;G3')
        self.assertEqual(row['mean_identity'], '100.000000;100.000000;100.000000')
        self.assertEqual(row['exon_reference_coordinates'], 'chr1:1000-1100;chr1:1000-1100;chr1:1000-1100')
        self.assertNotIn('alternatives_json', row)
        self.assertEqual(row['transcript_id'], 'A;B;C')
        self.assertEqual(row['gene_id'], 'G1;G2;G3')
        record = bed.convert(row, 'test')
        self.assertEqual(record.fields[9:12], ['1', '200,', '0,'])
        metadata = dict(zip(bed.EXTRA_COLUMNS, record.fields[12:]))
        self.assertNotIn('alternatives_json', metadata)
        self.assertEqual(metadata['matched_exon_query_coordinates'], '100-200;150-250;200-300')
        self.assertEqual(metadata['matched_exon_reference_coordinates'], row['exon_reference_coordinates'])

    def test_all_gene_evidence_and_differing_summaries_round_trip(self):
        a = make_call('A', gene='G1')
        b = replace(make_call('B', gene='G2', intervals=((150, 350),), expected=[1, 2],
                              transcript_length=300), protein_bonus=40)
        calls = []
        for call in (a, b):
            calls.append(replace(call, hits=[replace(h, interval_score=caller.MergedIntervalGeneScore(
                h.alignment.query_start, h.alignment.query_end, call.transcript_info.gene_id,
                h.score, h.alignment.exon_id, h.alignment.exon_length)) for h in call.hits]))
        selected = sorted(resolve(calls), key=lambda c: c.transcript_info.gene_id)
        self.assertEqual(len(selected), 2)
        individual = [dict(zip(caller.OUTPUT_HEADER, caller.call_to_output_row(c))) for c in selected]
        row = dict(zip(caller.OUTPUT_HEADER, caller.assignment_to_output_row(selected)))
        global_fields = {'query_id', 'query_contig', 'query_start', 'query_end', 'strand',
                         'tie_group_id', 'tie_count', 'assignment_status', 'pipeline_version'}
        for key in caller.OUTPUT_HEADER:
            if key in global_fields:
                continue
            values = row[key].split(';')
            if len(values) == 1:
                values *= 2
            self.assertEqual(values, [r[key] for r in individual], key)
        self.assertEqual(row['total_aligned_bases'], '200;100')
        self.assertEqual(row['expected_exons'], '2')
        self.assertEqual(row['found_exons'], '2;1')
        self.assertEqual(row['call_status'], 'complete;partial')
        self.assertEqual(row['fraction_expected_found'], '1.000000;0.500000')
        self.assertEqual(row['merged_interval_gene_scores'], '100.000000,100.000000;100.000000')
        self.assertEqual(row['transcript_exon_length'], '200;300')
        self.assertTrue(all(v for k,v in row.items() if k not in {'GENE_index','inserted_exon_numbers','inserted_exon_query_coordinates','insertion_run_unique_exons'}))
        self.assertEqual((row['query_start'], row['query_end']), ('100', '400'))
        record = bed.convert(row, 'sample')
        self.assertEqual(record.fields[9:12], ['1', '300,', '0,'])
        self.assertEqual(record.fields[4], '1000')

    def test_bed_rejects_misaligned_gene_exon_groups(self):
        calls = resolve([make_call('A'), make_call('B', gene='G2')])
        row = dict(zip(caller.OUTPUT_HEADER, caller.assignment_to_output_row(calls)))
        row['exon_reference_coordinates'] = 'chr1:1000-1100;chr1:1000-1100,chr1:2000-2100'
        with self.assertRaisesRegex(ValueError, 'inconsistent lengths or gene groups'):
            bed.convert(row, 'sample')

    def test_interleaved_ties_and_unique_calls_do_not_get_combined(self):
        left = resolve([make_call('A'), make_call('Z', gene='G2')])
        right = resolve([make_call('B', query='other'), make_call('Y', gene='G2', query='other')])
        unique = make_call('unique', query='third')
        groups = list(caller.iter_assignment_groups([left[0], right[0], unique,
                                                    left[1], right[1]]))
        self.assertEqual([len(g) for g in groups], [1, 2, 2])
        self.assertTrue(all(len({c.query_id for c in g}) == 1 for g in groups))

    def test_incomplete_tie_group_is_rejected(self):
        calls = resolve([make_call('A'), make_call('B', gene='G2')])
        with self.assertRaisesRegex(ValueError, 'incomplete tied'):
            list(caller.iter_assignment_groups(calls[:1]))
