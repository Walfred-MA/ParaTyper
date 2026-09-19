#!/usr/bin/env python3
"""Reproduce caller decisions and write a complete SMN1/SMN2 evidence audit.

Uses saved alignments, with the same eligibility, DP, bonuses and greedy
resolution as the production caller. No alignments or annotations are edited.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import html
import json
import math
from pathlib import Path

import call_genes_from_exon_alignments_v3 as c


def tsv(path, rows, fields=None):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def number(value):
    return f'{value:.6f}'


def bounds(call):
    return '-'.join(map(str, c.call_query_bounds(call)))


def table(rows, fields):
    header = '<tr>' + ''.join(f'<th>{html.escape(label)}</th>' for _, label in fields) + '</tr>'
    body = ''.join('<tr>' + ''.join(f'<td>{html.escape(str(row.get(key, "")))}</td>'
                                    for key, _ in fields) + '</tr>' for row in rows)
    return '<div class="scroll"><table><thead>' + header + '</thead><tbody>' + body + '</tbody></table></div>'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    parser.add_argument('--gff', type=Path, default=root / 'test/regression/data/SMN_CHM13_gene_filter.gff3')
    parser.add_argument('--alignments', type=Path, default=root / 'test/regression/data/SMN_CHM13_gene_filter.alignments.tsv')
    parser.add_argument('--eligible-exons', type=Path, default=root / 'test/regression/data/SMN_CHM13_gene_filter.exon_ids.tsv')
    parser.add_argument('--output', type=Path, default=root / 'results/smn_audit')
    parser.add_argument('--no-full-gene-transcripts', action='store_true',
                        help='audit only annotated transcripts, without synthetic full-gene candidates')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    ann = c.parse_gencode_gff3(str(args.gff))
    paralog_report = c.infer_report(str(args.eligible_exons))
    if paralog_report:
        c.apply_identical_paralogs(ann, paralog_report)
    if not args.no_full_gene_transcripts:
        c.write_full_gene_models(str(args.output / 'full_gene_models.tsv'), c.add_full_gene_transcripts(ann))
    eligible = c.load_eligible_exon_ids(str(args.eligible_exons))
    expected = c.truncate_expected_exons(ann, None, eligible)
    raw_alignments = c.read_alignments(str(args.alignments))
    alignments = c.dedup_same_exon_overlaps(c.filter_alignments(
        raw_alignments, ann, None, parse_query_suffix=False))
    smn_tids = {tid for tid, info in ann.transcripts.items() if info.gene_name in {'SMN1', 'SMN2'}}
    smn_exons = defaultdict(list)
    for tid in sorted(smn_tids):
        for exon in ann.exons_by_transcript[tid]:
            smn_exons[(exon.gene_name, exon.exon_id)].append(exon)

    # Build the same complete eligible hit associations as the production caller;
    # retain them here before DP so losing/unselected exons are also auditable.
    grouped = defaultdict(list)
    for aln in alignments:
        for tid, num, chrom, start, end in ann.exon_to_transcripts.get(aln.exon_id, []):
            if num in expected.get(tid, []):
                grouped[(aln.query_id, aln.strand, tid)].append(
                    c.TranscriptHit(num, aln, chrom, start, end))
    c.assign_merged_interval_gene_scores(grouped, ann)
    hit_index = {}
    block_hits = defaultdict(dict)
    for (query, strand, tid), hits in grouped.items():
        if tid not in smn_tids:
            continue
        gene = ann.transcripts[tid].gene_name
        for hit in hits:
            a = hit.alignment
            hit_index[(gene, a.exon_id, query, strand, a.query_start, a.query_end)] = hit
            block_hits[(query, strand, *hit.competition_interval, gene)][id(a)] = hit
    locus_keys = sorted({(q, s) for q, s, _, _, _ in block_hits})
    block_names = {}
    for query, strand in locus_keys:
        intervals = sorted({(s, e) for q, st, s, e, _ in block_hits if (q, st) == (query, strand)})
        for i, (s, e) in enumerate(intervals, 1):
            block_names[(query, strand, s, e)] = ('plus' if strand == '+' else 'minus') + f'.B{i:02d}'

    def blocks(call):
        return [block_names.get((call.query_id, call.strand, s.start, s.end), f'{s.start}-{s.end}')
                for s in call.merged_interval_scores]

    batches = []
    if args.no_full_gene_transcripts:
        raw_calls = c.build_transcript_calls(alignments, ann, expected, 10., 2., 1., 10)
        batches.append(dict(stage="transcript", parent="", calls=raw_calls, prefer_mane=True))
        reference = c.resolve_call_overlaps(raw_calls, True)
    else:
        reference, _count = c.build_hierarchical_calls(
            alignments, ann, expected, 10., 2., 1., 10, True, diagnostic_batches=batches)
        raw_calls = [call for batch in batches for call in batch['calls']]
    stage_by_call = {id(call): batch['stage'] for batch in batches for call in batch['calls']}
    chain_ids = {id(call): f'C{i:03d}' for i, call in enumerate(raw_calls, 1)}
    initial = {chain_ids[id(call)]: call for call in raw_calls}
    trace_events = []
    selection = {}
    selected_all = []
    decisions = defaultdict(list)
    groups = defaultdict(list)
    for batch_index, batch in enumerate(batches):
        for call in batch['calls']:
            groups[(batch_index, call.query_id, call.strand)].append(call)
    for (batch_index, query, strand), group in sorted(groups.items()):
        batch = batches[batch_index]
        trace = []
        selected_all.extend(c.resolve_transcript_overlaps_sparse(group, batch['prefer_mane'], trace=trace))
        owners = {}
        step = 0
        for event in trace:
            if event['event'] == 'isoform_tiebreak':
                loser, winner = event['call'], event['winner']
                cid, winner_id = chain_ids[id(loser)], chain_ids[id(winner)]
                reason = ('longer annotated transcript' if winner.transcript_info.exon_length > loser.transcript_info.exon_length
                          else 'stable transcript ID at equal length')
                decisions[cid].append(f'Same-gene score tie: prefer {winner_id} {winner.transcript_id} '
                    f'({winner.transcript_info.exon_length} bp) over {loser.transcript_id} '
                    f'({loser.transcript_info.exon_length} bp), by {reason}; defer to check residual evidence.')
                trace_events.append(dict(event='isoform_tiebreak', query=query, strand=strand, after_step=step,
                    chain_id=cid, preferred_chain=winner_id, reason=reason,
                    transcript_exon_length=loser.transcript_info.exon_length,
                    preferred_transcript_exon_length=winner.transcript_info.exon_length))
            elif event['event'] == 'select':
                step += 1
                for original, selected in zip(event['calls'], event['selected']):
                    cid = chain_ids[id(original)]
                    selection[cid] = (selected, step)
                    for hit in selected.hits:
                        owners.setdefault(hit.competition_interval, []).append(cid)
                    decisions[cid].append(f'Selected at step {step}; score {number(selected.weighted_score)}; '
                                          f'{selected.found_expected_count}/{selected.expected_count} expected exons; '
                                          f'{selected.tie_count} alternative(s).')
                    trace_events.append(dict(event='select', query=query, strand=strand, step=step,
                                             chain_id=cid, weighted_score=selected.weighted_score,
                                             intervals=blocks(selected), tie_group=selected.tie_group_id))
            else:
                call, remaining = event['call'], event['remaining']
                cid = chain_ids[id(call)]
                keep = {h.competition_interval for h in remaining.hits} if remaining else set()
                lost = {h.competition_interval for h in call.hits} - keep
                blocking = sorted({owner for interval in lost for owner in owners.get(interval, [])})
                reasons = []
                for owner in blocking:
                    winner = selection[owner][0]
                    reason = ('MANE priority' if c.transcript_call_rank_key(winner, True)[0] > c.transcript_call_rank_key(call, True)[0]
                              else 'equal-score isoform/gene competition' if math.isclose(
                                  winner.weighted_score, call.weighted_score, rel_tol=1e-12, abs_tol=1e-12)
                              else 'higher ranked score')
                    reasons.append(f'{owner} {winner.transcript_info.gene_name}/{winner.transcript_id} ({reason})')
                decisions[cid].append('Removed occupied blocks ' + ','.join(
                    block_names.get((query, strand, s, e), f'{s}-{e}') for s, e in sorted(lost)) +
                    '; occupied by ' + '; '.join(reasons) +
                    (f'; residual score {number(remaining.weighted_score)}.' if remaining else '; no evidence remains.'))
                trace_events.append(dict(event='trim', query=query, strand=strand, after_step=step,
                                         chain_id=cid, blocked_by=blocking, removed_intervals=sorted(lost),
                                         remaining_score=remaining.weighted_score if remaining else None))
                if remaining:
                    chain_ids[id(remaining)] = cid
        for event in trace_events:
            if 'stage' not in event:
                event['stage'] = batch['stage']
                event['parent'] = batch['parent']
    # Independently traced competitions must reproduce both production stages.
    def signature(call):
        return (call.transcript_id, call.query_id, call.strand, call.raw_score,
                call.weighted_score, call.tie_group_id,
                tuple((h.exon_number, h.alignment.exon_id, h.alignment.query_start,
                       h.alignment.query_end, h.is_insertion) for h in call.hits))
    assert sorted(map(signature, selected_all)) == sorted(map(signature, reference))
    selected_all = reference
    c.write_transcript_calls(selected_all, str(args.output / 'reproduced_assignments.tsv'),
                             'extended', False, False, False)
    (args.output / 'greedy_decisions.json').write_text(json.dumps(trace_events, indent=2) + '\n')

    candidate_rows = []
    for cid, call in initial.items():
        if call.transcript_id not in smn_tids:
            continue
        selected, step = selection.get(cid, (None, ''))
        retained = len(selected.hits) if selected else 0
        outcome = ('lost' if selected is None else 'selected residual' if retained < len(call.hits)
                   else 'selected intact')
        if selected and selected.tie_count > 1:
            outcome += ' (tied)'
        candidate_rows.append(dict(chain_id=cid, stage=stage_by_call[id(call)], gene=call.transcript_info.gene_name,
            transcript=call.transcript_id,
            transcript_type=call.transcript_info.transcript_type, MANE=int(call.transcript_info.is_mane),
            model_type=call.transcript_info.model_type,
            transcript_exon_length=call.transcript_info.exon_length,
            query=call.query_id, strand=call.strand, query_interval=bounds(call),
            expected_exons=call.expected_count, found_exons=call.found_expected_count,
            exon_numbers=','.join(map(str, call.found_exon_numbers)),
            exon_ids=','.join(h.alignment.exon_id for h in call.hits),
            scored_blocks=','.join(blocks(call)), raw_dp_score=number(call.raw_score),
            inserted_exons=call.inserted_exon_count, insertion_penalty=call.insertion_penalty,
            inserted_query_blocks=len(call.insertion_intervals),
            insertion_run_unique_exons=','.join(str(len(run)) for run in call.insertion_run_exons),
            protein_multiplier=call.protein_bonus if call.transcript_info.transcript_type == 'protein_coding' else 1,
            complete_multiplier=2 if call.expected_count and call.missing_count == 0 else 1,
            weighted_score=number(call.weighted_score), outcome=outcome, selected_step=step,
            selected_score=number(selected.weighted_score) if selected else '',
            selected_exons=retained, selected_blocks=','.join(blocks(selected)) if selected else '',
            tie_group=selected.tie_group_id if selected else '',
            explanation=' '.join(decisions[cid])))
    tsv(args.output / 'isoform_decisions.tsv', candidate_rows)

    block_rows = []
    for (query, strand, start, end, gene), entries in sorted(block_hits.items()):
        hits = list(entries.values())
        source = min((h.alignment for h in hits), key=c.interval_score_source_key)
        score = hits[0].interval_score
        assert (score.source_exon_id, score.source_exon_length, score.score) == (source.exon_id, source.exon_length, source.exon_score)
        block_rows.append(dict(block=block_names[(query, strand, start, end)], query=query, strand=strand,
            start=start, end=end, gene=gene, supplier_exon=source.exon_id, supplier_length=source.exon_length,
            supplier_identical_bases=source.identical_bases, supplier_score=number(score.score),
            former_max_score=number(max(h.alignment.exon_score for h in hits)),
            eligible_exon_hits=len(hits), exon_ids=','.join(sorted({h.alignment.exon_id for h in hits}))))
    tsv(args.output / 'merged_interval_scores.tsv', block_rows)

    original_rows = list(csv.DictReader(args.alignments.open(), delimiter='\t'))
    cigar_map = {(r['exon_id'], r['query_id'], r['strand'], int(r['query_start']), int(r['query_end'])): r
                 for r in original_rows}
    exon_rows = []
    for (gene, eid), annotations in sorted(smn_exons.items()):
        e = annotations[0]
        for query, strand in locus_keys:
            hits = [hit for key, hit in hit_index.items() if key[:4] == (gene, eid, query, strand)]
            for hit in hits or [None]:
                a = hit.alignment if hit else None
                call_members = [cid for cid, call in initial.items() if call.transcript_info.gene_name == gene
                                and a and any(h.alignment is a for h in call.hits)]
                selected_members = [cid for cid, (call, _) in selection.items()
                                    if call.transcript_info.gene_name == gene and a
                                    and any(h.alignment is a for h in call.hits)]
                raw = cigar_map.get((eid, query, strand, a.query_start, a.query_end), {}) if a else {}
                interval = hit.interval_score if hit else None
                supplied_members = [cid for cid, (call, _) in selection.items()
                    if interval and call.transcript_info.gene_name == gene
                    and call.query_id == query and call.strand == strand
                    and any(s.start == interval.start and s.end == interval.end and s.source_exon_id == eid
                            for s in call.merged_interval_scores)]
                exon_rows.append(dict(gene=gene, exon_id=eid, exon_id_full=e.exon_id_full,
                    reference=f'{e.chrom}:{e.start0}-{e.end0}:{e.strand}', annotated_length=e.length,
                    transcripts=','.join(sorted({x.transcript_id for x in annotations})),
                    MANE_exon_numbers=','.join(map(str, sorted({x.exon_number for x in annotations if x.is_mane}))),
                    eligible=int(eid in eligible), query=query, strand=strand,
                    alignment_status='aligned' if a else 'not in exon database' if eid not in eligible else 'no retained alignment',
                    query_start=a.query_start if a else '', query_end=a.query_end if a else '',
                    score_length=a.exon_length if a else '', identical_bases=a.identical_bases if a else '',
                    length_minus_identical=a.exon_length-a.identical_bases if a else '',
                    identity=a.percent_identity if a else '', coverage=a.exon_coverage if a else '',
                    individual_exon_score=number(a.exon_score) if a else '', AS=a.AS if a else '',
                    NM=raw.get('NM', ''), cigar=raw.get('cigar', ''),
                    merged_block=block_names[(query, strand, *hit.competition_interval)] if hit else '',
                    gene_interval_score=number(interval.score) if interval else '',
                    score_supplier=interval.source_exon_id if interval else '',
                    score_supplier_length=interval.source_exon_length if interval else '',
                    supplies_score=int(eid == interval.source_exon_id) if interval else '',
                    dp_chain_ids=','.join(call_members), selected_chain_ids=','.join(selected_members),
                    supplies_score_to_selected_chain_ids=','.join(supplied_members)))
    tsv(args.output / 'all_exons.tsv', exon_rows)

    mane_rows = []
    diagnostic_mane = c.build_calls_from_grouped_hits(
        {key: hits for key, hits in grouped.items() if ann.transcripts[key[2]].is_mane},
        ann, expected, 10., 2., 1., 10)
    for index, call in enumerate(diagnostic_mane, 1):
        row = next((r for r in candidate_rows if signature(initial[r['chain_id']])[:5] == signature(call)[:5]
                    and bounds(initial[r['chain_id']]) == bounds(call)), None)
        if row is None:
            row = dict(chain_id=f'D{index:03d}', gene=call.transcript_info.gene_name,
                       transcript=call.transcript_id, outcome='Diagnostic only; not admitted by gene stage')
        seen = set()
        for hit in sorted(call.hits, key=lambda h: h.exon_number):
            a, interval = hit.alignment, hit.interval_score
            contribution = 0 if hit.competition_interval in seen else hit.score
            seen.add(hit.competition_interval)
            mane_rows.append(dict(chain_id=row['chain_id'], gene=row['gene'], transcript=row['transcript'],
                query=call.query_id, strand=call.strand, exon_number=hit.exon_number, exon_id=a.exon_id,
                query_start=a.query_start, query_end=a.query_end, length=a.exon_length,
                identical_bases=a.identical_bases, individual_score=number(a.exon_score),
                block=block_names[(call.query_id, call.strand, *hit.competition_interval)],
                supplier=interval.source_exon_id, supplier_length=interval.source_exon_length,
                interval_score=number(hit.score), dp_contribution=number(contribution),
                reason='same block already counted' if contribution == 0 else 'first exon in this scored block',
                isoform_result=row['outcome']))
        assert abs(sum(float(r['dp_contribution']) for r in mane_rows if r['chain_id'] == row['chain_id']) - call.raw_score) < 1e-5
    tsv(args.output / 'mane_exon_comparison.tsv', mane_rows)

    inventory = []
    for tid in sorted(smn_tids):
        info = ann.transcripts[tid]
        candidates = [r for r in candidate_rows if r['transcript'] == tid]
        inventory.append(dict(gene=info.gene_name, transcript=tid,
            MANE=int(info.is_mane), transcript_type=info.transcript_type,
            model_type=info.model_type,
            transcript_exon_length=info.exon_length,
            annotated_exons=len({e.exon_number for e in ann.exons_by_transcript[tid]}), eligible_expected_exons=len(expected[tid]),
            candidate_chains=len(candidates), selected_chains=sum(r['outcome'] != 'lost' for r in candidates),
            chain_ids=','.join(r['chain_id'] for r in candidates),
            note='See isoform_decisions.tsv' if candidates else 'No candidate admitted inside a selected gene parent (or no eligible alignment)'))
    tsv(args.output / 'transcript_inventory.tsv', inventory)

    paths = [args.gff, args.alignments, args.eligible_exons, Path(c.__file__), Path(__file__)]
    if paralog_report:
        paths.append(paralog_report)
    manifest = dict(version=c.PIPELINE_VERSION, coordinate_system='0-based half-open',
        parameters=dict(prefer_mane=True, full_gene_transcripts=not args.no_full_gene_transcripts,
                        protein_bonus=10, gene_protein_bonus=1, complete_bonus=2, mane_bonus=1,
                        max_chains_per_transcript=10, score_supplier='longest reference exon per gene and merged interval',
                        length_tie_break='highest exon score, then stable exon ID and query coordinates',
                        isoform_tie_break='longest full annotated spliced transcript within gene; then transcript ID',
                        reported_ties='different genes only', gene_insertion_cost=c.GENE_INSERTION_COST,
                        max_gene_insertion_unique_exons=c.MAX_GENE_INSERTION_UNIQUE_EXONS,
                        insertion_count='unique full-gene exon numbers per consecutive insertion run',
                        calling='full genes first; real isoforms compete within each gene assignment'),
        inputs={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        counts=dict(input_alignments=len(raw_alignments), retained_alignments=len(alignments),
                    annotated_smn_transcripts=len(inventory), distinct_smn_exons=len(smn_exons),
                    exon_locus_rows=len(exon_rows), smn_candidate_chains=len(candidate_rows),
                    smn_selected_alternatives=sum(r['outcome'] != 'lost' for r in candidate_rows)))
    (args.output / 'provenance.json').write_text(json.dumps(manifest, indent=2) + '\n')

    leading = {}
    for row in candidate_rows:
        if not (row['MANE'] if args.no_full_gene_transcripts else row['model_type'] == 'full_gene'):
            continue
        key = (row['strand'], row['gene'])
        if key not in leading or c.transcript_call_rank_key(initial[row['chain_id']], True) > c.transcript_call_rank_key(initial[leading[key]['chain_id']], True):
            leading[key] = row
    main_candidates = list(leading.values())
    main_table = table(sorted(main_candidates, key=lambda r: (r['strand'], r['gene'])), [
        ('strand','Locus strand'), ('gene','Gene'), ('transcript','Candidate'), ('model_type','Model'), ('found_exons','Exons found'),
        ('raw_dp_score','DP score'), ('weighted_score','Weighted score'), ('outcome','Decision')])
    sections = []
    for strand in ('+', '-'):
        comparisons = []
        for num in range(1, 10):
            entries = {r['gene']:r for r in mane_rows if r['strand'] == strand and r['exon_number'] == num}
            if not {'SMN1', 'SMN2'} <= entries.keys():
                continue
            a, b = entries['SMN1'], entries['SMN2']
            comparisons.append(dict(exon=num, interval=f"{a['query_start']}-{a['query_end']}",
                a=f"{a['identical_bases']}/{a['length']} → {a['individual_score']}",
                b=f"{b['identical_bases']}/{b['length']} → {b['individual_score']}",
                block=a['block'], ga=a['interval_score'], gb=b['interval_score'],
                contribution=f"{a['dp_contribution']} / {b['dp_contribution']}"))
        sections.append(f'<h2>Every MANE exon: {html.escape(strand)} locus</h2>' +
            '<p>Individual scores use each MANE exon’s own alignment. DP uses the longest exon’s score from that gene’s entire merged block. '
            'An interval contributes only once; exon 9 adds zero because exon 8 already uses that block.</p>' + table(comparisons, [
            ('exon','GFF exon #'), ('interval','Query interval'), ('a','SMN1 identical/length → own score'),
            ('b','SMN2 identical/length → own score'), ('block','Merged block'),
            ('ga','SMN1 block score'), ('gb','SMN2 block score'), ('contribution','DP contribution SMN1 / SMN2')]))
    suppliers = table(block_rows, [('block','Block'), ('start','Start'), ('end','End'), ('gene','Gene'),
        ('supplier_exon','Longest exon supplying score'), ('supplier_length','Length'),
        ('supplier_identical_bases','Identical bases'), ('supplier_score','Current score'), ('former_max_score','Previous max score')])
    per_strand = {strand:{r['gene']:r for r in main_candidates if r['strand'] == strand} for strand in ('+', '-')}
    differences = []
    for strand, candidates in per_strand.items():
        winner = next((r for r in candidates.values() if r['outcome'] != 'lost'), None)
        loser = next((r for r in candidates.values() if r['outcome'] == 'lost'), None)
        if winner is None or loser is None:
            differences.append(f'<p>{html.escape(strand)} locus: see the candidate table and exact selection trace.</p>')
            continue
        different_blocks = []
        for row in block_rows:
            if row['strand'] != strand or row['gene'] != winner['gene'] or row['block'] not in winner['scored_blocks'].split(','):
                continue
            rival = next((r for r in block_rows if r['strand'] == strand and r['block'] == row['block'] and r['gene'] == loser['gene']), None)
            if rival is None:
                continue
            if row['supplier_score'] != rival['supplier_score']:
                different_blocks.append(f"{row['block']}: {winner['gene']} {row['supplier_score']} from {row['supplier_exon']} "
                    f"({row['supplier_identical_bases']}/{row['supplier_length']}) versus {loser['gene']} {rival['supplier_score']} "
                    f"from {rival['supplier_exon']} ({rival['supplier_identical_bases']}/{rival['supplier_length']})")
        differences.append('<p><b>' + html.escape(f"{strand} locus: {winner['gene']} wins by " +
            number(float(winner['weighted_score']) - float(loser['weighted_score'])) + ' weighted points.') + '</b> ' +
            html.escape('; '.join(different_blocks)) + '.</p>')
    file_links = ''.join(f'<li><a href="{name}">{label}</a></li>' for name, label in [
        ('all_exons.tsv','Every annotated SMN1/SMN2 exon, both loci: scores, coordinates, CIGAR and isoform membership'),
        ('mane_exon_comparison.tsv','36 MANE exon observations with actual DP contributions'),
        ('merged_interval_scores.tsv','All merged intervals and their longest-exon score suppliers'),
        ('isoform_decisions.tsv','Every generated SMN1/SMN2 isoform chain: score and exact selection/removal explanation'),
        ('transcript_inventory.tsv','All annotated SMN1/SMN2 transcripts, including those with no eligible chain'),
        ('greedy_decisions.json','Exact greedy selection and trimming trace for all genes'),
        ('reproduced_assignments.tsv','Final calls: one union record per tied competing assignment'),
        ('provenance.json','Parameters and SHA-256 hashes for inputs and caller')])
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>SMN1 / SMN2 exon audit</title>
<style>body{font:15px/1.55 system-ui,sans-serif;color:#172b3a;background:#f7fafc;margin:0}main{max-width:1450px;margin:auto;padding:35px}
h1{font-size:32px}h2{margin-top:36px}p{max-width:1100px}.card{padding:20px 26px;background:#e6f2f6;border-left:5px solid #126782}
table{border-collapse:collapse;width:100%;background:white;font-size:13px}th,td{padding:9px 11px;text-align:left;border-bottom:1px solid #dce4eb;vertical-align:top}
th{background:#15394a;color:white;position:sticky;top:0}tr:nth-child(even){background:#f0f5f8}.scroll{overflow:auto;max-height:640px;border:1px solid #dce4eb}
code{background:#e7edf1;padding:2px 4px}a{color:#086b8c}details{margin:16px 0}summary{cursor:pointer;font-weight:650}small{color:#526270}</style><main>
<h1>SMN1 / SMN2: exon scores and assignment decisions</h1>'''
    page += f'<p>ParaTyper {c.PIPELINE_VERSION} · CHM13 · saved exon alignments · coordinates are 0-based, half-open.</p>'
    page += '<div class="card"><b>Longest-exon scoring gives SMN1 at the + locus and SMN2 at the − locus.</b> '
    page += 'The table below includes the leading MANE or synthetic full-gene model for each gene and strand. '
    page += 'Within a gene, tied isoforms prefer the longest annotated transcript. Only ties between genes produce an ambiguous union record. '
    page += 'This reproduces the algorithm; it is not independent biological validation.</div>'
    page += main_table + '<h2>Why these genes win</h2>' + ''.join(differences)
    page += '<p>The former rule chose the highest normalized exon score within each gene and block. A short perfect alternative could therefore hide mismatches in longer exons. '
    page += 'The new rule first chooses the longest eligible aligned exon for that gene, then uses that exon’s normalized score. '
    page += 'At equal lengths it chooses the better score, then a stable exon ID. Different genes can still use different-length suppliers and different sequence coverage.</p>'
    page += '<h2>Score and selection rules</h2><ol><li>Individual exon score = <code>100 × (L − 4 × (L − identical_bases)) / L</code>. '
    page += 'L is the reference exon length. Missing/unmatched exon bases reduce this score. BLAST AS and percent identity are shown for auditing but are not the DP score.</li>'
    page += '<li>Merge overlapping eligible exon hits transitively on the same query and strand. Store a separate longest-exon score for every gene.</li>'
    page += '<li>DP chooses a chain with increasing exon numbers and nonoverlapping query hits in transcript orientation. Add each distinct merged interval once. '
    page += 'All isoforms of a gene reuse its block score, including a score supplied by an exon outside the chosen isoform.</li>'
    page += '<li>First compare only full-gene models, using matched exon scores minus 50 per unique full-gene exon number in each consecutive insertion run. Repeats within a run count once; runs over 20 unique exons are rejected. Then compete real isoforms inside each selected gene, pooling isoforms from tied genes and applying MANE priority. Gene-stage scores do not receive a protein-coding bonus. Only real transcript scores use the protein-coding multiplier (default ×10). Both stages retain the complete-model multiplier (default ×2). Parent and child rows share GENE_index. Synthetic models remain labeled full_gene, with MANE flag 0.</li>'
    page += '<li>Within each equal-rank competing group, select one isoform per gene using full annotated exon length (sum of exon lengths, excluding introns and deduplicating repeated annotation coordinates). '
    page += 'This includes exons absent from the database. Equal lengths use transcript ID. Length does not resolve cross-gene ties.</li>'
    page += '<li>Select the remaining equally ranked gene alternatives together. Consume the union of their merged intervals; trim remaining candidates, including losing isoforms, and recalculate their scores and completeness. Repeat.</li></ol>'
    page += ''.join(sections) + '<h2>Every merged block and its score supplier</h2>' + suppliers
    page += '<h2>Selected calls and residual evidence</h2><p>Selected models consume their supported merged blocks. '
    page += 'Any remaining blocks can support partial calls. These are residual exon evidence, not additional complete gene copies. '
    page += 'The current caller has no minimum fraction-of-transcript threshold, so one-exon partial assignments remain. '
    page += 'For example, an SMN1 residual label at the − locus does not override the main SMN2 assignment.</p>'
    page += table([r for r in candidate_rows if r['outcome'] != 'lost'], [('chain_id','Chain'), ('strand','Strand'), ('gene','Gene'),
        ('transcript','Transcript'), ('stage','Stage'), ('model_type','Model'), ('MANE','MANE'), ('weighted_score','Before competition'), ('selected_score','After competition'),
        ('selected_exons','Retained exons'), ('expected_exons','Expected'), ('outcome','Decision'), ('selected_blocks','Retained blocks')])
    page += '<h2>All isoform decisions</h2><p>Expand the complete table or download the TSV. A lost chain has no unoccupied exon evidence left; '
    page += 'a selected residual chain lost some blocks to an earlier choice. Step numbers are local to each query/strand.</p><details><summary>Show every generated SMN1/SMN2 chain</summary>'
    page += table(candidate_rows, [('chain_id','Chain'), ('gene','Gene'), ('transcript','Transcript'), ('strand','Strand'),
        ('MANE','MANE'), ('transcript_exon_length','Annotated exon length'), ('weighted_score','Initial score'),
        ('outcome','Outcome'), ('explanation','Exact decision trace')]) + '</details>'
    page += '<h2>Every annotated exon</h2><p>Each exon is listed at both loci. Blank scores mean no eligible alignment was available, '
    page += 'not a score of zero. “Selected exon chains” means this exon is part of the reported structure; “Supplies score to” can include other isoforms that reuse its gene/block score. '
    page += 'Chain IDs correspond to the decision table above; all transcript memberships and CIGAR strings are included in the downloadable exon TSV.</p>'
    page += '<details><summary>Show all exon/locus observations</summary><label>Find an exon, gene, block or chain: '
    page += '<input id="exon-search" type="search" placeholder="e.g. ENSE00002045198"></label><div id="exon-table">'
    page += table(exon_rows, [('gene','Gene'), ('exon_id','Exon ID'), ('strand','Locus strand'),
        ('MANE_exon_numbers','MANE exon #'), ('annotated_length','Length'), ('alignment_status','Alignment'),
        ('identical_bases','Identical'), ('individual_exon_score','Own score'), ('merged_block','Block'),
        ('gene_interval_score','Used block score'), ('score_supplier','Supplier exon'),
        ('selected_chain_ids','Selected exon chains'), ('supplies_score_to_selected_chain_ids','Supplies score to')]) + '</div></details>'
    page += '<h2>Coverage and downloadable evidence</h2><p>' + html.escape(
        f"The audit includes {len(inventory)} annotated SMN1/SMN2 transcripts and {len(smn_exons)} distinct gene/exon IDs; "
        f"{len(exon_rows)} exon/locus records include missing or database-ineligible exons. All {len(raw_alignments)} input alignments across all genes "
        f"were considered; {len(alignments)} survive caller filtering. {len(candidate_rows)} SMN1/SMN2 chains were generated.") + '</p><ul>' + file_links + '</ul>'
    page += '<p>Real transcript exon numbers come from the GFF; synthetic full-gene numbers label the gene-specific union blocks in strand order. Neither is asserted to match clinical exon nomenclature. '
    page += 'Union rows preserve each gene’s exon and alignment fields in the existing columns: semicolons separate genes, and commas separate exons within a gene. All fields follow the gene/transcript ID order. Equal summaries appear once; differing summaries use semicolons. Query bounds and BED geometry describe the union, while per-gene exon coordinates retain the separate structures. There is no expanded JSON list. '
    page += 'Competition is based on merged exon blocks on the same strand, not overlap of entire gene spans; opposite-strand antisense calls and calls within unoccupied introns can coexist.</p>'
    page += '<p><small>Regenerate from the project root: <code>python3 transcript/audit_smn_calls.py</code>. Existing earlier output directories are preserved.</small></p></main>'
    page += '''<script>document.getElementById('exon-search').addEventListener('input',function(){
const query=this.value.toLowerCase();document.querySelectorAll('#exon-table tbody tr').forEach(row=>{
row.hidden=!row.textContent.toLowerCase().includes(query);});});</script></html>'''
    (args.output / 'report.html').write_text(page)
    print(json.dumps(manifest['counts'], indent=2))
    print(args.output / 'report.html')


if __name__ == '__main__':
    main()
