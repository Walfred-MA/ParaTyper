"""Reference gene units connected by >100 bp of shared exonic sequence.

The builder and caller use the same coordinate-based rule. No sequences or
isoforms are discarded. Exons are unioned within each gene before comparing
genes, so annotation duplicates cannot inflate the overlap threshold.
"""
from __future__ import annotations

from collections import defaultdict
import csv
from dataclasses import replace


MIN_SHARED_EXON_BP = 100  # Strictly greater than this value.
MERGE_POLICY = "same_contig_strand_total_shared_exonic_gt100_v1"
REPORT_NAME = "shared_exon_genes.tsv"
REPORT_HEADER = [
    "merged_gene_id", "merged_gene_name", "reference_contig", "strand",
    "member_gene_ids", "member_gene_names", "max_pair_shared_exonic_bp",
]


def merge_shared_exon_rows(rows):
    """Return relabeled exon rows and one report row per merged locus.

    Accepts either builder RawExonRow or caller GffExon records. Components
    are local to a contig and strand; a gene's unrelated alternate locus is
    not relabeled just because its primary locus shares exons with a neighbor.
    Member ordering is deterministic by gene name, then unversioned gene ID.
    """
    rows = list(rows)
    intervals = defaultdict(set)
    names, full_ids = {}, {}
    for row in rows:
        if not row.gene_id or row.end0 <= row.start0:
            continue
        key = (row.gene_id, row.chrom, row.strand)
        intervals[key].add((row.start0, row.end0))
        names[key] = row.gene_name or row.gene_id
        full_ids[key] = row.gene_id_full or row.gene_id

    by_location = defaultdict(list)
    for key, spans in intervals.items():
        merged = []
        for start, end in sorted(spans):
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        for start, end in merged:
            by_location[key[1:]].append((start, end, key))

    parent = {key: key for key in intervals}
    size = {key: 1 for key in intervals}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    # Each gene's intervals are disjoint, so pairwise intersections can be
    # summed directly without counting any base twice across isoforms/exons.
    shared_bases = defaultdict(int)
    for location in sorted(by_location):
        active = []
        for start, end, key in sorted(by_location[location]):
            active = [(a, b, other) for a, b, other in active
                      if b > start]
            for _a, b, other in active:
                if other == key:
                    continue
                overlap = min(b, end) - start
                shared_bases[tuple(sorted((key, other)))] += overlap
            active.append((start, end, key))

    edges = []
    for (key, other), total in sorted(shared_bases.items()):
        if total <= MIN_SHARED_EXON_BP:
            continue
        left, right = find(key), find(other)
        if left != right:
            if size[left] < size[right]:
                left, right = right, left
            parent[right] = left
            size[left] += size[right]
        edges.append((key, total))

    components = defaultdict(list)
    for key in intervals:
        components[find(key)].append(key)
    max_overlap = defaultdict(int)
    for key, overlap in edges:
        root = find(key)
        max_overlap[root] = max(max_overlap[root], overlap)

    replacements, reports = {}, []
    for root, members in components.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda key: (names[key], key[0]))
        gene_id = "&".join(key[0] for key in members)
        gene_name = "&".join(names[key] for key in members)
        gene_id_full = "&".join(full_ids[key] for key in members)
        for key in members:
            replacements[key] = dict(gene_id=gene_id, gene_id_full=gene_id_full,
                                     gene_name=gene_name)
        reports.append(dict(
            merged_gene_id=gene_id, merged_gene_name=gene_name,
            reference_contig=root[1], strand=root[2],
            member_gene_ids=";".join(key[0] for key in members),
            member_gene_names=";".join(names[key] for key in members),
            max_pair_shared_exonic_bp=max_overlap[root],
        ))
    reports.sort(key=lambda row: (row["reference_contig"], row["strand"], row["merged_gene_id"]))
    return [replace(row, **replacements[(row.gene_id, row.chrom, row.strand)])
            if (row.gene_id, row.chrom, row.strand) in replacements else row
            for row in rows], reports


def write_report(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, REPORT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_report(path):
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != REPORT_HEADER:
            raise ValueError(f"{path}: incompatible shared-exon gene report")
        rows = list(reader)
    seen = set()
    for row in rows:
        ids, names = row["member_gene_ids"].split(";"), row["member_gene_names"].split(";")
        keys = {(gene, row["reference_contig"], row["strand"]) for gene in ids}
        if (len(ids) < 2 or len(ids) != len(names) or len(set(ids)) != len(ids)
                or not all(ids) or not all(names) or row["strand"] not in {"+", "-"}
                or not row["reference_contig"] or seen & keys
                or row["merged_gene_id"] != "&".join(ids)
                or row["merged_gene_name"] != "&".join(names)
                or int(row["max_pair_shared_exonic_bp"]) <= MIN_SHARED_EXON_BP):
            raise ValueError(f"{path}: invalid shared-exon gene group")
        seen.update(keys)
    return rows
