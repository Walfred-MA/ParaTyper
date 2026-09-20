#!/usr/bin/env python3
"""
build_exon_blastdb_v2.py

Build an exon-centric BLASTN database from a reference genome FASTA and a
GENCODE-style GFF3 annotation.

Genes with identical complete spliced MANE DNA sequence sets are merged before
exon extraction. Keep every isoform of the first gene in GFF3 order and rename
it <first_gene_name>merged; report the group in indenticalparalogs.tsv.

The exon FASTA header has 10 tab-separated fields:

    exon_id coordinates transcript_id transcript_index ifmane_transcript gene_id gene_name ifproteincoding left_anchor_length right_anchor_length

Coordinates are zero-based, right-open genomic coordinates with strand appended:

    chrom:start-end+
    chrom:start-end-

Each original exon shorter than 150 bp receives ceil((150 - length) / 2) bp
of flanking reference sequence on each side; longer exons have no anchors.
Actual oriented anchor lengths are stored because contig boundaries clip flanks.

Overlapping core exons of the same gene, contig and strand are unioned before
BLAST. Original exon boundaries and transcript associations remain in the alias
table, including their offsets within the union query. Identical anchored union
queries can also share one search without discarding their original aliases.

Outputs:
    <out>.exons.fa       exon FASTA used to build the BLAST database
    <out>.seq            BLAST ordinal-name map; first FASTA header token only
    <out>.exon_info.tsv  complete metadata used by the alignment converter
    <out>.exon_aliases.tsv  original retained exon metadata
    <out>.manifest.json  input identities, build settings and merge policy
    indenticalparalogs.tsv  merged gene groups, in the output prefix directory
    <out>.*              BLAST database files, unless --no-makeblastdb is used
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from identical_paralogs import (DATABASE_FORMAT, EXON_QUERY_MERGE_POLICY, MERGE_POLICY, REPORT_NAME,
                                file_identity, report_digest, write_report)
import shared_exon_genes

RC_TABLE = str.maketrans(
    "ACGTNacgtnRYKMSWBDHVrykmswbdhv",
    "TGCANtgcanYRMKSWVHDByrmkswvhdb",
)


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_attrs(attr_text: str) -> Dict[str, str]:
    """Parse GENCODE/GFF3 key=value attributes; tolerate simple GTF too."""
    attrs: Dict[str, str] = {}
    for field in attr_text.rstrip().split(";"):
        field = field.strip()
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
        elif " " in field:
            key, value = field.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def strip_version(identifier: str) -> str:
    """Strip Ensembl version suffix: ENSE000...1.4 -> ENSE000...1."""
    if "." not in identifier:
        return identifier
    head, tail = identifier.rsplit(".", 1)
    if tail.isdigit():
        return head
    return identifier


def attr_has_mane(attrs: Dict[str, str]) -> bool:
    for key, value in attrs.items():
        if "MANE" in key.upper() or "MANE" in value.upper():
            return True
    return False


def uniq_preserve_order(values: Iterable[str]) -> List[str]:
    seen = set()
    out = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


@dataclass
class RawExonRow:
    chrom: str
    start0: int
    end0: int
    strand: str
    exon_id_full: str
    exon_id: str
    transcript_id_full: str
    transcript_id: str
    exon_number: str
    gene_id_full: str
    gene_id: str
    gene_name: str
    transcript_type: str
    is_mane: bool

    @property
    def coordinates(self) -> str:
        return f"{self.chrom}:{self.start0}-{self.end0}{self.strand}"


@dataclass
class ExonRecord:
    chrom: str
    start0: int
    end0: int
    strand: str
    exon_id_full: str
    exon_id: str
    transcript_ids_full: List[str] = field(default_factory=list)
    transcript_ids: List[str] = field(default_factory=list)
    transcript_indices: List[str] = field(default_factory=list)
    gene_ids_full: List[str] = field(default_factory=list)
    gene_ids: List[str] = field(default_factory=list)
    gene_names: List[str] = field(default_factory=list)
    transcript_types: List[str] = field(default_factory=list)
    any_mane: bool = False

    @property
    def coordinates(self) -> str:
        return f"{self.chrom}:{self.start0}-{self.end0}{self.strand}"

    @property
    def transcript_id_field(self) -> str:
        return ",".join(uniq_preserve_order(self.transcript_ids_full))

    @property
    def transcript_index_field(self) -> str:
        # Same order as transcript_id_field for rows merged from the same exon.
        pairs = []
        seen = set()
        for tid, idx in zip(self.transcript_ids_full, self.transcript_indices):
            key = (tid, idx)
            if key in seen:
                continue
            seen.add(key)
            pairs.append(idx)
        return ",".join(pairs)

    @property
    def gene_id_field(self) -> str:
        return ",".join(uniq_preserve_order(self.gene_ids_full))

    @property
    def gene_name_field(self) -> str:
        return ",".join(uniq_preserve_order(self.gene_names))

    @property
    def ifproteincoding(self) -> str:
        return "1" if "protein_coding" in self.transcript_types else "0"

    def header_fields(self, left_anchor_length: int, right_anchor_length: int) -> List[str]:
        return [
            self.exon_id_full,
            self.coordinates,
            self.transcript_id_field,
            self.transcript_index_field,
            "1" if self.any_mane else "0",
            self.gene_id_field,
            self.gene_name_field,
            self.ifproteincoding,
            str(left_anchor_length),
            str(right_anchor_length),
        ]


@dataclass
class ExtractedExon:
    record: ExonRecord
    sequence: str
    unmasked: int
    left_anchor_length: int
    right_anchor_length: int

    @property
    def exon_length(self) -> int:
        return self.record.end0 - self.record.start0


@dataclass
class QueryExonAlias:
    exon: ExtractedExon
    core_start0: int
    core_end0: int
    anchor_start0: int
    anchor_end0: int
    overlap_merged: bool = False


@dataclass
class ExonQuery:
    target: ExtractedExon
    aliases: List[QueryExonAlias]


def genomic_query_span(exon: ExtractedExon) -> Tuple[int, int]:
    rec = exon.record
    left, right = exon.left_anchor_length, exon.right_anchor_length
    if rec.strand == "-":
        left, right = right, left
    return rec.start0 - left, rec.end0 + right


def exon_query_gene_key(record: ExonRecord) -> Tuple[str, Tuple[str, ...]]:
    """Prefer gene names for query unions; fall back to IDs for unnamed genes."""
    names = tuple(sorted({name.strip() for name in record.gene_names
                          if name.strip() not in {"", "."}}))
    if names:
        return "gene_name", names
    ids = tuple(sorted({gene.strip() for gene in record.gene_ids
                        if gene.strip() not in {"", "."}}))
    return "gene_id", ids


def merge_overlapping_exon_queries(
    exons: Sequence[ExtractedExon], merge_overlapping: bool = True
) -> List[ExonQuery]:
    """Union intersecting core intervals within a gene name/contig/strand.

    Unnamed genes fall back to gene ID. Shared reference-gene units (e.g.
    GA&GB) already have one name. Gene IDs remain in the original metadata.
    Flank overlap alone never joins separate exons. Each union keeps every
    original exon and its anchored interval in query-oriented coordinates.
    """
    by_gene = defaultdict(list)
    components = []
    for exon in exons:
        rec = exon.record
        gene_key = exon_query_gene_key(rec)
        if not gene_key[1] or not merge_overlapping:
            components.append([exon])
            continue
        key = (rec.chrom, rec.strand, gene_key)
        by_gene[key].append(exon)

    for key in sorted(by_gene):
        component = []
        end = -1
        for exon in sorted(by_gene[key], key=lambda x: (
            x.record.start0, x.record.end0, x.record.exon_id_full
        )):
            if component and exon.record.start0 >= end:
                components.append(component)
                component = []
                end = -1
            component.append(exon)
            end = max(end, exon.record.end0)
        if component:
            components.append(component)

    queries = []
    for component in components:
        if len(component) == 1:
            target = component[0]
        else:
            first = component[0].record
            start = min(x.record.start0 for x in component)
            end = max(x.record.end0 for x in component)
            ordered = sorted(component, key=lambda x: genomic_query_span(x))
            query_start = genomic_query_span(ordered[0])[0]
            cursor = query_start
            pieces = []
            for exon in ordered:
                left, right = genomic_query_span(exon)
                sequence = revcomp(exon.sequence) if first.strand == "-" else exon.sequence
                if left > cursor:
                    raise AssertionError("overlapping exon union has a reference gap")
                if right > cursor:
                    pieces.append(sequence[cursor - left:])
                    cursor = right
            sequence = "".join(pieces)
            left_anchor, right_anchor = start - query_start, cursor - end
            if first.strand == "-":
                sequence = revcomp(sequence)
                left_anchor, right_anchor = right_anchor, left_anchor
            identity = json.dumps([first.chrom, first.strand, start, end,
                                   exon_query_gene_key(first)], separators=(",", ":"))
            query_id = "PTEXON_" + hashlib.sha256(identity.encode()).hexdigest()
            fields = {
                name: [value for exon in component for value in getattr(exon.record, name)]
                for name in ("transcript_ids_full", "transcript_ids", "transcript_indices",
                             "gene_ids_full", "gene_ids", "gene_names", "transcript_types")
            }
            record = replace(first, start0=start, end0=end,
                             exon_id_full=query_id, exon_id=query_id,
                             any_mane=any(x.record.any_mane for x in component), **fields)
            target = ExtractedExon(record, sequence, unmasked_acgt_count(sequence),
                                   left_anchor, right_anchor)

        query_start, query_end = genomic_query_span(target)
        aliases = []
        for exon in component:
            rec = exon.record
            anchor_start, anchor_end = genomic_query_span(exon)
            if rec.strand == "-":
                core = (query_end - rec.end0, query_end - rec.start0)
                anchored = (query_end - anchor_end, query_end - anchor_start)
            else:
                core = (rec.start0 - query_start, rec.end0 - query_start)
                anchored = (anchor_start - query_start, anchor_end - query_start)
            if target.sequence[anchored[0]:anchored[1]] != exon.sequence:
                raise AssertionError("exon alias does not match its union query")
            aliases.append(QueryExonAlias(exon, *core, *anchored, len(component) > 1))
        queries.append(ExonQuery(target, aliases))
    return queries


def reciprocal_interval_overlap(a: ExonRecord, b: ExonRecord) -> float:
    overlap = max(0, min(a.end0, b.end0) - max(a.start0, b.start0))
    if overlap == 0:
        return 0.0
    return min(overlap / (a.end0 - a.start0), overlap / (b.end0 - b.start0))


def group_redundant_exons(
    exons: Sequence[ExtractedExon],
) -> List[List[ExtractedExon]]:
    """Group only anchored sequences with identical core-exon boundaries."""
    parent = list(range(len(exons)))
    rank = [0] * len(exons)

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        if rank[ra] == rank[rb]:
            rank[ra] += 1

    # The anchor context is part of the query.  Core-identical exons from
    # different loci are not interchangeable when their flanks differ.
    exact_representative: Dict[Tuple[str, int, int, int], int] = {}
    for i, exon in enumerate(exons):
        key = (
            exon.sequence,
            exon.left_anchor_length,
            exon.right_anchor_length,
            exon.exon_length,
        )
        prior = exact_representative.get(key)
        if prior is None:
            exact_representative[key] = i
        else:
            union(prior, i)

    grouped: Dict[int, List[ExtractedExon]] = defaultdict(list)
    for i, exon in enumerate(exons):
        grouped[find(i)].append(exon)

    groups = list(grouped.values())
    for group in groups:
        group.sort(
            key=lambda x: (
                -len(x.sequence),
                -x.unmasked,
                x.record.exon_id_full,
                x.record.coordinates,
            )
        )
    groups.sort(
        key=lambda g: (
            g[0].record.chrom,
            g[0].record.start0,
            g[0].record.end0,
            g[0].record.strand,
            g[0].record.exon_id_full,
        )
    )
    return groups


def parse_gencode_gff3(
    gff3_path: str,
    mane_by_gene: Optional[Dict[str, set[str]]] = None,
    gene_order: Optional[List[str]] = None,
) -> List[RawExonRow]:
    """Read exon rows from GENCODE-style GFF3, skipping comment lines."""
    mane_transcripts_full = set()
    mane_transcripts = set()
    seen_genes = set()

    with open_text(gff3_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            gid = strip_version(attrs.get("gene_id") or
                                (attrs.get("ID", "") if parts[2] == "gene" else attrs.get("Parent", "")))
            if gene_order is not None and gid and parts[2] in {"gene", "transcript"} and gid not in seen_genes:
                gene_order.append(gid)
                seen_genes.add(gid)
            if parts[2] != "transcript":
                continue
            transcript_id_full = attrs.get("transcript_id") or attrs.get("ID") or ""
            if transcript_id_full and attr_has_mane(attrs):
                mane_transcripts_full.add(transcript_id_full)
                mane_transcripts.add(strip_version(transcript_id_full))
                if mane_by_gene is not None and gid:
                    mane_by_gene.setdefault(gid, set()).add(strip_version(transcript_id_full))

    rows: List[RawExonRow] = []
    required_missing = defaultdict(int)
    with open_text(gff3_path) as handle:
        for line_num, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = parse_attrs(parts[8])
            exon_id_full = attrs.get("exon_id", "")
            transcript_id_full = attrs.get("transcript_id") or attrs.get("Parent", "")
            gene_id_full = attrs.get("gene_id", "")
            gene_name = attrs.get("gene_name", "")
            exon_number = attrs.get("exon_number", "")
            transcript_type = attrs.get("transcript_type") or attrs.get("gene_type", "")

            for key, value in [
                ("exon_id", exon_id_full),
                ("transcript_id", transcript_id_full),
                ("gene_id", gene_id_full),
                ("gene_name", gene_name),
                ("exon_number", exon_number),
                ("transcript_type", transcript_type),
            ]:
                if not value:
                    required_missing[key] += 1
            if not exon_id_full or not transcript_id_full:
                continue

            try:
                start0 = int(parts[3]) - 1
                end0 = int(parts[4])
            except ValueError:
                print(f"WARNING: skipping exon with non-integer coordinates at line {line_num}", file=sys.stderr)
                continue

            strand = parts[6]
            if strand not in {"+", "-"}:
                print(f"WARNING: skipping exon with unsupported strand {strand!r} at line {line_num}", file=sys.stderr)
                continue

            is_mane = (
                transcript_id_full in mane_transcripts_full
                or strip_version(transcript_id_full) in mane_transcripts
                or attr_has_mane(attrs)
            )
            gid = strip_version(gene_id_full)
            if gene_order is not None and gid and gid not in seen_genes:
                gene_order.append(gid)
                seen_genes.add(gid)
            if mane_by_gene is not None and is_mane and gid:
                mane_by_gene.setdefault(gid, set()).add(strip_version(transcript_id_full))
            rows.append(
                RawExonRow(
                    chrom=parts[0],
                    start0=start0,
                    end0=end0,
                    strand=strand,
                    exon_id_full=exon_id_full,
                    exon_id=strip_version(exon_id_full),
                    transcript_id_full=transcript_id_full,
                    transcript_id=strip_version(transcript_id_full),
                    exon_number=exon_number,
                    gene_id_full=gene_id_full,
                    gene_id=strip_version(gene_id_full),
                    gene_name=gene_name,
                    transcript_type=transcript_type,
                    is_mane=is_mane,
                )
            )

    if required_missing:
        msg = ", ".join(f"{k}:{v}" for k, v in sorted(required_missing.items()))
        print(f"WARNING: missing GFF3 attributes among exon rows: {msg}", file=sys.stderr)
    return rows


def collapse_exons(rows: Sequence[RawExonRow], per_transcript_records: bool = False) -> List[ExonRecord]:
    """Collapse shared exon_id rows unless per_transcript_records is requested."""
    records: Dict[Tuple[str, str, int, int, str, str], ExonRecord] = {}
    for row in rows:
        if per_transcript_records:
            key = (row.exon_id_full, row.transcript_id_full, row.start0, row.end0, row.strand, row.chrom)
        else:
            key = (row.exon_id_full, "", row.start0, row.end0, row.strand, row.chrom)
        rec = records.get(key)
        if rec is None:
            rec = ExonRecord(
                chrom=row.chrom,
                start0=row.start0,
                end0=row.end0,
                strand=row.strand,
                exon_id_full=row.exon_id_full,
                exon_id=row.exon_id,
            )
            records[key] = rec
        rec.transcript_ids_full.append(row.transcript_id_full)
        rec.transcript_ids.append(row.transcript_id)
        rec.transcript_indices.append(row.exon_number)
        rec.gene_ids_full.append(row.gene_id_full)
        rec.gene_ids.append(row.gene_id)
        rec.gene_names.append(row.gene_name)
        rec.transcript_types.append(row.transcript_type)
        rec.any_mane = rec.any_mane or row.is_mane

    out = list(records.values())
    out.sort(key=lambda r: (r.chrom, r.start0, r.end0, r.strand, r.exon_id_full, r.transcript_id_field))
    return out


def fasta_records(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def revcomp(seq: str) -> str:
    return seq.translate(RC_TABLE)[::-1]


def merge_identical_mane_genes(
    rows: Sequence[RawExonRow], genome_fasta: str,
    mane_by_gene: Dict[str, set[str]], gene_order: Sequence[str],
) -> Tuple[List[RawExonRow], List[Dict[str, str]]]:
    """Compare full spliced MANE DNA before anchoring or masking filters.

    Exact sequence sets must match for genes with multiple MANE transcripts.
    Case is ignored; incomplete or ambiguous DNA cannot prove identity. Keep
    every isoform of the first gene in GFF order and discard all other members.
    """
    by_transcript = defaultdict(list)
    names = {}
    for row in rows:
        names.setdefault(row.gene_id, row.gene_name or row.gene_id)
        if row.transcript_id in mane_by_gene.get(row.gene_id, set()):
            by_transcript[(row.gene_id, row.transcript_id)].append(row)
    by_chrom = defaultdict(set)
    for exons in by_transcript.values():
        for exon in exons:
            by_chrom[exon.chrom].add((exon.start0, exon.end0, exon.strand))
    core_sequences = {}
    if by_chrom:
        for chrom, sequence in fasta_records(genome_fasta):
            for start, end, strand in by_chrom.get(chrom, ()):
                if 0 <= start < end <= len(sequence):
                    core = sequence[start:end].upper()
                    core_sequences[(chrom, start, end, strand)] = revcomp(core) if strand == "-" else core
    transcript_sequences = {}
    invalid = 0
    for gene, tids in mane_by_gene.items():
        for tid in tids:
            exons = by_transcript.get((gene, tid), [])
            by_number = {}
            valid = bool(exons)
            for exon in exons:
                try:
                    number = int(exon.exon_number)
                except ValueError:
                    valid = False
                    break
                key = (exon.chrom, exon.start0, exon.end0, exon.strand)
                if number in by_number and by_number[number] != key:
                    valid = False
                by_number[number] = key
            numbers = sorted(by_number)
            if numbers != list(range(1, len(numbers) + 1)):
                valid = False
            keys = [by_number[n] for n in numbers]
            if len({(key[0], key[3]) for key in keys}) != 1 or any(key not in core_sequences for key in keys):
                valid = False
            for left, right in zip(keys, keys[1:]):
                if (left[3] == "+" and left[2] > right[1]) or (left[3] == "-" and right[2] > left[1]):
                    valid = False
            sequence = "".join(core_sequences[key] for key in keys) if valid else ""
            if not sequence or set(sequence) - set("ACGT"):
                invalid += 1
                continue
            transcript_sequences[(gene, tid)] = sequence
    by_sequence_set = {}
    for gene in gene_order:
        tids = mane_by_gene.get(gene, set())
        if not tids or any((gene, tid) not in transcript_sequences for tid in tids):
            continue
        # Actual strings, not hashes, determine equality.
        sequence_set = tuple(sorted({transcript_sequences[(gene, tid)] for tid in tids}))
        by_sequence_set.setdefault(sequence_set, []).append(gene)
    reports = []
    removed = set()
    renamed = {}
    for genes in by_sequence_set.values():
        if len(genes) < 2:
            continue
        first = genes[0]
        merged_name = names[first] + "merged"
        renamed[first] = merged_name
        removed.update(genes[1:])
        tids_by_gene = [sorted(mane_by_gene[gene]) for gene in genes]
        reports.append({
            "representative_gene_id": first,
            "representative_gene_name": names[first],
            "merged_gene_name": merged_name,
            "gene_ids": ";".join(genes),
            "gene_names": ";".join(names[g] for g in genes),
            "mane_transcript_ids": ";".join(",".join(tids) for tids in tids_by_gene),
            "mane_sequence_lengths": ";".join(",".join(str(len(transcript_sequences[(g, t)])) for t in tids)
                                                   for g, tids in zip(genes, tids_by_gene)),
            "mane_sequence_sha256": ";".join(",".join(hashlib.sha256(transcript_sequences[(g, t)].encode()).hexdigest()
                                                       for t in tids) for g, tids in zip(genes, tids_by_gene)),
        })
    if invalid:
        print(f"Excluded {invalid} incomplete/ambiguous MANE transcripts from exact-paralog merging", file=sys.stderr)
    retained = [replace(row, gene_name=renamed[row.gene_id]) if row.gene_id in renamed else row
                for row in rows if row.gene_id not in removed]
    print(f"Merged {len(removed)} identical MANE paralog genes into {len(reports)} representative groups", file=sys.stderr)
    return retained, reports


def unmasked_acgt_count(seq: str) -> int:
    return sum(1 for base in seq if base in "ACGT")


def wrap_fasta(seq: str, width: int = 80) -> Iterator[str]:
    for i in range(0, len(seq), width):
        yield seq[i : i + width]


def write_exon_fasta(
    genome_fasta: str,
    exon_records: Sequence[ExonRecord],
    out_fasta: str,
    out_seq: str,
    out_info: str,
    out_aliases: str,
    anchor_target_length: int,
    min_unmasked: int,
    merge_overlapping: bool = True,
) -> Tuple[int, int, int, int]:
    by_chrom: Dict[str, List[ExonRecord]] = defaultdict(list)
    for rec in exon_records:
        by_chrom[rec.chrom].append(rec)
    gene_bounds = {}
    for rec in exon_records:
        key = (exon_query_gene_key(rec), rec.chrom, rec.strand)
        start, end = gene_bounds.get(key, (rec.start0, rec.end0))
        gene_bounds[key] = (min(start, rec.start0), max(end, rec.end0))
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda r: (r.start0, r.end0, r.strand, r.exon_id_full))

    written = 0
    skipped_masked = 0
    skipped_out_of_range = 0
    chroms_seen = set()
    extracted: List[ExtractedExon] = []

    for chrom, chrom_seq in fasta_records(genome_fasta):
        if chrom not in by_chrom:
            continue
        chroms_seen.add(chrom)
        chrom_len = len(chrom_seq)
        for rec in by_chrom[chrom]:
            if rec.start0 < 0 or rec.end0 > chrom_len or rec.start0 >= rec.end0:
                skipped_out_of_range += 1
                print(
                    f"WARNING: skipping out-of-range exon {rec.exon_id_full} {rec.coordinates} for {chrom} length {chrom_len}",
                    file=sys.stderr,
                )
                continue
            anchor_size = max(0, (anchor_target_length - (rec.end0 - rec.start0) + 1) // 2)
            fragment_start = max(0, rec.start0 - anchor_size)
            fragment_end = min(chrom_len, rec.end0 + anchor_size)
            genomic_left_anchor = rec.start0 - fragment_start
            genomic_right_anchor = fragment_end - rec.end0
            seq = chrom_seq[fragment_start:fragment_end]
            if rec.strand == "-":
                # Report anchors in emitted/transcript FASTA orientation.
                left_anchor_length = genomic_right_anchor
                right_anchor_length = genomic_left_anchor
                seq = revcomp(seq)
            else:
                left_anchor_length = genomic_left_anchor
                right_anchor_length = genomic_right_anchor
            unmasked = unmasked_acgt_count(seq)
            if unmasked < min_unmasked:
                skipped_masked += 1
                continue
            extracted.append(
                ExtractedExon(
                    record=rec,
                    sequence=seq,
                    unmasked=unmasked,
                    left_anchor_length=left_anchor_length,
                    right_anchor_length=right_anchor_length,
                )
            )

    for chrom in sorted(set(by_chrom) - chroms_seen):
        print(f"WARNING: contig {chrom} appears in GFF3 but not in genome FASTA", file=sys.stderr)

    queries = merge_overlapping_exon_queries(extracted, merge_overlapping)
    query_aliases = {id(query.target): query.aliases for query in queries}
    groups = group_redundant_exons([query.target for query in queries])
    alias_count = sum(len(query.aliases) for query in queries)
    print(f"Merged {len(extracted)} retained exon records into {len(queries)} "
          "same-gene overlap queries before exact-sequence deduplication", file=sys.stderr)

    with open(out_fasta, "w", encoding="utf-8") as fa_out, \
        open(out_seq, "w", encoding="utf-8") as seq_out, \
        open(out_info, "w", encoding="utf-8") as info_out, \
        open(out_aliases, "w", encoding="utf-8") as alias_out:
        info_out.write(
            "exon_id_full\texon_id\tcoordinates\tchrom\tstart0\tend0\tstrand\t"
            "transcript_id_full\ttranscript_id\ttranscript_index\tifmane_transcript\t"
            "gene_id_full\tgene_id\tgene_name\tifproteincoding\tlength\t"
            "left_anchor_length\tright_anchor_length\tanchored_length\tunmasked_ACGT\n"
        )
        alias_out.write(
            "blast_exon_id_full\tblast_exon_id\texon_id_full\texon_id\tcoordinates\t"
            "chrom\tstart0\tend0\tstrand\tlength\tgene_id_full\tgene_id\tgene_name\t"
            "transcript_id_full\ttranscript_id\ttranscript_index\tifmane_transcript\t"
            "ifproteincoding\tleft_anchor_length\tright_anchor_length\t"
            "anchored_length\tmerge_reason\tquery_core_start0\tquery_core_end0\t"
            "query_anchor_start0\tquery_anchor_end0\tgene_start0\tgene_end0\n"
        )
        for group in groups:
            representative = group[0]
            rec = representative.record
            seq = representative.sequence
            header_fields = rec.header_fields(
                representative.left_anchor_length,
                representative.right_anchor_length,
            )
            header = "\t".join(header_fields)
            fa_out.write(f">{header}\n")
            for chunk in wrap_fasta(seq):
                fa_out.write(chunk + "\n")
            seq_out.write(header_fields[0].split()[0] + "\n")
            info_out.write(
                "\t".join(
                    [
                        rec.exon_id_full, rec.exon_id, rec.coordinates, rec.chrom,
                        str(rec.start0), str(rec.end0), rec.strand,
                        rec.transcript_id_field,
                        ",".join(uniq_preserve_order(rec.transcript_ids)),
                        rec.transcript_index_field, "1" if rec.any_mane else "0",
                        rec.gene_id_field, ",".join(uniq_preserve_order(rec.gene_ids)),
                        rec.gene_name_field, rec.ifproteincoding,
                        str(representative.exon_length),
                        str(representative.left_anchor_length),
                        str(representative.right_anchor_length), str(len(seq)),
                        str(representative.unmasked),
                    ]
                ) + "\n"
            )
            seen_aliases = set()
            for projection in (alias for target in group for alias in query_aliases[id(target)]):
                alias = projection.exon
                arec = alias.record
                alias_key = (arec.exon_id_full, arec.coordinates)
                if alias_key in seen_aliases:
                    continue
                seen_aliases.add(alias_key)
                if projection.overlap_merged:
                    reason = "same_gene_overlap"
                elif alias is representative:
                    reason = "representative"
                elif alias.sequence == representative.sequence:
                    reason = "exact_sequence"
                else:
                    raise AssertionError("non-identical anchored exons were grouped")
                alias_out.write(
                    "\t".join(
                        [
                            rec.exon_id_full, rec.exon_id,
                            arec.exon_id_full, arec.exon_id, arec.coordinates,
                            arec.chrom, str(arec.start0), str(arec.end0), arec.strand,
                            str(alias.exon_length), arec.gene_id_field,
                            ",".join(uniq_preserve_order(arec.gene_ids)),
                            arec.gene_name_field, arec.transcript_id_field,
                            ",".join(uniq_preserve_order(arec.transcript_ids)),
                            arec.transcript_index_field, "1" if arec.any_mane else "0",
                            arec.ifproteincoding,
                            str(alias.left_anchor_length),
                            str(alias.right_anchor_length), str(len(alias.sequence)),
                            reason, str(projection.core_start0), str(projection.core_end0),
                            str(projection.anchor_start0), str(projection.anchor_end0),
                            *(str(v) for v in gene_bounds[(exon_query_gene_key(arec), arec.chrom, arec.strand)]),
                        ]
                    ) + "\n"
                )
            written += 1
    return written, skipped_masked, skipped_out_of_range, alias_count


def run_makeblastdb(makeblastdb: str, exon_fasta: str, out_prefix: str) -> None:
    cmd = [
        makeblastdb,
        "-in",
        exon_fasta,
        "-dbtype",
        "nucl",
        "-blastdb_version",
        "5",
        "-out",
        out_prefix,
    ]
    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a BLASTN exon database from reference genome FASTA and GENCODE-style GFF3."
    )
    parser.add_argument("-f", "--genome", required=True, help="reference genome FASTA, optionally gzipped")
    parser.add_argument("-g", "--gff3", required=True, help="GENCODE-style GFF3 annotation, optionally gzipped")
    parser.add_argument("-o", "--out", required=True, help="output BLAST database prefix")
    parser.add_argument("--exon-fasta", default=None, help="output exon FASTA path [default: <out>.exons.fa]")
    parser.add_argument("--anchor-target-length", type=int, default=150, help="pad short exons to this minimum length, symmetrically; longer exons have no anchors [150]")
    parser.add_argument("--min-unmasked", type=int, default=50, help="minimum uppercase A/C/G/T bases in the full anchored sequence [50]")
    parser.add_argument("--makeblastdb", default="makeblastdb", help="path to makeblastdb [makeblastdb]")
    parser.add_argument("--no-makeblastdb", action="store_true", help="write FASTA/metadata only; do not run makeblastdb")
    parser.add_argument(
        "--merge-exon-overlap",
        type=float,
        default=99.0,
        help="legacy compatibility option (ignored); overlapping exons of the same gene are always unioned",
    )
    parser.add_argument(
        "--per-transcript-records",
        action="store_true",
        help="write one FASTA record per exon-transcript GFF3 row instead of collapsing shared exon IDs",
    )
    args = parser.parse_args()
    if args.anchor_target_length < 0:
        raise SystemExit("ERROR: --anchor-target-length cannot be negative")
    if args.min_unmasked < 0:
        raise SystemExit("ERROR: --min-unmasked cannot be negative")
    if not 0.0 < args.merge_exon_overlap <= 100.0:
        raise SystemExit("ERROR: --merge-exon-overlap must be greater than 0 and at most 100")

    exon_fasta = args.exon_fasta or f"{args.out}.exons.fa"
    seq_file = f"{args.out}.seq"
    info_file = f"{args.out}.exon_info.tsv"
    alias_file = f"{args.out}.exon_aliases.tsv"
    paralog_file = os.path.join(os.path.dirname(os.path.abspath(args.out)), REPORT_NAME)
    shared_gene_file = os.path.join(os.path.dirname(os.path.abspath(args.out)), shared_exon_genes.REPORT_NAME)

    # Output prefixes commonly include a new directory (for example,
    # exon_database/reference_exons).  Create every parent before opening the
    # FASTA/metadata files or invoking makeblastdb.
    for output_path in [exon_fasta, seq_file, info_file, alias_file, args.out]:
        parent = os.path.dirname(os.path.abspath(output_path))
        os.makedirs(parent, exist_ok=True)

    mane_by_gene: Dict[str, set[str]] = {}
    gene_order: List[str] = []
    raw_rows = parse_gencode_gff3(args.gff3, mane_by_gene, gene_order)
    if not raw_rows:
        raise SystemExit("ERROR: no exon records found. Check that the GFF3 is GENCODE-style and has feature type 'exon'.")
    # Shared reference genes must retain every member and isoform, even when
    # two members happen to have identical MANE sequences. Only other genes
    # participate in the older representative-only paralog collapse.
    provisional_rows, shared_groups = shared_exon_genes.merge_shared_exon_rows(raw_rows)
    del provisional_rows
    shared_members = {gene for group in shared_groups
                      for gene in group["member_gene_ids"].split(";")}
    mane_for_paralog_merge = {gene: tids for gene, tids in mane_by_gene.items()
                             if gene not in shared_members}
    raw_rows, paralog_groups = merge_identical_mane_genes(raw_rows, args.genome, mane_for_paralog_merge, gene_order)
    write_report(paralog_file, paralog_groups)
    raw_rows, shared_groups = shared_exon_genes.merge_shared_exon_rows(raw_rows)
    shared_exon_genes.write_report(shared_gene_file, shared_groups)
    print(f"Merged reference genes into {len(shared_groups)} shared-exon gene units", file=sys.stderr)
    exon_records = collapse_exons(raw_rows, per_transcript_records=args.per_transcript_records)

    print(f"Read {len(raw_rows)} GFF3 exon rows", file=sys.stderr)
    print(f"Prepared {len(exon_records)} exon FASTA records", file=sys.stderr)

    written, skipped_masked, skipped_oor, alias_count = write_exon_fasta(
        genome_fasta=args.genome,
        exon_records=exon_records,
        out_fasta=exon_fasta,
        out_seq=seq_file,
        out_info=info_file,
        out_aliases=alias_file,
        anchor_target_length=args.anchor_target_length,
        min_unmasked=args.min_unmasked,
    )
    print(f"Wrote {written} exon FASTA records to {exon_fasta}", file=sys.stderr)
    print(f"Padded short exons toward {args.anchor_target_length} bp with symmetric, contig-clipped anchors", file=sys.stderr)
    print(f"Skipped {skipped_masked} anchored exon records with uppercase A/C/G/T < {args.min_unmasked}", file=sys.stderr)
    print(f"Skipped {skipped_oor} out-of-range exon records", file=sys.stderr)
    print(f"Wrote BLAST ordinal name map to {seq_file}", file=sys.stderr)
    print(f"Wrote exon metadata to {info_file}", file=sys.stderr)
    print(f"Wrote {alias_count} original exon aliases to {alias_file}", file=sys.stderr)
    print(f"Collapsed {max(0, alias_count - written)} redundant exon targets", file=sys.stderr)

    if written == 0:
        raise SystemExit("ERROR: no exon sequences passed the filters; no BLAST database was built.")
    manifest = {
        "format": DATABASE_FORMAT,
        "reference": os.path.realpath(args.genome), "gff3": os.path.realpath(args.gff3),
        "reference_identity": file_identity(args.genome), "gff3_identity": file_identity(args.gff3),
        "anchor_target_length": args.anchor_target_length, "min_unmasked": args.min_unmasked,
        "merge_exon_overlap": args.merge_exon_overlap,
        "exon_query_merge_policy": EXON_QUERY_MERGE_POLICY,
        "identical_mane_policy": MERGE_POLICY,
        "identical_paralogs_sha256": report_digest(paralog_file),
        "shared_exon_gene_policy": shared_exon_genes.MERGE_POLICY,
        "shared_exon_genes_sha256": report_digest(shared_gene_file),
    }
    with open(f"{args.out}.manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Wrote identical MANE paralogs to {paralog_file}", file=sys.stderr)
    print(f"Wrote shared-exon gene units to {shared_gene_file}", file=sys.stderr)
    if not args.no_makeblastdb:
        run_makeblastdb(args.makeblastdb, exon_fasta, args.out)


if __name__ == "__main__":
    main()
