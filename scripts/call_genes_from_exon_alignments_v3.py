#!/usr/bin/env python3
"""
call_genes_from_exon_alignments_v3.py

Exon/transcript caller retaining gene ties after MANE priority and isoform selection.

It reads either:
  1. the extended, headered TSV from align_exon_blastdb.polished.py, or
  2. the original/legacy 13-column headerless TSV.

It parses a GENCODE-style GFF3, maps exon hits to transcripts, chains exon hits
in transcript exon order, filters consistently by exon coverage and identity,
and writes transcript/gene calls.

Coordinate convention:
  - GFF3 input is converted from 1-based closed to 0-based right-open.
  - --region and --truncate are interpreted as 0-based right-open.
  - Query FASTA names are expected to end with _start_end when absolute genomic
    output coordinates are desired, e.g. chr1_100000_120000.
  - Extended output stores per-exon query coordinates as start-end and reference
    coordinates as chrom:start-end; comma-separated entries correspond by index.
"""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
import csv
from functools import cached_property
import gzip
import hashlib
import heapq
from itertools import chain, groupby, islice
import json
import math
import multiprocessing as mp
import os
import pickle
import shutil
import sys
import tempfile
import time
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from identical_paralogs import infer_report, read_report
from shared_exon_genes import merge_shared_exon_rows, write_report as write_shared_exon_report

PIPELINE_VERSION = "3.9.0"
GENE_INSERTION_COST = 50.0
MAX_GENE_INSERTION_UNIQUE_EXONS = 20


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_attrs(attr_text: str) -> Dict[str, str]:
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


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return interval_overlap(a_start, a_end, b_start, b_end) > 0


class NonOverlappingIntervals:
    """Blocked sorted index for a set whose members never overlap each other."""

    BLOCK_SIZE = 512

    def __init__(self) -> None:
        self.blocks: List[List[Tuple[int, int]]] = []
        self.first_starts: List[int] = []

    def _block_index(self, start: int) -> int:
        if not self.blocks:
            return 0
        return max(0, bisect_right(self.first_starts, start) - 1)

    def _position(self, start: int) -> Tuple[int, int]:
        block_i = self._block_index(start)
        if not self.blocks:
            return block_i, 0
        item_i = bisect_left(self.blocks[block_i], (start, -1))
        return block_i, item_i

    def overlaps(self, start: int, end: int) -> bool:
        if not self.blocks:
            return False
        block_i, item_i = self._position(start)
        block = self.blocks[block_i]

        if item_i > 0:
            previous = block[item_i - 1]
        elif block_i > 0:
            previous = self.blocks[block_i - 1][-1]
        else:
            previous = None
        if previous is not None and previous[1] > start:
            return True

        if item_i < len(block):
            following = block[item_i]
        elif block_i + 1 < len(self.blocks):
            following = self.blocks[block_i + 1][0]
        else:
            following = None
        if following is not None and following[0] < end:
            return True
        return False

    def add(self, start: int, end: int) -> None:
        if not self.blocks:
            self.blocks.append([(start, end)])
            self.first_starts.append(start)
            return
        block_i, item_i = self._position(start)
        block = self.blocks[block_i]
        block.insert(item_i, (start, end))
        self.first_starts[block_i] = block[0][0]
        if len(block) > self.BLOCK_SIZE * 2:
            right = block[self.BLOCK_SIZE:]
            del block[self.BLOCK_SIZE:]
            self.blocks.insert(block_i + 1, right)
            self.first_starts.insert(block_i + 1, right[0][0])


@dataclass
class Region:
    chrom: str
    start: int
    end: int


def parse_region(text: str) -> Region:
    text = text.replace("/rc", "")
    if ":" not in text or "-" not in text:
        raise argparse.ArgumentTypeError(f"region must look like chrom:start-end, got {text!r}")
    chrom, rest = text.split(":", 1)
    start_s, end_s = rest.split("-", 1)
    start, end = int(start_s.replace(",", "")), int(end_s.replace(",", ""))
    if end < start:
        start, end = end, start
    return Region(chrom=chrom, start=start, end=end)


@dataclass
class QueryLocation:
    contig: str
    start: int
    end: int
    parsed: bool


def parse_query_location(query_id: str, parse_suffix: bool = True) -> QueryLocation:
    if not parse_suffix:
        return QueryLocation(contig=query_id, start=0, end=0, parsed=False)
    parts = query_id.split("_")
    if len(parts) >= 3:
        try:
            start = int(parts[-2])
            end = int(parts[-1])
            if end < start:
                start, end = end, start
            return QueryLocation(contig="_".join(parts[:-2]), start=start, end=end, parsed=True)
        except ValueError:
            pass
    return QueryLocation(contig=query_id, start=0, end=0, parsed=False)


@dataclass
class GffExon:
    exon_id: str
    exon_id_full: str
    transcript_id: str
    transcript_id_full: str
    exon_number: int
    gene_id: str
    gene_id_full: str
    gene_name: str
    transcript_type: str
    chrom: str
    start0: int
    end0: int
    strand: str
    is_mane: bool

    @property
    def length(self) -> int:
        return self.end0 - self.start0


@dataclass
class TranscriptInfo:
    transcript_id: str
    transcript_id_full: str
    gene_id: str
    gene_id_full: str
    gene_name: str
    transcript_type: str
    is_mane: bool = False
    exon_length: int = 0
    model_type: str = "transcript"
    source_gene_id: str = ""
    source_gene_name: str = ""

    @property
    def reported_gene_id(self) -> str:
        return self.source_gene_id if self.model_type == "transcript" and self.source_gene_id else self.gene_id

    @property
    def reported_gene_name(self) -> str:
        return self.source_gene_name if self.model_type == "transcript" and self.source_gene_name else self.gene_name


@dataclass
class Annotation:
    transcripts: Dict[str, TranscriptInfo]
    exons_by_transcript: Dict[str, List[GffExon]]
    exon_to_transcripts: Dict[
        str, List[Tuple[str, int, str, int, int]]
    ]
    shared_exon_gene_groups: List[dict] = field(default_factory=list)
    shared_exon_genes_applied: bool = False


def apply_identical_paralogs(annotation: Annotation, report_path) -> None:
    """Use only representative-gene isoforms, including for shared exon IDs."""
    shared_groups = apply_shared_exon_genes(annotation)
    shared_members = {gene for group in shared_groups
                      for gene in group["member_gene_ids"].split(";")}
    groups = read_report(report_path)
    removed = set()
    renamed = {}
    for group in groups:
        genes = group["gene_ids"].split(";")
        # Older reports may collapse genes now connected by reference overlap.
        # Keep all members of such groups under the new gene-unit policy.
        if shared_members.intersection(genes):
            continue
        removed.update(genes[1:])
        renamed[genes[0]] = group["merged_gene_name"]
    annotation.transcripts = {
        tid: replace(info, gene_name=renamed.get(info.gene_id, info.gene_name))
        for tid, info in annotation.transcripts.items() if info.gene_id not in removed
    }
    annotation.exons_by_transcript = {
        tid: [replace(exon, gene_name=renamed.get(exon.gene_id, exon.gene_name)) for exon in exons]
        for tid, exons in annotation.exons_by_transcript.items() if tid in annotation.transcripts
    }
    annotation.exon_to_transcripts = {
        eid: retained for eid, associations in annotation.exon_to_transcripts.items()
        if (retained := [item for item in associations if item[0] in annotation.transcripts])
    }


def apply_shared_exon_genes(annotation: Annotation) -> List[dict]:
    """Use reference overlap components as the gene identity in both stages."""
    if annotation.shared_exon_genes_applied:
        return annotation.shared_exon_gene_groups
    rows, groups = merge_shared_exon_rows(
        exon for exons in annotation.exons_by_transcript.values() for exon in exons)
    by_transcript = defaultdict(list)
    for exon in rows:
        by_transcript[exon.transcript_id].append(exon)
    for tid, exons in by_transcript.items():
        identities = {(e.gene_id, e.gene_id_full, e.gene_name) for e in exons}
        if len(identities) != 1:
            raise ValueError(f"transcript {tid} spans inconsistent reference gene units")
        gene_id, gene_id_full, gene_name = identities.pop()
        info = annotation.transcripts[tid]
        source = {}
        if gene_id != info.gene_id:
            source = dict(source_gene_id=info.source_gene_id or info.gene_id,
                          source_gene_name=info.source_gene_name or info.gene_name or info.gene_id)
        annotation.transcripts[tid] = replace(info,
            gene_id=gene_id, gene_id_full=gene_id_full, gene_name=gene_name, **source)
    annotation.exons_by_transcript = dict(by_transcript)
    annotation.shared_exon_gene_groups = groups
    annotation.shared_exon_genes_applied = True
    return groups


def add_full_gene_transcripts(annotation: Annotation) -> List[dict]:
    """Add exon-union candidates for shared-exon gene units, keeping all isoforms.

    Union reference exons within a gene ID, contig and strand. Original exon
    IDs still map the saved alignments to the synthetic block numbers. Keep
    original exon coordinates in those associations so alignment evidence is
    not presented as an alignment of the entire union block.
    """
    apply_shared_exon_genes(annotation)
    by_locus: Dict[Tuple[str, str, str], Dict[tuple, GffExon]] = defaultdict(dict)
    sources: Dict[Tuple[str, str, str], Dict[tuple, set[str]]] = defaultdict(lambda: defaultdict(set))
    for tid, exons in annotation.exons_by_transcript.items():
        if annotation.transcripts[tid].model_type == "full_gene":
            continue
        for exon in exons:
            if not exon.gene_id:
                continue
            locus = (exon.gene_id, exon.chrom, exon.strand)
            key = (exon.start0, exon.end0, exon.exon_id)
            by_locus[locus].setdefault(key, exon)
            sources[locus][key].add(tid)
    loci_by_gene: Dict[str, List[tuple]] = defaultdict(list)
    for locus in sorted(by_locus):
        loci_by_gene[locus[0]].append(locus)
    report = []
    for locus, unique in sorted(by_locus.items()):
        gene_id, chrom, strand = locus
        tid = gene_id + "_full_gene"
        if len(loci_by_gene[gene_id]) > 1:
            tid += "_" + str(loci_by_gene[gene_id].index(locus) + 1)
        if tid in annotation.transcripts:
            raise ValueError(f"full-gene transcript ID already exists: {tid}")
        ordered = sorted(unique.values(), key=lambda e: (e.start0, e.end0, e.exon_id))
        blocks: List[List[GffExon]] = []
        block_end = -1
        for exon in ordered:
            if not blocks or exon.start0 >= block_end:
                blocks.append([])
                block_end = exon.end0
            else:
                block_end = max(block_end, exon.end0)
            blocks[-1].append(exon)
        if strand == "-":
            blocks.reverse()
        first = ordered[0]
        biotypes = {e.transcript_type for e in ordered}
        biotype = "protein_coding" if "protein_coding" in biotypes else sorted(biotypes)[0]
        annotation.transcripts[tid] = TranscriptInfo(
            tid, tid, gene_id, first.gene_id_full, first.gene_name, biotype,
            is_mane=False,
            exon_length=sum(max(e.end0 for e in block) - min(e.start0 for e in block) for block in blocks),
            model_type="full_gene",
        )
        model_exons = []
        for number, block in enumerate(blocks, 1):
            transcript_ids = set()
            for exon in block:
                key = (exon.start0, exon.end0, exon.exon_id)
                transcript_ids.update(sources[locus][key])
                model_exons.append(replace(exon, transcript_id=tid, transcript_id_full=tid,
                                           exon_number=number, transcript_type=biotype, is_mane=False))
                association = (tid, number, chrom, exon.start0, exon.end0)
                annotation.exon_to_transcripts.setdefault(exon.exon_id, []).append(association)
            report.append({
                "transcript_id": tid, "model_type": "full_gene", "gene_id": gene_id,
                "gene_name": first.gene_name, "reference_contig": chrom, "strand": strand,
                "exon_number": number, "start": min(e.start0 for e in block),
                "end": max(e.end0 for e in block),
                "source_exon_ids": ",".join(sorted({e.exon_id for e in block})),
                "source_transcript_ids": ",".join(sorted(transcript_ids)),
            })
        annotation.exons_by_transcript[tid] = model_exons
    return report


def write_full_gene_models(path: str, rows: Sequence[dict]) -> None:
    fields = ["transcript_id", "model_type", "gene_id", "gene_name", "reference_contig",
              "strand", "exon_number", "start", "end", "source_exon_ids", "source_transcript_ids"]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_gencode_gff3(gff3_path: str) -> Annotation:
    mane_full = set()
    mane_norm = set()
    transcripts: Dict[str, TranscriptInfo] = {}

    with open_text(gff3_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = parse_attrs(parts[8])
            tid_full = attrs.get("transcript_id") or attrs.get("ID") or ""
            if not tid_full:
                continue
            tid = strip_version(tid_full)
            is_mane = attr_has_mane(attrs)
            if is_mane:
                mane_full.add(tid_full)
                mane_norm.add(tid)
            gid_full = attrs.get("gene_id") or attrs.get("Parent") or ""
            transcript_type = attrs.get("transcript_type") or attrs.get("gene_type", "")
            transcripts[tid] = TranscriptInfo(
                transcript_id=tid,
                transcript_id_full=tid_full,
                gene_id=strip_version(gid_full),
                gene_id_full=gid_full,
                gene_name=attrs.get("gene_name", ""),
                transcript_type=transcript_type,
                is_mane=is_mane,
            )

    exons_by_transcript: Dict[str, List[GffExon]] = defaultdict(list)
    exon_to_transcripts: Dict[
        str, List[Tuple[str, int, str, int, int]]
    ] = defaultdict(list)
    missing_attr = defaultdict(int)

    with open_text(gff3_path) as handle:
        for line_num, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = parse_attrs(parts[8])
            exon_full = attrs.get("exon_id", "")
            tid_full = attrs.get("transcript_id") or attrs.get("Parent", "")
            gid_full = attrs.get("gene_id", "")
            exon_num_s = attrs.get("exon_number", "")
            gene_name = attrs.get("gene_name", "")
            transcript_type = attrs.get("transcript_type") or attrs.get("gene_type", "")
            for key, value in [
                ("exon_id", exon_full),
                ("transcript_id", tid_full),
                ("gene_id", gid_full),
                ("exon_number", exon_num_s),
                ("gene_name", gene_name),
                ("transcript_type", transcript_type),
            ]:
                if not value:
                    missing_attr[key] += 1
            if not exon_full or not tid_full or not exon_num_s:
                continue
            try:
                exon_num = int(exon_num_s)
                start0 = int(parts[3]) - 1
                end0 = int(parts[4])
            except ValueError:
                print(f"WARNING: skipping malformed exon at GFF3 line {line_num}", file=sys.stderr)
                continue
            strand = parts[6]
            if strand not in {"+", "-"}:
                continue
            tid = strip_version(tid_full)
            gid = strip_version(gid_full)
            exon_id = strip_version(exon_full)
            is_mane = tid_full in mane_full or tid in mane_norm or attr_has_mane(attrs)
            exon = GffExon(
                exon_id=exon_id,
                exon_id_full=exon_full,
                transcript_id=tid,
                transcript_id_full=tid_full,
                exon_number=exon_num,
                gene_id=gid,
                gene_id_full=gid_full,
                gene_name=gene_name,
                transcript_type=transcript_type,
                chrom=parts[0],
                start0=start0,
                end0=end0,
                strand=strand,
                is_mane=is_mane,
            )
            exons_by_transcript[tid].append(exon)
            exon_to_transcripts[exon_id].append(
                (tid, exon_num, exon.chrom, exon.start0, exon.end0)
            )
            if tid not in transcripts:
                transcripts[tid] = TranscriptInfo(
                    transcript_id=tid,
                    transcript_id_full=tid_full,
                    gene_id=gid,
                    gene_id_full=gid_full,
                    gene_name=gene_name,
                    transcript_type=transcript_type,
                    is_mane=is_mane,
                )
            else:
                transcripts[tid].is_mane = transcripts[tid].is_mane or is_mane
                if not transcripts[tid].transcript_type:
                    transcripts[tid].transcript_type = transcript_type
                if not transcripts[tid].gene_name:
                    transcripts[tid].gene_name = gene_name

    for tid in exons_by_transcript:
        exons_by_transcript[tid].sort(key=lambda e: e.exon_number)
        # Full spliced transcript length, including database-ineligible exons.
        # Duplicate annotation rows must not inflate the isoform tie breaker.
        transcripts[tid].exon_length = sum(
            end - start for _chrom, start, end in {
                (e.chrom, e.start0, e.end0) for e in exons_by_transcript[tid]
            }
        )
    for exon_id in exon_to_transcripts:
        # Deduplicate because GENCODE can contain duplicate feature rows for related features.
        seen = set()
        deduped = []
        for item in exon_to_transcripts[exon_id]:
            transcript_exon = item[:2]
            if transcript_exon in seen:
                continue
            seen.add(transcript_exon)
            deduped.append(item)
        exon_to_transcripts[exon_id] = deduped

    if missing_attr:
        msg = ", ".join(f"{k}:{v}" for k, v in sorted(missing_attr.items()))
        print(f"WARNING: missing GFF3 attributes among exon rows: {msg}", file=sys.stderr)
    if not exons_by_transcript:
        raise SystemExit("ERROR: no exon rows were parsed from the GFF3")
    return Annotation(
        transcripts=transcripts,
        exons_by_transcript=dict(exons_by_transcript),
        exon_to_transcripts=dict(exon_to_transcripts),
    )


@dataclass
class ExonAlignment:
    query_id: str
    percent_identity: float
    query_start: int
    query_end: int
    strand: str
    exon_id: str
    exon_length: int
    exon_start: int
    exon_end: int
    aligned_exon_bases: int
    exon_coverage: float
    identical_bases: int
    AS: float

    @cached_property
    def exon_score(self) -> float:
        """Original caller's mismatch-penalized score, normalized per exon.

        A matching base contributes +1 and a non-identical or unaligned base
        contributes -3 (the original expression was length - 4 * errors).
        Keeping the score normalized prevents a longer paralogous exon from
        winning solely because it is longer.
        """
        if self.exon_length <= 0:
            return 0.0
        identical = min(max(self.identical_bases, 0), self.exon_length)
        quality_bases = self.exon_length - 4 * (self.exon_length - identical)
        return 100.0 * quality_bases / self.exon_length

    @cached_property
    def mismatch_adjusted_bases(self) -> int:
        if self.exon_length <= 0:
            return 0
        identical = min(max(self.identical_bases, 0), self.exon_length)
        return self.exon_length - 4 * (self.exon_length - identical)

    @property
    def qspan(self) -> int:
        return self.query_end - self.query_start


def safe_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def safe_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except Exception:
        return default


def looks_like_header(parts: List[str]) -> bool:
    lowered = {p.lower() for p in parts}
    return "query_id" in lowered or "percent_identity" in lowered or "exon_coverage" in lowered


def read_alignments(path: str) -> List[ExonAlignment]:
    alignments: List[ExonAlignment] = []
    with open_text(path) as handle:
        first = handle.readline()
        if not first:
            return []
        # Normal pipeline output is tab-delimited.  Some historical tables,
        # including SMN_align_CHM13.tsv, were written with runs of spaces.
        whitespace_delimited = "\t" not in first

        def split_row(line: str) -> List[str]:
            if whitespace_delimited:
                return line.split()
            return line.rstrip("\n").split("\t")

        first_parts = split_row(first)
        header: Optional[List[str]] = first_parts if looks_like_header(first_parts) else None
        lines: Iterable[str]
        if header is None:
            lines = chain([first], handle)
        else:
            lines = handle
            header_index = {name: i for i, name in enumerate(header)}

        for line_num, line in enumerate(lines, start=2 if header else 1):
            if not line.strip():
                continue
            parts = split_row(line)
            try:
                if header is not None:
                    def get(name: str, default: str = "") -> str:
                        i = header_index.get(name)
                        if i is None or i >= len(parts):
                            return default
                        return parts[i]

                    query_id = get("query_id")
                    pi = safe_float(get("percent_identity"))
                    qstart = safe_int(get("query_start"))
                    qend = safe_int(get("query_end"))
                    strand = get("strand") or "+"
                    exon_id = strip_version(get("exon_id") or get("exon_id_full"))
                    exon_length = safe_int(get("exon_length"))
                    estart = safe_int(get("exon_start"))
                    eend = safe_int(get("exon_end"))
                    aligned = safe_int(get("aligned_exon_bases"), abs(eend - estart))
                    cov_text = get("exon_coverage")
                    cov = safe_float(cov_text, 100.0 * aligned / exon_length if exon_length else 0.0)
                    identical = safe_int(get("identical_bases"), int(round(aligned * pi / 100.0)))
                    AS = safe_float(get("AS"))
                elif (
                    len(parts) >= 14
                    and strip_version(parts[5]).startswith("ENSE")
                    and strip_version(parts[6]) == strip_version(parts[5])
                ):
                    # Headerless extended table from align_exon_blastdb_v2.py.
                    query_id = parts[0]
                    pi = safe_float(parts[1])
                    qstart = safe_int(parts[2])
                    qend = safe_int(parts[3])
                    strand = parts[4]
                    exon_id = strip_version(parts[5])
                    exon_length = safe_int(parts[7])
                    estart = safe_int(parts[8])
                    eend = safe_int(parts[9])
                    aligned = safe_int(parts[10], abs(eend - estart))
                    cov = safe_float(
                        parts[11],
                        100.0 * aligned / exon_length if exon_length else 0.0,
                    )
                    identical = safe_int(
                        parts[12], int(round(aligned * pi / 100.0))
                    )
                    AS = safe_float(parts[13])
                else:
                    # Legacy 13-col table from align_exon_blastdb.polished.py --output-format legacy
                    # or the user's original format.
                    if len(parts) < 10:
                        print(f"WARNING: skipping short alignment row {line_num}", file=sys.stderr)
                        continue
                    query_id = parts[0]
                    pi = safe_float(parts[1], 100.0)
                    qstart = safe_int(parts[2])
                    qend = safe_int(parts[3])
                    strand = parts[4]
                    exon_id = strip_version(parts[5])
                    exon_length = safe_int(parts[6])
                    estart = safe_int(parts[7])
                    eend = safe_int(parts[8])
                    aligned = safe_int(parts[9], abs(eend - estart))
                    cov = safe_float(parts[11], 100.0 * aligned / exon_length if exon_length else 0.0) if len(parts) > 11 else (100.0 * aligned / exon_length if exon_length else 0.0)
                    identical = int(round(aligned * pi / 100.0))
                    AS = safe_float(parts[12] if len(parts) > 12 else parts[9])
            except Exception as exc:
                print(f"WARNING: skipping malformed alignment row {line_num}: {exc}", file=sys.stderr)
                continue

            if qend < qstart:
                qstart, qend = qend, qstart
            if eend < estart:
                estart, eend = eend, estart
            if not query_id or not exon_id:
                continue
            alignments.append(
                ExonAlignment(
                    query_id=query_id,
                    percent_identity=pi,
                    query_start=qstart,
                    query_end=qend,
                    strand=strand if strand in {"+", "-"} else "+",
                    exon_id=exon_id,
                    exon_length=exon_length,
                    exon_start=estart,
                    exon_end=eend,
                    aligned_exon_bases=aligned,
                    exon_coverage=cov,
                    identical_bases=identical,
                    AS=AS,
                )
            )
    return alignments


def region_allows_alignment(
    aln: ExonAlignment,
    region: Optional[Region],
    parse_query_suffix: bool = True,
) -> bool:
    if region is None:
        return True
    qloc = parse_query_location(aln.query_id, parse_suffix=parse_query_suffix)
    same_contig = qloc.contig == region.chrom or region.chrom in aln.query_id
    if not same_contig:
        return False
    abs_start = qloc.start + aln.query_start
    abs_end = qloc.start + aln.query_end
    return intervals_overlap(abs_start, abs_end, region.start, region.end)


def filter_alignments(
    alignments: Sequence[ExonAlignment],
    annotation: Annotation,
    region: Optional[Region],
    parse_query_suffix: bool = True,
) -> List[ExonAlignment]:
    kept = []
    dropped_no_exon = 0
    dropped_region = 0
    for aln in alignments:
        if aln.exon_id not in annotation.exon_to_transcripts:
            dropped_no_exon += 1
            continue
        if not region_allows_alignment(aln, region, parse_query_suffix=parse_query_suffix):
            dropped_region += 1
            continue
        kept.append(aln)
    print(f"Kept {len(kept)} exon alignments after caller filters", file=sys.stderr)
    print(f"Dropped {dropped_no_exon} not present in GFF3 exon_id map", file=sys.stderr)
    print(f"Dropped {dropped_region} outside --region", file=sys.stderr)
    return kept


def dedup_same_exon_overlaps(alignments: Sequence[ExonAlignment]) -> List[ExonAlignment]:
    """For each query/strand/exon, keep best non-overlapping HSPs."""
    grouped: Dict[Tuple[str, str, str], List[ExonAlignment]] = defaultdict(list)
    for aln in alignments:
        grouped[(aln.query_id, aln.strand, aln.exon_id)].append(aln)
    kept: List[ExonAlignment] = []
    for group in grouped.values():
        chosen: List[ExonAlignment] = []
        occupied = NonOverlappingIntervals()
        for aln in sorted(group, key=lambda a: (a.exon_score, a.aligned_exon_bases, a.AS), reverse=True):
            if occupied.overlaps(aln.query_start, aln.query_end):
                continue
            chosen.append(aln)
            occupied.add(aln.query_start, aln.query_end)
        kept.extend(chosen)
    kept.sort(key=lambda a: (a.query_id, a.strand, a.query_start, a.query_end, a.exon_id))
    return kept


@dataclass(frozen=True)
class MergedIntervalGeneScore:
    start: int
    end: int
    gene_id: str
    score: float
    source_exon_id: str = ""
    source_exon_length: int = 0


@dataclass
class TranscriptHit:
    exon_number: int
    alignment: ExonAlignment
    reference_chrom: str
    reference_start: int
    reference_end: int
    interval_score: Optional[MergedIntervalGeneScore] = None
    is_insertion: bool = False

    @property
    def score(self) -> float:
        if self.interval_score is not None:
            return self.interval_score.score
        return self.alignment.exon_score

    @property
    def competition_interval(self) -> Tuple[int, int]:
        if self.interval_score is not None:
            return self.interval_score.start, self.interval_score.end
        return self.alignment.query_start, self.alignment.query_end


@dataclass
class TranscriptCall:
    transcript_id: str
    transcript_info: TranscriptInfo
    hits: List[TranscriptHit]
    expected_exon_numbers: List[int]
    protein_bonus: float
    complete_bonus: float
    mane_bonus: float
    tie_group_id: str = ""
    tie_count: int = 1
    gene_index: str = ""
    gene_sort_score: float = 0.0
    gene_start: int = 0
    gene_end: int = 0

    @cached_property
    def found_exon_numbers(self) -> List[int]:
        return [h.exon_number for h in self.hits]

    @cached_property
    def found_expected_count(self) -> int:
        expected = set(self.expected_exon_numbers)
        return len({h.exon_number for h in self.hits if h.exon_number in expected and not h.is_insertion})

    @cached_property
    def expected_count(self) -> int:
        return len(self.expected_exon_numbers)

    @cached_property
    def missing_count(self) -> int:
        return max(0, self.expected_count - self.found_expected_count)

    @cached_property
    def raw_score(self) -> float:
        inserted = self.insertion_intervals
        return sum(s.score for s in self.merged_interval_scores if (s.start, s.end) not in inserted) + sum(
            h.score for h in self.hits if h.interval_score is None and not h.is_insertion
        ) - self.insertion_penalty

    @cached_property
    def insertion_intervals(self) -> set[Tuple[int, int]]:
        return {h.competition_interval for h in self.hits if h.is_insertion}

    @cached_property
    def insertion_run_exons(self) -> List[Tuple[int, ...]]:
        """Unique full-gene exon numbers in each consecutive insertion run."""
        runs = []
        current: set[int] = set()
        for hit in sorted(self.hits, key=lambda h: (*oriented_hit_interval(h, self.strand), h.exon_number)):
            if hit.is_insertion:
                current.add(hit.exon_number)
            elif current:
                runs.append(tuple(sorted(current)))
                current = set()
        if current:
            runs.append(tuple(sorted(current)))
        return runs

    @cached_property
    def inserted_exon_count(self) -> int:
        return sum(len(run) for run in self.insertion_run_exons)

    @cached_property
    def valid_insertion_runs(self) -> bool:
        return all(len(run) <= MAX_GENE_INSERTION_UNIQUE_EXONS for run in self.insertion_run_exons)

    @cached_property
    def insertion_penalty(self) -> float:
        return GENE_INSERTION_COST * self.inserted_exon_count

    @cached_property
    def merged_interval_scores(self) -> List[MergedIntervalGeneScore]:
        by_interval = {
            (h.interval_score.start, h.interval_score.end): h.interval_score
            for h in self.hits if h.interval_score is not None
        }
        return [by_interval[key] for key in sorted(by_interval)]

    @cached_property
    def total_aligned_bases(self) -> int:
        return sum(h.alignment.aligned_exon_bases for h in self.hits)

    @cached_property
    def mismatch_adjusted_bases(self) -> float:
        return sum(h.alignment.mismatch_adjusted_bases for h in self.hits)

    @cached_property
    def weighted_score(self) -> float:
        score = self.raw_score
        if self.transcript_info.transcript_type == "protein_coding":
            score *= self.protein_bonus
        if self.expected_count > 0 and self.found_expected_count >= self.expected_count:
            score *= self.complete_bonus
        if self.transcript_info.is_mane:
            score *= self.mane_bonus
        return score

    @cached_property
    def query_id(self) -> str:
        return self.hits[0].alignment.query_id if self.hits else ""

    @cached_property
    def strand(self) -> str:
        return self.hits[0].alignment.strand if self.hits else "+"

    def recalculated_with_hits(self, hits: List[TranscriptHit]) -> "TranscriptCall":
        return replace(self, hits=hits, tie_group_id="", tie_count=1)


def load_eligible_exon_ids(exon_info_path: Optional[str]) -> Optional[set[str]]:
    """Return version-stripped exon IDs present in the built BLAST database."""
    if exon_info_path is None:
        return None
    eligible: set[str] = set()
    with open_text(exon_info_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            exon_id_i = header.index("exon_id")
        except ValueError as exc:
            raise SystemExit(
                f"ERROR: {exon_info_path} is missing required column 'exon_id'"
            ) from exc
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if exon_id_i < len(parts) and parts[exon_id_i]:
                eligible.add(strip_version(parts[exon_id_i]))
    if not eligible:
        raise SystemExit(f"ERROR: no eligible exon IDs found in {exon_info_path}")
    return eligible


def truncate_expected_exons(
    annotation: Annotation,
    truncate: Optional[Region],
    eligible_exon_ids: Optional[set[str]] = None,
) -> Dict[str, List[int]]:
    expected: Dict[str, List[int]] = {}
    for tid, exons in annotation.exons_by_transcript.items():
        usable_exons = [
            e for e in exons
            if eligible_exon_ids is None or e.exon_id in eligible_exon_ids
        ]
        if truncate is None:
            nums = [e.exon_number for e in usable_exons]
        else:
            nums = [
                e.exon_number
                for e in usable_exons
                if e.chrom == truncate.chrom and intervals_overlap(e.start0, e.end0, truncate.start, truncate.end)
            ]
        expected[tid] = sorted(set(nums))
    return expected


class _FenwickCountStates:
    """Fenwick tree storing the best ended state for each chain length."""

    def __init__(self, size: int, scores: List[float]) -> None:
        self.tree: List[Dict[int, int]] = [dict() for _ in range(size + 1)]
        self.scores = scores

    def update(self, position: int, chain_count: int, state_index: int) -> None:
        while position < len(self.tree):
            current = self.tree[position].get(chain_count)
            if (
                current is None
                or self.scores[state_index] > self.scores[current]
                or (
                    self.scores[state_index] == self.scores[current]
                    and state_index < current
                )
            ):
                self.tree[position][chain_count] = state_index
            position += position & -position

    def query(self, position: int) -> List[int]:
        # For a fixed chain length, a lower-scoring state can never improve a
        # later chain. Retain the highest score and earliest original position.
        best_by_count: Dict[int, int] = {}
        while position > 0:
            for chain_count, state_index in self.tree[position].items():
                current = best_by_count.get(chain_count)
                if (
                    current is None
                    or self.scores[state_index] > self.scores[current]
                    or (
                        self.scores[state_index] == self.scores[current]
                        and state_index < current
                    )
                ):
                    best_by_count[chain_count] = state_index
            position -= position & -position
        return sorted(best_by_count.values())


def _best_chain_for_exon_hits(
    hits: Sequence[TranscriptHit], strand: str
) -> List[TranscriptHit]:
    """Sparse interval-sweep weighted exon-chain DP."""
    if not hits:
        return []

    if strand == "-":
        ordered = sorted(
            hits,
            key=lambda h: (
                -h.alignment.query_end,
                -h.alignment.query_start,
                h.exon_number,
            ),
        )

        def oriented_interval(hit: TranscriptHit) -> Tuple[int, int]:
            return -hit.alignment.query_end, -hit.alignment.query_start
    else:
        ordered = sorted(
            hits,
            key=lambda h: (
                h.alignment.query_start,
                h.alignment.query_end,
                h.exon_number,
            ),
        )

        def oriented_interval(hit: TranscriptHit) -> Tuple[int, int]:
            return hit.alignment.query_start, hit.alignment.query_end

    exon_numbers = sorted({hit.exon_number for hit in ordered})
    exon_rank = {number: i + 1 for i, number in enumerate(exon_numbers)}
    n = len(ordered)
    dp = [0.0] * n
    count = [0] * n
    prev = [-1] * n
    pending_by_end: List[Tuple[int, int]] = []
    fenwick = _FenwickCountStates(len(exon_numbers), dp)

    for i, hit in enumerate(ordered):
        start, end = oriented_interval(hit)
        while pending_by_end and pending_by_end[0][0] <= start:
            _prior_end, prior_i = heapq.heappop(pending_by_end)
            fenwick.update(
                exon_rank[ordered[prior_i].exon_number], count[prior_i], prior_i
            )

        dp[i] = hit.score
        count[i] = 1
        # Query strictly smaller exon numbers. Candidate states are restored to
        # their original scan order so the original floating-point/tie rule is
        # retained as closely as possible.
        for prior_i in fenwick.query(exon_rank[hit.exon_number] - 1):
            cand = dp[prior_i] + hit.score
            cand_count = count[prior_i] + 1
            if cand > dp[i] or (
                math.isclose(cand, dp[i]) and cand_count > count[i]
            ):
                dp[i] = cand
                count[i] = cand_count
                prev[i] = prior_i
        heapq.heappush(pending_by_end, (end, i))

    best_i = max(
        range(n),
        key=lambda i: (dp[i], count[i], ordered[i].alignment.aligned_exon_bases),
    )
    chain: List[TranscriptHit] = []
    i = best_i
    while i != -1:
        chain.append(ordered[i])
        i = prev[i]
    chain.reverse()
    return chain


def _best_chain_for_merged_hits(
    hits: Sequence[TranscriptHit], strand: str
) -> List[TranscriptHit]:
    """Chain isoform exons using each merged interval's gene score once.

    The global DP index contains states from earlier merged intervals. A local
    index handles exon transitions within the current interval at zero extra
    score. This preserves exon order/completeness without rewarding isoforms
    for subdividing the same scored interval into more exons.
    """
    if not hits or all(hit.interval_score is None for hit in hits):
        return _best_chain_for_exon_hits(hits, strand)
    if any(hit.interval_score is None for hit in hits):
        raise ValueError("mixed merged-interval and unscored exon hits")

    def oriented_interval(hit: TranscriptHit) -> Tuple[int, int]:
        a = hit.alignment
        return (-a.query_end, -a.query_start) if strand == "-" else (a.query_start, a.query_end)

    ordered = sorted(hits, key=lambda h: (*oriented_interval(h), h.exon_number))
    exon_numbers = sorted({h.exon_number for h in ordered})
    exon_rank = {number: i + 1 for i, number in enumerate(exon_numbers)}
    dp = [0.0] * len(ordered)
    count = [0] * len(ordered)
    previous = [-1] * len(ordered)
    global_index = _FenwickCountStates(len(exon_numbers), dp)
    local_index = _FenwickCountStates(len(exon_numbers), dp)
    pending: List[Tuple[int, int]] = []
    last_interval = None
    interval_states: List[int] = []

    for i, hit in enumerate(ordered):
        interval = hit.competition_interval
        if interval != last_interval:
            for prior in interval_states:
                global_index.update(exon_rank[ordered[prior].exon_number], count[prior], prior)
            local_index = _FenwickCountStates(len(exon_numbers), dp)
            pending = []
            interval_states = []
            last_interval = interval
        start, end = oriented_interval(hit)
        while pending and pending[0][0] <= start:
            _end, prior = heapq.heappop(pending)
            local_index.update(exon_rank[ordered[prior].exon_number], count[prior], prior)

        dp[i], count[i] = hit.score, 1
        for index, increment in ((global_index, hit.score), (local_index, 0.0)):
            for prior in index.query(exon_rank[hit.exon_number] - 1):
                score = dp[prior] + increment
                found = count[prior] + 1
                if score > dp[i] or (math.isclose(score, dp[i]) and found > count[i]):
                    dp[i], count[i], previous[i] = score, found, prior
        heapq.heappush(pending, (end, i))
        interval_states.append(i)

    best = max(range(len(ordered)), key=lambda i: (dp[i], count[i], ordered[i].alignment.aligned_exon_bases))
    chain: List[TranscriptHit] = []
    while best != -1:
        chain.append(ordered[best])
        best = previous[best]
    return list(reversed(chain))


def oriented_hit_interval(hit: TranscriptHit, strand: str) -> Tuple[int, int]:
    aln = hit.alignment
    return (-aln.query_end, -aln.query_start) if strand == "-" else (aln.query_start, aln.query_end)


def exon_copy_windows(
    hits: Sequence[TranscriptHit], strand: str,
) -> Dict[int, Tuple[float, float]]:
    """Bound links by disjoint alternative occurrences of either endpoint exon.

    A link p -> h cannot skip another occurrence of p or h that fits wholly
    between them. Overlapping hits are alternatives at one site, not barriers.
    Bounds use transcript-oriented query coordinates and impose no intron-size
    limit. Build them once from all hits so extracting one copy cannot erase
    the evidence separating copies in subsequent DP passes.
    """
    by_exon: Dict[int, List[TranscriptHit]] = defaultdict(list)
    for hit in hits:
        by_exon[hit.exon_number].append(hit)
    windows = {}
    for group in by_exon.values():
        intervals = [oriented_hit_interval(hit, strand) for hit in group]
        by_end = sorted((end, start) for start, end in intervals)
        ends = [end for end, _start in by_end]
        previous_starts = []
        maximum = -math.inf
        for _end, start in by_end:
            maximum = max(maximum, start)
            previous_starts.append(maximum)
        by_start = sorted(intervals)
        starts = [start for start, _end in by_start]
        following_ends = [math.inf] * len(by_start)
        minimum = math.inf
        for i in range(len(by_start) - 1, -1, -1):
            minimum = min(minimum, by_start[i][1])
            following_ends[i] = minimum
        for hit, (start, end) in zip(group, intervals):
            before = bisect_right(ends, start) - 1
            after = bisect_left(starts, end)
            windows[id(hit)] = (
                previous_starts[before] if before >= 0 else -math.inf,
                following_ends[after] if after < len(following_ends) else math.inf,
            )
    return windows


def best_chain_for_hits(
    hits: Sequence[TranscriptHit], strand: str,
    copy_windows: Optional[Dict[int, Tuple[float, float]]] = None,
) -> List[TranscriptHit]:
    """Find the best exon chain without jumping over supported exon copies."""
    if not hits:
        return []
    scored = [hit.interval_score is not None for hit in hits]
    if any(scored) and not all(scored):
        raise ValueError("mixed merged-interval and unscored exon hits")
    if copy_windows is None:
        copy_windows = exon_copy_windows(hits, strand)
    if all(copy_windows[id(hit)] == (-math.inf, math.inf) for hit in hits):
        return _best_chain_for_merged_hits(hits, strand)

    ordered = sorted(hits, key=lambda h: (*oriented_hit_interval(h, strand), h.exon_number))
    coordinates = [oriented_hit_interval(hit, strand) for hit in ordered]
    dp = [0.0] * len(ordered)
    count = [0] * len(ordered)
    previous = [-1] * len(ordered)
    pending: List[Tuple[int, int]] = []
    expiry: List[Tuple[float, int]] = []
    active: set[int] = set()
    for i, hit in enumerate(ordered):
        start, end = coordinates[i]
        while pending and pending[0][0] <= start:
            _end, prior = heapq.heappop(pending)
            active.add(prior)
            heapq.heappush(expiry, (copy_windows[id(ordered[prior])][1], prior))
        while expiry and expiry[0][0] <= start:
            _limit, prior = heapq.heappop(expiry)
            active.discard(prior)
        lower_bound = copy_windows[id(hit)][0]
        dp[i], count[i] = hit.score, 1
        for prior in sorted(active):
            predecessor = ordered[prior]
            if predecessor.exon_number >= hit.exon_number or coordinates[prior][1] <= lower_bound:
                continue
            same_interval = hit.interval_score is not None and predecessor.competition_interval == hit.competition_interval
            score = dp[prior] + (0.0 if same_interval else hit.score)
            found = count[prior] + 1
            if score > dp[i] or (math.isclose(score, dp[i]) and found > count[i]):
                dp[i], count[i], previous[i] = score, found, prior
        heapq.heappush(pending, (end, i))

    best = max(range(len(ordered)), key=lambda i: (dp[i], count[i], ordered[i].alignment.aligned_exon_bases))
    chain = []
    while best != -1:
        chain.append(ordered[best])
        best = previous[best]
    return list(reversed(chain))


class _GenePredecessorIndex:
    """Best DP state in a query-block range below a reference exon number.

    A segment tree over query blocks contains sparse Fenwick prefix maxima
    over reference exon numbers. Each insertion cost applies to a contiguous
    range of predecessors; at most 21 such ranges survive the unique-exon cap.
    This avoids rescanning arbitrarily long runs of repeated exon hits.
    """
    def __init__(self, blocks: Sequence[int], exons: Sequence[int], rank_key) -> None:
        self.size = 1
        while self.size <= max(blocks):
            self.size *= 2
        self.rank_key = rank_key
        coordinates = [set() for _ in range(2 * self.size)]
        for block, exon in zip(blocks, exons):
            coordinates[self.size + block].add(exon)
        for node in range(self.size - 1, 0, -1):
            coordinates[node] = coordinates[2 * node] | coordinates[2 * node + 1]
        self.coordinates = [sorted(values) for values in coordinates]
        self.trees = [[-1] * (len(values) + 1) for values in self.coordinates]

    def better(self, left: int, right: int) -> int:
        if left == -1:
            return right
        if right == -1:
            return left
        return right if self.rank_key(right) > self.rank_key(left) else left

    def update(self, block: int, exon: int, state: int) -> None:
        node = self.size + block
        while node:
            position = bisect_left(self.coordinates[node], exon) + 1
            tree = self.trees[node]
            while position < len(tree):
                tree[position] = self.better(tree[position], state)
                position += position & -position
            node //= 2

    def query(self, first_block: int, last_block: int, before_exon: int) -> int:
        left, right = first_block + self.size, last_block + self.size + 1
        best = -1
        while left < right:
            nodes = []
            if left & 1:
                nodes.append(left)
                left += 1
            if right & 1:
                right -= 1
                nodes.append(right)
            for node in nodes:
                position = bisect_left(self.coordinates[node], before_exon)
                while position:
                    best = self.better(best, self.trees[node][position])
                    position -= position & -position
            left //= 2
            right //= 2
        return best


def best_gene_chain_for_hits(
    hits: Sequence[TranscriptHit], strand: str,
    all_hits: Optional[Sequence[TranscriptHit]] = None,
) -> List[TranscriptHit]:
    """Local gene alignment with -50 per distinct exon per insertion run.

    Matches earn gene interval scores. A consecutive run of skipped query
    blocks costs 50 per unique reference full-gene exon number, using each
    block's deterministic longest-exon representative. More than 20 unique
    exons disallows that transition. A match ends the run; later runs count
    independently. Previously extracted blocks remain barriers.
    """
    if not hits:
        return []
    source = all_hits if all_hits is not None else hits
    intervals = sorted({h.competition_interval for h in source}, reverse=strand == "-")
    block_rank = {interval: i for i, interval in enumerate(intervals)}
    representatives = {}
    def representative_key(hit):
        return interval_score_source_key(hit.alignment), hit.exon_number
    for hit in source:
        block = block_rank[hit.competition_interval]
        current = representatives.get(block)
        if current is None or representative_key(hit) < representative_key(current):
            representatives[block] = hit
    ordered = sorted(hits, key=lambda h: (*oriented_hit_interval(h, strand), h.exon_number))
    blocks = [block_rank[h.competition_interval] for h in ordered]
    numbers = sorted({h.exon_number for h in ordered})
    ranks = {number: i + 1 for i, number in enumerate(numbers)}
    dp = [0.0] * len(ordered)
    count, inserted, previous = [0] * len(ordered), [0] * len(ordered), [-1] * len(ordered)
    def rank_key(i):
        return dp[i], count[i], -inserted[i], -i
    global_index = _GenePredecessorIndex(blocks, [h.exon_number for h in ordered], rank_key)
    local_tree = [-1] * (len(numbers) + 1)
    pending = []
    last_block = None
    segment_start = blocks[0]
    last_occurrence: OrderedDict[int, int] = OrderedDict()
    predecessor_ranges = []
    for i, hit in enumerate(ordered):
        block = blocks[i]
        if block != last_block:
            if last_block is not None:
                if block != last_block + 1:
                    # An extracted copy owns the missing block; never reinsert it.
                    segment_start = block
                    last_occurrence.clear()
                else:
                    number = representatives[last_block].exon_number
                    last_occurrence[number] = last_block
                    last_occurrence.move_to_end(number)
            # Gap distinct count = number of last occurrences after the
            # predecessor block. Consecutive predecessor ranges therefore
            # have constant costs 0, 50, ..., 1000; older ranges are dropped.
            recent = list(islice(reversed(last_occurrence.values()),
                                 MAX_GENE_INSERTION_UNIQUE_EXONS + 1))
            predecessor_ranges = []
            high = block - 1
            for unique_count, boundary in enumerate(recent + [segment_start]):
                if unique_count > MAX_GENE_INSERTION_UNIQUE_EXONS:
                    break
                low = max(segment_start, boundary)
                if low <= high:
                    predecessor_ranges.append((low, high, unique_count))
                high = boundary - 1
            local_tree = [-1] * (len(numbers) + 1)
            pending = []
            last_block = block
        start, end = oriented_hit_interval(hit, strand)
        while pending and pending[0][0] <= start:
            _end, prior = heapq.heappop(pending)
            position = ranks[ordered[prior].exon_number]
            while position < len(local_tree):
                local_tree[position] = global_index.better(local_tree[position], prior)
                position += position & -position
        dp[i], count[i] = hit.score, 1
        candidates = []
        for low, high, unique_count in predecessor_ranges:
            prior = global_index.query(low, high, hit.exon_number)
            if prior != -1:
                candidates.append((prior, hit.score, unique_count))
        position = ranks[hit.exon_number] - 1
        local_best = -1
        while position:
            local_best = global_index.better(local_best, local_tree[position])
            position -= position & -position
        if local_best != -1:
            candidates.append((local_best, 0.0, 0))
        for prior, reward, unique_count in candidates:
            score = dp[prior] + reward - GENE_INSERTION_COST * unique_count
            found, extra = count[prior] + 1, inserted[prior] + unique_count
            if (score, found, -extra) > (dp[i], count[i], -inserted[i]):
                dp[i], count[i], inserted[i], previous[i] = score, found, extra, prior
        global_index.update(block, hit.exon_number, i)
        heapq.heappush(pending, (end, i))
    best = max(range(len(ordered)), key=lambda i: (dp[i], count[i], -inserted[i],
               ordered[i].alignment.aligned_exon_bases))
    matched = []
    while best != -1:
        matched.append(best)
        best = previous[best]
    matched.reverse()
    chain = []
    for index, current in enumerate(matched):
        if index:
            for block in range(blocks[matched[index - 1]] + 1, blocks[current]):
                chain.append(replace(representatives[block], is_insertion=True))
        chain.append(ordered[current])
    return chain


def assign_merged_interval_gene_scores(
    grouped: Dict[Tuple[str, str, str], List[TranscriptHit]],
    annotation: Annotation,
) -> None:
    """Merge exon hits; each gene's longest eligible exon supplies its score.

    Isoforms of one gene reuse the same score object for each merged interval.
    No interval is assigned exclusively to a gene. Only candidate associations
    already eligible for transcript calling participate in this score table.
    """
    by_query: Dict[Tuple[str, str], Dict[int, Tuple[ExonAlignment, set[str]]]] = defaultdict(dict)
    gene_by_transcript: Dict[str, str] = {}
    for (query, strand, tid), hits in grouped.items():
        info = annotation.transcripts[tid]
        gene = info.gene_id or info.gene_name or "transcript:" + tid
        gene_by_transcript[tid] = gene
        evidence = by_query[(query, strand)]
        for hit in hits:
            key = id(hit.alignment)
            if key not in evidence:
                evidence[key] = (hit.alignment, set())
            evidence[key][1].add(gene)

    score_by_alignment: Dict[Tuple[int, str], MergedIntervalGeneScore] = {}
    for evidence in by_query.values():
        ordered = sorted(evidence.values(), key=lambda item: (
            item[0].query_start, item[0].query_end, item[0].exon_id
        ))
        block: List[Tuple[ExonAlignment, set[str]]] = []
        start, end = 0, -1

        def record_block() -> None:
            if not block:
                return
            best_by_gene: Dict[str, ExonAlignment] = {}
            for aln, genes in block:
                for gene in genes:
                    current = best_by_gene.get(gene)
                    if current is None or interval_score_source_key(aln) < interval_score_source_key(current):
                        best_by_gene[gene] = aln
            scores = {
                gene: MergedIntervalGeneScore(start, end, gene, aln.exon_score,
                                              aln.exon_id, aln.exon_length)
                for gene, aln in best_by_gene.items()
            }
            for aln, genes in block:
                for gene in genes:
                    score_by_alignment[(id(aln), gene)] = scores[gene]

        for aln, genes in ordered:
            if block and aln.query_start >= end:
                record_block()
                block = []
            if not block:
                start, end = aln.query_start, aln.query_end
            else:
                end = max(end, aln.query_end)
            block.append((aln, genes))
        record_block()

    for (_query, _strand, tid), hits in grouped.items():
        gene = gene_by_transcript[tid]
        for hit in hits:
            hit.interval_score = score_by_alignment[(id(hit.alignment), gene)]


def interval_score_source_key(alignment: ExonAlignment) -> tuple:
    # Full reference exon length determines the supplier, not query span or
    # similarity. Equal lengths prefer score, then stable identifiers/positions.
    return (-alignment.exon_length, -alignment.exon_score, alignment.exon_id,
            alignment.query_start, alignment.query_end)


def group_transcript_hits(
    alignments: Sequence[ExonAlignment],
    annotation: Annotation,
    expected_sets: Dict[str, set[int]],
) -> Dict[Tuple[str, str, str], List[TranscriptHit]]:
    # query/strand/transcript -> transcript hits
    grouped: Dict[Tuple[str, str, str], List[TranscriptHit]] = defaultdict(list)
    for aln in alignments:
        for (
            tid,
            exon_num,
            reference_chrom,
            reference_start,
            reference_end,
        ) in annotation.exon_to_transcripts.get(aln.exon_id, []):
            expected_nums = expected_sets.get(tid)
            if not expected_nums:
                continue
            if exon_num not in expected_nums:
                # When --truncate is used, only score exons from the truncation interval.
                continue
            grouped[(aln.query_id, aln.strand, tid)].append(
                TranscriptHit(
                    exon_number=exon_num,
                    alignment=aln,
                    reference_chrom=reference_chrom,
                    reference_start=reference_start,
                    reference_end=reference_end,
                )
            )

    return grouped


def build_transcript_calls(
    alignments: Sequence[ExonAlignment],
    annotation: Annotation,
    expected_by_transcript: Dict[str, List[int]],
    protein_bonus: float,
    complete_bonus: float,
    mane_bonus: float,
    max_chains_per_transcript: int,
    expected_sets: Optional[Dict[str, set[int]]] = None,
) -> List[TranscriptCall]:
    if expected_sets is None:
        expected_sets = {tid: set(nums) for tid, nums in expected_by_transcript.items() if nums}
    grouped = group_transcript_hits(alignments, annotation, expected_sets)
    assign_merged_interval_gene_scores(grouped, annotation)
    return build_calls_from_grouped_hits(grouped, annotation, expected_by_transcript,
                                        protein_bonus, complete_bonus, mane_bonus,
                                        max_chains_per_transcript)


def build_calls_from_grouped_hits(
    grouped: Dict[Tuple[str, str, str], List[TranscriptHit]],
    annotation: Annotation,
    expected_by_transcript: Dict[str, List[int]],
    protein_bonus: float,
    complete_bonus: float,
    mane_bonus: float,
    max_chains_per_transcript: int,
) -> List[TranscriptCall]:
    """Extract chains from hits with their original merged-interval gene scores."""
    calls: List[TranscriptCall] = []
    for (_query, strand, tid), hits in grouped.items():
        available = list(hits)
        info = annotation.transcripts[tid]
        copy_windows = exon_copy_windows(hits, strand) if info.model_type != "full_gene" else None
        chains_made = 0
        while available and chains_made < max_chains_per_transcript:
            info = annotation.transcripts[tid]
            if info.model_type == "full_gene":
                chain = best_gene_chain_for_hits(available, strand, all_hits=hits)
            else:
                chain = best_chain_for_hits(available, strand, copy_windows)
            if not chain:
                break
            if info.model_type != "full_gene":
                by_exon: Dict[int, TranscriptHit] = {}
                for hit in chain:
                    current = by_exon.get(hit.exon_number)
                    if current is None or hit.score > current.score:
                        by_exon[hit.exon_number] = hit
                chain = sorted(by_exon.values(), key=lambda h: h.exon_number)
            if not chain:
                break
            info = annotation.transcripts[tid]
            calls.append(
                TranscriptCall(
                    transcript_id=tid,
                    transcript_info=info,
                    hits=chain,
                    expected_exon_numbers=expected_by_transcript[tid],
                    protein_bonus=protein_bonus,
                    complete_bonus=complete_bonus,
                    mane_bonus=mane_bonus,
                )
            )
            used_intervals = {h.competition_interval for h in chain}
            available = [h for h in available if h.competition_interval not in used_intervals]
            chains_made += 1
    return calls


def call_query_bounds(call: TranscriptCall) -> Tuple[int, int]:
    return (
        min(h.alignment.query_start for h in call.hits),
        max(h.alignment.query_end for h in call.hits),
    )


def transcript_call_rank_key(
    call: TranscriptCall, prefer_mane: bool
) -> Tuple[int, float]:
    # Stage 1 passes prefer_mane=False and contains only full-gene models.
    # Stage 2 contains only real isoforms and preserves MANE priority.
    return (int(prefer_mane and call.transcript_info.is_mane), call.weighted_score)



class _TranscriptCallHeapItem:
    __slots__ = ("rank", "serial", "call")

    def __init__(
        self,
        rank: Tuple[int, float],
        serial: int,
        call: TranscriptCall,
    ) -> None:
        self.rank = rank
        self.serial = serial
        self.call = call

    def __lt__(self, other: "_TranscriptCallHeapItem") -> bool:
        # heapq is a min-heap; reverse rank comparison gives maximum priority.
        if self.rank != other.rank:
            return self.rank > other.rank
        return self.serial < other.serial


def label_tied_calls(calls: Sequence[TranscriptCall]) -> List[TranscriptCall]:
    """Label exon-overlap components within one equally ranked batch.

    Disjoint loci are separate copies, even if their scores are identical.
    Components describe connected ambiguity; members need not all overlap each
    other directly. A sweep joins each interval to the furthest-reaching prior
    overlapping interval, avoiding all-pairs comparisons.
    """
    parent = list(range(len(calls)))

    def root(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    intervals = sorted(
        (*hit.competition_interval, i)
        for i, call in enumerate(calls)
        for hit in call.hits
    )
    furthest_end, owner = -1, 0
    for start, end, i in intervals:
        if start < furthest_end:
            parent[root(i)] = root(owner)
        if start >= furthest_end or end > furthest_end:
            furthest_end, owner = end, i

    groups: Dict[int, List[int]] = defaultdict(list)
    for i in range(len(calls)):
        groups[root(i)].append(i)
    result = list(calls)
    for members in groups.values():
        group_id = ""
        if len(members) > 1:
            identities = sorted(
                (
                    calls[i].transcript_info.gene_id_full,
                    calls[i].transcript_info.transcript_id_full,
                    sorted(
                        (h.exon_number, h.alignment.exon_id,
                         h.alignment.query_start, h.alignment.query_end)
                        for h in calls[i].hits
                    ),
                )
                for i in members
            )
            first = calls[members[0]]
            payload = json.dumps(
                [first.query_id, first.strand, identities], separators=(",", ":")
            ).encode("utf-8")
            group_id = "tie_" + hashlib.sha256(payload).hexdigest()[:20]
        for i in members:
            result[i] = replace(
                calls[i], tie_group_id=group_id, tie_count=len(members)
            )
    return result


def choose_gene_isoforms(
    calls: Sequence[TranscriptCall], trace: Optional[List[dict]] = None,
) -> Tuple[List[TranscriptCall], List[TranscriptCall]]:
    """Keep one longest isoform per gene in each equal-rank overlap component.

    Length is only a within-gene tie breaker: it must not resolve gene ties or
    change MANE/score priority. Disjoint components remain independent copies.
    Losing isoforms are deferred so unoccupied residual blocks can be rescored.
    """
    by_component_gene: Dict[tuple, List[TranscriptCall]] = defaultdict(list)
    for index, (call, labeled) in enumerate(zip(calls, label_tied_calls(calls))):
        info = call.transcript_info
        component = labeled.tie_group_id or ("unique", index)
        gene = info.gene_id or info.gene_name or "transcript:" + call.transcript_id
        by_component_gene[(component, gene)].append(call)

    def isoform_key(call: TranscriptCall) -> tuple:
        return (-call.transcript_info.exon_length, call.transcript_id,
                tuple(sorted((h.exon_number, h.alignment.exon_id,
                              h.alignment.query_start, h.alignment.query_end) for h in call.hits)))

    kept, deferred = [], []
    for alternatives in by_component_gene.values():
        ordered = sorted(alternatives, key=isoform_key)
        winner = ordered[0]
        kept.append(winner)
        deferred.extend(ordered[1:])
        if trace is not None:
            for loser in ordered[1:]:
                trace.append({"event": "isoform_tiebreak", "call": loser, "winner": winner})
    return kept, deferred


def resolve_transcript_overlaps_sparse(
    calls: Sequence[TranscriptCall], prefer_mane: bool,
    trace: Optional[List[dict]] = None,
) -> List[TranscriptCall]:
    """Resolve competing ranks lazily, accepting a whole tied batch together.

    Ordinary ranks are upper bounds until a call reaches the heap top.
    Insertion-bearing calls are refreshed when occupancy changes, since
    removing a negative insertion can raise their score. Tied
    gene alternatives never trim each other. Same-gene isoform ties prefer the
    longer annotated transcript. Floating scores use a 1e-12 tolerance.
    """
    occupied = NonOverlappingIntervals()
    heap = [
        _TranscriptCallHeapItem(
            transcript_call_rank_key(call, prefer_mane), serial, call
        )
        for serial, call in enumerate(calls)
        if call.hits
    ]
    heapq.heapify(heap)
    next_serial = len(heap)
    selected: List[TranscriptCall] = []
    has_insertions = any(call.insertion_intervals for call in calls)
    insertion_dirty = False
    def prepare_top() -> None:
        nonlocal next_serial, insertion_dirty, has_insertions
        # Removing a negative insertion can increase rank. Refresh these
        # candidates once after occupancy changes; ordinary calls remain lazy.
        if insertion_dirty:
            insertion_dirty = False
            updated = []
            changed = False
            for item in heap:
                call = item.call
                if not call.insertion_intervals:
                    updated.append(item)
                    continue
                retained = [h for h in call.hits if not occupied.overlaps(*h.competition_interval)]
                if len(retained) == len(call.hits):
                    updated.append(item)
                    continue
                changed = True
                trimmed = call.recalculated_with_hits(retained) if any(not h.is_insertion for h in retained) else None
                if trimmed is not None and not trimmed.valid_insertion_runs:
                    trimmed = None
                if trimmed is not None:
                    updated.append(_TranscriptCallHeapItem(
                        transcript_call_rank_key(trimmed, prefer_mane), next_serial, trimmed))
                    next_serial += 1
                if trace is not None:
                    trace.append({"event": "trim", "call": call, "remaining": trimmed})
            if changed:
                heap[:] = updated
                heapq.heapify(heap)
            has_insertions = any(item.call.insertion_intervals for item in heap)
        while heap:
            call = heap[0].call
            trimmed_hits = [
                h for h in call.hits
                if not occupied.overlaps(*h.competition_interval)
            ]
            if len(trimmed_hits) == len(call.hits):
                return
            heapq.heappop(heap)
            trimmed = None
            if trimmed_hits:
                trimmed = call.recalculated_with_hits(trimmed_hits)
                heapq.heappush(heap, _TranscriptCallHeapItem(
                    transcript_call_rank_key(trimmed, prefer_mane), next_serial, trimmed
                ))
                next_serial += 1
            if trace is not None:
                trace.append({"event": "trim", "call": call, "remaining": trimmed})

    def same_score(a: float, b: float) -> bool:
        return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)

    while heap:
        prepare_top()
        if not heap:
            break
        best = heap[0].rank
        candidates = []
        # Collect all equal aggregate scores before consuming their intervals.
        while heap:
            prepare_top()
            if not heap or heap[0].rank[0] != best[0] or not same_score(heap[0].rank[1], best[1]):
                break
            candidates.append(heapq.heappop(heap))
        batch, deferred = choose_gene_isoforms([item.call for item in candidates], trace=trace)
        for call in deferred:
            heapq.heappush(heap, _TranscriptCallHeapItem(
                transcript_call_rank_key(call, prefer_mane), next_serial, call))
            next_serial += 1
        labeled = label_tied_calls(batch)
        selected.extend(labeled)
        if trace is not None:
            trace.append({"event": "select", "calls": batch, "selected": labeled})
        # The interval index requires disjoint entries; tied alternatives may
        # overlap or contain one another, so insert their union only once.
        for start, end in merge_intervals([
            h.competition_interval
            for call in batch for h in call.hits
        ]):
            occupied.add(start, end)
        insertion_dirty = has_insertions
    return selected


def resolve_call_overlaps(
    calls: Sequence[TranscriptCall],
    prefer_mane: bool = True,
) -> List[TranscriptCall]:
    """Resolve exon-chain competitors within one caller stage.

    Every candidate retains its gene and score. The hierarchical caller uses
    this separately for full genes and for each parent's real isoforms.
    """
    grouped: Dict[Tuple[str, str], List[TranscriptCall]] = defaultdict(list)
    for call in calls:
        if call.hits:
            grouped[(call.query_id, call.strand)].append(call)

    selected: List[TranscriptCall] = []
    for (query_id, strand), group_calls in grouped.items():
        large_group = len(group_calls) >= 10_000
        if large_group:
            print(
                f"Resolving large query/strand group {query_id} {strand}: "
                f"raw_transcript_calls={len(group_calls)}",
                file=sys.stderr,
            )
        group_selected = resolve_transcript_overlaps_sparse(group_calls, prefer_mane)
        selected.extend(group_selected)
        if large_group:
            print(
                f"Finished query/strand group {query_id} {strand}: "
                f"candidate_transcript_calls={len(group_calls)}, "
                f"selected_transcript_calls={len(group_selected)}",
                file=sys.stderr,
            )
    selected.sort(key=call_output_sort_key, reverse=True)
    return selected


def call_output_sort_key(call: TranscriptCall) -> tuple:
    if call.gene_index:
        return (call.gene_sort_score, call.query_id, -call.gene_start, -call.gene_end,
                call.gene_index, int(call.transcript_info.model_type == "full_gene"),
                call.weighted_score, call.transcript_id)
    return (
        call.weighted_score,
        call.query_id,
        call.transcript_id,
    )


def build_hierarchical_calls(
    alignments: Sequence[ExonAlignment], annotation: Annotation,
    expected_by_transcript: Dict[str, List[int]], protein_bonus: float,
    complete_bonus: float, mane_bonus: float, max_chains_per_transcript: int,
    prefer_mane: bool,
    diagnostic_batches: Optional[List[dict]] = None,
) -> Tuple[List[TranscriptCall], int]:
    """Lift full genes first, then compete real isoforms inside each parent.

    Interval scores are computed once across all eligible evidence. Children
    reuse those scores and only the parent's owned intervals and source genes;
    they cannot borrow evidence from another gene copy or an intronic gene.
    """
    expected_sets = {tid: set(nums) for tid, nums in expected_by_transcript.items() if nums}
    grouped = group_transcript_hits(alignments, annotation, expected_sets)
    assign_merged_interval_gene_scores(grouped, annotation)
    full_groups = {key: hits for key, hits in grouped.items()
                   if annotation.transcripts[key[2]].model_type == "full_gene"}
    # Gene identity must not depend on whether the reference annotation is
    # protein-coding: a partial coding paralog must not receive a 10x boost
    # over an exact complete pseudogene. Coding priority belongs to stage 2.
    raw_genes = build_calls_from_grouped_hits(
        full_groups, annotation, expected_by_transcript, 1.0,
        complete_bonus, mane_bonus, max_chains_per_transcript)
    if diagnostic_batches is not None:
        diagnostic_batches.append(dict(stage="gene", parent="", calls=raw_genes, prefer_mane=False))
    genes = resolve_call_overlaps(raw_genes, prefer_mane=False)
    raw_count = len(raw_genes)

    # Map each reference model to its own isoforms, including separate reference
    # loci carrying the same gene ID. Avoid scanning all alignments per parent.
    real_models = {}
    for tid in {key[2] for key in full_groups}:
        gene = annotation.transcripts[tid].gene_id
        members = set()
        for exon in annotation.exons_by_transcript[tid]:
            for child, _num, chrom, start, end in annotation.exon_to_transcripts[exon.exon_id]:
                info = annotation.transcripts[child]
                if (info.model_type == "transcript" and info.gene_id == gene
                        and (chrom, start, end) == (exon.chrom, exon.start0, exon.end0)):
                    members.add(child)
        real_models[tid] = members
    child_index = defaultdict(list)
    for key, hits in grouped.items():
        info = annotation.transcripts[key[2]]
        if info.model_type == "transcript":
            for hit in hits:
                child_index[(key[0], key[1], info.gene_id, hit.competition_interval)].append((key, hit))

    selected = []
    for parents in iter_assignment_groups(genes):
        intervals = {h.competition_interval for parent in parents for h in parent.hits}
        start = min(a for a, _b in intervals)
        end = max(b for _a, b in intervals)
        query, strand = parents[0].query_id, parents[0].strand
        identity = (query, strand, sorted(p.transcript_id for p in parents), sorted(intervals))
        parent_key = hashlib.sha256(repr(identity).encode()).hexdigest()[:24]
        metadata = dict(gene_index=parent_key, gene_sort_score=max(p.weighted_score for p in parents),
                        gene_start=start, gene_end=end)
        selected.extend(replace(parent, **metadata) for parent in parents)
        allowed_by_gene = defaultdict(set)
        for parent in parents:
            allowed_by_gene[parent.transcript_info.gene_id].update(real_models[parent.transcript_id])
        child_groups = defaultdict(list)
        for gene, allowed in allowed_by_gene.items():
            for interval in intervals:
                for key, hit in child_index.get((query, strand, gene, interval), ()):
                    if key[2] in allowed:
                        child_groups[key].append(hit)
        # Fresh DP inside each parent is necessary: trimming a globally chained
        # transcript could retain exon links that cross two different copies.
        raw_children = build_calls_from_grouped_hits(
            child_groups, annotation, expected_by_transcript, protein_bonus,
            complete_bonus, mane_bonus, max_chains_per_transcript)
        raw_count += len(raw_children)
        if diagnostic_batches is not None:
            diagnostic_batches.append(dict(stage="transcript", parent=parent_key,
                                           calls=raw_children, prefer_mane=prefer_mane))
        children = resolve_call_overlaps(raw_children, prefer_mane=prefer_mane)
        selected.extend(replace(child, **metadata) for child in children)
    selected.sort(key=call_output_sort_key, reverse=True)
    return selected, raw_count


def call_models(
    alignments, annotation, expected_by_transcript, protein_bonus, complete_bonus,
    mane_bonus, max_chains_per_transcript, prefer_mane, two_scale=False,
):
    if two_scale:
        return build_hierarchical_calls(
            alignments, annotation, expected_by_transcript, protein_bonus,
            complete_bonus, mane_bonus, max_chains_per_transcript, prefer_mane)
    raw = build_transcript_calls(alignments, annotation, expected_by_transcript,
                                protein_bonus, complete_bonus, mane_bonus,
                                max_chains_per_transcript)
    return resolve_call_overlaps(raw, prefer_mane=prefer_mane), len(raw)


# Forked workers inherit these large, read-only objects copy-on-write. This
# avoids reparsing or pickling the full GFF annotation for every contig task.
_WORKER_GROUPS: Dict[Tuple[str, str], List[ExonAlignment]] = {}
_WORKER_ANNOTATION: Optional[Annotation] = None
_WORKER_EXPECTED: Dict[str, List[int]] = {}
_WORKER_EXPECTED_SETS: Dict[str, set[int]] = {}
_WORKER_PARAMS: tuple = (10.0, 2.0, 1.0, 10, True, False)
_WORKER_SHARD_DIR = ""


def _call_query_strand_worker(key: Tuple[str, str]) -> Tuple[int, List[TranscriptCall]]:
    if _WORKER_ANNOTATION is None:
        raise RuntimeError("parallel caller worker was not initialized")
    selected, count = call_models(
        _WORKER_GROUPS[key], _WORKER_ANNOTATION, _WORKER_EXPECTED, *_WORKER_PARAMS)
    return count, selected


def _call_query_strand_shard_worker(
    task: Tuple[int, Tuple[str, str]],
) -> Tuple[int, int, str]:
    """Call one query/strand group and write results before returning to parent."""
    group_index, key = task
    raw_count, selected = _call_query_strand_worker(key)
    shard_path = os.path.join(_WORKER_SHARD_DIR, f"group_{group_index:09d}.pkl")
    with open(shard_path, "wb") as handle:
        for call_index, call in enumerate(selected):
            record = (call_output_sort_key(call), (group_index, call_index), call)
            pickle.dump(record, handle, protocol=pickle.HIGHEST_PROTOCOL)
    selected_count = len(selected)
    del selected
    return raw_count, selected_count, shard_path


class _ShardHeapItem:
    __slots__ = ("record", "source")

    def __init__(self, record, source: int) -> None:
        self.record = record
        self.source = source

    def __lt__(self, other: "_ShardHeapItem") -> bool:
        # Calls were historically globally sorted with reverse=True. Preserve
        # that order and Python sort stability using the original group/call
        # position as the secondary key.
        if self.record[0] != other.record[0]:
            return self.record[0] > other.record[0]
        return self.record[1] < other.record[1]


def _load_pickle_record(handle):
    try:
        return pickle.load(handle)
    except EOFError:
        return None


def _merge_shard_batch(input_paths: Sequence[str], output_path: str) -> None:
    handles = [open(path, "rb") for path in input_paths]
    heap: List[_ShardHeapItem] = []
    try:
        for source, handle in enumerate(handles):
            record = _load_pickle_record(handle)
            if record is not None:
                heapq.heappush(heap, _ShardHeapItem(record, source))
        with open(output_path, "wb") as output:
            while heap:
                item = heapq.heappop(heap)
                pickle.dump(item.record, output, protocol=pickle.HIGHEST_PROTOCOL)
                record = _load_pickle_record(handles[item.source])
                if record is not None:
                    heapq.heappush(heap, _ShardHeapItem(record, item.source))
    finally:
        for handle in handles:
            handle.close()


def merge_sorted_call_shards(
    shard_paths: Sequence[str],
    shard_dir: str,
    max_open_files: int = 64,
) -> str:
    """Disk-backed stable merge; memory holds one call per open shard."""
    paths = list(shard_paths)
    if not paths:
        empty_path = os.path.join(shard_dir, "calls.sorted.pkl")
        open(empty_path, "wb").close()
        return empty_path

    merge_round = 0
    while len(paths) > 1:
        next_paths: List[str] = []
        for batch_index in range(0, len(paths), max_open_files):
            batch = paths[batch_index : batch_index + max_open_files]
            if len(batch) == 1:
                next_paths.append(batch[0])
                continue
            merged = os.path.join(
                shard_dir,
                f"merge_{merge_round:03d}_{batch_index // max_open_files:06d}.pkl",
            )
            _merge_shard_batch(batch, merged)
            for path in batch:
                os.unlink(path)
            next_paths.append(merged)
        paths = next_paths
        merge_round += 1

    final_path = os.path.join(shard_dir, "calls.sorted.pkl")
    if paths[0] != final_path:
        os.replace(paths[0], final_path)
    return final_path


def iter_sorted_call_shard(path: str) -> Iterator[TranscriptCall]:
    with open(path, "rb") as handle:
        while True:
            record = _load_pickle_record(handle)
            if record is None:
                return
            yield record[2]


def call_transcripts_grouped(
    alignments: Sequence[ExonAlignment],
    annotation: Annotation,
    expected_by_transcript: Dict[str, List[int]],
    protein_bonus: float,
    complete_bonus: float,
    mane_bonus: float,
    max_chains_per_transcript: int,
    prefer_mane: bool,
    threads: int,
    two_scale: bool = False,
) -> Tuple[List[TranscriptCall], int]:
    """Call independently by query/strand, optionally using forked workers."""
    if threads <= 1:
        return call_models(alignments, annotation, expected_by_transcript, protein_bonus,
                           complete_bonus, mane_bonus, max_chains_per_transcript,
                           prefer_mane, two_scale)

    groups: Dict[Tuple[str, str], List[ExonAlignment]] = defaultdict(list)
    for aln in alignments:
        groups[(aln.query_id, aln.strand)].append(aln)
    if len(groups) <= 1:
        print(
            "Only one query/strand group is available; transcript calling remains single-process",
            file=sys.stderr,
        )
        return call_models(alignments, annotation, expected_by_transcript, protein_bonus,
                           complete_bonus, mane_bonus, max_chains_per_transcript,
                           prefer_mane, two_scale)

    if "fork" not in mp.get_all_start_methods():
        print(
            "WARNING: multiprocessing requires fork to share the large annotation; using one process",
            file=sys.stderr,
        )
        return call_models(alignments, annotation, expected_by_transcript, protein_bonus,
                           complete_bonus, mane_bonus, max_chains_per_transcript,
                           prefer_mane, two_scale)

    worker_count = min(max(1, threads), len(groups))
    print(
        f"Calling transcripts across {len(groups)} query/strand groups with {worker_count} processes",
        file=sys.stderr,
    )
    global _WORKER_GROUPS, _WORKER_ANNOTATION, _WORKER_EXPECTED, _WORKER_PARAMS
    _WORKER_GROUPS = dict(groups)
    _WORKER_ANNOTATION = annotation
    _WORKER_EXPECTED = expected_by_transcript
    _WORKER_PARAMS = (
        protein_bonus,
        complete_bonus,
        mane_bonus,
        max_chains_per_transcript,
        prefer_mane,
        two_scale,
    )
    try:
        context = mp.get_context("fork")
        with context.Pool(processes=worker_count) as pool:
            results = pool.map(_call_query_strand_worker, sorted(groups))
    finally:
        _WORKER_GROUPS = {}
        _WORKER_ANNOTATION = None
        _WORKER_EXPECTED = {}

    raw_count = sum(count for count, _calls in results)
    calls = [call for _count, selected in results for call in selected]
    calls.sort(key=call_output_sort_key, reverse=True)
    return calls, raw_count


def call_transcripts_grouped_to_disk(
    alignments: Sequence[ExonAlignment],
    annotation: Annotation,
    expected_by_transcript: Dict[str, List[int]],
    protein_bonus: float,
    complete_bonus: float,
    mane_bonus: float,
    max_chains_per_transcript: int,
    prefer_mane: bool,
    threads: int,
    shard_dir: str,
    two_scale: bool = False,
) -> Optional[Tuple[str, int, int]]:
    """Parallel caller whose workers return paths instead of large call lists."""
    if threads <= 1 or "fork" not in mp.get_all_start_methods():
        return None

    groups: Dict[Tuple[str, str], List[ExonAlignment]] = defaultdict(list)
    for aln in alignments:
        groups[(aln.query_id, aln.strand)].append(aln)
    if not groups:
        final_path = merge_sorted_call_shards([], shard_dir)
        return final_path, 0, 0

    worker_count = min(max(1, threads), len(groups))
    print(
        f"Calling {len(groups)} query/strand groups with {worker_count} processes; "
        f"completed groups write temporary result shards to {shard_dir}",
        file=sys.stderr,
    )
    global _WORKER_GROUPS, _WORKER_ANNOTATION, _WORKER_EXPECTED, _WORKER_EXPECTED_SETS
    global _WORKER_PARAMS, _WORKER_SHARD_DIR
    _WORKER_GROUPS = dict(groups)
    _WORKER_ANNOTATION = annotation
    _WORKER_EXPECTED = expected_by_transcript
    _WORKER_EXPECTED_SETS = {
        tid: set(nums) for tid, nums in expected_by_transcript.items() if nums
    }
    _WORKER_PARAMS = (
        protein_bonus,
        complete_bonus,
        mane_bonus,
        max_chains_per_transcript,
        prefer_mane,
        two_scale,
    )
    _WORKER_SHARD_DIR = shard_dir
    results: List[Tuple[int, int, str]] = []
    tasks = list(enumerate(sorted(groups)))
    try:
        context = mp.get_context("fork")
        # Recycle after every group so memory retained by Python's allocator is
        # returned to the operating system as soon as that shard is complete.
        with context.Pool(processes=worker_count, maxtasksperchild=1) as pool:
            for result in pool.imap_unordered(
                _call_query_strand_shard_worker, tasks, chunksize=1
            ):
                results.append(result)
    finally:
        _WORKER_GROUPS = {}
        _WORKER_ANNOTATION = None
        _WORKER_EXPECTED = {}
        _WORKER_EXPECTED_SETS = {}
        _WORKER_SHARD_DIR = ""

    raw_count = sum(raw for raw, _selected, _path in results)
    selected_count = sum(selected for _raw, selected, _path in results)
    sorted_path = merge_sorted_call_shards(
        [path for _raw, _selected, path in results], shard_dir
    )
    return sorted_path, raw_count, selected_count


def usable_temporary_directory(path: str) -> bool:
    return os.path.isdir(path) and os.access(path, os.W_OK | os.X_OK)


def choose_shard_parent(
    requested_tmp_dir: Optional[str],
    shard_storage: str,
    output_path: str,
) -> str:
    """Choose RAM, node-local scratch, or output storage for call shards."""
    if requested_tmp_dir:
        parent = os.path.abspath(os.path.expanduser(requested_tmp_dir))
        if not usable_temporary_directory(parent):
            raise SystemExit(
                f"ERROR: temporary directory is absent or not writable: {parent}"
            )
        print(f"Using explicitly requested transcript shard directory: {parent}", file=sys.stderr)
        return parent

    ram_parent = "/dev/shm"
    if shard_storage in {"auto", "ram"} and usable_temporary_directory(ram_parent):
        try:
            free_bytes = shutil.disk_usage(ram_parent).free
            free_text = f", free={free_bytes / (1024 ** 3):.2f} GiB"
        except OSError:
            free_text = ""
        print(
            f"Using RAM-backed transcript shard storage: {ram_parent}{free_text}",
            file=sys.stderr,
        )
        return ram_parent
    if shard_storage == "ram":
        raise SystemExit(
            "ERROR: --shard-storage ram requested, but /dev/shm is unavailable "
            "or not writable"
        )

    slurm_parent = os.environ.get("SLURM_TMPDIR", "")
    if slurm_parent:
        slurm_parent = os.path.abspath(os.path.expanduser(slurm_parent))
        if usable_temporary_directory(slurm_parent):
            print(
                f"Using SLURM node-local transcript shard storage: {slurm_parent}",
                file=sys.stderr,
            )
            return slurm_parent

    output_parent = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_parent, exist_ok=True)
    print(
        f"WARNING: using output filesystem for transcript shards: {output_parent}",
        file=sys.stderr,
    )
    return output_parent


def weighted_mean(values: Sequence[float], weights: Sequence[int]) -> float:
    total_w = sum(weights)
    if total_w == 0:
        return 0.0
    return sum(v * w for v, w in zip(values, weights)) / total_w


@dataclass(frozen=True)
class OverlapEvidence:
    query_start: int
    query_end: int
    candidates: Tuple[Tuple[str, str, str, str, int, int], ...]


class OverlapEvidenceIndex:
    """Compact interval index of every gene/transcript supported at a query exon."""

    def __init__(
        self,
        alignments: Sequence[ExonAlignment],
        annotation: Annotation,
    ) -> None:
        grouped: Dict[Tuple[str, str], List[OverlapEvidence]] = defaultdict(list)
        for alignment in alignments:
            candidates = set()
            for (
                transcript_id,
                _exon_number,
                reference_chrom,
                reference_start,
                reference_end,
            ) in annotation.exon_to_transcripts.get(alignment.exon_id, []):
                info = annotation.transcripts.get(transcript_id)
                if info is None:
                    continue
                gene_id = info.reported_gene_id or info.reported_gene_name
                if not gene_id:
                    continue
                candidates.add(
                    (
                        transcript_id,
                        gene_id,
                        info.reported_gene_name,
                        reference_chrom,
                        reference_start,
                        reference_end,
                    )
                )
            if candidates:
                grouped[(alignment.query_id, alignment.strand)].append(
                    OverlapEvidence(
                        query_start=alignment.query_start,
                        query_end=alignment.query_end,
                        candidates=tuple(sorted(candidates)),
                    )
                )

        self.groups: Dict[
            Tuple[str, str], Tuple[List[OverlapEvidence], List[int], List[int]]
        ] = {}
        for key, evidence in grouped.items():
            evidence.sort(key=lambda item: (item.query_start, item.query_end))
            starts = [item.query_start for item in evidence]
            prefix_max_end: List[int] = []
            maximum = -1
            for item in evidence:
                maximum = max(maximum, item.query_end)
                prefix_max_end.append(maximum)
            self.groups[key] = (evidence, starts, prefix_max_end)

    def overlapping(
        self,
        query_id: str,
        strand: str,
        start: int,
        end: int,
    ) -> List[OverlapEvidence]:
        indexed = self.groups.get((query_id, strand))
        if indexed is None:
            return []
        evidence, starts, prefix_max_end = indexed
        index = bisect_left(starts, end) - 1
        result: List[OverlapEvidence] = []
        while index >= 0 and prefix_max_end[index] > start:
            item = evidence[index]
            if item.query_end > start:
                result.append(item)
            index -= 1
        return result


OVERLAP_GROUP_HEADER = [
    "overlap_group_id",
    "query_id",
    "query_contig",
    "query_start",
    "query_end",
    "strand",
    "reported_transcript_id",
    "reported_gene_id",
    "reported_gene_name",
    "candidate_transcript_ids",
    "candidate_gene_ids",
    "candidate_gene_names",
    "overlap_regions",
    "candidate_reference_exons",
]


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def call_to_overlap_group_row(
    call: TranscriptCall,
    evidence_index: OverlapEvidenceIndex,
    group_number: int,
    parse_query_suffix: bool,
) -> List[str]:
    qloc = parse_query_location(call.query_id, parse_suffix=parse_query_suffix)
    call_intervals = [
        (hit.alignment.query_start, hit.alignment.query_end) for hit in call.hits
    ]
    candidates = set()
    overlap_intervals: List[Tuple[int, int]] = []
    for call_start, call_end in call_intervals:
        for evidence in evidence_index.overlapping(
            call.query_id, call.strand, call_start, call_end
        ):
            overlap_start = max(call_start, evidence.query_start)
            overlap_end = min(call_end, evidence.query_end)
            if overlap_end <= overlap_start:
                continue
            overlap_intervals.append(
                (qloc.start + overlap_start, qloc.start + overlap_end)
            )
            candidates.update(evidence.candidates)

    info = call.transcript_info
    reported_gene_id = info.reported_gene_id or info.reported_gene_name
    # A reported gene must remain a member even when no alternative alignment
    # survived the caller's initial filters.
    if not any(candidate[1] == reported_gene_id for candidate in candidates):
        for hit in call.hits:
            candidates.add(
                (
                    call.transcript_id,
                    reported_gene_id,
                    info.reported_gene_name,
                    hit.reference_chrom,
                    hit.reference_start,
                    hit.reference_end,
                )
            )

    candidates_sorted = sorted(candidates)
    gene_names = sorted({item[2] for item in candidates_sorted if item[2]})
    q_start, q_end = call_query_bounds(call)
    return [
        f"OG{group_number:09d}",
        call.query_id,
        qloc.contig,
        str(qloc.start + q_start),
        str(qloc.start + q_end),
        call.strand,
        call.transcript_id,
        reported_gene_id,
        info.reported_gene_name,
        ",".join(sorted({item[0] for item in candidates_sorted})),
        ",".join(sorted({item[1] for item in candidates_sorted})),
        ",".join(gene_names),
        ",".join(
            f"{start}-{end}" for start, end in merge_intervals(overlap_intervals)
        ),
        ",".join(
            sorted(
                {
                    f"{item[3]}:{item[4]}-{item[5]}"
                    for item in candidates_sorted
                }
            )
        ),
    ]


OUTPUT_HEADER = [
    "GENE_index",
    "transcript_id",
    "gene_id",
    "gene_name",
    "transcript_type",
    "ifmane_transcript",
    "weighted_score",
    "raw_exon_score",
    "total_aligned_bases",
    "mean_identity",
    "mean_exon_coverage",
    "expected_exons",
    "found_exons",
    "missing_exons",
    "fraction_expected_found",
    "found_exon_numbers",
    "query_id",
    "query_contig",
    "query_start",
    "query_end",
    "strand",
    "exon_query_coordinates",
    "exon_reference_coordinates",
    "exon_ids",
    "exon_coverages",
    "exon_identities",
    "alignment_AS",
    "call_status",
    "tie_group_id",
    "tie_count",
    "assignment_status",
    "pipeline_version",
    "merged_interval_coordinates",
    "merged_interval_gene_scores",
    "merged_interval_score_exon_ids",
    "merged_interval_score_exon_lengths",
    "transcript_exon_length",
    "model_type",
    "inserted_exons",
    "insertion_penalty",
    "inserted_exon_numbers",
    "inserted_exon_query_coordinates",
    "inserted_query_blocks",
    "insertion_run_unique_exons",
]


# Preserve gene/exon correspondence even when two alternatives have identical
# evidence. Other summary fields may appear once when their values all agree.
PAIRED_OUTPUT_FIELDS = (
    "transcript_id", "gene_id", "gene_name", "found_exon_numbers",
    "total_aligned_bases", "mean_identity", "mean_exon_coverage",
    "exon_query_coordinates", "exon_reference_coordinates", "exon_ids",
    "exon_coverages", "exon_identities", "alignment_AS",
    "merged_interval_coordinates", "merged_interval_gene_scores",
    "merged_interval_score_exon_ids", "merged_interval_score_exon_lengths",
    "model_type", "inserted_exon_numbers", "inserted_exon_query_coordinates",
    "insertion_run_unique_exons",
)


def call_to_output_row(call: TranscriptCall, parse_query_suffix: bool = True) -> List[str]:
    hits = sorted(call.hits, key=lambda h: min(h.alignment.query_start, h.alignment.query_end))
    qloc = parse_query_location(call.query_id, parse_suffix=parse_query_suffix)
    abs_starts = [qloc.start + h.alignment.query_start for h in hits]
    abs_ends = [qloc.start + h.alignment.query_end for h in hits]
    q_start = min(abs_starts) if abs_starts else 0
    q_end = max(abs_ends) if abs_ends else 0
    weights = [h.alignment.aligned_exon_bases for h in hits]
    mean_id = weighted_mean([h.alignment.percent_identity for h in hits], weights)
    mean_cov = sum(h.alignment.exon_coverage for h in hits) / len(hits) if hits else 0.0
    frac = call.found_expected_count / call.expected_count if call.expected_count else 0.0
    status = "complete" if call.missing_count == 0 else "partial"
    info = call.transcript_info
    if call.gene_index and info.model_type == "full_gene":
        q_start, q_end = qloc.start + call.gene_start, qloc.start + call.gene_end
    return [
        call.gene_index,
        info.gene_id if call.gene_index and info.model_type == "full_gene" else call.transcript_id,
        info.reported_gene_id,
        info.reported_gene_name,
        info.transcript_type,
        "1" if info.is_mane else "0",
        f"{call.weighted_score:.6f}",
        f"{call.raw_score:.6f}",
        str(call.total_aligned_bases),
        f"{mean_id:.6f}",
        f"{mean_cov:.6f}",
        str(call.expected_count),
        str(call.found_expected_count),
        str(call.missing_count),
        f"{frac:.6f}",
        ",".join(str(h.exon_number) for h in hits),
        call.query_id,
        qloc.contig,
        str(q_start),
        str(q_end),
        call.strand,
        ",".join(f"{start}-{end}" for start, end in zip(abs_starts, abs_ends)),
        ",".join(
            f"{h.reference_chrom}:{h.reference_start}-{h.reference_end}"
            for h in hits
        ),
        ",".join(h.alignment.exon_id for h in hits),
        ",".join(f"{h.alignment.exon_coverage:.3f}" for h in hits),
        ",".join(f"{h.alignment.percent_identity:.3f}" for h in hits),
        ",".join(f"{h.alignment.AS:.3f}" for h in hits),
        status,
        call.tie_group_id,
        str(call.tie_count),
        "tied" if call.tie_count > 1 else "unique",
        PIPELINE_VERSION,
        ",".join(
            f"{qloc.start + s.start}-{qloc.start + s.end}"
            for s in call.merged_interval_scores
        ),
        ",".join(f"{s.score:.6f}" for s in call.merged_interval_scores),
        ",".join(s.source_exon_id for s in call.merged_interval_scores),
        ",".join(str(s.source_exon_length) for s in call.merged_interval_scores),
        str(info.exon_length),
        info.model_type,
        str(call.inserted_exon_count),
        f"{call.insertion_penalty:.6f}",
        ",".join(str(h.exon_number) for h in hits if h.is_insertion),
        ",".join(f"{qloc.start + a}-{qloc.start + b}" for a, b in sorted(call.insertion_intervals)),
        str(len(call.insertion_intervals)),
        ",".join(str(len(run)) for run in call.insertion_run_exons),
    ]


def iter_assignment_groups(calls: Iterable[TranscriptCall]) -> Iterator[List[TranscriptCall]]:
    """Collect tied components even when their members interleave in shard order.

    Only incomplete tie groups are buffered; unique calls remain streaming.
    Resolution has already removed lower-ranked competing evidence.
    """
    pending: Dict[str, List[TranscriptCall]] = {}
    for call in calls:
        if call.tie_count == 1:
            yield [call]
            continue
        if not call.tie_group_id:
            raise ValueError("tied call has no tie_group_id")
        group = pending.setdefault(call.tie_group_id, [])
        if group and group[0].tie_count != call.tie_count:
            raise ValueError("inconsistent tie_count within assignment group")
        group.append(call)
        if len(group) == call.tie_count:
            del pending[call.tie_group_id]
            yield sorted(group, key=lambda c: (c.transcript_info.gene_id, c.transcript_id,
                                               call_query_bounds(c)))
    if pending:
        raise ValueError("incomplete tied assignment groups: " + ",".join(sorted(pending)))


def assignment_to_output_row(
    calls: Sequence[TranscriptCall], parse_query_suffix: bool = True
) -> List[str]:
    """One union record, preserving each gene's evidence in existing columns.

    Semicolons separate genes in the same order in every paired column; commas
    separate exons within a gene. Equal summary values appear once, differing
    values use semicolons. Query bounds describe the union, while coordinate
    arrays retain each gene's exon/score correspondence. No JSON duplication.
    """
    if len(calls) == 1:
        return call_to_output_row(calls[0], parse_query_suffix)
    calls = sorted(calls, key=lambda c: (c.transcript_info.gene_id, c.transcript_id,
                                        call_query_bounds(c)))
    alternatives = [dict(zip(OUTPUT_HEADER, call_to_output_row(c, parse_query_suffix)))
                    for c in calls]
    row = {key: alternatives[0][key] if len({a[key] for a in alternatives}) == 1
           else ";".join(a[key] for a in alternatives)
           for key in OUTPUT_HEADER}
    for key in PAIRED_OUTPUT_FIELDS:
        row[key] = ";".join(a[key] for a in alternatives)
    row["query_start"] = str(min(int(a["query_start"]) for a in alternatives))
    row["query_end"] = str(max(int(a["query_end"]) for a in alternatives))
    row["tie_count"] = str(len(calls))
    row["assignment_status"] = "tied"
    return [row[key] for key in OUTPUT_HEADER]


def call_to_legacy_row(
    call: TranscriptCall,
    truncate_used: bool,
    parse_query_suffix: bool = True,
) -> List[str]:
    if call.tie_count > 1 or call.gene_index:
        raise SystemExit(
            "ERROR: legacy output cannot represent tied assignments or gene/transcript hierarchy; use --output-format extended"
        )
    hits = sorted(call.hits, key=lambda h: min(h.alignment.query_start, h.alignment.query_end))
    info = call.transcript_info
    qloc = parse_query_location(call.query_id, parse_suffix=parse_query_suffix)
    abs_starts = [qloc.start + h.alignment.query_start for h in hits]
    abs_ends = [qloc.start + h.alignment.query_end for h in hits]
    q_start = min(abs_starts) if abs_starts else 0
    q_end = max(abs_ends) if abs_ends else 0
    simi = call.raw_score / call.expected_count if call.expected_count else 0.0
    transcript_id = call.transcript_id
    if truncate_used and call.expected_exon_numbers:
        transcript_id = f"{transcript_id}(:{min(call.expected_exon_numbers)}-{max(call.expected_exon_numbers)})"
    return [
        transcript_id,
        info.reported_gene_id,
        info.reported_gene_name,
        info.transcript_type,
        str(call.total_aligned_bases),
        f"{simi:.6f}",
        str(call.expected_count),
        str(call.missing_count),
        ",".join(str(h.exon_number) for h in hits),
        qloc.contig,
        str(q_start),
        str(q_end),
        ",".join(str(x) for x in abs_starts),
        ",".join(str(x) for x in abs_ends),
        ",".join(h.alignment.strand for h in hits),
        ",".join(h.alignment.exon_id for h in hits),
    ]


def is_pseudofragment_assignment(groups: Sequence[Sequence[TranscriptCall]]) -> bool:
    """Candidate gene-like fragment, classified from the whole selected parent.

    This is a structural screening label, not a pseudogene classification.
    Reference biotype is deliberately irrelevant. All tied gene alternatives
    must be partial with one matched reference block, and no selected child
    may be a complete MANE transcript. Keep a mixed tie in the main table.
    """
    parents = [call for group in groups for call in group
               if call.transcript_info.model_type == "full_gene"]
    return bool(parents) and all(
        call.missing_count > 0 and call.found_expected_count == 1 for call in parents
    ) and not any(
        call.transcript_info.model_type == "transcript"
        and call.transcript_info.is_mane and call.missing_count == 0
        for group in groups for call in group
    )


def iter_routed_assignments(calls: Iterable[TranscriptCall], split_processed: bool):
    """Buffer one parent at a time; parent-sorted worker shards stay streaming."""
    assignments = iter_assignment_groups(calls)
    if not split_processed:
        for group in assignments:
            yield group, False
        return
    for gene_index, groups in groupby(assignments, key=lambda group: group[0].gene_index):
        if not gene_index:
            # Transcript-only mode has no parent gene completeness to classify.
            for group in groups:
                yield group, False
            continue
        parent_groups = list(groups)
        processed = is_pseudofragment_assignment(parent_groups)
        for group in parent_groups:
            yield group, processed


def default_pseudofragments_path(output_path: str) -> str:
    suffix = ".transcript_calls.tsv"
    if output_path.endswith(suffix):
        return output_path[:-len(suffix)] + ".pseudofragments.tsv"
    return os.path.splitext(output_path)[0] + ".pseudofragments.tsv"


def write_transcript_calls(
    calls: Iterable[TranscriptCall],
    output_path: str,
    output_format: str,
    no_header: bool,
    truncate_used: bool,
    parse_query_suffix: bool,
    overlap_groups_output: Optional[str] = None,
    overlap_evidence_index: Optional[OverlapEvidenceIndex] = None,
    pseudofragments_output: Optional[str] = None,
) -> int:
    if pseudofragments_output is not None and os.path.realpath(pseudofragments_output) == os.path.realpath(output_path):
        raise ValueError("main and pseudofragments output paths must differ")
    written = 0
    processed_written = 0
    gene_indices: Dict[str, str] = {}
    processed_handle = None
    overlap_handle = None
    overlap_writer = None
    if overlap_groups_output is not None:
        if overlap_evidence_index is None:
            raise ValueError("overlap evidence is required with --overlap-groups-output")
        overlap_handle = open(
            overlap_groups_output, "w", encoding="utf-8", newline=""
        )
        overlap_writer = csv.writer(
            overlap_handle, delimiter="\t", lineterminator="\n"
        )
        overlap_writer.writerow(OVERLAP_GROUP_HEADER)

    try:
        out_handle = open(output_path, "w", encoding="utf-8", newline="")
    except Exception:
        if overlap_handle is not None:
            overlap_handle.close()
        raise
    try:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        processed_writer = None
        if pseudofragments_output is not None:
            processed_handle = open(pseudofragments_output, "w", encoding="utf-8", newline="")
            processed_writer = csv.writer(processed_handle, delimiter="\t", lineterminator="\n")
        if output_format == "extended" and not no_header:
            writer.writerow(OUTPUT_HEADER)
            if processed_writer is not None:
                processed_writer.writerow(OUTPUT_HEADER)
        for group, processed in iter_routed_assignments(calls, processed_writer is not None):
            call = group[0]
            target_writer = processed_writer if processed else writer
            if output_format == "extended":
                row = assignment_to_output_row(group, parse_query_suffix=parse_query_suffix)
                if call.gene_index:
                    row[0] = gene_indices.setdefault(call.gene_index, f"GENE_{len(gene_indices) + 1:06d}")
                target_writer.writerow(row)
            else:
                target_writer.writerow(
                    call_to_legacy_row(
                        call,
                        truncate_used=truncate_used,
                        parse_query_suffix=parse_query_suffix,
                    )
                )
            written += 1
            processed_written += int(processed)
            if overlap_writer is not None and overlap_evidence_index is not None:
                overlap_rows = [
                    call_to_overlap_group_row(
                        alternative,
                        evidence_index=overlap_evidence_index,
                        group_number=written,
                        parse_query_suffix=parse_query_suffix,
                    ) for alternative in group
                ]
                overlap_row = overlap_rows[0]
                if len(group) > 1:
                    overlap_row[3] = str(min(int(r[3]) for r in overlap_rows))
                    overlap_row[4] = str(max(int(r[4]) for r in overlap_rows))
                    for i in (6, 7, 8):
                        overlap_row[i] = ";".join(r[i] for r in overlap_rows)
                    for i in (9, 10, 11, 13):
                        overlap_row[i] = ",".join(sorted({v for r in overlap_rows for v in r[i].split(",") if v}))
                    overlap_row[12] = ",".join(f"{s}-{e}" for s, e in merge_intervals(
                        tuple(map(int, v.split("-"))) for r in overlap_rows for v in r[12].split(",") if v))
                overlap_writer.writerow(overlap_row)
    finally:
        out_handle.close()
        if processed_handle is not None:
            processed_handle.close()
        if overlap_handle is not None:
            overlap_handle.close()
    if pseudofragments_output is not None:
        print(f"Output rows: main={written - processed_written}, pseudofragments={processed_written}", file=sys.stderr)
    return written


def main() -> None:
    parser = argparse.ArgumentParser(description="Call transcripts/genes from exon BLAST alignment table.")
    parser.add_argument("--version", action="version", version=PIPELINE_VERSION)
    parser.add_argument("-i", "--input", required=True, help="exon alignment TSV from align_exon_blastdb_v2.py")
    parser.add_argument("-g", "--gff", required=True, help="GENCODE-style GFF3 annotation")
    parser.add_argument("-o", "--output", required=True, help="output transcript calls TSV")
    parser.add_argument("--pseudofragments-output", default=None,
                        help="candidate gene-like fragment table [derived from --output as .pseudofragments.tsv]")
    parser.add_argument(
        "--overlap-groups-output",
        default=None,
        help=(
            "optional TSV recording the alternative genes/transcripts whose "
            "exon evidence overlaps each final reported call"
        ),
    )
    parser.add_argument("-r", "--region", type=parse_region, default=None, help="optional query/genomic region, 0-based right-open chrom:start-end")
    parser.add_argument("-t", "--truncate", type=parse_region, default=None, help="optional expected-exon truncation region, 0-based right-open chrom:start-end")
    parser.add_argument("--min-exon-coverage", type=float, default=90.0, help="deprecated compatibility option; anchored-query coverage filtering is performed by align_exon_blastdb_v2.py")
    parser.add_argument("--min-identity", type=float, default=95.0, help="deprecated compatibility option; anchored-HSP identity filtering is performed by align_exon_blastdb_v2.py")
    parser.add_argument("--protein-bonus", type=float, default=10.0, help="protein-coding multiplier for real transcript selection only [10]")
    parser.add_argument("--complete-bonus", type=float, default=2.0, help="score multiplier when all expected exons are found [2]")
    parser.add_argument("--mane-bonus", type=float, default=1.0, help="score multiplier for MANE transcripts [1, disabled]")
    parser.add_argument(
        "--prefer-mane",
        action="store_true",
        default=True,
        help="prioritize MANE during transcript selection [default]",
    )
    parser.add_argument(
        "--no-prefer-mane", dest="prefer_mane", action="store_false",
        help="disable MANE priority during transcript selection",
    )
    parser.add_argument(
        "--eligible-exon-info",
        default=None,
        help="optional <db>.exon_info.tsv; only exons present in the database count as expected",
    )
    parser.add_argument(
        "--identical-paralogs", default=None,
        help="identical MANE paralog report [auto-detect indenticalparalogs.tsv beside eligible-exon-info]",
    )
    parser.add_argument("--full-gene-transcripts", dest="full_gene_transcripts", action="store_true",
                        help="lift full-gene exon unions first, then real isoforms within each gene [default]")
    parser.add_argument("--no-full-gene-transcripts", dest="full_gene_transcripts", action="store_false",
                        help="call only annotated transcripts")
    parser.set_defaults(full_gene_transcripts=True)
    parser.add_argument("--full-gene-models-output", default=None,
                        help="optional TSV of synthetic reference exon blocks and their source exons/isoforms")
    parser.add_argument("--shared-exon-genes-output", default=None,
                        help="optional TSV of reference genes merged by >100 bp of shared exonic sequence")
    parser.add_argument(
        "--no-gene-locus-resolution",
        action="store_true",
        help="deprecated compatibility flag; the gene prefilter has been removed",
    )
    parser.add_argument("--max-chains-per-transcript", type=int, default=10, help="maximum non-overlapping chains per transcript/query/strand [10]")
    parser.add_argument("--threads", type=int, default=1, help="parallel query/strand caller processes [1]")
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="explicit parent for temporary worker result shards; overrides --shard-storage",
    )
    parser.add_argument(
        "--shard-storage",
        choices=["auto", "ram", "disk"],
        default="auto",
        help=(
            "temporary transcript-call shard storage: auto prefers /dev/shm, "
            "ram requires /dev/shm, disk prefers $SLURM_TMPDIR [auto]"
        ),
    )
    parser.add_argument(
        "--output-format",
        choices=["extended", "legacy"],
        default="extended",
        help="extended headered output or legacy 16-column headerless output [extended]",
    )
    parser.add_argument("--no-header", action="store_true", help="omit header for extended output")
    parser.add_argument(
        "--query-coordinate-mode",
        choices=["header-suffix", "local"],
        default="header-suffix",
        help="interpret trailing _start_end in query IDs, or keep assembly-contig coordinates local [header-suffix]",
    )
    args = parser.parse_args()
    args.pseudofragments_output = args.pseudofragments_output or default_pseudofragments_path(args.output)
    if os.path.realpath(args.pseudofragments_output) == os.path.realpath(args.output):
        parser.error("main and pseudofragments output paths must differ")
    if args.full_gene_transcripts and args.output_format == "legacy":
        parser.error("two-stage gene/transcript calls require --output-format extended")
    if args.full_gene_models_output and not args.full_gene_transcripts:
        parser.error("--full-gene-models-output requires --full-gene-transcripts")

    parse_query_suffix = args.query_coordinate_mode == "header-suffix"

    stage_start = time.perf_counter()
    annotation = parse_gencode_gff3(args.gff)
    paralog_report = args.identical_paralogs or infer_report(args.eligible_exon_info)
    if paralog_report:
        apply_identical_paralogs(annotation, paralog_report)
        print(f"Applied identical MANE paralog representatives from {paralog_report}", file=sys.stderr)
    shared_groups = apply_shared_exon_genes(annotation)
    print(f"Merged reference genes into {len(shared_groups)} shared-exon gene units", file=sys.stderr)
    if args.shared_exon_genes_output:
        write_shared_exon_report(args.shared_exon_genes_output, shared_groups)
    if args.full_gene_transcripts:
        full_gene_models = add_full_gene_transcripts(annotation)
        model_count = len({row["transcript_id"] for row in full_gene_models})
        print(f"Added {model_count} full-gene models for stage 1; real isoforms are called within each selected gene", file=sys.stderr)
        if args.full_gene_models_output:
            write_full_gene_models(args.full_gene_models_output, full_gene_models)
    print(f"TIMING annotation: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)

    stage_start = time.perf_counter()
    alignments_raw = read_alignments(args.input)
    print(f"Read {len(alignments_raw)} alignment rows", file=sys.stderr)
    print(f"TIMING read_alignments: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)

    stage_start = time.perf_counter()
    alignments = filter_alignments(
        alignments_raw,
        annotation=annotation,
        region=args.region,
        parse_query_suffix=parse_query_suffix,
    )
    print(f"TIMING filter_alignments: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)
    del alignments_raw

    stage_start = time.perf_counter()
    alignments = dedup_same_exon_overlaps(alignments)
    print(f"Kept {len(alignments)} alignments after same-exon overlap deduplication", file=sys.stderr)
    print(f"TIMING deduplicate_alignments: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)

    overlap_evidence_index = None
    if args.overlap_groups_output is not None:
        stage_start = time.perf_counter()
        overlap_evidence_index = OverlapEvidenceIndex(alignments, annotation)
        print(
            f"TIMING overlap_evidence_index: "
            f"{time.perf_counter() - stage_start:.3f} seconds",
            file=sys.stderr,
        )

    stage_start = time.perf_counter()
    eligible_exon_ids = load_eligible_exon_ids(args.eligible_exon_info)
    expected_by_transcript = truncate_expected_exons(
        annotation,
        args.truncate,
        eligible_exon_ids=eligible_exon_ids,
    )
    print(f"TIMING expected_exons: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)

    stage_start = time.perf_counter()
    thread_count = max(1, args.threads)

    if thread_count > 1:
        temp_parent = choose_shard_parent(
            requested_tmp_dir=args.tmp_dir,
            shard_storage=args.shard_storage,
            output_path=args.output,
        )
        with tempfile.TemporaryDirectory(
            prefix="transcript_call_shards_", dir=temp_parent
        ) as shard_dir:
            disk_result = call_transcripts_grouped_to_disk(
                alignments=alignments,
                annotation=annotation,
                expected_by_transcript=expected_by_transcript,
                protein_bonus=args.protein_bonus,
                complete_bonus=args.complete_bonus,
                mane_bonus=args.mane_bonus,
                max_chains_per_transcript=args.max_chains_per_transcript,
                prefer_mane=args.prefer_mane,
                threads=thread_count,
                shard_dir=shard_dir,
                two_scale=args.full_gene_transcripts,
            )
            if disk_result is not None:
                sorted_shard, raw_call_count, selected_count = disk_result
                print(f"Built {raw_call_count} raw transcript chains", file=sys.stderr)
                print(f"Grouping {selected_count} selected alternatives into union assignments", file=sys.stderr)
                print(
                    f"TIMING chain_resolve_and_disk_merge: "
                    f"{time.perf_counter() - stage_start:.3f} seconds",
                    file=sys.stderr,
                )
                # The large parent-side inputs are no longer needed while the
                # disk-backed call stream is converted to TSV.
                del alignments, expected_by_transcript, annotation
                write_start = time.perf_counter()
                written = write_transcript_calls(
                    iter_sorted_call_shard(sorted_shard),
                    output_path=args.output,
                    output_format=args.output_format,
                    no_header=args.no_header,
                    truncate_used=args.truncate is not None,
                    parse_query_suffix=parse_query_suffix,
                    overlap_groups_output=args.overlap_groups_output,
                    overlap_evidence_index=overlap_evidence_index,
                    pseudofragments_output=args.pseudofragments_output,
                )
                print(f"Wrote {written} union assignments", file=sys.stderr)
                print(
                    f"TIMING write_output: {time.perf_counter() - write_start:.3f} seconds",
                    file=sys.stderr,
                )
                return

    calls, raw_call_count = call_transcripts_grouped(
        alignments=alignments,
        annotation=annotation,
        expected_by_transcript=expected_by_transcript,
        protein_bonus=args.protein_bonus,
        complete_bonus=args.complete_bonus,
        mane_bonus=args.mane_bonus,
        max_chains_per_transcript=args.max_chains_per_transcript,
        prefer_mane=args.prefer_mane,
        threads=1,
        two_scale=args.full_gene_transcripts,
    )
    print(f"Built {raw_call_count} raw transcript chains", file=sys.stderr)
    print(f"Selected {len(calls)} alternatives before union grouping", file=sys.stderr)
    print(f"TIMING chain_and_resolve: {time.perf_counter() - stage_start:.3f} seconds", file=sys.stderr)

    write_start = time.perf_counter()
    written = write_transcript_calls(
        calls,
        output_path=args.output,
        output_format=args.output_format,
        no_header=args.no_header,
        truncate_used=args.truncate is not None,
        parse_query_suffix=parse_query_suffix,
        overlap_groups_output=args.overlap_groups_output,
        overlap_evidence_index=overlap_evidence_index,
        pseudofragments_output=args.pseudofragments_output,
    )
    print(f"Wrote {written} union assignments", file=sys.stderr)
    print(f"TIMING write_output: {time.perf_counter() - write_start:.3f} seconds", file=sys.stderr)


if __name__ == "__main__":
    main()
