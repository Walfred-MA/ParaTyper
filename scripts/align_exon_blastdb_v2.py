#!/usr/bin/env python3
"""
align_exon_blastdb.polished.py

Align reference exons to an assembly and write an exon-alignment TSV for the
transcript caller.

The recommended ``--exons-as-query`` mode builds a temporary BLAST database
from the assembly, uses the representative exon FASTA as the BLAST query, and
parses explicit tabular coordinates.  This makes ``-max_target_seqs`` apply to
assembly loci per exon instead of limiting the total exon hits returned for a
large assembly contig.  The temporary database is placed in ``$SLURM_TMPDIR``
when available and is removed after the job.

The original assembly-query/SAM mode remains available for backward
compatibility.

The default output is a headered, explicit table.  BLAST selection uses the
complete anchored query:

    anchored_query_coverage >= 90%
    anchored_percent_identity > 95%

The selection values are not written to the result.  Output coordinates and
statistics describe only the core annotated exon after gap-aware removal of the
anchors.

A legacy 13-column output compatible with the original caller is available with
--output-format legacy.  In that legacy table, column 9 is aligned exon bases,
so the old `col9 >= 0.9 * exon_length` filter means exon coverage >= 90%.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from contextlib import closing
from dataclasses import dataclass, replace
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
BL_ORD_RE = re.compile(r"BL_ORD_ID(?::|\|)(\d+)")
DEFAULT_BLAST_QUERY_BATCH_BYTES = 100_000_000


def open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def strip_version(identifier: str) -> str:
    if "." not in identifier:
        return identifier
    head, tail = identifier.rsplit(".", 1)
    if tail.isdigit():
        return head
    return identifier


def read_fasta_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    name: Optional[str] = None
    length = 0
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].strip().split()[0]
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            lengths[name] = length
    return lengths


def load_seq_map(seq_path: str) -> Dict[int, str]:
    """Map zero-based BL_ORD_ID to sequence name from <db>.seq."""
    seq_map: Dict[int, str] = {}
    with open(seq_path, "r", encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            name = line.strip().split()[0]
            if name:
                seq_map[i] = name
    return seq_map


@dataclass
class ExonMeta:
    exon_id_full: str
    exon_id: str
    length: int
    coordinates: str = ""
    gene_name: str = ""
    gene_id: str = ""
    transcript_id: str = ""
    ifproteincoding: str = ""
    ifmane_transcript: str = ""
    chrom: str = ""
    start0: int = 0
    end0: int = 0
    strand: str = "+"
    left_anchor_length: int = 0
    right_anchor_length: int = 0
    anchored_length: int = 0


@dataclass
class ExonAlias:
    exon_id_full: str
    exon_id: str
    chrom: str
    start0: int
    end0: int
    strand: str
    length: int
    merge_reason: str = ""
    left_anchor_length: int = 0
    right_anchor_length: int = 0
    anchored_length: int = 0


def load_exon_info(info_path: str) -> Dict[str, ExonMeta]:
    """Load <db>.exon_info.tsv written by the build script."""
    meta: Dict[str, ExonMeta] = {}
    with open(info_path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for required in ["exon_id_full", "exon_id", "length"]:
            if required not in index:
                raise SystemExit(f"ERROR: {info_path} is missing required column {required!r}")
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            try:
                length = int(parts[index["length"]])
            except ValueError:
                continue
            exon_id_full = parts[index["exon_id_full"]]
            rec = ExonMeta(
                exon_id_full=exon_id_full,
                exon_id=parts[index["exon_id"]],
                length=length,
                coordinates=parts[index.get("coordinates", -1)] if "coordinates" in index else "",
                gene_name=parts[index.get("gene_name", -1)] if "gene_name" in index else "",
                gene_id=parts[index.get("gene_id", -1)] if "gene_id" in index else "",
                transcript_id=parts[index.get("transcript_id", -1)] if "transcript_id" in index else "",
                ifproteincoding=parts[index.get("ifproteincoding", -1)] if "ifproteincoding" in index else "",
                ifmane_transcript=parts[index.get("ifmane_transcript", -1)] if "ifmane_transcript" in index else "",
                chrom=parts[index.get("chrom", -1)] if "chrom" in index else "",
                start0=int(parts[index["start0"]]) if "start0" in index and parts[index["start0"]] else 0,
                end0=int(parts[index["end0"]]) if "end0" in index and parts[index["end0"]] else length,
                strand=parts[index.get("strand", -1)] if "strand" in index else "+",
                left_anchor_length=(
                    int(parts[index["left_anchor_length"]])
                    if "left_anchor_length" in index and parts[index["left_anchor_length"]]
                    else 0
                ),
                right_anchor_length=(
                    int(parts[index["right_anchor_length"]])
                    if "right_anchor_length" in index and parts[index["right_anchor_length"]]
                    else 0
                ),
                anchored_length=(
                    int(parts[index["anchored_length"]])
                    if "anchored_length" in index and parts[index["anchored_length"]]
                    else length
                ),
            )
            meta[exon_id_full] = rec
            # Also allow version-stripped lookup for custom DBs or fixed SAMs.
            meta.setdefault(rec.exon_id, rec)
    return meta


def load_exon_aliases(alias_path: str) -> Dict[str, List[ExonAlias]]:
    aliases: Dict[str, List[ExonAlias]] = defaultdict(list)
    with open(alias_path, "r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        required = [
            "blast_exon_id_full", "exon_id_full", "exon_id", "chrom",
            "start0", "end0", "strand", "length",
        ]
        for name in required:
            if name not in index:
                raise SystemExit(f"ERROR: {alias_path} is missing required column {name!r}")
        seen = set()
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            blast_id = parts[index["blast_exon_id_full"]]
            alias = ExonAlias(
                exon_id_full=parts[index["exon_id_full"]],
                exon_id=parts[index["exon_id"]],
                chrom=parts[index["chrom"]],
                start0=int(parts[index["start0"]]),
                end0=int(parts[index["end0"]]),
                strand=parts[index["strand"]],
                length=int(parts[index["length"]]),
                merge_reason=parts[index["merge_reason"]] if "merge_reason" in index else "",
                left_anchor_length=(
                    int(parts[index["left_anchor_length"]])
                    if "left_anchor_length" in index and parts[index["left_anchor_length"]]
                    else 0
                ),
                right_anchor_length=(
                    int(parts[index["right_anchor_length"]])
                    if "right_anchor_length" in index and parts[index["right_anchor_length"]]
                    else 0
                ),
                anchored_length=(
                    int(parts[index["anchored_length"]])
                    if "anchored_length" in index and parts[index["anchored_length"]]
                    else int(parts[index["length"]])
                ),
            )
            key = (blast_id, alias.exon_id_full, alias.chrom, alias.start0, alias.end0, alias.strand)
            if key in seen:
                continue
            seen.add(key)
            aliases[blast_id].append(alias)
            aliases.setdefault(strip_version(blast_id), aliases[blast_id])
    return dict(aliases)


def replace_blast_ord_id(text: str, seq_map: Dict[int, str]) -> str:
    # Depending on BLAST version/output format, local database identifiers are
    # rendered as BL_ORD_ID:0 or gnl|BL_ORD_ID|0.  The whole field represents
    # that ordinal, so return the original FASTA name rather than retaining a
    # synthetic ``gnl|`` prefix.
    match = BL_ORD_RE.search(text)
    if match:
        return seq_map.get(int(match.group(1)), text)
    return text


def parse_tags(fields: List[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for field in fields[11:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    if cigar == "*":
        return []
    ops = [(int(n), op) for n, op in CIGAR_RE.findall(cigar)]
    if not ops:
        raise ValueError(f"could not parse CIGAR: {cigar}")
    return ops


def leading_clip(ops: List[Tuple[int, str]]) -> int:
    total = 0
    for n, op in ops:
        if op in {"S", "H"}:
            total += n
        else:
            break
    return total


def trailing_clip(ops: List[Tuple[int, str]]) -> int:
    total = 0
    for n, op in reversed(ops):
        if op in {"S", "H"}:
            total += n
        else:
            break
    return total


def query_aligned_length(ops: List[Tuple[int, str]]) -> int:
    return sum(n for n, op in ops if op in {"M", "I", "=", "X"})


def target_consumed_length(ops: List[Tuple[int, str]]) -> int:
    return sum(n for n, op in ops if op in {"M", "D", "N", "=", "X"})


def alignment_columns(ops: List[Tuple[int, str]]) -> int:
    return sum(n for n, op in ops if op in {"M", "I", "D", "=", "X"})


def get_query_interval(qname: str, flag: int, ops: List[Tuple[int, str]], query_lengths: Dict[str, int]) -> Tuple[int, int]:
    q_aln = query_aligned_length(ops)
    left = leading_clip(ops)
    right = trailing_clip(ops)
    qlen = query_lengths.get(qname, left + q_aln + right)

    if flag & 16:
        # In SAM, leading soft clip is on the reverse-complemented query.
        qstart = right
        qend = qlen - left
    else:
        qstart = left
        qend = left + q_aln
    if qstart > qend:
        qstart, qend = qend, qstart
    return qstart, qend


def parse_float_tag(tags: Dict[str, str], key: str, default: Optional[float] = None) -> Optional[float]:
    value = tags.get(key)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def parse_int_tag(tags: Dict[str, str], key: str, default: Optional[int] = None) -> Optional[int]:
    value = tags.get(key)
    if value is None:
        return default
    try:
        return int(float(value))
    except ValueError:
        return default


@dataclass
class AlignmentRow:
    query_id: str
    percent_identity: float
    query_start: int
    query_end: int
    strand: str
    exon_id: str
    exon_id_full: str
    exon_length: int
    exon_start: int
    exon_end: int
    aligned_exon_bases: int
    exon_coverage: float
    identical_bases: int
    AS: float
    NM: int
    evalue: float
    cigar: str
    raw_target: str
    selection_percent_identity: float
    selection_coverage: float
    selection_AS: float

    def extended_fields(self) -> List[str]:
        return [
            self.query_id,
            f"{self.percent_identity:.6f}",
            str(self.query_start),
            str(self.query_end),
            self.strand,
            self.exon_id,
            self.exon_id_full,
            str(self.exon_length),
            str(self.exon_start),
            str(self.exon_end),
            str(self.aligned_exon_bases),
            f"{self.exon_coverage:.6f}",
            str(self.identical_bases),
            f"{self.AS:.6f}",
            str(self.NM),
            f"{self.evalue:.6g}",
            self.cigar,
            self.raw_target,
        ]

    def legacy_fields(self) -> List[str]:
        # Original script uses columns 0,2,3,4,5,6,7,8,9,12.
        return [
            self.query_id,                      # 0
            f"{self.percent_identity:.6f}",    # 1
            str(self.query_start),              # 2
            str(self.query_end),                # 3
            self.strand,                        # 4
            self.exon_id,                       # 5, normalized ENSE id
            str(self.exon_length),              # 6
            str(self.exon_start),               # 7
            str(self.exon_end),                 # 8
            str(self.aligned_exon_bases),       # 9, exon coverage bases
            f"{self.AS:.6f}",                  # 10
            f"{self.exon_coverage:.6f}",       # 11
            f"{self.AS:.6f}",                  # 12 placeholder expected by original
        ]


EXTENDED_HEADER = [
    "query_id",
    "percent_identity",
    "query_start",
    "query_end",
    "strand",
    "exon_id",
    "exon_id_full",
    "exon_length",
    "exon_start",
    "exon_end",
    "aligned_exon_bases",
    "exon_coverage",
    "identical_bases",
    "AS",
    "NM",
    "evalue",
    "cigar",
    "raw_target",
]


# BLAST tabular fields used by --exons-as-query.  qseqid is the representative
# reference exon and sseqid is the assembly contig.
BLAST_TABULAR_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gaps",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score",
    "qlen", "slen", "qseq", "sseq",
]

LEGACY_BLAST_TABULAR_FIELD_COUNT = len(BLAST_TABULAR_FIELDS) - 2


def aligned_strings_to_cigar(qseq: str, sseq: str) -> str:
    """Return a compact exon-core CIGAR from two gapped BLAST strings."""
    runs: List[Tuple[int, str]] = []
    for qbase, sbase in zip(qseq, sseq):
        if qbase == "-":
            op = "I"
        elif sbase == "-":
            op = "D"
        elif qbase.upper() == sbase.upper():
            op = "="
        else:
            op = "X"
        if runs and runs[-1][1] == op:
            runs[-1] = (runs[-1][0] + 1, op)
        else:
            runs.append((1, op))
    return "".join(f"{length}{op}" for length, op in runs)


def project_anchored_hsp_to_exon(
    qseq: str,
    sseq: str,
    qstart1: int,
    qend1: int,
    sstart1: int,
    send1: int,
    core_start0: int,
    core_end0: int,
) -> Optional[Tuple[int, int, int, int, int, float, int, float, int, str]]:
    """Project an anchored BLAST HSP onto its core exon, preserving gaps."""
    if len(qseq) != len(sseq) or not qseq:
        return None

    qstep = 1 if qend1 >= qstart1 else -1
    sstep = 1 if send1 >= sstart1 else -1
    qnext = qstart1 - 1
    snext = sstart1 - 1
    qcoords: List[Optional[int]] = []
    scoords: List[Optional[int]] = []

    for qbase, sbase in zip(qseq, sseq):
        if qbase == "-":
            qcoords.append(None)
        else:
            qcoords.append(qnext)
            qnext += qstep
        if sbase == "-":
            scoords.append(None)
        else:
            scoords.append(snext)
            snext += sstep

    core_columns = [
        i for i, coord in enumerate(qcoords)
        if coord is not None and core_start0 <= coord < core_end0
    ]
    if not core_columns:
        return None
    first_col = core_columns[0]
    last_col = core_columns[-1] + 1
    core_qseq = qseq[first_col:last_col]
    core_sseq = sseq[first_col:last_col]
    core_qcoords = [
        coord for coord in qcoords[first_col:last_col]
        if coord is not None and core_start0 <= coord < core_end0
    ]
    core_scoords = [coord for coord in scoords[first_col:last_col] if coord is not None]
    if not core_qcoords or not core_scoords:
        return None

    identical = sum(
        1
        for qbase, sbase in zip(core_qseq, core_sseq)
        if qbase != "-" and sbase != "-" and qbase.upper() == sbase.upper()
    )
    alignment_columns_count = len(core_qseq)
    percent_identity = (
        100.0 * identical / alignment_columns_count
        if alignment_columns_count else 0.0
    )
    nm = alignment_columns_count - identical
    # This is an exon-only mismatch-adjusted score.  The original full-HSP
    # BLAST score remains available internally for candidate selection.
    core_score = float(identical - 3 * nm)
    exon_start = min(core_qcoords) - core_start0
    exon_end = max(core_qcoords) + 1 - core_start0
    assembly_start = min(core_scoords)
    assembly_end = max(core_scoords) + 1
    aligned_exon_bases = len(core_qcoords)
    cigar = aligned_strings_to_cigar(core_qseq, core_sseq)
    return (
        exon_start,
        exon_end,
        assembly_start,
        assembly_end,
        aligned_exon_bases,
        percent_identity,
        identical,
        core_score,
        nm,
        cigar,
    )


def sam_line_to_alignment(
    line: str,
    seq_map: Dict[int, str],
    exon_meta: Dict[str, ExonMeta],
    query_lengths: Dict[str, int],
) -> Optional[AlignmentRow]:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return None
    try:
        flag = int(fields[1])
    except ValueError:
        return None
    if flag & 4:
        return None

    qname = fields[0]
    raw_rname = fields[2]
    rname = replace_blast_ord_id(raw_rname, seq_map)
    exon_id_full = rname.split("\t", 1)[0].split()[0]
    exon_id = strip_version(exon_id_full)
    meta = exon_meta.get(exon_id_full) or exon_meta.get(exon_id)
    if meta is None:
        # Custom DB fallback.  We can still parse, but exon length may be a lower bound.
        meta = ExonMeta(exon_id_full=exon_id_full, exon_id=exon_id, length=0)

    try:
        pos1 = int(fields[3])
    except ValueError:
        return None
    cigar = fields[5]
    if cigar == "*":
        return None

    tags = parse_tags(fields)
    ops = parse_cigar(cigar)
    qstart, qend = get_query_interval(qname, flag, ops, query_lengths)
    t_len = target_consumed_length(ops)
    aln_cols = alignment_columns(ops)
    sstart0 = pos1 - 1
    send0 = sstart0 + t_len
    strand = "-" if (flag & 16) else "+"

    as_score = parse_float_tag(tags, "AS", 0.0) or 0.0
    nm = parse_int_tag(tags, "NM", 0) or 0
    evalue = parse_float_tag(tags, "EV", 0.0) or 0.0
    pi = parse_float_tag(tags, "PI", None)
    if pi is None:
        if aln_cols > 0:
            pi = max(0.0, 100.0 * (aln_cols - nm) / aln_cols)
        else:
            pi = 0.0

    exon_len = meta.length if meta.length > 0 else max(send0, t_len)
    exon_coverage = 100.0 * t_len / exon_len if exon_len > 0 else 0.0
    identical_bases = int(round((pi / 100.0) * aln_cols)) if aln_cols > 0 else 0

    return AlignmentRow(
        query_id=qname,
        percent_identity=pi,
        query_start=qstart,
        query_end=qend,
        strand=strand,
        exon_id=meta.exon_id or exon_id,
        exon_id_full=meta.exon_id_full or exon_id_full,
        exon_length=exon_len,
        exon_start=sstart0,
        exon_end=send0,
        aligned_exon_bases=t_len,
        exon_coverage=exon_coverage,
        identical_bases=identical_bases,
        AS=as_score,
        NM=nm,
        evalue=evalue,
        cigar=cigar,
        raw_target=raw_rname,
        selection_percent_identity=pi,
        selection_coverage=exon_coverage,
        selection_AS=as_score,
    )


def tabular_line_to_alignment(
    line: str,
    assembly_seq_map: Dict[int, str],
    exon_meta: Dict[str, ExonMeta],
) -> Optional[AlignmentRow]:
    """Convert and trim one anchored-exon-query BLAST outfmt-6 HSP."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) not in {LEGACY_BLAST_TABULAR_FIELD_COUNT, len(BLAST_TABULAR_FIELDS)}:
        return None

    exon_id_full = fields[0].split()[0]
    exon_id = strip_version(exon_id_full)
    meta = exon_meta.get(exon_id_full) or exon_meta.get(exon_id)
    if meta is None:
        return None

    raw_subject = fields[1]
    assembly_id = replace_blast_ord_id(raw_subject, assembly_seq_map).split()[0]
    try:
        percent_identity = float(fields[2])
        alignment_length = int(fields[3])
        mismatches = int(fields[4])
        gaps = int(fields[5])
        qstart1 = int(fields[6])
        qend1 = int(fields[7])
        sstart1 = int(fields[8])
        send1 = int(fields[9])
        evalue = float(fields[10])
        raw_score = float(fields[12])
        qlen = int(fields[13])
    except ValueError:
        return None

    exon_length = meta.length if meta.length > 0 else qlen
    anchored_length = (
        meta.anchored_length
        if meta.anchored_length > 0
        else exon_length + meta.left_anchor_length + meta.right_anchor_length
    )
    selection_aligned_query_bases = abs(qend1 - qstart1) + 1
    selection_coverage = (
        100.0 * selection_aligned_query_bases / anchored_length
        if anchored_length > 0 else 0.0
    )

    if len(fields) == LEGACY_BLAST_TABULAR_FIELD_COUNT:
        if meta.left_anchor_length or meta.right_anchor_length:
            raise SystemExit(
                "ERROR: anchored exon databases require qseq and sseq in the "
                "BLAST tabular input; regenerate it with the current "
                "align_exon_blastdb_v2.py"
            )
        # Backward compatibility for unanchored saved BLAST tables.
        exon_start0 = min(qstart1, qend1) - 1
        exon_end0 = max(qstart1, qend1)
        assembly_start0 = min(sstart1, send1) - 1
        assembly_end0 = max(sstart1, send1)
        aligned_exon_bases = max(0, exon_end0 - exon_start0)
        exon_coverage = (
            100.0 * aligned_exon_bases / exon_length if exon_length > 0 else 0.0
        )
        identical_bases = int(round(percent_identity * alignment_length / 100.0))
        core_percent_identity = percent_identity
        core_score = raw_score
        core_nm = mismatches + gaps
        core_cigar = "blast-tabular"
    else:
        qseq = fields[15]
        sseq = fields[16]
        core_start0 = meta.left_anchor_length
        core_end0 = core_start0 + exon_length
        projection = project_anchored_hsp_to_exon(
            qseq=qseq,
            sseq=sseq,
            qstart1=qstart1,
            qend1=qend1,
            sstart1=sstart1,
            send1=send1,
            core_start0=core_start0,
            core_end0=core_end0,
        )
        if projection is None:
            return None
        (
            exon_start0,
            exon_end0,
            assembly_start0,
            assembly_end0,
            aligned_exon_bases,
            core_percent_identity,
            identical_bases,
            core_score,
            core_nm,
            core_cigar,
        ) = projection
        exon_coverage = (
            100.0 * aligned_exon_bases / exon_length if exon_length > 0 else 0.0
        )

    return AlignmentRow(
        query_id=assembly_id,
        percent_identity=core_percent_identity,
        query_start=assembly_start0,
        query_end=assembly_end0,
        # The anchored exon FASTA is transcript-oriented.  Equal query/subject
        # coordinate directions place the transcript on the assembly plus strand.
        strand="+" if (qend1 >= qstart1) == (send1 >= sstart1) else "-",
        exon_id=meta.exon_id or exon_id,
        exon_id_full=meta.exon_id_full or exon_id_full,
        exon_length=exon_length,
        exon_start=exon_start0,
        exon_end=exon_end0,
        aligned_exon_bases=aligned_exon_bases,
        exon_coverage=exon_coverage,
        identical_bases=identical_bases,
        AS=core_score,
        NM=core_nm,
        evalue=evalue,
        cigar=core_cigar,
        raw_target=raw_subject,
        selection_percent_identity=percent_identity,
        selection_coverage=selection_coverage,
        selection_AS=raw_score,
    )


def expand_alignment_aliases(
    row: AlignmentRow,
    representative: Optional[ExonMeta],
    alias_map: Dict[str, List[ExonAlias]],
) -> List[AlignmentRow]:
    aliases = alias_map.get(row.exon_id_full) or alias_map.get(row.exon_id)
    if not aliases or representative is None:
        return [row]

    expanded: List[AlignmentRow] = []
    for alias in aliases:
        aligned = 0
        alias_start = 0
        alias_end = 0
        same_locus = (
            representative.chrom
            and representative.chrom == alias.chrom
            and min(representative.end0, alias.end0) > max(representative.start0, alias.start0)
        )
        if same_locus:
            if representative.strand == "-":
                genomic_start = representative.end0 - row.exon_end
                genomic_end = representative.end0 - row.exon_start
            else:
                genomic_start = representative.start0 + row.exon_start
                genomic_end = representative.start0 + row.exon_end
            overlap_start = max(genomic_start, alias.start0)
            overlap_end = min(genomic_end, alias.end0)
            aligned = max(0, overlap_end - overlap_start)
            if aligned == 0:
                continue
            if alias.strand == "-":
                alias_start = alias.end0 - overlap_end
                alias_end = alias.end0 - overlap_start
            else:
                alias_start = overlap_start - alias.start0
                alias_end = overlap_end - alias.start0
        else:
            # Exact-sequence aliases can occur at different genomic loci.
            alias_start = min(row.exon_start, alias.length)
            alias_end = min(row.exon_end, alias.length)
            aligned = max(0, alias_end - alias_start)
            if aligned == 0:
                aligned = min(row.aligned_exon_bases, alias.length)
                alias_start = 0
                alias_end = aligned

        coverage = 100.0 * aligned / alias.length if alias.length else 0.0
        identical = (
            int(round(row.identical_bases * aligned / row.aligned_exon_bases))
            if row.aligned_exon_bases > 0 else 0
        )
        alias_alignment_strand = row.strand
        if (
            same_locus
            and representative.strand != alias.strand
            and alias.merge_reason != "exact_sequence"
        ):
            alias_alignment_strand = "-" if row.strand == "+" else "+"
        expanded.append(
            replace(
                row,
                strand=alias_alignment_strand,
                exon_id=alias.exon_id,
                exon_id_full=alias.exon_id_full,
                exon_length=alias.length,
                exon_start=alias_start,
                exon_end=alias_end,
                aligned_exon_bases=aligned,
                exon_coverage=coverage,
                identical_bases=identical,
            )
        )
    return expanded or [row]


def blast_command(args: argparse.Namespace) -> List[str]:
    cmd = [
        args.blastn,
        "-task",
        "megablast",
        "-query",
        args.query,
        "-db",
        args.db,
        "-outfmt",
        "17",
        "-word_size",
        str(args.word_size),
        "-num_threads",
        str(args.threads),
        "-evalue",
        args.evalue,
        "-dust",
        "yes",
        "-lcase_masking",
        "-perc_identity",
        str(args.blast_perc_identity),
    ]
    if args.qcov_hsp_perc is not None:
        cmd.extend(["-qcov_hsp_perc", str(args.qcov_hsp_perc)])
    cmd.extend(
        [
            "-max_target_seqs",
            str(args.max_target_seqs),
            "-out",
            "-",
        ]
    )
    return cmd


def exon_query_blast_command(
    args: argparse.Namespace, assembly_db: str, exon_fasta: Optional[str] = None
) -> List[str]:
    exon_fasta = exon_fasta or args.exon_fasta or f"{args.db}.exons.fa"
    cmd = [
        args.blastn,
        "-task", "megablast",
        "-query", exon_fasta,
        "-db", assembly_db,
        "-outfmt", "6 " + " ".join(BLAST_TABULAR_FIELDS),
        "-word_size", str(args.word_size),
        "-num_threads", str(args.threads),
        "-evalue", args.evalue,
        "-dust", "yes",
        "-lcase_masking",
        "-perc_identity", str(args.blast_perc_identity),
    ]
    if args.qcov_hsp_perc is not None:
        cmd.extend(["-qcov_hsp_perc", str(args.qcov_hsp_perc)])
    cmd.extend(["-max_target_seqs", str(args.max_target_seqs), "-out", "-"])
    return cmd


def assembly_ordinal_map(assembly_fasta: str) -> Dict[int, str]:
    """Return the FASTA record order used by makeblastdb BL_ORD_ID values."""
    names: Dict[int, str] = {}
    seen = set()
    with open_text(assembly_fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].strip().split()[0]
            if not name:
                raise SystemExit(f"ERROR: empty FASTA record name in {assembly_fasta}")
            if name in seen:
                raise SystemExit(
                    f"ERROR: duplicate FASTA record name {name!r} in {assembly_fasta}; "
                    "assembly contig names must be unique"
                )
            seen.add(name)
            names[len(names)] = name
    return names


def choose_temp_parent(args: argparse.Namespace) -> Optional[str]:
    requested = args.tmp_dir or os.environ.get("SLURM_TMPDIR")
    if requested:
        requested = os.path.abspath(os.path.expanduser(requested))
        if not os.path.isdir(requested):
            raise SystemExit(f"ERROR: temporary directory does not exist: {requested}")
        return requested
    output_parent = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_parent, exist_ok=True)
    return output_parent


def iter_query_records(path: str) -> Iterator[Tuple[str, str]]:
    """Read one FASTA record at a time, preserving identifiers and masking."""
    header: Optional[str] = None
    pieces: List[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    if not pieces:
                        raise ValueError(f"Empty exon query: {header}")
                    yield header, "".join(pieces)
                header = line.rstrip("\r\n")
                if not header[1:].strip():
                    raise ValueError(f"Empty FASTA header in {path}")
                pieces = []
            elif line.strip():
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                pieces.append("".join(line.split()))
        if header is not None:
            if not pieces:
                raise ValueError(f"Empty exon query: {header}")
            yield header, "".join(pieces)


def iter_exon_query_batches(
    exon_fasta: str, work_dir: str, batch_bytes: int
) -> Iterator[Tuple[str, int, int]]:
    """Yield one on-disk batch at a time; never split an exon across batches.

    The scratch FASTA is reused only after the consumer finishes a BLAST run.
    Only the exon ID is needed in the scratch header; the original descriptions
    and metadata remain in the reference files. An oversized single record is
    rejected rather than split or silently allowed to exceed the byte limit.
    Zero retains the original, unbatched query path.
    """
    if batch_bytes < 0:
        raise ValueError("BLAST query batch bytes must be nonnegative")
    if batch_bytes == 0:
        yield exon_fasta, 0, 0
        return

    batch_path = os.path.join(work_dir, "exon_query_batch.fa")
    count = bases = size = 0
    out = open(batch_path, "wb")
    try:
        with closing(iter_query_records(exon_fasta)) as records:
            for header, sequence in records:
                identifier = header[1:].split()[0]
                record = f">{identifier}\n{sequence}\n".encode("utf-8")
                if len(record) > batch_bytes:
                    raise ValueError(
                        f"Exon query {identifier} needs {len(record)} FASTA bytes, "
                        f"exceeding --blast-query-batch-bytes {batch_bytes}; "
                        "increase the limit to keep this exon intact"
                    )
                if count and size + len(record) > batch_bytes:
                    out.close()
                    yield batch_path, count, bases
                    out = open(batch_path, "wb")
                    count = bases = size = 0
                out.write(record)
                count += 1
                bases += len(sequence)
                size += len(record)
        if not count:
            raise ValueError(f"No exon queries found in {exon_fasta}")
        out.close()
        yield batch_path, count, bases
    finally:
        out.close()
        if os.path.exists(batch_path):
            os.unlink(batch_path)


def iter_exon_query_lines(args: argparse.Namespace) -> Iterator[str]:
    """Yield tabular HSPs with exons as query and the assembly as target."""
    if args.blast_tabular:
        with open_text(args.blast_tabular) as handle:
            yield from handle
        return

    exon_fasta = args.exon_fasta or f"{args.db}.exons.fa"
    if not os.path.isfile(exon_fasta):
        raise SystemExit(f"ERROR: reference exon FASTA does not exist: {exon_fasta}")

    temp_parent = choose_temp_parent(args)
    with tempfile.TemporaryDirectory(
        prefix="exon_blast_assembly_", dir=temp_parent
    ) as work_dir:
        assembly_db = os.path.join(work_dir, "assembly")
        make_cmd = [
            args.makeblastdb,
            "-in", args.query,
            "-dbtype", "nucl",
            "-blastdb_version", "5",
            "-out", assembly_db,
        ]
        print(f"Temporary assembly BLAST database: {work_dir}", file=sys.stderr)
        print("Running:", " ".join(make_cmd), file=sys.stderr)
        # Inherit stdout/stderr so makeblastdb diagnostics appear immediately in
        # the SLURM .out/.err stream rather than being hidden by Snakemake.
        subprocess.run(make_cmd, check=True)

        # All batches search the same complete database. Stream their evidence
        # through the same alias expansion/deduplication and call transcripts
        # only after every batch succeeds, preserving cross-gene competition.
        with closing(iter_exon_query_batches(
            exon_fasta, work_dir, args.blast_query_batch_bytes
        )) as batches:
            for batch_number, (query_path, count, bases) in enumerate(batches, 1):
                if args.blast_query_batch_bytes:
                    print(
                        f"BLAST query batch {batch_number}: {count} exons, "
                        f"{bases} bases, {os.path.getsize(query_path)} FASTA bytes "
                        f"(limit {args.blast_query_batch_bytes}), {args.threads} threads",
                        file=sys.stderr,
                    )
                else:
                    print("BLAST query batching disabled", file=sys.stderr)
                cmd = exon_query_blast_command(args, assembly_db, query_path)
                print("Running:", " ".join(cmd), file=sys.stderr)
                proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
                assert proc.stdout is not None
                try:
                    yield from proc.stdout
                    ret = proc.wait()
                    if ret != 0:
                        raise subprocess.CalledProcessError(ret, cmd)
                finally:
                    proc.stdout.close()
                    if proc.poll() is None:
                        proc.terminate()
                        proc.wait()


def iter_sam_lines(args: argparse.Namespace) -> Iterator[str]:
    if args.sam:
        with open_text(args.sam) as handle:
            for line in handle:
                yield line
        return

    cmd = blast_command(args)
    print("Running:", " ".join(cmd), file=sys.stderr)
    # Let BLAST diagnostics flow directly to the SLURM job's stderr.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        yield line
    ret = proc.wait()
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Project reference exons onto an assembly and emit a filtered exon-alignment TSV."
    )
    parser.add_argument("-q", "--query", required=True, help="input assembly FASTA")
    parser.add_argument("-d", "--db", required=True, help="reference exon prefix from build_exon_blastdb_v2.py")
    parser.add_argument("-o", "--output", required=True, help="output TSV")
    parser.add_argument("-t", "--threads", type=int, default=1, help="BLAST threads [1]")
    parser.add_argument(
        "--blast-query-batch-bytes", type=int, default=DEFAULT_BLAST_QUERY_BATCH_BYTES,
        help="maximum query FASTA bytes per sequential BLAST run; whole exons stay intact; 0 disables batching [100000000]",
    )
    parser.add_argument("--min-exon-coverage", type=float, default=90.0, help="minimum anchored-query coverage percentage [90]")
    parser.add_argument("--min-identity", type=float, default=95.0, help="anchored-HSP percent identity must be greater than this value [95]")
    parser.add_argument("--min-as", type=float, default=50.0, help="minimum BLAST raw alignment score; kept if score > this value [50]")
    parser.add_argument("--blast-perc-identity", type=float, default=95.0, help="BLAST -perc_identity [95]")
    parser.add_argument("--qcov-hsp-perc", type=float, default=None, help="optional BLAST -qcov_hsp_perc; omitted by default")
    parser.add_argument(
        "--no-qcov-hsp-perc",
        action="store_true",
        help="omit BLAST -qcov_hsp_perc; explicit exon-coverage filtering is still applied",
    )
    parser.add_argument("--word-size", type=int, default=19, help="BLAST -word_size [19]")
    parser.add_argument("--evalue", default="1e-30", help="BLAST -evalue [1e-30]")
    parser.add_argument("--max-target-seqs", type=int, default=100, help="BLAST -max_target_seqs; assembly loci per exon in recommended mode [100]")
    parser.add_argument("--blastn", default="blastn", help="path to blastn [blastn]")
    parser.add_argument("--makeblastdb", default="makeblastdb", help="path to makeblastdb [makeblastdb]")
    parser.add_argument(
        "--exons-as-query",
        action="store_true",
        help=(
            "recommended mode: query reference exons against a temporary BLAST "
            "database made from the assembly"
        ),
    )
    parser.add_argument(
        "--exon-fasta",
        default=None,
        help="representative reference exon FASTA [default: <db>.exons.fa]",
    )
    parser.add_argument(
        "--tmp-dir",
        default=None,
        help="parent for temporary assembly DB [default: $SLURM_TMPDIR, then output directory]",
    )
    parser.add_argument(
        "--blast-tabular",
        default=None,
        help="parse an existing exon-query BLAST tabular file instead of running BLAST",
    )
    parser.add_argument("--seq", default=None, help="name map file [default: <db>.seq]")
    parser.add_argument("--exon-info", default=None, help="exon metadata file [default: <db>.exon_info.tsv]")
    parser.add_argument(
        "--exon-aliases",
        default=None,
        help="merged-exon alias table [default: <db>.exon_aliases.tsv when present]",
    )
    parser.add_argument("--sam", default=None, help="parse existing BLAST -outfmt 17 SAM instead of running blastn")
    parser.add_argument(
        "--output-format",
        choices=["extended", "legacy"],
        default="extended",
        help="extended headered TSV for the polished caller, or legacy 13-column TSV for the original caller [extended]",
    )
    parser.add_argument("--no-header", action="store_true", help="do not write header for extended output")
    args = parser.parse_args()
    if args.blast_query_batch_bytes < 0:
        parser.error("--blast-query-batch-bytes must be nonnegative")
    if args.no_qcov_hsp_perc:
        args.qcov_hsp_perc = None
    if args.blast_tabular and not args.exons_as_query:
        raise SystemExit("ERROR: --blast-tabular requires --exons-as-query")
    if args.sam and args.exons_as_query:
        raise SystemExit("ERROR: --sam cannot be combined with --exons-as-query")

    output_parent = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_parent, exist_ok=True)

    seq_path = args.seq or f"{args.db}.seq"
    info_path = args.exon_info or f"{args.db}.exon_info.tsv"
    alias_path = args.exon_aliases or f"{args.db}.exon_aliases.tsv"
    exon_meta = load_exon_info(info_path)
    if not exon_meta:
        raise SystemExit(f"ERROR: could not load exon metadata from {info_path}")
    alias_map = load_exon_aliases(alias_path) if os.path.isfile(alias_path) else {}
    if alias_map:
        print(f"Loaded merged-exon aliases from {alias_path}", file=sys.stderr)
    if args.exons_as_query:
        seq_map = assembly_ordinal_map(args.query)
        if not seq_map:
            raise SystemExit(f"ERROR: no assembly FASTA records found in {args.query}")
        query_lengths: Dict[str, int] = {}
        alignment_lines = iter_exon_query_lines(args)
        record_label = "tabular BLAST HSP"
    else:
        if any(
            meta.left_anchor_length or meta.right_anchor_length
            for meta in exon_meta.values()
        ):
            raise SystemExit(
                "ERROR: anchored exon databases require --exons-as-query so "
                "core exon coordinates can be recovered gap-aware"
            )
        seq_map = load_seq_map(seq_path)
        if not seq_map:
            raise SystemExit(f"ERROR: could not load BLAST ordinal map from {seq_path}")
        query_lengths = read_fasta_lengths(args.query)
        if not query_lengths:
            raise SystemExit(f"ERROR: no query FASTA records found in {args.query}")
        alignment_lines = iter_sam_lines(args)
        record_label = "SAM alignment"

    seen = set()
    total = 0
    written = 0
    dropped_identity = 0
    dropped_coverage = 0
    dropped_as = 0

    with open(args.output, "w", encoding="utf-8") as out:
        if args.output_format == "extended" and not args.no_header:
            out.write("\t".join(EXTENDED_HEADER) + "\n")
        for line in alignment_lines:
            if not line.strip() or line.startswith("@"):
                continue
            total += 1
            if args.exons_as_query:
                row = tabular_line_to_alignment(line, seq_map, exon_meta)
            else:
                row = sam_line_to_alignment(line, seq_map, exon_meta, query_lengths)
            if row is None:
                continue
            if row.selection_AS <= args.min_as:
                dropped_as += 1
                continue
            if row.selection_percent_identity <= args.min_identity:
                dropped_identity += 1
                continue
            if row.selection_coverage < args.min_exon_coverage:
                dropped_coverage += 1
                continue
            representative = exon_meta.get(row.exon_id_full) or exon_meta.get(row.exon_id)
            for alias_row in expand_alignment_aliases(row, representative, alias_map):
                # Avoid duplicate rows caused by repeated SAM records or duplicated aliases.
                key = (
                    alias_row.query_id,
                    alias_row.query_start,
                    alias_row.query_end,
                    alias_row.strand,
                    alias_row.exon_id,
                    alias_row.exon_start,
                    alias_row.exon_end,
                    alias_row.cigar,
                )
                if key in seen:
                    continue
                seen.add(key)

                if args.output_format == "extended":
                    out.write("\t".join(alias_row.extended_fields()) + "\n")
                else:
                    out.write("\t".join(alias_row.legacy_fields()) + "\n")
                written += 1

    print(f"Read {total} {record_label} records", file=sys.stderr)
    print(f"Dropped {dropped_as} with anchored-HSP AS <= {args.min_as}", file=sys.stderr)
    print(f"Dropped {dropped_identity} with anchored-HSP identity <= {args.min_identity}", file=sys.stderr)
    print(f"Dropped {dropped_coverage} with anchored-query coverage < {args.min_exon_coverage}", file=sys.stderr)
    print(f"Wrote {written} filtered exon alignments to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
