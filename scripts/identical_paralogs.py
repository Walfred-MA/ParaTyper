"""Shared report format for genes collapsed by identical complete MANE DNA."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Dict, List, Sequence


REPORT_NAME = "indenticalparalogs.tsv"
MERGE_POLICY = "identical_complete_mane_sequence_sets_v1"
DATABASE_FORMAT = "anchored_exon_database_v2"
REPORT_HEADER = [
    "representative_gene_id", "representative_gene_name", "merged_gene_name",
    "gene_ids", "gene_names", "mane_transcript_ids", "mane_sequence_lengths",
    "mane_sequence_sha256",
]


def file_identity(path) -> dict:
    path = Path(path).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def report_digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_report(path: str, rows: Sequence[Dict[str, str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, REPORT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_report(path) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(REPORT_HEADER) <= set(reader.fieldnames or []):
            raise ValueError(f"{path}: incompatible identical-paralog report header")
        rows = list(reader)
    seen = set()
    for row in rows:
        ids = row["gene_ids"].split(";")
        names = row["gene_names"].split(";")
        if (len(ids) < 2 or len(ids) != len(names) or len(set(ids)) != len(ids)
                or not all(ids) or set(ids) & seen
                or row["representative_gene_id"] != ids[0]
                or row["representative_gene_name"] != names[0]
                or row["merged_gene_name"] != names[0] + "merged"):
            raise ValueError(f"{path}: invalid identical-paralog group {row.get('gene_ids')!r}")
        for key in ("mane_transcript_ids", "mane_sequence_lengths", "mane_sequence_sha256"):
            if len(row[key].split(";")) != len(ids):
                raise ValueError(f"{path}: {key} does not match the gene order")
        seen.update(ids)
    return rows


def infer_report(exon_info_path: str | None) -> Path | None:
    if exon_info_path:
        path = Path(exon_info_path).resolve().parent / REPORT_NAME
        if path.is_file():
            return path
    return None
