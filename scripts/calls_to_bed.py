#!/usr/bin/env python3
"""Convert old or new extended transcript calls to sorted BED12+ records."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from typing import Dict, List


EXTRA_COLUMNS = [
    "haplotype",
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_type",
    "ifmane_transcript",
    "call_status",
    "expected_exons",
    "found_exons",
    "missing_exons",
    "fraction_expected_found",
    "mean_identity",
    "mean_exon_coverage",
    "matched_exon_numbers",
    "matched_exon_ids",
    "matched_exon_query_coordinates",
    "matched_exon_reference_coordinates",
    "matched_exon_coverages",
    "matched_exon_identities",
    "tie_group_id",
    "tie_count",
    "assignment_status",
    "pipeline_version",
    "merged_interval_coordinates",
    "merged_interval_gene_scores",
    "model_type",
    "GENE_index",
    "inserted_exons",
    "insertion_penalty",
    "inserted_query_blocks",
    "insertion_run_unique_exons",
]


def split_csv(text: str) -> List[str]:
    return [] if not text else text.rstrip(",").split(",")


def as_int(text: str) -> int:
    return int(float(text))


def exon_query_intervals(row: Dict[str, str]) -> tuple[List[int], List[int]]:
    combined = row.get("exon_query_coordinates", "")
    if combined:
        starts: List[int] = []
        ends: List[int] = []
        for coordinate in (value for gene in combined.split(";") for value in split_csv(gene)):
            try:
                start_text, end_text = coordinate.split("-", 1)
            except ValueError as exc:
                raise ValueError(
                    f"invalid exon query coordinate {coordinate!r}; expected start-end"
                ) from exc
            starts.append(as_int(start_text))
            ends.append(as_int(end_text))
        return starts, ends

    starts = [as_int(x) for x in split_csv(row.get("exon_query_starts", ""))]
    ends = [as_int(x) for x in split_csv(row.get("exon_query_ends", ""))]
    return starts, ends


@dataclass
class BedRecord:
    chrom: str
    start: int
    end: int
    fields: List[str]


def convert(row: Dict[str, str], haplotype: str) -> BedRecord:
    starts, ends = exon_query_intervals(row)
    if not starts or len(starts) != len(ends):
        raise ValueError("missing or inconsistent exon query coordinates")
    reference_text = row.get("exon_reference_coordinates", "")
    if reference_text:
        query_text = row.get("exon_query_coordinates", "")
        reference_groups = [split_csv(gene) for gene in reference_text.split(";")]
        query_groups = [split_csv(gene) for gene in query_text.split(";")] if query_text else [starts]
        if len(reference_groups) != len(query_groups) or any(
            len(refs) != len(query) for refs, query in zip(reference_groups, query_groups)
        ):
            raise ValueError("exon query/reference coordinate arrays have inconsistent lengths or gene groups")

    # Alternative exon structures can overlap. BED geometry is their union;
    # metadata below retains the original per-gene coordinate arrays.
    grouped_exons = ";" in row.get("exon_query_coordinates", "")
    blocks = []
    for start, end in sorted((min(s, e), max(s, e)) for s, e in zip(starts, ends)):
        if grouped_exons and blocks and start <= blocks[-1][1]:
            blocks[-1] = (blocks[-1][0], max(blocks[-1][1], end))
        else:
            blocks.append((start, end))
    chrom_start = min(s for s, _ in blocks)
    chrom_end = max(e for _, e in blocks)
    block_sizes = [e - s for s, e in blocks]
    block_starts = [s - chrom_start for s, _ in blocks]

    fraction = max((float(value) for value in row.get("fraction_expected_found", "").split(";")
                    if value), default=0.0)
    score = max(0, min(1000, round(1000 * fraction)))
    mane = row.get("ifmane_transcript", "0")
    color = "0,102,204" if all(value == "1" for value in mane.split(";")) else "96,96,96"
    name = "|".join(
        [
            row.get("gene_name", "") or row.get("gene_id", "gene"),
            row.get("transcript_id", "transcript"),
            row.get("call_status", "partial"),
        ]
    )

    bed12 = [
        row["query_contig"],
        str(chrom_start),
        str(chrom_end),
        name,
        str(score),
        row["strand"],
        str(chrom_start),
        str(chrom_start),
        color,
        str(len(blocks)),
        ",".join(map(str, block_sizes)) + ",",
        ",".join(map(str, block_starts)) + ",",
    ]
    extras = [
        haplotype,
        row.get("gene_id", ""),
        row.get("gene_name", ""),
        row.get("transcript_id", ""),
        row.get("transcript_type", ""),
        mane,
        row.get("call_status", ""),
        row.get("expected_exons", ""),
        row.get("found_exons", ""),
        row.get("missing_exons", ""),
        row.get("fraction_expected_found", ""),
        row.get("mean_identity", ""),
        row.get("mean_exon_coverage", ""),
        row.get("found_exon_numbers", ""),
        row.get("exon_ids", ""),
        (
            row.get("exon_query_coordinates", "")
            if "exon_query_coordinates" in row
            else ",".join(map(str, starts))
        ),
        (
            reference_text
            if "exon_query_coordinates" in row
            else ",".join(map(str, ends))
        ),
        row.get("exon_coverages", ""),
        row.get("exon_identities", ""),
        row.get("tie_group_id", ""),
        row.get("tie_count", ""),
        row.get("assignment_status", ""),
        row.get("pipeline_version", ""),
        row.get("merged_interval_coordinates", ""),
        row.get("merged_interval_gene_scores", ""),
        row.get("model_type", "transcript"),
        row.get("GENE_index", ""),
        row.get("inserted_exons", "0"),
        row.get("insertion_penalty", "0.000000"),
        row.get("inserted_query_blocks", "0"),
        row.get("insertion_run_unique_exons", ""),
    ]
    return BedRecord(row["query_contig"], chrom_start, chrom_end, bed12 + extras)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--haplotype", required=True)
    args = parser.parse_args()

    records: List[BedRecord] = []
    with open(args.input, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "query_contig", "strand", "transcript_id", "call_status",
        }
        fieldnames = set(reader.fieldnames or [])
        missing = required - fieldnames
        if missing:
            raise SystemExit(f"ERROR: call table is missing columns: {', '.join(sorted(missing))}")
        has_old_coordinates = {
            "exon_query_starts", "exon_query_ends"
        } <= fieldnames
        has_new_coordinates = "exon_query_coordinates" in fieldnames
        if not has_old_coordinates and not has_new_coordinates:
            raise SystemExit(
                "ERROR: call table must contain exon_query_coordinates or both "
                "exon_query_starts and exon_query_ends"
            )
        for line_number, row in enumerate(reader, 2):
            try:
                records.append(convert(row, args.haplotype))
            except Exception as exc:
                raise SystemExit(f"ERROR: {args.input}:{line_number}: {exc}") from exc

    records.sort(key=lambda r: (r.chrom, r.start, r.end, r.fields[3]))
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for record in records:
            writer.writerow(record.fields)


if __name__ == "__main__":
    main()
