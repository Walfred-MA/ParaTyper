#!/usr/bin/env python3
"""Run the anchored-exon gene-annotation pipeline on assembly FASTAs.

The shared anchored-exon database is built once and reused.  For each named
assembly, exon alignments are kept below OUTPUT/temp/ and the final transcript
call table is written below OUTPUT/SAMPLE/.  A multi-record FASTA can instead
be aligned and called once, then split by its unique record names.

Version 3 reports exon-supported transcript assignments and preserves ties
after MANE priority. Coordinates describe the retained exon evidence.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import fcntl
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import threading
from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

from identical_paralogs import (DATABASE_FORMAT, MERGE_POLICY, REPORT_NAME,
                                file_identity, read_report, report_digest)
import shared_exon_genes


DATABASE_SUFFIXES = (
    ".exons.fa",
    ".seq",
    ".exon_info.tsv",
    ".exon_aliases.tsv",
    ".manifest.json",
)
SAMPLE_NAME_RE = re.compile(r"[A-Za-z0-9_.-]+")
PRINT_LOCK = threading.Lock()
PIPELINE_VERSION = "3.9.0"


@dataclass(frozen=True)
class AssemblyQuery:
    name: str
    fasta: Path


@dataclass(frozen=True)
class Scripts:
    build: Path
    align: Path
    call: Path


def log(message: str) -> None:
    with PRINT_LOCK:
        print(message, file=sys.stderr, flush=True)


def resolve_existing_file(path_text: str, label: str) -> Path:
    path = Path(path_text).expanduser().resolve()
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} does not exist or is not a file: {path}")
    return path


def read_queries(path: Path) -> List[AssemblyQuery]:
    """Read `sample_name assembly.fasta`, accepting shell-quoted paths."""
    result: List[AssemblyQuery] = []
    seen = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                fields = shlex.split(line)
            except ValueError as exc:
                raise SystemExit(f"ERROR: {path}:{line_number}: {exc}") from exc
            if len(fields) != 2:
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: expected "
                    "'sample_name assembly.fasta'"
                )
            name, fasta_text = fields
            if not SAMPLE_NAME_RE.fullmatch(name):
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: unsafe sample name {name!r}"
                )
            if name in seen:
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: duplicate sample {name!r}"
                )
            fasta = Path(fasta_text).expanduser()
            if not fasta.is_absolute():
                fasta = path.parent / fasta
            fasta = fasta.resolve()
            if not fasta.is_file():
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: assembly FASTA does not exist: "
                    f"{fasta}"
                )
            seen.add(name)
            result.append(AssemblyQuery(name=name, fasta=fasta))
    if not result:
        raise SystemExit(f"ERROR: no assembly FASTAs found in {path}")
    return result


def read_fasta_sample_names(path: Path) -> List[str]:
    """Return unique first-token FASTA record names and validate sequences."""
    names: List[str] = []
    seen = set()
    current_name = ""
    current_has_sequence = False
    saw_content = False
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            saw_content = True
            if line.startswith(">"):
                if current_name and not current_has_sequence:
                    raise SystemExit(
                        f"ERROR: {path}: FASTA record {current_name!r} has no sequence"
                    )
                header = line[1:].strip()
                if not header:
                    raise SystemExit(f"ERROR: {path}:{line_number}: empty FASTA header")
                name = header.split()[0]
                if not SAMPLE_NAME_RE.fullmatch(name):
                    raise SystemExit(
                        f"ERROR: {path}:{line_number}: unsafe FASTA sample name "
                        f"{name!r}; use only letters, numbers, ., _, and -"
                    )
                if name in seen:
                    raise SystemExit(
                        f"ERROR: {path}:{line_number}: duplicate FASTA sample "
                        f"name {name!r}"
                    )
                seen.add(name)
                names.append(name)
                current_name = name
                current_has_sequence = False
                continue
            if not current_name:
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: sequence data appears before "
                    "the first FASTA header"
                )
            if any(base.isspace() for base in line):
                raise SystemExit(
                    f"ERROR: {path}:{line_number}: whitespace inside FASTA sequence"
                )
            current_has_sequence = True
    if not saw_content or not names:
        raise SystemExit(f"ERROR: no FASTA records found in {path}")
    if current_name and not current_has_sequence:
        raise SystemExit(
            f"ERROR: {path}: FASTA record {current_name!r} has no sequence"
        )
    return names


def load_scripts(scripts_dir: Path) -> Scripts:
    paths = Scripts(
        build=scripts_dir / "build_exon_blastdb_v2.py",
        align=scripts_dir / "align_exon_blastdb_v2.py",
        call=scripts_dir / "call_genes_from_exon_alignments_v3.py",
    )
    required = [
        ("database builder", paths.build),
        ("exon aligner", paths.align),
        ("transcript caller", paths.call),
    ]
    for label, path in required:
        if not path.is_file():
            raise SystemExit(f"ERROR: {label} script is missing: {path}")
    return paths


def read_header(path: Path) -> List[str]:
    try:
        with path.open(encoding="utf-8") as handle:
            return handle.readline().rstrip("\n").split("\t")
    except (OSError, UnicodeError):
        return []


def validate_database(prefix: Path) -> Tuple[bool, str]:
    for suffix in DATABASE_SUFFIXES:
        path = Path(str(prefix) + suffix)
        if not path.is_file():
            return False, f"missing {path.name}"
        if path.stat().st_size == 0:
            return False, f"empty {path.name}"

    report = prefix.parent / REPORT_NAME
    try:
        read_report(report)
        manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
        if (manifest.get("format") != DATABASE_FORMAT
                or manifest.get("identical_mane_policy") != MERGE_POLICY
                or manifest.get("identical_paralogs_sha256") != report_digest(report)):
            return False, "missing or stale identical MANE paralog metadata"
    except (OSError, ValueError, KeyError, TypeError):
        return False, "missing or invalid identical MANE paralog report/manifest"

    shared_report = prefix.parent / shared_exon_genes.REPORT_NAME
    try:
        shared_exon_genes.read_report(shared_report)
        if (manifest.get("shared_exon_gene_policy") != shared_exon_genes.MERGE_POLICY
                or manifest.get("shared_exon_genes_sha256") != report_digest(shared_report)):
            return False, "missing or stale shared-exon gene metadata"
    except (OSError, ValueError, KeyError, TypeError):
        return False, "missing or invalid shared-exon gene report/manifest"

    info_path = Path(str(prefix) + ".exon_info.tsv")
    info_header = set(read_header(info_path))
    required_info = {
        "exon_id_full",
        "exon_id",
        "length",
        "left_anchor_length",
        "right_anchor_length",
        "anchored_length",
    }
    missing_info = required_info - info_header
    if missing_info:
        return False, "incompatible exon_info columns: " + ", ".join(
            sorted(missing_info)
        )

    alias_path = Path(str(prefix) + ".exon_aliases.tsv")
    alias_header = set(read_header(alias_path))
    required_aliases = {
        "blast_exon_id_full",
        "exon_id",
        "left_anchor_length",
        "right_anchor_length",
        "anchored_length",
    }
    missing_aliases = required_aliases - alias_header
    if missing_aliases:
        return False, "incompatible exon_aliases columns: " + ", ".join(
            sorted(missing_aliases)
        )

    fasta_path = Path(str(prefix) + ".exons.fa")
    try:
        with fasta_path.open(encoding="utf-8") as handle:
            first_header = next(
                (line.rstrip("\n") for line in handle if line.startswith(">")),
                "",
            )
    except (OSError, UnicodeError):
        first_header = ""
    if not first_header:
        return False, f"no FASTA records in {fasta_path.name}"
    if len(first_header[1:].split("\t")) < 10:
        return False, f"anchored metadata is absent from {fasta_path.name}"
    return True, "complete anchored-exon database"


def run_command(command: Sequence[str], label: str) -> None:
    log(f"[CMD] {label}: {shlex.join(list(command))}")
    subprocess.run(list(command), check=True)


def write_json_atomic(path: Path, value: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def ensure_database(
    prefix: Path,
    scripts: Scripts,
    python: str,
    reference: Path,
    gff3: Path,
    anchor_size: int,
    min_unmasked: int,
    merge_exon_overlap: float,
    force: bool,
) -> bool:
    """Ensure a current database, returning whether it was rebuilt."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(str(prefix) + ".build.lock")
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        valid, reason = validate_database(prefix)
        if valid:
            manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
            expected = {
                "reference_identity": file_identity(reference), "gff3_identity": file_identity(gff3),
                "anchor_size": anchor_size, "min_unmasked": min_unmasked,
                "merge_exon_overlap": merge_exon_overlap,
            }
            if any(manifest.get(key) != value for key, value in expected.items()):
                valid, reason = False, "reference, annotation, or database parameters changed"
        if valid and not force:
            log(f"[LOG] Reusing exon database {prefix} ({reason})")
            return False
        if force:
            log(f"[LOG] Rebuilding exon database {prefix} (--force-rebuild-database)")
        else:
            log(f"[LOG] Building exon database {prefix}: {reason}")

        with tempfile.TemporaryDirectory(
            prefix=f".{prefix.name}.build.",
            dir=prefix.parent,
        ) as temporary_dir:
            temporary_prefix = Path(temporary_dir) / prefix.name
            command = [
                python,
                str(scripts.build),
                "--genome",
                str(reference),
                "--gff3",
                str(gff3),
                "--out",
                str(temporary_prefix),
                "--anchor-size",
                str(anchor_size),
                "--min-unmasked",
                str(min_unmasked),
                "--merge-exon-overlap",
                str(merge_exon_overlap),
                "--no-makeblastdb",
            ]
            run_command(command, "build exon database")
            valid, reason = validate_database(temporary_prefix)
            if not valid:
                raise RuntimeError(
                    f"new exon database failed validation: {reason}"
                )
            os.replace(Path(temporary_dir) / REPORT_NAME, prefix.parent / REPORT_NAME)
            os.replace(Path(temporary_dir) / shared_exon_genes.REPORT_NAME,
                       prefix.parent / shared_exon_genes.REPORT_NAME)
            for suffix in DATABASE_SUFFIXES:
                source = Path(str(temporary_prefix) + suffix)
                destination = Path(str(prefix) + suffix)
                os.replace(source, destination)

        log(f"[DONE] Built exon database {prefix}")
        return True


def valid_call_table(path: Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    header = set(read_header(path))
    required = {
        "query_contig", "transcript_id", "call_status", "tie_group_id",
        "tie_count", "assignment_status", "pipeline_version", "model_type", "GENE_index", "inserted_exons", "insertion_penalty", "insertion_run_unique_exons",
    }
    if not required <= header or header & {
        "transcript_id_full", "gene_id_full", "alternatives_json"
    }:
        return False
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            first = next(csv.DictReader(handle, delimiter="\t"), None)
    except (OSError, UnicodeError, csv.Error):
        return False
    # Header-only tables are reusable only after validating the current schema.
    return first is None or first.get("pipeline_version") == PIPELINE_VERSION


def make_alignment_command(
    query_fasta: Path,
    output_path: Path,
    blast_temp: Path,
    args: argparse.Namespace,
    scripts: Scripts,
    database_prefix: Path,
) -> List[str]:
    return [
        args.python,
        str(scripts.align),
        "--query",
        str(query_fasta),
        "--db",
        str(database_prefix),
        "--exons-as-query",
        "--exon-fasta",
        str(database_prefix) + ".exons.fa",
        "--exon-info",
        str(database_prefix) + ".exon_info.tsv",
        "--exon-aliases",
        str(database_prefix) + ".exon_aliases.tsv",
        "--output",
        str(output_path),
        "--threads",
        str(args.blast_threads),
        "--tmp-dir",
        str(blast_temp),
        "--evalue",
        str(args.evalue),
        "--max-target-seqs",
        str(args.max_target_seqs),
        "--no-qcov-hsp-perc",
        "--min-exon-coverage",
        str(args.min_exon_coverage),
        "--min-identity",
        str(args.min_identity),
        "--blast-perc-identity",
        str(args.min_identity),
        "--min-as",
        str(args.min_alignment_score),
        "--word-size",
        str(args.word_size),
        "--blastn",
        args.blastn,
        "--makeblastdb",
        args.makeblastdb,
    ]


def make_caller_command(
    alignments: Path,
    output_path: Path,
    args: argparse.Namespace,
    scripts: Scripts,
    database_prefix: Path,
    pseudofragments_output: Path | None = None,
) -> List[str]:
    command = [
        args.python,
        str(scripts.call),
        "--input",
        str(alignments),
        "--gff",
        str(args.gff3),
        "--output",
        str(output_path),
        "--eligible-exon-info",
        str(database_prefix) + ".exon_aliases.tsv",
        "--query-coordinate-mode",
        args.query_coordinate_mode,
        "--protein-bonus",
        str(args.protein_bonus),
        "--complete-bonus",
        str(args.complete_bonus),
        "--max-chains-per-transcript",
        str(args.max_chains_per_transcript),
        "--threads",
        str(args.caller_threads),
        "--shard-storage",
        args.caller_shard_storage,
    ]
    command.append("--prefer-mane" if args.prefer_mane else "--no-prefer-mane")
    if pseudofragments_output is not None:
        command.extend(["--pseudofragments-output", str(pseudofragments_output)])
    command.append("--full-gene-transcripts" if getattr(args, "full_gene_transcripts", True)
                   else "--no-full-gene-transcripts")
    report = database_prefix.parent / REPORT_NAME
    if report.is_file():
        command.extend(["--identical-paralogs", str(report)])
    return command


def run_one_sample(
    query: AssemblyQuery,
    args: argparse.Namespace,
    scripts: Scripts,
    database_prefix: Path,
    output_root: Path,
    temp_root: Path,
) -> Path:
    sample_output = output_root / query.name
    sample_temp = temp_root / query.name
    sample_output.mkdir(parents=True, exist_ok=True)
    sample_temp.mkdir(parents=True, exist_ok=True)

    final_calls = sample_output / f"{query.name}.transcript_calls.tsv"
    final_processed = sample_output / f"{query.name}.pseudofragments.tsv"
    if not args.force_samples and all(valid_call_table(path) for path in (final_calls, final_processed)):
        log(f"[SKIP] {query.name}: both final call tables already exist")
        return final_calls

    alignments = sample_temp / f"{query.name}.exon_alignments.tsv"
    alignment_partial = sample_temp / f".{query.name}.exon_alignments.tsv.partial"
    calls_partial = sample_temp / f".{query.name}.transcript_calls.tsv.partial"
    processed_partial = sample_temp / f".{query.name}.pseudofragments.tsv.partial"
    for partial in (alignment_partial, calls_partial, processed_partial):
        if partial.exists():
            partial.unlink()

    # Own a separate scratch folder for this alignment. Remove it, including
    # any remaining assembly database files, before starting transcript calls.
    with tempfile.TemporaryDirectory(prefix="blast_tmp_", dir=sample_temp) as blast_temp:
        align_command = make_alignment_command(
            query_fasta=query.fasta,
            output_path=alignment_partial,
            blast_temp=Path(blast_temp),
            args=args,
            scripts=scripts,
            database_prefix=database_prefix,
        )
        run_command(align_command, f"{query.name} align exons")
    os.replace(alignment_partial, alignments)

    call_command = make_caller_command(
        alignments=alignments,
        output_path=calls_partial,
        args=args,
        scripts=scripts,
        database_prefix=database_prefix,
        pseudofragments_output=processed_partial,
    )
    run_command(call_command, f"{query.name} call transcripts")
    if not all(valid_call_table(path) for path in (calls_partial, processed_partial)):
        raise RuntimeError(
            f"caller did not produce two valid tables for {query.name}"
        )
    os.replace(processed_partial, final_processed)
    os.replace(calls_partial, final_calls)
    log(f"[DONE] {query.name} -> {final_calls}")
    log(f"[DONE] {query.name} pseudofragments -> {final_processed}")
    return final_calls


def split_combined_call_table(
    combined_calls: Path,
    sample_names: Sequence[str],
    output_root: Path,
    max_open_files: int = 64,
    table_suffix: str = "transcript_calls",
) -> List[Path]:
    """Stream one caller table into atomic per-sample output tables."""
    sample_set = set(sample_names)
    partial_paths: Dict[str, Path] = {}
    final_paths: Dict[str, Path] = {}
    open_outputs: "OrderedDict[str, Tuple[object, csv.DictWriter]]" = OrderedDict()

    with combined_calls.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        if "query_id" not in fieldnames:
            raise RuntimeError(
                f"combined caller table is missing query_id: {combined_calls}"
            )

        for sample_name in sample_names:
            sample_output = output_root / sample_name
            sample_output.mkdir(parents=True, exist_ok=True)
            final_path = sample_output / f"{sample_name}.{table_suffix}.tsv"
            partial_path = sample_output / f".{sample_name}.{table_suffix}.tsv.partial"
            if partial_path.exists():
                partial_path.unlink()
            with partial_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
            partial_paths[sample_name] = partial_path
            final_paths[sample_name] = final_path

        def writer_for(sample_name: str) -> csv.DictWriter:
            existing = open_outputs.pop(sample_name, None)
            if existing is not None:
                open_outputs[sample_name] = existing
                return existing[1]
            if len(open_outputs) >= max_open_files:
                _old_name, (old_handle, _old_writer) = open_outputs.popitem(
                    last=False
                )
                old_handle.close()
            handle = partial_paths[sample_name].open(
                "a", encoding="utf-8", newline=""
            )
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
            )
            open_outputs[sample_name] = (handle, writer)
            return writer

        try:
            for line_number, row in enumerate(reader, start=2):
                sample_name = row.get("query_id", "")
                if sample_name not in sample_set:
                    raise RuntimeError(
                        f"{combined_calls}:{line_number}: query_id "
                        f"{sample_name!r} is not a FASTA sample"
                    )
                writer_for(sample_name).writerow(row)
        finally:
            for handle, _writer in open_outputs.values():
                handle.close()
            open_outputs.clear()

    results: List[Path] = []
    for sample_name in sample_names:
        partial_path = partial_paths[sample_name]
        final_path = final_paths[sample_name]
        if not valid_call_table(partial_path):
            raise RuntimeError(
                f"invalid split call table for {sample_name}: {partial_path}"
            )
        os.replace(partial_path, final_path)
        results.append(final_path)
        log(f"[DONE] {sample_name} -> {final_path}")
    return results


def run_combined_fasta(
    query_fasta: Path,
    sample_names: Sequence[str],
    args: argparse.Namespace,
    scripts: Scripts,
    database_prefix: Path,
    output_root: Path,
    temp_root: Path,
) -> List[Path]:
    final_paths = [
        output_root / name / f"{name}.transcript_calls.tsv"
        for name in sample_names
    ]
    processed_paths = [output_root / name / f"{name}.pseudofragments.tsv" for name in sample_names]
    if not args.force_samples and all(valid_call_table(path) for path in final_paths + processed_paths):
        log(f"[SKIP] All {len(sample_names)} FASTA samples have both valid call tables")
        return final_paths

    combined_temp = temp_root / "combined_fasta"
    combined_temp.mkdir(parents=True, exist_ok=True)
    alignments = combined_temp / "combined.exon_alignments.tsv"
    combined_calls = combined_temp / "combined.transcript_calls.tsv"
    combined_processed = combined_temp / "combined.pseudofragments.tsv"
    alignment_partial = combined_temp / ".combined.exon_alignments.tsv.partial"
    calls_partial = combined_temp / ".combined.transcript_calls.tsv.partial"
    processed_partial = combined_temp / ".combined.pseudofragments.tsv.partial"
    for partial in (alignment_partial, calls_partial, processed_partial):
        if partial.exists():
            partial.unlink()

    with tempfile.TemporaryDirectory(prefix="blast_tmp_", dir=combined_temp) as blast_temp:
        align_command = make_alignment_command(
            query_fasta=query_fasta,
            output_path=alignment_partial,
            blast_temp=Path(blast_temp),
            args=args,
            scripts=scripts,
            database_prefix=database_prefix,
        )
        run_command(align_command, f"combined FASTA ({len(sample_names)} samples) align exons")
    os.replace(alignment_partial, alignments)

    call_command = make_caller_command(
        alignments=alignments,
        output_path=calls_partial,
        args=args,
        scripts=scripts,
        database_prefix=database_prefix,
        pseudofragments_output=processed_partial,
    )
    run_command(call_command, f"combined FASTA ({len(sample_names)} samples) call transcripts")
    if not all(valid_call_table(path) for path in (calls_partial, processed_partial)):
        raise RuntimeError("caller did not produce two valid combined tables")
    os.replace(processed_partial, combined_processed)
    os.replace(calls_partial, combined_calls)
    split_combined_call_table(
        combined_calls=combined_processed,
        sample_names=sample_names,
        output_root=output_root,
        table_suffix="pseudofragments",
    )
    split_calls = split_combined_call_table(
        combined_calls=combined_calls,
        sample_names=sample_names,
        output_root=output_root,
    )
    return split_calls


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be nonnegative")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def percentage(value: str) -> float:
    parsed = float(value)
    if not 0.0 <= parsed <= 100.0:
        raise argparse.ArgumentTypeError("must be between 0 and 100")
    return parsed


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(
        description=(
            "Annotate genes and transcripts on assembly FASTAs using the "
            "anchored-exon BLAST pipeline."
        )
    )
    parser.add_argument("--version", action="version", version=PIPELINE_VERSION)
    parser.add_argument("-r", "--reference", required=True, help="reference genome FASTA")
    parser.add_argument("-g", "--gff3", required=True, help="GENCODE-style GFF3 annotation")
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument("-q", "--query-list", help="two-column sample/assembly path file; one BLAST run per file")
    query_group.add_argument("--query-fasta", help="multi-record FASTA; first header token is the unique sample name and all records share one BLAST/caller run")
    parser.add_argument("-d", "--exon-database-dir", required=True, help="shared intermediate exon-database directory")
    parser.add_argument("-o", "--output", required=True, help="output root; one subdirectory per sample")
    parser.add_argument("--database-prefix", default="reference_exons", help="database filename prefix inside --exon-database-dir [reference_exons]")
    parser.add_argument("--scripts-dir", default=str(script_dir), help="directory containing the v3 exon/transcript pipeline scripts")
    parser.add_argument("--temp-subdir", default="temp", help="relative temporary directory below --output [temp]")
    parser.add_argument("--python", default=sys.executable, help="Python executable used for pipeline scripts [current Python]")

    parser.add_argument("--jobs", type=positive_int, default=1, help="assemblies processed concurrently [1]")
    parser.add_argument("--blast-threads", type=positive_int, default=16, help="BLAST threads per assembly [16]")
    parser.add_argument("--caller-threads", type=positive_int, default=16, help="caller processes per assembly [16]")
    parser.add_argument(
        "--caller-shard-storage",
        choices=("auto", "ram", "disk"),
        default="auto",
        help="caller shard storage; auto prefers RAM-backed /dev/shm [auto]",
    )

    parser.add_argument("--anchor-size", type=nonnegative_int, default=60, help="reference anchor bases on each exon side [60]")
    parser.add_argument("--min-unmasked", type=nonnegative_int, default=50, help="minimum unmasked bases in an anchored query [50]")
    parser.add_argument("--merge-exon-overlap", type=percentage, default=99.0, help="deprecated builder compatibility parameter [99]")
    parser.add_argument("--word-size", type=positive_int, default=19, help="BLAST word size [19]")
    parser.add_argument("--evalue", default="1e-30", help="BLAST E-value [1e-30]")
    parser.add_argument("--min-alignment-score", type=float, default=50.0, help="anchored-HSP raw score must be greater than this [50]")
    parser.add_argument("--min-exon-coverage", type=percentage, default=90.0, help="minimum anchored-query coverage [90]")
    parser.add_argument("--min-identity", type=percentage, default=95.0, help="anchored-HSP identity must be greater than this [95]")
    parser.add_argument("--max-target-seqs", type=positive_int, default=100, help="maximum assembly targets per anchored exon [100]")

    parser.add_argument("--protein-bonus", type=float, default=10.0, help="protein-coding multiplier for real transcript selection only [10]")
    parser.add_argument("--complete-bonus", type=float, default=2.0, help="complete-transcript multiplier [2]")
    parser.add_argument("--max-chains-per-transcript", type=positive_int, default=10, help="maximum chains per transcript [10]")
    parser.add_argument("--query-coordinate-mode", choices=("local", "header-suffix"), default="local", help="assembly coordinate interpretation [local]")
    parser.add_argument("--no-prefer-mane", dest="prefer_mane", action="store_false", help="disable MANE priority during transcript selection")
    parser.set_defaults(prefer_mane=True)
    parser.add_argument("--full-gene-transcripts", dest="full_gene_transcripts", action="store_true",
                        help="lift full-gene exon unions first, then isoforms within each gene [default]")
    parser.add_argument("--no-full-gene-transcripts", dest="full_gene_transcripts", action="store_false",
                        help="call only annotated transcripts")
    parser.set_defaults(full_gene_transcripts=True)

    parser.add_argument("--blastn", default="blastn", help="blastn executable [blastn]")
    parser.add_argument("--makeblastdb", default="makeblastdb", help="makeblastdb executable [makeblastdb]")
    parser.add_argument("--force-rebuild-database", action="store_true", help="rebuild even when a compatible database is present")
    parser.add_argument("--force-samples", action="store_true", help="rerun samples with valid final call tables")
    args = parser.parse_args()

    args.reference = resolve_existing_file(args.reference, "reference FASTA")
    args.gff3 = resolve_existing_file(args.gff3, "GFF3")
    if args.query_list:
        args.query_list = resolve_existing_file(args.query_list, "query list")
        args.query_fasta = None
    else:
        args.query_fasta = resolve_existing_file(args.query_fasta, "query FASTA")
        args.query_list = None
    args.scripts_dir = Path(args.scripts_dir).expanduser().resolve()
    args.exon_database_dir = Path(args.exon_database_dir).expanduser().resolve()
    args.output = Path(args.output).expanduser().resolve()
    temp_subdir = Path(args.temp_subdir)
    if temp_subdir.is_absolute() or ".." in temp_subdir.parts:
        parser.error("--temp-subdir must stay below --output")
    args.temp_subdir = temp_subdir
    if not SAMPLE_NAME_RE.fullmatch(args.database_prefix):
        parser.error("--database-prefix may contain only letters, numbers, ., _, and -")
    if args.min_identity >= 100.0:
        parser.error("--min-identity must be below 100 because selection is strict >")
    return args


def main() -> None:
    args = parse_args()
    scripts = load_scripts(args.scripts_dir)
    queries: List[AssemblyQuery] = []
    sample_names: List[str]
    if args.query_fasta is not None:
        sample_names = read_fasta_sample_names(args.query_fasta)
    else:
        queries = read_queries(args.query_list)
        sample_names = [query.name for query in queries]

    temp_top_level = args.temp_subdir.parts[0]
    conflicting_samples = [name for name in sample_names if name == temp_top_level]
    if conflicting_samples:
        raise SystemExit(
            f"ERROR: sample name {conflicting_samples[0]!r} conflicts with "
            f"temporary output directory {args.temp_subdir}"
        )
    output_root: Path = args.output
    temp_root = output_root / args.temp_subdir
    output_root.mkdir(parents=True, exist_ok=True)
    temp_root.mkdir(parents=True, exist_ok=True)
    database_prefix = args.exon_database_dir / args.database_prefix

    if args.query_fasta is not None:
        log(
            f"[LOG] Loaded {len(sample_names)} uniquely named sequences from "
            f"{args.query_fasta}"
        )
        if args.jobs != 1:
            log(
                "[LOG] --jobs is not used with --query-fasta; all sequences "
                "share one BLAST and one caller run"
            )
    else:
        log(f"[LOG] Loaded {len(queries)} assembly queries")

    database_rebuilt = ensure_database(
        prefix=database_prefix,
        scripts=scripts,
        python=args.python,
        reference=args.reference,
        gff3=args.gff3,
        anchor_size=args.anchor_size,
        min_unmasked=args.min_unmasked,
        merge_exon_overlap=args.merge_exon_overlap,
        force=args.force_rebuild_database,
    )
    if database_rebuilt:
        # Existing alignments and calls may contain genes removed by merging.
        args.force_samples = True
        log("[LOG] Refreshing sample results for the rebuilt exon database")

    if args.query_fasta is not None:
        try:
            run_combined_fasta(
                query_fasta=args.query_fasta,
                sample_names=sample_names,
                args=args,
                scripts=scripts,
                database_prefix=database_prefix,
                output_root=output_root,
                temp_root=temp_root,
            )
        except Exception as exc:
            raise SystemExit(f"ERROR: combined FASTA run failed: {exc}") from exc
        log(f"[DONE] Annotated {len(sample_names)} FASTA sequences")
        return

    failures: List[Tuple[str, Exception]] = []
    if args.jobs == 1 or len(queries) == 1:
        for query in queries:
            try:
                run_one_sample(
                    query,
                    args,
                    scripts,
                    database_prefix,
                    output_root,
                    temp_root,
                )
            except Exception as exc:
                failures.append((query.name, exc))
                log(f"[ERROR] {query.name}: {exc}")
    else:
        worker_count = min(args.jobs, len(queries))
        log(
            f"[LOG] Processing {len(queries)} assemblies with "
            f"{worker_count} concurrent jobs"
        )
        with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as pool:
            future_to_query = {
                pool.submit(
                    run_one_sample,
                    query,
                    args,
                    scripts,
                    database_prefix,
                    output_root,
                    temp_root,
                ): query
                for query in queries
            }
            for future in concurrent.futures.as_completed(future_to_query):
                query = future_to_query[future]
                try:
                    future.result()
                except Exception as exc:
                    failures.append((query.name, exc))
                    log(f"[ERROR] {query.name}: {exc}")

    if failures:
        names = ", ".join(name for name, _ in failures)
        raise SystemExit(f"ERROR: {len(failures)} sample(s) failed: {names}")
    log(f"[DONE] Annotated {len(queries)} assemblies")


if __name__ == "__main__":
    main()
