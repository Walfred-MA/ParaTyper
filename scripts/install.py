#!/usr/bin/env python3
"""Install external annotation-pipeline dependencies into the active environment."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import py_compile
import shlex
import shutil
import subprocess
import sys
from typing import List


CORE_EXECUTABLES = ("blastn", "makeblastdb", "blastdbcmd")
EXECUTABLE_PACKAGES = {
    "blastn": "blast",
    "makeblastdb": "blast",
    "blastdbcmd": "blast",
}
CORE_PIPELINE_SCRIPTS = (
    "identical_paralogs.py",
    "shared_exon_genes.py",
    "build_exon_blastdb_v2.py",
    "align_exon_blastdb_v2.py",
    "blast_gene_windows.py",
    "call_genes_from_exon_alignments_v3.py",
    "annotate_assemblies.py",
)


def missing_executables(required: tuple[str, ...]) -> List[str]:
    return [name for name in required if shutil.which(name) is None]


def check_python_scripts(script_dir: Path) -> None:
    scripts = CORE_PIPELINE_SCRIPTS
    for name in scripts:
        path = script_dir / name
        if not path.is_file():
            raise SystemExit(f"ERROR: required pipeline script is missing: {path}")
        try:
            py_compile.compile(str(path), doraise=True)
        except py_compile.PyCompileError as exc:
            raise SystemExit(f"ERROR: Python syntax check failed for {path}: {exc}") from exc


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Check the exon/transcript pipeline and install missing BLAST+ "
            "tools into the active environment."
        )
    )
    parser.add_argument("--mamba", default="mamba", help="mamba executable [mamba]")
    parser.add_argument("--dry-run", action="store_true", help="show the installation command without running it")
    args = parser.parse_args()

    if sys.version_info < (3, 9):
        raise SystemExit("ERROR: Python 3.9 or newer is required")
    script_dir = Path(__file__).resolve().parent
    check_python_scripts(script_dir)

    required = CORE_EXECUTABLES
    missing = missing_executables(required)
    if not missing:
        print("[OK] Python pipeline scripts compile", file=sys.stderr)
        print("[OK] BLAST+ is available (exon-based pipeline)", file=sys.stderr)
        return

    active_prefix = os.environ.get("CONDA_PREFIX", "")
    if not active_prefix:
        raise SystemExit(
            "ERROR: activate the target conda/mamba environment first; "
            "install.py never creates a new environment"
        )
    mamba = shutil.which(args.mamba)
    if mamba is None:
        raise SystemExit(
            f"ERROR: missing {', '.join(missing)} and mamba is not available"
        )

    packages = sorted({EXECUTABLE_PACKAGES[name] for name in missing})
    command = [
        mamba,
        "install",
        "--yes",
        "--channel",
        "conda-forge",
        "--channel",
        "bioconda",
        *packages,
    ]
    print(f"[CMD] {shlex.join(command)}", file=sys.stderr, flush=True)
    if args.dry_run:
        return
    subprocess.run(command, check=True)

    missing_after = missing_executables(required)
    if missing_after:
        raise SystemExit(
            "ERROR: installation completed but executables remain unavailable: "
            + ", ".join(missing_after)
        )
    print(
        f"[OK] Installed {', '.join(packages)} in active environment: {active_prefix}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
