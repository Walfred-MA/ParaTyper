"""Slurm orchestration checks that do not require a cluster."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT, runner


class SlurmRunnerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder = Path(temporary.name)
        self.reference = self.folder / "reference.fa"
        self.gff3 = self.folder / "genes.gff3"
        self.fasta = self.folder / "assembly.fa"
        self.query_list = self.folder / "queries.txt"
        self.reference.write_text(">chr1\nACGT\n")
        self.gff3.write_text("##gff-version 3\n")
        self.fasta.write_text(">contig\nACGT\n")
        self.query_list.write_text(f"sampleA {self.fasta}\n")
        self.output = self.folder / "output"
        self.output.mkdir()
        self.args = argparse.Namespace(
            reference=self.reference, gff3=self.gff3, query_list=self.query_list,
            query_fasta=None, exon_database_dir=self.folder / "database",
            database_prefix="reference_exons", output=self.output,
            scripts_dir=ROOT, temp_subdir=Path("temp"), python=sys.executable,
            jobs=1, slurm=2, slurm_command="--account=lab --partition=compute --mem=64G --time=24:00:00 --cpus-per-task=16",
            snakemake="snakemake", force_samples=False,
            force_rebuild_database=False, anchor_target_length=150,
            min_unmasked=50, merge_exon_overlap=99.0,
        )

    def test_quoted_sbatch_options_are_one_cli_value(self):
        argv = [
            "annotate_assemblies.py", "-r", str(self.reference), "-g", str(self.gff3),
            "-q", str(self.query_list), "-d", str(self.folder / "database"),
            "-o", str(self.output), "--slurm", "2", "--slurm-command",
            "--account=lab --partition=compute --mem=64G",
        ]
        with patch.object(sys, "argv", argv):
            args = runner.parse_args()
        self.assertEqual(args.slurm, 2)
        self.assertEqual(runner.slurm_command_options(args.slurm_command), [
            "sbatch", "--account=lab", "--partition=compute", "--mem=64G",
        ])

    def test_snakemake_owns_job_limit_and_output_lock(self):
        calls = []

        def fake_run(command, **kwargs):
            calls.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="6.15.1\n")

        with patch.object(runner, "current_database", return_value=(True, "current")), \
                patch.object(runner.subprocess, "run", side_effect=fake_run):
            runner.run_snakemake_slurm(
                self.args, [runner.AssemblyQuery("sampleA", self.fasta)],
                self.args.exon_database_dir / self.args.database_prefix,
            )
        command = calls[-1]
        self.assertEqual(command[command.index("--directory") + 1], str(self.output))
        self.assertEqual(command[command.index("--jobs") + 1], "2")
        self.assertIn("slurm_jobs=2", command)
        self.assertNotIn("--forcerun", command)
        config = json.loads(Path(command[command.index("--configfile") + 1]).read_text())
        self.assertTrue(config["args"]["force_samples"])
        self.assertEqual(config["queries"], {"sampleA": str(self.fasta)})
        self.assertEqual(config["slurm_command"], self.args.slurm_command)

    def test_rebuild_forces_database_and_sample_rules(self):
        calls = []

        def fake_run(command, **kwargs):
            calls.append(command)
            return subprocess.CompletedProcess(command, 0, stdout="6.15.1\n")

        with patch.object(runner, "current_database", return_value=(False, "missing")), \
                patch.object(runner.subprocess, "run", side_effect=fake_run):
            runner.run_snakemake_slurm(
                self.args, [runner.AssemblyQuery("sampleA", self.fasta)],
                self.args.exon_database_dir / self.args.database_prefix,
            )
        command = calls[-1]
        self.assertEqual(command[command.index("--forcerun") + 1:],
                         ["prepare_database", "annotate_sample"])

    def test_one_submission_uses_custom_sbatch_options_and_exit_file(self):
        spec = self.output / "spec.json"
        spec.write_text(json.dumps({
            "args": {"output": str(self.output), "python": sys.executable},
            "slurm_command": self.args.slurm_command,
        }))
        submitted = []

        def fake_run(command, **kwargs):
            submitted.append(command)
            script = Path(command[-1])
            (script.parent / (script.stem + ".exit")).write_text("0\n")
            return subprocess.CompletedProcess(command, 0, stdout="12345\n")

        with patch.object(runner.subprocess, "run", side_effect=fake_run):
            runner.run_slurm_task(spec, "sample", "sampleA")
        command = submitted[0]
        self.assertEqual(command[0], "sbatch")
        self.assertIn("--account=lab", command)
        self.assertIn("--cpus-per-task=16", command)
        self.assertIn("--parsable", command)
        self.assertIn("sampleA", Path(command[-1]).read_text())

    def test_failed_slurm_task_is_reported_to_snakemake(self):
        spec = self.output / "spec.json"
        spec.write_text(json.dumps({
            "args": {"output": str(self.output), "python": sys.executable},
            "slurm_command": "--mem=64G",
        }))

        def fake_run(command, **kwargs):
            script = Path(command[-1])
            (script.parent / (script.stem + ".exit")).write_text("9\n")
            return subprocess.CompletedProcess(command, 0, stdout="12345\n")

        with patch.object(runner.subprocess, "run", side_effect=fake_run):
            with self.assertRaisesRegex(RuntimeError, "Slurm job 12345 failed"):
                runner.run_slurm_task(spec, "sample", "sampleA")


if __name__ == "__main__":
    unittest.main()
