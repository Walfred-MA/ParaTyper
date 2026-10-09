"""Installer checks for the pinned Snakemake used by the Slurm workflow."""

import os
import subprocess
import sys
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT
import install as installer


class InstallTests(unittest.TestCase):
    def test_reads_exact_snakemake_version(self):
        with patch.object(installer.shutil, "which", return_value="/env/bin/snakemake"), \
                patch.object(installer.subprocess, "run", return_value=subprocess.CompletedProcess(
                    ["snakemake", "--version"], 0, stdout="6.15.1\n")), \
                patch.dict(os.environ, {"CONDA_PREFIX": ""}):
            self.assertEqual(installer.installed_snakemake_version(), "6.15.1")

    def test_snakemake_outside_active_environment_is_not_accepted(self):
        with patch.object(installer.shutil, "which", return_value="/other-env/bin/snakemake"), \
                patch.dict(os.environ, {"CONDA_PREFIX": "/nonexistent/paratyper-env"}):
            self.assertIsNone(installer.installed_snakemake_version())

    def test_missing_snakemake_installs_pinned_packages_into_active_env(self):
        versions = iter([None, "6.15.1"])
        with patch.object(installer, "check_python_scripts"), \
                patch.object(installer, "missing_executables", return_value=[]), \
                patch.object(installer, "installed_snakemake_version", side_effect=lambda: next(versions)), \
                patch.object(installer.shutil, "which", return_value="/usr/bin/mamba"), \
                patch.object(installer.subprocess, "run") as run, \
                patch.dict(os.environ, {"CONDA_PREFIX": "/tmp/paratyper-env"}), \
                patch.object(sys, "argv", ["install.py"]):
            installer.main()
        command = run.call_args.args[0]
        self.assertEqual(command[:4], ["/usr/bin/mamba", "install", "--prefix", "/tmp/paratyper-env"])
        self.assertIn("snakemake-minimal=6.15.1", command)
        self.assertIn("tabulate=0.8.10", command)

    def test_current_snakemake_needs_no_mamba_install(self):
        with patch.object(installer, "check_python_scripts"), \
                patch.object(installer, "missing_executables", return_value=[]), \
                patch.object(installer, "installed_snakemake_version", return_value="6.15.1"), \
                patch.object(installer.subprocess, "run") as run, \
                patch.object(sys, "argv", ["install.py"]):
            installer.main()
        run.assert_not_called()


if __name__ == "__main__":
    unittest.main()
