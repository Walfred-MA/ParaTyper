"""Bounded local jobs, isolated scratch files, and cancellation of owned tools."""
import argparse
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import patch

from test_exon_transcript_v3 import ROOT
import blast_gene_windows as w
from minimap_candidates import parse_candidate


class ParallelSearchTests(unittest.TestCase):
    def test_windows_of_one_gene_run_concurrently_with_independent_files_and_readers(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'queries.fa'
            source.write_bytes(b'012345678')
            genes = [w.Gene(0, 9)]
            barrier = threading.Barrier(3, timeout=5)
            lock = threading.Lock()
            active = peak = 0
            directories = set()

            def fake_local(args, api, gene, windows, lengths, assembly_db, work_dir, indexed, offsets, meta, processes):
                nonlocal active, peak
                self.assertIs(gene, genes[0])
                self.assertEqual(len(windows), 1)
                window_index = windows[0][1]
                with lock:
                    active += 1
                    peak = max(peak, active)
                    self.assertNotIn(work_dir, directories)
                    directories.add(work_dir)
                try:
                    scratch = Path(work_dir)/'gene_queries.fa'
                    scratch.write_text(str(window_index))
                    indexed.seek(window_index)
                    barrier.wait()
                    self.assertEqual(indexed.read(1), str(window_index).encode())
                    self.assertEqual(scratch.read_text(), str(window_index))
                    yield window_index
                finally:
                    with lock:
                        active -= 1

            with patch.object(w, 'local_rows', side_effect=fake_local):
                result = list(w.parallel_local_rows(
                    (w.WindowJob(i, 0, ('chr1', i, i+1)) for i in range(9)), 9,
                    argparse.Namespace(threads=3), None, genes, {}, 'db', tmp, str(source), {}, {}))
            self.assertEqual(sorted(result), list(range(9)))
            self.assertEqual(peak, 3)
            self.assertEqual(len(directories), 9)
            self.assertTrue(all(not Path(path).exists() for path in directories))
            self.assertEqual(list(Path(tmp).iterdir()), [source])

    def test_failure_stops_other_owned_processes_and_does_not_start_more_windows(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'queries.fa'; source.write_text('ACGT')
            started = threading.Event()
            launched = []
            visited = []

            def fake_local(args, api, gene, windows, lengths, assembly_db, work_dir, indexed, offsets, meta, processes):
                self.assertEqual(len(windows), 1)
                window_index = windows[0][1]
                visited.append(window_index)
                if window_index == 0:
                    self.assertTrue(started.wait(5))
                    raise ValueError('original window failure')
                proc = processes.start([sys.executable, '-c', 'import time; time.sleep(30)'])
                launched.append(proc)
                started.set()
                try:
                    proc.wait()
                    yield window_index
                finally:
                    processes.release(proc)

            start = time.monotonic()
            with patch.object(w, 'local_rows', side_effect=fake_local):
                with self.assertRaisesRegex(ValueError, 'original window failure'):
                    list(w.parallel_local_rows(
                        (w.WindowJob(i, 0, ('chr1', i, i+1)) for i in range(10)), 10, argparse.Namespace(threads=2),
                        None, [w.Gene(0, 10)], {}, 'db', tmp, str(source), {}, {}))
            self.assertLess(time.monotonic()-start, 5)
            self.assertEqual(set(visited), {0, 1})
            self.assertTrue(launched)
            self.assertTrue(all(proc.poll() is not None for proc in launched))

    def test_mapq_zero_and_partial_chains_are_valid_candidates(self):
        line = 'E1\t150\t50\t100\t-\tchr1\t10000\t500\t550\t40\t50\t0\ttp:A:S'
        self.assertEqual(parse_candidate(line), ('E1', 'chr1', 10000, 500, 550))
        with self.assertRaises(ValueError):
            parse_candidate('bad\trow')
        with self.assertRaises(ValueError):
            parse_candidate(line.replace('\t500\t550\t', '\t-1\t550\t'))


if __name__ == '__main__':
    unittest.main()
