"""Opt-in, network-touching smoke test: downloads a small real, accession-backed plasmid
cohort from NCBI and runs the full pipeline against it. Skipped by default -- PlasmiCord's
core analysis and its default test suite require no network access (see README.md,
docs/POSITIONING.md). Enable explicitly:

    PLASMICORD_REAL_COHORT_TEST=1 python -m unittest discover -s test -p test_real_cohort.py -v

or `make real-cohort-test`. This checks pipeline completeness on genuine sequence data, not
any specific biological outcome (no transmission truth is asserted or implied)."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


class RealCohortSmokeTest(unittest.TestCase):
    def setUp(self):
        if os.environ.get("PLASMICORD_REAL_COHORT_TEST") != "1":
            self.skipTest("Set PLASMICORD_REAL_COHORT_TEST=1 to run this network-touching test "
                          "(see `make real-cohort-test`)")
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)

    def tearDown(self):
        if hasattr(self, "tmp"):
            self.tmp.cleanup()

    def test_download_and_run_small_real_cohort(self):
        downloaded = self.p / "downloaded"
        download = subprocess.run([sys.executable, str(ROOT / "scripts" / "fetch_real_cohort.py"),
                                   "--out", str(downloaded)], cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(download.returncode, 0, f"stdout={download.stdout}\nstderr={download.stderr}")
        cohort_rows = (downloaded / "cohort.tsv").read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(cohort_rows) - 1, 3, "expected the default 3-isolate cohort")

        prepared = self.p / "prepared"
        prepare = subprocess.run([sys.executable, str(ROOT / "scripts" / "prepare_reference_validation.py"),
                                  "--cohort", str(downloaded / "cohort.tsv"),
                                  "--references", str(downloaded / "references"),
                                  "--out", str(prepared)], cwd=ROOT, capture_output=True, text=True, timeout=60)
        self.assertEqual(prepare.returncode, 0, f"stdout={prepare.stdout}\nstderr={prepare.stderr}")
        manifest_rows = (prepared / "manifest.tsv").read_text(encoding="utf-8").strip().splitlines()
        self.assertEqual(len(manifest_rows) - 1, 6, "expected 6 real plasmids across the 3 isolates")

        out = self.p / "results"
        # Real deployments cannot assume mash is installed; the exact k-mer engine is
        # explicitly designed for small sets (<=200 plasmids), which this cohort easily satisfies.
        run = subprocess.run([sys.executable, str(ROOT / "plasmicord.py"), "run",
                              "--manifest", str(prepared / "manifest.tsv"),
                              "--metadata", str(prepared / "metadata.tsv"),
                              "--out", str(out), "--mode", "precomputed", "--engine", "kmer",
                              "--threshold", "0.05", "--linkage", "complete", "--k", "21",
                              "--min-length", "200"], cwd=ROOT, capture_output=True, text=True, timeout=120)
        self.assertEqual(run.returncode, 0, f"stdout={run.stdout}\nstderr={run.stderr}")
        self.assertTrue((out / "REPORT.html").is_file())
        payload = json.loads((out / "report_data.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["counts"]["isolates"], 3)
        self.assertEqual(payload["counts"]["plasmids"], 6)
        # Pipeline-completeness check only: at least one dataset the pipeline can compute from
        # real sequence relationships (a sharing edge or a containment candidate) was evaluated.
        # No specific biological/transmission outcome is asserted.
        self.assertTrue(payload["edges"] or payload["containment"],
                        "expected at least one evaluated sharing edge or containment candidate on real sequence data")
        record = json.loads((out / "run_provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(record["status"], "complete")


if __name__ == "__main__":
    unittest.main()
