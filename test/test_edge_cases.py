"""Malformed-input and boundary-condition coverage not exercised elsewhere: missing manifest
columns, dangling file references, duplicate IDs across rows, non-UTF-8 input, an all-rejected
cohort, a missing `mash` binary, out-of-order CLI bounds, and a tampered provenance file."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "plasmicord.py"
sys.path.insert(0, str(ROOT))
from python.contracts import write_tsv, validate_manifest


def digest(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


SEQ = "ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC"


class EdgeCases(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *args, timeout=120):
        return subprocess.run([sys.executable, str(CLI), *map(str, args)],
                              cwd=ROOT, capture_output=True, text=True, timeout=timeout)

    def minimal_fixture(self, seq=SEQ, isolate="i1", plasmid="p1"):
        fasta = self.p / f"{plasmid}.fa"; fasta.write_text(f">{plasmid}\n{seq}\n")
        manifest = self.p / "manifest.tsv"
        write_tsv(manifest, ["isolate_id", "plasmid_id", "fasta_path"],
                  [dict(isolate_id=isolate, plasmid_id=plasmid, fasta_path=str(fasta))])
        meta = self.p / "metadata.tsv"
        write_tsv(meta, ["isolate_id"], [dict(isolate_id=isolate)])
        return manifest, meta

    # -- malformed manifest content ----------------------------------------
    def test_manifest_missing_required_column_entirely_is_rejected(self):
        manifest = self.p / "manifest.tsv"
        # No fasta_path column at all, not merely an empty value.
        write_tsv(manifest, ["isolate_id", "plasmid_id"], [dict(isolate_id="i1", plasmid_id="p1")])
        with self.assertRaisesRegex(ValueError, "missing required or duplicate columns"):
            validate_manifest(manifest, [dict(isolate_id="i1")], "precomputed", 1)

    def test_manifest_row_referencing_nonexistent_fasta_fails_cleanly_via_cli(self):
        manifest, meta = self.minimal_fixture()
        manifest.write_text(manifest.read_text().replace(str(self.p / "p1.fa"), str(self.p / "does_not_exist.fa")))
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ERROR", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_duplicate_plasmid_id_across_manifest_rows_is_rejected(self):
        fasta1 = self.p / "a.fa"; fasta1.write_text(">a\nACGTACGTACGTACGTACGT\n")
        fasta2 = self.p / "b.fa"; fasta2.write_text(">b\nTGCATGCATGCATGCATGCA\n")
        manifest = self.p / "manifest.tsv"
        write_tsv(manifest, ["isolate_id", "plasmid_id", "fasta_path"],
                  [dict(isolate_id="i1", plasmid_id="dup", fasta_path=str(fasta1)),
                   dict(isolate_id="i2", plasmid_id="dup", fasta_path=str(fasta2))])
        meta = [dict(isolate_id="i1"), dict(isolate_id="i2")]
        with self.assertRaisesRegex(ValueError, "duplicate"):
            validate_manifest(manifest, meta, "precomputed", 1)

    def test_non_utf8_binary_fasta_fails_cleanly_via_cli_not_with_a_traceback(self):
        manifest, meta = self.minimal_fixture()
        (self.p / "p1.fa").write_bytes(b"\xff\xfe\x00\x01>p1\n\x80\x81ACGT\n")
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ERROR", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    # -- boundary cohorts -----------------------------------------------
    def test_all_candidates_rejected_run_still_completes_with_zero_counts(self):
        manifest, meta = self.minimal_fixture(seq=SEQ)
        out = self.p / "out"
        # min-length far exceeds the only candidate's length: every candidate is rejected,
        # not merely absent -- a distinct code path from an empty manifest.
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", out,
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3",
                          "--min-length", str(len(SEQ) + 1000))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out / "REPORT.html").is_file())
        record = json.loads((out / "run_provenance.json").read_text(encoding="utf-8"))
        self.assertEqual(record["status"], "complete")
        self.assertEqual(record["n_plasmids"], 0)
        self.assertEqual(record["n_rejected"], 1)
        payload = json.loads((out / "report_data.json").read_text(encoding="utf-8"))
        self.assertEqual(payload["counts"]["plasmids"], 0)

    def test_single_isolate_single_plasmid_cohort_via_real_cli(self):
        manifest, meta = self.minimal_fixture()
        out = self.p / "out"
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", out,
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out / "REPORT.html").is_file())

    def test_mash_engine_without_binary_on_path_fails_cleanly(self):
        manifest, meta = self.minimal_fixture()
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "mash", "--threshold", "0.5", "--k", "3", "--min-length", "3")
        if result.returncode == 0:
            self.skipTest("A mash binary is on PATH in this environment; the absent-binary path cannot be exercised here")
        self.assertEqual(result.returncode, 1)
        self.assertIn("Mash is not on PATH", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    # -- CLI-level bounds validation (existing coverage is direct run() only) ----
    def test_containment_min_ratio_above_max_ratio_via_real_cli(self):
        manifest, meta = self.minimal_fixture()
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3",
                          "--containment-min-ratio", "0.9", "--containment-max-ratio", "0.5")
        self.assertEqual(result.returncode, 1)
        self.assertIn("containment-min-ratio", result.stderr)
        self.assertNotIn("Traceback", result.stderr)
        self.assertFalse((self.p / "out").exists())

    def test_negative_size_correction_cap_pct_via_real_cli(self):
        manifest, meta = self.minimal_fixture()
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3",
                          "--size-correction-cap-pct", "-10")
        self.assertEqual(result.returncode, 1)
        self.assertIn("size-correction-cap-pct", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    # -- tampered provenance before regeneration -----------------------------
    def test_tampered_run_provenance_before_report_regeneration_is_rejected(self):
        manifest, meta = self.minimal_fixture()
        out = self.p / "out"
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", out,
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "3")
        self.assertEqual(result.returncode, 0, result.stderr)
        provenance_path = out / "run_provenance.json"
        record = json.loads(provenance_path.read_text(encoding="utf-8"))
        # Corrupt one recorded checksum on a stable (non-regenerated) analysis file so the
        # regenerated report can no longer trust it; regenerated presentation files themselves
        # are deliberately excluded from this check (report_output.py::refresh_report's GENERATED
        # set), so target a source-of-truth file such as metadata.tsv instead.
        record["output_sha256"]["metadata.tsv"] = "0" * 64
        provenance_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        result = self.cli("report", "--results", out)
        self.assertEqual(result.returncode, 1)
        self.assertIn("ERROR", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
