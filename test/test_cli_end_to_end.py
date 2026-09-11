"""Real-subprocess CLI coverage: every plasmicord.py subcommand, invoked exactly as a user
would from a terminal (argparse parsing, exit codes and stderr text included), not just the
underlying Python function called directly with a hand-built argparse.Namespace."""
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
from python.contracts import write_tsv


def digest(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


class CliEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *args, timeout=120):
        return subprocess.run([sys.executable, str(CLI), *map(str, args)],
                              cwd=ROOT, capture_output=True, text=True, timeout=timeout)

    # -- top level -----------------------------------------------------
    def test_version_help_and_unknown_subcommand(self):
        result = self.cli("--version")
        self.assertEqual(result.returncode, 0)
        self.assertIn("plasmicord", result.stdout.lower())
        result = self.cli("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("usage", result.stdout.lower())
        result = self.cli("not-a-real-subcommand")
        self.assertEqual(result.returncode, 2)
        self.assertIn("usage", result.stderr.lower())

    def test_check_reports_expected_keys(self):
        result = self.cli("check")
        self.assertEqual(result.returncode, 0)
        payload = json.loads(result.stdout)
        for key in ("python", "standalone_core", "mash", "mob_recon", "bakta", "prokka",
                    "amrfinder", "mob_typer", "automatic_annotation"):
            self.assertIn(key, payload)

    # -- demo ------------------------------------------------------------
    def test_demo_end_to_end_and_rejects_nonempty_output(self):
        out = self.p / "demo_out"
        result = self.cli("demo", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out / "results" / "REPORT.html").is_file())
        self.assertTrue((out / "results" / "report_data.json").is_file())
        result = self.cli("demo", "--out", out)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("Traceback", result.stderr)
        self.assertIn("ERROR", result.stderr)

    # -- run: fixture helper ---------------------------------------------
    def build_run_fixture(self):
        """Three isolates exercising sharing, containment and typing in one real run:
        i1/i2 share an identical plasmid pA (-> sharing edge, multilayer edge); i3 carries
        an 80%-length prefix of pA (-> containment candidate at the default distance bound)."""
        seq_a = "ACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTACACGTTGCAACGTTCAGGATCCGATACCTAGCTGACTGGTAC"
        seq_a_short = seq_a[:int(len(seq_a) * 0.8)]
        p1 = self.p / "p1.fa"; p1.write_text(f">p1\n{seq_a}\n")
        p2 = self.p / "p2.fa"; p2.write_text(f">p2\n{seq_a}\n")
        p3 = self.p / "p3.fa"; p3.write_text(f">p3\n{seq_a_short}\n")
        manifest = self.p / "manifest.tsv"
        write_tsv(manifest, ["isolate_id", "plasmid_id", "fasta_path"],
                  [dict(isolate_id="i1", plasmid_id="p1", fasta_path=str(p1)),
                   dict(isolate_id="i2", plasmid_id="p2", fasta_path=str(p2)),
                   dict(isolate_id="i3", plasmid_id="p3", fasta_path=str(p3))])
        meta = self.p / "metadata.tsv"
        write_tsv(meta, ["isolate_id", "chromosomal_cluster"],
                  [dict(isolate_id="i1", chromosomal_cluster="CC1"),
                   dict(isolate_id="i2", chromosomal_cluster="CC1"),
                   dict(isolate_id="i3", chromosomal_cluster="CC2")])
        features = self.p / "features.tsv"
        write_tsv(features, ["plasmid_id", "feature_id", "start", "end", "strand", "amr_gene",
                             "drug_class", "hit_class", "annotation_confidence", "annotation_engine",
                             "database_name", "database_version"],
                  [dict(plasmid_id="p1", feature_id="f1", start=1, end=10, strand="+",
                        amr_gene="blaTEST", drug_class="beta-lactam", hit_class="strict",
                        annotation_confidence="high", annotation_engine="amrfinder",
                        database_name="amrfinderdb", database_version="2026")])
        quality_evidence = self.p / "quality_evidence.tsv"
        write_tsv(quality_evidence, ["plasmid_id", "sequence_sha256", "evidence_source", "copy_number"],
                  [dict(plasmid_id="p1", sequence_sha256=digest(seq_a), evidence_source="lab_estimate",
                        copy_number="3.5")])
        external_typing = self.p / "external_typing.tsv"
        write_tsv(external_typing, ["plasmid_id", "sequence_sha256", "evidence_source", "external_tool",
                                    "mob_primary_cluster_id"],
                  [dict(plasmid_id="p1", sequence_sha256=digest(seq_a), evidence_source="mob_suite_run",
                        external_tool="mob_suite", mob_primary_cluster_id="AA123")])
        return manifest, meta, features, quality_evidence, external_typing

    def test_run_kitchen_sink_populates_every_optional_dataset_in_report_data(self):
        manifest, meta, features, quality_evidence, external_typing = self.build_run_fixture()
        out = self.p / "results"
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--features", features,
                          "--quality-evidence", quality_evidence, "--external-typing", external_typing,
                          "--out", out, "--mode", "precomputed", "--engine", "kmer", "--threshold", "0.5",
                          "--linkage", "complete", "--k", "3", "--min-length", "10",
                          "--multilayer-network",
                          "--containment-min-ratio", "0.5", "--containment-max-ratio", "0.95",
                          "--containment-max-distance", "0.5",
                          "--size-correction-per-percent", "0.1", "--size-correction-cap-pct", "50")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out / "REPORT.html").is_file())
        payload = json.loads((out / "report_data.json").read_text(encoding="utf-8"))
        # Regression guard for the "silo" fix: these three optional datasets must reach the
        # payload through the real CLI path, not just a direct run() call.
        self.assertTrue(payload["containment"], "containment candidates must be non-empty for this fixture")
        self.assertTrue(payload["typing_crossreference"], "typing cross-reference must be non-empty for this fixture")
        self.assertTrue(payload["multilayer_edges"], "multilayer edges must be non-empty for this fixture")
        rule_ids = {f["rule_id"] for f in payload["dashboard"]["findings"]}
        self.assertIn("containment_candidates", rule_ids)
        self.assertIn("typing_coverage_external", rule_ids)
        self.assertIn("multilayer_network", rule_ids)

    def test_run_missing_required_threshold_exits_via_argparse(self):
        manifest, meta, *_ = self.build_run_fixture()
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", self.p / "out",
                          "--engine", "kmer")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--threshold", result.stderr)

    def test_run_rejects_features_together_with_annotation_config(self):
        manifest, meta, features, *_ = self.build_run_fixture()
        config = self.p / "annotation.json"; config.write_text("{}")
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--features", features,
                          "--annotation-config", config, "--out", self.p / "out", "--engine", "kmer",
                          "--threshold", "0.5", "--min-length", "10", "--k", "3")
        self.assertEqual(result.returncode, 1)
        self.assertIn("ERROR", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    # -- import: one real-CLI invocation per adapter, plus --samples batch mode ------------
    def test_import_generic_fasta_adapter(self):
        fasta = self.p / "candidate.fasta"; fasta.write_text(">contig1\nACGTACGTACGTACGTACGT\n")
        out = self.p / "imported"
        result = self.cli("import", "generic-fasta", "--input", fasta, "--isolate-id", "iso1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Imported 1 candidates", result.stdout)
        self.assertTrue((out / "manifest.tsv").is_file())

    def test_import_mob_recon_adapter(self):
        folder = self.p / "mob_recon_out"; folder.mkdir()
        (folder / "plasmid_AA.fasta").write_text(">a\nACGTACGTACGTACGTACGT\n")
        (folder / "contig_report.txt").write_text("contig_id\tmolecule_type\na\tplasmid\n")
        out = self.p / "imported"
        result = self.cli("import", "mob-recon", "--input", folder, "--isolate-id", "iso1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_import_flye_adapter_requires_record_ids(self):
        folder = self.p / "flye_out"; folder.mkdir()
        (folder / "assembly.fasta").write_text(">contig_1\nACGTACGTACGTACGTACGT\n")
        (folder / "assembly_info.txt").write_text("#seq_name\tlength\tcov.\tcirc.\ncontig_1\t20\t30\tY\n")
        out = self.p / "imported"
        result = self.cli("import", "flye", "--input", folder, "--isolate-id", "iso1",
                          "--record-ids", "contig_1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_import_unicycler_adapter(self):
        folder = self.p / "unicycler_out"; folder.mkdir()
        (folder / "assembly.fasta").write_text(">1 circular=true\nACGTACGTACGTACGTACGT\n")
        out = self.p / "imported"
        result = self.cli("import", "unicycler", "--input", folder, "--isolate-id", "iso1",
                          "--record-ids", "1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_import_plasbench_adapter(self):
        folder = self.p / "plasbench_out"; folder.mkdir()
        (folder / "selection_report.json").write_text(json.dumps(dict(selected_tool="mob_recon",
                                                                       selection_type="best_candidate")))
        (folder / "candidate.plasmid.fasta").write_text(">a\nACGTACGTACGTACGTACGT\n")
        out = self.p / "imported"
        result = self.cli("import", "plasbench", "--input", folder, "--isolate-id", "iso1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_import_tadrep_adapter(self):
        folder = self.p / "tadrep_out"; folder.mkdir()
        (folder / "sample1-refA-pseudo.fna").write_text(">refA\nACGTACGTACGTACGTACGT\n")
        header = "plasmid\tcontig\tcopy\tcoverage[%]\tidentity[%]\talignment length\tstrand\tstart\tend\treflen"
        (folder / "sample1-summary.tsv").write_text(header + "\nrefA\tcontig_1\t1\t95.0\t99.0\t20\t+\t1\t20\t20\n")
        out = self.p / "imported"
        result = self.cli("import", "tadrep", "--input", folder, "--isolate-id", "iso1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_import_samples_batch_mode(self):
        fasta1 = self.p / "c1.fasta"; fasta1.write_text(">contig1\nACGTACGTACGTACGTACGT\n")
        fasta2 = self.p / "c2.fasta"; fasta2.write_text(">contig2\nTGCATGCATGCATGCATGCA\n")
        samples = self.p / "samples.tsv"
        write_tsv(samples, ["isolate_id", "input_path"],
                  [dict(isolate_id="iso1", input_path=str(fasta1)),
                   dict(isolate_id="iso2", input_path=str(fasta2))])
        out = self.p / "imported"
        result = self.cli("import", "generic-fasta", "--samples", samples, "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Imported 2 candidates", result.stdout)

    # -- calibrate ---------------------------------------------------------
    def calibration_fixture(self):
        ids = list("ABCDEFGH")
        matrix = [[0 if i == j else .5 for j in range(8)] for i in range(8)]
        matrix[0][1] = matrix[1][0] = .005
        matrix[4][5] = matrix[5][4] = .02
        matrix_path = self.p / "matrix.tsv"
        matrix_path.write_text("item\t" + "\t".join(ids) + "\n" +
                               "".join(a + "\t" + "\t".join(map(str, row)) + "\n" for a, row in zip(ids, matrix)))
        index_path = self.p / "index.tsv"
        write_tsv(index_path, ["plasmid_id", "isolate_id", "sequence_sha256"],
                  [dict(plasmid_id=i, isolate_id=i, sequence_sha256=digest(i)) for i in ids])
        labels_path = self.p / "labels.tsv"
        rows = [dict(plasmid_a=a, plasmid_b=b, label=label, split=split, group_id=split,
                    evidence_source="independent laboratory evidence", evidence_type="experimental")
                for a, b, label, split in [("A", "B", "positive", "train"), ("C", "D", "negative", "train"),
                                           ("E", "F", "positive", "validation"), ("G", "H", "negative", "validation")]]
        write_tsv(labels_path, list(rows[0]), rows)
        return matrix_path, index_path, labels_path

    def test_calibrate_happy_path_via_real_cli(self):
        matrix_path, index_path, labels_path = self.calibration_fixture()
        out = self.p / "calibration"
        result = self.cli("calibrate", "--matrix", matrix_path, "--index", index_path,
                          "--labels", labels_path, "--cohort-kind", "real", "--engine", "mash",
                          "--thresholds", "0.01,0.03", "--min-groups", "1", "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads((out / "calibration.json").read_text(encoding="utf-8"))
        self.assertEqual(record["selected_threshold"], .01)

    def test_calibrate_leakage_rejection_via_real_cli_is_clean(self):
        matrix_path, index_path, labels_path = self.calibration_fixture()
        rows = [dict(plasmid_a="A", plasmid_b="B", label="positive", split="train", group_id="train",
                    evidence_source="lab", evidence_type="experimental"),
                dict(plasmid_a="C", plasmid_b="D", label="negative", split="train", group_id="train",
                    evidence_source="lab", evidence_type="experimental"),
                dict(plasmid_a="E", plasmid_b="F", label="positive", split="validation", group_id="train",
                    evidence_source="lab", evidence_type="experimental"),
                dict(plasmid_a="G", plasmid_b="H", label="negative", split="validation", group_id="validation",
                    evidence_source="lab", evidence_type="experimental")]
        write_tsv(labels_path, list(rows[0]), rows)
        out = self.p / "calibration"
        result = self.cli("calibrate", "--matrix", matrix_path, "--index", index_path,
                          "--labels", labels_path, "--cohort-kind", "real", "--engine", "mash",
                          "--out", out)
        self.assertEqual(result.returncode, 1)
        self.assertIn("group spans", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    # -- report / population-summary regeneration --------------------------
    def test_report_regeneration_and_calibration_attachment_via_real_cli(self):
        manifest, meta, *_ = self.build_run_fixture()
        out = self.p / "results"
        result = self.cli("run", "--manifest", manifest, "--metadata", meta, "--out", out,
                          "--engine", "kmer", "--threshold", "0.5", "--k", "3", "--min-length", "10")
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.cli("report", "--results", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.cli("population-summary", "--results", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((out / "population_summary.pu_level.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
