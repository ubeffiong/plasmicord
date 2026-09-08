"""Regression tests for scientific false links, validation and standalone outputs."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python.cli import run
from python.contracts import fasta_records, metadata, normalize_features, validate_manifest, write_tsv
from python.cluster_plasmids import read_matrix, cluster_complete
from python.kmer_distance import canonical_kmers, read_fasta_concat


class StandaloneTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def script(self, script, *args):
        return subprocess.run([sys.executable, str(ROOT / "python" / script), *map(str, args)], capture_output=True, text=True)

    def fixture(self, sequences=("ACGTACGTTGCA",), cluster="A"):
        meta = self.path / "metadata.tsv"
        write_tsv(meta, ["isolate_id", "chromosomal_cluster"], [{"isolate_id": "iso1", "chromosomal_cluster": cluster}, {"isolate_id": "iso2", "chromosomal_cluster": "B"}])
        rows = []
        for i, seq in enumerate(sequences):
            (self.path / f"p{i}.fa").write_text(f">c{i}\n{seq}\n")
            rows.append(dict(isolate_id=f"iso{min(i+1,2)}", plasmid_id=f"p{i}", fasta_path=f"p{i}.fa"))
        manifest = self.path / "manifest.tsv"
        write_tsv(manifest, ["isolate_id", "plasmid_id", "fasta_path"], rows)
        return argparse.Namespace(manifest=manifest, metadata=meta, features=None, mode="precomputed", engine="kmer", out=self.path / "out", threshold=.05, k=3, min_length=3, sketch_size=100, linkage="complete")

    def test_unknown_cluster_never_generates_crosslink(self):
        args = self.fixture(("ACGTACGTTGCA", "ACGTACGTTGCA"), "UNKNOWN")
        run(args)
        data = json.loads((args.out / "report_data.json").read_text())
        self.assertEqual(data["counts"]["sharing_pairs"], 1)
        self.assertEqual(data["counts"]["cross_cluster_pairs"], 0)
        self.assertEqual(data["counts"]["typed_isolates"], 1)
        self.assertIsNone(data["counts"]["observed_arg_plasmids"])
        self.assertEqual(data["plasmids"][1]["duplicate_of"], "p0")
        self.assertIn('REPORT.html', json.loads((args.out / "run_provenance.json").read_text())["output_sha256"])
        ET.parse(args.out / "network.graphml")
        with self.assertRaisesRegex(ValueError, "not empty"):
            run(args)

    def test_empty_and_singleton_runs(self):
        for seqs in ((), ("ACGTACGTTGCA",)):
            args = self.fixture(seqs)
            args.out = self.path / ("empty" if not seqs else "single")
            run(args)
            data = json.loads((args.out / "report_data.json").read_text())
            self.assertEqual(data["counts"]["units"], len(seqs))
            self.assertEqual(data["counts"]["isolates"], 2)

    def test_matrix_rejects_invalid_evidence(self):
        for values in ["item\tA\tA\nA\t0\t0\nA\t0\t0\n", "item\tA\tB\nA\t0\t.1\nB\t.2\t0\n", "item\tA\nA\tnan\n", "item\tA\nA\t.1\n", "item\tA\nA\t-1\n"]:
            p = self.path / "matrix.tsv"
            p.write_text(values)
            with self.assertRaises(SystemExit):
                read_matrix(p)

    def test_complete_linkage_ties_do_not_depend_on_input_order(self):
        original = [[0, .1, .3], [.1, 0, .1], [.3, .1, 0]]
        results = []
        for order in itertools.permutations(range(3)):
            ids = ["ABC"[i] for i in order]
            matrix = [[original[i][j] for j in order] for i in order]
            results.append({frozenset(ids[i] for i in c) for c in cluster_complete(ids, matrix, .1)})
        self.assertTrue(all(x == results[0] for x in results))

    def test_missing_mash_pairs_fail_instead_of_zero_fill(self):
        p = self.path / "mash.tsv"
        p.write_text("a.fa\ta.fa\t0\nb.fa\tb.fa\t0\n")
        result = self.script("mash_to_matrix.py", "--dist", p, "--out", self.path / "matrix.tsv")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing Mash pair", result.stderr)

    def test_kmers_never_cross_contig_boundaries_or_ambiguities(self):
        p = self.path / "seq.fa"
        p.write_text(">a\nAAAA\n>b\nCCCC\n")
        self.assertEqual(canonical_kmers(read_fasta_concat(p), 3), {"AAA", "CCC"})
        self.assertEqual(canonical_kmers("ARRA", 3), set())

    def test_invalid_fasta_and_duplicate_ids(self):
        for content in ("@raw\nACGT\n+\n!!!!\n", ">a\nACGT\n>a\nACGT\n", ">a\n", ">a\nACZ\n"):
            p = self.path / "bad.fa"
            p.write_text(content)
            with self.assertRaises(ValueError):
                fasta_records(p)

    def test_checksum_provenance_and_unknown_isolates(self):
        args = self.fixture()
        rows, seqs = validate_manifest(args.manifest, metadata(args.metadata), "longread", 3)
        self.assertEqual(rows[0]["quality_status"], "uncertain")
        rows[0]["sequence_sha256"] = "wrong"
        write_tsv(args.manifest, ["isolate_id", "plasmid_id", "fasta_path", "sequence_sha256"], rows)
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            validate_manifest(args.manifest, metadata(args.metadata), "precomputed", 3)

    def test_loose_arg_exclusion_and_coordinate_validation(self):
        args = self.fixture()
        rows, seqs = validate_manifest(args.manifest, metadata(args.metadata), "precomputed", 3)
        path = self.path / "features.tsv"
        feature = dict(plasmid_id="p0", feature_id="f1", start=1, end=4, strand="+", amr_gene="example", hit_class="loose", annotation_confidence="high", database_name="db", database_version="1", annotation_engine="caller")
        write_tsv(path, list(feature), [feature])
        normalized = normalize_features(path, rows, seqs, {"p0": "PU_0001"})
        self.assertEqual(normalized[0]["headline_eligible"], "false")
        feature["hit_class"] = "strict"
        write_tsv(path, list(feature), [feature])
        self.assertEqual(normalize_features(path, rows, seqs, {"p0": "PU_0001"})[0]["headline_eligible"], "true")
        feature["end"] = 100
        write_tsv(path, list(feature), [feature])
        with self.assertRaisesRegex(ValueError, "invalid coordinates"):
            normalize_features(path, rows, seqs, {})

    def test_metadata_html_is_escaped(self):
        args = self.fixture()
        write_tsv(args.metadata, ["isolate_id", "location"], [{"isolate_id": "iso1", "location": '</script><script>alert("x")</script>'}])
        run(args)
        content = (args.out / "REPORT.html").read_text(encoding="utf-8")
        self.assertNotIn('</script><script>alert("x")', content)
        self.assertIn('\\u003c/script', content)


if __name__ == "__main__":
    unittest.main()
