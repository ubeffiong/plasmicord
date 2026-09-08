"""Scientific reporting and portable-output regressions."""
import argparse
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from python import PROJECT_TITLE, PROJECT_TAGLINE
from python.contracts import write_tsv
from python.report_evidence import summary, catalog, sha256, result_files
from python.report_output import render_dashboard, refresh_report, attach_calibration
from python.cli import run

class DetailedReportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.out=self.root/"results";self.out.mkdir()
    def tearDown(self): self.tmp.cleanup()
    def payload(self):
        ps=[dict(plasmid_id=p,isolate_id=i,plasmid_unit="PU_0001",length=1000,n_contigs=1,quality_status="uncertain",quality_warnings="identity unresolved",source_tool="fixture") for p,i in [("p1","i1"),("p2","i2")]]
        meta=[dict(isolate_id=i,chromosomal_cluster="",date="",location="",organism="") for i in ["i1","i2"]]
        write_tsv(self.out/"metadata.tsv",list(meta[0]),meta)
        (self.out/"plasmid_matrix.tsv").write_text("item\tp1\tp2\np1\t0\t0.01\np2\t0.01\t0\n")
        return dict(provenance=dict(project_title=PROJECT_TITLE,project_tagline=PROJECT_TAGLINE,status="complete",
            framework_version="0.4.0",dataset_kind="synthetic demonstration",distance_engine="kmer",k=3,sketch_size=None,linkage="complete",threshold=.02,annotation_status="imported"),
            plasmids=ps,metadata=meta,features=[dict(plasmid_id="p2",isolate_id="i2",headline_eligible="true",amr_gene="ARGfixture",functional_category="amr")],
            annotation_status=[dict(plasmid_id="p1",stage="amr",status="complete",engine="amrfinder")],biological_quality=[],unit_functions=[],
            edges=[dict(source="i1",target="i2",weight=1,shared_units="PU_0001")],
            edge_evidence=[dict(source="i1",target="i2",plasmid_unit="PU_0001",direct_threshold_support="true",minimum_distance=.01,shared_args="")],
            crosslinks=[],safeguards=["Research evidence only."],contigs={"p1":[dict(id="c1",length=1000)],"p2":[dict(id="c2",length=1000)]},
            sensitivity=[dict(threshold=.02,n_units=1,n_singletons=0,linkage="complete")],
            counts=dict(isolates=2,plasmids=2,units=1,sharing_pairs=1,cross_cluster_pairs=0,cross_cluster_unit_links=0,typed_isolates=0,observed_arg_plasmids=1))
    def test_rules_distinguish_unknown_typing_imports_and_completed_zero_hits(self):
        result=summary(self.payload());findings={f["rule_id"]:f for f in result["findings"]}
        self.assertEqual(findings["typing_coverage"]["state"],"not_evaluated")
        self.assertEqual(findings["amr_coverage"]["state"],"warn")
        self.assertIn("1/2",findings["amr_coverage"]["interpretation"])
        self.assertIn("1 completed candidates",findings["amr_coverage"]["interpretation"])
        self.assertEqual(result["samples"][0]["amr_completed"],1)
        self.assertEqual(result["samples"][1]["amr_completed"],0)
        self.assertEqual(result["charts"]["functions"],[dict(label="amr",value=1)])
    def test_nested_inventory_escaping_previews_and_bundle_checksums(self):
        data=self.payload()
        nested=self.out/"annotation/p1/amr";nested.mkdir(parents=True)
        hostile="</script><script>alert(1)</script> __REPORT_DATA__"
        (nested/"a # log.txt").write_text(hostile+"x"*20000)
        (nested/"raw.bin").write_bytes(b"\0\xff\x00")
        render_dashboard(self.out,data)
        page=(self.out/"REPORT.html").read_text(encoding="utf-8")
        self.assertNotIn("</script><script>alert(1)</script>",page)
        self.assertIn("\\u003c/script",page)
        files={f["path"]:f for f in data["artifacts"]}
        f=files["annotation/p1/amr/a # log.txt"]
        self.assertIn("a%20%23%20log.txt",f["href"])
        self.assertTrue(f["preview_truncated"])
        self.assertIsNone(files["annotation/p1/amr/raw.bin"]["preview"])
        manifest=json.loads((self.out/"output_manifest.json").read_text())
        for row in manifest["files"]:
            self.assertEqual(sha256(self.out/row["path"]),row["sha256"])
        expected={p.relative_to(self.out).as_posix() for p in result_files(self.out)}-{"REPORT_BUNDLE.zip"}
        with zipfile.ZipFile(self.out/"REPORT_BUNDLE.zip") as archive:
            self.assertEqual(set(archive.namelist()),expected)
            self.assertEqual(archive.read("annotation/p1/amr/raw.bin"),b"\0\xff\x00")
        self.assertIsNone(files["run_provenance.json"]["sha256"])
    def test_file_catalog_does_not_follow_external_symlinks(self):
        outside=self.root/"private.txt";outside.write_text("outside")
        try: (self.out/"external.txt").symlink_to(outside)
        except OSError: self.skipTest("Symlink creation unavailable")
        self.assertNotIn("external.txt",{f["path"] for f in catalog(self.out)})
    def test_report_refresh_preserves_analysis_and_removes_stale_bundle_when_disabled(self):
        data=self.payload();render_dashboard(self.out,data)
        before=sha256(self.out/"plasmid_matrix.tsv")
        refresh_report(argparse.Namespace(results=self.out,calibration=None,no_bundle=True))
        self.assertEqual(before,sha256(self.out/"plasmid_matrix.tsv"))
        self.assertFalse((self.out/"REPORT_BUNDLE.zip").exists())
        result=json.loads((self.out/"report_data.json").read_text())
        self.assertEqual(result["report"]["renderer_version"],"0.4.0")
        self.assertFalse(result["report"]["bundle"])
    def test_refresh_rejects_changed_analysis_and_payload(self):
        data=self.payload();render_dashboard(self.out,data,bundle=False)
        (self.out/"metadata.tsv").write_text("changed")
        with self.assertRaisesRegex(ValueError,"changed or missing"):
            refresh_report(argparse.Namespace(results=self.out,calibration=None,no_bundle=True))
    def test_calibration_attachment_checks_cohort_and_linkage(self):
        data=self.payload();source=self.root/"calibration";source.mkdir()
        record=dict(matrix_sha256=sha256(self.out/"plasmid_matrix.tsv"),distance_engine="kmer",linkage="single",
                    status="insufficient_independent_labels",selected_threshold=None,target="plasmid_sequence_relatedness")
        (source/"calibration.json").write_text(json.dumps(record))
        for name in ("calibration_metrics.tsv","calibration_labels.tsv"):(source/name).write_text("split\tthreshold\n")
        with self.assertRaisesRegex(ValueError,"engine/linkage"):
            attach_calibration(self.out,source,data["provenance"])
        record["linkage"]="complete";(source/"calibration.json").write_text(json.dumps(record))
        attached=attach_calibration(self.out,source,data["provenance"])
        self.assertEqual(attached["record"]["status"],"insufficient_independent_labels")
        record["matrix_sha256"]="0"*64;(source/"calibration.json").write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError,"matrix does not match"):
            attach_calibration(self.out,source,data["provenance"])
    def test_annotation_failure_emits_diagnostic_dashboard_without_completed_counts(self):
        meta=self.root/"meta.tsv";write_tsv(meta,["isolate_id"],[dict(isolate_id="i1")])
        fasta=self.root/"p.fa";fasta.write_text(">p\nACGTACGTACGT\n")
        manifest=self.root/"manifest.tsv";write_tsv(manifest,["isolate_id","plasmid_id","fasta_path"],[dict(isolate_id="i1",plasmid_id="p",fasta_path=str(fasta))])
        args=argparse.Namespace(manifest=manifest,metadata=meta,features=None,mode="precomputed",engine="kmer",
            out=self.root/"failed",threshold=.01,k=3,min_length=3,sketch_size=100,linkage="complete",annotation_config=self.root/"config.json")
        with patch("python.annotation.annotate_candidates",side_effect=ValueError("annotation stage unavailable")):
            with self.assertRaisesRegex(ValueError,"annotation stage unavailable"):run(args)
        payload=json.loads((args.out/"report_data.json").read_text())
        self.assertEqual(payload["provenance"]["status"],"failed")
        self.assertIsNone(payload["counts"]["plasmids"])
        self.assertIsNone(payload["dashboard"]["samples"][0]["accepted_candidates"])
        self.assertEqual(payload["dashboard"]["samples"][0]["rejected_candidates"],0)
        self.assertEqual(payload["dashboard"]["findings"][0]["state"],"fail")
        self.assertTrue((args.out/"REPORT.html").exists())
        self.assertNotIn("report_error",json.loads((args.out/"run_provenance.json").read_text()))
if __name__=="__main__":unittest.main()
