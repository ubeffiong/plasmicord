#!/usr/bin/env python3
"""
Unit test: build_network.py + discordance.py on a hand-built outbreak scenario.

Scenario (5 isolates):
  Chromosomal clusters (from cgMLST/SNP, given in metadata):
      iso1, iso2  -> chrom cluster CC1
      iso3, iso4  -> chrom cluster CC2
      iso5        -> chrom cluster CC3   (plasmid-free)
  Plasmid units carried:
      iso1: PU_0001
      iso2: PU_0001, PU_0002
      iso3: PU_0001            <-- shares PU_0001 with CC1 isolates!
      iso4: PU_0002            <-- shares PU_0002 with iso2 (CC1) across clusters
      iso5: (none)

Expected network edges (share >=1 PU):
   iso1-iso2 (PU_0001)              weight 1
   iso1-iso3 (PU_0001)              weight 1   CROSS (CC1 vs CC2)
   iso2-iso3 (PU_0001)              weight 1   CROSS (CC1 vs CC2)
   iso2-iso4 (PU_0002)              weight 1   CROSS (CC1 vs CC2)
  (iso4 has PU_0002; iso2 has PU_0002 -> edge; iso3 has PU_0001 only)
   => 4 edges total. iso5 has no edges.

Expected cross-links (share PU but different chrom cluster):
   PU_0001: iso1(CC1)-iso3(CC2), iso2(CC1)-iso3(CC2)
   PU_0002: iso2(CC1)-iso4(CC2)
   => 3 cross-cluster pairs; 2 PUs cross clusters.
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
NET = os.path.join(HERE, "..", "python", "build_network.py")
DIS = os.path.join(HERE, "..", "python", "discordance.py")


def main():
    tmp = tempfile.mkdtemp(prefix="net_test_")

    index = os.path.join(tmp, "plasmid_index.tsv")
    with open(index, "w") as fh:
        fh.write("plasmid_id\tisolate_id\tlength\n")
        rows = [
            ("iso1__p1", "iso1"),
            ("iso2__p1", "iso2"),
            ("iso2__p2", "iso2"),
            ("iso3__p1", "iso3"),
            ("iso4__p1", "iso4"),
        ]
        for pid, iso in rows:
            fh.write(f"{pid}\t{iso}\t1000\n")

    clusters = os.path.join(tmp, "plasmid_clusters.tsv")
    with open(clusters, "w") as fh:
        fh.write("plasmid_id\tplasmid_unit\n")
        assign = {
            "iso1__p1": "PU_0001",
            "iso2__p1": "PU_0001",
            "iso2__p2": "PU_0002",
            "iso3__p1": "PU_0001",
            "iso4__p1": "PU_0002",
        }
        for pid, pu in assign.items():
            fh.write(f"{pid}\t{pu}\n")

    meta = os.path.join(tmp, "metadata.tsv")
    with open(meta, "w") as fh:
        fh.write("isolate_id\tchromosomal_cluster\tlocation\n")
        for iso, cc in [("iso1", "CC1"), ("iso2", "CC1"),
                        ("iso3", "CC2"), ("iso4", "CC2"), ("iso5", "CC3")]:
            fh.write(f"{iso}\t{cc}\tsite_{cc}\n")

    prefix = os.path.join(tmp, "net")
    subprocess.run([sys.executable, NET, "--index", index, "--clusters", clusters,
                    "--metadata", meta, "--out-prefix", prefix], check=True)

    # ---- check edges ----
    edges = set()
    with open(prefix + ".edges.tsv") as fh:
        fh.readline()
        for line in fh:
            s, t, w, shared = line.rstrip("\n").split("\t")
            edges.add((s, t, int(w)))
    print("edges:", sorted(edges))
    exp_edges = {("iso1", "iso2", 1), ("iso1", "iso3", 1),
                 ("iso2", "iso3", 1), ("iso2", "iso4", 1)}
    ok = (edges == exp_edges)
    print("  expected 4 edges ->", ok)

    # iso5 must still appear as a node in graphml (plasmid-free)
    with open(prefix + ".graphml") as fh:
        gml = fh.read()
    iso5_node = ('id="iso5"' in gml)
    print("  iso5 present as node ->", iso5_node)
    ok &= iso5_node

    # ---- discordance ----
    dprefix = os.path.join(tmp, "disc")
    subprocess.run([sys.executable, DIS, "--isolate-units", prefix + ".isolate_units.tsv",
                    "--metadata", meta, "--out-prefix", dprefix], check=True)
    crosspairs = set()
    with open(dprefix + ".crosslinks.tsv") as fh:
        fh.readline()
        for line in fh:
            pu, a, ca, b, cb = line.rstrip("\n").split("\t")
            crosspairs.add((pu, a, b))
    print("crosslinks:", sorted(crosspairs))
    exp_cross = {("PU_0001", "iso1", "iso3"),
                 ("PU_0001", "iso2", "iso3"),
                 ("PU_0002", "iso2", "iso4")}
    ok &= (crosspairs == exp_cross)
    print("  expected 3 cross-cluster pairs ->", crosspairs == exp_cross)

    with open(dprefix + ".summary.txt") as fh:
        summary = fh.read()
    print("\n--- summary.txt ---\n" + summary)
    ok &= ("PUs crossing chromosomal clusters .......... 2" in summary)

    print("NETWORK+DISCORDANCE TESTS PASSED" if ok else "NETWORK+DISCORDANCE TESTS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
