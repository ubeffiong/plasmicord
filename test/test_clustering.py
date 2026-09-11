#!/usr/bin/env python3
"""Unit test: clustering (single + complete linkage) on a hand-built matrix."""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
CLU = os.path.join(HERE, "..", "python", "cluster_plasmids.py")


def write_matrix(path, items, dist):
    with open(path, "w") as fh:
        fh.write("item\t" + "\t".join(items) + "\n")
        for i, it in enumerate(items):
            fh.write(it + "\t" + "\t".join(str(x) for x in dist[i]) + "\n")


def run(matrix, threshold, linkage, out, extra_args=()):
    subprocess.run([sys.executable, CLU, "--matrix", matrix,
                    "--threshold", str(threshold), "--linkage", linkage,
                    "--out", out, *extra_args], check=True)
    result = {}
    with open(out) as fh:
        fh.readline()
        for line in fh:
            pid, cid, *_ = line.rstrip("\n").split("\t")
            result[pid] = cid
    return result


def write_lengths(path, lengths):
    with open(path, "w") as fh:
        fh.write("plasmid_id\tlength\n")
        for pid, length in lengths.items():
            fh.write(f"{pid}\t{length}\n")


def margins(out):
    """plasmid_id -> (unit_max_internal_distance, unit_threshold_margin) as raw strings."""
    result = {}
    with open(out) as fh:
        fh.readline()
        for line in fh:
            pid, _cid, worst, margin = line.rstrip("\n").split("\t")
            result[pid] = (worst, margin)
    return result


def groups(res):
    """Return set of frozensets of items grouped by cluster id."""
    g = {}
    for item, cid in res.items():
        g.setdefault(cid, set()).add(item)
    return set(frozenset(v) for v in g.values())


def main():
    tmp = tempfile.mkdtemp(prefix="clu_test_")
    # 5 items. Distances designed so:
    #   A-B = 0.02, B-C = 0.03, A-C = 0.5 (A,B,C chain-link under single @0.05)
    #   D-E = 0.01 (own cluster). C-D = 0.4 (separate).
    items = ["A", "B", "C", "D", "E"]
    INF = 1.0
    d = [
        #  A     B     C     D     E
        [0.00, 0.02, 0.50, 0.60, 0.70],  # A
        [0.02, 0.00, 0.03, 0.55, 0.65],  # B
        [0.50, 0.03, 0.00, 0.40, 0.45],  # C
        [0.60, 0.55, 0.40, 0.00, 0.01],  # D
        [0.70, 0.65, 0.45, 0.01, 0.00],  # E
    ]
    m = os.path.join(tmp, "m.tsv")
    write_matrix(m, items, d)

    ok = True

    # SINGLE linkage @0.05: A-B (0.02) and B-C (0.03) link -> {A,B,C}; D-E -> {D,E}
    res = run(m, 0.05, "single", os.path.join(tmp, "s.tsv"))
    exp = {frozenset({"A", "B", "C"}), frozenset({"D", "E"})}
    got = groups(res)
    print("single@0.05 ->", got)
    ok &= (got == exp); print("  expected {A,B,C},{D,E} ->", got == exp)

    # Single linkage can chain beyond the threshold (A-C=0.50 > 0.05 within {A,B,C}):
    # the unit margin must go negative, exposing that diagnostic without touching clustering.
    m_single = margins(os.path.join(tmp, "s.tsv"))
    single_abc_margin = float(m_single["A"][1])
    ok &= single_abc_margin < 0
    print("  single {A,B,C} unit_threshold_margin < 0 ->", single_abc_margin < 0, f"({single_abc_margin})")
    single_de_margin = float(m_single["D"][1])
    ok &= abs(single_de_margin - 0.04) < 1e-6
    print("  single {D,E} unit_threshold_margin == 0.04 ->", abs(single_de_margin - 0.04) < 1e-6)

    # COMPLETE linkage @0.05: A-C is 0.50 > 0.05, so {A,B,C} cannot form under
    # complete linkage. A-B merge ok (0.02). C cannot join (needs C-A<=0.05).
    # => {A,B}, {C}, {D,E}
    res = run(m, 0.05, "complete", os.path.join(tmp, "c.tsv"))
    exp = {frozenset({"A", "B"}), frozenset({"C"}), frozenset({"D", "E"})}
    got = groups(res)
    print("complete@0.05 ->", got)
    ok &= (got == exp); print("  expected {A,B},{C},{D,E} ->", got == exp)

    # Complete linkage never produces a negative margin (algorithm invariant); a singleton
    # cluster (C) has no internal distance to report, so its margin fields are blank.
    m_complete = margins(os.path.join(tmp, "c.tsv"))
    complete_margins_nonneg = all(v != "" and float(v) >= 0 for _worst, v in m_complete.values() if v != "")
    ok &= complete_margins_nonneg
    print("  complete linkage unit_threshold_margin always >= 0 ->", complete_margins_nonneg)
    ok &= m_complete["C"] == ("", "")
    print("  complete {C} singleton margin blank ->", m_complete["C"] == ("", ""))

    # Threshold 0.0: everything singleton.
    res = run(m, 0.0, "single", os.path.join(tmp, "z.tsv"))
    got = groups(res)
    exp = {frozenset({x}) for x in items}
    print("single@0.0 -> all singletons:", got == exp)
    ok &= (got == exp)

    # Threshold 1.0: everything one cluster.
    res = run(m, 1.0, "single", os.path.join(tmp, "o.tsv"))
    got = groups(res)
    exp = {frozenset(items)}
    print("single@1.0 -> one cluster:", got == exp)
    ok &= (got == exp)

    # Determinism: same input -> identical cluster ids.
    r1 = run(m, 0.05, "single", os.path.join(tmp, "d1.tsv"))
    r2 = run(m, 0.05, "single", os.path.join(tmp, "d2.tsv"))
    ok &= (r1 == r2); print("deterministic ids:", r1 == r2)

    # Size-corrected threshold (opt-in via --lengths + --size-correction-per-percent):
    # F-G are two items whose raw distance (0.003) exceeds a base threshold of 0.001, but
    # whose 20% size difference loosens the effective threshold enough to merge them.
    items2 = ["F", "G"]
    d2 = [[0.00, 0.003], [0.003, 0.00]]
    m2 = os.path.join(tmp, "m2.tsv")
    write_matrix(m2, items2, d2)
    lengths_path = os.path.join(tmp, "lengths.tsv")
    write_lengths(lengths_path, {"F": 100, "G": 120})

    res = run(m2, 0.001, "single", os.path.join(tmp, "nocorr.tsv"))
    ok &= (groups(res) == {frozenset({"F"}), frozenset({"G"})})
    print("no size correction -> F,G stay separate:", groups(res) == {frozenset({"F"}), frozenset({"G"})})

    res = run(m2, 0.001, "single", os.path.join(tmp, "corr.tsv"),
              extra_args=["--lengths", lengths_path, "--size-correction-per-percent", "0.0003"])
    ok &= (groups(res) == {frozenset({"F", "G"})})
    print("size correction -> F,G merge:", groups(res) == {frozenset({"F", "G"})})

    # --size-correction-per-percent without --lengths is a clean CLI error, not a crash.
    failed_cleanly = subprocess.run([sys.executable, CLU, "--matrix", m2, "--threshold", "0.001",
                                     "--linkage", "single", "--out", os.path.join(tmp, "bad.tsv"),
                                     "--size-correction-per-percent", "0.0003"],
                                    capture_output=True, text=True).returncode != 0
    ok &= failed_cleanly
    print("missing --lengths with correction enabled -> clean CLI error:", failed_cleanly)

    # A zero/negative/non-numeric length must fail loudly, not silently fall back to the
    # uncorrected threshold for just that item's pairs (a previously-audited silent-fallback risk).
    bad_lengths_path = os.path.join(tmp, "bad_lengths.tsv")
    write_lengths(bad_lengths_path, {"F": 100, "G": 0})
    bad_length_failed = subprocess.run([sys.executable, CLU, "--matrix", m2, "--threshold", "0.001",
                                        "--linkage", "single", "--out", os.path.join(tmp, "bad2.tsv"),
                                        "--lengths", bad_lengths_path, "--size-correction-per-percent", "0.0003"],
                                       capture_output=True, text=True).returncode != 0
    ok &= bad_length_failed
    print("zero length in --lengths file -> clean CLI error, not silent fallback:", bad_length_failed)

    print("\nCLUSTERING TESTS PASSED" if ok else "\nCLUSTERING TESTS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
