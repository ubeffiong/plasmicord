# Recommendations for PlasBench (handoff note, not a PlasmiCord feature)

PlasBench is a separate project (reconstruction-quality benchmarking: "how trustworthy is
this reconstruction?"). PlasmiCord does not install, import, call or depend on it -- see
[project positioning](POSITIONING.md). This note exists only because the same competitor
research done for PlasmiCord's own roadmap (see [REVIEW.md](REVIEW.md)) surfaced items
relevant to PlasBench's mission instead. **Nothing here is implemented by PlasmiCord.**
It is a handoff for whoever maintains PlasBench to evaluate independently.

## Relevant lessons

### 1. Nearest-neighbor reconstruction-identity reporting (from MOB-suite)

**Source:** [phac-nml/mob-suite](https://github.com/phac-nml/mob-suite) (Apache-2.0).
**Why it matters for a reconstruction benchmark:** `mob_cluster` reports
`mash_nearest_neighbor`, `mash_neighbor_distance` and `mash_neighbor_identification` per
plasmid -- "how close is this reconstruction to the nearest known reference, and which
reference" is itself a trustworthiness signal distinct from PlasBench's existing
reference-truth comparison. Worth evaluating as an additional per-assembly metric.

### 2. Alignment-based containment/completeness check (from MobMess)

**Source:** [michaelkyu/MobMess](https://github.com/michaelkyu/MobMess) (GPL-3.0 --
note the copyleft license before reusing code directly; the *approach* is reusable
regardless).
**Why it matters:** MobMess uses MUMmer4 pairwise alignment (coverage ≥90%, identity ≥90%)
to decide whether one plasmid is a near-total subset of another. For a reconstruction
benchmark, an analogous alignment-based check against reference truth could give a
structural-completeness score more rigorous than a heuristic -- a concrete alternative or
complement to whatever completeness measure PlasBench currently uses.

### 3. ANI-based reference-network insertion with an explicit novel-vs-known decision (from COPLA)

**Source:** [santirdnd/COPLA](https://github.com/santirdnd/COPLA) (GPL-3.0).
**Why it matters:** COPLA computes ANI between a query and a reference plasmid network,
statistically tests cluster membership, and either assigns a known taxonomic unit or
flags the query as novel, with a confidence score. The "explicit novel-vs-known decision
with confidence" pattern maps directly onto PlasBench's "how trustworthy" framing --
worth considering for how PlasBench communicates uncertain/out-of-reference-set cases.

## Not relevant to PlasBench

Two repos from the same research are flagged here only to say they don't apply:
**phac-nml/plasmid_analysis** (population-surveillance aggregation) and **mge-cluster**
(reference-free clustering) are both closer to PlasmiCord's isolate/cohort-level mission
than to PlasBench's per-assembly reconstruction-quality mission. See PlasmiCord's own
[REVIEW.md](REVIEW.md) for what PlasmiCord took from them instead.

## How to use this note

Each item above is a recommendation to evaluate, not a description of work already done
in either project. If PlasBench's maintainer wants to pursue any of these, treat this
note as a starting pointer to the source project and its license, not as a specification.
