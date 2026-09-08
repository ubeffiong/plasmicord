# PlasmiCord: continuing beyond reconstruction benchmarking

**PlasmiCord: A Chromosome-Aware Plasmid Transmission and Surveillance Framework**

*Connecting plasmid sharing with chromosomal epidemiology.*

PlasBench evaluates reconstruction against reference truth and shows protein-coordinate recovery. PlasmiCord continues at cohort level: sequence-defined plasmid units, isolate sharing, chromosome-cluster discordance and functional cargo across available metadata.

| Question | PlasBench | PlasmiCord |
|---|---|---|
| Which method recovered a reference accurately? | Primary | Upstream context |
| Were genes and reference coordinates recovered? | Reconstruction evidence | Source context |
| Which isolates carry related plasmid candidates? | Outside main scope | Units and sharing network |
| Do units connect distinct chromosome clusters? | Outside main scope | Discordance evidence |
| What cargo occurs across isolates? | Supporting annotation | Main functional report |
| Is direct transmission proven? | No | No |

Lessons retained are inspectable coordinates, offline evidence reports, versioned provenance, explicit unevaluated states and reusable tables. No PlasBench application, import, execution or sibling directory is required.

The common boundary is candidate FASTA, isolate mapping and available provenance. Native adapters now support MOB-recon, Flye, Unicycler and selected PlasBench outputs. Curated references and other tools can supply the generic contract. Assembly and polishing remain upstream.

The separately buildable [plasmid-annotation-core](../shared_annotation/README.md) supplies versioned local annotation, normalized coordinates, an auditable cache and a PlasBench-compatible protein export. Both applications can consume this independent component. PlasBench itself has not been modified to use it automatically, and its recovery viewer remains separate.

This release adds biological evidence rules and public-reference threshold calibration. Next scientific steps include independently reviewed transmission cohorts, comparable chromosome typing, empirical quality-tier validation, structural confirmation, orthology/module analysis and surveillance-scale performance evaluation. See [the review](REVIEW.md).
