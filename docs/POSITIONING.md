# PlasmiCord: from PlasBench to plasmid-sharing investigation

**PlasmiCord: A Chromosome-Aware Plasmid Transmission and Surveillance Framework**

*Connecting plasmid sharing with chromosomal epidemiology.*

PlasBench v0.2.7, reviewed locally for this work, benchmarks reconstruction tools
against reference truth and offers candidate evaluation and protein-coordinate
recovery views. Project 2 continues the scientific programme at cohort level.

| Question | PlasBench | PlasmiCord |
|---|---|---|
| Which method reconstructed a reference most accurately? | Primary | Out of scope |
| Were ARGs/proteins and their reference context recovered? | Reconstruction evidence | Imported quality context only |
| Which isolates carry related plasmid candidates? | Out of scope | Sequence-defined units and network |
| Which units link distinct chromosome clusters? | Out of scope | Discordance evidence |
| What functional cargo occurs across isolates, places and dates? | Supporting annotation | Primary operational report |
| Is direct person-to-person transmission proven? | No | No |

Lessons carried forward: an offline evidence report, explicit unmeasured states,
versioned provenance, inspectable gene coordinates, reproducible synthetic fixtures,
downloadable tables and restrained interpretation. No PlasBench code, assets, datasets,
dependencies or sibling filesystem paths are required at runtime.

The stable boundary is **plasmid candidate FASTA + isolate mapping + available
provenance**, optionally accompanied by normalized annotations. MOB-recon, Platon,
plasmidSPAdes, Flye, Unicycler, curated references and future workflows can supply it.
Native directory adapters are a separate future convenience layer. Assembled ONT,
PacBio and hybrid candidates use this same boundary; raw reads require upstream QC,
assembly, polishing and plasmid identification.

The pasted design calls for a shared annotation engine/viewer across projects. This
release provides an independent normalized functional contract and Project 2 viewer;
it does not yet extract PlasBench's annotation code into a shared versioned package.
That should be a separately tested reusable dependency, not a runtime dependency on
the PlasBench application. PlasBench itself was read for comparison and not modified.

## Next-level scientific work

1. Common evidence-based biological gate: classification, contamination, circularity,
   read/assembly-graph support and source-specific quality evidence.
2. Versioned gene/AMR/replicon/mobility annotation profiles and cache keys incorporating
   sequence, engine version, database version and parameters.
3. Alignment-backed confirmation of prioritized links, containment/structural changes,
   comparable chromosome typing and epidemiological review.
4. Orthology, conserved/accessory cargo and module completeness from defined required
   components, with explicit annotation-coverage denominators.
5. Real-cohort validation, threshold calibration, performance measurement and external
   reproducibility before surveillance deployment.
