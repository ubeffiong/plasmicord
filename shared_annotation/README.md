# plasmid-annotation-core

This separately buildable Python package is the common functional annotation component used by PlasmiCord. It imports neither PlasmiCord nor PlasBench and can be installed by either project.

From the repository root:

```sh
python -m pip wheel --no-deps --wheel-dir dist ./shared_annotation
python -m pip install dist/plasmid_annotation_core-0.1.0-py3-none-any.whl
python -m plasmid_annotation --fasta candidate.fasta --config annotation.json --out annotation_output
```

The standalone command writes features.json, provenance.json and plasbench_proteins.tsv. The Python API exports `annotate`, `parse_gff`, `parse_amrfinder` and `to_plasbench`. See [annotation configuration](../docs/ANNOTATION.md).

The source lives in ../plasmid_annotation and is also bundled with the main project wheel. Install either distribution in an environment, rather than installing two separately versioned copies over the same package files. No network database lookup is performed during analysis. Tools and versioned databases are explicit external requirements.

PlasBench's application has not been changed to consume this library automatically. It can adopt the standalone wheel/API or consume the exported protein TSV; this is a reusable component and compatible boundary, not a dependency on either application. The existing PlasmiCord viewer remains project-specific.
