# Reading the final output

Start with **REPORT.html**. Check the dataset label, accepted/rejected candidates,
typing coverage and annotation status before interpreting the counts.

1. Inspect the network and cross-cluster table. A sharing pair means candidates from
   two isolates belong to the same sequence-defined unit. Unknown chromosome labels
   cannot establish chromosome discordance.
2. Select an edge. Inspect shared units, minimum distances and direct threshold
   support; cluster membership alone may reflect single-linkage chaining.
3. Filter on an eligible ARG, then open the functional cargo section. Inspect the
   original gene coordinates, source/database/version and candidate quality warnings.
4. Review locations and dates as context. Neither the collection order nor multiple
   links establishes direction, active movement or a direct epidemiological chain.
5. Compare threshold sensitivity. Confirm prioritized candidates with suitable
   sequence/assembly evidence and epidemiological review.

Keep **REPORT.html** for exploration, **REPORT.md** for a narrative, TSV/JSON for
analysis, and **network.graphml** for Cytoscape or Gephi. A browser-printed PDF is
a static snapshot; the HTML preserves interactive exploration.

Zero links does not prove no transmission. No candidate for an isolate does not prove
it is plasmid-free. Missing annotations are not gene absence. ARG presence does not
prove expression or phenotypic resistance. A metabolic annotation does not establish
a complete pathway. Read the report's safeguards alongside every result.
