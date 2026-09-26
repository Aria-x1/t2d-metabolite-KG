# T2D Metabolomics → Knowledge Graph Pipeline

A repeatable pipeline that takes statistically significant metabolite–disease
associations from metabolomics studies and turns them into
[Biolink Model](https://biolink.github.io/biolink-model/)-compliant triples,
ready for submission to the [ROBOKOP](https://robokop.renci.org/) biomedical
knowledge graph — with automatic duplicate-checking against the live graph.

The resulting triples are intended to feed a downstream drug-target discovery
pipeline (TARRAGON → SABLE), as part of an NCDRC pilot-grant project studying
Type 2 Diabetes metabolomics.

## Why this exists

Manually converting a table of "metabolite X is associated with disease Y"
findings into a knowledge-graph-ready format involves a lot of repetitive,
error-prone work: normalizing free-text names to standard ontology
identifiers (CHEBI for metabolites, MONDO for diseases), picking the correct
Biolink predicate for the type of statistical evidence at hand, and checking
whether the claim is actually new to the graph or a duplicate of something
already there. This pipeline automates all three steps and caches API calls
so it scales across many datasets (this project processes findings from
~80–100 studies).

## Pipeline

```
significant_metabolites_BH_vs_Holm.xlsx
        │  convert_holm_to_pipeline_input.py
        ▼
   pipeline_input.csv   (metabolite_name, disease_name, direction, p_value, evidence_type, ...)
        │  build_robokop_triples.py
        ▼
   robokop_triples.csv  (Biolink-compliant subject/predicate/object triples,
                          each flagged new_node / new_edge / possible_duplicate)
```

**`convert_holm_to_pipeline_input.py`** — reshapes a study-specific results
table (Holm-Bonferroni–corrected significant metabolites from an MTBLS8644
group comparison) into the pipeline's standard input format. This is the
"study-specific adapter" layer — each new dataset gets its own small
converter script like this one, feeding a common downstream pipeline.

**`build_robokop_triples.py`** — the general pipeline, dataset-agnostic:

1. **Entity normalization** — resolves free-text metabolite/disease names to
   canonical CURIEs via the NCATS Name Resolution and NodeNormalization APIs.
2. **Predicate assignment** — maps `(evidence_type, direction)` to a Biolink
   predicate, e.g. a group-comparison result becomes
   `biolink:associated_with_increased_likelihood_of`, while a continuous
   correlation becomes `biolink:positively_correlated_with`.
3. **Duplicate checking** — queries the live ROBOKOP graph (Automat Cypher
   endpoint) to check whether the subject node already exists and whether an
   edge already connects it to the object, so every row is labeled
   `new_node`, `new_edge`, `new_edge_type`, or `possible_duplicate`.
4. **Caching** — normalization and existence-check results are cached to
   disk (`curie_cache.json`, `existence_cache.json`) so re-running across many
   datasets doesn't repeatedly hit the same public APIs for the same entity.

## Data sources

- [MTBLS8644](https://www.ebi.ac.uk/metabolights/MTBLS8644) (MetaboLights) —
  cross-sectional cohort, 80 samples across four groups (Healthy / early
  prediabetes / late prediabetes / T2DM).
- [iPOP](https://www.nature.com/articles/s41591-019-0509-0) (Zhou et al.,
  *Nature Medicine* 2019) — longitudinal validation cohort, used as an
  independent external check on findings derived from MTBLS8644.

Raw data is not included in this repo — both datasets are publicly available
at the links above.

## Key methodological notes

- **Circular validation**: a signature's performance must be checked on data
  it wasn't derived from. Testing on the same cohort it was built on
  overstates accuracy and is called out explicitly wherever it applies.
- **CURIEs are always resolved via the normalization APIs**, never assigned
  by hand, to keep identifiers consistent with the rest of the ROBOKOP graph.
- **Predicate choice follows evidence type**, not just direction — a
  cross-sectional group comparison and a continuous correlation with a
  clinical variable (e.g. A1C) support different, specifically-worded biolink
  claims, and the pipeline keeps that distinction explicit rather than
  collapsing everything into a generic "associated_with".

## Usage

```bash
pip install -r requirements.txt

# 1. Convert a study-specific results table into the standard pipeline input
python scripts/convert_holm_to_pipeline_input.py significant_metabolites_BH_vs_Holm.xlsx pipeline_input.csv

# 2. Build Biolink-compliant, duplicate-checked triples
python scripts/build_robokop_triples.py pipeline_input.csv robokop_triples.csv
```

Both scripts call free, public, unauthenticated NCATS Translator / ROBOKOP
endpoints — run them from an environment with normal outbound internet
access.

## Status

This pipeline currently handles the metabolite–disease triples derived from
MTBLS8644 group comparisons. Ongoing work extends it to numeric-correlation
findings from the iPOP longitudinal cohort, and formalizes the delivery
format for submitting the resulting triple table to the ROBOKOP graph
maintainers.
