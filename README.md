# T2D Metabolomics → Knowledge Graph Pipeline

Turns statistically significant metabolite–disease associations into
[Biolink Model](https://biolink.github.io/biolink-model/)-compliant triples
for submission to [ROBOKOP](https://robokop.renci.org/), with automatic
duplicate-checking against the live graph. Part of an NCDRC pilot-grant
project studying Type 2 Diabetes metabolomics; feeds a downstream drug-target
discovery pipeline (TARRAGON → SABLE).

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



**`convert_holm_to_pipeline_input.py`** — study-specific adapter: reshapes an
MTBLS8644 results table into the pipeline's standard input format. Each new
dataset gets its own small adapter like this one.

**`build_robokop_triples.py`** — general, dataset-agnostic pipeline:

1. **Normalization** — free-text metabolite/disease names → canonical CURIEs
   via NCATS Name Resolution + NodeNormalization APIs.
2. **Predicate assignment** — `(evidence_type, direction)` → Biolink
   predicate (e.g. `associated_with_increased_likelihood_of` vs.
   `positively_correlated_with`).
3. **Duplicate checking** — queries the live ROBOKOP graph (Automat Cypher)
   to flag each row `new_node` / `new_edge` / `new_edge_type` /
   `possible_duplicate`.
4. **Caching** — normalization + existence checks cached locally
   (`curie_cache.json`, `existence_cache.json`) to scale across ~80–100
   datasets.

## Data sources

- [MTBLS8644](https://www.ebi.ac.uk/metabolights/MTBLS8644) — cross-sectional,
  80 samples, 4 groups (Healthy → T2DM).
- [iPOP](https://www.nature.com/articles/s41591-019-0509-0) (Zhou et al.,
  *Nature Medicine* 2019) — longitudinal external validation cohort.

Raw data not included — both are public at the links above.

## Usage

```bash
pip install -r requirements.txt

python scripts/convert_holm_to_pipeline_input.py significant_metabolites_BH_vs_Holm.xlsx pipeline_input.csv
python scripts/build_robokop_triples.py pipeline_input.csv robokop_triples.csv
```

## Status

**Actively in progress.** Handles MTBLS8644 group-comparison and iPOP
correlation findings so far. Currently: incorporating additional data
sources, integrating all findings into one deduplicated triple table, and
finalizing the delivery format to ROBOKOP maintainers.
