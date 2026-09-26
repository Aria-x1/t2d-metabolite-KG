"""
convert_holm_to_pipeline_input.py

One-off converter: reads Significant_Holm sheet from
significant_metabolites_BH_vs_Holm.xlsx and reshapes it into the input CSV
format expected by build_robokop_triples.py.

Column mapping:
    metabolite_identification -> metabolite_name
    comparison                -> determines disease_name / disease_curie
                                  (see COMPARISON_TO_DISEASE below)
    direction ("Up"/"Down")   -> direction ("increased"/"decreased")
    p_value_holm              -> p_value  (the Holm-corrected value, not raw p)
    (fixed)                   -> test_type, evidence_type, dataset_id

IMPORTANT: "Healthy vs Pre-diabetic" rows are mapped to Prediabetes
(MONDO:0006920), NOT Type 2 Diabetes Mellitus -- these two are biologically
distinct comparisons and treating them as equivalent would overstate the
evidence. Confirm this mapping is what you intend before running the full
batch through build_robokop_triples.py.

USAGE:
    python convert_holm_to_pipeline_input.py significant_metabolites_BH_vs_Holm.xlsx pipeline_input.csv
"""

import sys
import pandas as pd

COMPARISON_TO_DISEASE = {
    "Healthy vs T2DM": ("Type 2 Diabetes Mellitus", "MONDO:0005148"),
    "Pre-diabetic vs T2DM": ("Type 2 Diabetes Mellitus", "MONDO:0005148"),
    "Healthy vs Pre-diabetic": ("Prediabetes", "MONDO:0006920"),
}

DIRECTION_MAP = {"Up": "increased", "Down": "decreased"}

TEST_TYPE = "Mann-Whitney U"  # confirm this matches what was actually run
EVIDENCE_TYPE = "group_comparison"
DATASET_ID = "MTBLS8644"


def convert(input_xlsx, output_csv):
    df = pd.read_excel(input_xlsx, sheet_name="Significant_Holm")

    unknown_comparisons = set(df["comparison"]) - set(COMPARISON_TO_DISEASE)
    if unknown_comparisons:
        print(
            f"[error] Unrecognized comparison value(s): {unknown_comparisons}. "
            f"Add them to COMPARISON_TO_DISEASE before proceeding.",
            file=sys.stderr,
        )
        sys.exit(1)

    rows = []
    for _, row in df.iterrows():
        disease_name, disease_curie = COMPARISON_TO_DISEASE[row["comparison"]]
        rows.append({
            "metabolite_name": row["metabolite_identification"],
            "disease_name": disease_name,
            "disease_curie": disease_curie,
            "direction": DIRECTION_MAP[row["direction"]],
            "p_value": row["p_value_holm"],
            "test_type": TEST_TYPE,
            "evidence_type": EVIDENCE_TYPE,
            "dataset_id": DATASET_ID,
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(output_csv, index=False)
    print(f"Wrote {len(out_df)} rows to {output_csv}")
    print("\nBreakdown by disease object:")
    print(out_df["disease_name"].value_counts().to_string())


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_holm_to_pipeline_input.py input.xlsx output.csv")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
