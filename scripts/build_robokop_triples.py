"""
build_robokop_triples.py

Repeatable workflow: take a table of statistically significant metabolite-disease
associations and produce a Biolink-Model-compliant triple table ready for
downstream submission to ROBOKOP, with automatic duplicate-checking against
the live ROBOKOP KG.

Pipeline steps per row:
    1. Resolve the metabolite name to a standardized CURIE (NCATS NameResolution
       + NodeNormalization APIs).
    2. Assign a Biolink predicate based on evidence_type + direction:
         - group_comparison    -> associated_with_increased/decreased_likelihood_of
         - numeric_correlation -> positively/negatively_correlated_with
    3. Query the live ROBOKOP KG (Automat Cypher endpoint) to check:
         - does the subject node already exist in the graph?
         - does an edge between subject and object already exist, and with
           which predicate(s)?
    4. Flag rows as "new_node", "new_edge", or "possible_duplicate" so you can
       see at a glance which rows represent a genuinely new contribution.
    5. Cache both normalization lookups and existence checks locally so
       re-running across 80-100 datasets does not re-query the same
       metabolite or edge repeatedly.

INPUT CSV format (one row per significant metabolite-disease finding):
    metabolite_name, disease_name, disease_curie, direction, p_value,
    test_type, evidence_type, dataset_id

    - disease_curie is optional; if blank, disease_name is resolved the same
      way as metabolite_name.
    - direction: "increased" or "decreased".
    - evidence_type: "group_comparison" (e.g. Mann-Whitney U / trend test
      between disease status groups) or "numeric_correlation" (e.g. Pearson/
      Spearman correlation with a continuous variable like A1C). This
      determines which predicate family gets assigned -- see PREDICATE_MAP.

OUTPUT CSV format (Biolink-compliant triple table):
    subject_id, subject_name, subject_category, predicate,
    object_id, object_name, object_category,
    direction, p_value, evidence, knowledge_source, dataset_id,
    normalization_status, subject_node_exists, existing_edge_predicates,
    duplicate_check

USAGE:
    python build_robokop_triples.py input.csv output.csv

NOTE ON NETWORK ACCESS:
    This script calls public NCATS Translator services and the public
    ROBOKOP Automat Cypher endpoint:
        https://name-resolution-sri.renci.org
        https://nodenormalization-sri.renci.org
        https://automat-u24.apps.renci.org/robokopkg/1.3/query
    These are free, public, unauthenticated APIs. Run this script from an
    environment with normal internet access (your local machine), not from
    a sandboxed environment that restricts outbound domains.
"""

import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

NAME_RES_URL = "https://name-resolution-sri.renci.org/lookup"
NODE_NORM_URL = "https://nodenormalization-sri.renci.org/get_normalized_nodes"
ROBOKOP_CYPHER_URL = "https://automat-u24.apps.renci.org/robokopkg/1.3/query"

# Local caches so re-running the pipeline across many datasets does not
# repeatedly hit the same APIs for the same metabolite/disease/edge.
CURIE_CACHE_PATH = Path("curie_cache.json")
EXISTENCE_CACHE_PATH = Path("existence_cache.json")

# Predicate assignment: keyed by (evidence_type, direction).
# Extend this table as new study designs come up -- do not hardcode new
# logic elsewhere in the script.
PREDICATE_MAP = {
    ("group_comparison", "increased"): "biolink:associated_with_increased_likelihood_of",
    ("group_comparison", "decreased"): "biolink:associated_with_decreased_likelihood_of",
    ("numeric_correlation", "increased"): "biolink:positively_correlated_with",
    ("numeric_correlation", "decreased"): "biolink:negatively_correlated_with",
}
DEFAULT_PREDICATE = "biolink:associated_with"

KNOWLEDGE_SOURCE = "infores:aria-t2d-metabolomics-pipeline"


# --------------------------------------------------------------------------
# Caching helpers
# --------------------------------------------------------------------------

def load_json_cache(path):
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_json_cache(path, cache):
    path.write_text(json.dumps(cache, indent=2))


def http_get_json(url, params, timeout=15):
    query = urllib.parse.urlencode(params)
    req = urllib.request.Request(f"{url}?{query}", headers={"accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def http_post_json(url, payload, timeout=15):
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, headers={"Content-Type": "application/json", "Accept": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


# --------------------------------------------------------------------------
# Step 1: entity normalization
# --------------------------------------------------------------------------

def resolve_name_to_curie(name, cache):
    """
    Name Resolution API (free text -> candidate CURIE), then NodeNormalization
    API (candidate CURIE -> canonical preferred CURIE). Returns
    (curie, preferred_name, category, status); status is
    "resolved" / "unresolved" / "error".
    """
    if name in cache:
        entry = cache[name]
        return entry["curie"], entry["preferred_name"], entry["category"], entry["status"]

    curie, preferred_name, category, status = None, name, None, "unresolved"

    try:
        lookup_result = http_get_json(
            NAME_RES_URL, {"string": name, "autocomplete": "false", "limit": 5}
        )
        if isinstance(lookup_result, dict) and lookup_result:
            curie = next(iter(lookup_result.keys()))
        elif isinstance(lookup_result, list) and lookup_result:
            curie = lookup_result[0].get("curie")

        if curie:
            norm_result = http_get_json(NODE_NORM_URL, {"curie": curie})
            node_info = norm_result.get(curie)
            if node_info:
                curie = node_info["id"]["identifier"]
                preferred_name = node_info["id"].get("label", name)
                categories = node_info.get("type", [])
                category = categories[0] if categories else None
                status = "resolved"

    except Exception as exc:  # noqa: BLE001 - log and continue, don't crash the batch
        print(f"  [warn] normalization failed for '{name}': {exc}", file=sys.stderr)
        status = "error"

    cache[name] = {
        "curie": curie,
        "preferred_name": preferred_name,
        "category": category,
        "status": status,
    }
    return curie, preferred_name, category, status


# --------------------------------------------------------------------------
# Step 2: predicate assignment
# --------------------------------------------------------------------------

def assign_predicate(evidence_type, direction):
    key = ((evidence_type or "").strip().lower(), (direction or "").strip().lower())
    predicate = PREDICATE_MAP.get(key)
    if predicate is None:
        print(
            f"  [warn] no predicate rule for evidence_type='{evidence_type}', "
            f"direction='{direction}' -- falling back to '{DEFAULT_PREDICATE}'. "
            f"Add a rule to PREDICATE_MAP if this evidence type is expected.",
            file=sys.stderr,
        )
        predicate = DEFAULT_PREDICATE
    return predicate


# --------------------------------------------------------------------------
# Step 3: existence check against the live ROBOKOP KG
# --------------------------------------------------------------------------

def check_existence(subject_curie, object_curie, cache):
    """
    Queries the ROBOKOP Automat Cypher endpoint to check:
      - does the subject node exist at all?
      - does any edge already connect subject and object, and with what
        predicate(s)?
    Returns (node_exists, existing_predicates) where node_exists is
    True/False/None (None = check failed / could not verify) and
    existing_predicates is a list of predicate strings (empty if no edge).
    """
    cache_key = f"{subject_curie}|{object_curie}"
    if cache_key in cache:
        entry = cache[cache_key]
        return entry["node_exists"], entry["existing_predicates"]

    node_exists, existing_predicates = None, []

    if not subject_curie:
        cache[cache_key] = {"node_exists": None, "existing_predicates": []}
        return None, []

    try:
        node_query = {
            "query": f'MATCH (n {{id: "{subject_curie}"}}) RETURN n LIMIT 1'
        }
        node_result = http_post_json(ROBOKOP_CYPHER_URL, node_query)
        rows = node_result.get("results", [{}])[0].get("data", [])
        node_exists = len(rows) > 0

        if node_exists and object_curie:
            edge_query = {
                "query": (
                    f'MATCH (n1 {{id: "{subject_curie}"}})-[r]-(n2 {{id: "{object_curie}"}}) '
                    f"RETURN type(r) AS predicate"
                )
            }
            edge_result = http_post_json(ROBOKOP_CYPHER_URL, edge_query)
            edge_rows = edge_result.get("results", [{}])[0].get("data", [])
            existing_predicates = list(
                {row["row"][0] for row in edge_rows if row.get("row")}
            )

    except Exception as exc:  # noqa: BLE001
        print(f"  [warn] existence check failed for {subject_curie}: {exc}", file=sys.stderr)
        node_exists, existing_predicates = None, []

    cache[cache_key] = {"node_exists": node_exists, "existing_predicates": existing_predicates}
    return node_exists, existing_predicates


def classify_duplicate_risk(node_exists, existing_predicates, assigned_predicate):
    if node_exists is None:
        return "unverified (existence check failed)"
    if not node_exists:
        return "new_node"
    if not existing_predicates:
        return "new_edge"
    if assigned_predicate in existing_predicates:
        return "possible_duplicate (same predicate already present)"
    return "new_edge_type (node/edge exist, different predicate)"


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def build_triples(input_path, output_path):
    curie_cache = load_json_cache(CURIE_CACHE_PATH)
    existence_cache = load_json_cache(EXISTENCE_CACHE_PATH)
    rows_out = []
    unresolved = []

    with open(input_path, newline="", encoding="utf-8") as f_in:
        reader = csv.DictReader(f_in)
        for i, row in enumerate(reader, start=1):
            metabolite_name = row["metabolite_name"].strip()
            disease_name = row["disease_name"].strip()
            disease_curie_input = row.get("disease_curie", "").strip()

            print(f"[{i}] Resolving: {metabolite_name}")
            subj_curie, subj_name, subj_cat, subj_status = resolve_name_to_curie(
                metabolite_name, curie_cache
            )
            time.sleep(0.2)  # be polite to the public API

            if disease_curie_input:
                obj_curie, obj_name, obj_cat, obj_status = (
                    disease_curie_input, disease_name, "biolink:Disease", "provided",
                )
            else:
                obj_curie, obj_name, obj_cat, obj_status = resolve_name_to_curie(
                    disease_name, curie_cache
                )
                time.sleep(0.2)

            predicate = assign_predicate(row.get("evidence_type"), row.get("direction"))

            print("    Checking existence in ROBOKOP KG...")
            node_exists, existing_predicates = check_existence(
                subj_curie, obj_curie, existence_cache
            )
            time.sleep(0.2)
            duplicate_check = classify_duplicate_risk(node_exists, existing_predicates, predicate)

            triple_row = {
                "subject_id": subj_curie or "",
                "subject_name": subj_name,
                "subject_category": subj_cat or "",
                "predicate": predicate,
                "object_id": obj_curie or "",
                "object_name": obj_name,
                "object_category": obj_cat or "",
                "direction": row.get("direction", ""),
                "p_value": row.get("p_value", ""),
                "evidence": (
                    f"{row.get('evidence_type', 'unspecified')} study; "
                    f"test={row.get('test_type', 'unspecified')}; "
                    f"dataset={row.get('dataset_id', 'unspecified')}"
                ),
                "knowledge_source": KNOWLEDGE_SOURCE,
                "dataset_id": row.get("dataset_id", ""),
                "normalization_status": f"subject={subj_status}, object={obj_status}",
                "subject_node_exists": node_exists,
                "existing_edge_predicates": "; ".join(existing_predicates) if existing_predicates else "",
                "duplicate_check": duplicate_check,
            }
            rows_out.append(triple_row)

            if subj_status != "resolved":
                unresolved.append(metabolite_name)

    fieldnames = [
        "subject_id", "subject_name", "subject_category", "predicate",
        "object_id", "object_name", "object_category",
        "direction", "p_value", "evidence", "knowledge_source", "dataset_id",
        "normalization_status", "subject_node_exists", "existing_edge_predicates",
        "duplicate_check",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    save_json_cache(CURIE_CACHE_PATH, curie_cache)
    save_json_cache(EXISTENCE_CACHE_PATH, existence_cache)

    print(f"\nDone. Wrote {len(rows_out)} triples to {output_path}")
    if unresolved:
        print(
            f"\n[review needed] {len(unresolved)} metabolite name(s) did not "
            f"resolve automatically and need manual CURIE lookup:"
        )
        for name in unresolved:
            print(f"  - {name}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python build_robokop_triples.py input.csv output.csv")
        sys.exit(1)
    build_triples(sys.argv[1], sys.argv[2])
