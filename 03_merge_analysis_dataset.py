import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


PREFIXES = ["WS","DSM","LMG","ATCC","ICMP","CFBP","CCUG","JCM","NBRC","NCTC","CECT","CIP","KCTC","BCRC","MT"]


def strain_key(s):
    t = str(s).upper().replace("$","5").replace("OSM","DSM")
    for prefix in PREFIXES:
        m = re.search(rf"\b{re.escape(prefix)}\s*[-:]?\s*(\d+)\b", t)
        if m:
            return f"{prefix}{m.group(1)}"
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()
    root = Path(args.project_root)

    master = pd.read_csv(root/"preprint"/"results"/"AprX_master_table.tsv", sep="\t", dtype=str).fillna("")
    features = pd.read_csv(root/"preprint"/"results"/"auto"/"AprX_sequence_features.tsv", sep="\t", dtype=str).fillna("")
    phen = pd.read_csv(root/"preprint"/"results"/"auto"/"Figure2_phenotypes_auto.tsv", sep="\t", dtype=str).fillna("")

    master["strain_key"] = master["strain"].map(strain_key)
    phen["strain_key_norm"] = phen["strain"].map(strain_key)

    # Sequence features are keyed safely by assembly.
    fcols = [c for c in features.columns if c not in {"protein_id"}]
    merged = master.merge(features[fcols], on="assembly", how="left", suffixes=("","_seq"))

    pcols = [
        "strain_key_norm","strain","match_score","match_method",
        "day3_category","day3_score","day4_category","day4_score",
        "day7_category","day7_score"
    ]
    p = phen[[c for c in pcols if c in phen.columns]].copy()
    p = p.rename(columns={"strain":"published_strain"})

    # Never merge on an empty strain identifier.
    p = p[p["strain_key_norm"].astype(str).str.strip() != ""].copy()

    # If OCR produced more than one match for the same strain,
    # keep the highest-confidence match only.
    p["_match_score_num"] = pd.to_numeric(
        p["match_score"], errors="coerce"
    ).fillna(-1)

    p = (
        p.sort_values("_match_score_num", ascending=False)
         .drop_duplicates("strain_key_norm", keep="first")
         .drop(columns="_match_score_num")
    )

    expected_rows = len(merged)

    merged = merged.merge(
        p,
        left_on="strain_key",
        right_on="strain_key_norm",
        how="left"
    )

    if len(merged) != expected_rows:
        raise RuntimeError(
            f"Merge changed genome count from {expected_rows} to {len(merged)}. "
            "Check duplicated strain identifiers before analysis."
        )

    out = root/"preprint"/"results"/"auto"/"AprX_analysis_dataset.tsv"
    merged.to_csv(out, sep="\t", index=False)

    matched = merged["day7_category"].notna() & (merged["day7_category"] != "")
    print("Rows:", len(merged))
    print("Rows with matched Figure 2 phenotype:", int(matched.sum()))
    print("Phenotype match rate:", f"{matched.mean()*100:.1f}%")
    print("Saved:", out)


if __name__ == "__main__":
    main()
