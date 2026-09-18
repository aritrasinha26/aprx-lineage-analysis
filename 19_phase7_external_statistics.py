import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, chi2_contingency


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"
    outdir = root/"preprint"/"results"/"auto"

    df = pd.read_csv(
        ext/"external_integrated_serralysin_table.tsv",
        sep="\t", dtype=str
    ).fillna("")
    df = df[df["status"]=="candidate_serralysin"].copy()

    rows = []

    aprx = df["nearest_reference_class"]=="AprX-reference-proximal"
    apra = df["nearest_reference_class"]=="P.aeruginosa-AprA-reference-proximal"
    sv = df["SVMSY_aligned_pattern"]=="SVMSY"

    table = np.array([
        [(apra & sv).sum(), (apra & ~sv).sum()],
        [(aprx & sv).sum(), (aprx & ~sv).sum()]
    ])
    odds, p = fisher_exact(table)
    rows.append({
        "analysis":"external_SVMSY_enrichment_AprAprox_vs_AprXprox",
        "n":int(table.sum()),
        "effect_name":"odds_ratio",
        "effect":odds,
        "p_value":p,
        "notes":f"table_AprAprox_vs_AprXprox={table.tolist()}"
    })

    # Four-state motif distribution vs reference-proximity class
    tab = pd.crosstab(
        df["nearest_reference_class"],
        df["SVMSY_aligned_pattern"]
    )
    chi2, cp, dof, expected = chi2_contingency(tab)
    rows.append({
        "analysis":"external_marker_pattern_vs_reference_proximity_class",
        "n":len(df),
        "effect_name":"chi_square",
        "effect":chi2,
        "p_value":cp,
        "notes":f"dof={dof};patterns={','.join(tab.columns)}"
    })

    # Context enrichment between the two reference-proximity classes.
    for col in [
        "inhibitor_nearby",
        "T1SS_ABC_or_AprD_nearby",
        "HlyD_or_AprE_nearby",
        "TolC_or_AprF_nearby",
        "lipase_gene_nearby"
    ]:
        yes = df[col]=="YES"
        t = np.array([
            [(apra & yes).sum(), (apra & ~yes).sum()],
            [(aprx & yes).sum(), (aprx & ~yes).sum()]
        ])
        if t.sum() == 0:
            continue
        o, pv = fisher_exact(t)
        rows.append({
            "analysis":f"external_{col}_AprAprox_vs_AprXprox",
            "n":int(t.sum()),
            "effect_name":"odds_ratio",
            "effect":o,
            "p_value":pv,
            "notes":f"table={t.tolist()}"
        })

    res = pd.DataFrame(rows)
    res.to_csv(
        outdir/"phase7_external_statistics.tsv",
        sep="\t", index=False
    )

    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
