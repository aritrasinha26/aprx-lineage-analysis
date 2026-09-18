import argparse
from pathlib import Path
import math

import numpy as np
import pandas as pd
from Bio import Phylo
import statsmodels.api as sm


def brownian_covariance(tree, names):
    tree.root_at_midpoint()
    tips = {t.name: t for t in tree.get_terminals()}
    n = len(names)
    C = np.zeros((n,n), dtype=float)

    for i,a in enumerate(names):
        ta = tips[a]
        for j,b in enumerate(names[:i+1]):
            tb = tips[b]
            if i == j:
                v = tree.distance(tree.root, ta)
            else:
                mrca = tree.common_ancestor(ta, tb)
                v = tree.distance(tree.root, mrca)
            C[i,j] = C[j,i] = v

    # Numerical stabilization.
    scale = np.nanmedian(np.diag(C))
    nugget = max(scale * 1e-6, 1e-8)
    C += np.eye(n) * nugget
    return C


def fit_pgls(tree, frame, ycol, predictors, analysis_name):
    x = frame[[ycol] + predictors + ["assembly"]].copy()
    for c in [ycol] + predictors:
        x[c] = pd.to_numeric(x[c], errors="coerce")
    x = x.dropna()
    if len(x) < max(8, len(predictors)+4):
        return None

    names = x["assembly"].tolist()
    C = brownian_covariance(tree, names)
    X = sm.add_constant(x[predictors].astype(float), has_constant="add")
    y = x[ycol].astype(float)
    model = sm.GLS(y, X, sigma=C)
    result = model.fit()

    rows = []
    for term in result.params.index:
        rows.append({
            "analysis": analysis_name,
            "n": len(x),
            "term": term,
            "coefficient": result.params[term],
            "SE": result.bse[term],
            "t": result.tvalues[term],
            "p_value": result.pvalues[term],
            "AIC": result.aic,
            "R2_pseudo_note": "GLS; day-7 ordinal score treated as approximately continuous"
        })
    return rows


def fitch_transitions(tree, state_map):
    tree.root_at_midpoint()
    transitions = 0

    def post(clade):
        nonlocal transitions
        if clade.is_terminal():
            s = state_map.get(clade.name)
            return set() if s is None else {s}

        child_sets = [post(c) for c in clade.clades]
        child_sets = [s for s in child_sets if s]
        if not child_sets:
            return set()

        inter = set.intersection(*child_sets)
        if inter:
            return inter

        transitions += 1
        return set.union(*child_sets)

    post(tree.root)
    return transitions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    work = res/"core_phylogeny"

    treefile = work/"UBCG_core_tree.treefile"
    tree = Phylo.read(treefile, "newick")

    dataset = res/"AprX_analysis_dataset_with_operon.tsv"
    if not dataset.exists():
        dataset = res/"AprX_analysis_dataset.tsv"

    df = pd.read_csv(dataset, sep="\t", dtype=str).fillna("")
    tips = {t.name for t in tree.get_terminals()}
    df = df[df["assembly"].isin(tips)].copy()

    df["day7_score_num"] = pd.to_numeric(df["day7_score"], errors="coerce")
    df["AprX_distance_num"] = pd.to_numeric(df["AprX_distance"], errors="coerce")
    df["SVMSY_binary"] = np.where(df["SVMSY_exact"]=="YES",1,
                           np.where(df["SVMSY_exact"]=="NO",0,np.nan))

    aprx = df[df["classification"]=="AprX-reference-like"].copy()

    out_rows = []

    # Motif model.
    r = fit_pgls(
        Phylo.read(treefile, "newick"),
        aprx, "day7_score_num", ["SVMSY_binary"],
        "PGLS_day7_vs_SVMSY"
    )
    if r:
        out_rows.extend(r)

    # Distance model.
    r = fit_pgls(
        Phylo.read(treefile, "newick"),
        aprx, "day7_score_num", ["AprX_distance_num"],
        "PGLS_day7_vs_AprX_reference_distance"
    )
    if r:
        out_rows.extend(r)

    # Main operon contrast, using Type 1 vs Type 8 because both are common in this cohort.
    if "operon_type" in aprx.columns:
        aprx["operon1_vs_8"] = np.nan
        aprx.loc[aprx["operon_type"]=="1","operon1_vs_8"] = 1
        aprx.loc[aprx["operon_type"]=="8","operon1_vs_8"] = 0

        r = fit_pgls(
            Phylo.read(treefile, "newick"),
            aprx, "day7_score_num", ["operon1_vs_8"],
            "PGLS_day7_operon_type1_vs_type8"
        )
        if r:
            out_rows.extend(r)

    results = pd.DataFrame(out_rows)
    results.to_csv(res/"phylogenetically_controlled_tests.tsv", sep="\t", index=False)

    # How many evolutionary changes in the proposed marker are actually represented?
    marker_map = {
        r["assembly"]: int(r["SVMSY_binary"])
        for _,r in aprx.dropna(subset=["SVMSY_binary"]).iterrows()
    }
    marker_tree = Phylo.read(treefile, "newick")
    min_changes = fitch_transitions(marker_tree, marker_map)

    # Export tree annotation table for plotting/iTOL.
    cols = [
        "assembly","species","strain","isolation_source","classification",
        "SVMSY_region_aligned","SVMSY_exact","day7_category","day7_score"
    ]
    if "operon_type" in df.columns:
        cols.append("operon_type")
    cols = [c for c in cols if c in df.columns]
    df[cols].to_csv(res/"core_tree_metadata.tsv", sep="\t", index=False)

    summary = [
        "# Phase 3 phylogenetic-control summary",
        "",
        f"- Core tree tips represented in analysis dataset: {len(df)}",
        f"- AprX-reference-like genomes in core tree: {len(aprx)}",
        f"- Minimum Fitch changes required for the binary SVMSY state on the core tree: {min_changes}",
        "",
        "## PGLS results",
        ""
    ]

    if results.empty:
        summary.append("No PGLS models could be fitted.")
    else:
        for _,r in results.iterrows():
            summary.append(
                f"- {r['analysis']} / {r['term']}: "
                f"beta={r['coefficient']:.4g}, SE={r['SE']:.4g}, "
                f"p={r['p_value']:.4g}, n={int(r['n'])}"
            )

    summary += [
        "",
        "Interpretation note: day-7 phenotype is an ordinal 0–4 score. "
        "PGLS treats it as approximately continuous and should be used as a phylogenetic sensitivity analysis, "
        "not as the sole inferential model."
    ]

    (res/"phase3_phylogeny_summary.md").write_text("\n".join(summary)+"\n")

    print("\n".join(summary))
    print("\nSaved:")
    print(res/"phylogenetically_controlled_tests.tsv")
    print(res/"core_tree_metadata.tsv")
    print(res/"phase3_phylogeny_summary.md")


if __name__ == "__main__":
    main()
