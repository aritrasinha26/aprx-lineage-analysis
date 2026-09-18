import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.optimize import minimize_scalar
import statsmodels.api as sm


def brownian_covariance(tree, names):
    tree.root_at_midpoint()
    tips = {t.name: t for t in tree.get_terminals()}
    n = len(names)
    C = np.zeros((n, n), dtype=float)

    for i, a in enumerate(names):
        ta = tips[a]
        for j in range(i + 1):
            b = names[j]
            tb = tips[b]
            if a == b:
                v = tree.distance(tree.root, ta)
            else:
                mrca = tree.common_ancestor(ta, tb)
                v = tree.distance(tree.root, mrca)
            C[i, j] = C[j, i] = v
    return C


def to_correlation(C):
    d = np.diag(C).copy()
    if np.any(d <= 0):
        raise RuntimeError("Non-positive diagonal found in phylogenetic covariance.")
    s = np.sqrt(d)
    R = C / np.outer(s, s)
    np.fill_diagonal(R, 1.0)
    return R


def lambda_covariance(R, lam, nugget=1e-6):
    S = lam * R.copy()
    np.fill_diagonal(S, 1.0)
    S += np.eye(S.shape[0]) * nugget
    return S


def fit_lambda_pgls(treefile, frame, response, predictors, name):
    d = frame[["assembly", response] + predictors].copy()

    for c in [response] + predictors:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna().copy()

    if len(d) < 8:
        return [], {"analysis": name, "status": "too_few_rows", "n": len(d)}

    if any(d[c].nunique() < 2 for c in predictors):
        return [], {
            "analysis": name,
            "status": "predictor_has_no_variation",
            "n": len(d)
        }

    tree = Phylo.read(treefile, "newick")
    names = d["assembly"].tolist()

    C = brownian_covariance(tree, names)
    R = to_correlation(C)

    eig = np.linalg.eigvalsh(R)
    raw_condition = np.linalg.cond(R)

    X = sm.add_constant(d[predictors].astype(float), has_constant="add")
    y = d[response].astype(float)

    def objective(lam):
        try:
            S = lambda_covariance(R, float(lam))
            result = sm.GLS(y, X, sigma=S).fit()
            return np.inf if not np.isfinite(result.llf) else -result.llf
        except Exception:
            return np.inf

    opt = minimize_scalar(
        objective,
        bounds=(0.0, 0.999),
        method="bounded",
        options={"xatol": 1e-5}
    )

    lam = float(opt.x)
    S = lambda_covariance(R, lam)
    result = sm.GLS(y, X, sigma=S).fit()

    rows = []
    for term in result.params.index:
        rows.append({
            "analysis": name,
            "n": len(d),
            "lambda_ML": lam,
            "term": term,
            "coefficient": float(result.params[term]),
            "SE": float(result.bse[term]),
            "t": float(result.tvalues[term]),
            "p_value": float(result.pvalues[term]),
            "AIC": float(result.aic),
            "logLik": float(result.llf)
        })

    diag = {
        "analysis": name,
        "status": "ok",
        "n": len(d),
        "lambda_ML": lam,
        "raw_phylo_correlation_condition_number": raw_condition,
        "raw_min_eigenvalue": float(eig.min()),
        "raw_max_eigenvalue": float(eig.max()),
        "lambda_covariance_condition_number": float(np.linalg.cond(S)),
        "optimizer_success": bool(opt.success)
    }

    return rows, diag


def fitch_changes(treefile, state_map):
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    changes = 0

    def visit(clade):
        nonlocal changes
        if clade.is_terminal():
            state = state_map.get(clade.name)
            return set() if state is None else {state}

        child_sets = [visit(c) for c in clade.clades]
        child_sets = [x for x in child_sets if x]
        if not child_sets:
            return set()

        inter = set.intersection(*child_sets)
        if inter:
            return inter

        changes += 1
        return set.union(*child_sets)

    visit(tree.root)
    return changes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    treefile = res/"core_phylogeny"/"UBCG_core_tree.treefile"

    dataset = res/"AprX_analysis_dataset_with_operon.tsv"
    if not dataset.exists():
        dataset = res/"AprX_analysis_dataset.tsv"

    df = pd.read_csv(dataset, sep="\t", dtype=str).fillna("")

    df["day7_score_num"] = pd.to_numeric(df["day7_score"], errors="coerce")
    df["AprX_distance_num"] = pd.to_numeric(df["AprX_distance"], errors="coerce")
    df["SVMSY_binary"] = np.where(
        df["SVMSY_exact"] == "YES", 1,
        np.where(df["SVMSY_exact"] == "NO", 0, np.nan)
    )

    aprx = df[df["classification"] == "AprX-reference-like"].copy()

    if "operon_type" in aprx.columns:
        aprx["operon1_vs_8"] = np.nan
        aprx.loc[aprx["operon_type"] == "1", "operon1_vs_8"] = 1
        aprx.loc[aprx["operon_type"] == "8", "operon1_vs_8"] = 0

    models = [
        ("lambdaPGLS_day7_vs_SVMSY", ["SVMSY_binary"]),
        ("lambdaPGLS_day7_vs_AprX_distance", ["AprX_distance_num"])
    ]

    if "operon1_vs_8" in aprx.columns:
        models.append(
            ("lambdaPGLS_day7_operon1_vs_8", ["operon1_vs_8"])
        )

    result_rows = []
    diag_rows = []

    for name, predictors in models:
        rows, diag = fit_lambda_pgls(
            treefile, aprx, "day7_score_num", predictors, name
        )
        result_rows.extend(rows)
        diag_rows.append(diag)

    results = pd.DataFrame(result_rows)
    diagnostics = pd.DataFrame(diag_rows)

    results.to_csv(res/"lambda_pgls_results.tsv", sep="\t", index=False)
    diagnostics.to_csv(res/"lambda_pgls_diagnostics.tsv", sep="\t", index=False)

    transition_rows = []

    sv = aprx.dropna(subset=["SVMSY_binary"])
    sv_map = dict(zip(sv["assembly"], sv["SVMSY_binary"].astype(int)))
    transition_rows.append({
        "trait": "SVMSY_binary",
        "n_tips": len(sv_map),
        "minimum_Fitch_changes": fitch_changes(treefile, sv_map)
    })

    if "operon_type" in aprx.columns:
        op = aprx[aprx["operon_type"].isin(["1", "8"])].copy()
        op_map = dict(zip(op["assembly"], op["operon_type"]))
        transition_rows.append({
            "trait": "operon_type_1_vs_8",
            "n_tips": len(op_map),
            "minimum_Fitch_changes": fitch_changes(treefile, op_map)
        })

    ph = aprx.dropna(subset=["day7_score_num"]).copy()
    ph["high_activity"] = (ph["day7_score_num"] >= 3).astype(int)
    ph_map = dict(zip(ph["assembly"], ph["high_activity"]))
    transition_rows.append({
        "trait": "high_day7_activity",
        "n_tips": len(ph_map),
        "minimum_Fitch_changes": fitch_changes(treefile, ph_map)
    })

    transitions = pd.DataFrame(transition_rows)
    transitions.to_csv(
        res/"phylogenetic_trait_transition_counts.tsv",
        sep="\t",
        index=False
    )

    summary = []
    summary.append("# Phase 4 phylogenetic diagnostic summary")
    summary.append("")
    summary.append(
        "The previous fixed-Brownian PGLS produced extremely large standard errors. "
        "This phase diagnoses that covariance and re-fits PGLS with ML-estimated "
        "Pagel-like lambda on a normalized phylogenetic correlation matrix."
    )
    summary.append("")
    summary.append("## Covariance diagnostics")
    summary.append("")

    for _, r in diagnostics.iterrows():
        summary.append(
            f"- {r['analysis']}: n={int(r['n'])}; status={r['status']}; "
            f"lambda={r.get('lambda_ML', np.nan):.4g}; "
            f"raw condition number={r.get('raw_phylo_correlation_condition_number', np.nan):.4g}; "
            f"lambda-adjusted condition number={r.get('lambda_covariance_condition_number', np.nan):.4g}"
        )

    summary.append("")
    summary.append("## Lambda-PGLS")
    summary.append("")

    if results.empty:
        summary.append("- No model was successfully fitted.")
    else:
        for _, r in results.iterrows():
            if r["term"] == "const":
                continue
            summary.append(
                f"- {r['analysis']}: lambda={r['lambda_ML']:.4g}; "
                f"beta={r['coefficient']:.4g}; SE={r['SE']:.4g}; "
                f"p={r['p_value']:.4g}; n={int(r['n'])}"
            )

    summary.append("")
    summary.append("## Independent phylogenetic replication")
    summary.append("")

    for _, r in transitions.iterrows():
        summary.append(
            f"- {r['trait']}: {int(r['minimum_Fitch_changes'])} minimum changes "
            f"across {int(r['n_tips'])} observed tips."
        )

    summary.append("")
    summary.append(
        "Interpretation rule: a trait with only a few transitions has little independent "
        "evolutionary replication, even if a pooled strain-level association appears strong."
    )

    out = res/"phase4_phylo_diagnostic_summary.md"
    out.write_text("\n".join(summary) + "\n")

    print("\n".join(summary))
    print("\nSaved:")
    print(res/"lambda_pgls_results.tsv")
    print(res/"lambda_pgls_diagnostics.tsv")
    print(res/"phylogenetic_trait_transition_counts.tsv")
    print(out)


if __name__ == "__main__":
    main()
