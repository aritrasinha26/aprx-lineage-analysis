import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import mannwhitneyu
from scipy.optimize import minimize_scalar
import statsmodels.api as sm


def brownian_covariance(tree, names):
    tree.root_at_midpoint()
    tips = {t.name:t for t in tree.get_terminals()}
    n=len(names)
    C=np.zeros((n,n),float)
    for i,a in enumerate(names):
        ta=tips[a]
        for j in range(i+1):
            b=names[j]; tb=tips[b]
            if a==b:
                v=tree.distance(tree.root,ta)
            else:
                m=tree.common_ancestor(ta,tb)
                v=tree.distance(tree.root,m)
            C[i,j]=C[j,i]=v
    d=np.diag(C)
    s=np.sqrt(d)
    R=C/np.outer(s,s)
    np.fill_diagonal(R,1.0)
    return R


def lambda_sigma(R, lam):
    S=lam*R.copy()
    np.fill_diagonal(S,1.0)
    S+=np.eye(S.shape[0])*1e-6
    return S


def fit_lambda_pgls(treefile, d):
    tree=Phylo.read(treefile,"newick")
    R=brownian_covariance(tree,d["assembly"].tolist())
    X=sm.add_constant(d[["type1_vs_type8"]].astype(float),has_constant="add")
    y=d["day7_score_num"].astype(float)

    def obj(lam):
        try:
            r=sm.GLS(y,X,sigma=lambda_sigma(R,float(lam))).fit()
            return -r.llf if np.isfinite(r.llf) else np.inf
        except Exception:
            return np.inf

    opt=minimize_scalar(obj,bounds=(0,0.999),method="bounded")
    lam=float(opt.x)
    fit=sm.GLS(y,X,sigma=lambda_sigma(R,lam)).fit()
    return lam, fit


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()

    root=Path(args.project_root)
    res=root/"preprint"/"results"/"auto"
    opdir=res/"operon_homology"

    df=pd.read_csv(
        opdir/"pilot_reconstructed_operon_classes.tsv",
        sep="\t",dtype=str
    ).fillna("")
    df["day7_score_num"]=pd.to_numeric(df["day7_score"],errors="coerce")

    comp=df[
        df["reconstructed_operon_class"].isin([
            "Type1-like_local_prtAB",
            "Type8-like_distal_prtAB"
        ])
        & df["day7_score_num"].notna()
    ].copy()

    g1=comp.loc[
        comp["reconstructed_operon_class"]=="Type1-like_local_prtAB",
        "day7_score_num"
    ]
    g8=comp.loc[
        comp["reconstructed_operon_class"]=="Type8-like_distal_prtAB",
        "day7_score_num"
    ]

    rows=[]

    if len(g1)>0 and len(g8)>0:
        u,p=mannwhitneyu(g1,g8,alternative="two-sided")
        rows.append({
            "analysis":"reconstructed_type1_vs_type8_unadjusted",
            "n":len(g1)+len(g8),
            "effect_name":"mean_day7_difference_type1_minus_type8",
            "effect":g1.mean()-g8.mean(),
            "p_value":p,
            "notes":f"n_type1={len(g1)};n_type8={len(g8)};U={u}"
        })

        comp["type1_vs_type8"]=(
            comp["reconstructed_operon_class"]=="Type1-like_local_prtAB"
        ).astype(int)

        treefile=res/"core_phylogeny"/"UBCG_core_tree.treefile"
        lam,fit=fit_lambda_pgls(treefile,comp)

        rows.append({
            "analysis":"reconstructed_type1_vs_type8_lambdaPGLS",
            "n":len(comp),
            "effect_name":"beta_type1_vs_type8",
            "effect":fit.params["type1_vs_type8"],
            "p_value":fit.pvalues["type1_vs_type8"],
            "notes":f"lambda={lam};SE={fit.bse['type1_vs_type8']};AIC={fit.aic}"
        })

    pd.DataFrame(rows).to_csv(
        opdir/"corrected_operon_phenotype_tests.tsv",
        sep="\t",index=False
    )

    print("Corrected reconstructed-operon phenotype comparison:")
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nDay-7 summaries:")
    print(
        comp.groupby("reconstructed_operon_class")["day7_score_num"]
        .agg(["count","mean","median","min","max"])
        .to_string()
    )


if __name__=="__main__":
    main()
