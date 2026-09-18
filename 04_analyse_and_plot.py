import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr, kruskal
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


def fisher_or(a,b,c,d):
    table = np.array([[a,b],[c,d]], dtype=int)
    odds, p = fisher_exact(table)
    return odds, p


def species_stratified_permutation(df, marker_col, score_col, n_iter, seed):
    x = df[[marker_col, score_col, "species"]].dropna().copy()
    x = x[x[marker_col].isin(["YES","NO"])]
    if x.empty or x[marker_col].nunique() < 2:
        return np.nan, np.nan, 0, 0

    x["_m"] = (x[marker_col] == "YES").astype(int)
    x["_y"] = pd.to_numeric(x[score_col], errors="coerce")
    x = x.dropna(subset=["_y"])

    def effect(frame):
        yes = frame.loc[frame["_m"]==1, "_y"]
        no = frame.loc[frame["_m"]==0, "_y"]
        if len(yes)==0 or len(no)==0:
            return np.nan
        return yes.mean() - no.mean()

    obs = effect(x)
    variable_species = [
        sp for sp,g in x.groupby("species") if g["_m"].nunique() > 1
    ]
    informative_n = int(x[x["species"].isin(variable_species)].shape[0])
    if not variable_species:
        return obs, np.nan, 0, informative_n

    rng = np.random.default_rng(seed)
    vals = []
    base = x.copy()
    for _ in range(n_iter):
        perm = base.copy()
        for sp, idx in perm.groupby("species").groups.items():
            idx = list(idx)
            if len(set(perm.loc[idx,"_m"])) > 1:
                perm.loc[idx,"_m"] = rng.permutation(perm.loc[idx,"_m"].values)
        vals.append(effect(perm))
    vals = np.asarray(vals, dtype=float)
    p = (1 + np.sum(np.abs(vals) >= abs(obs))) / (1 + len(vals))
    return obs, p, len(variable_species), informative_n


def site_association(root, data, min_count):
    aln_path = root/"preprint"/"results"/"serralysin_reference_alignment.faa"
    aln = {r.id:str(r.seq).upper() for r in SeqIO.parse(aln_path,"fasta")}
    candidates = {k:v for k,v in aln.items() if k.startswith("GCF_")}

    # phenotype keyed by assembly
    dd = data.copy()
    dd["day7_score_num"] = to_num(dd["day7_score"])
    pheno = dict(zip(dd["assembly"], dd["day7_score_num"]))

    cy_id = next((k for k in aln if "REF_AprX_CY091" in k), None)
    cy_seq = aln.get(cy_id, "")
    cypos = []
    p=0
    for aa in cy_seq:
        if aa != "-":
            p += 1
            cypos.append(p)
        else:
            cypos.append("")

    rows=[]
    L = len(next(iter(candidates.values())))
    for col in range(L):
        groups=defaultdict(list)
        for sid, seq in candidates.items():
            assembly=sid.split("|")[0]
            y=pheno.get(assembly, np.nan)
            aa=seq[col]
            if aa in "-X?" or pd.isna(y):
                continue
            groups[aa].append(float(y))
        groups={aa:v for aa,v in groups.items() if len(v)>=min_count}
        if len(groups)<2:
            continue
        try:
            stat,pv=kruskal(*groups.values())
        except ValueError:
            continue
        rows.append({
            "alignment_column_1based":col+1,
            "CY091_position":cypos[col] if cypos else "",
            "n_groups_tested":len(groups),
            "groups":";".join(f"{aa}:{len(v)}" for aa,v in sorted(groups.items())),
            "kruskal_H":stat,
            "p_value":pv
        })
    out=pd.DataFrame(rows)
    if not out.empty:
        out["FDR_BH"]=multipletests(out["p_value"], method="fdr_bh")[1]
        out=out.sort_values(["FDR_BH","p_value"])
    else:
        out = pd.DataFrame(columns=[
            "alignment_column_1based","CY091_position","n_groups_tested",
            "groups","kruskal_H","p_value","FDR_BH"
        ])
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args=ap.parse_args()
    root=Path(args.project_root)
    cfg=json.loads((root/"preprint"/"automation"/"config.json").read_text())

    resdir=root/"preprint"/"results"/"auto"
    figdir=root/"preprint"/"figures"/"auto"
    figdir.mkdir(parents=True, exist_ok=True)

    df=pd.read_csv(resdir/"AprX_analysis_dataset.tsv", sep="\t", dtype=str).fillna("")
    for c in ["day3_score","day4_score","day7_score","AprX_distance","percent_identity_to_CY091"]:
        if c in df:
            df[c+"_num"]=to_num(df[c])

    results=[]

    # 1. Raw milk association with AprX-reference-like classification.
    raw=df["isolation_source"].str.lower().eq("raw milk")
    aprx=df["classification"].eq("AprX-reference-like")
    a=int((raw & aprx).sum()); b=int((raw & ~aprx).sum())
    c=int((~raw & aprx).sum()); d=int((~raw & ~aprx).sum())
    odds,p=fisher_or(a,b,c,d)
    results.append({
        "analysis":"raw_milk_vs_AprX_reference_like",
        "n":len(df),
        "effect":odds,
        "effect_name":"odds_ratio",
        "p_value":p,
        "notes":f"table=[[{a},{b}],[{c},{d}]]"
    })

    # 2. SVMSY vs phenotype.
    sub=df[(df["classification"]=="AprX-reference-like") & df["day7_score_num"].notna()].copy()
    if "SVMSY_exact" in sub and sub["SVMSY_exact"].isin(["YES","NO"]).any():
        yes=sub.loc[sub["SVMSY_exact"]=="YES","day7_score_num"].dropna()
        no=sub.loc[sub["SVMSY_exact"]=="NO","day7_score_num"].dropna()
        if len(yes)>0 and len(no)>0:
            u,pv=mannwhitneyu(yes,no,alternative="two-sided")
            results.append({
                "analysis":"SVMSY_vs_day7_ordinal_unadjusted",
                "n":len(yes)+len(no),
                "effect":yes.mean()-no.mean(),
                "effect_name":"mean_score_difference_YES_minus_NO",
                "p_value":pv,
                "notes":f"n_yes={len(yes)};n_no={len(no)};Mann-Whitney_U={u}"
            })

            obs,pp,nsp,ninf=species_stratified_permutation(
                sub,"SVMSY_exact","day7_score_num",
                cfg["permutation_iterations"],cfg["random_seed"]
            )
            results.append({
                "analysis":"SVMSY_vs_day7_species_stratified_permutation",
                "n":len(yes)+len(no),
                "effect":obs,
                "effect_name":"mean_score_difference_YES_minus_NO",
                "p_value":pp,
                "notes":f"species_with_within_species_marker_variation={nsp};informative_rows={ninf}"
            })

            # Strong/very strong vs <= moderate
            sub2=sub[sub["SVMSY_exact"].isin(["YES","NO"])].copy()
            sub2["high"]=sub2["day7_score_num"]>=3
            aa=int(((sub2["SVMSY_exact"]=="YES") & sub2["high"]).sum())
            bb=int(((sub2["SVMSY_exact"]=="YES") & ~sub2["high"]).sum())
            cc=int(((sub2["SVMSY_exact"]=="NO") & sub2["high"]).sum())
            dd=int(((sub2["SVMSY_exact"]=="NO") & ~sub2["high"]).sum())
            if min(aa+bb,cc+dd)>0:
                oo,fp=fisher_or(aa,bb,cc,dd)
                results.append({
                    "analysis":"SVMSY_vs_high_day7_activity",
                    "n":len(sub2),
                    "effect":oo,
                    "effect_name":"odds_ratio",
                    "p_value":fp,
                    "notes":f"high=strong_or_very_strong;table=[[{aa},{bb}],[{cc},{dd}]]"
                })

    # 3. AprX phylogenetic distance vs phenotype.
    if "AprX_distance_num" in sub:
        x=sub[["AprX_distance_num","day7_score_num"]].dropna()
        if len(x)>=5 and x["AprX_distance_num"].nunique()>1:
            rho,pv=spearmanr(x["AprX_distance_num"],x["day7_score_num"])
            results.append({
                "analysis":"nearest_AprX_reference_distance_vs_day7",
                "n":len(x),
                "effect":rho,
                "effect_name":"Spearman_rho",
                "p_value":pv,
                "notes":"descriptive association; phylogenetic non-independence remains"
            })

    stats=pd.DataFrame(results)
    stats.to_csv(resdir/"statistical_results.tsv", sep="\t", index=False)

    sites=site_association(root,df,cfg["site_test_min_residue_count"])
    sites.to_csv(resdir/"AprX_site_association_day7.tsv", sep="\t", index=False)

    # Figure 1: species distribution by classification.
    counts=(df.groupby(["species","classification"]).size()
              .unstack(fill_value=0))
    totals=counts.sum(axis=1)
    counts=counts.loc[totals.sort_values(ascending=True).index]
    fig,ax=plt.subplots(figsize=(9,max(5,0.28*len(counts))))
    counts.plot(kind="barh",stacked=True,ax=ax)
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("")
    ax.set_title("AprX/serralysin classification across the 87-genome cohort")
    fig.tight_layout()
    fig.savefig(figdir/"Fig_species_by_classification.png",dpi=300)
    plt.close(fig)

    # Figure 2: source vs classification, top sources.
    src=(df.assign(source=df["isolation_source"].replace({"":"missing"}))
           .groupby(["source","classification"]).size().unstack(fill_value=0))
    src=src.loc[src.sum(axis=1).sort_values(ascending=True).index]
    fig,ax=plt.subplots(figsize=(9,max(4,0.32*len(src))))
    src.plot(kind="barh",stacked=True,ax=ax)
    ax.set_xlabel("Number of genomes")
    ax.set_ylabel("")
    ax.set_title("Isolation source and serralysin classification")
    fig.tight_layout()
    fig.savefig(figdir/"Fig_source_by_classification.png",dpi=300)
    plt.close(fig)

    # Figure 3: SVMSY marker vs day-7 phenotype.
    if "SVMSY_exact" in sub and sub["SVMSY_exact"].isin(["YES","NO"]).sum()>0:
        pdat=sub[sub["SVMSY_exact"].isin(["YES","NO"])].copy()
        groups=[pdat.loc[pdat["SVMSY_exact"]==g,"day7_score_num"].dropna().values for g in ["NO","YES"]]
        if all(len(g)>0 for g in groups):
            fig,ax=plt.subplots(figsize=(5.5,4.5))
            ax.boxplot(groups,tick_labels=["SVMSY absent","SVMSY present"])
            rng=np.random.default_rng(cfg["random_seed"])
            for i,g in enumerate(groups,1):
                jitter=rng.normal(i,0.035,len(g))
                ax.scatter(jitter,g,s=18,alpha=0.7)
            ax.set_ylabel("Day-7 proteolytic activity score (0–4)")
            ax.set_title("Published Maier phenotype versus SVMSY marker")
            ax.set_yticks([0,1,2,3,4],["non","weak","moderate","strong","very strong"])
            fig.tight_layout()
            fig.savefig(figdir/"Fig_SVMSY_vs_day7_phenotype.png",dpi=300)
            plt.close(fig)

    # Figure 4: nearest AprX-reference distance vs phenotype.
    if "AprX_distance_num" in sub:
        x=sub[["AprX_distance_num","day7_score_num"]].dropna()
        if len(x)>=5:
            fig,ax=plt.subplots(figsize=(5.5,4.5))
            ax.scatter(x["AprX_distance_num"],x["day7_score_num"],alpha=0.75)
            ax.set_xlabel("Phylogenetic distance to nearest AprX reference")
            ax.set_ylabel("Day-7 proteolytic activity score (0–4)")
            ax.set_yticks([0,1,2,3,4])
            fig.tight_layout()
            fig.savefig(figdir/"Fig_AprX_distance_vs_day7.png",dpi=300)
            plt.close(fig)

    print("Statistical tests:")
    print(stats.to_string(index=False) if not stats.empty else "No tests could be run.")
    print(f"Exploratory variable-site tests: {len(sites)}")
    if not sites.empty:
        print("Lowest FDR site tests:")
        print(sites.head(10).to_string(index=False))


if __name__=="__main__":
    main()
