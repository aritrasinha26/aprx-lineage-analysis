import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from Bio import SeqIO


def num(s):
    return pd.to_numeric(s, errors="coerce")


def eta_squared(groups, y):
    valid = np.array([g not in {"-", "X", "?"} for g in groups]) & np.isfinite(y)
    groups = np.asarray(groups)[valid]
    y = np.asarray(y, dtype=float)[valid]
    if len(y) < 5 or len(set(groups)) < 2:
        return np.nan
    grand = y.mean()
    ss_total = np.sum((y-grand)**2)
    if ss_total <= 0:
        return 0.0
    ss_between = 0.0
    for aa in set(groups):
        vals = y[groups == aa]
        ss_between += len(vals) * (vals.mean()-grand)**2
    return ss_between/ss_total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--permutations", type=int, default=2000)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    figs = root/"preprint"/"figures"/"auto"
    figs.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(res/"AprX_analysis_dataset.tsv", sep="\t", dtype=str).fillna("")
    df["day7_score_num"] = num(df["day7_score"])
    aprx = df[(df["classification"]=="AprX-reference-like") & df["day7_score_num"].notna()].copy()

    motif_col = "SVMSY_region_aligned"
    if motif_col not in aprx.columns:
        raise RuntimeError(f"{motif_col} not found in analysis dataset.")

    # Species x motif x phenotype summary
    species_rows = []
    for sp, g in aprx.groupby("species"):
        motifs = sorted(x for x in g[motif_col].unique() if x)
        species_rows.append({
            "species": sp,
            "n": len(g),
            "motif_patterns": ";".join(motifs),
            "n_motif_patterns": len(motifs),
            "day7_mean": g["day7_score_num"].mean(),
            "day7_median": g["day7_score_num"].median(),
            "day7_min": g["day7_score_num"].min(),
            "day7_max": g["day7_score_num"].max()
        })
    species_summary = pd.DataFrame(species_rows).sort_values(["n","species"], ascending=[False,True])
    species_summary.to_csv(res/"species_motif_phenotype_summary.tsv", sep="\t", index=False)

    # Motif variant phenotype summary
    motif_summary = (aprx.groupby(motif_col)
        .agg(n=("assembly","size"),
             day7_mean=("day7_score_num","mean"),
             day7_median=("day7_score_num","median"),
             day7_min=("day7_score_num","min"),
             day7_max=("day7_score_num","max"),
             n_species=("species","nunique"))
        .reset_index()
        .sort_values("n", ascending=False))
    motif_summary.to_csv(res/"motif_variant_phenotype_summary.tsv", sep="\t", index=False)

    tests = []
    # Primary direct comparison if both common variants exist
    if {"SVMSY","SLMSY"}.issubset(set(aprx[motif_col])):
        a = aprx.loc[aprx[motif_col]=="SVMSY","day7_score_num"].dropna()
        b = aprx.loc[aprx[motif_col]=="SLMSY","day7_score_num"].dropna()
        if len(a) and len(b):
            u,p = mannwhitneyu(a,b,alternative="two-sided")
            tests.append({
                "analysis":"SVMSY_vs_SLMSY_day7_unadjusted",
                "n":len(a)+len(b),
                "effect":a.mean()-b.mean(),
                "effect_name":"mean_day7_difference_SVMSY_minus_SLMSY",
                "p_value":p,
                "notes":f"n_SVMSY={len(a)};n_SLMSY={len(b)};U={u}"
            })

    usable = [pat for pat,n in aprx[motif_col].value_counts().items() if pat and n>=3]
    if len(usable)>=2:
        groups = [aprx.loc[aprx[motif_col]==pat,"day7_score_num"].dropna().values for pat in usable]
        H,p = kruskal(*groups)
        tests.append({
            "analysis":"motif_region_variants_vs_day7_unadjusted",
            "n":sum(map(len,groups)),
            "effect":H,
            "effect_name":"Kruskal_H",
            "p_value":p,
            "notes":"groups="+";".join(f"{pat}:{len(g)}" for pat,g in zip(usable,groups))
        })

    variable_species = species_summary.loc[species_summary["n_motif_patterns"]>1,"species"].tolist()
    tests.append({
        "analysis":"within_species_motif_variation_check",
        "n":len(aprx),
        "effect":len(variable_species),
        "effect_name":"species_with_more_than_one_motif_pattern",
        "p_value":np.nan,
        "notes":";".join(variable_species) if variable_species else "none"
    })
    pd.DataFrame(tests).to_csv(res/"lineage_aware_marker_tests.tsv", sep="\t", index=False)

    # Alignment-wide screen restricted to sites that vary within at least one species.
    aln_path = root/"preprint"/"results"/"serralysin_reference_alignment.faa"
    aln = {r.id:str(r.seq).upper() for r in SeqIO.parse(aln_path,"fasta") if r.id.startswith("GCF_")}
    amap = dict(zip(aprx["assembly"], aprx.index))
    rows = []
    ids = []
    species = []
    y = []
    seqs = []
    for sid, seq in aln.items():
        assembly = sid.split("|")[0]
        if assembly not in amap:
            continue
        r = aprx.loc[amap[assembly]]
        ids.append(sid); seqs.append(seq); species.append(r["species"]); y.append(float(r["day7_score_num"]))
    species = np.asarray(species, dtype=object)
    y = np.asarray(y, dtype=float)

    # Species-centred phenotype removes between-species mean differences.
    yres = y.copy()
    for sp in set(species):
        ix = np.where(species==sp)[0]
        yres[ix] = y[ix] - y[ix].mean()

    rng = np.random.default_rng(20260915)
    permutations = []
    for _ in range(args.permutations):
        yp = yres.copy()
        for sp in set(species):
            ix = np.where(species==sp)[0]
            if len(ix)>1:
                yp[ix] = rng.permutation(yp[ix])
        permutations.append(yp)
    permutations = np.asarray(permutations)

    if seqs:
        L = len(seqs[0])
        for col in range(L):
            aa = np.asarray([s[col] for s in seqs], dtype=object)
            valid = np.array([x not in {"-","X","?"} for x in aa])
            if valid.mean() < 0.90:
                continue
            counts = Counter(aa[valid])
            common = [x for x,n in counts.items() if n>=5]
            if len(common)<2:
                continue

            # Must vary within at least one species, otherwise it is pure lineage signal here.
            within_var = False
            variable_sps = []
            for sp in set(species):
                ix = np.where((species==sp) & valid)[0]
                vals = set(aa[ix])
                if len(vals)>=2:
                    within_var = True
                    variable_sps.append(sp)
            if not within_var:
                continue

            obs = eta_squared(aa, yres)
            if not np.isfinite(obs):
                continue
            null = np.array([eta_squared(aa, yp) for yp in permutations])
            null = null[np.isfinite(null)]
            p = (1 + np.sum(null >= obs)) / (1 + len(null))
            rows.append({
                "alignment_column_1based":col+1,
                "allele_counts":";".join(f"{k}:{v}" for k,v in counts.most_common()),
                "species_with_within_species_variation":";".join(sorted(variable_sps)),
                "n_variable_species":len(variable_sps),
                "eta_squared_species_centered":obs,
                "permutation_p":p
            })

    site = pd.DataFrame(rows)
    if not site.empty:
        site["FDR_BH"] = multipletests(site["permutation_p"], method="fdr_bh")[1]
        site = site.sort_values(["FDR_BH","permutation_p"])
    else:
        site = pd.DataFrame(columns=[
            "alignment_column_1based","allele_counts","species_with_within_species_variation",
            "n_variable_species","eta_squared_species_centered","permutation_p","FDR_BH"
        ])
    site.to_csv(res/"within_species_variable_site_tests.tsv", sep="\t", index=False)

    # Figure: phenotype by motif variant
    order = motif_summary[motif_col].tolist()
    plot_groups = []
    labels = []
    for pat in order:
        vals = aprx.loc[aprx[motif_col]==pat,"day7_score_num"].dropna().values
        if len(vals):
            plot_groups.append(vals); labels.append(f"{pat}\n(n={len(vals)})")
    if plot_groups:
        fig,ax = plt.subplots(figsize=(6,4.5))
        ax.boxplot(plot_groups, tick_labels=labels)
        rng2=np.random.default_rng(123)
        for i,g in enumerate(plot_groups,1):
            ax.scatter(rng2.normal(i,0.035,len(g)),g,s=18,alpha=0.7)
        ax.set_ylabel("Day-7 proteolytic activity score (0–4)")
        ax.set_yticks([0,1,2,3,4],["non","weak","moderate","strong","very strong"])
        ax.set_title("AprX motif-region variant versus published phenotype")
        fig.tight_layout()
        fig.savefig(figs/"Fig_motif_variant_vs_day7.png",dpi=300)
        plt.close(fig)

    print("\nSpecies with >1 motif-region variant:", len(variable_species))
    print(species_summary.to_string(index=False))
    print("\nMotif summary:")
    print(motif_summary.to_string(index=False))
    print("\nMarker tests:")
    print(pd.DataFrame(tests).to_string(index=False))
    print("\nWithin-species-variable alignment sites tested:", len(site))
    if not site.empty:
        print("Sites with FDR < 0.05:", int((site["FDR_BH"]<0.05).sum()))
        print(site.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
