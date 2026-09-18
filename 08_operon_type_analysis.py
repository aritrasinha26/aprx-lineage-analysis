import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import kruskal


def strain_key(s):
    t=str(s).upper().replace("$","5").replace("OSM","DSM")
    for prefix in ["WS","DSM","LMG","ATCC","ICMP","CFBP","CCUG","JCM","NBRC","NCTC","CECT","CIP","KCTC","BCRC","MT"]:
        m=re.search(rf"\b{prefix}\s*[-:]?\s*(\d+)\b",t)
        if m:
            return f"{prefix}{m.group(1)}"
    return ""


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()
    root=Path(args.project_root)
    res=root/"preprint"/"results"/"auto"
    figs=root/"preprint"/"figures"/"auto"

    df=pd.read_csv(res/"AprX_analysis_dataset.tsv",sep="\t",dtype=str).fillna("")
    op=pd.read_csv(res/"Figure2_operon_types_auto.tsv",sep="\t",dtype=str).fillna("")
    op["key"]=op["strain"].map(strain_key)
    op=op[op["key"]!=""].drop_duplicates("key")
    df["key"]=df["strain"].map(strain_key)
    merged=df.merge(op[["key","operon_type"]],on="key",how="left")
    merged["day7_score_num"]=pd.to_numeric(merged["day7_score"],errors="coerce")
    merged.to_csv(res/"AprX_analysis_dataset_with_operon.tsv",sep="\t",index=False)

    x=merged[(merged["operon_type"].fillna("")!="") & merged["day7_score_num"].notna()].copy()
    summary=(x.groupby("operon_type")
        .agg(n=("assembly","size"),
             day7_mean=("day7_score_num","mean"),
             day7_median=("day7_score_num","median"),
             day7_min=("day7_score_num","min"),
             day7_max=("day7_score_num","max"))
        .reset_index())
    summary["type_num"]=pd.to_numeric(summary["operon_type"],errors="coerce")
    summary=summary.sort_values("type_num").drop(columns="type_num")
    summary.to_csv(res/"operon_type_phenotype_summary.tsv",sep="\t",index=False)

    common=summary.loc[summary["n"]>=3,"operon_type"].tolist()
    test=""
    if len(common)>=2:
        groups=[x.loc[x["operon_type"]==t,"day7_score_num"].values for t in common]
        H,p=kruskal(*groups)
        test=f"Kruskal-Wallis across types with n>=3: H={H:.4f}, p={p:.6g}, groups={common}"

    if common:
        plot_groups=[x.loc[x["operon_type"]==t,"day7_score_num"].values for t in common]
        fig,ax=plt.subplots(figsize=(max(6,0.7*len(common)),4.5))
        ax.boxplot(plot_groups,tick_labels=[f"Type {t}\n(n={len(g)})" for t,g in zip(common,plot_groups)])
        rng=np.random.default_rng(42)
        for i,g in enumerate(plot_groups,1):
            ax.scatter(rng.normal(i,0.035,len(g)),g,s=18,alpha=0.7)
        ax.set_ylabel("Day-7 proteolytic activity score (0–4)")
        ax.set_yticks([0,1,2,3,4],["non","weak","moderate","strong","very strong"])
        ax.set_title("Published proteolytic phenotype by aprX-lipA2 operon type")
        fig.tight_layout()
        fig.savefig(figs/"Fig_operon_type_vs_day7.png",dpi=300)
        plt.close(fig)

    print(summary.to_string(index=False))
    if test:
        print("\n"+test)
    print("\nRows with operon type and phenotype:",len(x))


if __name__=="__main__":
    main()
