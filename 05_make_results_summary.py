import argparse
from pathlib import Path

import pandas as pd


def fmt(x):
    try:
        v=float(x)
        if pd.isna(v):
            return "NA"
        return f"{v:.4g}"
    except Exception:
        return str(x)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args=ap.parse_args()
    root=Path(args.project_root)
    res=root/"preprint"/"results"/"auto"

    df=pd.read_csv(res/"AprX_analysis_dataset.tsv",sep="\t",dtype=str).fillna("")
    stats=pd.read_csv(res/"statistical_results.tsv",sep="\t",dtype=str).fillna("")
    sites=pd.read_csv(res/"AprX_site_association_day7.tsv",sep="\t",dtype=str).fillna("")

    lines=[]
    lines.append("# Automated AprX preprint analysis summary")
    lines.append("")
    lines.append("This file is generated directly from the current project data. It is a results/QC summary, not manuscript prose.")
    lines.append("")
    lines.append("## Cohort")
    lines.append("")
    lines.append(f"- Total genomes in master dataset: {len(df)}")
    for k,v in df["classification"].value_counts().items():
        lines.append(f"- {k}: {v}")
    matched=(df["day7_category"]!="").sum() if "day7_category" in df else 0
    lines.append(f"- Genomes matched to a Figure 2 phenotype row: {matched}")
    if "SVMSY_exact" in df:
        sub=df[df["classification"]=="AprX-reference-like"]
        lines.append(f"- AprX-reference-like proteins carrying exact SVMSY: {(sub['SVMSY_exact']=='YES').sum()} / {len(sub)}")

    lines.append("")
    lines.append("## Statistical analyses")
    lines.append("")
    if stats.empty:
        lines.append("No statistical analyses were available.")
    else:
        for _,r in stats.iterrows():
            lines.append(
                f"- **{r['analysis']}**: n={r['n']}; {r['effect_name']}={fmt(r['effect'])}; "
                f"p={fmt(r['p_value'])}. {r['notes']}"
            )

    lines.append("")
    lines.append("## Exploratory residue association screen")
    lines.append("")
    lines.append(f"- Variable alignment sites passing the minimum-count rule and tested against day-7 phenotype: {len(sites)}")
    if not sites.empty:
        q=pd.to_numeric(sites["FDR_BH"],errors="coerce")
        lines.append(f"- Sites with FDR < 0.05: {(q<0.05).sum()}")
        lines.append("")
        lines.append("Top ten site tests:")
        lines.append("")
        lines.append("|Alignment col.|CY091 pos.|Groups|p|FDR|")
        lines.append("|---:|---:|---|---:|---:|")
        for _,r in sites.head(10).iterrows():
            lines.append(
                f"|{r['alignment_column_1based']}|{r['CY091_position']}|{r['groups']}|"
                f"{fmt(r['p_value'])}|{fmt(r['FDR_BH'])}|"
            )

    lines.append("")
    lines.append("## Interpretation safeguards")
    lines.append("")
    lines.append("- Figure 2 phenotype extraction must be visually QC-checked once using `Figure2_phenotype_QC.png`.")
    lines.append("- Exact SVMSY presence is evaluated as a proposed marker. Association does not establish causality.")
    lines.append("- Species/lineage can confound residue-phenotype associations. The species-stratified permutation result is therefore more informative than the unadjusted marker test when within-species marker variation exists.")
    lines.append("- Site-wise association tests are exploratory and FDR-corrected.")
    lines.append("- Absence of a convincing AprX-like hit in an assembly is not equivalent to proof of biological gene absence.")

    out=res/"results_summary.md"
    out.write_text("\n".join(lines)+"\n")
    print("Saved:",out)


if __name__=="__main__":
    main()
