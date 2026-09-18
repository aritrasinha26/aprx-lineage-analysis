import argparse
from pathlib import Path

import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()
    root=Path(args.project_root)

    outdir=root/"preprint"/"results"/"auto"/"operon_homology"
    s=pd.read_csv(outdir/"operon_homology_genome_summary.tsv",sep="\t",dtype=str).fillna("")

    genes=[
        "aprI","aprD","aprE","aprF","prtA_like","prtB_like","lipA1","lipA2"
    ]

    rows=[]
    for cohort,g in s.groupby("cohort"):
        n=len(g)
        for gene in genes:
            present=int((g[gene+"_present"]=="YES").sum())
            near=int((g[gene+"_within_50kb"]=="YES").sum())
            rows.append({
                "cohort":cohort,
                "gene":gene,
                "n_genomes":n,
                "present_anywhere":present,
                "present_anywhere_pct":100*present/n,
                "within_50kb_of_serralysin":near,
                "within_50kb_pct":100*near/n
            })

    out=pd.DataFrame(rows)
    out.to_csv(outdir/"operon_homology_frequency_summary.tsv",sep="\t",index=False)

    lines=["# Phase 8 homology-based operon reconstruction summary",""]
    for cohort,g in out.groupby("cohort"):
        lines.append(f"## {cohort.capitalize()} cohort")
        lines.append("")
        for _,r in g.iterrows():
            lines.append(
                f"- {r['gene']}: present anywhere {int(r['present_anywhere'])}/{int(r['n_genomes'])} "
                f"({r['present_anywhere_pct']:.1f}%); within 50 kb of serralysin "
                f"{int(r['within_50kb_of_serralysin'])}/{int(r['n_genomes'])} "
                f"({r['within_50kb_pct']:.1f}%)."
            )
        lines.append("")

    lines.append(
        "The homology calls use the nine Pseudomonas protegens CHA0 reference proteins "
        "AGL85002.1–AGL85010.1 with best-hit, coverage, identity, E-value and reciprocal-best-reference filtering. "
        "They should be treated as automated high-confidence homology calls rather than manually curated orthology assignments."
    )

    summary=outdir/"phase8_operon_homology_summary.md"
    summary.write_text("\n".join(lines)+"\n")
    print("\n".join(lines))
    print("\nSaved:",summary)

if __name__=="__main__":
    main()
