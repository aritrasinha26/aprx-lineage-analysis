import argparse
from pathlib import Path
import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()
    root=Path(args.project_root)
    opdir=root/"preprint"/"results"/"auto"/"operon_homology"

    arch=pd.read_csv(opdir/"operon_architectures_final.tsv",sep="\t",dtype=str).fillna("")
    val=pd.read_csv(opdir/"type1_type8_targeted_validation.tsv",sep="\t",dtype=str).fillna("")
    freq=[]

    lines=["# Final computational-analysis freeze checkpoint",""]
    lines.append("## Eligible homology-based operon cohorts")
    lines.append("")
    for cohort,g in arch.groupby("cohort"):
        lines.append(f"- {cohort}: {len(g)} serralysin-positive genomes included.")
        for gene,label in [
            ("aprI","AprI"),("aprD","AprD"),("aprE","AprE"),("aprF","AprF"),
            ("prtA_like","PrtA-like"),("prtB_like","PrtB-like"),("lipA1","LipA1"),("lipA2","LipA2")
        ]:
            n=(g[gene+"_within50kb"]=="YES").sum()
            lines.append(f"  - {label} within 50 kb: {n}/{len(g)} ({100*n/len(g):.1f}%).")
        lines.append("")

    lines.append("## Pilot published-operon validation")
    lines.append("")
    if val.empty:
        lines.append("- No targeted type 1/type 8 validation rows were available.")
    else:
        for _,r in val.iterrows():
            lines.append(
                f"- Published type {r['published_operon_type']}: "
                f"{r['matches_expected_relative_architecture']}/{r['n']} "
                f"({100*float(r['fraction']):.1f}%) matched the expected automated relative architecture."
            )

    lines.append("")
    lines.append("## Freeze rule")
    lines.append("")
    lines.append(
        "The large-scale exploratory analysis is complete. Subsequent changes should be limited to "
        "manual QC, figure presentation, manuscript wording, or correction of a demonstrable pipeline error. "
        "New exploratory analyses should not be added unless they address a specific reviewer-style weakness."
    )

    out=root/"preprint"/"results"/"auto"/"FINAL_ANALYSIS_FREEZE.md"
    out.write_text("\n".join(lines)+"\n")
    print("\n".join(lines))
    print("\nSaved:",out)

if __name__=="__main__":
    main()
