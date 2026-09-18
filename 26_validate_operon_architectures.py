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
    arch=arch[arch["cohort"]=="pilot"].copy()

    pilot=pd.read_csv(
        root/"preprint"/"results"/"auto"/"AprX_analysis_dataset_with_operon.tsv",
        sep="\t",dtype=str
    ).fillna("")

    merged=pilot.merge(
        arch,
        left_on="assembly",
        right_on="accession",
        how="left"
    ).fillna("")

    usable=merged[
        (merged["operon_type"].astype(str).str.strip()!="") &
        (merged["architecture_100kb"].astype(str).str.strip()!="")
    ].copy()

    rows=[]
    for typ,g in usable.groupby("operon_type"):
        vc=g["architecture_100kb"].value_counts()
        dominant=vc.index[0]
        rows.append({
            "published_operon_type":typ,
            "n_with_reconstructed_architecture":len(g),
            "n_distinct_reconstructed_architectures":len(vc),
            "dominant_reconstructed_architecture":dominant,
            "dominant_n":int(vc.iloc[0]),
            "dominant_fraction":float(vc.iloc[0]/len(g))
        })

    summary=pd.DataFrame(rows)
    if not summary.empty:
        summary["type_num"]=pd.to_numeric(summary["published_operon_type"],errors="coerce")
        summary=summary.sort_values("type_num").drop(columns="type_num")
    summary.to_csv(opdir/"published_type_vs_reconstructed_architecture.tsv",sep="\t",index=False)

    # Targeted checks for the two major types used in phenotype analysis.
    checks=[]
    for typ in ["1","8"]:
        g=usable[usable["operon_type"]==typ].copy()
        if g.empty:
            continue
        if typ=="1":
            # Expected relative order for published type 1:
            # aprXIDEF prtAB lipA2
            def ok(s):
                tokens=["aprX","aprI","aprD","aprE","aprF","prtA-like","prtB-like","lipA2"]
                pos=[s.find(t) for t in tokens]
                return all(p>=0 for p in pos) and pos==sorted(pos) and "|" not in s
        else:
            # Published type 8: aprXIDEF lipA2 | prtAB.
            def ok(s):
                core=["aprX","aprI","aprD","aprE","aprF","lipA2"]
                pos=[s.find(t) for t in core]
                core_ok=all(p>=0 for p in pos) and pos==sorted(pos)
                ab_ok=("prtA-like" in s and "prtB-like" in s and s.find("prtA-like")<s.find("prtB-like"))
                separated="|" in s
                return core_ok and ab_ok and separated
        vals=g["architecture_100kb"].map(ok)
        checks.append({
            "published_operon_type":typ,
            "n":len(g),
            "matches_expected_relative_architecture":int(vals.sum()),
            "fraction":float(vals.mean())
        })

    checkdf=pd.DataFrame(checks)
    checkdf.to_csv(opdir/"type1_type8_targeted_validation.tsv",sep="\t",index=False)

    print("Published operon type vs reconstructed architecture:")
    print(summary.to_string(index=False))
    print("\nTargeted type 1 / type 8 checks:")
    print(checkdf.to_string(index=False) if not checkdf.empty else "No usable type 1/8 rows.")

if __name__=="__main__":
    main()
