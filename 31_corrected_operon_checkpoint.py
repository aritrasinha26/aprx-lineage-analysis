import argparse
from pathlib import Path
import pandas as pd


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()

    root=Path(args.project_root)
    res=root/"preprint"/"results"/"auto"
    opdir=res/"operon_homology"

    cls=pd.read_csv(
        opdir/"pilot_reconstructed_operon_classes.tsv",
        sep="\t",dtype=str
    ).fillna("")
    testfile = opdir/"corrected_operon_phenotype_tests.tsv"
    try:
        tests=pd.read_csv(testfile,sep="\t",dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        tests=pd.DataFrame()

    lines=[
        "# Corrected operon-analysis checkpoint",
        "",
        "The OCR-derived individual operon-type assignments are no longer used as the primary inferential variable. "
        "A targeted QC showed row-level OCR misassignment for some published type-8 labels. "
        "Operon classes below are reconstructed directly from homology and genomic position.",
        "",
        "## Reconstructed pilot classes",
        ""
    ]

    for k,v in cls["reconstructed_operon_class"].value_counts().items():
        lines.append(f"- {k}: {v}")

    lines += ["","## Corrected phenotype tests",""]

    for _,r in tests.iterrows():
        lines.append(
            f"- {r['analysis']}: n={r['n']}; {r['effect_name']}={r['effect']}; "
            f"p={r['p_value']}. {r['notes']}"
        )

    lines += [
        "",
        "## Interpretation",
        "",
        "Type-8-like reconstruction allows PrtA-like and PrtB-like homologues to be distal from the local aprX/aprI/aprDEF/lipA2 cluster, including on a different assembly contig. "
        "AprI uses a gene-specific rescue rule because the short inhibitor homologues in P. lundensis were annotated as AprI/Inh proteins but failed the earlier global E-value threshold.",
        "",
        "The previously reported OCR-based type-1 versus type-8 statistical result should be replaced by the reconstructed-operon result above."
    ]

    out=res/"CORRECTED_OPERON_ANALYSIS.md"
    out.write_text("\n".join(lines)+"\n")
    print("\n".join(lines))
    print("\nSaved:",out)


if __name__=="__main__":
    main()
