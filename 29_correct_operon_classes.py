import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def yes(v):
    return str(v).strip().upper() == "YES"


def rescue_apri(hit):
    """Gene-specific rescue for short/divergent AprI homologues."""
    try:
        qcov = float(hit["qcov"])
        pident = float(hit["pident"])
        evalue = float(hit["evalue"])
    except Exception:
        return False

    annot = str(hit.get("annotation_product","")).lower()
    recip = str(hit.get("reciprocal_best_same_gene","")).upper() == "YES"

    annotation_support = (
        "apri" in annot
        or "inh family metalloprotease inhibitor" in annot
        or "protease inhibitor inh" in annot
    )

    return (
        qcov >= 70
        and pident >= 25
        and evalue <= 1e-10
        and (recip or annotation_support)
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    opdir = res/"operon_homology"

    pilot = pd.read_csv(
        res/"AprX_analysis_dataset.tsv",
        sep="\t", dtype=str
    ).fillna("")

    summary = pd.read_csv(
        opdir/"operon_homology_genome_summary.tsv",
        sep="\t", dtype=str
    ).fillna("")
    summary = summary[summary["cohort"]=="pilot"].copy()

    hits = pd.read_csv(
        opdir/"operon_homology_all_hits.tsv",
        sep="\t", dtype=str
    ).fillna("")
    hits = hits[hits["cohort"]=="pilot"].copy()

    # Restrict to the 76 serralysin-positive pilot genomes used in the sequence analysis.
    eligible = pilot[
        pilot["classification"].isin([
            "AprX-reference-like",
            "P.aeruginosa-AprA-like"
        ])
    ].copy()

    merged = eligible.merge(
        summary,
        left_on="assembly",
        right_on="accession",
        how="left",
        suffixes=("","_op")
    ).fillna("")

    apri_hits = hits[hits["gene_label"]=="aprI"].copy()
    apri_rescue = {}
    for _, r in apri_hits.iterrows():
        apri_rescue[r["accession"]] = rescue_apri(r)

    rows = []

    for _, r in merged.iterrows():
        acc = r["assembly"]

        # Local = high-confidence homologue within 50 kb of serralysin.
        local = {}
        present = {}

        for gene in [
            "aprI","aprD","aprE","aprF",
            "prtA_like","prtB_like","lipA1","lipA2"
        ]:
            present[gene] = yes(r.get(gene+"_present",""))
            local[gene] = yes(r.get(gene+"_within_50kb",""))

        # Rescue AprI only, because short inhibitors were clearly being lost by the global E-value threshold.
        apri_rescued = False
        if not present["aprI"] and apri_rescue.get(acc, False):
            present["aprI"] = True

            ah = hits[
                (hits["accession"]==acc) &
                (hits["gene_label"]=="aprI")
            ]

            xh = hits[
                (hits["accession"]==acc) &
                (hits["gene_label"]=="aprX_like_serralysin") &
                (hits["high_confidence_homologue"]=="YES")
            ]

            if not ah.empty and not xh.empty:
                a = ah.iloc[0]
                x = xh.iloc[0]

                try:
                    amid = (float(a["start"]) + float(a["end"])) / 2
                    xmid = (float(x["start"]) + float(x["end"])) / 2

                    local["aprI"] = (
                        a["contig"] == x["contig"] and
                        abs(amid - xmid) <= 50000
                    )
                except Exception:
                    local["aprI"] = False

            apri_rescued = True

        core_local = all(local[g] for g in ["aprI","aprD","aprE","aprF"])
        lipA2_local = local["lipA2"]
        prta_local = local["prtA_like"]
        prtb_local = local["prtB_like"]
        prta_present = present["prtA_like"]
        prtb_present = present["prtB_like"]

        if core_local and lipA2_local and prta_local and prtb_local:
            cls = "Type1-like_local_prtAB"
        elif (
            core_local and lipA2_local
            and prta_present and prtb_present
            and not (prta_local and prtb_local)
        ):
            cls = "Type8-like_distal_prtAB"
        elif (
            core_local and lipA2_local
            and not prta_present and not prtb_present
        ):
            cls = "Type2-like_no_prtAB"
        else:
            cls = "Other_or_incomplete"

        rows.append({
            "assembly": acc,
            "species": r["species"],
            "strain": r["strain"],
            "classification": r["classification"],
            "day7_category": r["day7_category"],
            "day7_score": r["day7_score"],
            "SVMSY_region_aligned": r.get("SVMSY_region_aligned",""),
            "reconstructed_operon_class": cls,
            "aprI_rescued_gene_specific": "YES" if apri_rescued else "NO",
            "aprI_local": "YES" if local["aprI"] else "NO",
            "aprD_local": "YES" if local["aprD"] else "NO",
            "aprE_local": "YES" if local["aprE"] else "NO",
            "aprF_local": "YES" if local["aprF"] else "NO",
            "prtA_present": "YES" if prta_present else "NO",
            "prtA_local": "YES" if prta_local else "NO",
            "prtB_present": "YES" if prtb_present else "NO",
            "prtB_local": "YES" if prtb_local else "NO",
            "lipA2_local": "YES" if lipA2_local else "NO",
        })

    out = pd.DataFrame(rows)
    out.to_csv(
        opdir/"pilot_reconstructed_operon_classes.tsv",
        sep="\t", index=False
    )

    print("Reconstructed pilot operon classes:")
    print(out["reconstructed_operon_class"].value_counts().to_string())

    print("\nAprI gene-specific rescues:")
    print((out["aprI_rescued_gene_specific"]=="YES").sum())

    print("\nClass by species:")
    print(pd.crosstab(out["species"], out["reconstructed_operon_class"]).to_string())


if __name__ == "__main__":
    main()
