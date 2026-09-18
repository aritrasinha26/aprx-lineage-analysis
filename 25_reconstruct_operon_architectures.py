import argparse
from pathlib import Path
import numpy as np
import pandas as pd


GENE_ORDER = [
    "aprX_like_serralysin","aprI","aprD","aprE","aprF",
    "prtA_like","prtB_like","lipA1","lipA2"
]

DISPLAY = {
    "aprX_like_serralysin":"aprX",
    "aprI":"aprI",
    "aprD":"aprD",
    "aprE":"aprE",
    "aprF":"aprF",
    "prtA_like":"prtA-like",
    "prtB_like":"prtB-like",
    "lipA1":"lipA1",
    "lipA2":"lipA2",
}

def midpoint(row):
    return (float(row["start"]) + float(row["end"])) / 2

def oriented_gap(a, b, strand):
    # Absolute intergenic distance between adjacent features in genomic coordinates.
    a1,a2 = sorted([float(a["start"]), float(a["end"])])
    b1,b2 = sorted([float(b["start"]), float(b["end"])])
    if b1 >= a2:
        return b1 - a2
    if a1 >= b2:
        return a1 - b2
    return 0.0

def architecture_for_genome(g, gap_separator_bp=10000, window_bp=100000):
    g = g[g["high_confidence_homologue"]=="YES"].copy()
    g["start_num"] = pd.to_numeric(g["start"], errors="coerce")
    g["end_num"] = pd.to_numeric(g["end"], errors="coerce")
    g = g.dropna(subset=["start_num","end_num"])

    anchor = g[g["gene_label"]=="aprX_like_serralysin"]
    if anchor.empty:
        return None
    anchor = anchor.iloc[0]
    contig = anchor["contig"]
    strand = anchor["strand"]
    ac = (anchor["start_num"] + anchor["end_num"]) / 2

    same = g[g["contig"]==contig].copy()
    same["mid"] = (same["start_num"] + same["end_num"]) / 2
    same["distance_from_anchor"] = same["mid"] - ac
    same = same[same["distance_from_anchor"].abs() <= window_bp].copy()

    ascending = strand != "-"
    same = same.sort_values("mid", ascending=ascending)

    toks = []
    details = []
    prev = None
    for _, r in same.iterrows():
        gene = r["gene_label"]
        if gene not in DISPLAY:
            continue
        label = DISPLAY[gene]
        if prev is not None:
            gap = oriented_gap(prev, r, strand)
            toks.append(" | " if gap > gap_separator_bp else "-")
        toks.append(label)
        details.append(
            f"{label}:{r['contig']}:{int(r['start_num'])}-{int(r['end_num'])}:{r['strand']}"
        )
        prev = r

    distal = g[
        (g["contig"] != contig) |
        (((g["start_num"]+g["end_num"])/2 - ac).abs() > window_bp)
    ]["gene_label"].map(DISPLAY).dropna().tolist()

    return {
        "anchor_contig": contig,
        "anchor_strand": strand,
        "architecture_100kb": "".join(toks),
        "distal_homologues": ",".join(distal),
        "feature_details_100kb": " ; ".join(details),
        "n_features_100kb": len(details),
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--gap-separator-bp", type=int, default=10000)
    args = ap.parse_args()

    root = Path(args.project_root)
    opdir = root/"preprint"/"results"/"auto"/"operon_homology"
    hits = pd.read_csv(opdir/"operon_homology_all_hits.tsv", sep="\t", dtype=str).fillna("")

    pilot = pd.read_csv(
        root/"preprint"/"results"/"auto"/"AprX_analysis_dataset.tsv",
        sep="\t", dtype=str
    ).fillna("")
    pilot_keep = set(
        pilot.loc[
            pilot["classification"].isin(["AprX-reference-like","P.aeruginosa-AprA-like"]),
            "assembly"
        ]
    )

    external = pd.read_csv(
        root/"preprint"/"external_validation"/"external_integrated_serralysin_table.tsv",
        sep="\t", dtype=str
    ).fillna("")
    external_keep = set(
        external.loc[external["status"]=="candidate_serralysin","accession"]
    )

    eligible = {
        "pilot": pilot_keep,
        "external": external_keep
    }

    rows = []
    for (cohort, accession), g in hits.groupby(["cohort","accession"]):
        if cohort not in eligible or accession not in eligible[cohort]:
            continue
        rec = architecture_for_genome(
            g,
            gap_separator_bp=args.gap_separator_bp,
            window_bp=100000
        )
        if rec is None:
            continue
        rec.update({"cohort":cohort, "accession":accession})
        # Presence and within-50kb flags from the original Phase 8 summary.
        anchor_row = g[g["gene_label"]=="aprX_like_serralysin"]
        if not anchor_row.empty:
            anchor_row = anchor_row.iloc[0]
            ac = (float(anchor_row["start"]) + float(anchor_row["end"])) / 2 if anchor_row["start"] and anchor_row["end"] else None
            contig = anchor_row["contig"]
            for gene in GENE_ORDER[1:]:
                q = g[(g["gene_label"]==gene) & (g["high_confidence_homologue"]=="YES")]
                present = not q.empty
                within50 = False
                if present and ac is not None:
                    qr = q.iloc[0]
                    if qr["contig"]==contig and qr["start"] and qr["end"]:
                        qc = (float(qr["start"])+float(qr["end"]))/2
                        within50 = abs(qc-ac) <= 50000
                rec[gene+"_present"] = "YES" if present else "NO"
                rec[gene+"_within50kb"] = "YES" if within50 else "NO"
        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv(opdir/"operon_architectures_final.tsv", sep="\t", index=False)

    print("Eligible reconstructed architectures:")
    print(out["cohort"].value_counts().to_string())
    print("\nTop pilot architectures:")
    print(out[out["cohort"]=="pilot"]["architecture_100kb"].value_counts().head(15).to_string())
    print("\nTop external architectures:")
    print(out[out["cohort"]=="external"]["architecture_100kb"].value_counts().head(15).to_string())

if __name__=="__main__":
    main()
