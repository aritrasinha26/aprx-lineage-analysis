import argparse
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


def parse_attrs(s):
    d = {}
    for item in s.strip().split(";"):
        if "=" in item:
            k,v = item.split("=",1)
            d[k] = unquote(v)
    return d


def find_gff(base):
    cands = list(base.glob("*.gff")) + list(base.glob("*.gff3"))
    if not cands:
        cands = list(base.rglob("*.gff")) + list(base.rglob("*.gff3"))
    return cands[0] if cands else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--window-bp", type=int, default=20000)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"
    scan = pd.read_csv(
        ext/"external_type_aprx_marker_table.tsv",
        sep="\t", dtype=str
    ).fillna("")
    scan = scan[scan["status"]=="candidate_serralysin"].copy()

    data_root = ext/"type_refseq_genomes"/"ncbi_dataset"/"data"

    rows = []

    for i, r in scan.iterrows():
        acc = r["accession"]
        pid = r["protein_id"]
        gff = find_gff(data_root/acc)

        base = {
            "accession": acc,
            "protein_id": pid,
            "gff_found": "YES" if gff else "NO"
        }

        if gff is None:
            base["candidate_found_in_gff"] = "NO"
            rows.append(base)
            continue

        features = []
        candidate = None

        with open(gff, errors="replace") as h:
            for line in h:
                if not line or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) != 9:
                    continue
                seqid, source, ftype, start, end, score, strand, phase, attr = f
                if ftype not in {"CDS","gene"}:
                    continue
                a = parse_attrs(attr)
                rec = {
                    "seqid": seqid,
                    "type": ftype,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attrs": a,
                    "gene": a.get("gene",""),
                    "product": a.get("product",""),
                    "protein_id": a.get("protein_id",""),
                    "locus_tag": a.get("locus_tag","")
                }
                features.append(rec)
                if ftype == "CDS" and (
                    rec["protein_id"] == pid or
                    pid in attr
                ):
                    candidate = rec

        if candidate is None:
            base["candidate_found_in_gff"] = "NO"
            rows.append(base)
            continue

        center = (candidate["start"] + candidate["end"]) // 2
        nearby = [
            x for x in features
            if x["type"]=="CDS"
            and x["seqid"]==candidate["seqid"]
            and x["end"] >= center-args.window_bp
            and x["start"] <= center+args.window_bp
        ]
        nearby = sorted(nearby, key=lambda x:x["start"])

        def blob(x):
            return f"{x['gene']} {x['product']} {x['protein_id']}".lower()

        inhibitor = []
        abc = []
        hlyd = []
        tolc = []
        prt = []
        lip = []

        for x in nearby:
            b = blob(x)
            if re.search(r"\bapri\b|proteinase inhibitor|protease inhibitor|alkaline.*inhibitor", b):
                inhibitor.append(x)
            if (
                "type i secretion" in b and
                ("atp" in b or "abc" in b or "permease" in b)
            ) or re.search(r"\baprd\b", b):
                abc.append(x)
            if "hlyd" in b or "membrane fusion" in b or re.search(r"\bapre\b", b):
                hlyd.append(x)
            if "tolc" in b or re.search(r"\baprf\b", b):
                tolc.append(x)
            if re.search(r"\bprt[ab]\b", b):
                prt.append(x)
            if re.search(r"\blip[a-z0-9]*\b", b) or "lipase" in b:
                lip.append(x)

        def yes(xs):
            return "YES" if xs else "NO"

        neigh_text = " | ".join(
            f"{x['locus_tag']}:{x['gene']}:{x['product']}"
            for x in nearby
        )

        base.update({
            "candidate_found_in_gff": "YES",
            "contig": candidate["seqid"],
            "candidate_start": candidate["start"],
            "candidate_end": candidate["end"],
            "candidate_strand": candidate["strand"],
            "candidate_gene": candidate["gene"],
            "candidate_product": candidate["product"],
            "candidate_locus_tag": candidate["locus_tag"],
            "inhibitor_nearby": yes(inhibitor),
            "T1SS_ABC_or_AprD_nearby": yes(abc),
            "HlyD_or_AprE_nearby": yes(hlyd),
            "TolC_or_AprF_nearby": yes(tolc),
            "prtAB_nearby": yes(prt),
            "lipase_gene_nearby": yes(lip),
            "n_nearby_CDS": len(nearby),
            "neighborhood_text": neigh_text
        })
        rows.append(base)

        if len(rows) % 50 == 0:
            print(f"Processed {len(rows)}/{len(scan)} candidates")

    out = pd.DataFrame(rows)
    out.to_csv(ext/"external_genomic_context.tsv", sep="\t", index=False)

    print("\nExternal candidates with GFF locus found:",
          (out["candidate_found_in_gff"]=="YES").sum(), "/", len(out))
    for c in [
        "inhibitor_nearby","T1SS_ABC_or_AprD_nearby",
        "HlyD_or_AprE_nearby","TolC_or_AprF_nearby",
        "prtAB_nearby","lipase_gene_nearby"
    ]:
        if c in out:
            print(c, "YES:", (out[c]=="YES").sum())


if __name__ == "__main__":
    main()
