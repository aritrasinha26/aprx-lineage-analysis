import argparse
import csv
import re
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


ZINC_RE = re.compile(r"HE..H..G..H")


def blast_best(query, subject):
    cmd = [
        "blastp",
        "-query", str(query),
        "-subject", str(subject),
        "-max_target_seqs", "1",
        "-outfmt", "6 sseqid pident length qlen slen qcovs evalue bitscore"
    ]
    p = subprocess.run(cmd, text=True, capture_output=True, check=True)
    line = p.stdout.strip().splitlines()
    if not line:
        return None
    f = line[0].split("\t")
    return {
        "protein_id": f[0],
        "pident": float(f[1]),
        "alignment_length": int(f[2]),
        "query_length": int(f[3]),
        "subject_length": int(f[4]),
        "qcov": float(f[5]),
        "evalue": float(f[6]),
        "bitscore": float(f[7])
    }


def find_protein(proteome, protein_id):
    for rec in SeqIO.parse(proteome, "fasta"):
        if rec.id == protein_id:
            return str(rec.seq).upper()
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"
    unpacked = ext/"type_refseq_genomes"/"ncbi_dataset"/"data"
    meta = pd.read_csv(
        ext/"pseudomonas_type_refseq_selected.tsv",
        sep="\t", dtype=str
    ).fillna("")

    query = root/"preprint"/"references"/"CY091_AprX_clean.faa"
    if not query.exists():
        query = root/"preprint"/"references"/"CY091_AprX_AAC38255.2.faa"
    if not query.exists():
        raise FileNotFoundError("CY091 AprX reference FASTA not found.")

    rows = []
    candidate_records = []

    for i, r in meta.iterrows():
        acc = r["accession"]
        d = unpacked/acc
        proteome = d/"protein.faa"

        if not proteome.exists():
            # NCBI packages can occasionally nest/rename files; find a protein FASTA.
            found = list(d.glob("*protein.faa"))
            if found:
                proteome = found[0]

        if not proteome.exists():
            rows.append({
                "accession": acc,
                "organism": r.get("organism",""),
                "strain": r.get("strain",""),
                "status": "protein_fasta_missing"
            })
            continue

        hit = blast_best(query, proteome)
        if hit is None:
            rows.append({
                "accession": acc,
                "organism": r.get("organism",""),
                "strain": r.get("strain",""),
                "status": "no_blast_hit"
            })
            continue

        seq = find_protein(proteome, hit["protein_id"])
        motif = ZINC_RE.search(seq) if seq else None

        convincing = (
            hit["qcov"] >= 90 and
            hit["evalue"] <= 1e-50 and
            430 <= len(seq) <= 550 and
            motif is not None
        )

        row = {
            "accession": acc,
            "organism": r.get("organism",""),
            "strain": r.get("strain",""),
            **hit,
            "protein_length": len(seq) if seq else "",
            "zinc_motif": motif.group(0) if motif else "",
            "zinc_motif_present": "YES" if motif else "NO",
            "status": "candidate_serralysin" if convincing else "no_convincing_AprX_like_candidate"
        }
        rows.append(row)

        if convincing:
            candidate_records.append(
                (f"{acc}|{hit['protein_id']}|{r.get('organism','').replace(' ','_')}", seq)
            )

        if (i+1) % 25 == 0 or i+1 == len(meta):
            print(f"[{i+1}/{len(meta)}] scanned")

    out = pd.DataFrame(rows)
    out.to_csv(ext/"external_type_aprx_scan.tsv", sep="\t", index=False)

    with open(ext/"external_type_serralysin_candidates.faa", "w") as h:
        for name, seq in candidate_records:
            h.write(f">{name}\n{seq}\n")

    print("\nExternal type genomes:", len(out))
    print("Convincing AprX-like/serralysin candidates:",
          (out["status"]=="candidate_serralysin").sum())
    print("No convincing candidate:",
          (out["status"]=="no_convincing_AprX_like_candidate").sum())
    print("Candidate FASTA:", ext/"external_type_serralysin_candidates.faa")


if __name__ == "__main__":
    main()
