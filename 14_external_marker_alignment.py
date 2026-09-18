import argparse
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO


MARKER = "SVMSY"


def find_marker_anchor(fasta_path):
    for rec in SeqIO.parse(fasta_path, "fasta"):
        seq = str(rec.seq).upper()
        if MARKER in seq:
            return rec.id, seq
    return None, None


def motif_columns(aln_seq, motif=MARKER):
    ungapped = aln_seq.replace("-", "")
    pos = ungapped.find(motif)
    if pos < 0:
        return None

    target_positions = set(range(pos, pos + len(motif)))
    cols = []
    u = -1

    for col, aa in enumerate(aln_seq):
        if aa != "-":
            u += 1
        if aa != "-" and u in target_positions:
            cols.append(col)

    return cols if len(cols) == len(motif) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"

    external_candidates = ext/"external_type_serralysin_candidates.faa"
    pilot_candidates = root/"preprint"/"results"/"serralysin_candidates.faa"

    if not external_candidates.exists():
        raise FileNotFoundError(external_candidates)
    if not pilot_candidates.exists():
        raise FileNotFoundError(pilot_candidates)

    anchor_id, anchor_seq = find_marker_anchor(pilot_candidates)

    if anchor_seq is None:
        raise RuntimeError(
            "No exact SVMSY-bearing sequence was found in the pilot candidate FASTA."
        )

    print("Using pilot marker anchor:", anchor_id)
    print("Anchor contains exact SVMSY:", MARKER in anchor_seq)

    refs = [
        root/"preprint"/"references"/"CY091_AprX_clean.faa",
        root/"preprint"/"references"/"A506_AprX.faa",
        root/"preprint"/"references"/"TSS_AprX.faa",
        root/"preprint"/"references"/"B52_AprX.faa",
        root/"preprint"/"references"/"PAO1_AprA_clean.faa",
        root/"preprint"/"references"/"PA14_AprA.faa"
    ]

    panel = ext/"external_reference_panel.faa"

    with open(panel, "w") as out:
        # Dedicated exact-SVMSY marker anchor.
        out.write(f">MARKER_ANCHOR_SVMSY|{anchor_id}\n{anchor_seq}\n")

        # External candidates.
        for rec in SeqIO.parse(external_candidates, "fasta"):
            out.write(f">{rec.id}\n{str(rec.seq)}\n")

        # Curated references.
        for p in refs:
            if not p.exists():
                continue
            text = p.read_text()
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")

    aln = ext/"external_reference_alignment.faa"

    print("Running MAFFT on external candidate/reference panel...")
    with open(aln, "w") as oh, open(ext/"mafft_external.log", "w") as eh:
        subprocess.run(
            ["mafft", "--auto", str(panel)],
            stdout=oh,
            stderr=eh,
            check=True
        )

    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln, "fasta")}

    anchor_name = next(
        (k for k in seqs if k.startswith("MARKER_ANCHOR_SVMSY|")),
        None
    )

    if anchor_name is None:
        raise RuntimeError("Marker anchor is missing from the external alignment.")

    cols = motif_columns(seqs[anchor_name])

    if cols is None:
        raise RuntimeError(
            "Could not map SVMSY even from the dedicated exact-SVMSY anchor."
        )

    print("SVMSY alignment columns (1-based):", [c + 1 for c in cols])

    rows = []

    for sid, aseq in seqs.items():
        if not sid.startswith("GCF_"):
            continue

        fields = sid.split("|")
        accession = fields[0]
        protein_id = fields[1] if len(fields) > 1 else ""

        pattern = "".join(
            aseq[c] if c < len(aseq) else "-"
            for c in cols
        )

        rows.append({
            "accession": accession,
            "protein_id": protein_id,
            "SVMSY_aligned_pattern": pattern,
            "exact_SVMSY": "YES" if pattern == MARKER else "NO"
        })

    pat = pd.DataFrame(rows)

    if pat.empty:
        raise RuntimeError("No external GCF candidate sequences were recovered from the alignment.")

    scan = pd.read_csv(
        ext/"external_type_aprx_scan.tsv",
        sep="\t",
        dtype=str
    ).fillna("")

    merged = scan.merge(
        pat,
        on=["accession", "protein_id"],
        how="left"
    )

    merged.to_csv(
        ext/"external_type_aprx_marker_table.tsv",
        sep="\t",
        index=False
    )

    cand = merged[
        merged["status"] == "candidate_serralysin"
    ].copy()

    summary = (
        cand.assign(
            SVMSY_aligned_pattern=cand["SVMSY_aligned_pattern"].replace("", "unresolved")
        )
        .groupby("SVMSY_aligned_pattern")
        .agg(
            n=("accession", "size"),
            species=("organism", "nunique")
        )
        .reset_index()
        .sort_values(["n", "SVMSY_aligned_pattern"], ascending=[False, True])
    )

    summary.to_csv(
        ext/"external_marker_pattern_summary.tsv",
        sep="\t",
        index=False
    )

    unresolved = (
        cand["SVMSY_aligned_pattern"].fillna("").eq("").sum()
        if "SVMSY_aligned_pattern" in cand.columns
        else len(cand)
    )

    print("\nExternal candidates:", len(cand))
    print("Resolved marker patterns:", len(cand) - unresolved)
    print("Unresolved marker patterns:", unresolved)
    print("\nExternal motif patterns:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
