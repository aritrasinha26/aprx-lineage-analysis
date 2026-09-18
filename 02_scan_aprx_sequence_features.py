import argparse
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


ZINC_RE = re.compile(r"HE..H..G..H")


def shannon_entropy(chars):
    chars = [c for c in chars if c not in "-X?"]
    if not chars:
        return 0.0
    n = len(chars)
    counts = Counter(chars)
    return -sum((v/n) * math.log2(v/n) for v in counts.values())


def aligned_columns_for_motif(aln_seq, motif):
    ungapped = aln_seq.replace("-", "")
    pos = ungapped.find(motif)
    if pos < 0:
        return None
    target_ungapped_positions = set(range(pos, pos + len(motif)))
    cols = []
    ungapped_i = -1
    for col, aa in enumerate(aln_seq):
        if aa != "-":
            ungapped_i += 1
        if ungapped_i in target_ungapped_positions and aa != "-":
            cols.append(col)
    return cols if len(cols) == len(motif) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    cfg = json.loads((root/"preprint"/"automation"/"config.json").read_text())
    marker = cfg["marker_motif"]

    fasta = root/"preprint"/"results"/"serralysin_candidates.faa"
    aln_fasta = root/"preprint"/"results"/"serralysin_reference_alignment.faa"
    outdir = root/"preprint"/"results"/"auto"
    outdir.mkdir(parents=True, exist_ok=True)

    seqs = {r.id: str(r.seq).upper() for r in SeqIO.parse(fasta, "fasta")}
    if not seqs:
        raise RuntimeError(f"No sequences found in {fasta}")

    rows = []
    for sid, seq in seqs.items():
        assembly = sid.split("|")[0]
        protein_id = sid.split("|")[1] if "|" in sid else ""
        z = ZINC_RE.search(seq)
        rows.append({
            "sequence": sid,
            "assembly": assembly,
            "protein_id": protein_id,
            "protein_length": len(seq),
            "SVMSY_exact": "YES" if marker in seq else "NO",
            "SVMSY_start_1based": seq.find(marker) + 1 if marker in seq else "",
            "zinc_motif": z.group(0) if z else "",
            "zinc_motif_start_1based": z.start()+1 if z else "",
            "zinc_motif_present": "YES" if z else "NO",
        })

    feat = pd.DataFrame(rows)

    aln = {r.id: str(r.seq).upper() for r in SeqIO.parse(aln_fasta, "fasta")}
    candidate_aln = {k:v for k,v in aln.items() if k.startswith("GCF_")}

    motif_cols = None
    motif_anchor = None
    for sid, aseq in candidate_aln.items():
        cols = aligned_columns_for_motif(aseq, marker)
        if cols is not None:
            motif_cols = cols
            motif_anchor = sid
            break

    if motif_cols is not None:
        patterns = {}
        for sid, aseq in candidate_aln.items():
            patterns[sid] = "".join(aseq[c] for c in motif_cols)
        feat["SVMSY_region_aligned"] = feat["sequence"].map(patterns)
    else:
        feat["SVMSY_region_aligned"] = ""
        print("WARNING: no exact SVMSY-bearing sequence was found in the alignment.")

    feat.to_csv(outdir/"AprX_sequence_features.tsv", sep="\t", index=False)

    # Alignment-wide variable sites among candidate proteins.
    ids = list(candidate_aln)
    aln_len = len(next(iter(candidate_aln.values())))
    cy_id = next((k for k in aln if "REF_AprX_CY091" in k), None)
    cy_seq = aln.get(cy_id, "")
    cy_pos = 0
    cy_positions = []
    for aa in cy_seq:
        if aa != "-":
            cy_pos += 1
            cy_positions.append(cy_pos)
        else:
            cy_positions.append("")

    var_rows = []
    for col in range(aln_len):
        chars = [candidate_aln[s][col] for s in ids]
        nongap = [c for c in chars if c not in "-X?"]
        occupancy = len(nongap) / len(chars)
        counts = Counter(nongap)
        if occupancy < 0.90 or len(counts) < 2:
            continue
        common = counts.most_common()
        minor_count = len(nongap) - common[0][1]
        var_rows.append({
            "alignment_column_1based": col+1,
            "CY091_position": cy_positions[col] if cy_positions else "",
            "consensus_residue": common[0][0],
            "consensus_count": common[0][1],
            "minor_count": minor_count,
            "n_residues": len(nongap),
            "n_alleles": len(counts),
            "entropy_bits": shannon_entropy(chars),
            "allele_counts": ";".join(f"{aa}:{n}" for aa,n in common)
        })
    var_df = pd.DataFrame(var_rows)
    if not var_df.empty:
        var_df = var_df.sort_values(["entropy_bits","minor_count"], ascending=False)
    else:
        var_df = pd.DataFrame(columns=[
            "alignment_column_1based","CY091_position","consensus_residue",
            "consensus_count","minor_count","n_residues","n_alleles",
            "entropy_bits","allele_counts"
        ])
    var_df.to_csv(outdir/"AprX_alignment_variable_positions.tsv", sep="\t", index=False)

    print(f"Candidate proteins: {len(feat)}")
    print("Exact SVMSY:", (feat["SVMSY_exact"]=="YES").sum())
    if motif_cols is not None:
        print("SVMSY alignment anchor:", motif_anchor)
        print("SVMSY alignment columns (1-based):", [x+1 for x in motif_cols])
        print("SVMSY-region variants:")
        print(feat["SVMSY_region_aligned"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
