import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO, Phylo


APR_X_REFS = [
    "REF_AprX_CY091",
    "REF_AprX_A506",
    "REF_AprX_TSS",
    "REF_AprX_B52",
]

APR_A_REFS = [
    "REF_AprA_PAO1",
    "REF_AprA_PA14",
]


def run(cmd, stdout=None, stderr=None):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def iqtree_exe():
    for x in ("iqtree2", "iqtree"):
        p = shutil.which(x)
        if p:
            return p
    raise RuntimeError("iqtree2/iqtree not found in PATH")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"

    candidates = ext/"external_type_serralysin_candidates.faa"
    refs = [
        root/"preprint"/"references"/"CY091_AprX_clean.faa",
        root/"preprint"/"references"/"A506_AprX.faa",
        root/"preprint"/"references"/"TSS_AprX.faa",
        root/"preprint"/"references"/"B52_AprX.faa",
        root/"preprint"/"references"/"PAO1_AprA_clean.faa",
        root/"preprint"/"references"/"PA14_AprA.faa",
    ]

    panel = ext/"external_phylogeny_panel.faa"
    with open(panel, "w") as out:
        for rec in SeqIO.parse(candidates, "fasta"):
            out.write(f">{rec.id}\n{str(rec.seq)}\n")
        for p in refs:
            if not p.exists():
                print("WARNING missing reference:", p)
                continue
            for rec in SeqIO.parse(p, "fasta"):
                out.write(f">{rec.id}\n{str(rec.seq)}\n")

    nseq = sum(1 for _ in SeqIO.parse(panel, "fasta"))
    print("Sequences in external phylogeny panel:", nseq)

    aln = ext/"external_phylogeny_alignment.faa"
    if not aln.exists() or aln.stat().st_size == 0:
        with open(aln, "w") as oh, open(ext/"mafft_external_phylogeny.log", "w") as eh:
            run(["mafft", "--auto", str(panel)], stdout=oh, stderr=eh)

    prefix = ext/"external_serralysin_tree"
    treefile = Path(str(prefix)+".treefile")

    if not treefile.exists():
        run([
            iqtree_exe(),
            "-s", str(aln),
            "-m", "MFP",
            "-B", "1000",
            "-alrt", "1000",
            "-T", "AUTO",
            "-pre", str(prefix)
        ])

    tree = Phylo.read(treefile, "newick")
    terminals = {t.name: t for t in tree.get_terminals()}

    def resolve(prefixes):
        out = {}
        for pref in prefixes:
            name = next((n for n in terminals if n.startswith(pref)), None)
            if name is None:
                raise RuntimeError(f"Reference {pref} not found in tree")
            out[pref] = name
        return out

    aprx = resolve(APR_X_REFS)
    apra = resolve(APR_A_REFS)

    rows = []
    for name, tip in terminals.items():
        if not name.startswith("GCF_"):
            continue

        fields = name.split("|")
        accession = fields[0]
        protein_id = fields[1] if len(fields) > 1 else ""

        dx = {k: tree.distance(tip, terminals[v]) for k,v in aprx.items()}
        da = {k: tree.distance(tip, terminals[v]) for k,v in apra.items()}

        nearest_x = min(dx, key=dx.get)
        nearest_a = min(da, key=da.get)
        min_x = dx[nearest_x]
        min_a = da[nearest_a]

        rows.append({
            "accession": accession,
            "protein_id": protein_id,
            "nearest_AprX_reference": nearest_x,
            "AprX_distance": min_x,
            "nearest_AprA_reference": nearest_a,
            "AprA_distance": min_a,
            "AprA_minus_AprX_distance": min_a - min_x,
            "nearest_reference_class":
                "AprX-reference-proximal" if min_x < min_a
                else "P.aeruginosa-AprA-reference-proximal"
        })

    out = pd.DataFrame(rows)
    out.to_csv(ext/"external_reference_proximity.tsv", sep="\t", index=False)

    print("\nNearest-reference classes:")
    print(out["nearest_reference_class"].value_counts().to_string())
    print("\nTree:", treefile)
    print("Distances:", ext/"external_reference_proximity.tsv")


if __name__ == "__main__":
    main()
