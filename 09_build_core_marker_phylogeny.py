import argparse
import csv
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def which_iqtree():
    for cmd in ("iqtree2", "iqtree"):
        p = shutil.which(cmd)
        if p:
            return p
    raise RuntimeError("Neither iqtree2 nor iqtree is in PATH.")


def run(cmd, stdout=None, stderr=None):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout, stderr=stderr)


def parse_tblout(path, evalue_cutoff):
    hits = defaultdict(list)
    with open(path) as h:
        for line in h:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            # HMMER --tblout:
            # target name, accession, query name, accession, E-value, score, bias, ...
            if len(f) < 7:
                continue
            target = f[0]
            query = f[2]
            try:
                evalue = float(f[4])
                score = float(f[5])
            except ValueError:
                continue
            if evalue <= evalue_cutoff:
                hits[query].append((target, evalue, score))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--evalue", type=float, default=1e-10)
    ap.add_argument("--min-presence", type=float, default=0.95)
    ap.add_argument("--max-duplicate-fraction", type=float, default=0.05)
    ap.add_argument("--threads", default="AUTO")
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    work = root/"preprint"/"results"/"auto"/"core_phylogeny"
    logs = root/"preprint"/"logs"/"auto"/"core_phylogeny"
    work.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    master = pd.read_csv(root/"preprint"/"results"/"AprX_master_table.tsv", sep="\t", dtype=str).fillna("")
    assemblies = master["assembly"].tolist()

    hmm = root/"preprint"/"tools"/"bcgTree"/"data"/"ubcg.hmm"
    if not hmm.exists():
        raise FileNotFoundError(f"{hmm} not found. Run setup_phase3.sh first.")

    hmm_dir = work/"hmmsearch"
    hmm_dir.mkdir(exist_ok=True)

    all_hits = {}
    presence = defaultdict(int)
    duplicate = defaultdict(int)

    for i, assembly in enumerate(assemblies, 1):
        proteome = root/"preprint"/"genomes"/"AprX_pilot_genomes"/"ncbi_dataset"/"data"/assembly/"protein.faa"
        if not proteome.exists():
            raise FileNotFoundError(proteome)

        tbl = hmm_dir/f"{assembly}.tbl"
        log = logs/f"{assembly}.hmmsearch.log"

        if not tbl.exists() or tbl.stat().st_size == 0:
            with open(log, "w") as lg:
                run([
                    "hmmsearch", "--noali", "--tblout", str(tbl),
                    str(hmm), str(proteome)
                ], stdout=lg, stderr=subprocess.STDOUT)

        hits = parse_tblout(tbl, args.evalue)
        all_hits[assembly] = hits

        for marker, hs in hits.items():
            presence[marker] += 1
            if len(hs) > 1:
                duplicate[marker] += 1

        print(f"[{i}/{len(assemblies)}] {assembly}: {len(hits)} marker models detected")

    n = len(assemblies)
    markers = sorted(set(presence))
    qc_rows = []
    selected = []

    for marker in markers:
        pres_frac = presence[marker] / n
        dup_frac = duplicate[marker] / n
        keep = pres_frac >= args.min_presence and dup_frac <= args.max_duplicate_fraction
        qc_rows.append({
            "marker": marker,
            "genomes_present": presence[marker],
            "presence_fraction": pres_frac,
            "genomes_with_multiple_hits": duplicate[marker],
            "duplicate_fraction": dup_frac,
            "selected": "YES" if keep else "NO"
        })
        if keep:
            selected.append(marker)

    qc = pd.DataFrame(qc_rows).sort_values(
        ["selected","presence_fraction","marker"],
        ascending=[False,False,True]
    )
    qc.to_csv(work/"core_marker_QC.tsv", sep="\t", index=False)

    if len(selected) < 40:
        raise RuntimeError(
            f"Only {len(selected)} markers passed QC. Stop and inspect core_marker_QC.tsv."
        )

    print(f"\nSelected {len(selected)} core markers for concatenation.")

    # Load each proteome only once.
    proteomes = {}
    for assembly in assemblies:
        proteome = root/"preprint"/"genomes"/"AprX_pilot_genomes"/"ncbi_dataset"/"data"/assembly/"protein.faa"
        proteomes[assembly] = {r.id: str(r.seq) for r in SeqIO.parse(proteome, "fasta")}

    marker_root = work/"markers"
    marker_root.mkdir(exist_ok=True)
    align_root = work/"alignments"
    align_root.mkdir(exist_ok=True)

    alignment_lengths = {}
    best_hit_rows = []

    for marker in selected:
        faa = marker_root/f"{marker}.faa"

        with open(faa, "w") as out:
            for assembly in assemblies:
                hs = all_hits[assembly].get(marker, [])
                if not hs:
                    continue
                hs = sorted(hs, key=lambda x: (-x[2], x[1]))
                target, evalue, score = hs[0]
                seq = proteomes[assembly].get(target)
                if not seq:
                    continue
                out.write(f">{assembly}\n{seq}\n")
                best_hit_rows.append({
                    "assembly": assembly,
                    "marker": marker,
                    "protein_id": target,
                    "evalue": evalue,
                    "bitscore": score,
                    "n_significant_hits": len(hs)
                })

        aln = align_root/f"{marker}.faa"
        if not aln.exists() or aln.stat().st_size == 0:
            with open(aln, "w") as oh, open(logs/f"{marker}.mafft.log", "w") as eh:
                run(["mafft", "--auto", str(faa)], stdout=oh, stderr=eh)

        recs = list(SeqIO.parse(aln, "fasta"))
        if not recs:
            raise RuntimeError(f"No alignment sequences for {marker}")
        alignment_lengths[marker] = len(recs[0].seq)

    pd.DataFrame(best_hit_rows).to_csv(
        work/"core_marker_best_hits.tsv", sep="\t", index=False
    )

    # Concatenate, allowing missing markers as gaps.
    concat = {a: [] for a in assemblies}
    partitions = []
    start = 1

    for marker in selected:
        aln = align_root/f"{marker}.faa"
        d = {r.id: str(r.seq) for r in SeqIO.parse(aln, "fasta")}
        L = alignment_lengths[marker]

        for assembly in assemblies:
            concat[assembly].append(d.get(assembly, "-"*L))

        end = start + L - 1
        partitions.append((marker, start, end))
        start = end + 1

    concat_path = work/"UBCG_core_concatenated.faa"
    with open(concat_path, "w") as out:
        for assembly in assemblies:
            out.write(f">{assembly}\n{''.join(concat[assembly])}\n")

    with open(work/"UBCG_partitions.nex", "w") as out:
        out.write("#nexus\nbegin sets;\n")
        for marker, s, e in partitions:
            safe = "".join(c if c.isalnum() or c=="_" else "_" for c in marker)
            out.write(f"  charset {safe} = {s}-{e};\n")
        out.write("end;\n")

    pd.DataFrame(partitions, columns=["marker","start","end"]).to_csv(
        work/"UBCG_partitions.tsv", sep="\t", index=False
    )

    iqtree = which_iqtree()
    prefix = work/"UBCG_core_tree"

    if not Path(str(prefix)+".treefile").exists():
        run([
            iqtree,
            "-s", str(concat_path),
            "-m", "MFP",
            "-B", "1000",
            "-alrt", "1000",
            "-T", str(args.threads),
            "-pre", str(prefix)
        ])

    print("\nCore phylogeny complete.")
    print("Markers retained:", len(selected))
    print("Concatenated AA alignment length:", start-1)
    print("Tree:", str(prefix)+".treefile")
    print("Marker QC:", work/"core_marker_QC.tsv")


if __name__ == "__main__":
    main()
