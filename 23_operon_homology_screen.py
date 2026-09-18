import argparse
import subprocess
from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote

import pandas as pd
from Bio import SeqIO


MIN_QCOV = 60.0
MIN_PIDENT = 25.0
MAX_EVALUE = 1e-20


def parse_attrs(s):
    d = {}
    for item in s.strip().split(";"):
        if "=" in item:
            k,v = item.split("=",1)
            d[k] = unquote(v)
    return d


def find_gff(folder):
    hits = list(folder.glob("*.gff")) + list(folder.glob("*.gff3"))
    if not hits:
        hits = list(folder.rglob("*.gff")) + list(folder.rglob("*.gff3"))
    return hits[0] if hits else None


def load_gff_protein_map(gff):
    m = {}
    if gff is None:
        return m
    with open(gff, errors="replace") as h:
        for line in h:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue
            a = parse_attrs(f[8])
            pid = a.get("protein_id","")
            if not pid:
                continue
            m[pid] = {
                "contig": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
                "gene": a.get("gene",""),
                "product": a.get("product",""),
                "locus_tag": a.get("locus_tag","")
            }
    return m


def blast_panel(panel, proteome):
    cmd = [
        "blastp",
        "-query", str(panel),
        "-subject", str(proteome),
        "-max_target_seqs", "1",
        "-outfmt", "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore"
    ]
    p = subprocess.run(cmd, text=True, capture_output=True, check=True)
    rows = {}
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) != 9:
            continue
        rows[f[0]] = {
            "qseqid": f[0],
            "protein_id": f[1],
            "pident": float(f[2]),
            "alignment_length": int(f[3]),
            "query_length": int(f[4]),
            "subject_length": int(f[5]),
            "qcov": float(f[6]),
            "evalue": float(f[7]),
            "bitscore": float(f[8])
        }
    return rows


def reciprocal_best(hit_seq, panel, tmp_faa):
    tmp_faa.write_text(f">hit\n{hit_seq}\n")
    cmd = [
        "blastp",
        "-query", str(tmp_faa),
        "-subject", str(panel),
        "-max_target_seqs", "1",
        "-outfmt", "6 sseqid evalue bitscore"
    ]
    p = subprocess.run(cmd, text=True, capture_output=True, check=True)
    lines = p.stdout.strip().splitlines()
    return lines[0].split("\t")[0] if lines else ""


def process_cohort(cohort_name, accessions, data_root, panel, outdir):
    all_rows = []
    summary_rows = []
    tmp = outdir/"_tmp_hit.faa"

    for i, acc in enumerate(accessions,1):
        folder = data_root/acc
        proteome = folder/"protein.faa"
        if not proteome.exists():
            found = list(folder.glob("*protein.faa")) + list(folder.rglob("*protein.faa"))
            if found:
                proteome = found[0]
        if not proteome.exists():
            continue

        seqs = {r.id:str(r.seq).upper() for r in SeqIO.parse(proteome,"fasta")}
        gff = find_gff(folder)
        pos = load_gff_protein_map(gff)
        hits = blast_panel(panel, proteome)

        row = {"cohort":cohort_name, "accession":acc}
        accepted_positions = {}

        for qid, h in hits.items():
            gene = qid.split("|")[0].replace("REF_","")
            seq = seqs.get(h["protein_id"],"")
            rb = reciprocal_best(seq, panel, tmp) if seq else ""
            reciprocal_ok = (rb == qid)
            high_conf = (
                h["qcov"] >= MIN_QCOV and
                h["pident"] >= MIN_PIDENT and
                h["evalue"] <= MAX_EVALUE and
                reciprocal_ok
            )
            info = pos.get(h["protein_id"], {})
            if high_conf and info:
                accepted_positions[gene] = info

            all_rows.append({
                "cohort": cohort_name,
                "accession": acc,
                "gene_label": gene,
                **h,
                "reciprocal_best_reference": rb,
                "reciprocal_best_same_gene": "YES" if reciprocal_ok else "NO",
                "high_confidence_homologue": "YES" if high_conf else "NO",
                "contig": info.get("contig",""),
                "start": info.get("start",""),
                "end": info.get("end",""),
                "strand": info.get("strand",""),
                "annotation_gene": info.get("gene",""),
                "annotation_product": info.get("product",""),
                "locus_tag": info.get("locus_tag","")
            })

            row[gene+"_present"] = "YES" if high_conf else "NO"

        # Relative genomic position to the serralysin reference hit.
        anchor = accepted_positions.get("aprX_like_serralysin")
        for gene in [
            "aprI","aprD","aprE","aprF","prtA_like","prtB_like","lipA1","lipA2"
        ]:
            near = "NO"
            same_contig = "NO"
            distance = ""
            if anchor and gene in accepted_positions:
                other = accepted_positions[gene]
                if other["contig"] == anchor["contig"]:
                    same_contig = "YES"
                    ac = (anchor["start"]+anchor["end"])//2
                    oc = (other["start"]+other["end"])//2
                    distance = oc-ac
                    if abs(distance) <= 50000:
                        near = "YES"
            row[gene+"_same_contig_as_serralysin"] = same_contig
            row[gene+"_within_50kb"] = near
            row[gene+"_distance_bp"] = distance

        present = [
            g for g in [
                "aprX_like_serralysin","aprI","aprD","aprE","aprF",
                "prtA_like","prtB_like","lipA1","lipA2"
            ] if row.get(g+"_present")=="YES"
        ]
        row["high_confidence_gene_set"] = ",".join(present)
        summary_rows.append(row)

        if i % 25 == 0 or i == len(accessions):
            print(f"{cohort_name}: [{i}/{len(accessions)}]")

    if tmp.exists():
        tmp.unlink()

    return pd.DataFrame(all_rows), pd.DataFrame(summary_rows)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()
    root=Path(args.project_root)

    panel=root/"preprint"/"references"/"operon_homology"/"Maier_operon_reference_panel.faa"
    outdir=root/"preprint"/"results"/"auto"/"operon_homology"
    outdir.mkdir(parents=True,exist_ok=True)

    pilot_meta=pd.read_csv(
        root/"preprint"/"results"/"AprX_master_table.tsv",
        sep="\t",dtype=str
    ).fillna("")
    pilot_acc=pilot_meta["assembly"].tolist()
    pilot_root=root/"preprint"/"genomes"/"AprX_pilot_genomes"/"ncbi_dataset"/"data"

    ext_meta=pd.read_csv(
        root/"preprint"/"external_validation"/"pseudomonas_type_refseq_selected.tsv",
        sep="\t",dtype=str
    ).fillna("")
    ext_acc=ext_meta["accession"].tolist()
    ext_root=root/"preprint"/"external_validation"/"type_refseq_genomes"/"ncbi_dataset"/"data"

    p_hits,p_sum=process_cohort("pilot",pilot_acc,pilot_root,panel,outdir)
    e_hits,e_sum=process_cohort("external",ext_acc,ext_root,panel,outdir)

    hits=pd.concat([p_hits,e_hits],ignore_index=True)
    summ=pd.concat([p_sum,e_sum],ignore_index=True)

    hits.to_csv(outdir/"operon_homology_all_hits.tsv",sep="\t",index=False)
    summ.to_csv(outdir/"operon_homology_genome_summary.tsv",sep="\t",index=False)

    print("\nHigh-confidence homologue frequencies:")
    genes=[
        "aprX_like_serralysin","aprI","aprD","aprE","aprF",
        "prtA_like","prtB_like","lipA1","lipA2"
    ]
    for cohort,g in summ.groupby("cohort"):
        print("\n",cohort,"n=",len(g))
        for gene in genes:
            c=gene+"_present"
            print(gene, int((g[c]=="YES").sum()), "/", len(g))

if __name__=="__main__":
    main()
