import argparse
from pathlib import Path
import pandas as pd


LEVEL_RANK = {
    "Complete Genome": 4,
    "Chromosome": 3,
    "Scaffold": 2,
    "Contig": 1
}

REFSEQ_RANK = {
    "reference genome": 3,
    "representative genome": 2,
    "na": 1,
    "": 1
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"
    inp = ext/"pseudomonas_type_refseq_metadata.tsv"

    df = pd.read_csv(inp, sep="\t", dtype=str).fillna("")
    rename = {
        "Assembly Accession": "accession",
        "Organism Name": "organism",
        "Organism Infraspecific Names Strain": "strain",
        "Assembly Level": "assembly_level",
        "Assembly Refseq Category": "refseq_category",
        "Assembly Stats Number of Contigs": "contigs",
        "Assembly Stats Contig N50": "contig_n50",
        "Type Material Label": "type_material_label",
        "Type Material Display Text": "type_material"
    }
    df = df.rename(columns={k:v for k,v in rename.items() if k in df.columns})

    required = ["accession","organism"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Required columns missing from metadata TSV: {missing}")

    for c in ["strain","assembly_level","refseq_category","contigs","contig_n50",
              "type_material_label","type_material"]:
        if c not in df.columns:
            df[c] = ""

    df = df[df["accession"].str.startswith("GCF_")].copy()
    df = df[df["organism"].str.strip() != ""].copy()

    df["level_rank"] = df["assembly_level"].map(LEVEL_RANK).fillna(0)
    df["refseq_rank"] = (
        df["refseq_category"].str.lower().map(REFSEQ_RANK).fillna(1)
    )
    df["contigs_num"] = pd.to_numeric(df["contigs"], errors="coerce").fillna(10**9)
    df["n50_num"] = pd.to_numeric(df["contig_n50"], errors="coerce").fillna(0)

    # One current type-material RefSeq assembly per taxonomic organism name.
    # Preference: RefSeq category, assembly level, fewer contigs, larger N50.
    df = df.sort_values(
        ["organism","refseq_rank","level_rank","contigs_num","n50_num"],
        ascending=[True,False,False,True,False]
    )
    selected = df.drop_duplicates("organism", keep="first").copy()

    selected.to_csv(ext/"pseudomonas_type_refseq_selected.tsv", sep="\t", index=False)
    selected["accession"].to_csv(
        ext/"pseudomonas_type_refseq_accessions.txt",
        index=False, header=False
    )

    print("Type-material RefSeq records returned:", len(df))
    print("Unique organism names selected:", len(selected))
    print("\nAssembly levels:")
    print(selected["assembly_level"].value_counts(dropna=False).to_string())
    print("\nSaved:")
    print(ext/"pseudomonas_type_refseq_selected.tsv")
    print(ext/"pseudomonas_type_refseq_accessions.txt")


if __name__ == "__main__":
    main()
