import argparse
import urllib.request
from pathlib import Path
from Bio import SeqIO
from io import StringIO

REFS = [
    ("aprX_like_serralysin", "AGL85010.1"),
    ("aprI", "AGL85009.1"),
    ("aprD", "AGL85008.1"),
    ("aprE", "AGL85007.1"),
    ("aprF", "AGL85006.1"),
    ("lipA1", "AGL85005.1"),
    ("prtA_like", "AGL85004.1"),
    ("prtB_like", "AGL85003.1"),
    ("lipA2", "AGL85002.1"),
]

def fetch_fasta(accession):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={accession}&rettype=fasta&retmode=text"
    )
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    outdir = root/"preprint"/"references"/"operon_homology"
    outdir.mkdir(parents=True, exist_ok=True)

    panel = outdir/"Maier_operon_reference_panel.faa"
    meta = outdir/"Maier_operon_reference_mapping.tsv"

    mapping_rows = []
    with open(panel, "w") as out:
        for gene, acc in REFS:
            print("Fetching", gene, acc)
            txt = fetch_fasta(acc)
            recs = list(SeqIO.parse(StringIO(txt), "fasta"))
            if len(recs) != 1:
                raise RuntimeError(f"Expected one record for {acc}, got {len(recs)}")
            rec = recs[0]
            seq = str(rec.seq)
            out.write(f">REF_{gene}|{acc}\n{seq}\n")
            mapping_rows.append((gene, acc, rec.description, len(seq)))

    with open(meta, "w") as h:
        h.write("gene_label\taccession\tNCBI_description\tprotein_length\n")
        for row in mapping_rows:
            h.write("\t".join(map(str,row))+"\n")

    print("\nSaved:", panel)
    print("Saved:", meta)

if __name__ == "__main__":
    main()
