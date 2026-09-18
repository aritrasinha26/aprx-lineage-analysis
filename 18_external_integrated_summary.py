import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo


def fitch_changes(treefile, state_map):
    tree = Phylo.read(treefile, "newick")
    tree.root_at_midpoint()
    changes = 0

    def walk(clade):
        nonlocal changes
        if clade.is_terminal():
            s = state_map.get(clade.name)
            return set() if s is None else {s}

        cs = [walk(c) for c in clade.clades]
        cs = [x for x in cs if x]
        if not cs:
            return set()
        inter = set.intersection(*cs)
        if inter:
            return inter
        changes += 1
        return set.union(*cs)

    walk(tree.root)
    return changes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"

    marker = pd.read_csv(
        ext/"external_type_aprx_marker_table.tsv",
        sep="\t", dtype=str
    ).fillna("")
    prox = pd.read_csv(
        ext/"external_reference_proximity.tsv",
        sep="\t", dtype=str
    ).fillna("")
    context = pd.read_csv(
        ext/"external_genomic_context.tsv",
        sep="\t", dtype=str
    ).fillna("")

    df = marker.merge(
        prox,
        on=["accession","protein_id"],
        how="left"
    ).merge(
        context,
        on=["accession","protein_id"],
        how="left"
    )

    df.to_csv(ext/"external_integrated_serralysin_table.tsv", sep="\t", index=False)

    cand = df[df["status"]=="candidate_serralysin"].copy()

    cross = (
        cand.groupby(["nearest_reference_class","SVMSY_aligned_pattern"])
        .size()
        .reset_index(name="n")
        .sort_values(["nearest_reference_class","n"], ascending=[True,False])
    )
    cross.to_csv(ext/"external_marker_by_reference_class.tsv", sep="\t", index=False)

    ctx_cols = [
        "inhibitor_nearby","T1SS_ABC_or_AprD_nearby",
        "HlyD_or_AprE_nearby","TolC_or_AprF_nearby",
        "prtAB_nearby","lipase_gene_nearby"
    ]

    ctx = []
    for grp, g in cand.groupby("nearest_reference_class"):
        row = {"nearest_reference_class":grp, "n":len(g)}
        for c in ctx_cols:
            if c in g:
                row[c+"_YES"] = int((g[c]=="YES").sum())
        ctx.append(row)
    pd.DataFrame(ctx).to_csv(
        ext/"external_context_by_reference_class.tsv",
        sep="\t", index=False
    )

    # Fitch changes for the 4-state marker pattern on the serralysin tree.
    treefile = ext/"external_serralysin_tree.treefile"
    tree = Phylo.read(treefile, "newick")
    tip_names = {t.name for t in tree.get_terminals()}

    state_map = {}
    for _, r in cand.iterrows():
        prefix = f"{r['accession']}|{r['protein_id']}|"
        name = next((n for n in tip_names if n.startswith(prefix)), None)
        if name and r["SVMSY_aligned_pattern"]:
            state_map[name] = r["SVMSY_aligned_pattern"]

    changes = fitch_changes(treefile, state_map)

    summary = []
    summary.append("# Phase 6 external serralysin classification summary")
    summary.append("")
    summary.append(f"- External candidate serralysins: {len(cand)}")
    summary.append("")
    summary.append("## Reference proximity")
    summary.append("")
    for k,v in cand["nearest_reference_class"].value_counts().items():
        summary.append(f"- {k}: {v}")

    summary.append("")
    summary.append("## Marker pattern by reference-proximity class")
    summary.append("")
    for _, r in cross.iterrows():
        summary.append(
            f"- {r['nearest_reference_class']} / {r['SVMSY_aligned_pattern']}: {int(r['n'])}"
        )

    summary.append("")
    summary.append("## External marker phylogenetic repetition")
    summary.append("")
    summary.append(
        f"- Minimum Fitch changes for the multi-state five-residue marker pattern: {changes} "
        f"across {len(state_map)} candidate tips."
    )

    summary.append("")
    summary.append("## Genomic-context recovery")
    summary.append("")
    found = int((cand["candidate_found_in_gff"]=="YES").sum())
    summary.append(f"- Candidate locus recovered from GFF: {found}/{len(cand)}")
    for c in ctx_cols:
        if c in cand:
            summary.append(f"- {c} YES: {int((cand[c]=='YES').sum())}/{len(cand)}")

    summary.append("")
    summary.append(
        "Interpretation: nearest-reference class is an evolutionary proximity descriptor, "
        "not a definitive functional gene name. External type genomes provide taxonomic and "
        "genomic-context validation but not independent spoilage phenotype validation."
    )

    out = ext/"phase6_external_summary.md"
    out.write_text("\n".join(summary)+"\n")
    print("\n".join(summary))
    print("\nSaved:", out)


if __name__ == "__main__":
    main()
