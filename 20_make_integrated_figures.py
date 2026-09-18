import argparse
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo


def tip_y_positions(tree):
    terminals = tree.get_terminals()
    return {t.name:i for i,t in enumerate(terminals)}


def plot_tree_with_tracks(treefile, meta, outfile, title,
                          species_col=None,
                          motif_col=None,
                          phenotype_col=None,
                          operon_col=None):
    tree = Phylo.read(treefile, "newick")
    terminals = tree.get_terminals()
    ypos = tip_y_positions(tree)

    maxdist = max(tree.distance(tree.root, t) for t in terminals)
    n = len(terminals)

    fig_h = max(8, n*0.085)
    fig, ax = plt.subplots(figsize=(15, fig_h))
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False, label_func=lambda x: None)

    ax.set_title(title)
    ax.set_xlabel("Branch length")
    ax.set_ylabel("")

    meta = meta.set_index("tip_name", drop=False)

    # track positions
    x0 = maxdist * 1.04
    dx = maxdist * 0.055

    motif_vals = sorted(v for v in meta[motif_col].unique() if v) if motif_col else []
    motif_map = {v:i for i,v in enumerate(motif_vals)}

    phenotype_order = {"non":0,"weak":1,"moderate":2,"strong":3,"very_strong":4}

    for tip in terminals:
        name = tip.name
        if name not in meta.index:
            continue
        r = meta.loc[name]
        y = ypos[name]

        if motif_col:
            v = r[motif_col]
            if v:
                ax.scatter(x0, y, s=18, marker="s", c=[motif_map.get(v,0)], cmap="tab10",
                           vmin=0, vmax=max(1,len(motif_vals)-1))

        if phenotype_col:
            v = r[phenotype_col]
            if v in phenotype_order:
                ax.scatter(x0+dx, y, s=18, marker="o", c=[phenotype_order[v]],
                           cmap="viridis", vmin=0, vmax=4)

        if operon_col:
            v = str(r[operon_col]).strip()
            if v:
                try:
                    num = int(float(v))
                    ax.scatter(x0+2*dx, y, s=18, marker="^", c=[num],
                               cmap="plasma", vmin=1, vmax=22)
                except Exception:
                    pass

    ax.text(x0, n+1, "Marker", rotation=90, va="bottom", ha="center")
    if phenotype_col:
        ax.text(x0+dx, n+1, "Day-7 phenotype", rotation=90, va="bottom", ha="center")
    if operon_col:
        ax.text(x0+2*dx, n+1, "Operon type", rotation=90, va="bottom", ha="center")

    ax.set_xlim(left=0, right=x0 + (3.2*dx if operon_col else 2.2*dx))
    fig.tight_layout()
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root/"preprint"/"results"/"auto"
    ext = root/"preprint"/"external_validation"
    figs = root/"preprint"/"figures"/"final"
    figs.mkdir(parents=True, exist_ok=True)

    # Pilot core tree metadata
    core_meta = pd.read_csv(
        res/"core_tree_metadata.tsv",
        sep="\t", dtype=str
    ).fillna("")
    core_meta["tip_name"] = core_meta["assembly"]

    pilot_tree = res/"core_phylogeny"/"UBCG_core_tree.treefile"
    plot_tree_with_tracks(
        pilot_tree,
        core_meta,
        figs/"Fig_core_tree_marker_phenotype_operon.png",
        "Core-genome phylogeny of the 87-genome pilot cohort",
        motif_col="SVMSY_region_aligned",
        phenotype_col="day7_category",
        operon_col="operon_type" if "operon_type" in core_meta.columns else None
    )

    # External tree metadata
    ex = pd.read_csv(
        ext/"external_integrated_serralysin_table.tsv",
        sep="\t", dtype=str
    ).fillna("")
    ex = ex[ex["status"]=="candidate_serralysin"].copy()

    # Reconstruct actual tree tip names from the FASTA convention.
    ex["tip_name"] = (
        ex["accession"] + "|" + ex["protein_id"] + "|" +
        ex["organism"].str.replace(" ","_",regex=False)
    )

    external_tree = ext/"external_serralysin_tree.treefile"
    plot_tree_with_tracks(
        external_tree,
        ex,
        figs/"Fig_external_serralysin_tree_marker.png",
        "External type-genome serralysin phylogeny",
        motif_col="SVMSY_aligned_pattern"
    )

    # Reference class x marker stacked bar
    cross = pd.crosstab(
        ex["nearest_reference_class"],
        ex["SVMSY_aligned_pattern"]
    )
    order = [x for x in ["SVMSY","SLMSY","SIMSY","TVMSY"] if x in cross.columns]
    cross = cross[order]
    frac = cross.div(cross.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(8,4.5))
    frac.plot(kind="barh", stacked=True, ax=ax)
    ax.set_xlabel("Proportion of candidate serralysins")
    ax.set_ylabel("")
    ax.set_xlim(0,1)
    ax.set_title("Five-residue marker pattern by reference-proximity class")
    ax.legend(title="Marker pattern", bbox_to_anchor=(1.02,1), loc="upper left")
    fig.tight_layout()
    fig.savefig(figs/"Fig_external_marker_by_reference_class.png", dpi=300)
    plt.close(fig)

    # Context conservation
    ctx_cols = [
        "inhibitor_nearby",
        "T1SS_ABC_or_AprD_nearby",
        "HlyD_or_AprE_nearby",
        "TolC_or_AprF_nearby",
        "lipase_gene_nearby"
    ]
    ctx_labels = [
        "Protease inhibitor",
        "T1SS ABC/AprD",
        "HlyD/AprE",
        "TolC/AprF",
        "Lipase"
    ]

    vals=[]
    for c in ctx_cols:
        vals.append((ex[c]=="YES").mean()*100)

    fig,ax=plt.subplots(figsize=(7,4.5))
    ax.barh(ctx_labels,vals)
    ax.set_xlabel("Candidates with feature nearby (%)")
    ax.set_xlim(0,100)
    ax.set_title("Conservation of serralysin genomic context")
    for y,v in enumerate(vals):
        ax.text(v+1,y,f"{v:.1f}%",va="center")
    fig.tight_layout()
    fig.savefig(figs/"Fig_external_genomic_context_summary.png",dpi=300)
    plt.close(fig)

    print("Created:")
    for p in sorted(figs.glob("*.png")):
        print(p)


if __name__ == "__main__":
    main()
