import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo


MOTIF_COLORS = {
    "SVMSY":"#1f77b4",
    "SLMSY":"#ff7f0e",
    "SIMSY":"#2ca02c",
    "TVMSY":"#d62728",
}
PHENO_COLORS = {
    "non":"#2b7bba",
    "weak":"#7fb2bd",
    "moderate":"#f6f0b8",
    "strong":"#ef9d6e",
    "very_strong":"#d62728",
}

def node_positions(tree, cladogram=False):
    terminals=tree.get_terminals()
    y={t: i for i,t in enumerate(terminals)}
    x={}
    def setx(clade,depth=0,dist=0.0):
        x[clade]=depth if cladogram else dist
        for child in clade.clades:
            bl=child.branch_length or 0.0
            setx(child,depth+1,dist+bl)
    setx(tree.root)
    def sety(clade):
        if clade in y:
            return y[clade]
        vals=[sety(c) for c in clade.clades]
        y[clade]=sum(vals)/len(vals)
        return y[clade]
    sety(tree.root)
    return x,y,terminals

def draw_rect_tree(ax,tree,cladogram=False):
    x,y,terms=node_positions(tree,cladogram)
    for clade in tree.find_clades(order="preorder"):
        if not clade.clades:
            continue
        ys=[y[c] for c in clade.clades]
        ax.plot([x[clade],x[clade]],[min(ys),max(ys)],color="black",lw=0.6)
        for child in clade.clades:
            ax.plot([x[clade],x[child]],[y[child],y[child]],color="black",lw=0.6)
    ax.invert_yaxis()
    return x,y,terms

def save_both(fig, stem):
    fig.savefig(stem.with_suffix(".png"),dpi=600,bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"),bbox_inches="tight")

def core_tree_figure(root,figdir):
    tree=Phylo.read(root/"preprint"/"results"/"auto"/"core_phylogeny"/"UBCG_core_tree.treefile","newick")
    meta=pd.read_csv(root/"preprint"/"results"/"auto"/"core_tree_metadata.tsv",sep="\t",dtype=str).fillna("")
    meta=meta.set_index("assembly")

    fig,ax=plt.subplots(figsize=(12,18))
    x,y,terms=draw_rect_tree(ax,tree,cladogram=False)
    xmax=max(x.values())
    dx=xmax*0.055
    xs=[xmax+dx,xmax+2*dx,xmax+3*dx]

    # categorical operon colors
    ops=sorted({meta.loc[t.name,"operon_type"] for t in terms if t.name in meta.index and meta.loc[t.name,"operon_type"]})
    cmap=plt.get_cmap("tab20")
    opcolors={op:cmap(i%20) for i,op in enumerate(ops)}

    for t in terms:
        if t.name not in meta.index:
            continue
        r=meta.loc[t.name]
        yy=y[t]
        m=r.get("SVMSY_region_aligned","")
        if m:
            ax.scatter(xs[0],yy,s=18,marker="s",color=MOTIF_COLORS.get(m,"0.6"),edgecolors="none")
        p=r.get("day7_category","")
        if p in PHENO_COLORS:
            ax.scatter(xs[1],yy,s=18,marker="o",color=PHENO_COLORS[p],edgecolors="none")
        op=r.get("operon_type","")
        if op:
            ax.scatter(xs[2],yy,s=20,marker="^",color=opcolors.get(op,"0.6"),edgecolors="none")

    top=-2
    for xx,label in zip(xs,["Marker","Day-7 phenotype","Operon type"]):
        ax.text(xx,top,label,rotation=45,ha="left",va="bottom",fontsize=9)

    from matplotlib.lines import Line2D
    motif_handles=[Line2D([0],[0],marker="s",linestyle="",color=MOTIF_COLORS[k],label=k) for k in MOTIF_COLORS if k in set(meta["SVMSY_region_aligned"])]
    ph_handles=[Line2D([0],[0],marker="o",linestyle="",color=PHENO_COLORS[k],label=k) for k in PHENO_COLORS]
    op_handles=[Line2D([0],[0],marker="^",linestyle="",color=opcolors[k],label=f"Type {k}") for k in ops]

    leg1=ax.legend(handles=motif_handles,title="Five-residue state",loc="upper left",bbox_to_anchor=(1.01,1.0),frameon=False)
    ax.add_artist(leg1)
    leg2=ax.legend(handles=ph_handles,title="Day-7 phenotype",loc="upper left",bbox_to_anchor=(1.01,0.72),frameon=False)
    ax.add_artist(leg2)
    ax.legend(handles=op_handles,title="Published operon type",loc="upper left",bbox_to_anchor=(1.01,0.42),frameon=False,ncol=2)

    ax.set_title("Core-genome phylogeny of the 87-genome pilot cohort")
    ax.set_xlabel("Branch length")
    ax.set_yticks([])
    ax.spines[["top","right","left"]].set_visible(False)
    ax.set_xlim(0,xs[-1]+dx)
    fig.tight_layout()
    save_both(fig,figdir/"Fig1_core_tree_integrated")
    plt.close(fig)

def external_tree_figure(root,figdir):
    tree=Phylo.read(root/"preprint"/"external_validation"/"external_serralysin_tree.treefile","newick")
    meta=pd.read_csv(root/"preprint"/"external_validation"/"external_integrated_serralysin_table.tsv",sep="\t",dtype=str).fillna("")
    meta=meta[meta["status"]=="candidate_serralysin"].copy()
    meta["tip_name"]=meta["accession"]+"|"+meta["protein_id"]+"|"+meta["organism"].str.replace(" ","_",regex=False)
    meta=meta.set_index("tip_name")

    fig,ax=plt.subplots(figsize=(10,18))
    x,y,terms=draw_rect_tree(ax,tree,cladogram=True)
    xmax=max(x.values())
    xs=[xmax+1.0,xmax+2.0]
    class_colors={
        "AprX-reference-proximal":"#4c78a8",
        "P.aeruginosa-AprA-reference-proximal":"#e45756",
    }

    ref_tips=[]
    for t in terms:
        yy=y[t]
        if t.name in meta.index:
            r=meta.loc[t.name]
            c=r.get("nearest_reference_class","")
            m=r.get("SVMSY_aligned_pattern","")
            if c:
                ax.scatter(xs[0],yy,s=13,marker="s",color=class_colors.get(c,"0.6"),edgecolors="none")
            if m:
                ax.scatter(xs[1],yy,s=13,marker="s",color=MOTIF_COLORS.get(m,"0.6"),edgecolors="none")
        elif t.name.startswith("REF_"):
            ref_tips.append(t)
            ax.scatter(x[t],yy,s=30,marker="*",color="black",zorder=5)
            ax.text(xmax+0.2,yy,t.name.split("|")[0],fontsize=7,va="center",ha="left")

    ax.text(xs[0],-3,"Reference proximity",rotation=45,ha="left",va="bottom",fontsize=9)
    ax.text(xs[1],-3,"Marker state",rotation=45,ha="left",va="bottom",fontsize=9)

    from matplotlib.lines import Line2D
    cls_handles=[Line2D([0],[0],marker="s",linestyle="",color=v,label=k.replace("P.aeruginosa-","")) for k,v in class_colors.items()]
    motif_handles=[Line2D([0],[0],marker="s",linestyle="",color=MOTIF_COLORS[k],label=k) for k in MOTIF_COLORS if k in set(meta["SVMSY_aligned_pattern"])]
    ref_handle=[Line2D([0],[0],marker="*",linestyle="",color="black",label="Curated reference")]

    l1=ax.legend(handles=cls_handles,title="Nearest-reference class",loc="upper left",bbox_to_anchor=(1.01,1.0),frameon=False)
    ax.add_artist(l1)
    l2=ax.legend(handles=motif_handles,title="Five-residue state",loc="upper left",bbox_to_anchor=(1.01,0.78),frameon=False)
    ax.add_artist(l2)
    ax.legend(handles=ref_handle,loc="upper left",bbox_to_anchor=(1.01,0.58),frameon=False)

    ax.set_title("Reference-aware serralysin phylogeny across 231 Pseudomonas type genomes")
    ax.set_xlabel("Cladogram depth")
    ax.set_yticks([])
    ax.spines[["top","right","left"]].set_visible(False)
    ax.set_xlim(0,xs[-1]+1.0)
    fig.tight_layout()
    save_both(fig,figdir/"Fig3_external_serralysin_cladogram")
    plt.close(fig)

def marker_bar(root,figdir):
    df=pd.read_csv(root/"preprint"/"external_validation"/"external_integrated_serralysin_table.tsv",sep="\t",dtype=str).fillna("")
    df=df[df["status"]=="candidate_serralysin"].copy()
    order=["AprX-reference-proximal","P.aeruginosa-AprA-reference-proximal"]
    patterns=["SVMSY","SLMSY","SIMSY","TVMSY"]
    counts=pd.crosstab(df["nearest_reference_class"],df["SVMSY_aligned_pattern"]).reindex(index=order,fill_value=0)
    counts=counts.reindex(columns=patterns,fill_value=0)
    frac=counts.div(counts.sum(axis=1),axis=0)

    fig,ax=plt.subplots(figsize=(8,4))
    left=np.zeros(len(order))
    for p in patterns:
        vals=frac[p].values
        bars=ax.barh(range(len(order)),vals,left=left,label=p,color=MOTIF_COLORS[p])
        for i,(v,l) in enumerate(zip(vals,left)):
            if v>=0.025:
                ax.text(l+v/2,i,f"{100*v:.1f}%",ha="center",va="center",fontsize=9)
        left+=vals
    labels=[
        f"AprX-proximal (n={counts.loc[order[0]].sum()})",
        f"AprA-proximal (n={counts.loc[order[1]].sum()})"
    ]
    ax.set_yticks(range(2),labels)
    ax.set_xlim(0,1)
    ax.set_xticks(np.linspace(0,1,6),[f"{int(x*100)}%" for x in np.linspace(0,1,6)])
    ax.set_xlabel("Candidate serralysins")
    ax.set_title("Five-residue marker pattern by reference-proximity class")
    ax.legend(title="Marker pattern",bbox_to_anchor=(1.02,1),loc="upper left",frameon=False)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    save_both(fig,figdir/"Fig4A_external_marker_distribution")
    plt.close(fig)

def homology_context(root,figdir,cohort):
    opdir=root/"preprint"/"results"/"auto"/"operon_homology"
    arch=pd.read_csv(opdir/"operon_architectures_final.tsv",sep="\t",dtype=str).fillna("")
    g=arch[arch["cohort"]==cohort].copy()
    genes=["aprI","aprD","aprE","aprF","prtA_like","prtB_like","lipA1","lipA2"]
    labels=["AprI","AprD","AprE","AprF","PrtA-like","PrtB-like","LipA1","LipA2"]
    vals=[]
    nums=[]
    for gene in genes:
        n=(g[gene+"_within50kb"]=="YES").sum()
        vals.append(100*n/len(g))
        nums.append(n)
    fig,ax=plt.subplots(figsize=(7,4.8))
    y=np.arange(len(labels))
    ax.barh(y,vals)
    ax.set_yticks(y,labels)
    ax.invert_yaxis()
    ax.set_xlim(0,108)
    ax.set_xlabel("Eligible serralysin genomes with homologue within 50 kb (%)")
    title="Pilot cohort" if cohort=="pilot" else "External type-genome cohort"
    ax.set_title(f"Homology-based serralysin locus context: {title}")
    for i,(v,n) in enumerate(zip(vals,nums)):
        ax.text(v+1,i,f"{n}/{len(g)} ({v:.1f}%)",va="center",fontsize=8)
    ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout()
    save_both(fig,figdir/f"Fig4B_{cohort}_homology_context")
    plt.close(fig)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()
    root=Path(args.project_root)
    figdir=root/"preprint"/"figures"/"publication"
    figdir.mkdir(parents=True,exist_ok=True)

    core_tree_figure(root,figdir)
    external_tree_figure(root,figdir)
    marker_bar(root,figdir)
    homology_context(root,figdir,"pilot")
    homology_context(root,figdir,"external")

    print("Publication figures written to:",figdir)
    for p in sorted(figdir.glob("*")):
        print(p)

if __name__=="__main__":
    main()
