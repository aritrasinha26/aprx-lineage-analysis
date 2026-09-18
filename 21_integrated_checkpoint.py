import argparse
from pathlib import Path

import pandas as pd


def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project-root",required=True)
    args=ap.parse_args()

    root=Path(args.project_root)
    res=root/"preprint"/"results"/"auto"
    ext=root/"preprint"/"external_validation"

    pilot=pd.read_csv(res/"AprX_analysis_dataset_with_operon.tsv",sep="\t",dtype=str).fillna("")
    if pilot.empty:
        pilot=pd.read_csv(res/"AprX_analysis_dataset.tsv",sep="\t",dtype=str).fillna("")

    extdf=pd.read_csv(ext/"external_integrated_serralysin_table.tsv",sep="\t",dtype=str).fillna("")
    extcand=extdf[extdf["status"]=="candidate_serralysin"].copy()

    phase4=(res/"phase4_phylo_diagnostic_summary.md").read_text() if (res/"phase4_phylo_diagnostic_summary.md").exists() else ""

    stats=pd.read_csv(res/"phase7_external_statistics.tsv",sep="\t")
    row=stats[stats["analysis"]=="external_SVMSY_enrichment_AprAprox_vs_AprXprox"].iloc[0]

    aprx=extcand[extcand["nearest_reference_class"]=="AprX-reference-proximal"]
    apra=extcand[extcand["nearest_reference_class"]=="P.aeruginosa-AprA-reference-proximal"]

    svx=(aprx["SVMSY_aligned_pattern"]=="SVMSY").sum()
    sva=(apra["SVMSY_aligned_pattern"]=="SVMSY").sum()

    lines=[]
    lines.append("# Integrated analysis checkpoint")
    lines.append("")
    lines.append("## Pilot phenotype cohort")
    lines.append("")
    lines.append(f"- Genomes: {len(pilot)}")
    lines.append(f"- AprX-reference-like: {(pilot['classification']=='AprX-reference-like').sum()}")
    lines.append(f"- Day-7 phenotype available: {pilot['day7_category'].isin(['non','weak','moderate','strong','very_strong']).sum()}")
    lines.append("- Within-species-variable AprX site screen: 11 testable sites; 0 FDR-significant.")
    lines.append("- The five-residue motif state showed no within-species variation in the pilot cohort.")
    lines.append("- Phylogenetically regularized analysis did not support an independent SVMSY–phenotype association.")
    lines.append("- Operon type 1 versus type 8 retained a large phenotype difference in the lambda-PGLS sensitivity analysis, but the trait had few independent transitions.")
    lines.append("")
    lines.append("## External type-genome cohort")
    lines.append("")
    lines.append(f"- Type-material genomes screened: {len(extdf)}")
    lines.append(f"- Full-length serralysin candidates: {len(extcand)}")
    lines.append(f"- AprX-reference-proximal: {len(aprx)}")
    lines.append(f"- P. aeruginosa AprA-reference-proximal: {len(apra)}")
    lines.append(f"- Exact SVMSY in AprX-reference-proximal group: {svx}/{len(aprx)} ({100*svx/len(aprx):.1f}%)")
    lines.append(f"- Exact SVMSY in AprA-reference-proximal group: {sva}/{len(apra)} ({100*sva/len(apra):.1f}%)")
    lines.append(f"- Fisher exact comparison, AprA-proximal versus AprX-proximal: OR={row['effect']:.3g}, p={row['p_value']:.3g}.")
    lines.append("- Therefore, SVMSY is not specific to the AprX-reference-proximal group in the broad external serralysin dataset.")
    lines.append("")
    lines.append("## Current biological interpretation")
    lines.append("")
    lines.append(
        "The combined analysis supports a model in which AprX/serralysin sequence states and proteolytic phenotype are strongly structured by bacterial lineage. "
        "The five-residue SVMSY motif is common but neither universal nor AprX-specific in the external reference-aware dataset, while the pilot phenotype cohort does not provide evidence that it independently predicts higher proteolysis after phylogenetic structure is considered. "
        "By contrast, broader genomic organization, particularly operon architecture, is more strongly associated with phenotype in the pilot cohort."
    )
    lines.append("")
    lines.append("This conclusion concerns the datasets analysed here. It does not by itself exclude a lineage-specific predictive use of SVMSY in another defined strain population.")

    out=res/"integrated_analysis_checkpoint.md"
    out.write_text("\n".join(lines)+"\n")
    print("\n".join(lines))
    print("\nSaved:",out)


if __name__=="__main__":
    main()
