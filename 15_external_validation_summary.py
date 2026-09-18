import argparse
from pathlib import Path
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    ext = root/"preprint"/"external_validation"

    external = pd.read_csv(
        ext/"external_type_aprx_marker_table.tsv",
        sep="\t",
        dtype=str
    ).fillna("")

    pilot = pd.read_csv(
        root/"preprint"/"results"/"auto"/"AprX_analysis_dataset.tsv",
        sep="\t",
        dtype=str
    ).fillna("")

    cand = external[external["status"]=="candidate_serralysin"].copy()

    lines = [
        "# External taxonomic validation summary",
        "",
        "This external cohort consists of current annotated RefSeq Pseudomonas assemblies flagged by NCBI as type material, with one selected assembly per organism name.",
        "",
        "## External cohort",
        "",
        f"- Selected type-material genomes: {len(external)}",
        f"- Genomes with a convincing full-length serralysin candidate under the current sequence criteria: {len(cand)}",
        f"- Genomes without a convincing candidate under the current sequence criteria: {len(external)-len(cand)}",
        "",
        "## Marker-region diversity in external candidates",
        ""
    ]

    if not cand.empty:
        vc = cand["SVMSY_aligned_pattern"].replace("", "unresolved").value_counts()
        for pattern, n in vc.items():
            lines.append(f"- {pattern}: {n}")

        exact = int((cand["SVMSY_aligned_pattern"]=="SVMSY").sum())
        lines.append("")
        lines.append(
            f"- Exact SVMSY among external serralysin candidates: {exact}/{len(cand)} "
            f"({100*exact/len(cand):.1f}%)"
        )

    lines += [
        "",
        "## Pilot cohort comparison",
        ""
    ]

    p = pilot[pilot["classification"]=="AprX-reference-like"].copy()

    if "SVMSY_region_aligned" in p.columns:
        pvc = p["SVMSY_region_aligned"].replace("", "unresolved").value_counts()
        for pattern, n in pvc.items():
            lines.append(f"- Pilot {pattern}: {n}")

    lines += [
        "",
        "Interpretation: this phase tests taxonomic prevalence and sequence diversity of the proposed marker region. It does not independently validate spoilage phenotype because the external type-genome cohort lacks harmonized proteolytic measurements."
    ]

    out = ext/"external_validation_summary.md"
    out.write_text("\n".join(lines)+"\n")

    print("\n".join(lines))
    print("\nSaved:", out)


if __name__ == "__main__":
    main()
