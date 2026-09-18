#!/usr/bin/env python3

from pathlib import Path
import csv
import math

VALIDATED_CALLS = Path(
    "preprint/external_validation/metturn_audit/"
    "svmsy_corrected_calls.tsv"
)

INTEGRATED = Path(
    "preprint/external_validation/"
    "external_integrated_serralysin_table.tsv"
)

OUTDIR = Path(
    "preprint/external_validation/metturn_audit"
)
OUTDIR.mkdir(parents=True, exist_ok=True)

VALIDATED_TABLE = OUTDIR / (
    "external_integrated_serralysin_table_validated.tsv"
)

STATS_OUT = OUTDIR / (
    "validated_svmsy_statistics.tsv"
)

REPORT_OUT = OUTDIR / (
    "validated_svmsy_statistics_report.txt"
)


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(
            csv.DictReader(f, delimiter="\t")
        )


def write_tsv(path, rows, fieldnames):
    with open(
        path,
        "w",
        newline="",
        encoding="utf-8"
    ) as f:
        w = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


def clean(x):
    if x is None:
        return ""
    return str(x).strip()


def boolify(x):
    return clean(x).lower() in {
        "true",
        "yes",
        "1"
    }


def log_comb(n, k):
    if k < 0 or k > n:
        return float("-inf")

    return (
        math.lgamma(n + 1)
        - math.lgamma(k + 1)
        - math.lgamma(n - k + 1)
    )


def hypergeom_prob(x, row1, col1, n):
    """
    Probability for cell a=x given fixed margins.
    """

    col2 = n - col1

    if x < 0:
        return 0.0

    if row1 - x < 0:
        return 0.0

    if x > col1:
        return 0.0

    if row1 - x > col2:
        return 0.0

    logp = (
        log_comb(col1, x)
        + log_comb(col2, row1 - x)
        - log_comb(n, row1)
    )

    return math.exp(logp)


def fisher_exact_two_sided(a, b, c, d):
    """
    Two-sided Fisher exact test using the same general
    probability-ordering definition used by SciPy.

    Table:
        [[a, b],
         [c, d]]
    """

    row1 = a + b
    row2 = c + d
    col1 = a + c
    col2 = b + d
    n = row1 + row2

    xmin = max(0, row1 - col2)
    xmax = min(row1, col1)

    p_obs = hypergeom_prob(
        a,
        row1,
        col1,
        n
    )

    tol = 1e-12

    p_two = 0.0

    for x in range(xmin, xmax + 1):
        p = hypergeom_prob(
            x,
            row1,
            col1,
            n
        )

        if p <= p_obs + tol:
            p_two += p

    p_two = min(1.0, p_two)

    if b == 0 or c == 0:
        if a * d > 0:
            odds = float("inf")
        elif a * d == 0:
            odds = float("nan")
        else:
            odds = 0.0
    else:
        odds = (a * d) / (b * c)

    return odds, p_two


def odds_ratio_ci(a, b, c, d):
    """
    Approximate Wald CI on log odds ratio.

    Adds 0.5 only if a zero cell is present.
    """

    cells = [a, b, c, d]

    corrected = any(x == 0 for x in cells)

    if corrected:
        a2 = a + 0.5
        b2 = b + 0.5
        c2 = c + 0.5
        d2 = d + 0.5
    else:
        a2, b2, c2, d2 = (
            float(a),
            float(b),
            float(c),
            float(d)
        )

    odds = (
        a2 * d2
    ) / (
        b2 * c2
    )

    se = math.sqrt(
        1/a2
        + 1/b2
        + 1/c2
        + 1/d2
    )

    z = 1.959963984540054

    lower = math.exp(
        math.log(odds) - z * se
    )

    upper = math.exp(
        math.log(odds) + z * se
    )

    return odds, lower, upper, corrected


print("=" * 70)
print("VALIDATED SVMSY STATISTICAL REANALYSIS")
print("=" * 70)

if not VALIDATED_CALLS.exists():
    raise SystemExit(
        f"Missing: {VALIDATED_CALLS}"
    )

if not INTEGRATED.exists():
    raise SystemExit(
        f"Missing: {INTEGRATED}"
    )


calls = read_tsv(VALIDATED_CALLS)
integrated = read_tsv(INTEGRATED)

print()
print(
    "Validated call rows:",
    len(calls)
)

print(
    "Integrated table rows:",
    len(integrated)
)


if not calls:
    raise SystemExit(
        "Validated call table is empty."
    )

if not integrated:
    raise SystemExit(
        "Integrated table is empty."
    )


integrated_fields = list(
    integrated[0].keys()
)

call_fields = list(
    calls[0].keys()
)


assembly_col = None

for candidate in [
    "accession",
    "assembly"
]:
    if candidate in integrated_fields:
        assembly_col = candidate
        break

if assembly_col is None:
    raise SystemExit(
        "Could not identify accession/assembly "
        "column in integrated table."
    )

if "protein_id" not in integrated_fields:
    raise SystemExit(
        "protein_id column not found in "
        "integrated table."
    )


required_call_cols = [
    "assembly",
    "protein",
    "local_state5",
    "local_state6",
    "corrected_exact_SVMSY",
    "positional_status"
]

for c in required_call_cols:
    if c not in call_fields:
        raise SystemExit(
            f"Validated table missing column: {c}"
        )


# ---------------------------------------------------------
# Build validated lookup
# ---------------------------------------------------------

lookup = {}

for r in calls:

    key = (
        clean(r["assembly"]),
        clean(r["protein"])
    )

    if key in lookup:
        raise SystemExit(
            f"Duplicate validation key: {key}"
        )

    lookup[key] = r


# ---------------------------------------------------------
# Join validation onto integrated table
# ---------------------------------------------------------

merged = []

matched = 0
unmatched = []

for row in integrated:

    new = dict(row)

    key = (
        clean(row.get(assembly_col)),
        clean(row.get("protein_id"))
    )

    call = lookup.get(key)

    if call is not None:
        matched += 1

        new["SVMSY_original_alignment_pattern"] = clean(
            row.get("SVMSY_aligned_pattern")
        )

        new["SVMSY_validated_pattern"] = clean(
            call.get("local_state5")
        )

        new[
            "SVMSY_validated_six_residue_state"
        ] = clean(
            call.get("local_state6")
        )

        new["SVMSY_validated_exact"] = (
            "YES"
            if boolify(
                call.get(
                    "corrected_exact_SVMSY"
                )
            )
            else "NO"
        )

        new[
            "SVMSY_positional_validation"
        ] = clean(
            call.get("positional_status")
        )

        new[
            "SVMSY_assignment_method"
        ] = (
            "independent_zinc_anchored_"
            "catalytic_region_mapping"
        )

        new[
            "SVMSY_zinc_motif"
        ] = clean(
            call.get("zinc_motif")
        )

        new[
            "SVMSY_zinc_offset"
        ] = clean(
            call.get(
                "local_offset_from_zinc"
            )
        )

        new[
            "SVMSY_catalytic_identity_anchor"
        ] = clean(
            call.get(
                "catalytic_region_identity_to_anchor"
            )
        )

    else:
        unmatched.append(key)

        new[
            "SVMSY_original_alignment_pattern"
        ] = clean(
            row.get("SVMSY_aligned_pattern")
        )

        new["SVMSY_validated_pattern"] = ""
        new[
            "SVMSY_validated_six_residue_state"
        ] = ""
        new["SVMSY_validated_exact"] = ""
        new[
            "SVMSY_positional_validation"
        ] = ""
        new["SVMSY_assignment_method"] = ""
        new["SVMSY_zinc_motif"] = ""
        new["SVMSY_zinc_offset"] = ""
        new[
            "SVMSY_catalytic_identity_anchor"
        ] = ""

    merged.append(new)


print()
print(
    "Matched to validated calls:",
    matched
)

print(
    "Integrated rows unmatched:",
    len(unmatched)
)


if matched != 231:
    print()
    print(
        "WARNING: expected 231 validated "
        "candidate matches."
    )


# ---------------------------------------------------------
# Compare old versus validated state
# ---------------------------------------------------------

validated_rows = [
    r for r in merged
    if clean(
        r.get("SVMSY_validated_pattern")
    )
]

disagreements = []

for r in validated_rows:

    old = clean(
        r.get(
            "SVMSY_original_alignment_pattern"
        )
    )

    new = clean(
        r.get(
            "SVMSY_validated_pattern"
        )
    )

    if old != new:
        disagreements.append(r)


print()
print("STATE ASSIGNMENT AUDIT")

print(
    "Compared:",
    len(validated_rows)
)

print(
    "Exact agreements:",
    len(validated_rows)
    - len(disagreements)
)

print(
    "Disagreements:",
    len(disagreements)
)


if disagreements:
    print()

    for r in disagreements:
        print(
            clean(r.get(assembly_col)),
            clean(r.get("protein_id")),
            clean(
                r.get(
                    "SVMSY_original_alignment_pattern"
                )
            ),
            "->",
            clean(
                r.get(
                    "SVMSY_validated_pattern"
                )
            )
        )


# ---------------------------------------------------------
# Write validated integrated table
# ---------------------------------------------------------

extra_fields = [
    "SVMSY_original_alignment_pattern",
    "SVMSY_validated_pattern",
    "SVMSY_validated_six_residue_state",
    "SVMSY_validated_exact",
    "SVMSY_positional_validation",
    "SVMSY_assignment_method",
    "SVMSY_zinc_motif",
    "SVMSY_zinc_offset",
    "SVMSY_catalytic_identity_anchor",
]

out_fields = list(integrated_fields)

for x in extra_fields:
    if x not in out_fields:
        out_fields.append(x)

write_tsv(
    VALIDATED_TABLE,
    merged,
    out_fields
)


# ---------------------------------------------------------
# Find reference proximity classes
# ---------------------------------------------------------

class_col = "nearest_reference_class"

if class_col not in integrated_fields:
    raise SystemExit(
        "nearest_reference_class not found."
    )


classes = sorted(
    {
        clean(r.get(class_col))
        for r in validated_rows
        if clean(r.get(class_col))
    }
)


print()
print("REFERENCE CLASSES")

for c in classes:
    print(c)


aprx_classes = [
    x for x in classes
    if (
        "aprx" in x.lower()
        and "reference" in x.lower()
        and "proximal" in x.lower()
    )
]

apra_classes = [
    x for x in classes
    if (
        "apra" in x.lower()
        and "reference" in x.lower()
        and "proximal" in x.lower()
    )
]


if len(aprx_classes) != 1:
    raise SystemExit(
        "Could not uniquely identify "
        "AprX-reference-proximal class: "
        + repr(aprx_classes)
    )

if len(apra_classes) != 1:
    raise SystemExit(
        "Could not uniquely identify "
        "AprA-reference-proximal class: "
        + repr(apra_classes)
    )


aprx_class = aprx_classes[0]
apra_class = apra_classes[0]


# ---------------------------------------------------------
# Validated contingency table
# ---------------------------------------------------------

apra_rows = [
    r for r in validated_rows
    if clean(r.get(class_col))
    == apra_class
]

aprx_rows = [
    r for r in validated_rows
    if clean(r.get(class_col))
    == aprx_class
]


A_sv = sum(
    clean(
        r.get(
            "SVMSY_validated_pattern"
        )
    ) == "SVMSY"
    for r in apra_rows
)

A_non = len(apra_rows) - A_sv

X_sv = sum(
    clean(
        r.get(
            "SVMSY_validated_pattern"
        )
    ) == "SVMSY"
    for r in aprx_rows
)

X_non = len(aprx_rows) - X_sv


validated_table = [
    [A_sv, A_non],
    [X_sv, X_non]
]


validated_or, validated_p = (
    fisher_exact_two_sided(
        A_sv,
        A_non,
        X_sv,
        X_non
    )
)


ci_or, ci_low, ci_high, corrected = (
    odds_ratio_ci(
        A_sv,
        A_non,
        X_sv,
        X_non
    )
)


if (
    not math.isnan(validated_or)
    and not math.isinf(validated_or)
    and validated_or != 0
):
    reciprocal_or = 1 / validated_or
else:
    reciprocal_or = float("nan")


# ---------------------------------------------------------
# Original contingency table on same rows
# ---------------------------------------------------------

OA_sv = sum(
    clean(
        r.get(
            "SVMSY_original_alignment_pattern"
        )
    ) == "SVMSY"
    for r in apra_rows
)

OA_non = len(apra_rows) - OA_sv

OX_sv = sum(
    clean(
        r.get(
            "SVMSY_original_alignment_pattern"
        )
    ) == "SVMSY"
    for r in aprx_rows
)

OX_non = len(aprx_rows) - OX_sv


original_table = [
    [OA_sv, OA_non],
    [OX_sv, OX_non]
]


original_or, original_p = (
    fisher_exact_two_sided(
        OA_sv,
        OA_non,
        OX_sv,
        OX_non
    )
)


print()
print("=" * 70)
print("VALIDATED FISHER TEST")
print("=" * 70)

print()
print("AprA-proximal class:")
print(apra_class)

print()
print("AprX-proximal class:")
print(aprx_class)

print()
print(
    "AprA-proximal:",
    "SVMSY =",
    A_sv,
    "; non-SVMSY =",
    A_non
)

print(
    "AprX-proximal:",
    "SVMSY =",
    X_sv,
    "; non-SVMSY =",
    X_non
)

print()
print(
    "Validated contingency table:",
    validated_table
)

print(
    "OR AprA vs AprX:",
    validated_or
)

print(
    "OR AprX vs AprA:",
    reciprocal_or
)

print(
    "Two-sided Fisher p:",
    validated_p
)

print(
    "Approximate 95% CI AprA vs AprX:",
    ci_low,
    "to",
    ci_high
)

print(
    "Zero-cell correction for CI:",
    corrected
)


print()
print("=" * 70)
print("ORIGINAL VERSUS VALIDATED")
print("=" * 70)

print(
    "Original table:",
    original_table
)

print(
    "Validated table:",
    validated_table
)

print(
    "Original OR:",
    original_or
)

print(
    "Validated OR:",
    validated_or
)

print(
    "Original p:",
    original_p
)

print(
    "Validated p:",
    validated_p
)

same_table = (
    original_table == validated_table
)

same_or = (
    original_or == validated_or
)

same_p = (
    abs(
        original_p
        - validated_p
    ) < 1e-15
)

print(
    "Contingency table unchanged:",
    same_table
)

print(
    "OR unchanged:",
    same_or
)

print(
    "p unchanged:",
    same_p
)


# ---------------------------------------------------------
# Save statistics
# ---------------------------------------------------------

stats_rows = [{
    "analysis":
        "validated_external_SVMSY_"
        "enrichment_AprAprox_vs_AprXprox",

    "assignment_method":
        "independent_zinc_anchored_"
        "catalytic_region_mapping",

    "AprA_class":
        apra_class,

    "AprX_class":
        aprx_class,

    "AprA_SVMSY":
        A_sv,

    "AprA_non_SVMSY":
        A_non,

    "AprX_SVMSY":
        X_sv,

    "AprX_non_SVMSY":
        X_non,

    "OR_AprA_vs_AprX":
        validated_or,

    "OR_AprX_vs_AprA":
        reciprocal_or,

    "OR_95CI_lower_AprA_vs_AprX":
        ci_low,

    "OR_95CI_upper_AprA_vs_AprX":
        ci_high,

    "Fisher_two_sided_p":
        validated_p,

    "original_OR":
        original_or,

    "original_p":
        original_p,

    "state_disagreements":
        len(disagreements),

    "contingency_table_unchanged":
        same_table,

    "OR_unchanged":
        same_or,

    "p_unchanged":
        same_p,
}]


write_tsv(
    STATS_OUT,
    stats_rows,
    list(stats_rows[0].keys())
)


# ---------------------------------------------------------
# Final report
# ---------------------------------------------------------

report = []

report.append(
    "VALIDATED SVMSY STATISTICAL REANALYSIS"
)

report.append("=" * 60)

report.append("")
report.append(
    f"Validated rows: {len(validated_rows)}"
)

report.append(
    f"State disagreements: "
    f"{len(disagreements)}"
)

report.append("")
report.append(
    f"AprA-proximal class: {apra_class}"
)

report.append(
    f"AprX-proximal class: {aprx_class}"
)

report.append("")
report.append(
    f"Validated table: {validated_table}"
)

report.append(
    f"Validated OR AprA vs AprX: "
    f"{validated_or}"
)

report.append(
    f"Validated OR AprX vs AprA: "
    f"{reciprocal_or}"
)

report.append(
    f"Validated Fisher p: "
    f"{validated_p}"
)

report.append(
    f"Approximate 95% CI: "
    f"{ci_low} to {ci_high}"
)

report.append("")
report.append(
    f"Original table: {original_table}"
)

report.append(
    f"Original OR: {original_or}"
)

report.append(
    f"Original p: {original_p}"
)

report.append("")
report.append(
    f"Contingency table unchanged: "
    f"{same_table}"
)

report.append(
    f"OR unchanged: {same_or}"
)

report.append(
    f"p unchanged: {same_p}"
)

report.append("")

if (
    len(disagreements) == 0
    and same_table
    and same_or
    and same_p
):
    report.append(
        "CONCLUSION: Independent zinc-anchored "
        "catalytic-region mapping reproduces the "
        "original SVMSY classification and statistical "
        "result exactly. The reported association is "
        "therefore not caused by transferring global "
        "alignment columns 258-262 from a single "
        "anchor sequence."
    )
else:
    report.append(
        "CONCLUSION: The independent positional "
        "validation changes at least part of the "
        "original statistical analysis. Use the "
        "validated result."
    )


REPORT_OUT.write_text(
    "\n".join(report) + "\n",
    encoding="utf-8"
)


print()
print("=" * 70)
print("FILES WRITTEN")
print("=" * 70)

print(VALIDATED_TABLE)
print(STATS_OUT)
print(REPORT_OUT)
