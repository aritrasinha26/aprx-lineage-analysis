#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import re
import math

ROOT = Path("preprint")
EV = ROOT / "external_validation"

CANDIDATES = EV / "external_type_serralysin_candidates.faa"
GLOBAL_ALIGNMENT = EV / "external_reference_alignment.faa"
REFERENCE_PANEL = EV / "external_reference_panel.faa"

OUTDIR = EV / "metturn_audit"
OUTDIR.mkdir(parents=True, exist_ok=True)

ANCHOR_ID = "WP_170037105.1"

# Use the core M10 zinc-binding motif for POSITIONAL anchoring.
# Do not require the following Pro because that would incorrectly
# reject genuine sequence variants.
ZINC_RE = re.compile(r"HE..H..G..H")

ANCHOR_STATE5 = "SVMSY"
ANCHOR_STATE6 = "SVMSYW"

# Catalytic-region extraction relative to independently detected zinc motif.
UPSTREAM = 25
DOWNSTREAM = 120


def read_fasta(path):
    records = []
    header = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)

        if header is not None:
            records.append((header, "".join(seq)))

    return records


def write_fasta(records, path, width=80):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + "\n")


def parse_header(header):
    parts = header.split("|")

    assembly = ""
    protein = ""
    species = ""

    for part in parts:
        part = part.strip()

        if (
            part.startswith(("GCF_", "GCA_"))
            and not assembly
        ):
            assembly = part
            continue

        # RefSeq/GenBank protein accessions.
        # Explicitly exclude genome assembly accessions.
        if (
            not protein
            and not part.startswith(("GCF_", "GCA_"))
            and re.match(
                r"^[A-Z]{1,4}_\d+\.\d+$",
                part
            )
        ):
            protein = part

    if len(parts) >= 3:
        species = parts[-1].strip()

    return assembly, protein, species


def aln_col_to_resnum(aln):
    result = []
    n = 0

    for aa in aln:
        if aa == "-":
            result.append(None)
        else:
            n += 1
            result.append(n)

    return result


def needleman_wunsch(a, b, match=2, mismatch=-1, gap=-2):
    """
    Simple global protein alignment.

    We use it only on short catalytic regions that have already been
    independently anchored by the zinc-binding motif. This therefore
    provides an independent check on the original whole-protein MSA.
    """

    n = len(a)
    m = len(b)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = "L"

    for i in range(1, n + 1):
        ai = a[i - 1]

        for j in range(1, m + 1):
            bj = b[j - 1]

            diag = score[i - 1][j - 1] + (
                match if ai == bj else mismatch
            )
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap

            best = max(diag, up, left)
            score[i][j] = best

            # Prefer residue-residue correspondence when tied.
            if best == diag:
                trace[i][j] = "D"
            elif best == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    aa = []
    bb = []

    i = n
    j = m

    while i > 0 or j > 0:
        t = trace[i][j]

        if i > 0 and j > 0 and t == "D":
            aa.append(a[i - 1])
            bb.append(b[j - 1])
            i -= 1
            j -= 1

        elif i > 0 and (j == 0 or t == "U"):
            aa.append(a[i - 1])
            bb.append("-")
            i -= 1

        else:
            aa.append("-")
            bb.append(b[j - 1])
            j -= 1

    return "".join(reversed(aa)), "".join(reversed(bb))


def anchor_to_candidate_map(aligned_anchor, aligned_candidate):
    """
    Map each ungapped anchor-region residue index to the corresponding
    ungapped candidate-region residue index.
    """

    mapping = {}

    ia = 0
    ib = 0

    for a, b in zip(aligned_anchor, aligned_candidate):

        anchor_index = None
        candidate_index = None

        if a != "-":
            anchor_index = ia

        if b != "-":
            candidate_index = ib

        if anchor_index is not None:
            mapping[anchor_index] = candidate_index

        if a != "-":
            ia += 1

        if b != "-":
            ib += 1

    return mapping


def alignment_identity(a, b):
    matches = 0
    compared = 0

    for x, y in zip(a, b):
        if x == "-" or y == "-":
            continue

        compared += 1

        if x == y:
            matches += 1

    if compared == 0:
        return None

    return matches / compared


def contiguous_positions(pos):
    if not pos:
        return False

    if any(x is None for x in pos):
        return False

    return pos == list(range(pos[0], pos[0] + len(pos)))


def candidate_global_match(header, raw_seq, aligned_records, seq_index):
    """
    Prefer a header-specific match. Fall back to exact ungapped sequence.
    """

    assembly, protein, species = parse_header(header)

    direct = []

    for ah, aln in aligned_records:
        if assembly and assembly not in ah:
            continue

        if protein and protein not in ah:
            continue

        if aln.replace("-", "") == raw_seq:
            direct.append((ah, aln))

    if len(direct) == 1:
        return direct[0], "header_and_sequence"

    matches = seq_index.get(raw_seq, [])

    if matches:
        return matches[0], "sequence_only"

    return None, "not_found"


def extract_region(seq, zinc_start0):
    start = max(0, zinc_start0 - UPSTREAM)
    end = min(len(seq), zinc_start0 + DOWNSTREAM)

    return start, end, seq[start:end]


def map_anchor_targets(
    anchor_target_full0,
    anchor_region_start,
    candidate_region_start,
    candidate_seq,
    mapping
):
    chars = []
    positions = []

    for full_anchor_pos0 in anchor_target_full0:

        anchor_region_index = (
            full_anchor_pos0 - anchor_region_start
        )

        candidate_region_index = mapping.get(
            anchor_region_index
        )

        if candidate_region_index is None:
            chars.append("-")
            positions.append(None)
            continue

        full_candidate_pos0 = (
            candidate_region_start
            + candidate_region_index
        )

        chars.append(
            candidate_seq[full_candidate_pos0]
        )

        positions.append(
            full_candidate_pos0 + 1
        )

    return "".join(chars), positions


def safe_join(values):
    return ",".join(
        "" if x is None else str(x)
        for x in values
    )


def scan_analysis_sources(root):
    keywords = [
        "SVMSY",
        "fisher_exact",
        "odds ratio",
        "odds_ratio",
        "oddsratio",
        "13.0",
        "proteolytic",
        "spoilage",
    ]

    suffixes = {
        ".py", ".R", ".r", ".md",
        ".txt", ".tsv", ".csv"
    }

    hits = []

    for p in root.rglob("*"):

        if not p.is_file():
            continue

        if p.suffix not in suffixes:
            continue

        try:
            if p.stat().st_size > 20_000_000:
                continue
        except OSError:
            continue

        try:
            text = p.read_text(errors="ignore")
        except Exception:
            continue

        lines = text.splitlines()

        for lineno, line in enumerate(lines, 1):
            low = line.lower()

            matched = [
                k for k in keywords
                if k.lower() in low
            ]

            if matched:
                hits.append({
                    "path": str(p),
                    "line": lineno,
                    "keywords": ",".join(matched),
                    "text": line[:500]
                })

    return hits


def scan_possible_phenotype_tables(root):
    phenotype_terms = [
        "phenotype",
        "spoil",
        "proteol",
        "activity",
        "casein",
        "milk",
        "class",
        "positive",
        "negative",
    ]

    id_terms = [
        "assembly",
        "accession",
        "protein",
        "genome",
        "strain",
        "isolate",
    ]

    results = []

    for p in root.rglob("*"):

        if not p.is_file():
            continue

        if p.suffix.lower() not in {".tsv", ".csv"}:
            continue

        try:
            if p.stat().st_size > 50_000_000:
                continue
        except OSError:
            continue

        delim = "\t" if p.suffix.lower() == ".tsv" else ","

        try:
            with open(p, newline="", errors="ignore") as f:
                reader = csv.reader(f, delimiter=delim)
                header = next(reader, [])
        except Exception:
            continue

        lower = [x.lower() for x in header]

        phenotype_cols = [
            header[i]
            for i, x in enumerate(lower)
            if any(t in x for t in phenotype_terms)
        ]

        id_cols = [
            header[i]
            for i, x in enumerate(lower)
            if any(t in x for t in id_terms)
        ]

        if phenotype_cols and id_cols:
            results.append({
                "path": str(p),
                "id_columns": ";".join(id_cols),
                "phenotype_like_columns":
                    ";".join(phenotype_cols),
                "all_columns":
                    ";".join(header)
            })

    return results


candidates = read_fasta(CANDIDATES)
global_aln = read_fasta(GLOBAL_ALIGNMENT)

if len(candidates) == 0:
    raise SystemExit("No candidate sequences found.")

# ------------------------------------------------------------
# Establish anchor
# ------------------------------------------------------------

anchor_hits = [
    (h, s)
    for h, s in global_aln
    if ANCHOR_ID in h
]

if len(anchor_hits) != 1:
    raise SystemExit(
        f"Expected exactly one anchor, found "
        f"{len(anchor_hits)}"
    )

anchor_header, anchor_aligned = anchor_hits[0]
anchor_seq = anchor_aligned.replace("-", "")

sv5_start0 = anchor_seq.find(ANCHOR_STATE5)
sv6_start0 = anchor_seq.find(ANCHOR_STATE6)

if sv5_start0 == -1:
    raise SystemExit(
        f"{ANCHOR_STATE5} not found in anchor."
    )

if sv6_start0 == -1:
    raise SystemExit(
        f"{ANCHOR_STATE6} not found in anchor."
    )

anchor_zinc_hits = list(ZINC_RE.finditer(anchor_seq))

if len(anchor_zinc_hits) != 1:
    raise SystemExit(
        "Anchor does not contain exactly one "
        "HEXXHXXGXXH motif."
    )

anchor_zinc = anchor_zinc_hits[0]
anchor_zinc_start0 = anchor_zinc.start()

anchor_region_start, anchor_region_end, anchor_region = (
    extract_region(
        anchor_seq,
        anchor_zinc_start0
    )
)

anchor_target5_full0 = list(
    range(sv5_start0, sv5_start0 + 5)
)

anchor_target6_full0 = list(
    range(sv6_start0, sv6_start0 + 6)
)

# ------------------------------------------------------------
# Recover original whole-protein MSA columns
# ------------------------------------------------------------

anchor_colmap = aln_col_to_resnum(anchor_aligned)

global_cols5 = []
global_cols6 = []

for full0 in anchor_target5_full0:
    residue1 = full0 + 1

    hits = [
        i
        for i, r in enumerate(anchor_colmap)
        if r == residue1
    ]

    if len(hits) != 1:
        raise SystemExit(
            f"Cannot uniquely map anchor residue "
            f"{residue1}."
        )

    global_cols5.append(hits[0])

for full0 in anchor_target6_full0:
    residue1 = full0 + 1

    hits = [
        i
        for i, r in enumerate(anchor_colmap)
        if r == residue1
    ]

    if len(hits) != 1:
        raise SystemExit(
            f"Cannot uniquely map anchor residue "
            f"{residue1}."
        )

    global_cols6.append(hits[0])

# ------------------------------------------------------------
# Index global alignment
# ------------------------------------------------------------

seq_index = defaultdict(list)

for h, aln in global_aln:
    seq_index[aln.replace("-", "")].append(
        (h, aln)
    )

rows = []
raw_by_header = {}

# ------------------------------------------------------------
# Main audit
# ------------------------------------------------------------

for raw_header, raw_seq in candidates:

    raw_by_header[raw_header] = raw_seq

    assembly, protein, species = parse_header(
        raw_header
    )

    global_match, global_match_method = (
        candidate_global_match(
            raw_header,
            raw_seq,
            global_aln,
            seq_index
        )
    )

    if global_match is None:
        global_header = ""
        global_seq = ""
        global_state5 = ""
        global_state6 = ""
        global_positions5 = []
        global_positions6 = []
        global_contiguous5 = False
        global_contiguous6 = False

    else:
        global_header, global_seq = global_match
        cmap = aln_col_to_resnum(global_seq)

        global_state5 = "".join(
            global_seq[c]
            for c in global_cols5
        )

        global_state6 = "".join(
            global_seq[c]
            for c in global_cols6
        )

        global_positions5 = [
            cmap[c]
            for c in global_cols5
        ]

        global_positions6 = [
            cmap[c]
            for c in global_cols6
        ]

        global_contiguous5 = contiguous_positions(
            global_positions5
        )

        global_contiguous6 = contiguous_positions(
            global_positions6
        )

    # Independent zinc-site identification.
    zinc_hits = list(ZINC_RE.finditer(raw_seq))
    zinc_count = len(zinc_hits)

    strict_following_pro = False
    zinc_seq = ""
    zinc_start1 = None
    zinc_third_h1 = None

    local_state5 = ""
    local_state6 = ""
    local_positions5 = []
    local_positions6 = []
    local_contiguous5 = False
    local_contiguous6 = False
    pairwise_identity = None

    if zinc_count == 1:

        z = zinc_hits[0]

        zinc_seq = z.group()
        zinc_start1 = z.start() + 1
        zinc_third_h1 = z.start() + 11

        if z.end() < len(raw_seq):
            strict_following_pro = (
                raw_seq[z.end()] == "P"
            )

        cand_region_start, cand_region_end, cand_region = (
            extract_region(
                raw_seq,
                z.start()
            )
        )

        aligned_anchor_region, aligned_cand_region = (
            needleman_wunsch(
                anchor_region,
                cand_region
            )
        )

        pairwise_identity = alignment_identity(
            aligned_anchor_region,
            aligned_cand_region
        )

        amap = anchor_to_candidate_map(
            aligned_anchor_region,
            aligned_cand_region
        )

        local_state5, local_positions5 = (
            map_anchor_targets(
                anchor_target5_full0,
                anchor_region_start,
                cand_region_start,
                raw_seq,
                amap
            )
        )

        local_state6, local_positions6 = (
            map_anchor_targets(
                anchor_target6_full0,
                anchor_region_start,
                cand_region_start,
                raw_seq,
                amap
            )
        )

        local_contiguous5 = contiguous_positions(
            local_positions5
        )

        local_contiguous6 = contiguous_positions(
            local_positions6
        )

    global_local_same5 = (
        global_state5 == local_state5
        and global_positions5 == local_positions5
        and local_contiguous5
    )

    global_local_same6 = (
        global_state6 == local_state6
        and global_positions6 == local_positions6
        and local_contiguous6
    )

    corrected_exact_svmsy = (
        local_contiguous5
        and local_state5 == "SVMSY"
    )

    original_exact_svmsy = (
        global_contiguous5
        and global_state5 == "SVMSY"
    )

    central_met_conserved = (
        len(local_state6) == 6
        and local_state6[2] == "M"
    )

    terminal_w_conserved = (
        len(local_state6) == 6
        and local_state6[5] == "W"
    )

    if (
        zinc_count == 1
        and local_contiguous5
        and global_local_same5
    ):
        if central_met_conserved:
            positional_status = "VALIDATED"
        else:
            positional_status = (
                "VALIDATED_POSITION_NONCANONICAL_METTURN"
            )

    elif (
        zinc_count == 1
        and local_contiguous5
        and not global_local_same5
    ):
        positional_status = (
            "REASSIGNED_BY_CATALYTIC_ALIGNMENT"
        )

    else:
        positional_status = "UNCERTAIN"

    local_offset = None

    if (
        zinc_start1 is not None
        and local_positions5
        and local_positions5[0] is not None
    ):
        local_offset = (
            local_positions5[0] - zinc_start1
        )

    global_offset = None

    if (
        zinc_start1 is not None
        and global_positions5
        and global_positions5[0] is not None
    ):
        global_offset = (
            global_positions5[0] - zinc_start1
        )

    rows.append({
        "candidate_header": raw_header,
        "assembly": assembly,
        "protein": protein,
        "species": species,
        "protein_length": len(raw_seq),

        "zinc_motif_count": zinc_count,
        "zinc_motif": zinc_seq,
        "zinc_start_residue": zinc_start1,
        "zinc_third_H_residue": zinc_third_h1,
        "proline_immediately_after_zinc_motif":
            strict_following_pro,

        "global_match_method": global_match_method,
        "global_alignment_header": global_header,

        "global_state5": global_state5,
        "global_state6": global_state6,
        "global_positions5":
            safe_join(global_positions5),
        "global_positions6":
            safe_join(global_positions6),
        "global_contiguous5":
            global_contiguous5,
        "global_contiguous6":
            global_contiguous6,
        "global_offset_from_zinc":
            global_offset,

        "local_state5": local_state5,
        "local_state6": local_state6,
        "local_positions5":
            safe_join(local_positions5),
        "local_positions6":
            safe_join(local_positions6),
        "local_contiguous5":
            local_contiguous5,
        "local_contiguous6":
            local_contiguous6,
        "local_offset_from_zinc":
            local_offset,

        "catalytic_region_identity_to_anchor":
            (
                round(pairwise_identity, 5)
                if pairwise_identity is not None
                else ""
            ),

        "global_local_same5":
            global_local_same5,
        "global_local_same6":
            global_local_same6,

        "central_met_conserved":
            central_met_conserved,
        "terminal_w_conserved":
            terminal_w_conserved,

        "original_exact_SVMSY":
            original_exact_svmsy,
        "corrected_exact_SVMSY":
            corrected_exact_svmsy,

        "state_changed":
            (
                original_exact_svmsy
                != corrected_exact_svmsy
            ),

        "positional_status":
            positional_status,
    })

# ------------------------------------------------------------
# Write authoritative corrected calls
# ------------------------------------------------------------

calls_path = (
    OUTDIR / "svmsy_corrected_calls.tsv"
)

with open(calls_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(rows[0].keys()),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

# ------------------------------------------------------------
# Write problem sequences separately
# ------------------------------------------------------------

flagged = [
    r for r in rows
    if r["positional_status"] != "VALIDATED"
]

flagged_path = (
    OUTDIR / "svmsy_flagged_sequences.tsv"
)

with open(flagged_path, "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(rows[0].keys()),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(flagged)

# ------------------------------------------------------------
# Build catalytic-region FASTA
# ------------------------------------------------------------

cat_records = [
    (
        "ANCHOR|"
        + anchor_header.replace(" ", "_"),
        anchor_region
    )
]

for i, (header, seq) in enumerate(candidates, 1):
    z = list(ZINC_RE.finditer(seq))

    if len(z) != 1:
        continue

    start, end, region = extract_region(
        seq,
        z[0].start()
    )

    cat_records.append(
        (
            f"C{i:04d}|{header.replace(' ', '_')}",
            region
        )
    )

write_fasta(
    cat_records,
    OUTDIR / "catalytic_regions.faa"
)

# ------------------------------------------------------------
# Prepare structural validation panel
# ------------------------------------------------------------

selected_headers = {}
eligible = [
    r for r in rows
    if r["positional_status"].startswith(
        "VALIDATED"
    )
]

# Include the most divergent representative of every
# observed six-residue state.
by_state6 = defaultdict(list)

for r in eligible:
    by_state6[r["local_state6"]].append(r)

for state, group in by_state6.items():

    def identity_value(x):
        v = x[
            "catalytic_region_identity_to_anchor"
        ]
        return float(v) if v != "" else 999

    chosen = min(group, key=identity_value)

    selected_headers[
        chosen["candidate_header"]
    ] = (
        f"representative_state_{state}"
    )

# Include overall most divergent candidate.
if eligible:
    divergent = min(
        eligible,
        key=lambda x: (
            float(
                x[
                    "catalytic_region_identity_to_anchor"
                ]
            )
            if x[
                "catalytic_region_identity_to_anchor"
            ] != ""
            else 999
        )
    )

    selected_headers[
        divergent["candidate_header"]
    ] = "lowest_catalytic_identity"

# Include smallest and largest zinc-to-state offsets.
offset_rows = [
    r for r in eligible
    if r["local_offset_from_zinc"] is not None
]

if offset_rows:

    min_offset = min(
        offset_rows,
        key=lambda x: int(
            x["local_offset_from_zinc"]
        )
    )

    max_offset = max(
        offset_rows,
        key=lambda x: int(
            x["local_offset_from_zinc"]
        )
    )

    selected_headers[
        min_offset["candidate_header"]
    ] = (
        selected_headers.get(
            min_offset["candidate_header"],
            ""
        )
        + ";minimum_offset"
    ).strip(";")

    selected_headers[
        max_offset["candidate_header"]
    ] = (
        selected_headers.get(
            max_offset["candidate_header"],
            ""
        )
        + ";maximum_offset"
    ).strip(";")

# Include every sequence that required reassignment
# or remains uncertain.
for r in flagged:
    selected_headers[
        r["candidate_header"]
    ] = (
        selected_headers.get(
            r["candidate_header"],
            ""
        )
        + ";flagged"
    ).strip(";")

structure_records = [
    (
        "ANCHOR|"
        + anchor_header.replace(" ", "_"),
        anchor_seq
    )
]

structure_meta = [{
    "header": anchor_header,
    "reason": "anchor_exact_SVMSY"
}]

for header, reason in selected_headers.items():
    if header in raw_by_header:
        structure_records.append(
            (
                header.replace(" ", "_"),
                raw_by_header[header]
            )
        )

        structure_meta.append({
            "header": header,
            "reason": reason
        })

# Add AprA references if they are present locally.
if REFERENCE_PANEL.exists():

    refs = read_fasta(REFERENCE_PANEL)

    apra_added = 0

    for h, seq in refs:
        upper = h.upper()

        if (
            "APRA" in upper
            or "PAO1" in upper
            or "PA14" in upper
        ):
            structure_records.append(
                (
                    "REFERENCE|"
                    + h.replace(" ", "_"),
                    seq
                )
            )

            structure_meta.append({
                "header": h,
                "reason": "AprA_reference"
            })

            apra_added += 1

            if apra_added >= 2:
                break

write_fasta(
    structure_records,
    OUTDIR / "structure_validation_panel.faa"
)

with open(
    OUTDIR / "structure_validation_panel.tsv",
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=["header", "reason"],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(structure_meta)

# ------------------------------------------------------------
# Search project for current phenotype / OR analysis
# ------------------------------------------------------------

source_hits = scan_analysis_sources(ROOT)

with open(
    OUTDIR / "analysis_source_hits.tsv",
    "w",
    newline=""
) as f:

    if source_hits:
        writer = csv.DictWriter(
            f,
            fieldnames=list(
                source_hits[0].keys()
            ),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(source_hits)

possible_tables = scan_possible_phenotype_tables(
    ROOT
)

with open(
    OUTDIR / "possible_phenotype_tables.tsv",
    "w",
    newline=""
) as f:

    fieldnames = [
        "path",
        "id_columns",
        "phenotype_like_columns",
        "all_columns"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(possible_tables)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

status_counts = Counter(
    r["positional_status"]
    for r in rows
)

state5_global_counts = Counter(
    r["global_state5"]
    for r in rows
)

state5_local_counts = Counter(
    r["local_state5"]
    for r in rows
)

state6_local_counts = Counter(
    r["local_state6"]
    for r in rows
)

offset_counts = Counter(
    r["local_offset_from_zinc"]
    for r in rows
    if r["local_offset_from_zinc"] is not None
)

changed = [
    r for r in rows
    if r["state_changed"]
]

disagreements = [
    r for r in rows
    if not r["global_local_same5"]
]

unique_zinc = sum(
    r["zinc_motif_count"] == 1
    for r in rows
)

strict_p = sum(
    bool(
        r["proline_immediately_after_zinc_motif"]
    )
    for r in rows
)

global_contig = sum(
    bool(r["global_contiguous5"])
    for r in rows
)

local_contig = sum(
    bool(r["local_contiguous5"])
    for r in rows
)

met_conserved = sum(
    bool(r["central_met_conserved"])
    for r in rows
)

terminal_w = sum(
    bool(r["terminal_w_conserved"])
    for r in rows
)

original_svmsy = sum(
    bool(r["original_exact_SVMSY"])
    for r in rows
)

corrected_svmsy = sum(
    bool(r["corrected_exact_SVMSY"])
    for r in rows
)

report = []

report.append(
    "SVMSY / SERRALYSIN MET-TURN POSITION AUDIT"
)

report.append("=" * 60)

report.append("")
report.append(f"Candidates: {len(rows)}")
report.append(
    f"Anchor: {anchor_header}"
)

report.append(
    f"Anchor zinc motif: "
    f"{anchor_zinc.group()} "
    f"at residue {anchor_zinc.start()+1}"
)

report.append(
    f"Anchor SVMSY: residues "
    f"{sv5_start0+1}-{sv5_start0+5}"
)

report.append(
    "Original SVMSY global-alignment columns: "
    + ",".join(
        str(x + 1)
        for x in global_cols5
    )
)

report.append("")
report.append(
    "INDEPENDENT CATALYTIC-LANDMARK CHECK"
)

report.append(
    f"Exactly one HEXXHXXGXXH zinc motif: "
    f"{unique_zinc}/{len(rows)}"
)

report.append(
    f"Proline immediately after zinc motif: "
    f"{strict_p}/{len(rows)}"
)

report.append(
    f"Original columns map to five contiguous "
    f"raw residues: "
    f"{global_contig}/{len(rows)}"
)

report.append(
    f"Independent catalytic-region mapping gives "
    f"five contiguous residues: "
    f"{local_contig}/{len(rows)}"
)

report.append(
    f"Central Met at homologous Met-turn position: "
    f"{met_conserved}/{len(rows)}"
)

report.append(
    f"Terminal W at homologous +5 position: "
    f"{terminal_w}/{len(rows)}"
)

report.append("")
report.append(
    "GLOBAL WHOLE-PROTEIN MSA VERSUS "
    "INDEPENDENT CATALYTIC-REGION ALIGNMENT"
)

report.append(
    f"Exact five-residue positional agreement: "
    f"{len(rows)-len(disagreements)}/{len(rows)}"
)

report.append(
    f"Disagreements: {len(disagreements)}"
)

report.append(
    f"Original exact-SVMSY calls: "
    f"{original_svmsy}"
)

report.append(
    f"Corrected exact-SVMSY calls: "
    f"{corrected_svmsy}"
)

report.append(
    f"SVMSY binary states changed: "
    f"{len(changed)}"
)

report.append("")
report.append("POSITIONAL STATUS")

for status, n in sorted(
    status_counts.items()
):
    report.append(
        f"{status}: {n}"
    )

report.append("")
report.append(
    "CORRECTED SIX-RESIDUE STATES"
)

for state, n in state6_local_counts.most_common():
    report.append(
        f"{state}: {n}"
    )

report.append("")
report.append(
    "ZINC-TO-MET-TURN START OFFSETS"
)

for offset, n in sorted(
    offset_counts.items()
):
    report.append(
        f"{offset}: {n}"
    )

report.append("")
report.append("FLAGGED / DISAGREEING SEQUENCES")

if not flagged and not disagreements:
    report.append("None.")
else:
    already = set()

    for r in flagged + disagreements:

        key = r["candidate_header"]

        if key in already:
            continue

        already.add(key)

        report.append(
            f"{key}\t"
            f"status={r['positional_status']}\t"
            f"global={r['global_state5']}\t"
            f"local={r['local_state5']}\t"
            f"zinc={r['zinc_motif']}\t"
            f"offset={r['local_offset_from_zinc']}"
        )

report.append("")
report.append(
    "CURRENT ANALYSIS / PHENOTYPE SOURCE SEARCH"
)

report.append(
    f"Relevant source-code/text hits: "
    f"{len(source_hits)}"
)

report.append(
    f"Possible phenotype tables: "
    f"{len(possible_tables)}"
)

for item in possible_tables[:30]:
    report.append(
        f"{item['path']}\t"
        f"ID=[{item['id_columns']}]\t"
        f"phenotype=["
        f"{item['phenotype_like_columns']}]"
    )

report.append("")
report.append("OUTPUTS")

for p in [
    calls_path,
    flagged_path,
    OUTDIR / "catalytic_regions.faa",
    OUTDIR / "structure_validation_panel.faa",
    OUTDIR / "structure_validation_panel.tsv",
    OUTDIR / "analysis_source_hits.tsv",
    OUTDIR / "possible_phenotype_tables.tsv",
]:
    report.append(str(p))

report.append("")
report.append("INTERPRETATION")

if (
    len(changed) == 0
    and len(disagreements) == 0
    and local_contig == len(rows)
):
    report.append(
        "The original exact-SVMSY binary calls are "
        "fully reproduced by an independent "
        "zinc-anchored catalytic-region alignment. "
        "The whole-protein MSA therefore did not "
        "change the state assignment."
    )

elif len(changed) == 0:
    report.append(
        "The exact-SVMSY binary classification is "
        "unchanged, although some positional details "
        "require review. Use corrected calls and the "
        "flag table for manuscript reporting."
    )

else:
    report.append(
        "One or more original exact-SVMSY calls change "
        "when the site is independently mapped through "
        "the zinc-anchored catalytic region. The "
        "corrected calls should replace the original "
        "binary state before the odds ratio is "
        "recalculated."
    )

report.append(
    "The structure_validation_panel.faa file contains "
    "the representative proteins needed for the final "
    "3D structural check against a serralysin "
    "reference structure."
)

report_path = OUTDIR / "audit_report.txt"
report_path.write_text(
    "\n".join(report) + "\n"
)

print("\n".join(report))
print()
print(
    f"Full report written to: {report_path}"
)
