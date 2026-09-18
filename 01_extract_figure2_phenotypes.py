import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from docx import Document
from rapidfuzz import fuzz, process
import pytesseract


def norm_text(s):
    s = str(s)
    replacements = {
        "$": "5", "OSM": "DSM", "D5M": "DSM",
        "P,": "P.", "sp,": "sp.", "  ": " "
    }
    for a, b in replacements.items():
        s = s.replace(a, b)
    s = re.sub(r"\s+", " ", s).strip()
    s = s.replace("(T)", " T").replace("™", "")
    return s


def norm_match(s):
    s = norm_text(s).lower()
    s = s.replace("pseudomonas", "p.")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


PREFIXES = [
    "WS", "DSM", "LMG", "ATCC", "ICMP", "CFBP", "CCUG", "JCM",
    "NBRC", "NCTC", "CECT", "CIP", "KCTC", "BCRC", "MT"
]


def strain_key(s):
    t = norm_text(s).upper()
    for prefix in PREFIXES:
        m = re.search(rf"\b{re.escape(prefix)}\s*[-:]?\s*(\d+)\b", t)
        if m:
            return f"{prefix}{m.group(1)}"
    # WS can be OCRed with punctuation around it.
    m = re.search(r"\bWS\s*(\d+)\b", t)
    if m:
        return f"WS{m.group(1)}"
    return ""


def table2_strains(docx_path):
    doc = Document(docx_path)
    if len(doc.tables) < 2:
        raise RuntimeError("Expected at least two tables in the Maier DOCX.")
    table = doc.tables[1]
    vals = []
    for row in table.rows[1:]:
        if not row.cells:
            continue
        s = norm_text(row.cells[0].text)
        if s:
            vals.append(s)
    return vals


def ocr_lines(image, cfg):
    w, h = image.size
    rw, rh = cfg["figure2_reference_width"], cfg["figure2_reference_height"]
    x0, y0, x1, y1 = cfg["figure2_ocr_crop"]
    x0, x1 = int(x0 * w / rw), int(x1 * w / rw)
    y0, y1 = int(y0 * h / rh), int(y1 * h / rh)

    crop = image.crop((x0, y0, x1, y1))
    scale = 4
    up = crop.resize((crop.width * scale, crop.height * scale)).convert("L")

    df = pytesseract.image_to_data(
        up, config="--psm 6", output_type=pytesseract.Output.DATAFRAME
    )
    df = df.dropna(subset=["text"])
    df = df[df["conf"] >= cfg["ocr_min_confidence"]]

    lines = []
    if df.empty:
        return lines

    for _, g in df.groupby(["block_num", "par_num", "line_num"], sort=False):
        g = g.sort_values("left")
        text = norm_text(" ".join(g["text"].astype(str)))
        if not text:
            continue
        top = g["top"].min() / scale + y0
        bottom = (g["top"] + g["height"]).max() / scale + y0
        y = (top + bottom) / 2
        conf = float(g["conf"].mean())
        lines.append({"ocr_text": text, "y": y, "ocr_conf": conf})
    return lines


def match_lines(lines, expected, cfg):
    expected_norm = {s: norm_match(s) for s in expected}
    key_to_expected = {}
    for s in expected:
        k = strain_key(s)
        if k:
            key_to_expected.setdefault(k, []).append(s)

    used = set()
    out = []

    for line in sorted(lines, key=lambda x: x["y"]):
        txt = line["ocr_text"]
        k = strain_key(txt)

        matched = None
        score = 0.0
        method = ""

        if k and len(key_to_expected.get(k, [])) == 1:
            cand = key_to_expected[k][0]
            if cand not in used:
                matched = cand
                score = 100.0
                method = "identifier"

        if matched is None:
            q = norm_match(txt)
            choices = {s: v for s, v in expected_norm.items() if s not in used}
            if choices and q:
                best = process.extractOne(
                    q, choices, scorer=fuzz.token_set_ratio, score_cutoff=cfg["fuzzy_match_min_score"]
                )
                if best is not None:
                    # rapidfuzz with mapping returns (value, score, key)
                    _, score, cand = best
                    matched = cand
                    method = "fuzzy"

        # Only retain lines that plausibly represent one of the expected phenotype strains.
        if matched:
            used.add(matched)
            row = dict(line)
            row.update({
                "strain": matched,
                "strain_key": strain_key(matched),
                "match_score": float(score),
                "match_method": method
            })
            out.append(row)

    unmatched = [s for s in expected if s not in used]
    return out, unmatched


def sample_rgb(arr, x, y, radius=2):
    h, w, _ = arr.shape
    x, y = int(round(x)), int(round(y))
    patch = arr[max(0, y-radius):min(h, y+radius+1),
                max(0, x-radius):min(w, x+radius+1), :]
    if patch.size == 0:
        return np.array([255, 255, 255], dtype=float)
    return np.median(patch.reshape(-1, 3), axis=0)


def classify_color(rgb, palette):
    rgb = np.asarray(rgb, dtype=float)
    # White / near-white means no phenotype cell is shown.
    if rgb.mean() > 245 and (rgb.max() - rgb.min()) < 20:
        return "missing", np.nan
    d = {name: float(np.linalg.norm(rgb - np.asarray(col, dtype=float)))
         for name, col in palette.items()}
    best = min(d, key=d.get)
    vals = sorted(d.values())
    margin = vals[1] - vals[0] if len(vals) > 1 else np.nan
    # Very distant from all palette colors is treated as uncertain rather than forced.
    if d[best] > 95:
        return "uncertain", margin
    return best, margin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    auto = root / "preprint" / "automation"
    cfg = json.loads((auto / "config.json").read_text())

    image_path = root / "preprint" / "published_data" / "Maier2020_Figure2.png"
    docx_path = root / "preprint" / "published_data" / "Maier2020_Table1.docx"
    outdir = root / "preprint" / "results" / "auto"
    figdir = root / "preprint" / "figures" / "auto"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    if not image_path.exists():
        raise FileNotFoundError(image_path)
    if not docx_path.exists():
        raise FileNotFoundError(docx_path)

    image = Image.open(image_path).convert("RGB")
    arr = np.asarray(image)
    expected = table2_strains(docx_path)
    lines = ocr_lines(image, cfg)
    matched, unmatched = match_lines(lines, expected, cfg)

    w, h = image.size
    xs = [
        x * w / cfg["figure2_reference_width"]
        for x in cfg["figure2_heatmap_x_centers"]
    ]
    palette = cfg["phenotype_palette"]
    score_map = cfg["phenotype_scores"]

    rows = []
    for m in matched:
        row = dict(m)
        for day, x in zip([3, 4, 7], xs):
            rgb = sample_rgb(arr, x, m["y"])
            label, margin = classify_color(rgb, palette)
            row[f"day{day}_category"] = label
            row[f"day{day}_score"] = score_map.get(label, np.nan)
            row[f"day{day}_rgb"] = ",".join(str(int(v)) for v in rgb)
            row[f"day{day}_color_margin"] = margin
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "Figure2_phenotypes_auto.tsv", sep="\t", index=False)

    pd.DataFrame({"strain": unmatched, "strain_key": [strain_key(x) for x in unmatched]}).to_csv(
        outdir / "Figure2_unmatched_strains.tsv", sep="\t", index=False
    )
    pd.DataFrame(lines).to_csv(outdir / "Figure2_all_OCR_lines.tsv", sep="\t", index=False)

    # QC overlay
    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)
    for r in rows:
        y = int(round(r["y"]))
        for day, x in zip([3,4,7], xs):
            x = int(round(x))
            draw.ellipse((x-3, y-3, x+3, y+3), outline="black", width=1)
        # small row marker at left edge of OCR crop
        draw.line((270*w/cfg["figure2_reference_width"], y,
                   278*w/cfg["figure2_reference_width"], y), fill="black", width=1)
    overlay.save(figdir / "Figure2_phenotype_QC.png", dpi=(300,300))

    print(f"Expected strains from Supplementary Table 2: {len(expected)}")
    print(f"Matched to Figure 2 OCR: {len(df)}")
    print(f"Unmatched: {len(unmatched)}")
    if not df.empty:
        for day in [3,4,7]:
            print(f"Day {day}:")
            print(df[f"day{day}_category"].value_counts(dropna=False).to_string())
    print("QC overlay:", figdir / "Figure2_phenotype_QC.png")


if __name__ == "__main__":
    main()
