import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import pytesseract


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    args = ap.parse_args()

    root = Path(args.project_root)
    res = root / "preprint" / "results" / "auto"
    figs = root / "preprint" / "figures" / "auto"

    phen = pd.read_csv(
        res / "Figure2_phenotypes_auto.tsv",
        sep="\t",
        dtype=str
    ).fillna("")

    image_path = root / "preprint" / "published_data" / "Maier2020_Figure2.png"
    image = Image.open(image_path).convert("RGB")

    # Calibrated against the original Maier Figure 2 (974 x 1700 px).
    # The printed operon-type numbers occupy approximately x = 855–875.
    ref_w, ref_h = 974, 1700
    sx = image.width / ref_w
    sy = image.height / ref_h

    x0 = int(round(855 * sx))
    x1 = int(round(875 * sx))
    y0 = int(round(15 * sy))
    y1 = int(round(1695 * sy))

    crop = image.crop((x0, y0, x1, y1)).convert("L")

    scale = 12
    up = crop.resize((crop.width * scale, crop.height * scale))

    arr = np.asarray(up)
    # Keep black type-number text; suppress pale operon arrows and grey rules.
    arr = np.where(arr < 175, 0, 255).astype(np.uint8)
    proc = Image.fromarray(arr)

    ocr = pytesseract.image_to_data(
        proc,
        config="--psm 6 -c tessedit_char_whitelist=0123456789",
        output_type=pytesseract.Output.DATAFRAME
    )

    ocr = ocr.dropna(subset=["text"])
    tokens = []

    for _, r in ocr.iterrows():
        txt = str(r["text"]).strip()
        nums = re.findall(r"\d{1,2}", txt)

        if not nums:
            continue

        value = int(nums[0])

        # Maier Figure 2 uses operon types 1–22.
        if not 1 <= value <= 22:
            continue

        y_center_up = float(r["top"]) + float(r["height"]) / 2
        y_center = y0 + y_center_up / scale

        tokens.append({
            "operon_type": str(value),
            "type_y": y_center,
            "ocr_conf": float(r["conf"]),
            "ocr_raw": txt
        })

    tokens = pd.DataFrame(tokens)

    if tokens.empty:
        raise RuntimeError(
            "No operon type numbers were detected. "
            "Check Figure2_operon_type_column_processed.png."
        )

    # Save processed crop for transparent QC.
    proc.save(figs / "Figure2_operon_type_column_processed.png")

    # Assign each detected type number to the nearest strain row.
    strain_rows = phen.copy()
    strain_rows["y_num"] = pd.to_numeric(strain_rows["y"], errors="coerce")

    assigned = []
    used_token_indices = set()

    for _, row in strain_rows.iterrows():
        y = float(row["y_num"])

        distances = (tokens["type_y"] - y).abs()

        if len(distances) == 0:
            continue

        idx = distances.idxmin()
        distance = float(distances.loc[idx])

        # Typical row spacing is only a few pixels. 4.5 px is strict enough
        # to avoid assigning a type number from a neighbouring strain.
        if distance <= 4.5 and idx not in used_token_indices:
            t = tokens.loc[idx]
            used_token_indices.add(idx)

            assigned.append({
                "strain": row["strain"],
                "strain_key": row.get("strain_key", ""),
                "figure_y": y,
                "operon_type": t["operon_type"],
                "type_y": t["type_y"],
                "y_distance": distance,
                "ocr_conf": t["ocr_conf"],
                "ocr_raw": t["ocr_raw"]
            })

    out = pd.DataFrame(assigned)

    out.to_csv(
        res / "Figure2_operon_types_auto.tsv",
        sep="\t",
        index=False
    )

    # Also save unassigned OCR tokens for QC.
    unused = tokens.loc[
        [i for i in tokens.index if i not in used_token_indices]
    ].copy()

    unused.to_csv(
        res / "Figure2_operon_type_unassigned_tokens.tsv",
        sep="\t",
        index=False
    )

    overlay = image.copy()
    draw = ImageDraw.Draw(overlay)

    # Show OCR column.
    draw.rectangle(
        (x0, y0, x1, y1),
        outline="black",
        width=1
    )

    # Mark assigned numbers at their row positions.
    for _, r in out.iterrows():
        y = int(round(float(r["figure_y"])))
        draw.line((x0 - 8, y, x0 - 2, y), fill="black", width=1)
        draw.rectangle(
            (x0, y - 4, x1, y + 4),
            outline="black",
            width=1
        )

    overlay.save(
        figs / "Figure2_operon_type_QC.png",
        dpi=(300, 300)
    )

    print("Detected numeric OCR tokens:", len(tokens))
    print("Assigned operon types to strain rows:", len(out))
    print("Unassigned numeric tokens:", len(unused))

    print("\nType counts:")
    if not out.empty:
        counts = out["operon_type"].value_counts()
        counts.index = counts.index.astype(int)
        counts = counts.sort_index()
        print(counts.to_string())
    else:
        print("None")

    print("\nMedian strain/type y-distance:",
          out["y_distance"].median() if not out.empty else "NA")

    print("\nQC files:")
    print(figs / "Figure2_operon_type_QC.png")
    print(figs / "Figure2_operon_type_column_processed.png")


if __name__ == "__main__":
    main()
