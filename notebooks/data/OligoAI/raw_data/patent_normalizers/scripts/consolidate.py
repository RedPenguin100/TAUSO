#!/usr/bin/env python
"""Merge the per-agent batch_*.csv files into the master annotations + fetch log + summary.

Reconciles the two schemas in play:
  - batch_1..15 : 12 columns (no `fetch_method`)  -> fetch_method derived from `notes`
  - batch_16..  : 13 columns (with `fetch_method`)

Outputs (all in the parent patent_normalizers/ folder):
  patent_normalizer_annotations.csv   master, one row per patent, 13 columns
  fetch_log.csv                       provenance: patent_id, target_gene, n_rows, fetch_method, source_url, confidence
  SUMMARY.txt                         tallies incl. the ACTB count

  /home/michael/miniforge3/envs/tauso/bin/python consolidate.py
"""
import glob, os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
DATASET = os.path.normpath(os.path.join(BASE, "..", "aso_inhibitions_with_canonical_gene.csv.gz"))
COLS = ["patent_id", "target_gene", "cell_line", "n_rows", "normalizer_gene", "normalization_basis",
        "assay_method", "used_actb", "confidence", "fetch_method", "source_url", "supporting_quote", "notes"]


def derive_fetch(notes):
    n = str(notes).lower()
    if any(k in n for k in ("image pdf", "image-only", "pdftoppm", "uspto image", "rendered")):
        return "uspto_image_pdf_ocr"
    if any(k in n for k in ("granted sibling", "granted patent of the identical", "identical application family", "granted sibling rather")):
        return "granted_sibling_html"
    return "google_patents_html_grep"


def main():
    frames = []
    for p in sorted(glob.glob(os.path.join(BASE, "batch_*.csv")),
                    key=lambda x: int(os.path.basename(x).split("_")[1].split(".")[0])):
        d = pd.read_csv(p)
        if "fetch_method" not in d.columns:
            d["fetch_method"] = d["notes"].apply(derive_fetch)
        frames.append(d)
    df = pd.concat(frames, ignore_index=True).reindex(columns=COLS)
    df = df.drop_duplicates("patent_id", keep="first").sort_values("n_rows", ascending=False).reset_index(drop=True)

    df.to_csv(os.path.join(BASE, "patent_normalizer_annotations.csv"), index=False)
    df[["patent_id", "target_gene", "n_rows", "fetch_method", "source_url", "confidence"]].to_csv(
        os.path.join(BASE, "fetch_log.csv"), index=False)

    total_rows = int(pd.read_csv(DATASET, usecols=["patent_id"]).shape[0]) if os.path.exists(DATASET) else 173191
    cov = int(df["n_rows"].sum())
    actb = df[df["used_actb"].astype(str).str.lower() == "yes"]
    lines = []
    lines.append(f"patents annotated: {len(df)} / 343")
    lines.append(f"rows covered: {cov} / {total_rows} = {cov / total_rows * 100:.1f}%")
    lines.append("")
    for col in ["normalization_basis", "assay_method", "confidence", "fetch_method", "used_actb"]:
        lines.append(f"== {col} ==")
        lines.append(df[col].value_counts(dropna=False).to_string())
        lines.append("")
    lines.append(f"### ACTB-as-normalizer: {len(actb)} patent(s)")
    for _, r in actb.iterrows():
        lines.append(f"  {r['patent_id']} {r['target_gene']} (n={r['n_rows']})")
    lines.append("")
    hk = df[df["normalization_basis"] == "housekeeping_gene"]
    lines.append(f"### housekeeping_gene: {len(hk)} patent(s), {int(hk['n_rows'].sum())} rows")
    for _, r in hk.iterrows():
        lines.append(f"  {r['patent_id']} {r['target_gene']} -> {r['normalizer_gene']} (n={r['n_rows']}, {r['confidence']})")

    out = "\n".join(lines)
    with open(os.path.join(BASE, "SUMMARY.txt"), "w") as fh:
        fh.write(out + "\n")
    print(out)


if __name__ == "__main__":
    main()
