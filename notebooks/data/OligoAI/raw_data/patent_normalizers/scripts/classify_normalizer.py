#!/usr/bin/env python
"""Programmatically classify the normalizer for EVERY cached patent, to reach 100% coverage.

Manual agent curation (113 patents) stays authoritative. For the remaining ~230 patents
(overwhelmingly the homogeneous early-2000s Isis family) this extracts the normalization
sentence(s) straight from the cached HTML and classifies the basis. It then:
  1. validates itself against the 113 manual labels (prints agreement),
  2. builds a full-coverage master = manual rows + auto rows for the rest,
  3. tallies ACTB-as-normalizer across all 343.

Heuristic: find "<...> (m)RNA levels were normalized/adjusted to/according to/using <X>"
and "normalized ... by quantifying total RNA using RiboGreen", classify <X>, and prefer
in-vitro matches (downgrade windows mentioning mouse / in vivo / transgenic / Northern).

Outputs (parent folder):
  classified_normalizers.csv              auto result for every cached patent (+ evidence sentence)
  patent_normalizer_annotations_full.csv  113 manual + auto for the rest (col `source`)
  SUMMARY_full.txt                         full-coverage tallies incl. ACTB

  /home/michael/miniforge3/envs/tauso/bin/python classify_normalizer.py
"""
import csv, glob, gzip, os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
CACHE = os.path.join(BASE, "html_cache")
MASTER = os.path.join(BASE, "patent_normalizer_annotations.csv")
DATASET = os.path.normpath(os.path.join(BASE, "..", "aso_inhibitions_with_canonical_gene.csv.gz"))

TAG_RE = re.compile(r"<[^>]+>")
WS_RE = re.compile(r"\s+")
EX_RE = re.compile(r"example\s+\d+", re.IGNORECASE)
# normalization sentence: capture the object of "normalized/adjusted to/according to/using ..."
NORM_RE = re.compile(
    r"(?:(m?RNA|transcript)\s+levels?\s+(?:were|are|was)\s+)?"
    r"(?:normali[sz]ed|adjusted)\s+(?:to|according to|using|by quantifying)\s+([^.;:]{3,140})",
    re.IGNORECASE)
CONTEXT_BAD = re.compile(r"\b(mouse|murine|in vivo|transgenic|northern|western|protein|liver of|spinal|tissue|animals?)\b", re.I)


def strip_html(html):
    t = TAG_RE.sub(" ", html)
    for a, b in [("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">"), ("&#39;", "'"), ("&quot;", '"'), ("&nbsp;", " ")]:
        t = t.replace(a, b)
    return WS_RE.sub(" ", t).strip()


def classify_object(obj):
    o = obj.lower()
    if "either" in o and ("gapdh" in o or "ribogreen" in o or "total rna" in o):
        return "either_or"        # the ambiguous Isis boilerplate
    if "ribogreen" in o or "total rna" in o:
        return "total_RNA"
    if "gusb" in o or "glucuronidase" in o:
        return "GUSB"
    if re.search(r"\bactb\b|β-?actin|b-actin|beta-actin", o):
        return "ACTB"
    if re.search(r"\bgapdh\b|\bg3pdh\b", o):
        return "GAPDH"
    if "cyclophilin" in o:
        return "cyclophilin"
    return None


def nearest_example(ex_spans, off):
    best = None
    for m in ex_spans:
        if m.start() <= off:
            best = m
        else:
            break
    return best.group(0).strip().title() if best else ""


def analyze(text):
    """Return (basis, normalizer, confidence, actb_flag, evidence_sentence, example)."""
    ex_spans = list(EX_RE.finditer(text))
    counts = {}          # type -> count (in-vitro-ish)
    counts_all = {}
    first_ev = {}        # type -> (sentence, example)
    actb_flag = False
    actb_ev = None
    for m in NORM_RE.finditer(text):
        obj = m.group(2).strip()
        kind = classify_object(obj)
        if not kind:
            continue
        window = text[max(0, m.start() - 160): m.end() + 40]
        invitro = not CONTEXT_BAD.search(window)
        counts_all[kind] = counts_all.get(kind, 0) + 1
        if invitro:
            counts[kind] = counts.get(kind, 0) + 1
            first_ev.setdefault(kind, (window.strip(), nearest_example(ex_spans, m.start())))
        if kind == "ACTB":
            # ACTB as an mRNA normalizer (any context) — record, flag only if in-vitro mRNA
            if actb_ev is None:
                actb_ev = (window.strip(), nearest_example(ex_spans, m.start()), invitro)
            if invitro:
                actb_flag = True

    has_ribo = "ribogreen" in text.lower()

    # decide basis, preferring committed in-vitro target-gene normalization
    def ev(kind):
        return first_ev.get(kind, ("", ""))

    if counts.get("total_RNA", 0) > 0:
        basis, norm, conf = "total_RNA_ribogreen", "total RNA", "high"
        sent, exh = ev("total_RNA")
    elif counts.get("GUSB", 0) > 0:
        basis, norm, conf = "housekeeping_gene", "GUSB", "high"; sent, exh = ev("GUSB")
    elif counts.get("ACTB", 0) > 0 and counts.get("ACTB", 0) >= counts.get("GAPDH", 0):
        basis, norm, conf = "housekeeping_gene", "ACTB", "medium"; sent, exh = ev("ACTB")
    elif counts.get("either_or", 0) > 0 or counts_all.get("either_or", 0) > 0:
        basis, norm, conf = "unknown", "GAPDH or total RNA (boilerplate)", "medium"; sent, exh = ev("either_or") if ev("either_or")[0] else ("", "")
    elif counts.get("GAPDH", 0) > 0:
        basis, norm, conf = "housekeeping_gene", "GAPDH", "medium"; sent, exh = ev("GAPDH")
    elif has_ribo:
        basis, norm, conf = "total_RNA_ribogreen", "total RNA", "medium"; sent, exh = "", ""
    elif counts.get("cyclophilin", 0) > 0:
        basis, norm, conf = "housekeeping_gene", "cyclophilin A", "low"; sent, exh = ev("cyclophilin")
    else:
        basis, norm, conf = "unknown", "unknown", "low"; sent, exh = "", ""
    return basis, norm, conf, actb_flag, sent[:300], exh, actb_ev


def main():
    vc = pd.read_csv(DATASET, usecols=["patent_id"])["patent_id"].value_counts()
    meta = pd.read_csv(DATASET, usecols=["patent_id", "target_gene", "cell_line"])
    gene = meta.groupby("patent_id")["target_gene"].agg(lambda s: s.value_counts().index[0] if len(s.value_counts()) else "")
    cell = meta.groupby("patent_id")["cell_line"].agg(lambda s: s.value_counts().index[0] if len(s.dropna().value_counts()) else "")

    rows = []
    for p in sorted(glob.glob(os.path.join(CACHE, "*.html.gz"))):
        pid = os.path.basename(p)[:-8]
        with gzip.open(p, "rb") as fh:
            text = strip_html(fh.read().decode("utf-8", errors="replace"))
        basis, norm, conf, actb_flag, sent, exh, actb_ev = analyze(text)
        rows.append(dict(patent_id=pid, n_rows=int(vc.get(pid, 0)),
                         normalizer_auto=norm, basis_auto=basis, confidence_auto=conf,
                         actb_as_normalizer=actb_flag,
                         actb_evidence=(actb_ev[0][:200] if actb_ev else ""),
                         example_heading=exh, evidence_sentence=sent))
    auto = pd.DataFrame(rows).sort_values("n_rows", ascending=False)
    auto.to_csv(os.path.join(BASE, "classified_normalizers.csv"), index=False)

    # ---- validate against manual ----
    man = pd.read_csv(MASTER)
    j = man.merge(auto, on="patent_id", how="inner")
    agree = (j["normalization_basis"] == j["basis_auto"]).mean()
    print(f"VALIDATION on {len(j)} manually-curated patents:")
    print(f"  basis agreement: {agree*100:.1f}%")
    print("  confusion (manual basis -> auto basis):")
    print(j.groupby(["normalization_basis", "basis_auto"]).size().to_string())
    man_actb = set(man[man["used_actb"].astype(str).str.lower() == "yes"]["patent_id"])
    auto_actb = set(auto[auto["actb_as_normalizer"]]["patent_id"])
    print(f"\n  manual ACTB-yes: {sorted(man_actb)}")
    print(f"  auto  ACTB-flag (intersection w/ manual set): {sorted(auto_actb & set(man['patent_id']))}")

    # ---- build full-coverage master ----
    man2 = man.copy(); man2["source"] = "manual_curation"
    have = set(man["patent_id"])
    add = auto[~auto["patent_id"].isin(have)].copy()
    add["target_gene"] = add["patent_id"].map(gene)
    add["cell_line"] = add["patent_id"].map(cell)
    add = add.rename(columns={"normalizer_auto": "normalizer_gene", "basis_auto": "normalization_basis",
                              "confidence_auto": "confidence"})
    add["assay_method"] = "RT-qPCR"
    add["used_actb"] = add["actb_as_normalizer"].map({True: "yes", False: "no"})
    add["fetch_method"] = "google_patents_html_grep"
    add["source_url"] = "https://patents.google.com/patent/" + add["patent_id"] + "/en"
    add["supporting_quote"] = add["evidence_sentence"]
    add["notes"] = "auto-classified from cached HTML; " + add["example_heading"].fillna("")
    add["source"] = "auto_classified"
    cols = list(man2.columns)
    full = pd.concat([man2, add.reindex(columns=cols)], ignore_index=True).sort_values("n_rows", ascending=False)
    full.to_csv(os.path.join(BASE, "patent_normalizer_annotations_full.csv"), index=False)

    total = int(vc.sum())
    cov = int(full["n_rows"].sum())
    actb_full = full[full["used_actb"].astype(str).str.lower() == "yes"]
    lines = [f"FULL COVERAGE: {len(full)} / 343 patents | rows {cov}/{total} = {cov/total*100:.1f}%",
             f"  (manual: {len(man2)}, auto-classified: {len(add)})", ""]
    for c in ["normalization_basis", "confidence", "source", "used_actb"]:
        lines += [f"== {c} ==", full[c].value_counts(dropna=False).to_string(), ""]
    lines.append(f"### ACTB-as-normalizer (full): {len(actb_full)} patent(s)")
    for _, r in actb_full.iterrows():
        lines.append(f"  {r['patent_id']} {r['target_gene']} (n={r['n_rows']}, {r['source']})")
    out = "\n".join(lines)
    with open(os.path.join(BASE, "SUMMARY_full.txt"), "w") as fh:
        fh.write(out + "\n")
    print("\n" + out)


if __name__ == "__main__":
    main()
