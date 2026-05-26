#!/usr/bin/env python
"""Tie each normalizer finding to the EXACT place in the source patent.

For every patent in the master annotations that has a cached HTML page, this:
  1. strips HTML to plain text,
  2. locates the recorded `supporting_quote` inside that text (verbatim, then a
     whitespace/punctuation-normalized fuzzy match on the first ~14 words),
  3. records the char offset and the nearest preceding "Example N" heading,
  4. independently scans the text for every normalization-signal keyword and
     logs each hit's offset + example heading (the raw audit trail).

Outputs:
  evidence.csv                  one row per patent: where the quote was found + keyword-hit counts
  evidence/<PATENT_ID>.md       human-readable: located quote in context + full keyword-hit list

Re-run cost is ~zero (pure local parsing of the cache). Run AFTER fetch_patents.py.
  /home/michael/miniforge3/envs/tauso/bin/python extract_evidence.py
"""
import csv, glob, gzip, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
CACHE = os.path.join(BASE, "html_cache")
MASTER = os.path.join(BASE, "patent_normalizer_annotations.csv")
EVID_DIR = os.path.join(BASE, "evidence")
EVID_CSV = os.path.join(BASE, "evidence.csv")

# normalization-signal keywords -> regex (case-insensitive, on stripped text)
SIGNALS = [
    ("RIBOGREEN",       r"ribogreen"),
    ("total_RNA_norm",  r"adjusted according to total rna|normali[sz]ed?[^.]{0,40}total rna|total rna content|quantifying total rna"),
    ("GAPDH",           r"\bgapdh\b|\bg3pdh\b"),
    ("cyclophilin",     r"cyclophilin"),
    ("GUSB",            r"\bgusb\b|glucuronidase"),
    ("ACTB_bactin",     r"\bactb\b|β-?actin|b-actin|beta-actin"),
    ("bDNA_QuantiGene", r"\bbdna\b|quantigene|branched dna"),
    ("housekeeping",    r"housekeeping"),
]

TAG_RE = re.compile(r"<[^>]+>")
WS_RE = re.compile(r"\s+")
EX_RE = re.compile(r"example\s+\d+", re.IGNORECASE)


def strip_html(html):
    text = TAG_RE.sub(" ", html)
    text = (text.replace("&amp;", "&").replace("&lt;", "<").replace("&gt;", ">")
                .replace("&#39;", "'").replace("&quot;", '"').replace("&nbsp;", " "))
    return WS_RE.sub(" ", text).strip()


def norm(s):
    """lowercase, strip non-alphanumerics to spaces, collapse — for fuzzy matching."""
    return WS_RE.sub(" ", re.sub(r"[^a-z0-9 ]", " ", str(s).lower())).strip()


def nearest_example(ex_spans, offset):
    """Return the text of the closest 'Example N' heading at or before `offset`."""
    best = None
    for m in ex_spans:
        if m.start() <= offset:
            best = m
        else:
            break
    return best.group(0).strip().title() if best else ""


def load_html(pid):
    p = os.path.join(CACHE, f"{pid}.html.gz")
    if not os.path.exists(p):
        return None
    with gzip.open(p, "rb") as fh:
        return fh.read().decode("utf-8", errors="replace")


def locate_quote(text, ntext, ex_spans, quote):
    """Find quote (verbatim then fuzzy). Return (offset, how) or (-1, '')."""
    if not quote or quote.strip().lower() in ("", "nan"):
        return -1, ""
    i = text.find(quote)
    if i >= 0:
        return i, "verbatim"
    nq = norm(quote)
    if nq:
        j = ntext.find(nq)
        if j >= 0:
            return j, "normalized"
        words = nq.split()
        if len(words) >= 6:
            frag = " ".join(words[:14])
            k = ntext.find(frag)
            if k >= 0:
                return k, "fuzzy-first14"
    return -1, ""


def main():
    import pandas as pd
    os.makedirs(EVID_DIR, exist_ok=True)
    df = pd.read_csv(MASTER)
    rows = []
    located = missing_cache = quote_unfound = 0

    for _, r in df.iterrows():
        pid = r["patent_id"]
        html = load_html(pid)
        if html is None:
            missing_cache += 1
            rows.append(dict(patent_id=pid, normalizer_gene=r.get("normalizer_gene"),
                             normalization_basis=r.get("normalization_basis"), used_actb=r.get("used_actb"),
                             cache_file="", quote_found="no_cache", char_offset=-1, example_heading="",
                             keyword_hits="", context_snippet=""))
            continue

        text = strip_html(html)
        ntext = norm(text)
        ex_spans = list(EX_RE.finditer(text))

        # 1) locate the recorded supporting quote
        off, how = locate_quote(text, ntext, ex_spans, r.get("supporting_quote", ""))
        if off >= 0:
            located += 1
            # map offset back to original text for the example heading + context
            ex = nearest_example(ex_spans, off if how == "verbatim" else _map_to_text(text, ntext, off))
            base = off if how == "verbatim" else _map_to_text(text, ntext, off)
            context = text[max(0, base - 220): base + 220]
        else:
            quote_unfound += 1
            ex, context = "", ""

        # 2) independent keyword scan (audit trail)
        hit_lines, hit_counts = [], []
        low = text.lower()
        for label, pat in SIGNALS:
            ms = list(re.finditer(pat, low))
            if ms:
                hit_counts.append(f"{label}:{len(ms)}")
                for m in ms[:8]:
                    hit_lines.append((label, m.start(), nearest_example(ex_spans, m.start()),
                                      text[max(0, m.start() - 90): m.start() + 90]))

        rows.append(dict(patent_id=pid, normalizer_gene=r.get("normalizer_gene"),
                         normalization_basis=r.get("normalization_basis"), used_actb=r.get("used_actb"),
                         cache_file=f"html_cache/{pid}.html.gz",
                         quote_found=how if off >= 0 else "NOT_FOUND",
                         char_offset=off, example_heading=ex,
                         keyword_hits="; ".join(hit_counts), context_snippet=context))

        # per-patent evidence file
        with open(os.path.join(EVID_DIR, f"{pid}.md"), "w") as fh:
            fh.write(f"# {pid} — {r.get('target_gene')} ({r.get('cell_line')})\n\n")
            fh.write(f"- normalizer: **{r.get('normalizer_gene')}** (basis: {r.get('normalization_basis')}, used_actb: {r.get('used_actb')}, confidence: {r.get('confidence')})\n")
            fh.write(f"- source_url: {r.get('source_url')}\n")
            fh.write(f"- fetch_method: {r.get('fetch_method')}\n")
            fh.write(f"- cached source: html_cache/{pid}.html.gz\n\n")
            fh.write(f"## Recorded supporting quote ({how or 'NOT FOUND in cache'})\n\n")
            fh.write(f"> {r.get('supporting_quote')}\n\n")
            if context:
                fh.write(f"**Found near `{ex or 'n/a'}` (char offset {off}):**\n\n…{context}…\n\n")
            fh.write("## All normalization-keyword hits in this patent (offset · example · context)\n\n")
            for label, o, exh, ctx in hit_lines:
                fh.write(f"- `{label}` @ {o} · {exh or '-'} · …{ctx.strip()}…\n")
            if not hit_lines:
                fh.write("- (none — see source_url; may be an image-only patent read via PDF/OCR)\n")

    with open(EVID_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["patent_id", "normalizer_gene", "normalization_basis",
                                           "used_actb", "cache_file", "quote_found", "char_offset",
                                           "example_heading", "keyword_hits", "context_snippet"])
        w.writeheader()
        w.writerows(rows)

    print(f"patents: {len(df)} | quote located: {located} | quote not found: {quote_unfound} | no cache: {missing_cache}")
    print(f"wrote {EVID_CSV} and {len(df) - missing_cache} files in {EVID_DIR}/")


def _map_to_text(text, ntext, noff):
    """Approximate: map an offset in the normalized text back to the original text.
    norm() only substitutes chars 1:1 with spaces and lowercases (no length change
    except leading-strip), so offsets align closely; we nudge to the nearest space."""
    o = min(noff, len(text) - 1)
    # walk back to a word boundary for a clean context window
    while o > 0 and text[o] not in " .":
        o -= 1
    return o


if __name__ == "__main__":
    main()
