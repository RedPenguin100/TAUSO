# Normalizer-extraction pipeline (scripts)

Reproducible, **resumable** pipeline that turns the per-agent `batch_*.csv` curations
into a master table, caches each source patent's HTML, and ties every finding to the
exact place it was found. Designed so that if the session is shut down and the task is
re-asked, the cached HTML + master CSV mean it finishes in seconds, not hours.

All scripts use the tauso env python:
`/home/michael/miniforge3/envs/tauso/bin/python`

## End-to-end

```bash
cd notebooks/data/OligoAI/raw_data/patent_normalizers
PY=/home/michael/miniforge3/envs/tauso/bin/python

$PY scripts/consolidate.py        # batch_*.csv  -> patent_normalizer_annotations.csv, fetch_log.csv, SUMMARY.txt
$PY scripts/fetch_patents.py      # download full-text HTML -> html_cache/<id>.html.gz  (SKIPS already-cached)
$PY scripts/extract_evidence.py   # locate each quote in its cached source -> evidence.csv, evidence/<id>.md
```

## What each script does

| script | reads | writes | re-run cost |
|---|---|---|---|
| `consolidate.py` | `batch_*.csv` (+ dataset for coverage %) | `patent_normalizer_annotations.csv`, `fetch_log.csv`, `SUMMARY.txt` | instant |
| `fetch_patents.py` | master (or `--top N` / `--ids`) | `html_cache/<id>.html.gz` | **instant if cached** (network only for new ids) |
| `extract_evidence.py` | master + `html_cache/` | `evidence.csv`, `evidence/<id>.md` | instant (local parse) |

## Resuming / extending to more patents

The whole point of the cache is cheap re-runs. To add patents beyond the current 113:

1. Run new curation agents (see `../../../../.claude/skills/aso-patent-normalizer/SKILL.md`)
   writing `batch_<N>.csv` here, **or** just fetch raw HTML first:
   `$PY scripts/fetch_patents.py --top 343`   (caches every patent; skips the 113 already done)
2. `$PY scripts/consolidate.py` to fold new batches into the master.
3. `$PY scripts/extract_evidence.py` to regenerate the evidence map.

`fetch_patents.py` only downloads ids whose `html_cache/<id>.html.gz` is missing, so
re-running it is safe and fast. `--force` re-downloads.

## Why `curl`/urllib of the HTML and not WebFetch

Plain `WebFetch` of `patents.google.com/patent/<id>/en` truncates **before** the
Examples section, which is exactly where the normalization sentence lives. The full
HTML (1–7 MB) contains the OCR'd full text including Examples, so we fetch and grep
that. Image-only scanned patents (4 of 113) have no usable text even in the HTML;
those were read via the USPTO image PDF + `pdftoppm` OCR and are marked
`fetch_method=uspto_image_pdf_ocr` (their quote may not resolve in `evidence.csv`).

## Provenance columns

- `fetch_log.csv` — how each patent was obtained: `patent_id, target_gene, n_rows, fetch_method, source_url, confidence`.
- `evidence.csv` — where the finding lives in the source: `patent_id, normalizer_gene, normalization_basis, used_actb, cache_file, quote_found, char_offset, example_heading, keyword_hits, context_snippet`.
- `evidence/<id>.md` — per-patent: the recorded quote shown in its surrounding context near the located `Example N`, plus every normalization-keyword hit with offset + example (the raw audit trail behind the call).
