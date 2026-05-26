#!/usr/bin/env python
"""Fetch Google Patents full-text HTML into a gzip cache. RESUMABLE.

This is the slow, network-bound step. Caching the raw HTML here means a later
re-run (or a fresh session re-asked to do this task) skips every already-fetched
patent and finishes in seconds instead of re-downloading hundreds of MB.

Cache layout:  <base>/html_cache/<PATENT_ID>.html.gz   (one file per patent)
A patent is considered "done" if its cache file exists and is > MIN_BYTES.

Usage (run with the tauso env python):
  /home/michael/miniforge3/envs/tauso/bin/python fetch_patents.py            # ids from master annotations
  ... fetch_patents.py --top 150                                             # top-N patents by row count
  ... fetch_patents.py --ids US20230002771A1,US20180312845A1                 # explicit ids
  ... fetch_patents.py --force                                              # re-download even if cached

Source of patent IDs (priority): --ids  >  --top N (from dataset)  >  master annotations CSV.
"""
import argparse, gzip, os, sys, time, urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)                       # .../patent_normalizers
CACHE = os.path.join(BASE, "html_cache")
MASTER = os.path.join(BASE, "patent_normalizer_annotations.csv")
DATASET = os.path.normpath(os.path.join(BASE, "..", "aso_inhibitions_with_canonical_gene.csv.gz"))
URL = "https://patents.google.com/patent/{pid}/en"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124 Safari/537.36"
MIN_BYTES = 50_000          # a real full-text page is ~1-7 MB; smaller => failed/blocked
DELAY = 0.7                 # politeness between requests (seconds)


def cache_path(pid):
    return os.path.join(CACHE, f"{pid}.html.gz")


def is_cached(pid):
    p = cache_path(pid)
    return os.path.exists(p) and os.path.getsize(p) > 2000  # gzip of a real page is >2KB


def fetch_one(pid, retries=3):
    """Download one patent's HTML, return raw bytes (decompressed) or None."""
    url = URL.format(pid=pid)
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            if len(data) < MIN_BYTES:
                raise ValueError(f"page too small ({len(data)} bytes) - likely blocked/empty")
            return data
        except (urllib.error.URLError, ValueError, TimeoutError) as e:
            print(f"  [{pid}] attempt {attempt}/{retries} failed: {e}", file=sys.stderr)
            time.sleep(2 * attempt)
    return None


def load_ids(args):
    if args.ids:
        return [s.strip() for s in args.ids.split(",") if s.strip()]
    if args.top:
        import pandas as pd
        df = pd.read_csv(DATASET, usecols=["patent_id"])
        return list(df["patent_id"].value_counts().head(args.top).index)
    import pandas as pd
    return list(pd.read_csv(MASTER, usecols=["patent_id"])["patent_id"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids")
    ap.add_argument("--top", type=int)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    os.makedirs(CACHE, exist_ok=True)
    ids = load_ids(args)
    print(f"{len(ids)} patents requested; cache dir = {CACHE}")

    done = fetched = failed = 0
    for i, pid in enumerate(ids, 1):
        if not args.force and is_cached(pid):
            done += 1
            continue
        print(f"[{i}/{len(ids)}] fetching {pid} ...")
        data = fetch_one(pid)
        if data is None:
            failed += 1
            print(f"  !! FAILED {pid} (try the USPTO image-PDF fallback by hand)", file=sys.stderr)
            continue
        with gzip.open(cache_path(pid), "wb") as fh:
            fh.write(data)
        fetched += 1
        time.sleep(DELAY)

    print(f"\nDONE. already-cached: {done} | newly-fetched: {fetched} | failed: {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
