# Packing the feature cache

`notebooks/features/pack_features.py` builds the wide feature-cache parquet that
`tauso setup-features` downloads from Zenodo. Normal users never run it -- they just
`setup-features`. Only use it to publish a new feature-cache version.

## Usage

```bash
python -m notebooks.features.pack_features --run oligo \
    --source notebooks/saved_features_oligo \
    --out feature_pack_build
```

Input: a directory of per-feature CSVs keyed on `index_<run>`.
Output: `<run>_features.parquet` (+ its md5, printed to stdout).

After uploading the parquet to Zenodo, set `ZENODO_FEATURES_RECORD` and the matching
`FEATURE_CACHE_FILES[<run>]` entry in `src/tauso/populate/feature_cache.py` to the new
record id + md5. Add a new key (e.g. `oligo_v2`) to ship a new run alongside.
