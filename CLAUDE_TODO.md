# Claude TODO

Low-priority tasks to pick up when there's nothing more pressing.

---

## Eternal SHA256 for DepMap files

`setup-depmap` currently does no integrity verification — it just checks if the
file exists. `setup-mrna-halflife` shows the right pattern: a hardcoded
`EXPECTED_SHA256` constant in source that is checked against the downloaded file.

The blocker is that DepMap doesn't publish stable, versioned SHA256 hashes
alongside their data files (their API has an inconsistent `checksum` column).
The work here is:

1. Decide where to store the canonical hashes — options:
   - Hardcoded constants in `cli.py` (like `setup_mrna_halflife` does)
   - A small committed `depmap_hashes.json` in the repo
   - A pinned entry in `CLAUDE.md` or similar

2. After deciding, add the known-good SHA256 for each DepMap 25Q3 file:
   - `Model.csv`
   - `OmicsProfiles.csv`
   - `OmicsExpressionTPMLogp1HumanAllGenesStranded.parquet` (post-conversion)

3. Wire them into `setup-depmap` the same way `should_skip_download` works.

This would prevent silent data drift when DepMap updates files within a quarter
without changing the release name (which is what caused the regression test
failures on the `feat/scavenger-receptor-expression` branch).
