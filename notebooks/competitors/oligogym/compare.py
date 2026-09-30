"""One table: the OligoGym baselines and the shipped TAUSO model on the same held-out test rows.

Reads what `run_baselines.py` and `score_tauso.py` wrote and prints them ranked by within-experiment
Spearman, with the gap to TAUSO alongside each baseline.

Run:  python notebooks/competitors/oligogym/compare.py
"""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS_DIR = HERE / "results"
METRICS = ("exp_med", "exp_mean", "gxc_med", "gxc_mean", "top5", "p5", "p10", "globP")


def main():
    # Every baseline file except the smoke run; a later successful entry replaces an earlier
    # failed one for the same combination, so a re-run of one featurizer folds in cleanly.
    baselines = {}
    for path in sorted(RESULTS_DIR.glob("oligogym_baselines*.json")):
        if path.stem.endswith("_smoke"):
            continue
        for entry in json.loads(path.read_text()):
            key = (entry["featurizer"], entry["model"], entry["variant"],
                   tuple(entry.get("extra", [])))
            if "error" not in entry or key not in baselines:
                baselines[key] = entry

    # Reference rows: the shipped TAUSO model, and OligoAI as the published middle rung.
    references = []
    for name in ("tauso_reference.json", "competitors_reference.json"):
        path = RESULTS_DIR / name
        if path.exists():
            references.extend(json.loads(path.read_text()))
    if not references:
        raise SystemExit("no reference rows; run score_tauso.py / score_oligoai.py first")

    rows = [(entry["model"], entry, entry["kind"]) for entry in references]
    for (featurizer, model, variant, extra), entry in baselines.items():
        label = f"{featurizer}/{model}/{variant}" + ("+cov" if extra else "")
        if "error" in entry:
            print(f"skipping failed combination {label}: {entry['error']}")
            continue
        kind = "sequence+chemistry+covariates" if extra else "sequence+chemistry"
        rows.append((label, entry, kind))
    rows.sort(key=lambda row: -row[1]["exp_med"])

    reference = max(entry["exp_med"] for _, entry, _ in rows)
    width = max(len(label) for label, _, _ in rows) + 2

    lines = [
        f"OligoGym baselines vs the shipped TAUSO model -- held-out test of the frozen OligoAI split",
        f"test rows: {references[0]['test_rows']}",
        "",
        "exp_med / exp_mean: within-experiment (custom_id) Spearman, median and mean",
        "gxc_*:              same, grouped by cohort (gene x cell line)",
        "top5:               mean actual inhibition of the top 5% predicted, median over experiments",
        "p5 / p10:           precision@5 / @10 within an experiment",
        "globP:              Pearson over all test rows at once",
        "",
        " " * width + " ".join(f"{m:>9}" for m in METRICS) + "   vs TAUSO   features",
    ]
    for label, scores, kind in rows:
        delta = scores["exp_med"] - reference
        gap = "  reference" if delta == 0 else f"{delta:+10.3f}"
        lines.append(f"{label:{width}s}" + " ".join(f"{scores[m]:>9.3f}" for m in METRICS)
                     + gap + f"   {kind}")

    text = "\n".join(lines)
    print(text)
    (RESULTS_DIR / "comparison.txt").write_text(text + "\n")
    print(f"\nsaved -> {RESULTS_DIR / 'comparison.txt'}")


if __name__ == "__main__":
    main()
