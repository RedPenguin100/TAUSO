"""Generate MALAT1 2'-MOE gapmer ASOs and rank them with OligoAI.

Tiles the MALAT1 sense transcript into 20-mer 5-10-5 (MOE / DNA / MOE) full-PS gapmers,
formats every candidate into OligoAI's inhibition-model schema, runs OligoAI inference with
a trained checkpoint, and writes the candidates ranked by predicted knockdown.

OligoAI is a scorer, not a generator, so "design" here means enumerate-then-score: one
candidate per transcript position, each scored for percent target knockdown.

Run (from this directory):
    python design_malat1_oligoai.py --ckpt /path/to/OligoAI_11_09_25.ckpt

The inference step runs inside the `oligo_5090_hybrid` conda env (torch cu128 + flash-attn
2.8.4, required for RTX 5090 / Blackwell). Everything else is plain pandas and runs anywhere.
"""

import argparse
import subprocess
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DEFAULT_FASTA = HERE / "MALAT1.fasta"
DEFAULT_OLIGOAI_REPO = Path("/home/michael/career/tauso_article/OligoAI-fork")
DEFAULT_CKPT = Path("/mnt/c/Users/micha/Downloads/OligoAI_11_09_25.ckpt")
DEFAULT_ENV = "oligo_5090_hybrid"

# --- 2'-MOE "vanilla" gapmer design ---
ASO_LEN = 20
WING = 5  # MOE nucleotides on each end -> 5-10-5 MOE/DNA/MOE
FLANK = 50  # nt of sense context on each side of the binding site (OligoAI uses +/-50)
SUGAR_MODS = ["MOE"] * WING + ["DNA"] * (ASO_LEN - 2 * WING) + ["MOE"] * WING
BACKBONE_MODS = ["PS"] * (ASO_LEN - 1) + ["<PAD>"]  # full PS: 19 linkages, terminal pad

# Lipofection screening dose the model conditions on (nM); in-domain mode of OligoAI's data.
DEFAULT_DOSE = 100.0
DEFAULT_DELIVERY = "Lipofection"

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}


def revcomp(seq: str) -> str:
    return "".join(_COMP[b] for b in reversed(seq))


def read_fasta(path: Path) -> str:
    bases = [ln.strip() for ln in Path(path).read_text().splitlines() if ln and not ln.startswith(">")]
    return "".join(bases).upper()


def build_candidates(transcript: str, delivery: str, dose: float, limit: int | None = None) -> pd.DataFrame:
    """One 5-10-5 MOE gapmer per transcript position, in OligoAI's input schema."""
    rows = []
    last = len(transcript) - ASO_LEN
    for start in range(0, last + 1):
        site = transcript[start : start + ASO_LEN]  # sense target window (DNA)
        if set(site) - set("ACGT"):
            continue
        aso = revcomp(site)  # ASO, 5'->3', DNA letters
        up = transcript[max(0, start - FLANK) : start]
        down = transcript[start + ASO_LEN : start + ASO_LEN + FLANK]
        rna_context = (up + site + down).replace("T", "U")  # sense context, RNA letters
        rows.append(
            {
                "aso_sequence_5_to_3": aso,
                "inhibition_percent": 0.0,  # unknown; placeholder so OligoAI retains the row
                "chemistry": "2'MOE 5-10-5 gapmer, full PS",
                "custom_id": "MALAT1_2moe_design",
                "target_mrna": "MALAT1",
                "target_gene": "MALAT1",
                "cell_line": "",
                "cell_line_species": "human",
                "dosage": dose,
                "cells_per_well": "",
                "transfection_method": delivery,
                "steric_blocking": False,
                "rna_context": rna_context,
                "sugar_mods": SUGAR_MODS,
                "backbone_mods": BACKBONE_MODS,
                "split": "test",
                "target_start": start,  # 0-based position on the MALAT1 transcript
            }
        )
        if limit is not None and len(rows) >= limit:
            break
    return pd.DataFrame(rows)


def run_oligoai(input_csv: Path, pred_csv: Path, repo: Path, ckpt: Path, env: str, batch_size: int) -> None:
    cmd = [
        "conda", "run", "--no-capture-output", "-n", env,
        "python", str(repo / "run_inference.py"), str(input_csv),
        "--model_checkpoint", str(ckpt),
        "--device", "auto",
        "--batch_size", str(batch_size),
        "--output_path", str(pred_csv),
    ]
    print("Running OligoAI inference:\n  " + " ".join(cmd))
    subprocess.run(cmd, cwd=repo, check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    p.add_argument("--oligoai-repo", type=Path, default=DEFAULT_OLIGOAI_REPO)
    p.add_argument("--ckpt", type=Path, default=DEFAULT_CKPT, help="trained OligoAI checkpoint (.ckpt)")
    p.add_argument("--env", default=DEFAULT_ENV, help="conda env with OligoAI's deps")
    p.add_argument("--delivery", default=DEFAULT_DELIVERY, choices=["Lipofection", "Gymnosis", "Electroporation", "Other"])
    p.add_argument("--dose", type=float, default=DEFAULT_DOSE, help="dose the model conditions on (nM)")
    p.add_argument("--batch-size", type=int, default=128)
    p.add_argument("--limit", type=int, default=None, help="tile only the first N positions (quick test)")
    p.add_argument("--out-dir", type=Path, default=HERE)
    p.add_argument("--no-run-inference", action="store_true", help="build input only; skip OligoAI")
    args = p.parse_args()

    transcript = read_fasta(args.fasta)
    df = build_candidates(transcript, args.delivery, args.dose, args.limit)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    input_csv = args.out_dir / "MALAT1_2moe_oligoai_input.csv"
    df.to_csv(input_csv, index=False)
    print(f"MALAT1: {len(transcript)} nt -> {len(df)} candidate 2'-MOE gapmers "
          f"({args.delivery} @ {args.dose} nM) -> {input_csv.name}")

    if args.no_run_inference:
        return

    pred_csv = args.out_dir / "MALAT1_2moe_oligoai_input.with_predictions.csv"
    run_oligoai(input_csv, pred_csv, args.oligoai_repo, args.ckpt, args.env, args.batch_size)

    pred = pd.read_csv(pred_csv)
    ranked = pred.sort_values("predicted_inhibition_percent", ascending=False).reset_index(drop=True)
    ranked.insert(0, "rank", range(1, len(ranked) + 1))
    cols = ["rank", "target_start", "aso_sequence_5_to_3", "predicted_inhibition_percent",
            "chemistry", "transfection_method", "dosage"]
    final_csv = args.out_dir / "MALAT1_2moe_oligoai_ranked.csv"
    ranked[cols].to_csv(final_csv, index=False)
    print(f"\nRanked {len(ranked)} candidates -> {final_csv.name}")
    print(ranked[["rank", "target_start", "aso_sequence_5_to_3", "predicted_inhibition_percent"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
