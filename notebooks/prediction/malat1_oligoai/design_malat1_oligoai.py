"""Design MALAT1 2'-MOE gapmer ASOs with TAUSO and rank them with the OligoAI model.

TAUSO tiles the transcript into 20-mer 5-10-5 (MOE/DNA/MOE) full-PS gapmers; each candidate is
formatted for OligoAI's inhibition model, scored, and written ranked by predicted knockdown.
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from tauso.aso_generation import default_config, get_initial_data
from tauso.common.modifications import transform_linkage_to_oligo, transform_pattern_to_oligo
from tauso.data.consts import ASO_SEQUENCE
from tauso.populate.calculators.cache import AssetCache

HERE = Path(__file__).resolve().parent

GENE = "MALAT1"
ASO_LEN = 20
FLANK = 50  # nt of sense context on each side of the site (OligoAI's +/-50 window)
SUGAR_MODS = transform_pattern_to_oligo(default_config().standard_chemical_pattern)
BACKBONE_MODS = transform_linkage_to_oligo("else", ASO_LEN)  # "else" => every linkage is PS
CHEMISTRY = "2'MOE 5-10-5 gapmer, full PS"

INPUT_COLS = ["aso_sequence_5_to_3", "sugar_mods", "backbone_mods", "rna_context",
              "dosage", "transfection_method", "inhibition_percent", "custom_id", "split", "target_start"]
OUTPUT_COLS = ["rank", "target_start", "aso_sequence_5_to_3", "predicted_inhibition_percent",
               "chemistry", "transfection_method", "dosage"]


def build_candidates(delivery, dose, genome="GRCh38"):
    transcript = str(AssetCache(genome=genome).get_full_gene_data()[GENE].full_mrna)
    candidates = get_initial_data(transcript, aso_sizes=[ASO_LEN], canonical_name=GENE)
    rows = []
    for start, aso in enumerate(candidates[ASO_SEQUENCE]):
        up = transcript[max(0, start - FLANK):start]
        down = transcript[start + ASO_LEN:start + ASO_LEN + FLANK]
        context = (up + transcript[start:start + ASO_LEN] + down).replace("T", "U")
        rows.append({
            "aso_sequence_5_to_3": aso,
            "sugar_mods": SUGAR_MODS,
            "backbone_mods": BACKBONE_MODS,
            "rna_context": context,
            "dosage": dose,
            "transfection_method": delivery,
            "inhibition_percent": 0.0,  # placeholder label; OligoAI drops rows with a null label
            "custom_id": f"{GENE}_2moe_design",
            "split": "test",
            "target_start": start,
        })
    return pd.DataFrame(rows)


def score_with_oligoai(input_csv, pred_csv, repo, ckpt, env, batch_size):
    subprocess.run(
        ["conda", "run", "--no-capture-output", "-n", env,
         "python", str(Path(repo) / "run_inference.py"), str(input_csv),
         "--model_checkpoint", str(ckpt), "--device", "auto",
         "--batch_size", str(batch_size), "--output_path", str(pred_csv)],
        cwd=repo, check=True,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--delivery", default="Lipofection",
                   choices=["Lipofection", "Gymnosis", "Electroporation", "Other"])
    p.add_argument("--dose", type=float, default=100.0, help="dose OligoAI conditions on (nM)")
    p.add_argument("--out-dir", type=Path, default=HERE)
    p.add_argument("--oligoai-repo", default=os.environ.get("OLIGOAI_REPO"),
                   help="path to the OligoAI checkout (run_inference.py); or set OLIGOAI_REPO")
    p.add_argument("--ckpt", default=os.environ.get("OLIGOAI_CKPT"),
                   help="trained OligoAI checkpoint; or set OLIGOAI_CKPT")
    p.add_argument("--conda-env", default=os.environ.get("OLIGOAI_ENV", "oligo_5090_hybrid"),
                   help="conda env holding OligoAI's dependencies")
    p.add_argument("--batch-size", type=int, default=128)
    args = p.parse_args()

    candidates = build_candidates(args.delivery, args.dose)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not (args.oligoai_repo and args.ckpt):
        input_csv = args.out_dir / f"{GENE}_2moe_oligoai_input.csv"
        candidates[INPUT_COLS].to_csv(input_csv, index=False)
        print(f"Built {len(candidates)} candidates -> {input_csv}")
        print("Pass --oligoai-repo and --ckpt (or set OLIGOAI_REPO / OLIGOAI_CKPT) to score and rank.")
        return

    with tempfile.TemporaryDirectory() as work:
        input_csv, pred_csv = Path(work) / "input.csv", Path(work) / "pred.csv"
        candidates[INPUT_COLS].to_csv(input_csv, index=False)
        score_with_oligoai(input_csv, pred_csv, args.oligoai_repo, args.ckpt, args.conda_env, args.batch_size)
        pred = pd.read_csv(pred_csv)

    ranked = pred.sort_values("predicted_inhibition_percent", ascending=False, ignore_index=True)
    ranked.insert(0, "rank", range(1, len(ranked) + 1))
    ranked["chemistry"] = CHEMISTRY
    ranked_csv = args.out_dir / f"{GENE}_2moe_oligoai_ranked.csv"
    ranked[OUTPUT_COLS].to_csv(ranked_csv, index=False)
    print(f"Ranked {len(ranked)} candidates -> {ranked_csv}")
    print(ranked[OUTPUT_COLS[:4]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
