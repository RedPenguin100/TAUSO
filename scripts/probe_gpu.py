"""
GPU probe script — run this as a SLURM job to discover what GPU hardware
is actually visible and whether XGBoost + cupy work on it.

Usage (after finding a GPU partition via `sinfo -o "%P %G"`):
    sbatch -p <gpu-partition> --gres=gpu:1 --wrap="
        source /tamir2/kovaliov/miniconda3/etc/profile.d/conda.sh
        conda activate tauso
        python /tamir2/kovaliov/working_directory/aso_prediction/TAUSO/scripts/probe_gpu.py
    "

Or interactively:
    srun -p <gpu-partition> --gres=gpu:1 --pty bash
    python scripts/probe_gpu.py
"""
import subprocess
import sys


def section(title):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


# ---------------------------------------------------------------------------
# 1. nvidia-smi: raw hardware view
# ---------------------------------------------------------------------------
section("nvidia-smi")
try:
    out = subprocess.check_output(["nvidia-smi"], text=True)
    print(out)
except FileNotFoundError:
    print("nvidia-smi not found — no NVIDIA driver visible from this node.")
    sys.exit(1)
except subprocess.CalledProcessError as e:
    print(f"nvidia-smi failed: {e}")
    sys.exit(1)

# GPU list in parseable form
section("GPU list (nvidia-smi -L)")
try:
    out = subprocess.check_output(["nvidia-smi", "-L"], text=True)
    print(out)
except Exception as e:
    print(f"Failed: {e}")

# CUDA_VISIBLE_DEVICES (SLURM sets this to restrict which GPUs a job sees)
import os
section("CUDA_VISIBLE_DEVICES")
print(os.environ.get("CUDA_VISIBLE_DEVICES", "(not set — all GPUs visible)"))

# ---------------------------------------------------------------------------
# 2. cupy
# ---------------------------------------------------------------------------
section("cupy")
try:
    import cupy as cp
    print(f"cupy version     : {cp.__version__}")
    print(f"CUDA runtime     : {cp.cuda.runtime.runtimeGetVersion()}")
    n = cp.cuda.runtime.getDeviceCount()
    print(f"Devices visible  : {n}")
    for i in range(n):
        props = cp.cuda.runtime.getDeviceProperties(i)
        name = props["name"].decode() if isinstance(props["name"], bytes) else props["name"]
        mem_gb = props["totalGlobalMem"] / 1024**3
        print(f"  [{i}] {name}  ({mem_gb:.1f} GB)")
    # Smoke test: small allocation + matmul
    a = cp.random.rand(1000, 1000, dtype=cp.float32)
    b = cp.random.rand(1000, 1000, dtype=cp.float32)
    _ = cp.dot(a, b)
    cp.cuda.Stream.null.synchronize()
    print("cupy matmul smoke test: PASSED")
except ImportError:
    print("cupy not installed — GPU arrays will not be used (XGBoost will still use GPU internally).")
except Exception as e:
    print(f"cupy error: {e}")

# ---------------------------------------------------------------------------
# 3. XGBoost GPU
# ---------------------------------------------------------------------------
section("XGBoost GPU")
try:
    import numpy as np
    import xgboost as xgb
    print(f"xgboost version  : {xgb.__version__}")

    rng = np.random.default_rng(0)
    X = rng.standard_normal((5000, 20)).astype(np.float32)
    y = rng.standard_normal(5000).astype(np.float32)

    dtrain = xgb.DMatrix(X, label=y)
    params = {"tree_method": "hist", "device": "cuda", "objective": "reg:squarederror", "seed": 0}
    try:
        bst = xgb.train(params, dtrain, num_boost_round=10, verbose_eval=False)
        preds = bst.predict(dtrain)
        print(f"XGBoost GPU train: PASSED  (pred mean={preds.mean():.4f})")
    except Exception as e:
        print(f"XGBoost GPU train: FAILED\n  {e}")

    # Also test QuantileDMatrix (used in RFE)
    try:
        dq = xgb.QuantileDMatrix(X, label=y)
        bst2 = xgb.train(params, dq, num_boost_round=10, verbose_eval=False)
        print("XGBoost QuantileDMatrix + GPU: PASSED")
    except Exception as e:
        print(f"XGBoost QuantileDMatrix + GPU: FAILED\n  {e}")

except Exception as e:
    print(f"XGBoost import/setup error: {e}")

# ---------------------------------------------------------------------------
# 4. Summary
# ---------------------------------------------------------------------------
section("Summary")
print("If all sections above show PASSED, you can run train_model.py with --device cuda.")
print("If cupy is missing but XGBoost GPU passed, GPU tree building works but without cupy array speedup.")
print("If XGBoost GPU failed, check CUDA version compatibility with your xgboost build.")
