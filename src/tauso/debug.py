"""Diagnostics for slow or memory-hungry runs. Nothing here is on a hot path."""

import functools
import gc
import logging
import os

import psutil

logger = logging.getLogger(__name__)


def log_dataframe_memory(df, label=""):
    """Log a DataFrame's total memory and its heaviest columns.

    The populate pipeline accumulates one feature per column on a single frame, so
    the per-column breakdown is what identifies a heavy feature family. `deep=True`
    is needed to size object columns rather than just their pointers.
    """
    usage = df.memory_usage(deep=True) / (1024**2)
    heaviest = ", ".join(f"{col} {mb:.1f}MB" for col, mb in usage.sort_values(ascending=False).head(5).items())
    logger.info("[MEM] %s | %.1fMB over %d columns | heaviest: %s", label, usage.sum(), df.shape[1], heaviest)


def log_memory_usage(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        process = psutil.Process(os.getpid())
        mem_before = process.memory_info().rss / (1024**2)

        # Force GC to get a clean baseline
        gc.collect()

        result = func(*args, **kwargs)

        gc.collect()
        mem_after = process.memory_info().rss / (1024**2)
        logger.debug(
            f"--- [MEM] {func.__name__} | Before: {mem_before:.1f}MB | After: {mem_after:.1f}MB | Change: {mem_after - mem_before:.1f}MB ---"
        )
        return result

    return wrapper
