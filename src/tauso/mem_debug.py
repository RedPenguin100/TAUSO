"""
Lightweight memory instrumentation for the off-target / RIsearch pipeline.

All probes are opt-in via TAUSO_MEM_DEBUG (set to "1" to enable). When
disabled the helpers are near-zero cost (just an env check + early return),
so it's safe to leave the calls in production code.

What we track:
  - Process-tree RSS: parent + all subprocess descendants (RIsearch is a
    separate process, so its memory doesn't show up in the parent's RSS).
  - Python heap peak (via tracemalloc) optionally, for narrowing the
    parent-side contribution.
  - Sizes of intermediate buffers (RIsearch output bytes, DataFrame
    memory_usage, dict entry count).

Log lines all carry [mem pid=… tid=… label=…] so they can be parsed
out of large log files. RSS values are reported in MB.
"""

from __future__ import annotations

import logging
import os
import threading
import time
import tracemalloc
from contextlib import contextmanager

import psutil

logger = logging.getLogger(__name__)


def _enabled() -> bool:
    return os.environ.get("TAUSO_MEM_DEBUG", "").lower() in ("1", "true", "yes", "on")


_MAIN_PROC = psutil.Process(os.getpid())


def total_rss_mb() -> float:
    """RSS of parent + all live descendant processes in MB.

    RIsearch runs as a subprocess so its working set is NOT counted in the
    Python parent's RSS. Walk children to get the real memory cost.
    """
    rss = 0
    try:
        rss += _MAIN_PROC.memory_info().rss
    except psutil.Error:
        return 0.0
    try:
        for child in _MAIN_PROC.children(recursive=True):
            try:
                rss += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
    except psutil.Error:
        pass
    return rss / 1024**2


def parent_rss_mb() -> float:
    """RSS of the Python parent process only (excludes child subprocesses)."""
    try:
        return _MAIN_PROC.memory_info().rss / 1024**2
    except psutil.Error:
        return 0.0


def log_mem(label: str, level: int = logging.INFO) -> None:
    """One-shot snapshot: parent RSS, tree RSS, optional tracemalloc peak."""
    if not _enabled():
        return
    parent = parent_rss_mb()
    tree = total_rss_mb()
    children = tree - parent
    extra = ""
    if tracemalloc.is_tracing():
        cur, peak = tracemalloc.get_traced_memory()
        extra = f" heap={cur / 1024**2:.0f}M peak={peak / 1024**2:.0f}M"
    logger.log(
        level,
        "[mem pid=%d tid=%d] %s | parent=%.0fM tree=%.0fM children=%.0fM%s",
        os.getpid(),
        threading.get_native_id(),
        label,
        parent,
        tree,
        children,
        extra,
    )


@contextmanager
def track(label: str, level: int = logging.INFO):
    """Context manager: logs delta RSS and elapsed time over the block.

    Use around the suspect hot spot:
        with track("risearch_call"):
            result = subprocess.check_output(...)
    """
    if not _enabled():
        yield
        return
    parent_before = parent_rss_mb()
    tree_before = total_rss_mb()
    t0 = time.perf_counter()
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        parent_after = parent_rss_mb()
        tree_after = total_rss_mb()
        logger.log(
            level,
            "[mem pid=%d tid=%d] %s | %.2fs | parent %.0f→%.0f (Δ%+.0f)M "
            "tree %.0f→%.0f (Δ%+.0f)M",
            os.getpid(),
            threading.get_native_id(),
            label,
            dt,
            parent_before,
            parent_after,
            parent_after - parent_before,
            tree_before,
            tree_after,
            tree_after - tree_before,
        )


def log_obj_size(label: str, obj, level: int = logging.INFO) -> None:
    """Best-effort size report for an intermediate buffer / DataFrame / dict."""
    if not _enabled():
        return
    size_mb: float = 0.0
    extra = ""
    try:
        if isinstance(obj, (str, bytes, bytearray)):
            size_mb = len(obj) / 1024**2
            extra = f" len={len(obj):,}"
        elif hasattr(obj, "memory_usage"):
            mu = obj.memory_usage(deep=True)
            size_mb = float(mu.sum() if hasattr(mu, "sum") else mu) / 1024**2
            if hasattr(obj, "shape"):
                extra = f" shape={obj.shape}"
        elif isinstance(obj, dict):
            size_mb = len(obj) * 64 / 1024**2  # rough: 64 bytes per entry
            extra = f" n_entries={len(obj):,} (estimate)"
        elif isinstance(obj, (list, tuple)):
            size_mb = len(obj) * 8 / 1024**2  # rough
            extra = f" n_entries={len(obj):,} (estimate)"
    except Exception as e:
        extra = f" (sizing failed: {e})"
    logger.log(
        level,
        "[mem pid=%d tid=%d] obj %s | %.1fM%s",
        os.getpid(),
        threading.get_native_id(),
        label,
        size_mb,
        extra,
    )


def start_tracemalloc() -> None:
    """Optional: enable Python heap tracking. Adds overhead — only when needed."""
    if _enabled() and not tracemalloc.is_tracing():
        tracemalloc.start(25)


def dump_top_allocators(label: str, n: int = 10, level: int = logging.INFO) -> None:
    """Dump the top-N tracemalloc allocation sites by current size."""
    if not _enabled() or not tracemalloc.is_tracing():
        return
    snapshot = tracemalloc.take_snapshot()
    stats = snapshot.statistics("lineno")[:n]
    logger.log(level, "[mem pid=%d tid=%d] top allocators @ %s:", os.getpid(), threading.get_native_id(), label)
    for s in stats:
        logger.log(
            level,
            "    %.1fM | %s",
            s.size / 1024**2,
            "; ".join(str(f) for f in s.traceback.format()[-3:]),
        )
