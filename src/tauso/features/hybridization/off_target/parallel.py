"""Shared thread-pool dispatch for the RIsearch scoring paths.

Both the general/specific orchestration (populate_off_target) and the per-gene
orchestration (gene_chunk_scoring) submit their RIsearch tasks through here, so there is
a single dispatch implementation rather than one per caller.
"""

from concurrent.futures import ThreadPoolExecutor


def run_tasks_parallel(tasks, fn, n_jobs):
    """Run ``[(key, *args), ...]`` as ``fn(*args)`` on up to ``n_jobs`` threads.

    Returns ``{key: result}``. Runs serially when ``n_jobs == 1`` or there is a single task.
    """
    results = {}
    if n_jobs > 1 and len(tasks) > 1:
        with ThreadPoolExecutor(max_workers=min(n_jobs, len(tasks))) as pool:
            future_to_key = {pool.submit(fn, *args): key for key, *args in tasks}
            for fut, key in future_to_key.items():
                results[key] = fut.result()
    else:
        for key, *args in tasks:
            results[key] = fn(*args)
    return results
