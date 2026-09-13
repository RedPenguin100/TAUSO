"""Shared thread-pool dispatch for the RIsearch scoring paths.

Both the general/specific orchestration (populate_off_target) and the per-gene
orchestration (gene_chunk_scoring) submit their RIsearch tasks through here, so there is
a single dispatch implementation rather than one per caller.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed


def run_tasks_parallel(tasks, fn, n_jobs):
    """Run ``[(key, *args), ...]`` as ``fn(*args)`` on up to ``n_jobs`` threads.

    Returns ``{key: result}``. Runs serially when ``n_jobs == 1`` or there is a single task.

    A task that fails is raised as soon as it fails: the tasks not yet started are
    cancelled and only the ones already running are waited for, so a bad gene early in
    a long scan is reported early.
    """
    results = {}
    if n_jobs > 1 and len(tasks) > 1:
        pool = ThreadPoolExecutor(max_workers=min(n_jobs, len(tasks)))
        try:
            key_of = {pool.submit(fn, *args): key for key, *args in tasks}
            for fut in as_completed(key_of):
                results[key_of[fut]] = fut.result()
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
    else:
        for key, *args in tasks:
            results[key] = fn(*args)
    return results
