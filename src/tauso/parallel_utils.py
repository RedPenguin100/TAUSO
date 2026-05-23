import logging
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

logger = logging.getLogger(__name__)


def init_pandarallel(n_jobs, progress_bar=False, verbose=0, use_memory_fs=True):
    """
    Initialize pandarallel when n_jobs > 1.
    Returns True if parallel mode was activated, False for serial (n_jobs <= 1 or ImportError).
    """
    if n_jobs <= 1:
        return False
    try:
        from pandarallel import pandarallel

        pandarallel.initialize(
            nb_workers=n_jobs, progress_bar=progress_bar, verbose=verbose, use_memory_fs=use_memory_fs
        )
        return True
    except ImportError:
        logger.warning("pandarallel not found; falling back to single-core apply")
        return False


def make_apply_fn(obj, n_jobs, progress_bar=False, verbose=0, use_memory_fs=True):
    """
    Return obj.parallel_apply (pandarallel) or obj.apply (standard pandas).
    Works with Series, DataFrame, and GroupBy objects.
    Initializes pandarallel once; pass the returned callable into the hot loop.
    """
    use_parallel = init_pandarallel(n_jobs, progress_bar=progress_bar, verbose=verbose, use_memory_fs=use_memory_fs)
    return obj.parallel_apply if use_parallel else obj.apply


def run_parallel_tasks(fn, args_list, n_jobs, executor="thread"):
    """
    Run fn(*args) for each args in args_list, returning results in submission order.

    When n_jobs <= 1 or there is only one task, runs serially with no pool overhead.
    executor: "thread" (default) — good for subprocess-heavy work (GIL released).
              "process" — true multiprocessing; args/results must be picklable.
    """
    if not args_list:
        return []
    if n_jobs <= 1 or len(args_list) == 1:
        return [fn(*args) for args in args_list]

    workers = min(n_jobs, len(args_list))
    cls = ThreadPoolExecutor if executor == "thread" else ProcessPoolExecutor
    with cls(max_workers=workers) as pool:
        futures = [pool.submit(fn, *args) for args in args_list]
        return [f.result() for f in futures]


def run_keyed_parallel_tasks(tasks, fn, n_jobs, executor="thread"):
    """
    Run fn(*args) for each (key, *args) in tasks, returning {key: result}.

    Uses as_completed internally so faster tasks don't block slower ones.
    When n_jobs <= 1 or there is only one task, runs serially.
    executor: "thread" (default) or "process" (args/results must be picklable).

    Note: closures are not picklable — use executor="thread" for closure-based
    functions.
    """
    results: dict = {}
    if not tasks:
        return results
    if n_jobs <= 1 or len(tasks) == 1:
        for key, *args in tasks:
            results[key] = fn(*args)
        return results

    workers = min(n_jobs, len(tasks))
    cls = ThreadPoolExecutor if executor == "thread" else ProcessPoolExecutor
    with cls(max_workers=workers) as pool:
        future_to_key = {pool.submit(fn, *args): key for key, *args in tasks}
        for fut in as_completed(future_to_key):
            results[future_to_key[fut]] = fut.result()
    return results
