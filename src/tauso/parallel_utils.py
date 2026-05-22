import logging

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

        pandarallel.initialize(nb_workers=n_jobs, progress_bar=progress_bar, verbose=verbose, use_memory_fs=use_memory_fs)
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
