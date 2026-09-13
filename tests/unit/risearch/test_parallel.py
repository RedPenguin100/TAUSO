"""run_tasks_parallel: results keyed by task, and a failure that does not wait for the queue."""

import time

import pytest

from tauso.features.hybridization.off_target.parallel import run_tasks_parallel


def double(x):
    return x * 2


@pytest.mark.parametrize("n_jobs", [1, 4])
def test_results_are_keyed_by_task(n_jobs):
    tasks = [("a", 1), ("b", 2), ("c", 3)]

    assert run_tasks_parallel(tasks, double, n_jobs) == {"a": 2, "b": 4, "c": 6}


def test_a_failure_is_raised_without_running_the_rest():
    started = []

    def work(i):
        started.append(i)
        if i == 0:
            raise RuntimeError("bad gene")
        time.sleep(0.05)
        return i

    with pytest.raises(RuntimeError, match="bad gene"):
        run_tasks_parallel([(i, i) for i in range(50)], work, n_jobs=2)

    assert len(started) < 50
