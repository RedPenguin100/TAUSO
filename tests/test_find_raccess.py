import pytest

import subprocess
from pathlib import Path

from tauso._raccess.core import find_raccess


def test_sanity():
    exe = find_raccess()
    p = Path(exe)
    print(f"raccess path: {p}")

    # basic sanity about the path
    assert p.is_file()
    assert p.is_absolute()

    # run with no args, just to see that it executes at all
    cp = subprocess.run(
        [exe],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    # it *should* fail because we didn't give seqfile/outfile, but that means it ran
    assert cp.returncode != 127  # 127 would be "command not found" if we used a shell
    assert "cannot open file" in cp.stderr or "NO FILE" in cp.stderr
