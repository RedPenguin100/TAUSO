# src/tauso/_raccess/cli.py
from importlib.resources import files
import subprocess
import sys


def main():
    script = files("tauso._raccess") / "install_raccess.sh"

    # forward ALL user args to the shell script
    cmd = ["bash", str(script), *sys.argv[1:]]
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)
