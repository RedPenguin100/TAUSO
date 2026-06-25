"""Availability checks for the external competitor tools, used to skip integration tests."""
import shutil
import subprocess


def tool_missing(tool):
    """True if `tool` is not on PATH."""
    return shutil.which(tool) is None


def pfred_container_running(name="pfred"):
    """True if the named PFRED Docker container is up (lightweight check, no submodule import)."""
    try:
        r = subprocess.run(["docker", "inspect", "-f", "{{.State.Running}}", name], capture_output=True, text=True)
        return r.stdout.strip() == "true"
    except Exception:
        return False
