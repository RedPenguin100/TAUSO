# CLAUDE.md

## Environment

All commands must be run using the `tauso_claude` mamba environment. Mamba must be initialized first using this exact sequence — no other method:

```bash
export MAMBA_EXE="/home/michael/miniforge3/bin/mamba" && \
export MAMBA_ROOT_PREFIX="/home/michael/miniforge3" && \
source "/home/michael/miniforge3/etc/profile.d/mamba.sh" && \
mamba run -n tauso_claude <command>
```

Never use `python`, `python3`, `pytest`, or any other executable directly. Always initialize mamba as above and prefix commands with `mamba run -n tauso_claude`.

## Filesystem Access

Do not access `/home/michael` directly without explicit user permission. The mamba init lines above are the only exception — they are copied verbatim from the user-provided alias and must not be changed.

## Imports inside the package

Always use relative imports (e.g. `from ..hybridization.fast_hybridization import ...`, `from ...data.consts import ...`) when writing code inside the `src/tauso` package. Never use absolute `from tauso...` imports inside package files.
