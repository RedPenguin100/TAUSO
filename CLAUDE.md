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

## pandas groupby.apply and pandarallel

`include_groups=False` was added to `groupby.apply` in pandas 2.x to suppress a FutureWarning about grouping columns being passed to the function. **pandarallel's `parallel_apply` does not accept this kwarg** and will raise `TypeError` at runtime.

When using `groupby.apply` alongside `pandarallel.parallel_apply`:
1. Check the pandas version before using `include_groups`: `pd.__version__ >= "2.0"`.
2. Or — preferred — drop the grouping column from the DataFrame before the groupby so neither path needs `include_groups`:
   ```python
   groups = df["group_col"]
   df_no_group = df.drop(columns=["group_col"])
   df_no_group.groupby(groups).apply(fn)          # pandas, no include_groups needed
   df_no_group.groupby(groups).parallel_apply(fn) # pandarallel, same
   ```
