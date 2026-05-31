# Running on GPU

```bash
USE_GPU=1 ./amber_pipe_new_KLKB1_v3.sh <pdb> <n_threads> <nodes>
```

- Runs one `pmemd.cuda` process per GPU, pinned to `CUDA_VISIBLE_DEVICES=0`.
  `n_threads`/`nodes` apply only to the CPU path.
- Plain `pmemd.cuda` is the `_SPFP` build (fast). Do not use `_DPFP`.
- NPT stages require `barostat=2` (Monte-Carlo); `barostat=1` (Berendsen)
  aborts `pmemd.cuda` at step 1 on these OPC systems.
- On the cluster, build/binary is per-machine: use the cluster's own
  `pmemd.cuda` (the local one is compiled for the RTX 5090's arch), submit to
  a GPU partition with `--gres=gpu:<model>:1 --ntasks=1`, and set `USE_GPU=1`.

Throughput (KLKB1, ~55.7k atoms): ~346 ns/day on one RTX 5090 vs ~3.4 ns/day
on 32 CPU ranks.
