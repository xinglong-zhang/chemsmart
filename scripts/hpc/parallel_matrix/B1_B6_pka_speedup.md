# B1 vs B6 pKa batch wall-time speedup

Live comparison of serial (`%1`) vs parallel (`%5`) array submission for
`pka_scale.csv` (5 molecules) on CUHK CHPC.

| | **B1** (serial) | **B6** (parallel) |
|--|--|--|
| Run dir | `…/parallel_matrix/runs/B1_20260729_160355` | `…/parallel_matrix/runs/B6_20260729_172207` |
| Slurm job | `2011482` | `2012047` |
| Array | `1-5%1` | `1-5%5` |
| Throttle | M=1 | M=5 |
| Tasks | 5/5 COMPLETED (exit 0) | 5/5 COMPLETED (exit 0) |

## Per-task sacct

### B1 — serial (`1-5%1`)

| Task | State | Elapsed | TotalCPU | Node |
|------|-------|---------|----------|------|
| 2011482_1 | COMPLETED | 00:10:15 | 05:11:20 | chpc-cn073 |
| 2011482_2 | COMPLETED | 00:01:12 | 30:01 | chpc-cn073 |
| 2011482_3 | COMPLETED | 00:01:02 | 25:47 | chpc-cn073 |
| 2011482_4 | COMPLETED | 00:04:40 | 02:21:23 | chpc-cn073 |
| 2011482_5 | COMPLETED | 00:01:37 | 46:30 | chpc-cn073 |

Tasks ran back-to-back (16:04:25 → 16:23:12).

### B6 — parallel (`1-5%5`)

| Task | State | Elapsed | TotalCPU | Node |
|------|-------|---------|----------|------|
| 2012047_1 | COMPLETED | 00:09:17 | 04:37:12 | chpc-cn072 |
| 2012047_2 | COMPLETED | 00:01:19 | 28:33 | chpc-cn072 |
| 2012047_3 | COMPLETED | 00:00:55 | 23:50 | chpc-cn073 |
| 2012047_4 | COMPLETED | 00:04:16 | 02:09:05 | chpc-cn073 |
| 2012047_5 | COMPLETED | 00:01:30 | 43:09 | chpc-cn073 |

All tasks started together at 17:22:49.

## Speedup

| Metric | B1 serial | B6 parallel | B1 / B6 |
|--------|-----------|-------------|---------|
| **Array wall** (`max End − min Start`) | **1127 s (18.8 min)** | **557 s (9.3 min)** | **2.02×** |
| Σ task Elapsed | 1126 s | 1037 s | 1.09 |
| Σ TotalCPU (sacct) | 33301 s (9.25 CPU-h) | 30109 s (8.36 CPU-h) | 1.11 |
| Σ Gaussian `Job cpu time` | 33284 s (554.7 min) | 30092 s (501.5 min) | 1.11 |
| max task Elapsed | 615 s | 557 s | — |

- **Wall speedup ≈ 2.0×** — user wait time roughly halves under `%5`.
- Serial wall ≈ Σ Elapsed (efficiency ~1.0×) → classic `%1` chaining.
- Parallel wall ≈ max(task); Σ Elapsed / wall ≈ **1.86×** effective overlap.
- Chemistry cost is similar (TotalCPU / Gaussian CPU ratio ~1.1); small differences are run-to-run variance.

## Why not ~5×

Molecule **1a** dominates (~10 min vs ~1–4 min for others). The imbalance limit is
about Σ / max ≈ **1.8×**. Observing **2.0×** is at/near that bound: parallelization
works; the pKa CSV is uneven.

## Gaussian Job cpu time by molecule

| Molecule | B1 (s) | B6 (s) |
|----------|--------|--------|
| 1a | 18675 | 16629 |
| 2a | 1797 | 1709 |
| 3a | 1544 | 1426 |
| 4a | 8481 | 7741 |
| 5a | 2787 | 2587 |
| **sum** | **33284** | **30092** |
