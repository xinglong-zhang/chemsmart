# Crest / QRC nestable wall-time speedup

Goal: compare serial nestable parents (`--no-run-in-parallel`) vs parallel
array children (`--run-in-parallel -M …`) for **crest** and **qrc**, same style
as [`B1_B6_pka_speedup.md`](B1_B6_pka_speedup.md).

## Status (after live speedup runs)

| Job | Serial live | Parallel live | Fair speedup pair? | Wall speedup |
|-----|-------------|---------------|--------------------|--------------|
| **ORCA qrc** | yes (`2012440`) | yes (`2014355`, `1-2%2`) | **Yes** | **2.09×** |
| **Gaussian qrc** | yes (`2014502`) | yes (`2014514`, `1-2%2`) | **Yes** | **1.75×** (noisy) |
| **Gaussian crest** | yes (`2014625`) | yes (`2014626`, `1-5%5`) | **Yes** | **3.6×** (too short; c4/c5 fail) |
| **dias** (nestable proxy) | yes (parent N1) | yes (parent N2/N4) | **Yes** | **1.95×** |

Run roots:

- ORCA qrc serial: `…/runs/speedup_nestable_20260729_180037/orca_qrc_serial`
- All other pairs: `…/runs/speedup_nestable_20260729_214808/{orca_qrc_par,g16_qrc_serial,g16_qrc_par,crest_serial,crest_par}`

---

## Live matched pairs

### 1. ORCA qrc (recommended — minutes of work)

| | **Serial** | **Parallel** |
|--|--|--|
| Run dir | `…/speedup_nestable_20260729_180037/orca_qrc_serial` | `…/speedup_nestable_20260729_214808/orca_qrc_par` |
| Slurm job | `2012440` (non-array) | `2014355` (`1-2%2`) |
| Label | `speed_orca_qrc_serial` | `speed_orca_qrc_par` |
| Input | `minimum/inputs/orca_freq.out` | same |
| State | COMPLETED (exit 0) | 2/2 COMPLETED (exit 0) |

**Serial** — one parent, 2 children (`f_opt`, `r_opt`) chained in-process.

| JobID | State | Elapsed | TotalCPU | Node |
|-------|-------|---------|----------|------|
| 2012440 | COMPLETED | 00:14:41 | 05:24:54 | chpc-cn078 |

**Parallel** — array tasks start together at 21:51:01.

| Task | State | Elapsed | TotalCPU | Node |
|------|-------|---------|----------|------|
| 2014355_1 | COMPLETED | 00:05:36 | 01:57:52 | chpc-cn073 |
| 2014355_2 | COMPLETED | 00:07:02 | 01:58:55 | chpc-cn074 |

**Speedup**

| Metric | Serial | Parallel | Serial / Parallel |
|--------|--------|----------|-------------------|
| **Job / array wall** | **881 s (14.7 min)** | **422 s (7.0 min)** | **2.09×** |
| Σ task Elapsed | 881 s | 758 s | 1.16 |
| Σ TotalCPU | 19494 s | 14207 s | 1.37 |
| max task Elapsed | — | 422 s | — |

- Wall speedup ≈ **2.1×** — consistent with B1/B6 pKa (**2.02×**) and DIAS N1/N2 (**1.95×**).
- Two children of similar length (~6–7 min each); parallel wall ≈ max(child).
- Chemistry cost ratio Σ TotalCPU ≈ 1.37 reflects run-to-run variance, not lost parallelism.

---

### 2. Gaussian qrc (minimum fixture — ~30 s)

| | **Serial** | **Parallel** |
|--|--|--|
| Run dir | `…/speedup_nestable_20260729_214808/g16_qrc_serial` | `…/g16_qrc_par` |
| Slurm job | `2014502` (non-array) | `2014514` (`1-2%2`) |
| Label | `speed_g16_qrc_serial` | `speed_g16_qrc_par` |
| Input | `minimum/inputs/ts_freq.log` | same |
| State | COMPLETED (exit 0) | 2/2 COMPLETED (exit 0) |

| Task | Serial | Parallel |
|------|--------|----------|
| JobID | 2014502 | 2014514_1, 2014514_2 |
| Elapsed | 00:00:28 | 00:00:15, 00:00:16 |
| TotalCPU | 09:59 | 05:00, 04:25 |

**Speedup**

| Metric | Serial | Parallel | Serial / Parallel |
|--------|--------|----------|-------------------|
| **Job / array wall** | **28 s** | **16 s** | **1.75×** |
| Σ task Elapsed | 28 s | 31 s | 0.90 |
| Σ TotalCPU | 599 s | 566 s | 1.06 |

- Wiring is correct (non-array serial vs `1-2%2` parallel, both exit 0).
- Wall ratio is **not meaningful** at this scale: Slurm/Gaussian startup dominates.
- Confirms the earlier min N2 estimate (~2× from Σ Elapsed) was optimistic for real serial parents.

---

### 3. Gaussian crest (`five_mols.xyz` — ~30 s total)

| | **Serial** | **Parallel** |
|--|--|--|
| Run dir | `…/speedup_nestable_20260729_214808/crest_serial` | `…/crest_par` |
| Slurm job | `2014625` (non-array) | `2014626` (`1-5%5`) |
| Label | `speed_crest_serial_opt` | `speed_crest_par_opt_array` |
| Input | `minimum/inputs/five_mols.xyz` | same |
| State | COMPLETED (exit 0) | **3/5 COMPLETED**, 2/5 FAILED (exit 1) |

| Task | Serial (in-process) | Parallel |
|------|---------------------|----------|
| c1–c3 | chained, exit 0 | 2014626_1–3 COMPLETED (8–9 s each) |
| c4–c5 | chained, Gaussian L701 error | 2014626_4–5 FAILED (L701 error) |

**Speedup**

| Metric | Serial | Parallel | Serial / Parallel |
|--------|--------|----------|-------------------|
| **Job / array wall** | **32 s** | **9 s** | **3.56×** |
| Σ task Elapsed | 32 s | 44 s | 0.73 |
| Σ TotalCPU | 649 s | 542 s | 1.20 |

- Parallel overlap works (5 tasks on 3 nodes from 22:25:54).
- **Do not treat 3.6× as a production speedup**: total wall is ~30 s; c4/c5 fail with the same
  Gaussian `L701.exe` error in **both** serial and parallel (fixture/chemistry issue, not parallelization).
- Parent serial job reports exit 0 despite c4/c5 errors; array tasks 4/5 correctly report FAILED.

---

## Best existing nestable pair: parent DIAS (N1 vs N2)

Same expand-to-array path crest/qrc use (`--child-index` + shared runscript).

| | **N1** (serial parent) | **N2** (parallel `%3`) |
|--|--|--|
| Run dir | `…/runs/N1_20260728_174622` | `…/runs/N2_20260728_213805` |
| Slurm job | `2007586` | `2009622` |
| Submit | non-array `chemsmart_sub_ts.sh` | `#SBATCH --array=1-3%3` |
| Tasks | 1 parent, 3 children in-process | 3 array tasks |
| State | COMPLETED | 3/3 COMPLETED |

### Per-task / parent sacct

**N1 serial**

| JobID | State | Elapsed | TotalCPU |
|-------|-------|---------|----------|
| 2007586 | COMPLETED | 00:40:25 | 21:10:22 |

**N2 parallel**

| JobID | State | Elapsed | TotalCPU | Node |
|-------|-------|---------|----------|------|
| 2009622_1 | COMPLETED | 00:20:46 | 10:52:20 | chpc-cn082 |
| 2009622_2 | COMPLETED | 00:00:49 | 11:37 | chpc-cn082 |
| 2009622_3 | COMPLETED | 00:16:21 | 08:31:19 | chpc-cn072 |

### Speedup

| Metric | N1 serial | N2 parallel | N1 / N2 |
|--------|-----------|-------------|---------|
| **Array / job wall** | **2425 s (40.4 min)** | **1246 s (20.8 min)** | **1.95×** |
| Σ task Elapsed | 2425 s | 2276 s | 1.07 |
| Σ TotalCPU | 76222 s | 70516 s | 1.08 |
| max child Elapsed | — | 1246 s | — |

- Wall speedup ≈ **2.0×**.
- Parallel efficiency Σ/wall ≈ **1.83×** (3 children; one is ~49 s, two are ~16–21 min).
- N3 (`1-3%1`) wall ≈ 2241 s ≈ serial (array still one-at-a-time).
- N4 (`1-3%2`) wall ≈ 1240 s ≈ N2 (**1.96×** vs N1).

### Gaussian child elapsed (from logs)

| Child | N1 elapsed (s) | N2 elapsed (s) |
|-------|----------------|----------------|
| `ts_p1` | 1320 | 1234 |
| `ts_p1_f1` | 37 | 37 |
| `ts_p1_f2` | 1063 | 968 |

---

## Minimum matrix runs (pre-speedup estimates)

### Minimum crest (serial only) — superseded by live pair above

| | **min N1 gaussian crest** |
|--|--|
| Run dir | `…/runs/min_N1_gaussian_crest_20260729_140331` |
| Job | `2011388` (non-array) |
| Wall | **27 s** |

Per-child Gaussian elapsed ~4 s each; c4/c5 did not fully Normal-terminate.
Live crest pair confirms the same c4/c5 L701 failures.

### Minimum qrc (parallel only) — estimate vs live

**Gaussian qrc** — live serial/parallel pair above replaces the Σ-Elapsed estimate.

**ORCA qrc** — min N2 parallel-only estimate was **~1.74×** (wall 539 s, Σ Elapsed 938 s).
Live matched pair gives **2.09×** (881 s vs 422 s); the estimate was conservative.

| | **min N2 orca qrc** |
|--|--|
| Run dir | `…/runs/min_N2_orca_qrc_20260729_140354` |
| Job | `2011392` |
| Parallel wall | 539 s |
| Σ task Elapsed | 938 s |
| Estimated speedup (no serial job) | ~1.74× |

---

## One-off sandbox commands (copy-paste)

Use **distinct labels** per side so the duplicate-submit guard and shared scratch
do not collide. Run **serial first**, wait until it leaves the queue, then
parallel. Prefer **ORCA qrc** (minutes of work); gaussian qrc / crest on minimum
fixtures are short.

```bash
ssh sandbox
source ~/miniconda3/etc/profile.d/conda.sh && conda activate chemsmart
export CHEMSMART_CONFIG_DIR=~/.chemsmart

MIN=/project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix/minimum/inputs
SPEED=/project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix/runs/speedup_nestable_$(date +%Y%m%d_%H%M%S)
mkdir -p "$SPEED"/{orca_qrc_serial,orca_qrc_par,g16_qrc_serial,g16_qrc_par,crest_serial,crest_par}
```

### 1. ORCA qrc (recommended)

```bash
# --- serial: one parent, 2 children chained ---
cd "$SPEED/orca_qrc_serial"
ln -sfn "$MIN/orca_freq.out" .
chemsmart sub -s cu_batch --no-run-in-parallel -R \
  orca -p test -f orca_freq.out -c 0 -m 1 -l speed_orca_qrc_serial \
  qrc

# wait until job leaves queue, then:

# --- parallel: array 1-2%2 ---
cd "$SPEED/orca_qrc_par"
ln -sfn "$MIN/orca_freq.out" .
chemsmart sub -s cu_batch --run-in-parallel -M 2 -R \
  orca -p test -f orca_freq.out -c 0 -m 1 -l speed_orca_qrc_par \
  qrc
```

Expect: serial non-array submit; parallel `#SBATCH --array=1-2%2`.

### 2. Gaussian qrc

```bash
cd "$SPEED/g16_qrc_serial"
ln -sfn "$MIN/ts_freq.log" .
chemsmart sub -s cu_batch --no-run-in-parallel -R \
  gaussian -p test -f ts_freq.log -c 0 -m 1 -l speed_g16_qrc_serial \
  qrc

cd "$SPEED/g16_qrc_par"
ln -sfn "$MIN/ts_freq.log" .
chemsmart sub -s cu_batch --run-in-parallel -M 2 -R \
  gaussian -p test -f ts_freq.log -c 0 -m 1 -l speed_g16_qrc_par \
  qrc
```

Note: minimum `ts_freq.log` finishes in ~30 s/child — fine for wiring, weak for wall speedup.

### 3. Gaussian crest

```bash
cd "$SPEED/crest_serial"
ln -sfn "$MIN/five_mols.xyz" .
chemsmart sub -s cu_batch --no-run-in-parallel -R \
  gaussian -p test -f five_mols.xyz -c 0 -m 1 -l speed_crest_serial \
  crest -j opt

cd "$SPEED/crest_par"
ln -sfn "$MIN/five_mols.xyz" .
chemsmart sub -s cu_batch --run-in-parallel -M 5 -R \
  gaussian -p test -f five_mols.xyz -c 0 -m 1 -l speed_crest_par \
  crest -j opt
```

Note: `five_mols` crest is ~30 s total — use a larger xyz/log if you want a
meaningful crest wall ratio. `-M 5` matches five crest children from this input.

### After jobs finish

```bash
sacct -j 2012440,2014355,2014502,2014514,2014625,2014626 --parsable2 \
  -o JobID,State,ExitCode,Elapsed,TotalCPU,Start,End
```

Compare **array/job wall** (`max End − min Start`) and Σ TotalCPU as in
[`B1_B6_pka_speedup.md`](B1_B6_pka_speedup.md).

---

## Summary

| Comparison | Wall speedup | Confidence | Notes |
|------------|--------------|------------|-------|
| **ORCA qrc serial vs parallel** | **2.09×** | **High** | Matched live pair; ~15 min chemistry |
| DIAS N1 vs N2 (nestable proxy) | **1.95×** | High | Matched live pair; longer DIAS path |
| B1 vs B6 pKa (batch array) | **2.02×** | High | Related batch mechanism |
| Gaussian qrc serial vs parallel | **1.75×** | Low | Correct wiring; ~30 s jobs, overhead-dominated |
| Gaussian crest serial vs parallel | **3.56×** | Low | ~30 s jobs; c4/c5 L701 fail in both modes |
| ORCA qrc parallel-only estimate (min N2) | ~1.74× | Superseded | Live serial run confirms ~2× |

**Conclusion:** Nestable parallelization delivers roughly **2× wall-time speedup** when
child jobs run for several minutes (ORCA qrc, DIAS, pKa batch). Short Gaussian fixtures
(qrc, crest on `five_mols`) validate submit/array wiring but are too fast for reliable
wall ratios. Crest c4/c5 need a better input or method fix (Gaussian L701) before using
crest for speedup benchmarking.
