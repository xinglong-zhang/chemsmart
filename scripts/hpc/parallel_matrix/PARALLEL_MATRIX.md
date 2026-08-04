# Parallel matrix HPC tests

Scripts live in:

```text
/project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix/
```

They exercise `chemsmart sub` batch / nestable parallel options on CUHK CHPC.

## Prerequisites

1. Work on **sandbox** (login node kills `chemsmart`):

   ```bash
   ssh sandbox
   source ~/miniconda3/etc/profile.d/conda.sh && conda activate chemsmart
   export CHEMSMART_CONFIG_DIR=~/.chemsmart
   cd /project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix
   chmod +x *.sh
   ```

2. Inputs (already expected under project paths; each run dir only links the
   files for that suite):

   | Suite | Input |
   |-------|-------------|
   | **B*** / **E2–E5** | `$PROJECT/chemsmart/test_pka/gaussian/pka_scale.csv` (+ `1a.xyz`…`5a.xyz`) |
   | **E1** | `$PROJECT/chemsmart/test_batch/water.com` |
   | **N*** | `$PROJECT/chemsmart/test_batch/model_R_c5_ts_irc_flatr_flat_opt.log` (DIAS) |
   | **minimum/** | `parallel_matrix/minimum/inputs/` (tiny xyz/com/csv + freq logs; DIAS still from test_batch) |

3. Default mode is **dry-run** (`DRY_RUN=1` → `chemsmart sub --test`). No jobs are queued unless you set `DRY_RUN=0`.

## Quick start

```bash
# Script-generation checks only (recommended first)
DRY_RUN=1 ./run_all.sh

# One suite
DRY_RUN=1 ./run_B.sh
DRY_RUN=1 ./run_E.sh
DRY_RUN=1 ./run_N_dias.sh

# Subset of cases
CASES=B1,B6 DRY_RUN=1 ./run_B.sh
CASES=N1,N2 DRY_RUN=1 ./run_N_dias.sh

# Real submit + wait (uses cluster time; waits between cases)
DRY_RUN=0 WAIT=1 CASES=B6 ./run_B.sh
DRY_RUN=0 WAIT=1 CASES=B1 ./run_B.sh
DRY_RUN=0 WAIT=1 CASES=E3,E4,E5 ./run_E.sh
```

Results:

- Per-case workdirs: `parallel_matrix/runs/<CASE>_<timestamp>/`
- Logs / summaries: `parallel_matrix/logs/`

## Case catalogue

### B-series — BatchJob (pKa CSV, 5 rows)

| Case | Flags | Expected `#SBATCH --array` |
|------|--------|----------------------------|
| **B1** | `--no-run-in-parallel` (default serial) | `1-5%1` |
| **B2** | `--no-run-in-parallel -M 5` | still `1-5%1` (`-M` ignored when serial) |
| **B3** | `--run-in-parallel` | `1-5` or `1-5%<cap>` |
| **B4** | `--run-in-parallel -M 1` | `1-5%1` |
| **B5** | `--run-in-parallel -M 2` | `1-5%2` |
| **B6** | `--run-in-parallel -M 5` | `1-5%5` |
| **B7** | `--run-in-parallel -M 10` | `1-5%5` (clamped to N) |

Also checks that only **one** shared `chemsmart_run_array_<array_log_stem>.py`
is written (not legacy `chemsmart_run_array.py` or per-task `*_1.py`…`*_5.py`).
Shared stdout is `<array_log_stem>.out` (with matching `.err` / `.loglock`).

Equivalent manual command (B6):

```bash
cd /project/xlzhang/jingyi/chemsmart/test_pka/gaussian
chemsmart sub -s cu_batch --test --run-in-parallel -M 5 \
  gaussian -p pka_wat -f pka_scale.csv \
  pka --scheme direct batch -R
```

### E-series — controls / edges

| Case | What it checks | Notes |
|------|----------------|-------|
| **E1** | Single-molecule `sp` → non-array submit | baseline |
| **E2** | `-R` rewrites shared runscript | dry-run OK |
| **E3** | Wall time B1 vs B6 | needs `DRY_RUN=0` |
| **E4** | `squeue` shows array tasks | needs `DRY_RUN=0` |
| **E5** | flock / headed shared log / `.loglock` | dry-run checks script wiring; live log needs real run |

### N-series — nestable DIAS (TS mode → 3 children)

Uses:

```bash
dias -i 1-10 -n 20 --mode ts
# on model_R_c5_ts_irc_flatr_flat_opt.log (defaults: parent 2/2, f1 0/1, f2 2/2)
```

| Case | Flags | Expected |
|------|--------|----------|
| **N1** | `--no-run-in-parallel` | **one** non-array `chemsmart_sub_ts*.sh` |
| **N2** | `--run-in-parallel` | array `1-3` (or `%M`), shared runscript with `--child-index` |
| **N3** | `--run-in-parallel -M 1` | `1-3%1` |
| **N4** | `--run-in-parallel -M 2` | `1-3%2` |

Equivalent manual commands:

```bash
cd /project/xlzhang/jingyi/chemsmart/test_batch

# N1 — single parent
chemsmart sub -s cu_batch --test --no-run-in-parallel gaussian -p test \
  -f model_R_c5_ts_irc_flatr_flat_opt.log --charge 2 --multiplicity 2 \
  dias -i 1-10 -n 20 --mode ts -c1 0 -m1 1 -c2 2 -m2 2 -R

# N4 — array of 3 children, 2 concurrent
chemsmart sub -s cu_batch --test --run-in-parallel -M 2 gaussian -p test \
  -f model_R_c5_ts_irc_flatr_flat_opt.log --charge 2 --multiplicity 2 \
  dias -i 1-10 -n 20 --mode ts -c1 0 -m1 1 -c2 2 -m2 2 -R
```

## Environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `DRY_RUN` | `1` | `1` → `--test` (no sbatch) |
| `WAIT` | `0` | `1` → after each live case, wait for its job(s) to leave the queue (B, E, N) |
| `SERVER` | `cu_batch` | `-s` server yaml |
| `GAUSSIAN_PROJECT_BATCH` | `pka_wat` | project for pKa batch |
| `GAUSSIAN_PROJECT_OPT` | `test` | project for E1 |
| `GAUSSIAN_PROJECT_DIAS` | `test` | project for N-series |
| `CASES` | suite default | e.g. `B1,B6` |
| `DIAS_FRAGMENTS` | `1-10` | DIAS `-i` |
| `DIAS_EVERY_N` | `20` | DIAS `-n` |
| `DIAS_MODE` | `ts` | `ts` keeps child count at 3 |
| `FORCE_RERUN` | `1` | pass `-R` |

## Interpreting flags

- **Default / `--no-run-in-parallel`**: BatchJob arrays use `%1`; nestable jobs stay one parent.
- **`--run-in-parallel`**: enables concurrent array tasks; nestable parents expand to arrays.
- **`-M` / `--max-tasks`**: only meaningful with `--run-in-parallel` (ignored under serial policy — see B2).
- Nestable dry-runs use colliding labels (`-l ts` / `-l irc` / `-l opt`) and assert `--child-index`
  appears **after** the job token in `TASK_CLI` (guards CLI rewrite anchoring).

## Suggested order

1. `DRY_RUN=1 ./run_all.sh` — script / array-line checks  
2. `DRY_RUN=1 CASES=N1,N2,N3,N4 ./run_N_dias.sh`  
3. `DRY_RUN=0 WAIT=1 CASES=B6 ./run_B.sh` then `CASES=B1` — live batch (one case at a time; `WAIT=1` now waits between cases)  
4. `DRY_RUN=0 WAIT=1 CASES=E3,E4,E5 ./run_E.sh` — queue + shared-log checks  
5. `DRY_RUN=1 ./minimum/run_minimum.sh` — covering map over all top-level job types  

Live B/N cases share batch labels (`pka_scale_pka_batch`, `ts`). With `WAIT=1`, later cases wait for earlier jobs to leave the queue so the duplicate-submit guard and shared scratch stay consistent. Prefer one live case per invocation when comparing serial vs parallel wall time.

## Minimum set (cases × job types)

See [`minimum/README.md`](minimum/README.md).

Covers **B1, B2, B6, O2, N1, N2, E3** and all top-level Gaussian (18) + ORCA (10)
commands via a **28-row covering map** (one case per job type), using tiny
inputs under `minimum/inputs/`. Default `DRY_RUN=1`. Optional live gate:

```bash
cd minimum
DRY_RUN=1 ./run_minimum.sh
DRY_RUN=0 WAIT=1 CASES=E3 ./run_minimum.sh   # five_mols SP serial vs parallel wall time
```

## Cleanup

```bash
rm -rf /project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix/runs
```
