# Minimum parallel matrix

Covering suite for batch / nestable parallel behavior across **all** top-level
Gaussian (18) and ORCA (10) job types, using the seven matrix cases
**B1, B2, B6, O2, N1, N2, E3**.

Not a full cartesian product: [`JOB_CASE_MAP.tsv`](JOB_CASE_MAP.tsv) assigns
**exactly one case per job type** (28 rows).

## Layout

```text
minimum/
  inputs/           # tiny fixtures (xyz/com/inp/csv + small freq logs)
  JOB_CASE_MAP.tsv  # covering assignment
  run_minimum.sh    # driver
  README.md
```

On HPC after sync:

```text
$PROJECT/chemsmart/test_batch/parallel_matrix/minimum/
```

## Prerequisites

Same as the parent suite ([`../PARALLEL_MATRIX.md`](../PARALLEL_MATRIX.md)):
run on **sandbox**, conda env `chemsmart`, `CHEMSMART_CONFIG_DIR` set.

DIAS (`N1`) still needs the existing IRC/opt log under
`$PROJECT/chemsmart/test_batch/` (linked by kind `dias`). All other inputs are
under `minimum/inputs/`.

## Quick start

```bash
ssh sandbox
source ~/miniconda3/etc/profile.d/conda.sh && conda activate chemsmart
export CHEMSMART_CONFIG_DIR=~/.chemsmart
cd /project/xlzhang/jingyi/chemsmart/test_batch/parallel_matrix/minimum
chmod +x run_minimum.sh

# Script-generation checks only (default)
DRY_RUN=1 ./run_minimum.sh

# Subset
CASES=B1,B6 DRY_RUN=1 ./run_minimum.sh
JOBS=gaussian/sp,orca/sp DRY_RUN=1 ./run_minimum.sh

# Live wall-time gate: five_mols gaussian SP serial vs parallel
DRY_RUN=0 WAIT=1 CASES=E3 ./run_minimum.sh
```

Default `DRY_RUN=1` passes `--test` (no `sbatch`). Live compute is intended for
the **E3 timing** gate only.

## Case policies (from the map)

| Case | Flags | Typical expect |
|------|--------|----------------|
| **B1** | `--no-run-in-parallel` | `#SBATCH --array=1-N%1` |
| **B2** | `--no-run-in-parallel -M 5` | still `%1` (`-M` ignored when serial) |
| **B6** | `--run-in-parallel -M 5` | `1-5%5` |
| **O2** | `--run-in-parallel -M 2` on 2 mols | array + distinct `--label` in `TASK_CLI` |
| **N1** | `--no-run-in-parallel` | non-array nestable parent |
| **N2** | `--run-in-parallel -M 2` | array + `--child-index` in shared runscript |
| **E3** | singles (non-array) + optional live SP timing | non-array submit; timing needs `DRY_RUN=0 WAIT=1` |

Canonical owners: `gaussian/sp`→B1, `gaussian/opt`→B2, `orca/sp`→B6,
`gaussian/td`→O2, `gaussian/dias`→N1, `gaussian/qrc`→N2; remaining job types
follow the same case policies (see the TSV).

## E3 timing

When `CASES` includes `E3` and `DRY_RUN=0 WAIT=1`, the runner times:

- B1-shaped: `gaussian … five_mols.xyz -i 1,2,3,4,5 sp` with `--no-run-in-parallel`
- B6-shaped: same SP with `--run-in-parallel -M 5`

This uses HF/`test` project settings on tiny H₂ structures — not full pKa.

## Results

- Workdirs: `parallel_matrix/runs/min_<case>_<program>_<job>_<timestamp>/`
- Summary: `parallel_matrix/logs/minimum_summary_<timestamp>.tsv`

## Environment

| Variable | Default | Meaning |
|----------|---------|---------|
| `DRY_RUN` | `1` | `1` → `--test` |
| `WAIT` | `0` | `1` → wait for E3 timing jobs |
| `CASES` | `B1,B2,B6,O2,N1,N2,E3` | case filter |
| `JOBS` | (all) | e.g. `gaussian/sp,orca/opt` |
| `MINIMUM_INPUTS` | `./inputs` | fixture directory |
| `GAUSSIAN_PROJECT_MIN` | `test` (`GAUSSIAN_PROJECT_OPT`) | project for most gaussian rows |
| `ORCA_PROJECT` | `test` | project for ORCA rows |
| `SERVER` | `cu_batch` | `-s` server yaml |
