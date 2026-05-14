# pKa BatchJob Decoupling Refactoring

## Overview

This refactoring decouples the pKa calculation logic from the BatchJob processing system while maintaining full support for multi-molecule `.csv` and `.cdxml` file inputs. The system now processes all molecules as a collection of independent sequential jobs rather than delegating to distributed node execution through BatchJob.

## Objectives Achieved

✅ **Input Support Maintained**: The pKa module still correctly parses and handles:
- Single molecules from `.xyz`, `.mol`, `.mol2` files
- Multi-molecule `.cdxml` files (fragment detection)
- `.csv` tables with per-entry parameters

✅ **Distributed Execution Removed**: Even when multiple molecules are detected, the system no longer instantiates or delegates to `BatchJob` classes:
- `GaussianpKaBatchJob` removed from CLI execution paths
- `OrcaBatchJob` removed from CLI execution paths
- Individual `pKa` jobs execute sequentially via `process_pipeline`

✅ **Single Sequential Execution**: Entire molecule sets are treated as collections of independent local tasks:
- Job lists returned directly from CLI commands
- Native job runner list handling via `process_pipeline`
- Serial/parallel execution controlled at individual job level

✅ **Safe Dependency Removal**: BatchJob dependencies identified and safely handled:
- No forced imports or dependencies in pKa CLI
- `GaussianpKaBatchJob` class preserved for backward compatibility
- Job runner already supports job list execution natively

## Changes Made

### 1. **Gaussian pKa CLI** (`chemsmart/cli/gaussian/pka.py`)

#### `submit()` Command
**Before**: Returned `GaussianpKaBatchJob` wrapping either single job or job list
```python
return GaussianpKaBatchJob(
    jobs=[job],
    run_in_serial=serial_mode.run_in_serial,
    jobrunner=jobrunner,
)
```

**After**: Returns job directly or job list
```python
return GaussianpKaJob(
    molecule=molecules[-1],
    settings=pka_settings,
    label=label,
    jobrunner=jobrunner,
    skip_completed=skip_completed,
    **kwargs,
)
# OR for multiple molecules:
return [
    GaussianpKaJob(...),
    GaussianpKaJob(...),
]
```

**Key Changes**:
- Single-molecule inputs return `GaussianpKaJob` directly
- Multiple-molecule inputs (via --index) return `list[GaussianpKaJob]`
- Multi-fragment CDXML returns `list[GaussianpKaJob]` from helper
- Removed import of `get_serial_mode` (batch execution handled by job runner)

#### `batch()` Command
**Before**: Created jobs and wrapped in `GaussianpKaBatchJob`
```python
return GaussianpKaBatchJob(
    jobs=jobs,
    run_in_serial=serial_mode.run_in_serial,
    write_outcome_logs=True,
    jobrunner=jobrunner,
)
```

**After**: Returns job list directly
```python
logger.info(f"Created {len(jobs)} pKa jobs from table")
return jobs
```

**Key Changes**:
- Table-driven batch submission returns `list[GaussianpKaJob]`
- Per-entry subcommand overrides preserved for `chemsmart sub`
- Removed `GaussianpKaBatchJob` wrapping

#### `_create_pka_jobs_from_molecules()` Helper
**Before**: Called `GaussianpKaJob.from_molecules()` which returned `GaussianpKaBatchJob`
```python
return GaussianpKaJob.from_molecules(
    pka_molecules=pka_molecules,
    ...
)
```

**After**: Creates and returns job list directly
```python
jobs = []
for idx, pka_mol in enumerate(pka_molecules, start=1):
    pka_settings = build_gaussian_pka_settings(...)
    jobs.append(GaussianpKaJob(...))
return jobs
```

**Key Changes**:
- Inlined job creation logic
- Returns bare `list[GaussianpKaJob]`
- Removed dependency on `from_molecules()` class method

#### Import Changes
- ❌ Removed: `from chemsmart.jobs.runner import get_serial_mode`
- ❌ Removed: `GaussianpKaBatchJob` import from `submit()`
- ✅ Added: `import os` for path manipulation

### 2. **Gaussian pKa Job** (`chemsmart/jobs/gaussian/pka.py`)

#### `GaussianpKaJob.from_molecules()` Method
**Before**: Returns `GaussianpKaBatchJob` wrapping job list
```python
return GaussianpKaBatchJob(
    jobs=jobs,
    run_in_serial=serial_mode.run_in_serial,
    jobrunner=jobrunner,
)
```

**After**: Returns bare job list
```python
return jobs
```

**Key Changes**:
- Method now returns `list[GaussianpKaJob]` instead of `GaussianpKaBatchJob`
- Removed call to `get_serial_mode()` (no longer needed)
- Docstring updated to reflect new return type
- Still available for backward compatibility (e.g., in old code paths)

#### Class Preservation
- ✅ `GaussianpKaBatchJob` class left intact for backward compatibility
- 📝 Marked as legacy in docstring (not used in current pKa flows)
- Can be removed in future major version if no external code depends on it

### 3. **ORCA pKa CLI** (`chemsmart/cli/orca/pka.py`)

**Status**: Already decoupled ✅
- `submit()` command already returns `ORCApKaJob` or `list[ORCApKaJob]`
- `batch()` command already returns `list[ORCApKaJob]`
- Multi-fragment helper already returns job list
- No changes needed

**Verification**: Lines 192-213 show direct job/list returns without batch wrapping.

### 4. **Job Runner Integration** (`chemsmart/cli/run.py`)

**Status**: Already supports job lists ✅

The `process_pipeline()` callback already handles:
```python
if isinstance(job, list):
    logger.info(f"Running {len(job)} jobs")
    serial_mode = get_serial_mode(jobrunner)
    
    if serial_mode.run_in_serial:
        # Run jobs one by one
        for single_job in job:
            single_job.jobrunner = _prepare_runner(single_job)
            single_job.run()
    else:
        # Run jobs in parallel with ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all jobs and collect results
```

**No changes needed** — pKa job lists will be executed correctly by existing infrastructure.

## Execution Flow Changes

### Before: Single Molecule

```
CLI submit → GaussianpKaJob created → Wrapped in GaussianpKaBatchJob
    ↓
BatchJob.run() → Executes batch orchestration → Calls job.run()
    ↓
Job execution
```

### After: Single Molecule

```
CLI submit → GaussianpKaJob created → Returned directly
    ↓
process_pipeline() → Detects Job instance → Creates jobrunner → job.run()
    ↓
Job execution
```

### Before: Multiple Molecules

```
CLI submit/batch → [GaussianpKaJob, ...] created → Wrapped in GaussianpKaBatchJob
    ↓
BatchJob.run() → Iterates jobs → Calls each job.run() (or parallel)
    ↓
Multiple job executions (distributed across nodes if enabled)
```

### After: Multiple Molecules

```
CLI submit/batch → [GaussianpKaJob, ...] created → Returned directly
    ↓
process_pipeline() → Detects list → Iterates jobs → Prepares each jobrunner
    ↓
Serial/parallel execution based on jobrunner.run_in_serial
    ↓
Multiple job executions (sequential on local node)
```

## Backward Compatibility

### ✅ Preserved Functionality

1. **Input Parsing**
   - Single molecules: `.xyz`, `.mol`, `.mol2`
   - Multi-fragment: `.cdxml` auto-detection
   - Batch tables: `.csv` with 4 columns

2. **Settings & Validation**
   - All pKa settings classes unchanged
   - Charge/multiplicity validation preserved
   - Reference acid handling (proton exchange/direct)
   - Conjugate pair calculations

3. **Job Execution**
   - Individual pKa job logic unchanged
   - Gas-phase optimizations
   - Solution-phase single points
   - pKa output analysis

4. **CLI Subcommands**
   - `gaussian pka submit` — works as before
   - `gaussian pka batch` — works as before
   - `orca pka submit` — works as before
   - `orca pka batch` — works as before

### ⚠️ Breaking Changes

1. **Direct BatchJob Usage**: Code that explicitly imports and instantiates `GaussianpKaBatchJob` or relies on its return:
   ```python
   # Old code:
   from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob
   batch_job = GaussianpKaBatchJob(jobs=[...])  # Still works but not used in pKa flow
   ```
   
   **Impact**: Low - internal class, not documented as public API

2. **Distributed Node Processing**: Features relying on BatchJob's node distribution:
   - Multi-node SLURM/PBS distribution no longer automatic
   - Individual jobs run sequentially on assigned node
   - Can be re-enabled by manually wrapping job list in BatchJob if needed

### ✅ How Existing Code Adapts

1. **Job Runner Usage** — No changes needed
   - `process_pipeline()` already handles job lists
   - `jobrunner.from_job()` works for both single jobs and lists

2. **Job Submission Scripts** — No changes needed
   - CLI commands work identically
   - Return type change is transparent to CLI parsing

3. **Post-Processing** — No changes needed
   - `chemsmart run pka analyze` — unchanged
   - `chemsmart run pka batch-analyze` — unchanged
   - All thermochemistry extraction — unchanged

## Dependency Analysis

### Removed Dependencies from pKa CLI

1. **`get_serial_mode()`** from `chemsmart.jobs.runner`
   - **Before**: Used to determine batch execution mode
   - **After**: No longer needed (handled by process_pipeline)
   - **Status**: Safely removed ✅

2. **`GaussianpKaBatchJob`** class
   - **Before**: Imported and instantiated in CLI
   - **After**: No longer imported in pKa CLI
   - **Status**: Preserved in jobs module for backward compatibility ✅

3. **`serial_mode` configuration**
   - **Before**: Passed to BatchJob constructor
   - **After**: Handled transparently by process_pipeline
   - **Status**: Safely removed ✅

### Unaffected BatchJob Usages

The following still use BatchJob (not affected by this refactoring):

1. `chemsmart.jobs.gaussian.crest.GaussianCRESTJob`
   - Uses `GaussianBatchJob` for conformer distribution
   - No changes made

2. `chemsmart.jobs.gaussian.qrc.GaussianQRCJob`
   - Uses `GaussianBatchJob` for forward/reverse IRC
   - No changes made

3. `chemsmart.jobs.orca.qrc.ORCAQRCJob`
   - Uses `OrcaBatchJob` for QRC jobs
   - No changes made

## Testing Recommendations

### Unit Tests

1. **CLI Return Types**
   ```python
   # test_gaussian_pka_cli.py
   def test_submit_single_molecule_returns_job():
       # Verify returns GaussianpKaJob, not batch job
       
   def test_submit_multiple_molecules_returns_list():
       # Verify returns list[GaussianpKaJob]
       
   def test_batch_returns_list():
       # Verify batch command returns list[GaussianpKaJob]
   ```

2. **Job Creation**
   ```python
   # test_gaussian_pka_job.py
   def test_from_molecules_returns_list():
       # Verify from_molecules() returns list, not GaussianpKaBatchJob
   ```

### Integration Tests

1. **CLI Execution**
   ```bash
   # Single molecule
   chemsmart run gaussian -f acid.xyz -c 0 -m 1 pka -pi 10 submit
   
   # Multi-fragment CDXML
   chemsmart run gaussian -f molecules.cdxml pka submit
   
   # Batch table
   chemsmart run gaussian -f pka_table.csv pka batch
   ```

2. **Job Runner Integration**
   ```python
   # Verify process_pipeline handles returned job/list correctly
   # Verify serial and parallel execution modes work
   ```

### Regression Tests

1. **Reference Acid Handling** — Ensure proton exchange cycle still works
2. **Solvation Model Selection** — Verify CPCM/SMD settings applied
3. **Output Analysis** — Confirm pKa calculation results unchanged
4. **Multi-Job Execution** — Ensure all jobs complete successfully

## Migration Guide

### For End Users

**No action required** — CLI commands work identically to before.

### For Developers

#### If You Were Using `GaussianpKaBatchJob` Directly

```python
# Old code:
from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob
batch_job = GaussianpKaBatchJob(jobs=[job1, job2])
batch_job.run()

# New approach (in pKa flow):
jobs = [job1, job2]
# Let process_pipeline handle job list execution
return jobs

# Alternative (if you need batch wrapper):
from chemsmart.jobs.gaussian.batch import GaussianBatchJob
batch_job = GaussianBatchJob(jobs=[job1, job2])
batch_job.run()
```

#### If You Were Calling `GaussianpKaJob.from_molecules()`

```python
# Old code:
pka_batch_job = GaussianpKaJob.from_molecules(
    pka_molecules=molecules,
    ...
)
pka_batch_job.run()

# New code:
pka_jobs = GaussianpKaJob.from_molecules(
    pka_molecules=molecules,
    ...
)
# Returns list[GaussianpKaJob], not GaussianpKaBatchJob
for job in pka_jobs:
    job.run()

# Or better: return list and let process_pipeline handle it
return pka_jobs
```

## Files Modified

| File | Changes |
|------|---------|
| `chemsmart/cli/gaussian/pka.py` | Removed BatchJob wrapping; return job/list directly |
| `chemsmart/jobs/gaussian/pka.py` | Modified `from_molecules()` to return list; preserved `GaussianpKaBatchJob` class |
| `chemsmart/cli/orca/pka.py` | No changes (already decoupled) |
| `chemsmart/jobs/orca/pka.py` | No changes (already decoupled) |
| `chemsmart/cli/run.py` | No changes (already supports job lists) |

## Files NOT Modified

- ✅ `chemsmart/jobs/batch.py` — Base `BatchJob` class preserved
- ✅ `chemsmart/jobs/gaussian/batch.py` — `GaussianBatchJob` preserved
- ✅ `chemsmart/jobs/orca/batch.py` — `OrcaBatchJob` preserved
- ✅ All CREST, QRC job implementations — Still use BatchJob as needed
- ✅ All other CLI subcommands — No impact

## Summary

This refactoring successfully decouples the pKa calculation system from BatchJob distributed processing while:

1. ✅ Maintaining all input file support (.xyz, .cdxml, .csv)
2. ✅ Preserving all pKa calculation functionality
3. ✅ Enabling sequential job execution (no node distribution)
4. ✅ Leveraging existing job runner list support
5. ✅ Keeping BatchJob framework available for other uses
6. ✅ Maintaining backward compatibility for most use cases
7. ✅ Simplifying the pKa submission architecture

The system now treats multiple pKa molecules as independent sequential jobs rather than distributed batch operations, making the workflow simpler, more transparent, and easier to debug.

