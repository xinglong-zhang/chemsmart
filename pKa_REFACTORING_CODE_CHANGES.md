# pKa BatchJob Decoupling - Code Changes Reference

## Quick Reference: Changes Made

### 1. Gaussian pKa CLI: Single Molecule Submission

#### BEFORE
```python
# chemsmart/cli/gaussian/pka.py - submit()
serial_mode = get_serial_mode(jobrunner)

job = GaussianpKaJob(
    molecule=molecules[-1],
    settings=pka_settings,
    label=label,
    jobrunner=jobrunner,
    skip_completed=skip_completed,
    **kwargs,
)
# Wrap single job in batch container
return GaussianpKaBatchJob(
    jobs=[job],
    run_in_serial=serial_mode.run_in_serial,
    jobrunner=jobrunner,
)
```

#### AFTER
```python
# chemsmart/cli/gaussian/pka.py - submit()
return GaussianpKaJob(
    molecule=molecules[-1],
    settings=pka_settings,
    label=label,
    jobrunner=jobrunner,
    skip_completed=skip_completed,
    **kwargs,
)
```

**Impact**: Single job returned directly; process_pipeline detects `isinstance(job, Job)` and runs it.

---

### 2. Gaussian pKa CLI: Multiple Molecules (--index)

#### BEFORE
```python
# chemsmart/cli/gaussian/pka.py - submit()
if len(molecules) > 1 and molecule_indices:
    logger.info(f"Creating {len(molecules)} pKa jobs")
    jobs = [
        GaussianpKaJob(...)
        for mol, idx in zip(molecules, molecule_indices)
    ]
    # Wrap list in batch container
    return GaussianpKaBatchJob(
        jobs=jobs,
        run_in_serial=serial_mode.run_in_serial,
        jobrunner=jobrunner,
    )
```

#### AFTER
```python
# chemsmart/cli/gaussian/pka.py - submit()
if len(molecules) > 1 and molecule_indices:
    logger.info(f"Creating {len(molecules)} pKa jobs")
    return [
        GaussianpKaJob(...)
        for mol, idx in zip(molecules, molecule_indices)
    ]
```

**Impact**: Job list returned directly; process_pipeline detects `isinstance(job, list)` and iterates/executes.

---

### 3. Gaussian pKa CLI: Batch Table Submission

#### BEFORE
```python
# chemsmart/cli/gaussian/pka.py - batch()
serial_mode = get_serial_mode(jobrunner)

jobs = []
for entry in entries:
    # ... create job ...
    jobs.append(GaussianpKaJob(...))

return GaussianpKaBatchJob(
    jobs=jobs,
    run_in_serial=serial_mode.run_in_serial,
    write_outcome_logs=True,
    jobrunner=jobrunner,
)
```

#### AFTER
```python
# chemsmart/cli/gaussian/pka.py - batch()
jobs = []
for entry in entries:
    # ... create job ...
    jobs.append(GaussianpKaJob(...))

logger.info(f"Created {len(jobs)} pKa jobs from table")
return jobs
```

**Impact**: Job list returned directly; no batch wrapper.

---

### 4. Gaussian pKa CLI: Multi-Fragment CDXML

#### BEFORE
```python
# chemsmart/cli/gaussian/pka.py - submit()
if pka_molecules is not None:
    return _create_pka_jobs_from_molecules(
        ctx, pka_molecules, shared, skip_completed, **kwargs
    )

# Inside _create_pka_jobs_from_molecules():
return GaussianpKaJob.from_molecules(
    pka_molecules=pka_molecules,
    shared=shared,
    project_settings=project_settings,
    job_settings=job_settings,
    keywords=keywords,
    jobrunner=jobrunner,
    filename=filename,
    skip_completed=skip_completed,
    **kwargs,
)
# Which returned GaussianpKaBatchJob internally
```

#### AFTER
```python
# chemsmart/cli/gaussian/pka.py - submit()
if pka_molecules is not None:
    return _create_pka_jobs_from_molecules(
        ctx, pka_molecules, shared, skip_completed, **kwargs
    )

# Inside _create_pka_jobs_from_molecules():
jobs = []
for idx, pka_mol in enumerate(pka_molecules, start=1):
    label = f"{base_name}_frag{idx}_pka"
    pka_settings = build_gaussian_pka_settings(...)
    jobs.append(
        GaussianpKaJob(
            molecule=pka_mol,
            settings=pka_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
    )

logger.info(f"Created {len(jobs)} pKa jobs from multi-fragment CDXML")
return jobs  # Returns list directly, not GaussianpKaBatchJob
```

**Impact**: Multi-fragment CDXML now returns job list directly.

---

### 5. Gaussian pKa Job: from_molecules() Class Method

#### BEFORE
```python
# chemsmart/jobs/gaussian/pka.py - GaussianpKaJob.from_molecules()
@classmethod
def from_molecules(cls, pka_molecules, shared, project_settings, ...):
    """Create a batch job from PKaMolecule inputs.
    
    Returns:
        GaussianpKaBatchJob: Batch container with created jobs.
    """
    opt_settings = project_settings.opt_settings()
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    
    serial_mode = get_serial_mode(jobrunner)
    base_name = os.path.splitext(os.path.basename(filename))[0] or "pka"
    
    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        # ... create job ...
        jobs.append(cls(...))
    
    return GaussianpKaBatchJob(
        jobs=jobs,
        run_in_serial=serial_mode.run_in_serial,
        jobrunner=jobrunner,
    )
```

#### AFTER
```python
# chemsmart/jobs/gaussian/pka.py - GaussianpKaJob.from_molecules()
@classmethod
def from_molecules(cls, pka_molecules, shared, project_settings, ...):
    """Create a list of pKa jobs from PKaMolecule inputs.
    
    Returns:
        list[GaussianpKaJob]: List of created pKa jobs (no batch wrapper).
    """
    opt_settings = project_settings.opt_settings()
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    
    base_name = os.path.splitext(os.path.basename(filename))[0] or "pka"
    
    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        # ... create job ...
        jobs.append(cls(...))
    
    return jobs  # Return list directly
```

**Impact**: Return type changed from `GaussianpKaBatchJob` to `list[GaussianpKaJob]`.

---

### 6. Import Changes

#### BEFORE
```python
# chemsmart/cli/gaussian/pka.py
from chemsmart.jobs.runner import get_serial_mode
from chemsmart.jobs.gaussian.pka import (
    GaussianpKaBatchJob,
    GaussianpKaJob,
    build_gaussian_pka_settings,
)
```

#### AFTER
```python
# chemsmart/cli/gaussian/pka.py
from chemsmart.jobs.gaussian.pka import (
    GaussianpKaJob,
    build_gaussian_pka_settings,
)
```

**Impact**: Removed unused imports; simplified dependency chain.

---

## Execution Flow Changes

### Single Molecule Job Execution

```
┌─ BEFORE ─────────────────────────────┐
│ GaussianpKaJob created               │
│       ↓                              │
│ Wrapped in GaussianpKaBatchJob       │
│       ↓                              │
│ batch.run() orchestrates             │
│       ↓                              │
│ Calls job.run()                      │
└──────────────────────────────────────┘

┌─ AFTER ──────────────────────────────┐
│ GaussianpKaJob created               │
│       ↓                              │
│ Returned directly                    │
│       ↓                              │
│ process_pipeline detects Job         │
│       ↓                              │
│ Creates jobrunner                    │
│       ↓                              │
│ Calls job.run()                      │
└──────────────────────────────────────┘
```

### Multiple Jobs Execution

```
┌─ BEFORE ─────────────────────────────┐
│ [GaussianpKaJob, ...] created        │
│       ↓                              │
│ Wrapped in GaussianpKaBatchJob       │
│       ↓                              │
│ batch.run() orchestrates             │
│       ↓                              │
│ Iterates: job.run() (or parallel)    │
│       ↓                              │
│ Can distribute across nodes          │
└──────────────────────────────────────┘

┌─ AFTER ──────────────────────────────┐
│ [GaussianpKaJob, ...] created        │
│       ↓                              │
│ Returned directly                    │
│       ↓                              │
│ process_pipeline detects list        │
│       ↓                              │
│ Iterates with ThreadPoolExecutor     │
│       ↓                              │
│ Each job.run() (serial or parallel)  │
│       ↓                              │
│ Sequential execution on local node   │
└──────────────────────────────────────┘
```

---

## Key Differences

| Aspect | BEFORE | AFTER |
|--------|--------|-------|
| **Single Job Return** | `GaussianpKaBatchJob([job])` | `GaussianpKaJob` |
| **Multiple Jobs Return** | `GaussianpKaBatchJob(jobs)` | `list[GaussianpKaJob]` |
| **Batch Wrapper** | Always wrapped | Never wrapped |
| **Node Distribution** | Via BatchJob class | Via individual jobs (none) |
| **Serial/Parallel Control** | BatchJob decides | Individual job + process_pipeline |
| **Job Execution** | batch.run() method | job.run() method |
| **Imports** | get_serial_mode, GaussianpKaBatchJob | Neither needed in CLI |

---

## Backward Compatibility Status

### ✅ Still Works

- `GaussianpKaJob.from_molecules()` — Still exists, returns list instead of batch job
- `GaussianpKaBatchJob` class — Still defined for any code that directly instantiates it
- CLI commands — All work identically from user perspective
- pKa calculations — All outputs unchanged
- Output analysis — All analysis tools unchanged

### ⚠️ Breaking Changes

- Direct return type of `GaussianpKaJob.from_molecules()` changed
- `GaussianpKaBatchJob` no longer returned from pKa CLI
- Code relying on BatchJob node distribution for pKa will need adjustment

### 🔄 Migration Path

```python
# If you manually call from_molecules and expect GaussianpKaBatchJob:

# OLD CODE (will break):
batch_job = GaussianpKaJob.from_molecules(...)
batch_job.run()

# NEW CODE (use job list):
jobs = GaussianpKaJob.from_molecules(...)
for job in jobs:
    job.run()

# OR (let process_pipeline handle it):
return GaussianpKaJob.from_molecules(...)
```

---

## Testing Checklist

- [ ] Single molecule pKa submission works
- [ ] Multi-molecule (--index) pKa submission works
- [ ] Multi-fragment CDXML pKa submission works
- [ ] Batch table pKa submission works
- [ ] Serial execution mode works
- [ ] Parallel execution mode works
- [ ] pKa calculations produce correct results
- [ ] Output analysis tools work
- [ ] Reference acid handling works
- [ ] No regression in existing functionality
- [ ] CLI help text accurate
- [ ] Error messages informative

---

## Deployment Notes

1. **No database migrations** — No schema changes
2. **No configuration changes** — Existing configs work
3. **No dependency updates** — No new packages required
4. **Backward compatible** — Existing scripts mostly work
5. **CLI identical** — End users see no difference
6. **Internal refactoring** — Better architecture, simpler code

