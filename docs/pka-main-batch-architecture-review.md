# Architectural Integrity Review: pKa Refactor

## Scope and review lens

This review evaluates the current refactor state for:

- JobRunner lifecycle and copy safety
- shared phase-runner behavior across pKa phases
- consistency of `get_serial_mode(...)` as serial/parallel source of truth
- batch scalability and submission behavior under large job arrays
- containment of pKa-specific `fail_fast` semantics

Code reviewed includes:

- `chemsmart/cli/run.py`
- `chemsmart/cli/gaussian/pka.py`
- `chemsmart/cli/orca/pka.py`
- `chemsmart/jobs/runner.py`
- `chemsmart/jobs/job.py`
- `chemsmart/jobs/batch.py`
- `chemsmart/jobs/gaussian/pka.py`
- `chemsmart/jobs/orca/pka.py`

---

## Findings (ordered by severity)

### 1) High — Batch default policy now biases to serial unless caller explicitly overrides(FIXED)

- **Evidence:** `chemsmart/jobs/batch.py:53-69` (`BatchJob.__init__`)
- **What changed logically:** `self.run_in_serial` is now `run_in_serial OR runner.run_in_serial` (via `get_serial_mode`).
- **Risk:** because constructor default is `run_in_serial=True`, any call site that does not pass `run_in_serial` explicitly remains serial even when runner policy is parallel. This can silently suppress intended parallel behavior in batch subclasses that rely on defaults.
- **Impact area:** scalability and policy consistency.
- **Recommendation:** make precedence explicit and non-ambiguous:
  - either change constructor default to `run_in_serial=None` and resolve from runner when `None`,
  - or keep default but enforce explicit passing at all `BatchJob` call sites and document that requirement.

### 2) High — Parallel pKa path does not gate SP phase on optimization failures (FIXED)

- **Evidence:**
  - `chemsmart/jobs/gaussian/pka.py:530-573` (`GaussianpKaJob._run_parallel`)
  - `chemsmart/jobs/orca/pka.py:608-651` (`ORCApKaJob._run_parallel`)
- **What happens:** after collecting opt phase failures, code still proceeds to create and run SP jobs.
- **Risk:** SP jobs may run on fallback/unoptimized geometries after failed opt jobs; this breaks phase safety expectation and can waste resources at scale.
- **Recommendation:** gate SP phase with `if all_failures: raise` (or phase-level policy flag) before SP dispatch.

### 3) Medium — Error handling in batch parallel execution hides failures from callers (FIXED)

- **Evidence:**
  - `_submit_job` catches and logs exceptions internally (`chemsmart/jobs/batch.py:142-154`)
  - `_run_jobs_in_parallel` also catches future exceptions, but futures usually do not fail because `_submit_job` swallows exceptions (`chemsmart/jobs/batch.py:120-135`)
- **Risk:** parent workflow cannot reliably detect batch failure without polling output/state; this can produce false-positive orchestration success.
- **Recommendation:** aggregate per-job failures and return/raise structured summary from `BatchJob.run` for upstream phase control.

### 4) Medium — `get_serial_mode` is mostly centralized, but not yet complete repo-wide

- **Evidence:** direct boolean branch remains outside reviewed pKa paths in `chemsmart/cli/sub.py:241` (`not jobrunner.run_in_serial`).
- **Risk:** policy drift if serial-mode semantics evolve (e.g., tri-state, environment overrides) and non-central callers remain.
- **Recommendation:** migrate remaining direct checks to `get_serial_mode` or formally scope helper to pKa + run pipeline only.

### 5) Medium — Potential oversubscription hotspots in large batches

- **Evidence:**
  - `chemsmart/cli/run.py:138` uses `ThreadPoolExecutor(max_workers=len(job))`
  - `chemsmart/jobs/batch.py:120` uses `ThreadPoolExecutor()` default worker count
- **Risk:** large batch lists can create too many concurrent submitters while each job itself may be multi-core; this can bottleneck scheduler APIs and local head nodes.
- **Recommendation:** cap worker count using a policy-aware bound (`min(len(jobs), configured_max_submitters)`) and optionally tie to runner/server capabilities.

### 6) Low — Phase-runner abstraction is centralized but parallel branch still bypasses it

- **Evidence:**
  - shared abstraction present in `run_phase_jobs` -> `Job._execute_phase_jobs` (`chemsmart/jobs/runner.py:45-66`, `chemsmart/jobs/job.py:127-167`)
  - both engine `_run_parallel` implementations use custom loops (`chemsmart/jobs/gaussian/pka.py:458-591`, `chemsmart/jobs/orca/pka.py:556-670`)
- **Risk:** behavior differences can reappear between serial and parallel paths (phase gating, reporting, stop semantics).
- **Recommendation:** keep custom resource splitting, but move phase transition decisions (opt -> sp gating) through a shared phase contract.

---

## JobRunner lifecycle trace (integrity check)

1. **Instantiation (CLI root):**
   - `chemsmart/cli/run.py:65-74` creates base `JobRunner` with CLI flags.
2. **Type specialization:**
   - `process_pipeline` creates engine-specific runner via `from_job(...)` for each job (`chemsmart/cli/run.py:116-127`, `160-170`).
3. **Batch propagation:**
   - `BatchJob._build_jobrunner` uses `Job._propagate_runner(...)` copy semantics (`chemsmart/jobs/batch.py:155-173`).
4. **Phase propagation:**
   - `Job._execute_phase_jobs` copies runner per child (`chemsmart/jobs/job.py:152-156`).
5. **Parallel worker propagation:**
   - pKa parallel workers also call `_propagate_runner` per job (`chemsmart/jobs/gaussian/pka.py:399-425`, `chemsmart/jobs/orca/pka.py:497-523`).

**Conclusion:** runner copy discipline is generally sound; no major shared-reference mutation risk found in the main pKa execution flow.

---

## `get_serial_mode` consistency check

### In-scope paths (good)

- CLI list execution: `chemsmart/cli/run.py:114-154`
- Gaussian pKa CLI: `chemsmart/cli/gaussian/pka.py` (submit/batch/analyze/batch-analyze/thermo)
- ORCA pKa CLI: `chemsmart/cli/orca/pka.py:151-153`, `241-243`
- pKa job classes: `chemsmart/jobs/gaussian/pka.py`, `chemsmart/jobs/orca/pka.py`
- BatchJob constructor integration: `chemsmart/jobs/batch.py:65-69`

### Drift observed

- Non-pKa path still performs local boolean math in `chemsmart/cli/sub.py:241`.

---

## Fail-fast encapsulation check

- `get_serial_mode(...)` contains only serial-mode state (`chemsmart/jobs/runner.py:31-42`) and does **not** include fail-fast.
- pKa-specific stop behavior remains phase-level (`stop_on_incomplete=True` in pKa phase calls, serial-gated in `Job._execute_phase_jobs`).

**Conclusion:** fail-fast remains correctly encapsulated and has not leaked into shared serial-mode utility.

---

## Methods currently deviating from the intended pattern

1. `chemsmart/jobs/gaussian/pka.py:458-591` — `_run_parallel` (custom phase transition path bypasses shared phase contract)
2. `chemsmart/jobs/orca/pka.py:556-670` — `_run_parallel` (same deviation)
3. `chemsmart/jobs/batch.py:142-154` — `_submit_job` (swallows exceptions; weak upstream failure propagation)
4. `chemsmart/cli/sub.py:241` — local serial-mode boolean logic outside `get_serial_mode`

---

## Recommended next refactor slice

1. Normalize `BatchJob` serial precedence (`None` sentinel or explicit-callsite contract).
2. Add opt->SP gate in both pKa `_run_parallel` methods.
3. Introduce batch failure aggregation (return/raise summary from batch run).
4. Replace remaining direct `run_in_serial` branching with `get_serial_mode` where policy uniformity is required.

---

## Overall structural assessment

The refactor has strong architectural direction: runner propagation is centralized and copy-safe, phase orchestration for serial path is shared, and serial-mode policy is mostly unified behind `get_serial_mode`. Remaining risks are concentrated in parallel phase gating and batch failure surfacing, not in core object lifecycle correctness.
