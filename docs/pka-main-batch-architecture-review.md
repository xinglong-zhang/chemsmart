# Architectural Integrity Review: pKa Refactor (Current State)

## Scope and review lens

This review evaluates the current refactor state for:

- JobRunner lifecycle and copy safety
- shared phase-runner behavior across pKa phases
- consistency of `get_serial_mode(...)` as serial/parallel source of truth
- batch scalability and submission behavior under large job arrays
- containment of pKa-specific fail/stop semantics

Code reviewed includes:

- `chemsmart/cli/run.py`
- `chemsmart/cli/sub.py`
- `chemsmart/cli/gaussian/pka.py`
- `chemsmart/cli/orca/pka.py`
- `chemsmart/jobs/runner.py`
- `chemsmart/jobs/job.py`
- `chemsmart/jobs/batch.py`
- `chemsmart/jobs/gaussian/pka.py`
- `chemsmart/jobs/orca/pka.py`
- `tests/test_runner.py`

---

## Status snapshot

- Batch serial precedence fix: **FIXED**
- Parallel Opt -> SP gating fix: **FIXED**
- Batch failure surfacing to callers: **FIXED**
- `get_serial_mode` centralization in reviewed paths: **FIXED**
- Submitter oversubscription capping: **FIXED**
- Shared phase transition contract for parallel/serial: **PARTIALLY FIXED** (transition contract shared; worker-dispatch loops still engine-local)

---

## Findings (ordered by severity)

### 1) Medium - Multi-node coordinator exceptions are now converted into explicit batch failures
(FIXED)
- **Evidence:** `BatchJob._run_multi_node` now appends a synthetic failed outcome when a node future raises (`chemsmart/jobs/batch.py`).
- **What happens now:** node-level coordination errors become structured failure entries and are included in aggregate batch outcome handling.
- **Impact:** upstream callers reliably see node-level failures through the same failure channel as job-level failures.
- **Follow-up:** keep regression coverage for node-future exception paths to prevent drift.

### 2) Low - Parallel worker execution remains engine-specific even though phase transition policy is centralized

- **Evidence:** both `GaussianpKaJob._run_parallel` and `ORCApKaJob._run_parallel` keep custom `ThreadPoolExecutor` loops and per-engine worker functions (`chemsmart/jobs/gaussian/pka.py`, `chemsmart/jobs/orca/pka.py`).
- **What improved:** phase transition/gating is now shared through `decide_phase_transition(...)` in both serial and parallel paths.
- **Residual risk:** behavior drift can reappear in future if one engine changes worker collection/reporting semantics and the other does not.
- **Recommendation:** keep resource-splitting engine-local, but consider extracting a shared parallel phase collector helper (future map -> successes/failures -> transition decision).

---

## Previously flagged issues now resolved

### 1) Batch default serial bias due to ambiguous constructor precedence - FIXED

- `BatchJob.__init__` now uses `run_in_serial: Optional[bool] = None`.
- Resolution is explicit:
  - `None` -> inherit from `get_serial_mode(jobrunner).run_in_serial`
  - explicit `True`/`False` -> call-site override wins
- Covered by tests in `tests/test_runner.py` (`test_batch_serial_mode_when_unset`, `test_batch_serial_mode_keeps_explicit_false`, `test_batch_serial_mode_keeps_explicit_true`).

### 2) Parallel pKa Opt failure not gating SP dispatch - FIXED

- Both engines now gate via `decide_phase_transition(phase_name="Optimization", failures=phase_failures)` before SP job creation/dispatch.
- If optimization failures exist, `_run_parallel` raises and does not proceed to SP.
- Covered by tests:
  - `test_gaussian_pka_parallel_stops_before_sp_on_opt_failure`
  - `test_orca_pka_parallel_stops_before_sp_on_opt_failure`

### 3) Batch failures hidden from callers - FIXED

- Batch execution now aggregates per-job outcomes and raises `BatchExecutionError` from `BatchJob._run` when any child reports failure.
- Optional outcome logs are controlled by `write_outcome_logs` (default `False`), avoiding global side effects.
- Covered by tests (`test_batch_writes_success_and_failed_logs`, `test_batch_run_raises_with_failed_job_summary`).

### 4) `get_serial_mode` centralization incomplete - FIXED in reviewed paths

- Reviewed run/sub/pKa paths now route serial-policy decisions through `get_serial_mode(...)`.
- `chemsmart/cli/sub.py` no longer uses direct `not jobrunner.run_in_serial` logic.

### 5) Oversubscription hotspots in list/batch dispatch - FIXED

- `chemsmart/cli/run.py` and `chemsmart/jobs/batch.py` now cap workers through `get_submitter_worker_count(...)`.
- Cap derives from policy in `get_configured_max_submitters(...)` (`CHEMSMART_MAX_SUBMITTERS`, runner/server caps, cores, CPU fallback).
- Covered by tests in `tests/test_runner.py` for policy cap and env override.

---

## JobRunner lifecycle trace (integrity check)

1. **Instantiation (CLI root):**
   - `run`/`sub` create base `JobRunner` using CLI options (`chemsmart/cli/run.py`, `chemsmart/cli/sub.py`).
2. **Type specialization:**
   - pipeline resolves engine-specific runners via `from_job(...)` (`chemsmart/cli/run.py`).
3. **Batch propagation:**
   - `BatchJob._build_jobrunner` uses `Job._propagate_runner(...)` copy semantics (`chemsmart/jobs/batch.py`).
4. **Phase propagation:**
   - `run_phase_jobs(...)` delegates to `Job._execute_phase_jobs(...)`, which copies runner per child (`chemsmart/jobs/runner.py`, `chemsmart/jobs/job.py`).
5. **Parallel worker propagation:**
   - pKa parallel workers also propagate patched runner copies per worker job (`chemsmart/jobs/gaussian/pka.py`, `chemsmart/jobs/orca/pka.py`).

**Conclusion:** runner copy discipline remains sound in the reviewed pKa and batch execution paths.

---

## `get_serial_mode` consistency check

### In-scope paths (good)

- CLI list execution: `chemsmart/cli/run.py`
- Submission path: `chemsmart/cli/sub.py`
- Gaussian pKa CLI: `chemsmart/cli/gaussian/pka.py`
- ORCA pKa CLI: `chemsmart/cli/orca/pka.py`
- pKa job classes: `chemsmart/jobs/gaussian/pka.py`, `chemsmart/jobs/orca/pka.py`
- Batch constructor precedence: `chemsmart/jobs/batch.py`

### Drift observed

- No direct drift found in reviewed run/sub/pKa paths.

---

## pKa-specific stop semantics encapsulation

- `get_serial_mode(...)` remains intentionally limited to serial-mode state only (`chemsmart/jobs/runner.py`).
- Phase gating and failure decisions remain in pKa orchestration via `decide_phase_transition(...)` and phase-specific completion checks.

**Conclusion:** pKa-specific stop behavior is still encapsulated outside shared serial-mode utility.

---

## Methods with residual deviation from full unification

1. `chemsmart/jobs/gaussian/pka.py` - `_run_parallel` keeps engine-local parallel dispatch loop (intentional for resource splitting).
2. `chemsmart/jobs/orca/pka.py` - `_run_parallel` keeps engine-local parallel dispatch loop (same rationale).

---

## Recommended next refactor slice

1. Extract a shared parallel phase collector utility for both engines while retaining engine-specific resource partitioning.
2. Add a short operator-facing doc note for `CHEMSMART_MAX_SUBMITTERS` and submitter-cap resolution order.

---

## Overall structural assessment

The architecture is materially improved versus the earlier review baseline. Core lifecycle safety, serial-policy centralization, phase-transition gating, failure surfacing, and worker capping are now in place and test-backed. Remaining gaps are focused and low risk, primarily around optional further deduplication of parallel worker loops.
