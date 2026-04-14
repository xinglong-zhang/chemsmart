# Architectural Review: `pka` Branch vs `main`

## Scope

This review compares batch job processing architecture between `pka` and `main`, with focus on:

- `run_in_serial` behavior
- jobrunner propagation
- DRY refactor opportunities in job orchestration

Per request, Gaussian-vs-ORCA chemistry behavior differences are not counted as duplication unless they duplicate engine-agnostic orchestration logic.

---

## Executive Differences

### 1) `run_in_serial` support changed from implicit to explicit (but is unevenly applied)

In `main`, list execution was serial by default in `chemsmart/cli/run.py`, but there was no runner-level serial flag.

In `pka`, `run_in_serial` is added and plumbed through:

- `chemsmart/cli/run.py`
- `chemsmart/jobs/runner.py`
- `chemsmart/cli/jobrunner.py`

However, generic list execution in `chemsmart/cli/run.py` still executes lists serially in practice; the flag mostly changes logging there.

### 2) Strong new abstraction: shared batch orchestrator

`pka` introduces `chemsmart/jobs/batch.py`, which centralizes:

- serial/parallel child submission
- child-level fault tolerance
- SLURM/PBS node discovery
- multi-node chunking and dispatch
- runner copy/propagation via `Job._propagate_runner`

This is a major architectural improvement over `main`, where no shared batch layer exists.

### 3) Runner propagation became centralized

`pka` adds `Job._propagate_runner(...)` in `chemsmart/jobs/job.py`, creating one canonical copy-and-override mechanism.

This is good DRY progress, but some pKa flows still repeat execution-policy patterns around this helper.

### 4) pKa policy handling remains duplicated across layers

Execution policy (`serial`/`parallel`/`fail-fast`) is currently split between:

- CLI wrappers (`chemsmart/cli/gaussian/pka.py`, `chemsmart/cli/orca/pka.py`)
- job classes (`chemsmart/jobs/gaussian/pka.py`, `chemsmart/jobs/orca/pka.py`)
- batch infrastructure (`chemsmart/jobs/batch.py`)

This causes drift and repeated logic.

---

## DRY Refactor Targets (Specific Methods/Lines)

### A. Generic run pipeline

- `chemsmart/cli/run.py:111-136` (`process_pipeline` list branch)
  - Currently iterates list serially regardless of effective policy.
  - Should delegate list orchestration to a shared policy executor.

### B. Gaussian pKa CLI hard-coded batch serial mode

- `chemsmart/cli/gaussian/pka.py`
  - Serial hardcoding appears at: `231-235`, `246-248`, `391-393`, `445-447`, `456-458`, `558-560`, `569-571`, `726-728`.
  - Replace with centralized execution-policy resolver.

### C. ORCA pKa CLI policy derivation duplicated in two places

- `chemsmart/cli/orca/pka.py:151-153`, `242-244`
  - `parallel = not jobrunner.run_in_serial` repeated.
  - Move to shared helper used by both engines.

### D. ORCA pKa serial stop rules duplicated per phase

- `chemsmart/jobs/orca/pka.py`
  - `_run_opt_jobs` (`380-390`)
  - `_run_ref_opt_jobs` (`392-405`)
  - `_run_sp_jobs` (`407-418`)
  - `_run_ref_sp_jobs` (`420-434`)
  - plus stage gating in `_run` (`684-712`)

All repeat the same pattern: propagate runner, run, stop-on-incomplete if serial.

### E. Gaussian pKa phase loops are structurally similar

- `chemsmart/jobs/gaussian/pka.py`
  - `_run_opt_jobs` / `_run_ref_opt_jobs` / `_run_sp_jobs` / `_run_ref_sp_jobs`
  - Similar loop+propagation orchestration can share one engine-agnostic phase executor.

---

## Optimal Design (Unified Batch Runner)

### Design goals

1. One orchestration primitive for all phase execution.
2. One place to apply runner propagation and optional resource overrides.
3. One policy model: `serial`, `parallel`, `fail_fast`.
4. Keep engine-specific behavior only in job creation/settings and optional node command adaptation.

### Example implementation sketch

```python
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable, Sequence

from chemsmart.jobs.job import Job


@dataclass(frozen=True)
class PhaseSpec:
    name: str
    jobs: Sequence[Job]
    parallel: bool
    fail_fast: bool = True
    refresh_before: Callable[[], None] | None = None


class UnifiedBatchExecutor:
    """Engine-agnostic executor for pKa and other multi-phase workflows."""

    def __init__(self, parent_runner, configure_runner_for_node=None):
        self.parent_runner = parent_runner
        self.configure_runner_for_node = configure_runner_for_node

    def _attach_runner(self, job: Job, *, num_cores=None, mem_gb=None, node=None):
        runner = Job._propagate_runner(
            self.parent_runner,
            job,
            num_cores=num_cores,
            mem_gb=mem_gb,
        )
        if runner and node and self.configure_runner_for_node:
            job.jobrunner = self.configure_runner_for_node(runner, node, job)

    def _run_one(self, job: Job):
        self._attach_runner(job)
        job.run()

    def run_phase(self, spec: PhaseSpec) -> list[Exception]:
        if spec.refresh_before:
            spec.refresh_before()

        errors: list[Exception] = []

        if not spec.parallel:
            for job in spec.jobs:
                self._run_one(job)
                if spec.fail_fast and not job.is_complete():
                    break
            return errors

        max_workers = max(1, len(spec.jobs))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(self._run_one, job) for job in spec.jobs]
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    errors.append(exc)

        return errors
```

---

## Recommended Migration Sequence

1. Add a shared execution-policy resolver (single function) used by both Gaussian and ORCA pKa CLI.
2. Introduce `run_phase(...)` in shared batch infrastructure (or adjacent orchestration module).
3. Migrate `ORCApKaJob` phase loops first (highest repetition, cleanest payoff).
4. Migrate `GaussianpKaJob` to same phase API.
5. Update `process_pipeline` list handling to use unified collection execution behavior.
6. Keep engine-specific node pinning in engine batch subclasses only.

---

## Net Outcome

The `pka` branch already delivers the right foundation (`BatchJob`, `Job._propagate_runner`, runner serial flag). The remaining work is consolidation of policy decisions and phase-loop orchestration into one reusable execution layer so behavior is consistent, DRY, and easier to evolve.
