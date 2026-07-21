# Bugs found while adding test coverage

Found while writing direct unit tests for job classes (see commits on
`improve_coverage`). Each bug is documented with a reproducing test in the
listed test file. Not fixed here — tracked for a separate branch.

---

## 1. `thermochemistry boltzmann` CLI command always crashes

**Files:** `chemsmart/cli/thermochemistry/boltzmann.py:40-46`,
`chemsmart/jobs/thermochemistry/job.py:73`
**Test:** `tests/test_boltzmann_job_unit.py::TestBoltzmannJobFromFiles::test_from_files_requires_filename_workaround`

`BoltzmannAverageThermochemistryJob.__init__` forwards `**kwargs` to
`ThermochemistryJob.__init__`, whose parent unconditionally does:

```python
# chemsmart/jobs/thermochemistry/job.py:73
if not filename.endswith((".log", ".out")):
```

`filename` defaults to `None`, so this raises
`AttributeError: 'NoneType' object has no attribute 'endswith'` unless a
`filename` kwarg is explicitly supplied.

The CLI command that constructs this job never supplies one:

```python
# chemsmart/cli/thermochemistry/boltzmann.py:40-46
boltzmann_thermochemistry = BoltzmannAverageThermochemistryJob(
    files=files,
    energy_type=energy_type_for_weighting,
    outputfile=outputfile,
    settings=job_settings.copy(),
    skip_completed=skip_completed,
)
```

**Impact:** `chemsmart sub thermochemistry boltzmann ...` cannot succeed on
any real invocation today.

**Reproduce:**
```python
from chemsmart.jobs.thermochemistry.boltzmann import BoltzmannAverageThermochemistryJob
BoltzmannAverageThermochemistryJob(files=["a.log", "b.log"])
# AttributeError: 'NoneType' object has no attribute 'endswith'
```

**Suggested direction:** either guard the extension check in
`ThermochemistryJob.__init__` for `filename is None`, or have the CLI/
`BoltzmannAverageThermochemistryJob` pass a real `filename` (e.g. the first
file in `files`).

**Related dead code:** even when a `filename` workaround is supplied,
`BoltzmannAverageThermochemistryJob.__init__`'s own label-generation block
(common prefix of `files`, suffixed `_boltzmann_avg_by_{energy_type}`) never
runs, because `ThermochemistryJob.__init__` already sets `self.label` from
`filename` first. See
`tests/test_boltzmann_job_unit.py::TestBoltzmannJobConstruction::test_label_actually_derives_from_filename_not_files`.

---

## 2. `GaussianJob.from_jobtype("g16", ..., jobrunner=<provided>)` wrongly raises

**File:** `chemsmart/jobs/gaussian/job.py:410-435`
**Test:** `tests/test_gaussian_job_base_unit.py::TestGaussianJobFactories::test_from_jobtype_g16_with_explicit_jobrunner_is_buggy`

```python
elif jobtype.lower() == "g16":
    # Create jobrunner if not provided
    if jobrunner is None:
        jobrunner = JobRunner.from_job(...)

        return GaussianGeneralJob(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
    else:
        raise ValueError(f"Invalid job type: {jobtype}")
```

The `return GaussianGeneralJob(...)` and the `else: raise ValueError(...)`
are both nested inside `if jobrunner is None:` — not at the level of the
`elif jobtype.lower() == "g16":` block. So whenever a caller already has a
jobrunner in hand (a common pattern — see how parent jobs propagate their
runner to children elsewhere in the codebase) and calls
`from_jobtype("g16", ..., jobrunner=<runner>)`, it falls into the `else`
branch and raises `ValueError: Invalid job type: g16`, even though "g16" is
a valid, recognized job type.

**Reproduce:**
```python
from unittest.mock import MagicMock
from chemsmart.jobs.gaussian.job import GaussianJob
GaussianJob.from_jobtype("g16", molecule=some_molecule, settings=some_settings, jobrunner=MagicMock())
# ValueError: Invalid job type: g16
```

**Suggested direction:** dedent the `return GaussianGeneralJob(...)` (and
drop/relocate the `else: raise ValueError`) so it executes regardless of
whether a jobrunner was already supplied.

---

## 3. `GaussianJob.from_jobtype(<unrecognized type>, ...)` silently returns `None`

**File:** `chemsmart/jobs/gaussian/job.py:357-435`
**Test:** `tests/test_gaussian_job_base_unit.py::TestGaussianJobFactories::test_from_jobtype_invalid_silently_returns_none`

Closely related to bug #2. The only `raise ValueError("Invalid job type...")`
in this method lives inside the `elif jobtype.lower() == "g16":` branch (see
above). There is no top-level `else` for jobtypes that are not
`"opt"`/`"com"`/`"g16"` at all, so calling e.g.
`GaussianJob.from_jobtype("bogus", ...)` falls through every branch and the
function implicitly returns `None` instead of raising.

Contrast with the equivalent `ORCAJob.from_jobtype` in
`chemsmart/jobs/orca/job.py:292-403`, which correctly has a top-level
`else: raise ValueError(f"Invalid job type: {jobtype}")` and behaves as
expected (see `tests/test_orca_job_base_unit.py::TestORCAJobFactories::test_from_jobtype_invalid_raises`).

**Reproduce:**
```python
from chemsmart.jobs.gaussian.job import GaussianJob
result = GaussianJob.from_jobtype("bogus", molecule=some_molecule, settings=some_settings)
print(result)  # None, no exception
```

**Suggested direction:** add a top-level `else: raise ValueError(f"Invalid job type: {jobtype}")` at the end of the `if`/`elif` chain (mirroring
`ORCAJob.from_jobtype`), and fix bug #2's indentation so the "g16" branch's
own `raise` (if kept) doesn't shadow the real error for valid jobtypes.

---

## 4. `GaussianComJob.from_filename(...)` always crashes

**File:** `chemsmart/jobs/gaussian/job.py:463-528`
**Test:** `tests/test_gaussian_job_base_unit.py::TestGaussianComJob::test_from_filename_crashes_on_none_molecule`

```python
@classmethod
def from_filename(
    cls, filename, settings=None, label=None, jobrunner=None, **kwargs
):
    ...
    return cls(
        molecule=None,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
```

`from_filename` only reads route/settings info from the `.com` file (via
`Gaussian16Input`/`GaussianJobSettings.from_filepath`) — it never extracts
molecular coordinates — and always constructs the job with `molecule=None`.
But the parent `GaussianJob.__init__` unconditionally requires
`isinstance(molecule, Molecule)`:

```python
# chemsmart/jobs/gaussian/job.py:77-80
if not isinstance(molecule, Molecule):
    raise ValueError(
        f"Molecule must be instance of Molecule for {self}, but is {molecule} instead!"
    )
```

So `GaussianComJob.from_filename(...)` cannot currently succeed at all.

**Note:** the actual `chemsmart sub gaussian ... com` CLI command
(`chemsmart/cli/gaussian/com.py`) does **not** use this factory method — it
builds `GaussianComJob` via the normal constructor with a real `Molecule`
obtained from the CLI's own file-loading logic, so the CLI path is
unaffected. This bug only affects direct use of the
`GaussianComJob.from_filename` classmethod.

**Reproduce:**
```python
from chemsmart.jobs.gaussian.job import GaussianComJob
GaussianComJob.from_filename(filename="some_input.com")
# ValueError: Molecule must be instance of Molecule for ..., but is None instead!
```

**Suggested direction:** either read the molecule from the `.com` file (e.g.
via `Molecule.from_filepath(filename)`) and pass it through, or relax the
parent's molecule-type check for job types that are known to not need one
(mirrors how `ThermochemistryJob` allows `molecule=None`).
