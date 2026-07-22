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

---

## 5. `SLFSubmitter._write_scheduler_options` always crashes (missing `Server.num_nodes`)

**File:** `chemsmart/settings/submitters.py:894-907`
**Test:** `tests/test_submitters_unit.py::TestSLFSubmitterBug::test_scheduler_options_crashes_on_missing_num_nodes`

```python
def _write_scheduler_options(self, f):
    ...
    f.write(f"#BSUB -nnodes {self.server.num_nodes}\n")
```

`Server` (`chemsmart/settings/server.py`) has no `num_nodes` attribute or
property anywhere — only `num_cores`, `num_gpus`, `num_hours`, etc. Any
attempt to write an LSF/SLF submission script (i.e. any job submitted with
`SCHEDULER="SLF"`) crashes with
`AttributeError: 'Server' object has no attribute 'num_nodes'`.

**Reproduce:**
```python
from io import StringIO
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import SLFSubmitter
server = Server("s", SCHEDULER="SLF", NUM_CORES=8, MEM_GB=24, NUM_GPUS=0, NUM_HOURS=24)
job = type("J", (), {"label": "j1"})()
SLFSubmitter(job=job, server=server)._write_scheduler_options(StringIO())
# AttributeError: 'Server' object has no attribute 'num_nodes'
```

**Suggested direction:** add a `num_nodes` property to `Server` (likely
reading a `NUM_NODES` kwarg with a sensible default of `1`, matching how
`num_cores`/`num_gpus` are exposed), or have `SLFSubmitter` fall back to a
hardcoded `1` node like the other submitters implicitly do.

**Related fragility (not confirmed reachable):** a few lines above,
`project_number` is only assigned inside `if user_settings is not None:`,
then read unconditionally on the next line. `user_settings` is a
module-level singleton that's always truthy in practice, so this isn't
currently reachable, but it would raise `UnboundLocalError` if that ever
changed.

---

## 6. `FUGAKUSubmitter._write_scheduler_options` always crashes (undefined `self.project`)

**File:** `chemsmart/settings/submitters.py:955-977`
**Test:** `tests/test_submitters_unit.py::TestFUGAKUSubmitterBug::test_scheduler_options_crashes_on_missing_project_attr`

```python
def _write_scheduler_options(self, f):
    if user_settings is not None:
        f.write(f'#PJM -L rscgrp={user_settings.data["RSCGRP"]}\n')
    f.write("#PJM -L node=1\n")
    f.write(f"#PJM -L elapse={self.server.num_hours}\n")
    f.write(f"#PJM --mpi proc={self.server.num_cores}\n")
    f.write(f"#PJM -g {self.project}\n")
    ...
```

`self.project` is never assigned anywhere in `FUGAKUSubmitter` or the base
`Submitter` class. Any attempt to write a FUGAKU submission script crashes
with `AttributeError: 'FUGAKUSubmitter' object has no attribute 'project'`
(after first requiring `user_settings.data["RSCGRP"]` to be set, which
raises its own `KeyError` if missing).

**Reproduce:**
```python
from io import StringIO
from unittest.mock import patch
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import FUGAKUSubmitter, user_settings
server = Server("s", SCHEDULER="FUGAKU", NUM_CORES=8, MEM_GB=24, NUM_HOURS=24)
job = type("J", (), {"label": "j1"})()
with patch.object(user_settings, "data", {"RSCGRP": "small"}):
    FUGAKUSubmitter(job=job, server=server)._write_scheduler_options(StringIO())
# AttributeError: 'FUGAKUSubmitter' object has no attribute 'project'
```

**Suggested direction:** likely meant to read a project/group ID from user
settings (mirrors `PROJECT` used by `PBSSubmitter`/`SLURMSubmitter`) —
probably should be `user_settings.data["PROJECT"]` or a dedicated
`RSCGRP`-adjacent key, not `self.project`.

---

## 7. `Server.register()` always crashes for a fresh instance

**File:** `chemsmart/settings/server.py:296-310`,
`chemsmart/utils/mixins.py:1483-1501`
**Test:** `tests/test_server_class_unit.py::TestServerRegisterBug::test_register_crashes_on_unrelated_registry_entries`

```python
def register(self):
    # if server already in registry, pass
    if self in Server._REGISTRY:
        return self
    Server._REGISTRY.append(self)
    return self
```

`RegistryMixin`'s metaclass (`RegistryMeta`, in `chemsmart/utils/mixins.py`)
initializes `_REGISTRY` as a class attribute the *first* time any
`RegistryMixin`-based class is defined, then every subsequent
`RegistryMixin` subclass — `Server`, `Executable`, `Submitter`, and all of
*their* subclasses — inherits that same single list via normal Python
attribute lookup (`hasattr(cls, "_REGISTRY")` is `True` for all of them, so
a fresh list is never created per-hierarchy). The metaclass populates it
with **classes** (not instances) across every hierarchy:

```pycon
>>> from chemsmart.settings.server import Server
>>> [x.__name__ for x in Server._REGISTRY]
['Executable', 'GaussianExecutable', 'ORCAExecutable', 'NCIPLOTExecutable',
 'Submitter', 'PBSSubmitter', 'SLURMSubmitter', 'SLFSubmitter',
 'FUGAKUSubmitter', 'Server', 'YamlServerSettings', 'SLURMServer',
 'PBSServer', 'LSFServer', 'SGE_Server']
```

`Server.register()` then does `if self in Server._REGISTRY`, appending a
**Server instance** to a list of classes belonging to unrelated
hierarchies. Python's `in` short-circuits to `True` only on `is`
(identity); for every other element it falls back to `Server.__eq__`
(`self.name == other.name`). The first non-identical element it compares
against — e.g. the `Executable` class — has no `.name` attribute at all,
so this raises `AttributeError: type object 'Executable' has no attribute
'name'` for any freshly constructed `Server`, before `register()` can ever
succeed.

**Reproduce:**
```python
from chemsmart.settings.server import Server
Server("myserver").register()
# AttributeError: type object 'Executable' has no attribute 'name'
```

**Suggested direction:** give `Server` (or `RegistryMixin` subclasses in
general) their own per-hierarchy registry rather than sharing one global
list across unrelated class families — e.g. initialize `_REGISTRY` keyed
by the immediate root class, or give `Server` a dedicated
`_INSTANCE_REGISTRY` list separate from the metaclass's class-level
`_REGISTRY`, since the two are conceptually different (registered *classes*
for dispatch vs. registered *instances* for caching/dedup).

---

## 8. `SLURMServer()` always crashes (wrong kwarg name to parent `__init__`)

**File:** `chemsmart/settings/server.py:865-872`
**Test:** `tests/test_yaml_server_settings_unit.py::TestSLURMServerBug::test_construction_raises_type_error`

```python
class SLURMServer(YamlServerSettings):
    NAME = "SLURM"
    SCHEDULER_TYPE = "SLURM"

    def __init__(self, **kwargs):
        super().__init__(filename=f"{self.NAME}.yaml", **kwargs)
```

`YamlServerSettings.__init__(self, name, **kwargs)` (and `Server.__init__`
above it) take a positional/`name` argument — there is no `filename`
parameter anywhere in the chain. Compare with the sibling classes, which
pass it correctly:

```python
class PBSServer(YamlServerSettings):
    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)   # positional — correct
```

So `SLURMServer(...)` always raises
`TypeError: YamlServerSettings.__init__() missing 1 required positional argument: 'name'`,
while `PBSServer`/`LSFServer`/`SGE_Server` construct fine.

**Reproduce:**
```python
from chemsmart.settings.server import SLURMServer
SLURMServer(NUM_CORES=8)
# TypeError: YamlServerSettings.__init__() missing 1 required positional argument: 'name'
```

**Impact:** Currently limited — `Server.from_scheduler_type()` never calls
`SLURMServer()` directly; it dispatches through
`server_cls.from_servername(scheduler_type)`, which loads
`YamlServerSettings` via `ServerSettingsManager` instead. So SLURM
autodetection itself doesn't hit this path today. But `SLURMServer` is
public API and any direct instantiation (e.g. from a future caller, or a
test) crashes.

**Suggested direction:** change to `super().__init__(self.NAME, **kwargs)`,
matching `PBSServer`/`LSFServer`/`SGE_Server`.

---

## 9. `GaussianDIASJob._sample_molecules` duplicates the endpoint molecule

**File:** `chemsmart/jobs/gaussian/dias.py:183-200`
**Test:** `tests/test_gaussian_dias_job_unit.py::TestSampleMolecules::test_samples_every_n_points_and_appends_last_again`

```python
def _sample_molecules(self, molecules):
    filtered_molecules = molecules[0 :: self.every_n_points]
    if (self.num_molecules - 1) / self.every_n_points != 0:
        filtered_molecules.append(molecules[-1])
    return filtered_molecules
```

The guard uses true division (`/`) instead of modulo (`%`). The evident
intent — based on the surrounding docstring ("ensuring the last molecule
is always included") — was to append the final molecule **only when the
slice `[0::every_n_points]` doesn't already land on it**, i.e. `(num_molecules
- 1) % every_n_points != 0`. Instead, `(num_molecules - 1) / every_n_points`
is a float that is `!= 0` for essentially any real trajectory (it's only
`0` when `num_molecules == 1`), so the last molecule is **unconditionally
appended a second time** regardless of whether the slice already included
it.

**Concrete example:** 5 molecules, `every_n_points=2`. The slice
`molecules[0::2]` is `[m0, m2, m4]`, which already includes the endpoint
`m4`. The buggy guard still evaluates `(5-1)/2 = 2.0 != 0` → `True`, so
`m4` is appended again, producing `[m0, m2, m4, m4]` — one point processed
(and later run as a full Gaussian job) twice, and the reported "sampled
count" is off by one from what a caller would expect.

**Impact:** Every IRC-mode DI-AS job (`GaussianDIASJob.all_molecules_jobs`
/ `.fragment1_jobs` / `.fragment2_jobs` when `mode="irc"`) submits one
redundant duplicate calculation for the final trajectory point, wasting
compute and duplicating a row in the DI-AS energy analysis output.

**Reproduce:**
```python
job._sample_molecules([m0, m1, m2, m3, m4])  # every_n_points=2
# [m0, m2, m4, m4]  -- m4 present twice
```

**Suggested direction:** change the guard to modulo:
`if (self.num_molecules - 1) % self.every_n_points != 0:`.
