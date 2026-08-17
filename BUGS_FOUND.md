# Bugs found while adding test coverage

Found while writing direct unit tests for job classes (see commits on
`improve_coverage`). Each bug is documented with a reproducing test in the
listed test file. Not fixed here — tracked for a separate branch.

---

## 1. `thermochemistry boltzmann` CLI command always crashes

**Files:** `chemsmart/cli/thermochemistry/boltzmann.py:40-46`,
`chemsmart/jobs/thermochemistry/job.py:72-73`
**Test:** `tests/test_boltzmann_job_unit.py::TestBoltzmannJobFromFiles::test_from_files_requires_filename_workaround`

**Update:** `ThermochemistryJob.__init__` now explicitly guards against
`filename is None` (upstream `main` change, merged into this branch after
this entry was first written):

```python
# chemsmart/jobs/thermochemistry/job.py:72-73
if filename is None:
    raise ValueError("'filename' must be provided.")
```

This is exactly the fix suggested below, so the confusing
`AttributeError: 'NoneType' object has no attribute 'endswith'` no longer
occurs — the failure is now a clear `ValueError: 'filename' must be
provided.` raised earlier, directly from `__init__`.

**Fixed:** the CLI's `boltzmann()` command now passes
`filename=files[0] if files else None` to
`BoltzmannAverageThermochemistryJob`, so
`chemsmart sub thermochemistry boltzmann -f a.log -f b.log ...` no longer
crashes. See `tests/test_thermochemistry_cli.py::TestThermochemistryBoltzmannCommand`
for CLI-level regression tests covering this fix.

**Impact (historical):** `chemsmart sub thermochemistry boltzmann ...`
could not succeed on any real invocation prior to this fix.

**Reproduce (historical, prior to fix):**
```python
from chemsmart.jobs.thermochemistry.boltzmann import BoltzmannAverageThermochemistryJob
BoltzmannAverageThermochemistryJob(files=["a.log", "b.log"])
# ValueError: 'filename' must be provided.
```

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

---

## 10. `GaussianDIASLogFolder` DI-AS analysis crashes on normal filenames

**Files:** `chemsmart/analysis/dias.py:439`,
`chemsmart/utils/repattern.py:135`
**Test:** `tests/test_dias_analysis_unit.py::TestGaussianDIASBugFullMoleculeGroupIndex::test_full_molecule_scan_crashes_without_dias_in_filename`

`GaussianDIASLogFolder._get_all_files_full_molecule` matches candidate
filenames against
`gaussian_dias_filename_point_without_fragment_without_reactant`:

```python
# chemsmart/utils/repattern.py:135
gaussian_dias_filename_point_without_fragment_without_reactant = (
    r"(?:(?!.*_r[12](?:_|\.)).*dias_p(\d+)(?:_((?:(?!f\d).)+))?\.log)"
    r"|(?:(?!.*_r[12](?:_|\.)).*_p(\d+)(?:_((?:(?!f\d).)+))?\.log)"
)
```

This is a two-alternative regex. The **first** alternative requires the
literal substring `dias_p` in the filename and puts the point number in
capture group 1. The **second** (fallback, plain `_p<digits>`)
alternative puts the point number in group 3 instead, leaving group 1
as `None`. But the calling code always reads group 1:

```python
# chemsmart/analysis/dias.py:439
match = full_molecule_pattern.match(filename)
if match:
    number_after_p = int(match.group(1))  # None when 2nd alt matched
```

**Concrete example:** a full-molecule output named `rxn_p1.log` (no
"dias" substring) matches via the second alternative —
`match.groups() == (None, None, '1', None)` — so
`int(match.group(1))` raises
`TypeError: int() argument must be a string, a bytes-like object or a
real number, not 'NoneType'`.

**Impact:** Real `GaussianDIASJob` output files are named
`f"{label}_p{i}.log"` (see `chemsmart/jobs/gaussian/dias.py:316`), with
no "dias" substring unless the user's job `label` happens to contain
the word "dias". So `GaussianDIASLogFolder` — used for the
post-processing/plotting analysis step, separate from the job
submission itself — crashes on essentially any normally-named DI-AS
job folder. (`ORCADIASOutFolder`'s equivalent pattern,
`orca_dias_filename_point_without_fragment`, has only one alternative
and does not have this bug.)

**Reproduce:**
```python
import re
from chemsmart.utils.repattern import (
    gaussian_dias_filename_point_without_fragment_without_reactant as p,
)
re.compile(p).match("rxn_p1.log").groups()
# (None, None, '1', None)  -- group(1) is None
```

**Suggested direction:** either read whichever group is non-`None`
(e.g. `next(g for g in match.groups() if g is not None)`), or collapse
the regex to a single alternative with an optional non-capturing
`dias_` prefix so the point number always lands in the same group.

---

## 11. `GaussianUVVISJob.__init__` always crashes

**File:** `chemsmart/jobs/gaussian/uvvis.py:59`
**Test:** `tests/test_gaussian_uvvis_job_unit.py::TestGaussianUVVISJobInitCrash::test_init_always_raises_typeerror`

```python
def __init__(self, folder, atoms, settings):
    super().__init__(folder=folder, atoms=atoms, settings=settings)
```

`GaussianUVVISJob` extends `GaussianJob`, whose `__init__` signature is
`__init__(self, molecule, settings=None, label=None, jobrunner=None,
**kwargs)` — it takes `molecule`, not `atoms`, and `molecule` has no
default. Passing `atoms=` instead means `molecule` is never supplied, so
every construction fails immediately:

```
TypeError: GaussianJob.__init__() missing 1 required positional argument: 'molecule'
```

**Impact:** `GaussianUVVISJob` cannot be instantiated at all today. It is
already flagged `# TODO: incomplete job` in the source and is not wired
into any CLI command (`grep` finds no caller outside
`chemsmart/jobs/gaussian/runner.py`'s `TYPE` registry and the package
`__init__.py` export), so this doesn't affect any currently-reachable
user workflow — but the class is public API and crashes on the most
basic possible use.

**Suggested direction:** rename the `atoms` parameter/kwarg to
`molecule` to match `GaussianJob.__init__`.

---

## 12. `ORCAQMMMJobSettings` intermediate-layer validation has a silent gap

**File:** `chemsmart/jobs/orca/settings.py:1898-1965`
**Test:** `tests/test_orca_qmmm_neb_settings_unit.py::TestQMMMIntermediateValidation::test_intermediate_params_for_non_qm2_jobtype_without_low_level_is_silently_accepted`

`_validate_intermediate_parameters` is meant to reject intermediate-level
(QM2) parameters supplied for a job type that doesn't use a QM2 layer
(anything other than `QM/QM2` or `QM/QM2/MM`):

```python
if not requires_qm2 and has_intermediate_params:
    provided_params = [...]
    if self.low_level_method is not None:
        raise ValueError(
            f"Job type '{self.jobtype}' does not require QM2 (intermediate) layer, but "
            ...
        )
```

The `provided_params` list is built unconditionally but the `raise` is
additionally gated on `self.low_level_method is not None`. So a plain
`QMMM` job (or any non-QM2 job type) that specifies
`intermediate_level_functional`/`intermediate_level_basis`/
`intermediate_level_method` **without** also setting `low_level_method`
passes construction silently, even though the error message's own logic
("does not require QM2 layer, but intermediate-level parameters were
provided") says it should be rejected.

**Reproduce:**
```python
from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings
s = ORCAQMMMJobSettings(
    jobtype="QMMM",
    intermediate_level_functional="B3LYP",
    intermediate_level_basis="def2-SVP",
    charge_total=0,
    mult_total=1,
)
s.intermediate_level_functional  # "B3LYP" -- no error raised
```

**Impact:** A caller who misconfigures a `QMMM` job with a leftover/typo'd
intermediate-level layer (but no `low_level_method`) gets no validation
error; the unused intermediate parameters are simply ignored by
`_get_level_of_theory_string`/`_write_qmmm_block` rather than flagged.

**Suggested direction:** drop the `self.low_level_method is not None`
condition from the `raise` so the check fires whenever
`not requires_qm2 and has_intermediate_params`, matching the error
message's stated intent.

---

## 13. `GaussianCrestJob` empty/non-list `molecules` bypasses its own validation

**File:** `chemsmart/jobs/gaussian/crest.py:81-82`
**Test:** `tests/test_gaussian_crest_job_unit.py::TestGaussianCrestJobConstruction::test_empty_molecules_list_bypasses_validation_and_crashes`

```python
if not isinstance(molecules, list) and len(molecules) == 0:
    raise ValueError("Molecules must be a list of Molecule objects.")
```

The guard uses `and` where the intent (per the error message, which
complains about *either* condition) is clearly `or`. Consequences:

* An **empty list** — `not isinstance([], list)` is `False` — makes the
  whole condition `False`, so the intended `ValueError` never fires.
  Execution instead reaches `molecules[0]` a few lines later and crashes
  with an unrelated `IndexError: list index out of range`.
* A **non-list, non-empty** iterable (e.g. a tuple) — `not isinstance(x,
  list)` is `True` but `len(x) == 0` is `False` — also makes the
  condition `False`, so passing a tuple instead of a list is silently
  accepted rather than rejected as the message promises.

**Reproduce:**
```python
from chemsmart.jobs.gaussian.crest import GaussianCrestJob
GaussianCrestJob(molecules=[])
# IndexError: list index out of range  (not the intended ValueError)
```

**Suggested direction:** change `and` to `or`:
`if not isinstance(molecules, list) or len(molecules) == 0:`.

---

## 14. `GaussianLinkJob.backup_files` is unreachable dead code

**File:** `chemsmart/jobs/gaussian/link.py:183-193`
**Test:** `tests/test_gaussian_link_job_unit.py::TestBackupFilesIsDeadCode::test_backup_files_is_not_the_method_invoked_by_backup`

`Job.backup(**kwargs)` (the only call site for backups in the codebase)
invokes the private, underscore-prefixed hook:

```python
# chemsmart/jobs/job.py:183-187
def backup(self, **kwargs):
    self._backup_files(**kwargs)
```

`GaussianLinkJob` defines an IRC-aware override named `backup_files`
(no leading underscore) instead of `_backup_files`:

```python
# chemsmart/jobs/gaussian/link.py:183
def backup_files(self, backup_chk=False):
    if self._is_irc_job():
        irc_jobs = self._get_irc_jobs()
        for job in irc_jobs:
            self.backup_file(job.inputfile)
            ...
```

Because the method name doesn't match the hook `Job.backup()` actually
calls, `GaussianLinkJob` silently inherits `GaussianJob._backup_files`
unchanged — confirmed directly: `GaussianLinkJob._backup_files is
GaussianJob._backup_files`. The IRC-subjob-aware backup logic in
`backup_files` is never invoked by the normal job lifecycle; it only
runs if a caller happens to spell out `job.backup_files(...)` explicitly
instead of `job.backup(...)`.

**Impact:** Calling `.backup()` on an IRC-mode `GaussianLinkJob` backs up
only the link job's own (largely unused, since IRC link jobs delegate to
forward/reverse subjobs) files, not the forward/reverse subjobs'
`.com`/`.log`/`.chk` files that the override was clearly written to
protect.

**Suggested direction:** rename `backup_files` to `_backup_files` (and
match the base class's `backup_chk` keyword-argument handling via
`**kwargs` if needed for signature compatibility with `Job.backup`).

---

## 15. `NCIPLOTJob(settings=None)` always crashes

**File:** `chemsmart/jobs/nciplot/job.py:67-112`
**Test:** `tests/test_nciplot_job_unit.py::TestNCIPLOTJobSettingsDefaultCrashes::test_settings_none_always_crashes`

`settings` defaults to `None` in the signature, and is only
conditionally validated:

```python
def __init__(self, filenames=None, molecule=None, settings=None, ...):
    ...
    if settings is not None and not isinstance(settings, NCIPLOTJobSettings):
        raise ValueError(...)
    ...
    self.settings = settings.copy()  # unconditional
```

The `isinstance` check is skipped when `settings is None` (as the
default value is meant to allow), but `self.settings = settings.copy()`
a few lines later runs unconditionally regardless of whether `settings`
was ever supplied, so relying on the documented default always crashes:

```
AttributeError: 'NoneType' object has no attribute 'copy'
```

**Reproduce:**
```python
from chemsmart.jobs.nciplot.job import NCIPLOTJob
NCIPLOTJob(filenames=["a.xyz"], settings=None)
# AttributeError: 'NoneType' object has no attribute 'copy'
```

**Impact:** The CLI (`chemsmart/cli/nciplot/nciplot.py`) always builds
and passes an explicit `NCIPLOTJobSettings` instance, so this isn't
reachable through the current CLI entry point — but any other caller
(scripts, notebooks, tests) that constructs `NCIPLOTJob` without
supplying `settings` explicitly hits this immediately.

**Suggested direction:** default to `settings = settings or
NCIPLOTJobSettings()` before validating/copying, matching the pattern
other job classes use for optional settings.

---

## 16. `GaussianNCIJob.__init__` ignores its `jobrunner` argument

**File:** `chemsmart/jobs/gaussian/nci.py:32-51`
**Test:** `tests/test_gaussian_small_job_wrappers_unit.py::TestGaussianNCIJob::test_jobrunner_argument_is_ignored_and_class_used_instead`

```python
from chemsmart.jobs.runner import JobRunner

class GaussianNCIJob(GaussianJob):
    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=JobRunner,   # <-- the class, not the `jobrunner` param
            **kwargs,
        )
```

`__init__` accepts a `jobrunner` parameter but never uses it — it passes
the literal `JobRunner` class object to the parent constructor instead.
Every `GaussianNCIJob` therefore gets `self.jobrunner is JobRunner` (the
abstract base class itself, not an instance), regardless of what the
caller supplied.

**Reproduce:**
```python
job = GaussianNCIJob(molecule=mol, settings=settings, label="x", jobrunner=my_runner)
job.jobrunner is my_runner   # False
job.jobrunner is JobRunner   # True -- the class itself
```

**Impact:** Any caller-supplied jobrunner (e.g. a configured
`GaussianJobRunner` instance from `JobRunner.from_job`) is silently
discarded. Calling `.run()` on a `GaussianNCIJob` would invoke methods
on the `JobRunner` class rather than a real runner instance, which is
very unlikely to work correctly for actual job submission.

**Suggested direction:** pass through the `jobrunner` parameter like
every other Gaussian job class does: `jobrunner=jobrunner`.

---

## 17. `EnergyGrouper._merge_groups_to_target`'s "no connections" fallback is unreachable dead code

**File:** `chemsmart/jobs/grouper/energy.py:428-458`
**Test:** `tests/test_energy_grouper_unit.py::TestMergeGroupsToTargetDirect::test_merges_into_the_most_connected_group`
(see the docstring note; the fallback itself has no reachable test because it cannot execute)

```python
best_merge_idx = -1
best_connection_count = -1

for i, target_indices in enumerate(index_groups):
    if i == min_idx:
        continue
    connection_count = 0
    for src_idx in index_groups[min_idx]:
        for tgt_idx in target_indices:
            if adj_matrix[src_idx, tgt_idx]:
                connection_count += 1
    if connection_count > best_connection_count:
        best_connection_count = connection_count
        best_merge_idx = i

if best_merge_idx >= 0:
    ...
else:
    # No connections found, merge into the largest group
    ...
```

`best_connection_count` starts at `-1`, but `connection_count` starts at
`0` for every candidate group. Since `0 > -1` is always true, the very
first candidate considered always sets `best_merge_idx` to a real index
(`>= 0`). The `while len(groups) > self.num_groups` loop guarantees at
least 2 groups exist whenever this code runs, so there is always at
least one candidate group to consider — meaning `best_merge_idx` can
never remain `-1`. The `else` branch ("no connections found, merge into
the largest group") is therefore unreachable in practice, even when the
adjacency matrix is fully disconnected (verified directly: with an
all-`False` matrix, `best_merge_idx` still ends up `0`, not `-1`).

**Impact:** Purely a latent dead-code / defensive-code bug — behavior
is unaffected today because the "most connected" branch already merges
into *some* group (arbitrarily the first one iterated, when all
connection counts are tied at 0), which happens to coincide with
reasonable behavior. But the fallback was clearly intended to have a
distinct effect ("merge into the *largest* group" instead of the first
tied one) and never fires.

**Suggested direction:** initialize `best_connection_count = -1` is fine,
but change the comparison to only prefer strictly-positive connection
counts (e.g. track `best_merge_idx` separately for `connection_count >
0` cases), so ties-at-zero genuinely fall through to the "merge into
largest" branch.

**Also affects:** `chemsmart/jobs/grouper/rmsd.py`'s
`RMSDGrouper._merge_groups_to_target` (around line 382-442) has the
exact same `best_connection_count = -1` / `connection_count = 0` pattern
and is subject to the identical dead-code fallback. See
`tests/test_rmsd_grouper_base_unit.py::TestMergeGroupsToTarget::test_merge_reduces_group_count_and_prefers_connected_group`.
`TorsionFingerprintGrouper._merge_groups_to_target` in
`chemsmart/jobs/grouper/tfd.py` is unaffected — it merges purely by
group size and has no connection-counting fallback.

---

## 18. `chemsmart database export` CLI: unsupported output extension raises a raw `ValueError` instead of a clean CLI error

**File:** `chemsmart/cli/database/export.py:141-216`
**Test:** `tests/test_database_cli.py::TestExportCommandValidation::test_unsupported_output_extension_skips_cli_validation`

The `export` command validates `-o/--output`'s extension against
`{.json, .csv, .xyz, .extxyz}` implicitly, via an `if ext in
(".json",".csv"): ... elif ext in (".xyz",".extxyz"): ...` chain (lines
141-190) that raises friendly `click.UsageError`s for bad option
combinations. But if `ext` matches *neither* branch (e.g. `-o
out.txt`), both branches are skipped entirely — no error is raised
there. The only place that actually rejects an unsupported extension is
`DatabaseExporter.__init__` → `_infer_format()`, which raises a plain
`ValueError`. That constructor call (`DatabaseExporter(...)` at line
197) is *not* wrapped in a try/except — only the later
`exporter.export()` call is (lines 210-213). So an unsupported output
extension surfaces as an uncaught `ValueError` with a Python traceback,
rather than the clean `click.UsageError`/`ClickException` messages used
everywhere else in this command.

**Reproduce:**
```
chemsmart run database export -f my.db --rid abc123 -o out.txt
# Traceback (most recent call last):
#   ...
# ValueError: Unsupported output format '.txt'. Supported extensions: ...
```

**Impact:** Minor UX inconsistency — a mistyped output extension gives
a stack trace instead of a one-line usage error like every other
invalid-input case in this command.

**Suggested direction:** wrap the `DatabaseExporter(...)` construction
in the same `try/except ValueError: raise click.ClickException(...)`
used around `exporter.export()`, or validate `ext` against
`chemsmart.database.export.SUPPORTED_FORMATS` up front alongside the
existing `if`/`elif` chain.

---

## 19. `FileConverter._convert_all_files`'s per-file "unsupported type" branch is unreachable dead code

**File:** `chemsmart/io/converter.py:100-190`
**Test:** `tests/test_converter_unit.py` (see `TestConvertAllFilesDirectoryDispatch`; no test targets the dead branch itself since it cannot execute)

`_convert_all_files` validates `type` twice against the exact same set
of eight values (`log`, `com`, `gjf`, `out`, `inp`, `xyz`, `sdf`, `pdb`,
plus `cdxml`/`cdx`):

1. Once at the top (lines 100-146) to decide which folder-listing
   helper to call, with an `else: raise ValueError(...)` for anything
   else.
2. Again per-file inside the loop (lines 152-190) to decide which file
   class to instantiate, with an identical `else: raise
   ValueError(f"File type {type} is not supported.")` at line 190.

Since `type` is a local parameter that isn't reassigned between the two
checks, and both `elif` chains list precisely the same values, reaching
the loop at all already proves `type` matched one of the first chain's
branches — so the second `else` at line 190 can never execute.

**Impact:** None today — purely redundant defensive code with no
behavioral effect, since the outer validation always catches an
unsupported type before the loop is ever entered.

**Suggested direction:** harmless to leave, but could be simplified by
dropping the second validation (or converting it to an `assert` /
`AssertionError` documenting the invariant instead of a duplicated
user-facing `ValueError`).

---

## 20. `GaussianpKaJobSettings.build_gaussian_pka_settings` crashes when `opt_settings` carries solvent info not overridden by `shared`

**File:** `chemsmart/jobs/gaussian/settings.py:1187-1252`
**Test:** `tests/test_gaussian_settings_route_strings_unit.py::TestBuildGaussianPkaSettings::test_solvent_settings_from_opt_settings_crashes`

```python
opt_kwargs = {
    key: value
    for key, value in vars(opt_settings).items()
    if key in gs_params and value is not None and key not in pka_kwargs
}
...
solvent_model = _first_non_none(
    pka_kwargs.get("solvent_model"),
    getattr(opt_settings, "solvent_model", None),
    ...,
    "SMD",
)
...
pka_kwargs["solvent_model"] = solvent_model
pka_kwargs["solvent_id"] = solvent_id

return cls(proton_index=proton_index, **pka_kwargs, **opt_kwargs)
```

`opt_kwargs` is built by copying every non-`None` attribute off
`opt_settings` that isn't already a key in `pka_kwargs` **at that
point** — but `solvent_model`/`solvent_id` aren't added to `pka_kwargs`
until several lines later. So when `opt_settings.solvent_model` (or
`.solvent_id`) is set and `shared` doesn't specify its own value, the
resolved fallback ends up in `pka_kwargs["solvent_model"]` *and* the
verbatim `opt_settings` value ends up in `opt_kwargs["solvent_model"]`
— both dicts carrying the same key when unpacked into the same `cls(...)`
call.

**Reproduce:**
```python
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings, GaussianpKaJobSettings,
)
opt_settings = GaussianJobSettings(
    functional="b3lyp", basis="sto-3g",
    solvent_model="PCM", solvent_id="dmso",
)
GaussianpKaJobSettings.build_gaussian_pka_settings(
    proton_index=2, shared={}, opt_settings=opt_settings,
)
# TypeError: ...got multiple values for keyword argument 'solvent_model'
```

**Impact:** Any pKa CLI invocation where the optimization step's
settings already specify a solvent (the common case — pKa calculations
are almost always run in solution) and the pKa-specific CLI flags don't
redundantly repeat `--solvent-model`/`--solvent-id` will crash instead
of inheriting the optimization step's solvent.

**Suggested direction:** move the `pka_kwargs["solvent_model"] = ...` /
`pka_kwargs["solvent_id"] = ...` assignments to *before* `opt_kwargs` is
computed, so the `key not in pka_kwargs` filter correctly excludes them
from `opt_kwargs` once resolved.

---

## 21. `GaussianQMMMJobSettings._get_charge_and_multiplicity`'s per-sub-level fill branches are dead code

**File:** `chemsmart/jobs/gaussian/settings.py:2791-2944`
**Test:** `tests/test_gaussian_qmmm_settings_unit.py::TestChargeAndMultiplicityErrors`
(see the docstring notes on `test_three_layer_fills_from_intermediate_when_only_real_and_int_set`;
the unreachable branches themselves have no test since they cannot execute)

`_get_charge_and_multiplicity` builds a flat list of
charge/multiplicity pairs for each ONIOM sub-level, then picks one of
several "fill in the blanks" branches depending on which suffix of the
list is entirely `None`:

```python
# 3-layer case
charge_and_multiplicity_list = [
    real_low_charge, real_low_multiplicity,      # [0:2]
    int_med_charge, int_med_multiplicity,        # [2:4]
    int_low_charge, int_low_multiplicity,        # [4:6]  <- always == [2:4]
    model_high_charge, model_high_multiplicity,  # [6:8]
    model_med_charge, model_med_multiplicity,    # [8:10] <- always == [6:8]
    model_low_charge, model_low_multiplicity,    # [10:12] <- always == [6:8]
]
if all(v is None for v in list[2:]): ...      # reachable
elif all(v is None for v in list[4:]): ...    # DEAD
elif all(v is None for v in list[6:]): ...    # reachable
elif all(v is None for v in list[8:]): ...    # DEAD
elif all(v is None for v in list[10:]): ...   # DEAD
elif all(v is not None for v in list): pass   # reachable
else: raise ValueError(...)                   # reachable
```

Just above this method, `int_low_charge`/`int_low_multiplicity` are
assigned the *exact same* `self.charge_intermediate`/`self.mult_intermediate`
values as `int_med_charge`/`int_med_multiplicity` — there is no
independent "intermediate, low level-of-theory" constructor parameter,
only one `charge_intermediate`/`mult_intermediate` (aliased from legacy
`int_charge`/`int_multiplicity`) pair. Likewise `model_high_charge` ==
`model_med_charge` == `model_low_charge` (and the multiplicity
equivalents), all sourced from the single `charge_high`/`mult_high`
pair (aliased from legacy `model_charge`/`model_multiplicity`).

Because of this, list positions `[2:4]` and `[4:6]` are always
identical, and positions `[6:8]`, `[8:10]`, `[10:12]` are always
identical. That makes it impossible for `list[4:]` to be "all None"
without `list[2:]` *also* being all None (in which case the first,
higher-priority branch already fired) — and likewise for `list[8:]`
and `list[10:]` relative to `list[6:]`. The exact same redundancy
exists in the 2-layer variant (`model_high_charge` == `model_low_charge`
always, making its `list[4:]`-all-None branch at line 2862 equally
unreachable).

**Reproduce (informal):** no input can make `list[4:]`, `list[8:]`, or
`list[10:]` all-`None` while the branch immediately above it was
`False`, because the values being checked are literally aliases of
values already checked by that branch. Verified directly: constructing
`GaussianQMMMJobSettings` with only `model_charge`/`model_multiplicity`
set (no `int_charge`/`int_multiplicity`) hits the `list[6:]`-all-None
branch first, not `list[8:]`, since positions 8-11 being None forces
positions 6-7 (the identical `model_high` values) to also be None.

**Impact:** None today — these branches were presumably intended to
support independently specifying charge/multiplicity for each of the
three model-system sub-levels (high/medium/low), matching the detailed
docstring on `_get_charge_and_multiplicity`, but the constructor never
actually exposes separate parameters for `int_low_*` vs `int_med_*` or
for `model_med_*`/`model_low_*` vs `model_high_*`, so those branches
can never be exercised as designed.

**Suggested direction:** either add the missing independent
constructor parameters (e.g. `charge_intermediate_low`,
`charge_model_medium`, `charge_model_low`, etc.) so each sub-level can
truly be set independently, or — if the simplification is intentional
— remove the dead branches and update the docstring to reflect that
only "real", "intermediate", and "model" (not six independent
sub-levels) can be specified.

## 22. `BasePreprocessor._get_max_bonding_capacity`'s tuple-handling branch is unreachable dead code

**Location:** `chemsmart/jobs/iterate/iterate.py`, `_get_max_bonding_capacity` (lines 58-71).

```python
default_valence = periodic_table.GetDefaultValence(atomic_num)

if isinstance(default_valence, tuple):
    return max(default_valence)
return default_valence
```

The comment above this code claims `GetDefaultValence` "returns a
tuple of possible valences" for some elements, but RDKit's
`GetDefaultValence` (unlike `GetValenceList`) always returns a single
`int` — never a tuple — for every element in the periodic table.

**Reproduce (informal):** iterated `Chem.GetPeriodicTable().GetDefaultValence(z)`
for every atomic number `z` from 1 to 99 and confirmed the return type
is `int` in every case; `isinstance(default_valence, tuple)` is `False`
for all of them, so the `return max(default_valence)` line can never
execute.

**Impact:** None observed — the `return default_valence` fallback is
correct and always taken, so `_has_available_bonding_position` computes
the right max-bonding-capacity regardless. This is purely dead code
left over from confusing `GetDefaultValence` with `GetValenceList`
(which does return a tuple/list of allowed valences).

**Suggested direction:** remove the `isinstance`/`tuple` branch (and
the stale comment) and just return `default_valence` directly — or,
if the intent was genuinely to allow for multi-valence elements, switch
to `GetValenceList` and take its max instead.

## 23. `SkeletonPreprocessor._dfs_collect_branch`'s `node == excluded` guard is unreachable dead code

**Location:** `chemsmart/jobs/iterate/iterate.py`, `_dfs_collect_branch`
(lines 343-380).

```python
while stack:
    node = stack.pop()
    if node in visited or node == excluded:
        continue
    visited.add(node)
    branch_atoms.append(node)

    for neighbor in graph.neighbors(node):
        if neighbor not in visited and neighbor != excluded:
            stack.append(neighbor)
```

`excluded` can only ever end up on `stack` in one of two ways: as the
initial `start` value, or via the neighbor-expansion loop. The only
caller, `_find_non_skeleton_branches`, always invokes this with
`start=neighbor` and `excluded=self.link_index`, where `neighbor` is by
construction a graph-neighbor of `link_index` and therefore never
equal to it (no self-loops). The neighbor-expansion loop itself already
filters with `neighbor != excluded` before pushing onto the stack. So
`excluded` is never pushed onto `stack` by either path, meaning the
`node == excluded` half of the guard on line 370 can never be `True` —
`node in visited` is the only condition that can ever fire `continue`.

**Reproduce (informal):** constructed a 3-membered carbon ring and
called `SkeletonPreprocessor._find_non_skeleton_branches()` with the
ring's link atom — even in this "worst case" topology (a cycle
containing the link atom, where DFS revisits nodes on the way back
around the ring), the link atom index never reaches line 370 with
`node in visited` being `False`, confirming the excluded-check never
independently triggers a `continue`.

**Impact:** None — the neighbor-expansion filter already guarantees
correctness; this is purely redundant/dead defensive code.

**Suggested direction:** drop the `or node == excluded` clause from
line 370 (and, if desired, the `excluded` parameter entirely, since
nothing pushes it onto the stack) since the neighbor-expansion filter
already fully excludes it.

## 24. `compute_pka_thermochemistry`'s nested `get_species_thermo`'s `filepath is None` guard is unreachable dead code

**Location:** `chemsmart/cli/pka.py`, `compute_pka_thermochemistry`
(lines 273-360), nested helper `get_species_thermo` (lines 306-350).

```python
def get_species_thermo(filepath, name):
    if filepath is None:
        return None
    thermo = Thermochemistry(filename=filepath, **thermo_kwargs)
    ...

if ha_file is not None:
    results["HA"] = get_species_thermo(ha_file, "HA")
if a_file is not None:
    results["A"] = get_species_thermo(a_file, "A-")
if href_file is not None:
    results["HRef"] = get_species_thermo(href_file, "HRef")
if ref_file is not None:
    results["Ref"] = get_species_thermo(ref_file, "Ref-")
```

`get_species_thermo` is a purely local nested function with exactly
four call sites, and every one of them is already guarded by an
`if X_file is not None:` check before the call. So `filepath` can never
be `None` inside `get_species_thermo` — the function's own
`if filepath is None: return None` guard (line 307-308) can never be
taken.

**Reproduce (informal):** grepped all call sites of
`get_species_thermo` within `compute_pka_thermochemistry`; each of the
four calls passes `ha_file`/`a_file`/`href_file`/`ref_file` directly
from inside a block that already tested that same value `is not None`.

**Impact:** None — the outer guards already produce the correct
"omit this species from the results dict" behavior; this is purely
redundant/dead defensive code inside the nested helper.

**Suggested direction:** drop the `if filepath is None: return None`
guard from `get_species_thermo`, since every caller already ensures
`filepath` is not `None` before invoking it.

## 25. `gaussian` CLI group crashes computing the default label for a filename-less (PubChem-only) job with no `-l`/`-a`

**Location:** `chemsmart/cli/gaussian/gaussian.py`, the `gaussian()`
group callback's label-resolution block (lines ~782-794).

```python
if label is None and append_label is None:
    label = os.path.splitext(os.path.basename(filename))[0]
    if is_chemsmart_db:
        if structure_id is not None:
            label = f"{label}_SID-{structure_id}"
        elif record_id is not None:
            label = f"{label}_RID-{record_id}"
        elif record_index is not None:
            label = f"{label}_RI-{record_index}"
    if filename:
        label = os.path.splitext(os.path.basename(filename))[0]
    else:
        label = "output"
    if ctx.invoked_subcommand:
        label = f"{label}_{ctx.invoked_subcommand}"
```

Two problems in this one block:

1. **Crash when `filename` is `None`.** The very first line
   unconditionally calls `os.path.basename(filename)` — even though the
   `if filename: ... else: label = "output"` guard a few lines down
   exists specifically to handle `filename is None` (e.g. a
   `--pubchem`-only job with no `-f`). Since `os.path.basename(None)`
   raises `TypeError` immediately, the `else: label = "output"`
   fallback (and the `ctx.invoked_subcommand` suffixing after it) is
   never reached in that case — the whole command crashes instead.
2. **Dead SID/RID/RI suffixing for chemsmart-db default labels.** Even
   when `filename` *is* given, the `is_chemsmart_db` suffix computed on
   lines 2-8 (`label = f"{label}_SID-..."` etc.) is immediately
   discarded: the very next line unconditionally recomputes
   `label = os.path.splitext(os.path.basename(filename))[0]` from
   scratch when `filename` is truthy, wiping out whatever suffix was
   just added. The analogous `append_label is not None` branch just
   above this block does **not** have this problem (it appends
   `_{append_label}` after the suffix, so the suffix survives) — this
   block is clearly meant to mirror it but recomputes the base label a
   second time by mistake.

**Reproduce:**
```
chemsmart run gaussian --pubchem 222 -c 0 -m 1 opt
# TypeError: expected str, bytes or os.PathLike object, not NoneType
```
(also reproduced directly via `tests/test_gaussian_cli.py::TestGaussianCLIGroupValidation::test_pubchem_only_without_label_crashes`)

**Impact:** Any PubChem-only Gaussian job submission (no `-f`) that
doesn't also pass `-l`/`--label` or `-a`/`--append-label` crashes
outright instead of falling back to the documented `"output"` label.
Additionally, for chemsmart-database inputs selected by `--sid`/`--rid`/
`--ri` with no explicit `-l`/`-a`, the resulting default label never
actually contains the `_SID-`/`_RID-`/`_RI-` suffix it appears to
compute, silently losing the disambiguating suffix between multiple
structures/records pulled from the same database file.

**Suggested direction:** guard the whole block with `if filename:`
before computing `label`, mirroring the `append_label` branch above
it, e.g.:
```python
if label is None and append_label is None:
    if filename:
        label = os.path.splitext(os.path.basename(filename))[0]
        if is_chemsmart_db:
            if structure_id is not None:
                label = f"{label}_SID-{structure_id}"
            elif record_id is not None:
                label = f"{label}_RID-{record_id}"
            elif record_index is not None:
                label = f"{label}_RI-{record_index}"
    else:
        label = "output"
    if ctx.invoked_subcommand:
        label = f"{label}_{ctx.invoked_subcommand}"
```

## 26. `gaussian` CLI: `click_gaussian_qmmm_options` and the group-level `qmmm` conversion block are both unreachable dead code

**Location:** `chemsmart/cli/gaussian/gaussian.py`,
`click_gaussian_qmmm_options` (lines 361-490) and the
`ctx.invoked_subcommand == "qmmm"` block inside the `gaussian()` group
callback (lines ~822-859).

`click_gaussian_qmmm_options` is a click-options decorator function
defining `-hx/-hb/-hf/-mx/...` QMMM-layer options, but it is never
applied as a decorator anywhere in the codebase (confirmed via
project-wide grep for `click_gaussian_qmmm_options`) — the actual
`qmmm` subcommand's options are defined independently in
`chemsmart/cli/gaussian/qmmm.py`. This helper is simply orphaned,
presumably superseded when `qmmm` was refactored into a subcommand
nested under each jobtype (`opt`, `ts`, `sp`, `scan`, `qrc`, `modred`
all call `create_qmmm_subcommand(<jobtype>)` in their own modules)
rather than being a sibling top-level command under `gaussian` itself.

Relatedly, the `gaussian()` group callback contains:
```python
try:
    if ctx.invoked_subcommand == "qmmm":
        ...  # convert molecules to QMMMMolecule
except Exception as exc:
    ...
```
But since `qmmm` is *always* registered as a child of a jobtype
subcommand (e.g. `opt qmmm`, `ts qmmm`) rather than a direct child of
`gaussian`, `ctx.invoked_subcommand` at the `gaussian` group's own
callback is always one of `opt`/`ts`/`sp`/`scan`/`qrc`/`modred`/etc. —
**never** `"qmmm"` itself (Click's `ctx.invoked_subcommand` only
reflects the immediate child command). So this condition can never be
`True`, and the whole `QMMMMolecule` conversion body is unreachable
group-level code; QMMM molecule conversion, if it happens at all, must
be (and is) handled independently within `qmmm.py`'s own subcommand
callback.

**Reproduce (informal):** grepped the entire codebase for
`create_qmmm_subcommand(gaussian)` (the top-level group) — zero
matches; `qmmm` is only ever attached via `create_qmmm_subcommand` to
`opt`, `ts`, `sp`, `scan`, `qrc`, and `modred`, confirming
`ctx.invoked_subcommand` at the `gaussian` group level can never equal
`"qmmm"`.

**Impact:** None — dead code only; `chemsmart/cli/gaussian/qmmm.py`
already fully owns QMMM option definitions and any necessary molecule
conversion for its own subcommand scope.

**Suggested direction:** remove `click_gaussian_qmmm_options` (lines
361-490) entirely, and remove the `ctx.invoked_subcommand == "qmmm"`
try/except block from the `gaussian()` group callback, since it can
never execute.

## 27. `chemsmart/cli/gaussian/qmmm.py`'s `_populate_charge_and_multiplicity_on_settings` is an unused orphaned duplicate

**Location:** `chemsmart/cli/gaussian/qmmm.py`, lines 385-413.

This module defines
`_populate_charge_and_multiplicity_on_settings(qs)` (charge/multiplicity
fallback resolution from `charge_intermediate`/`charge_high`/
`charge_total`), but nothing in `qmmm()` (the actual QMMM subcommand
callback defined just above it in the same file) ever calls it, and no
other module imports it from here. `chemsmart/cli/orca/qmmm.py` has its
own **separately defined** copy of a function with the exact same name
(line 460 in that file) which *is* actually called (line 418) — this
Gaussian-side copy is simply an unused leftover, presumably from
factoring the ORCA version out of a shared original or vice versa,
never wired up on the Gaussian side.

**Reproduce (informal):** grepped the whole codebase for
`_populate_charge_and_multiplicity_on_settings` — the only call site is
inside `chemsmart/cli/orca/qmmm.py`, calling *that file's own*
same-named function, never this one.

**Impact:** None today — the Gaussian `qmmm()` callback resolves
charge/multiplicity through its own inline CLI-option assignments
(`charge_total`/`mult_total` etc. from lines ~294-305), so nothing is
missing in practice; this is simply 29 lines of dead code.

**Suggested direction:** delete
`_populate_charge_and_multiplicity_on_settings` from
`chemsmart/cli/gaussian/qmmm.py`, or — if the intent was for the
Gaussian `qmmm()` callback to also use it instead of duplicating the
same fallback logic inline — wire it in and remove the redundant
inline assignments.

## 28. ORCA `qmmm` subcommand's `-h`/`--high-level-h-bond-length` option is unusable for any real value

**Location:** `chemsmart/cli/orca/qmmm.py`, the `-h` option
declaration (`type=dict`, lines ~152-157) and its later use
(`ast.literal_eval(high_level_h_bond_length)` at lines 439-442).

```python
@click.option(
    "-h",
    "--high-level-h-bond-length",
    type=dict,
    help="Custom high-level-H bond lengths",
)
...
if high_level_h_bond_length is not None:
    high_level_h_bond_length_dict = ast.literal_eval(
        high_level_h_bond_length
    )
    molecule.scale_factors = high_level_h_bond_length_dict
```

`type=dict` in Click just wraps the builtin `dict` callable as the
option's converter — it calls `dict(<raw CLI string>)` on whatever the
user typed. Python's `dict()` constructor only accepts a mapping or an
iterable of key-value pairs; calling it on an arbitrary string (e.g.
`"{1: 1.1}"`, the exact kind of value the later
`ast.literal_eval(...)` call is clearly meant to parse) always raises,
so Click rejects the option with `Invalid value for '-h': ...` before
the callback body ever runs. The only string that doesn't immediately
error is `""` (`dict("")` → `{}`), but that's an empty dict — not a
string — so the subsequent `ast.literal_eval({})` would itself raise
`TypeError: literal_eval() ... expected string`.

**Reproduce:**
```
chemsmart run orca -p <project> -f mol.xyz opt qmmm -h "{1: 1.1}"
# Error: Invalid value for '-h' / '--high-level-h-bond-length': {1: 1.1}
```
(also reproduced via
`tests/test_orca_qmmm_cli.py::TestOrcaQmmmSubcommand::test_high_level_h_bond_length_option_is_unusable`)

**Impact:** The `-h/--high-level-h-bond-length` CLI option can never
be used to actually set custom high-level-H bond lengths — every
attempt to pass a real value fails at argument parsing, making this
documented feature completely inaccessible from the CLI (the
underlying `ORCAQMMMJobSettings.high_level_h_bond_length` attribute and
`molecule.scale_factors` plumbing are otherwise intact and would work
if the value ever reached them).

**Suggested direction:** change the option to `type=str` (matching how
`-sf`/`--scale-factors` and `-ba`/`--bonded-atoms` are declared
elsewhere in this same file, both of which are later parsed with
`ast.literal_eval`/similar) so the raw string reaches the existing
`ast.literal_eval` call intact.

## 29. ORCA `ts` subcommand: `-j scants` (job type) is silently overridden back to "optts"

**Location:** `chemsmart/cli/orca/ts.py`, lines 184-193.

```python
jobtype_normalized = (jobtype or "").lower()
cli_tssearch_type = tssearch_type.lower() if tssearch_type else None
effective_tssearch_type = ts_settings.tssearch_type or "optts"

if jobtype_normalized == "scants":
    effective_tssearch_type = "scants"
if cli_tssearch_type is not None:
    effective_tssearch_type = cli_tssearch_type

ts_settings.tssearch_type = effective_tssearch_type
```

The `-ts`/`--tssearch-type` click option is declared with
`default="optts"` (never `None` unless something very unusual
happens), so `cli_tssearch_type` is **never** `None` for any real CLI
invocation. That means the second `if` (line 190) always fires and
always overwrites whatever the first `if` (line 188-189, the
`-j scants` jobtype-inference branch) just set — `-j scants` alone
(without also explicitly passing `-ts scants`) has no effect at all;
`tssearch_type` silently stays `"optts"`.

**Reproduce:**
```
chemsmart run orca -p <project> -f mol.xyz ts -j scants
# tssearch_type ends up "optts", not "scants" -- the "-j scants" is ignored.
```
(also reproduced via
`tests/test_orca_cli.py::TestORCACLITsSubcommand::test_jobtype_scants_alone_is_overridden_by_tssearch_type_default`)

**Impact:** The documented convenience of triggering ScanTS mode via
`-j scants` (mentioned in the `ScanTS (--tssearch-type scants or
-j scants) requires ...` error message elsewhere in this same file)
doesn't actually work — users must always use the explicit
`-ts/--tssearch-type scants` flag; `-j scants` on its own silently
does nothing, which will confuse anyone following the error message's
own suggestion.

**Suggested direction:** give `-ts/--tssearch-type` a `default=None`
(matching the "only override if explicitly given" pattern used for
every other option in this same function) and apply the "optts"
fallback separately, so `-j scants` isn't unconditionally overridden
when the user hasn't explicitly set `-ts`.

## 30. `orca` CLI group: default job label always doubles the subcommand-name suffix (and crashes for filename-less jobs)

**Location:** `chemsmart/cli/orca/orca.py`, the `orca()` group callback's
label-resolution block (lines ~753-768).

```python
if label is None and append_label is None:
    label = os.path.splitext(os.path.basename(filename))[0]
    if filename:
        label = os.path.splitext(os.path.basename(filename))[0]
    else:
        label = "output"
    if ctx.invoked_subcommand:
        label = f"{label}_{ctx.invoked_subcommand}"
    if is_chemsmart_db:
        if structure_id is not None:
            label = f"{label}_SID-{structure_id}"
        elif record_id is not None:
            label = f"{label}_RID-{record_id}"
        elif record_index is not None:
            label = f"{label}_RI-{record_index}"
    label = f"{label}_{ctx.invoked_subcommand}"
```

Two problems, mirroring bug #25 in the Gaussian CLI's equivalent block
but slightly worse here:

1. **Crash when `filename` is `None`.** Exactly like #25: the first
   line unconditionally calls `os.path.basename(filename)` before the
   `if filename: ... else: label = "output"` guard a few lines down,
   so a `--pubchem`-only job with no `-f`/`-l`/`-a` crashes with
   `TypeError` instead of falling back to `"output"`.
2. **The subcommand-name suffix is *always* appended twice.** Line
   `if ctx.invoked_subcommand: label = f"{label}_{ctx.invoked_subcommand}"`
   conditionally appends the subcommand name — but the very last line
   of the block unconditionally does the exact same append again,
   regardless of whether the conditional branch above it ran. For any
   normal invocation (a subcommand is always given), this means every
   default-labeled ORCA job's label ends up as
   `<basename>_<subcommand>_<subcommand>` instead of
   `<basename>_<subcommand>`.

**Reproduce:**
```python
from unittest.mock import MagicMock, patch
from click.testing import CliRunner
from chemsmart.cli.orca.orca import orca

with patch("chemsmart.jobs.orca.opt.ORCAOptJob") as mock_job_cls:
    mock_job_cls.return_value = MagicMock()
    CliRunner().invoke(
        orca,
        ["-p", "gas_solv", "-f", "mol.xyz", "-c", "0", "-m", "1", "opt"],
    )
    print(mock_job_cls.call_args[1]["label"])
    # "mol_opt_opt" -- "_opt" appears twice
```
(also reproduced via
`tests/test_orca_cli.py::TestORCACLIGroupValidation::test_default_label_doubles_subcommand_suffix`
and `...::test_pubchem_only_without_label_crashes`)

**Impact:** Every ORCA job submitted without an explicit `-l`/`-a`
label gets a doubled subcommand suffix in its output filenames/labels
(e.g. `mol_opt_opt.inp` instead of `mol_opt.inp`), which is cosmetically
wrong and could cause confusion or collide unexpectedly with
differently-named files. Separately, any PubChem-only ORCA job
submission with no explicit label crashes outright instead of using
`"output"` as documented.

**Suggested direction:** remove the unconditional final
`label = f"{label}_{ctx.invoked_subcommand}"` line (the conditional
`if ctx.invoked_subcommand:` block above it already handles this
correctly), and guard the whole block with `if filename:` before
computing the basename, mirroring the fix suggested for bug #25.

## 31. `chemsmart/cli/gaussian/link.py`'s `jobtype is None` label branch is unreachable dead code

**Location:** `chemsmart/cli/gaussian/link.py`, lines 136-139.

```python
if jobtype is None:
    label = label
else:
    label += f"_{jobtype}"
    ...
```

`link()` calls `get_setting_from_jobtype_for_gaussian(project_settings,
jobtype, ...)` near the top of the function (line 79), and that helper
(`chemsmart/utils/cli.py:393-394`) does:
```python
if jobtype is None:
    raise ValueError("Jobtype must be provided for Crest and Link job.")
```
So by the time execution reaches line 136, `jobtype` can never be
`None` — the function would already have raised `ValueError` before
this point. The `if jobtype is None: label = label` branch (a no-op
assignment even if it were reachable) can never execute.

**Reproduce (informal):** `tests/test_gaussian_cli.py::TestGaussianCLILinkCommand::test_link_requires_jobtype`
already demonstrates that `link` with no `-j` raises `ValueError`
before any label logic runs, confirming `jobtype` is always non-`None`
by line 136.

**Impact:** None — purely dead code left over, presumably, from before
the `get_setting_from_jobtype_for_gaussian` validation was added (or
tightened) to make `jobtype` mandatory.

**Suggested direction:** remove the `if jobtype is None: label = label`
branch and keep only the `else` body's logic unconditionally, since
`jobtype` is guaranteed non-`None` at this point.

## 32. `AtomsChargeMultiplicity.from_atoms`'s bare-`FixAtoms` branch is unreachable through any public ASE API

**Location:** `chemsmart/io/molecules/atoms.py`, `from_atoms`
(lines 131-148).

```python
if atoms.constraints:
    if isinstance(atoms.constraints, list):
        for i, constraint in enumerate(atoms.constraints):
            if isinstance(constraint, FixAtoms):
                ...
    elif isinstance(atoms.constraints, FixAtoms):
        indices = FixAtoms.todict(atoms.constraints)["kwargs"]["indices"]
        ...
```

The `elif isinstance(atoms.constraints, FixAtoms):` branch assumes
`atoms.constraints` can sometimes be a bare (non-list) constraint
object. But ASE's own `Atoms.constraints` property setter (which is
`set_constraint`) always normalizes to a list:

```python
# ase/atoms.py, Atoms.set_constraint
if constraint is None:
    self._constraints = []
elif isinstance(constraint, list):
    self._constraints = constraint
elif isinstance(constraint, tuple):
    self._constraints = list(constraint)
else:
    self._constraints = [constraint]   # <-- always wraps a single constraint
```

So `atoms.constraints = FixAtoms(...)` (or `atoms.set_constraint(...)`)
always ends up as `[FixAtoms(...)]`, never a bare `FixAtoms`. There is
no supported way to make `atoms.constraints` a bare `FixAtoms` instance
through ASE's public API — only by writing the private
`atoms._constraints` attribute directly, bypassing the property
entirely.

**Reproduce (informal):**
```python
from ase import Atoms
from ase.constraints import FixAtoms
a = Atoms("Ar2", positions=[(0, 0, 0), (3.5, 0, 0)])
a.constraints = FixAtoms(indices=[0])
print(type(a.constraints))  # <class 'list'>, not FixAtoms
```
`tests/test_atoms_charge_multiplicity_unit.py::TestFromAtomsSingleFixAtomsConstraint::test_single_fixatoms_constraint_not_wrapped_in_list`
only reaches this branch by assigning the private `_constraints`
attribute directly (`simple_ase_atoms._constraints = FixAtoms(...)`),
confirming the public-API path can't.

**Impact:** None — the `isinstance(atoms.constraints, list)` branch
above it already handles every constraint (including a lone `FixAtoms`
wrapped in a single-element list) that any real ASE `Atoms` object can
carry, so no real input is mishandled. This is purely defensive dead
code for an ASE internal representation that doesn't occur in
practice.

**Suggested direction:** remove the `elif isinstance(atoms.constraints,
FixAtoms):` branch, since `atoms.constraints` is always either `[]` or
a `list` for any `Atoms` object constructed or mutated through ASE's
public API.

## 33. `GaussianInputWriter._append_gen_genecp_basis`'s trailing-newline check is always true

**Location:** `chemsmart/jobs/gaussian/writer.py`, `_append_gen_genecp_basis`
(lines 487-491).

```python
f.write(genecp_section.string)
# Check that the last line of genecp_section.string is empty,
# if not, add an empty line
if genecp_section.string_list[-1] != "\n":
    f.write("\n")
```

`genecp_section.string_list` is `genecp_section.string.split("\n")`
(see `chemsmart/io/gaussian/gengenecp.py`, `GenGenECPSection.string_list`).
`str.split("\n")` never returns an element that is literally `"\n"` —
a trailing newline in the source string produces an empty string `""`
as the last element, not `"\n"`, and a string with no trailing newline
produces its last real line (never `"\n"` either, since the line
content itself never contains the delimiter). So
`string_list[-1] != "\n"` is true unconditionally, regardless of
whether `genecp_section.string` actually ends with a blank line, and
`f.write("\n")` always executes.

**Reproduce (informal):**
```python
s1 = "C 0\nsto-3g\n****\n"
print(s1.split("\n")[-1] != "\n")   # True (last element is "")
s2 = "C 0\nsto-3g\n****"
print(s2.split("\n")[-1] != "\n")   # True (last element is "****")
```
Neither case can make the comparison false.

**Impact:** Low — an extra blank line is always appended after the
genecp section, whether or not one was already present. Since Gaussian
input files generally tolerate extra blank lines between sections,
this has not manifested as a job-breaking bug, but the guard does not
do what its comment says.

**Suggested direction:** the intended check was likely
`genecp_section.string_list[-1] != ""` (i.e., whether the string ends
with a newline) rather than comparing against the literal `"\n"`.

## 34. `GaussianInputWriter._write_route_section`'s QMMM branch is unreachable through the normal write pipeline

**Location:** `chemsmart/jobs/gaussian/writer.py`, `_write_route_section`
(lines 176-179), dispatched from `_write_all` (lines 100-108).

```python
# _write_all:
if isinstance(self.settings, GaussianQMMMJobSettings):
    self._write_route_section_qmmm(f)
    ...
else:
    self._write_route_section(f)
    ...

# _write_route_section:
route_string = self.settings.route_string
if isinstance(self.settings, GaussianQMMMJobSettings):
    route_string = self.settings._route_string
```

`_write_all` only calls `_write_route_section` from the `else` branch,
i.e. only when `self.settings` is *not* a `GaussianQMMMJobSettings`
instance (QMMM jobs are routed to the dedicated
`_write_route_section_qmmm` instead). So by the time
`_write_route_section` runs, `isinstance(self.settings,
GaussianQMMMJobSettings)` is guaranteed `False` — the check at line
178 can never be `True` through the real write pipeline.

**Reproduce (informal):**
`tests/test_GaussianWriter.py::TestGaussianInputWriter::test_write_qmmm_job`
(and `test_write_qmmm_input_from_logfile`) already exercise a full
QMMM job write end-to-end and only ever go through
`_write_route_section_qmmm`, never `_write_route_section` — confirming
the two code paths are mutually exclusive by construction.

**Impact:** None — purely dead code, presumably left over from before
`_write_route_section_qmmm` was split out into its own method.

**Suggested direction:** remove the `isinstance(self.settings,
GaussianQMMMJobSettings)` check and the `self.settings._route_string`
fallback from `_write_route_section`, since that method is now only
ever called for non-QMMM settings.

## 35. `update_irc_label` always appends `_flat` even when `direction is None`, contradicting its own docstring and causing a double `_flat` suffix downstream

**Location:** `chemsmart/utils/cli.py`, `update_irc_label` (lines 725-756).

```python
def update_irc_label(label, direction, flat_irc):
    """
    ...
    Appends 'f' for forward direction, 'r' for reverse direction,
    and '_flat' if flat_irc is True and a direction is specified.

    When direction is None (both forward and reverse IRC sub-jobs will be
    created), the '_flat' suffix is not added here because the sub-jobs
    created by ``_ircf_job()`` / ``_ircr_job()`` append the direction and
    '_flat' suffix themselves, avoiding double-application.
    ...
    """
    if direction is not None:
        if direction.lower() == "forward":
            label += "f"
        elif direction.lower() == "reverse":
            label += "r"
        else:
            raise ValueError(...)
    if flat_irc and not label.endswith("_flat"):
        label += "_flat"
    return label
```

The docstring explicitly promises that when `direction is None`, the
`_flat` suffix is deferred to `_ircf_job()`/`_ircr_job()` to avoid
double-application. But the `if flat_irc and not label.endswith(...)`
check at the end is unconditional — it is not nested inside `if
direction is not None:` — so it runs regardless of `direction`,
directly contradicting the docstring.

This was introduced by commit `ba929ad9` ("update irc flat labels
(#649)", merged into this branch from `main`), which unindented what
had previously been a `direction is not None`-guarded block (added by
the earlier commit `ac67b1b2`, "Fix IRC flat job label: don't add
_flat when direction is None") back out to always run, while leaving
the docstring's claim about `direction is None` unchanged.

The real downstream consequence (from `chemsmart/jobs/gaussian/irc.py`,
`_ircf_job`/`_ircr_job`, lines 96-123 and 125-155+): when
`self.settings.direction is None`, these methods derive the sub-job
label from `self.label` (which was already built by
`update_irc_label` at the CLI layer) by checking
`if label.endswith("_irc"): label += "f"/"r" else: label += "_ircf"/"_ircr"`,
then unconditionally append `"_flat"` again if `self.settings.flat_irc`
is true — with no "already ends with _flat" guard the way
`update_irc_label` has. So for a `direction=None`, `flat_irc=True` job:
1. `update_irc_label` now produces e.g. `"mol_flat"` (bug: should be `"mol"`).
2. `_ircf_job` sees `"mol_flat"` doesn't end with `"_irc"`, so appends
   `"_ircf"` → `"mol_flat_ircf"`.
3. Since `settings.flat_irc` is true, appends `"_flat"` again →
   `"mol_flat_ircf_flat"` — a doubled, out-of-order `_flat` suffix.

**Reproduce (informal):**
```python
from chemsmart.utils.cli import update_irc_label
update_irc_label("mol", None, True)  # returns "mol_flat", not "mol"
```
`tests/test_utils_cli.py::TestUpdateIrcLabel::test_none_direction_unchanged`
now fails against this behavior and has been updated to assert the
actual (buggy) output while referencing this entry.

**Impact:** Medium — any IRC job run without an explicit `--direction`
(i.e. both forward and reverse sub-jobs auto-created) combined with
`--flat-irc`/`flat_irc=True` gets a malformed, doubled `_flat` suffix
in both sub-job labels, e.g. `mol_flat_ircf_flat` /
`mol_flat_ircr_flat` instead of the intended `mol_ircf_flat` /
`mol_ircr_flat`.

**Suggested direction:** re-nest the `if flat_irc and not
label.endswith("_flat"): label += "_flat"` check inside the `if
direction is not None:` block in `update_irc_label`, restoring the
behavior described in its own docstring (and originally introduced by
`ac67b1b2`).

## 36. `ThermochemistryJob.compute_thermochemistry`'s docstring claims a `ValueError` that can no longer be raised

**Location:** `chemsmart/jobs/thermochemistry/job.py`,
`compute_thermochemistry` (lines 273-288).

```python
def compute_thermochemistry(self):
    """
    ...
    Raises:
        ValueError: If no input file is provided
        Exception: If calculation fails during processing
    """
    # Set default output file if not specified
    if self.settings.outputfile is None:
        self.settings.outputfile = self.outputfile
    ...
```

The docstring claims `compute_thermochemistry` raises `ValueError` when
no input file is provided, but there is no such check anywhere in the
method body — it goes straight to `self.outputfile`, which requires
`self.label`/`self.folder` to be set. `ThermochemistryJob.__init__` now
unconditionally validates `filename is not None` (see the update to
entry #1 above), so under any legitimate construction path,
`self.filename` can never be `None` by the time
`compute_thermochemistry` runs — the check the docstring describes
would be genuinely unreachable/redundant if it still existed, and
apparently was removed for that reason without updating the docstring.

The only way to observe a missing-filename-like failure here is to
bypass `__init__` entirely (e.g. `ThermochemistryJob.__new__(...)`),
which then fails with an unrelated `AttributeError` (`'NoneType'
object has no attribute 'label'`/similar) from the `self.outputfile`
property, not the documented `ValueError`.

**Reproduce (informal):**
```python
job = ThermochemistryJob.__new__(ThermochemistryJob)
job.filename = None
job.settings = ThermochemistryJobSettings()
job.compute_thermochemistry()
# AttributeError: 'ThermochemistryJob' object has no attribute 'label'
# (not the docstring's claimed ValueError)
```
`tests/test_thermochemistry_job_unit.py::TestThermochemistryJobComputeAndShow::test_compute_thermochemistry_requires_filename`
was rewritten to instead test the real, reachable guard — that
`ThermochemistryJob(filename=None)` raises `ValueError` from
`__init__` — since that is the only place this validation can actually
occur through legitimate construction.

**Impact:** None — purely a stale docstring describing behavior that
either was removed as redundant (now handled at `__init__` time) or
never existed after a refactor.

**Suggested direction:** drop the `Raises: ValueError: If no input
file is provided` line from `compute_thermochemistry`'s docstring,
since that validation now lives entirely in `__init__`.

## 37. `chemsmart/cli/mol/align.py`'s `align()` command has several defensive checks that are unreachable safeguards

**Location:** `chemsmart/cli/mol/align.py`, `align()` (lines 49-238).

Several checks in this function are dead code — each is guaranteed to
never trigger given the validation performed earlier in the same
function (or, in one case, by the CLI's option semantics):

1. **Lines 136-141** (inner `else` inside `if directory:` / `if
   filetype:`):
   ```python
   if directory:
       ...
       if filetype:
           ...
       else:
           # This should not happen due to
           # validation above, but keep as safeguard
           raise click.BadParameter(
               "Directory specified but no filetype provided. ..."
           )
   ```
   The outer guard at line 61 (`if directory and not filetype: raise
   click.BadParameter(...)`) already ensures that whenever `directory`
   is truthy inside this block, `filetype` is also truthy — the inner
   `else` can never execute. The comment ("should not happen ... keep
   as safeguard") already acknowledges this.

2. **Lines 181-185** (final `else` of the `if directory: ... elif
   filenames: ... else:` chain):
   ```python
   else:
       # This should not happen due to validation above, but keep as
       # safeguard
       raise click.BadParameter(
           "No input files specified. This should have been caught
           earlier."
       )
   ```
   The very first check in the function (line 56) already raises
   whenever both `filenames` and `directory` are falsy, so reaching
   this branch would require `directory` falsy and `filenames` falsy
   simultaneously — already excluded.

3. **Lines 68-78**, the final implicit `else` (no branch taken, i.e.
   `index` stays `None`): reaching `if index is None:` with none of
   the three `if`/`elif` conditions true requires `filenames` falsy
   *and* `not (directory and filetype)`. Combined with line 56's
   guard (which requires `directory` truthy whenever `filenames` is
   falsy) and line 61's guard (which requires `filetype` truthy
   whenever `directory` is truthy), this combination is impossible.

4. **Lines 187-189**:
   ```python
   if not isinstance(molecules, list):
       molecules = list(molecules) if molecules else []
   ```
   `molecules` is initialized as `[]` (line 79) and only ever mutated
   via `.extend(...)`, so it is always a `list` by construction —
   this coercion never fires.

5. **Line 208** (`base_label = "molecules"` fallback): by the point
   label generation runs, either the `directory` branch has set
   `base_file_for_label = matched_files[0]` (line 135, only reached
   when `matched_files` is non-empty, since empty is caught at line
   122) or the `filenames` branch has set `base_file_for_label =
   filenames[0]` (line 179, unconditional whenever `filenames` is
   truthy) — the final `else` chain guarantees one of these two
   branches always runs before label generation, so
   `base_file_for_label` is always truthy here.

6. **Line 148** (`if isinstance(filenames, str): filenames =
   [filenames]`): every code path in `chemsmart/cli/mol/mol.py` that
   sets `ctx.obj["filenames"]` for the `align` subcommand assigns
   either the tuple produced by click's `multiple=True` `-f` option
   or a `list` from `glob.glob(...)` — never a bare string. This
   coercion is unreachable through the real CLI.

**Reproduce (informal):** `tests/test_mol_cli.py::TestMolCLIAlignCommand`
now covers every legitimately reachable branch in `align()` (directory
without filetype via the `-p`/program group path, directory with
filetype, single- and multi-file label generation for both the
`<=2` and `>2` structure cases, explicit `-l` label, not-enough-
molecules, no-files-found-for-filetype, and out-of-range index errors
from both the single-file and per-file-helper code paths); coverage
tops out at 91% with only the branches listed above remaining
uncovered, confirming they are unreachable.

**Impact:** None — all are defensive checks whose own comments (in
cases 1 and 2) already flag them as belt-and-suspenders code that
"should not happen."

**Suggested direction:** remove the six dead branches above; if
desired, replace the multi-condition duplication (case 3) with a
single unconditional fallback, and drop the `isinstance(..., str)`
coercion (case 6) since it does not correspond to any real call path.

## 38. `mol_qmmm` CLI group is entirely orphaned — defined but never reachable

**Location:** `chemsmart/cli/mol/mol.py`, `mol_qmmm` (lines 646-805) and
`mol_qmmm_process_pipeline` (lines 798-805).

`mol_qmmm` is a fully-implemented `@click.group(cls=MyGroup)` that
mirrors the regular `mol` group's molecule-loading logic but converts
everything to `QMMMMolecule`. However:

1. **No subcommand is ever registered on it.** Every subcommand file in
   `chemsmart/cli/mol/` (`align.py`, `nci.py`, `mo.py`, `movie.py`,
   `spin.py`, `irc.py`, `visualize.py`) decorates its command with
   `@mol.command(...)` — none use `@mol_qmmm.command(...)`. A `grep`
   across the whole `chemsmart/` tree confirms zero occurrences of
   `@mol_qmmm.command`, `mol_qmmm.add_command`, or
   `add_command(mol_qmmm)`.
2. **It is never attached to any parent CLI group** — nothing imports
   `mol_qmmm` anywhere else in the codebase (only `mol.py` itself
   defines it), so it is not reachable as `chemsmart run mol-qmmm ...`
   or under any other entry point.
3. Even if it *were* invoked directly (e.g. in a test via
   `CliRunner().invoke(mol_qmmm, [...])`), Click refuses to run it:
   since it is a `MultiCommand` with an empty `commands` dict and no
   subcommand token can ever resolve to one, Click's argument parsing
   raises `Error: Missing command.` (exit code 2) before the group's
   own callback body ever executes — confirmed empirically:
   ```python
   from click.testing import CliRunner
   from chemsmart.cli.mol.mol import mol_qmmm
   print(mol_qmmm.commands)  # {}
   CliRunner().invoke(mol_qmmm, ["-f", "some.xyz"]).output
   # "Error: Missing command."
   ```
   The only way to exercise the callback body at all is to bypass
   Click's command routing entirely and invoke the raw function via
   `click.Context(mol_qmmm)` + `ctx.invoke(mol_qmmm.callback, ...)`,
   which is not a real invocation path available to any user.

**Impact:** None currently (dead code), but it represents ~125 lines of
unreachable, unmaintained duplicate logic that will silently drift out
of sync with the `mol` group it was cloned from (e.g., it does not
handle `-i`/`--si` aliasing at all, unlike `mol`'s equivalent check at
lines 407-414, and its append/default label logic does not merge in
the chemsmart-db-specific suffixes that `mol`'s does).

**Reproduce:** see the code excerpt above; also
`tests/test_mol_cli.py::TestMolQmmmGroupDirectInvocation` exercises the
callback via the `ctx.invoke(mol_qmmm.callback, ...)` bypass technique
to get direct-unit coverage on its logic despite it being unreachable
through the real CLI.

**Suggested direction:** either finish wiring `mol_qmmm` up (attach
QMMM-aware subcommands and register the group somewhere reachable), or
remove it entirely if QMMM-aware molecule loading for `mol` subcommands
is meant to happen some other way (e.g. via the `--qmmm`
flag/`ctx.obj["qmmm"]` mechanism already used elsewhere in the
codebase).

## 39. `boltzmann` CLI command's `outputfile` parameter is vestigial — always `None`

**Location:** `chemsmart/cli/thermochemistry/boltzmann.py`, `boltzmann()`
(function signature and body).

```python
def boltzmann(
    ctx,
    skip_completed,
    energy_type_for_weighting="gibbs",
    outputfile=None,
):
    ...
    boltzmann_thermochemistry = BoltzmannAverageThermochemistryJob(
        files=files,
        energy_type=energy_type_for_weighting,
        filename=files[0] if files else None,
        outputfile=outputfile,
        ...
    )
```

`boltzmann()` has an `outputfile=None` parameter, but there is no
`@click.option` anywhere in this file (or inherited via
`click_job_options`) that supplies a value for it — `-o`/`--outputfile`
is only defined on the parent `thermochemistry` group itself (and
consumed there to build each per-file job's `job_settings.outputfile`
and to set `ctx.obj["outputfile"]`). Consequently `outputfile` in
`boltzmann()` is **always** `None` regardless of what the user passes;
attempting `chemsmart run thermochemistry -f a.log -o out.dat
boltzmann` even fails outright, since `-o` is not a recognized option
at the `boltzmann` subcommand scope:

```
$ chemsmart run thermochemistry -f a.log -T 298.15 boltzmann -o out.dat
Error: No such option: -o
```

There is currently no way to specify a custom output path for the
Boltzmann-averaged result specifically — the job falls back entirely
to its own `outputfile` property (`{label}.dat`).

**Reproduce:** `tests/test_thermochemistry_cli.py::TestThermochemistryBoltzmannCommand::test_outputfile_parameter_is_always_none`
confirms `outputfile` in the constructed job's kwargs is `None` even
when the group-level `-o` is supplied before the `boltzmann` subcommand
token (since it's consumed by the group, not forwarded).

**Impact:** Low — users cannot redirect the Boltzmann-averaged output
file via CLI; the default `{label}.dat` naming is always used instead.

**Suggested direction:** either read `ctx.obj.get("outputfile")` inside
`boltzmann()` and pass that through, or drop the dead `outputfile`
parameter from the function signature entirely if this was never
intended to be configurable per-command.

## 40. `add_lines_in_yaml_files`'s `except yaml.YAMLError` clause is unreachable

**Location:** `chemsmart/cli/config.py`, `add_lines_in_yaml_files`
(lines 411-456).

```python
for yaml_file in target_directory.glob("*.yaml"):
    try:
        with open(yaml_file, "r") as f:
            yaml_content = f.readlines()  # Read file line by line
        ...
        with open(yaml_file, "w") as f:
            for line in updated_content:
                f.write(line)
    except yaml.YAMLError as e:
        logger.info(f"Error reading {yaml_file}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error while processing {yaml_file}: {e}")
```

Despite the `.yaml` extension and the `import yaml` at the top of the
module, this function never actually parses YAML — it reads the file
as a plain list of lines (`f.readlines()`) and writes lines back
verbatim (`f.write(line)`), treating the file as arbitrary text so it
can insert lines after a matching marker. No `yaml.safe_load`,
`yaml.load`, or any other YAML-parsing call exists anywhere in this
function's body, so `yaml.YAMLError` can never actually be raised here
— the specific `except yaml.YAMLError` branch is dead code, and any
real failure (a missing file disappearing mid-scan, a permissions
error, etc.) falls through to the generic `except Exception` clause
instead.

**Reproduce:** `tests/test_config.py::TestAddLinesInYamlFilesFunction::test_malformed_yaml_is_logged_not_raised`
demonstrates that a simulated I/O failure (`OSError` from a patched
`open`) is caught by the generic `except Exception` clause; there is
no way to make the `except yaml.YAMLError` clause fire since nothing
in the `try` block can raise that exception type.

**Impact:** None — purely dead exception-handling code, presumably
left over from an earlier version of this function that used a real
YAML parser, or added defensively by analogy with the `.yaml` file
extension without checking whether it was actually reachable.

**Suggested direction:** remove the `except yaml.YAMLError` branch
(and the now-unused `import yaml` if nothing else in the module needs
it) since the function never parses YAML — it only rewrites lines of
text within `.yaml`-named files.

## 41. `get_thermochemistry.py`'s `except TypeError` handler crashes with `UnboundLocalError` when the constructor's frequency-missing warning is actually needed

**Location:** `chemsmart/scripts/get_thermochemistry.py`, `entry_point`
(lines 365-493).

```python
for file in filenames:
    try:
        thermochemistry = Thermochemistry(file, ...)
        structure = os.path.splitext(os.path.basename(file))[0]
        energy = thermochemistry.electronic_energy * unit_conversion
        zero_point_energy = thermochemistry.zero_point_energy * unit_conversion
        ...
    except TypeError:
        log(
            "{:2} {:39} {:13.6f} {:<50}\n".format(
                " ×",
                structure,
                energy,
                "  Warning! Frequency information not found ...",
            )
        )
        continue
```

The `except TypeError` handler's own message ("Warning! Frequency
information not found") makes clear its intended purpose: when a
structure lacks frequency data, some `thermochemistry.*` property
computation returns `None`, and multiplying `None * unit_conversion`
raises `TypeError`. But by the time that first arithmetic line
(`energy = thermochemistry.electronic_energy * unit_conversion`)
raises, `energy` itself has **not yet been assigned** — the handler
then tries to format `energy` into the warning message and crashes
with `UnboundLocalError: local variable 'energy' referenced before
assignment`, masking the original TypeError entirely and aborting the
whole script (since `UnboundLocalError` is not itself caught by
anything here).

**Reproduce:**
```python
from unittest.mock import MagicMock, patch
from click.testing import CliRunner
from chemsmart.scripts.get_thermochemistry import entry_point

runner = CliRunner()
with runner.isolated_filesystem():
    open("mol.log", "w").write("dummy")
    with patch("chemsmart.scripts.get_thermochemistry.Thermochemistry") as mock_cls:
        thermo = MagicMock()
        thermo.electronic_energy = None  # simulates missing frequency data
        mock_cls.return_value = thermo
        result = runner.invoke(entry_point, ["-f", "mol.log"], catch_exceptions=True)
        print(result.exception)
        # UnboundLocalError: local variable 'energy' referenced before assignment
```
See `tests/test_scripts_analysis_cli.py::TestGetThermochemistryScript::test_missing_frequency_data_crashes_instead_of_warning`.

**Impact:** Medium — any real output file missing frequency
information (the exact case this handler exists to report gracefully)
instead crashes the entire batch run with a confusing
`UnboundLocalError`, rather than logging a per-structure warning and
continuing to the next file as intended.

**Suggested direction:** initialize `energy = None` (and any other
values referenced in the `except` block) before the `try`, or restrict
the message to `structure` alone / use a placeholder string when
`energy` was never computed.

## 42. `pka batch`'s proton-exchange reference validation catches the wrong exception type, so a failed CDXML reference-proton auto-detect crashes instead of producing a helpful `UsageError`

**Location:** `chemsmart/cli/gaussian/pka.py`, `batch()` (lines
266-278).

```python
elif shared["reference_proton_index"] is None:
    ref = shared["reference"]
    if ref.endswith((".cdx", ".cdxml")):
        try:
            shared["reference_proton_index"] = (
                PKaCDXFile.resolve_reference_proton(
                    ref,
                    None,
                    shared["reference_color_code"],
                )
            )
        except click.UsageError:
            missing.append("-rpi/--reference-proton-index")
```

`PKaCDXFile.resolve_reference_proton` (`chemsmart/io/file.py`, line
582) explicitly documents and raises `ValueError` on parse/color-code
failure, never `click.UsageError`. The `except click.UsageError:`
clause here can therefore never catch a real failure from this call —
when reference-proton auto-detection genuinely fails (e.g. the
reference CDXML has no uniquely coloured proton), the `ValueError`
propagates uncaught out of `batch()` instead of being converted into
the intended "missing -rpi/--reference-proton-index" `UsageError`
alongside the other reference-option checks.

**Reproduce:**
```python
from unittest.mock import patch
import click

with patch(
    "chemsmart.cli.gaussian.pka.PKaCDXFile.resolve_reference_proton",
    side_effect=ValueError("no coloured proton found"),
):
    ...  # invoking `pka batch` with scheme="proton exchange" and a
         # .cdxml reference lacking reference_proton_index raises
         # ValueError instead of click.UsageError
```
See `tests/test_pka.py::TestPKa::test_batch_proton_exchange_cdxml_reference_resolve_failure_raises_usage_error`,
which patches the exception type raised to demonstrate both what the
code currently does (crash) and what the `except` clause was clearly
meant to handle.

**Impact:** Low-medium — only affects the specific combination of
`pka batch` + `scheme="proton exchange"` + a `.cdxml`/`.cdx` reference
file with no `--reference-proton-index` given + auto-detection
failing. Users hit a raw `ValueError` traceback instead of the
friendly, actionable `UsageError` listing all missing reference
options.

**Suggested direction:** change the `except` clause to catch
`ValueError` (matching what `resolve_reference_proton` actually
raises), or have `resolve_reference_proton` raise `click.UsageError`
directly if that is the intended contract across its other callers.

## 43. `Gaussian16Input.constrained_atoms` setter always crashes with `RecursionError` instead of setting anything

**Location:** `chemsmart/io/gaussian/input.py`, lines 270-281.

```python
@property
def constrained_atoms(self):
    """
    Get atoms with coordinate constraints.
    """
    return self.coordinate_block.constrained_atoms

@constrained_atoms.setter
def constrained_atoms(self, value):
    """
    Set atoms with coordinate constraints.
    """
    self.constrained_atoms = value
```

The setter assigns to `self.constrained_atoms`, which is the same
property it is defining -- this re-invokes the setter itself,
infinitely, until Python's recursion limit is hit. Any attempt to do
`some_gaussian_input.constrained_atoms = [...]` crashes with
`RecursionError: maximum recursion depth exceeded` instead of storing
the value anywhere (there is no backing `_constrained_atoms`
attribute for it to write to).

**Reproduce:**
```python
from chemsmart.io.gaussian.input import Gaussian16Input

gi = Gaussian16Input(filename="some_valid_input.com")
gi.constrained_atoms = [1, 2]
# RecursionError: maximum recursion depth exceeded
```
See `tests/test_GaussianIO.py::TestGaussian16InputDirectPropertyCoverage::test_constrained_atoms_setter_recurses_infinitely`.

**Impact:** Low today -- grep shows no production code anywhere
assigns to `.constrained_atoms` on a `Gaussian16Input`/
`Gaussian16QMMMInput` instance, so this is currently dead code. But
the setter exists and is part of the public property surface; anyone
who reaches for it (a very natural thing to try, given the getter)
will hit an opaque crash instead of the constraint being recorded.

**Suggested direction:** either remove the setter entirely (the
value is derived from `coordinate_block`, so a plain read-only
property may be all that's intended), or have it write to a real
backing attribute (e.g. `self._constrained_atoms = value`) that the
getter falls back to.

## 44. `Gaussian16Input.oniom_charge`/`oniom_multiplicity` crash with `AttributeError` on any non-QMMM input, because the `RecursionError` guard around `self.partition` catches the wrong exception type

**Location:** `chemsmart/io/gaussian/input.py`, lines 140-147 and
337-363.

```python
@property
def oniom_charge(self):
    oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
    return oniom_charge
...
def _get_oniom_charge_and_multiplicity(self, use_partition=True):
    ...
    if use_partition:
        try:
            partition_len = len(self.partition)
        except RecursionError:
            partition_len = None
```

`oniom_charge`/`oniom_multiplicity` call
`_get_oniom_charge_and_multiplicity()` with the default
`use_partition=True`, which accesses `self.partition`. But
`partition` is only defined on the `Gaussian16QMMMInput` subclass
(line 488) -- the base `Gaussian16Input` class has no such attribute
at all. Accessing it raises `AttributeError`, not `RecursionError`,
so the `except RecursionError` guard never catches it, and the
`AttributeError` propagates uncaught out of `oniom_charge`/
`oniom_multiplicity` for any plain (non-QMMM) `Gaussian16Input`
instance.

Contrast this with `charge`/`multiplicity` (lines 105-137), which
call the same helper with `use_partition=False` explicitly and so
never hit this path -- only the `oniom_charge`/`oniom_multiplicity`
properties are affected.

**Reproduce:**
```python
from chemsmart.io.gaussian.input import Gaussian16Input

gi = Gaussian16Input(filename="some_valid_non_qmmm_input.com")
gi.oniom_charge
# AttributeError: 'Gaussian16Input' object has no attribute 'partition'
```
See `tests/test_GaussianIO.py::TestGaussian16InputDirectPropertyCoverage::test_oniom_charge_crashes_on_base_class_non_qmmm_input`.

**Impact:** Low-medium -- only affects code that calls
`.oniom_charge`/`.oniom_multiplicity` on a plain `Gaussian16Input`
(rather than `Gaussian16QMMMInput`) instance. Since ONIOM-specific
data only makes sense for QMMM inputs, this may never happen in
practice, but nothing prevents a caller from doing so, and the result
is a confusing `AttributeError` rather than a clean `None`/empty dict.

**Suggested direction:** either catch `AttributeError` as well (or
instead), or guard with `getattr(self, "partition", None)` so the
base class degenerates gracefully instead of crashing.

## 45. `Gaussian16QMMMInput.model_charge`/`model_multiplicity` crash with `KeyError` for real 2-layer ONIOM systems

**Location:** `chemsmart/io/gaussian/input.py`, lines 463-466,
483-486, 505-528.

```python
@property
def model_charge(self):
    oniom_charge, _ = self._get_oniom_charge_and_multiplicity()
    return int(oniom_charge["model_charge"])
...
def _get_oniom_charge_and_multiplicity(self, use_partition=True):
    ...
    if len(self.partition) == 2:
        print(charge_multiplicity_list)
        charge_multiplicity_list = charge_multiplicity_list[0:3]
        full_line = 6
```

For a 2-layer ONIOM system, `Gaussian16QMMMInput`'s own override of
`_get_oniom_charge_and_multiplicity` slices
`charge_multiplicity_list[0:3]`, keeping only
`["charge_total", "real_multiplicity", "int_charge"]` -- it never
includes `"model_charge"`/`"model_multiplicity"` in the resulting
`oniom_charge`/`oniom_multiplicity` dicts for this case. But the
`model_charge`/`model_multiplicity` properties unconditionally index
`oniom_charge["model_charge"]`/`oniom_multiplicity["model_multiplicity"]`,
so calling either on a 2-layer system raises `KeyError`.

Compare with the base class's own `_get_oniom_charge_and_multiplicity`
(lines 337-381), which for the equivalent 2-layer case keeps
`charge_multiplicity_list[0:2] + charge_multiplicity_list[4:6]` --
i.e. `charge_total`/`real_multiplicity`/`model_charge`/
`model_multiplicity` (dropping the intermediate layer, keeping the
model layer) -- the logically consistent choice for a 2-layer system,
where "model" is the inner/high-level layer. The QMMM subclass's
override uses different, inconsistent slicing that drops
`model_charge`/`model_multiplicity` instead of `int_charge`/
`int_multiplicity`.

There's also a stray `print(charge_multiplicity_list)` debug
statement left in this branch (line 526).

**Reproduce:**
```python
from chemsmart.io.gaussian.input import Gaussian16QMMMInput

gi = Gaussian16QMMMInput(filename="a_2layer_oniom_input.com")
gi.model_charge
# KeyError: 'model_charge'
```
See `tests/test_GaussianIO.py::TestGaussian16InputDirectPropertyCoverage::test_qmmm_2layer_model_charge_crashes_with_keyerror`.

**Impact:** Medium -- any 2-layer ONIOM/QMMM job that calls
`.model_charge`/`.model_multiplicity` (e.g. to report or reuse the
inner-layer charge/multiplicity) crashes outright. 3-layer systems are
unaffected (they keep the full 6-entry list).

**Suggested direction:** align the QMMM subclass's 2-layer slicing
with the base class's `[0:2] + [4:6]` (or whatever the intended
semantics are for a 2-layer system), and remove the leftover debug
`print`.

## 46. `Gaussian16Input.charge`/`.multiplicity` return a `str` instead of `int` when falling back to oniom-line parsing

**Location:** `chemsmart/io/gaussian/input.py`, lines 105-137.

```python
@property
def charge(self):
    charge_multiplicity = self._get_charge_and_multiplicity()
    if charge_multiplicity is not None:
        charge, _ = charge_multiplicity
        return charge  # int, from _get_charge_and_multiplicity

    oniom_charge, _ = self._get_oniom_charge_and_multiplicity(
        use_partition=False
    )
    if oniom_charge:
        return oniom_charge.get("charge_total")  # str, unconverted
    return None
```

When the normal single "charge multiplicity" line is found,
`_get_charge_and_multiplicity()` explicitly does `int(line_elements[0])`
/ `int(line_elements[1])` (lines 333-334), so `.charge`/`.multiplicity`
return `int`. But when that line isn't found (e.g. a combined
multi-layer ONIOM charge/mult line with more than two numbers) and the
code falls back to `_get_oniom_charge_and_multiplicity`, the values in
`oniom_charge`/`oniom_multiplicity` are taken directly from
`line.split()` with no `int()` conversion (see `_get_oniom_charge_and_multiplicity`,
lines 322-381) -- so `.charge`/`.multiplicity` silently return a `str`
in this case. Any caller that assumes an `int` (e.g. arithmetic like
`pka_settings.charge - 1` in `chemsmart/cli/gaussian/pka.py`) would
crash with `TypeError: unsupported operand type(s) for -: 'str' and
'int'` if it ever received a Gaussian16Input parsed from such a file.

**Reproduce:**
```python
from chemsmart.io.gaussian.input import Gaussian16Input

gi = Gaussian16Input(filename="a_3layer_oniom_input_with_combined_charge_mult_line.com")
gi.charge  # '0' (str), not 0 (int)
gi.charge - 1  # TypeError: unsupported operand type(s) for -: 'str' and 'int'
```
See `tests/test_GaussianIO.py::TestGaussian16InputDirectPropertyCoverage::test_charge_and_multiplicity_fall_back_to_oniom_parsing`.

**Impact:** Low-medium -- only affects `Gaussian16Input` (not the
`Gaussian16QMMMInput` subclass) parsing a file whose charge/multiplicity
line doesn't match the plain two-integer pattern. Silent type
inconsistency is worse than an upfront crash since it can propagate
before failing far from the actual cause.

**Suggested direction:** wrap the oniom fallback values in `int(...)`
before returning, matching the normal-path behavior.

## 47. `ORCAQMMMInput.qm_force_field` always crashes with `KeyError` -- `_get_qmmm_block`'s `ORCAFFFilename` check is compared against an already-lowercased line

**Location:** `chemsmart/io/orca/input.py`, `_get_qmmm_block` (lines
490-511), used by the `qm_force_field` property (line 287-289).

```python
def _get_qmmm_block(self):
    block = {}
    for line in self.contents:
        line = line.lower()
        if "qmatoms" in line:
            ...
        if "qm2atoms" in line:
            ...
        if "optregion_fixedatoms" in line:
            ...
        if "ORCAFFFilename" in line:
            force_field = line.split()[-1]
            block["force field"] = force_field
    return block
```

`line` is reassigned to `line.lower()` at the top of the loop, but the
`ORCAFFFilename` check still uses the original mixed-case spelling.
Since a lowercased string can never contain the substring
`"ORCAFFFilename"` (only `"orcafffilename"`), this branch can never
match, so `block["force field"]` is never set -- for any input file,
regardless of whether an `ORCAFFFilename` directive is actually
present. `qm_force_field` therefore always raises `KeyError:
'force field'`.

This was already noticed once: `tests/test_ORCAIO.py`'s
`test_orca_qmmm_input` has a commented-out
`# assert orca_inp.qm_force_field` line, suggesting a previous author
hit this and quietly worked around it rather than filing it.

**Reproduce:**
```python
from chemsmart.io.orca.input import ORCAQMMMInput

qi = ORCAQMMMInput(filename="any_qmmm_input_with_ORCAFFFilename.inp")
qi.qm_force_field
# KeyError: 'force field'
```
See `tests/test_ORCAIO.py::TestORCAQMMMInputDirectPropertyCoverage::test_qm_force_field_always_crashes_with_keyerror`.

**Impact:** Medium -- any QM/MM job that specifies a force-field file
via `ORCAFFFilename` and later calls `.qm_force_field` (e.g. to report
or re-emit the setting) crashes outright; the property can never
succeed as written.

**Suggested direction:** compare against the lowercase spelling
(`"orcafffilename"`), matching the other checks in this same loop.

## 48. `PyMOLJob._backup_files(backup_chk=True)` always crashes with `AttributeError`, since no PyMOL job defines a `chkfile` property

**Location:** `chemsmart/jobs/mol/job.py`, `_backup_files` (lines
199-214).

```python
def _backup_files(self, backup_chk=False, **kwargs):
    folder = self._create_backup_folder_name()
    self.backup_file(self.inputfile, folder=folder, **kwargs)
    self.backup_file(self.outputfile, folder=folder, **kwargs)
    if backup_chk:
        self.backup_file(self.chkfile, folder=folder, **kwargs)
```

`PyMOLJob` and its subclasses never define a `chkfile` property or
attribute anywhere (checkpoint files are a Gaussian-specific concept;
PyMOL jobs only have `.xyz`/`.pse`/`.err`/log files). Passing
`backup_chk=True` to `_backup_files` (or to `Job.backup(**kwargs)`,
which forwards to it) unconditionally accesses `self.chkfile`, which
raises `AttributeError` before `backup_file` is ever called.

**Reproduce:**
```python
from chemsmart.jobs.mol.job import PyMOLJob

job = PyMOLJob(molecule=some_molecule, label="test")
job._backup_files(backup_chk=True)
# AttributeError: 'PyMOLJob' object has no attribute 'chkfile'
```
See `tests/test_pymol_job_base_unit.py::TestPyMOLJobBackupFiles::test_backup_chk_crashes_since_pymol_jobs_have_no_chkfile`.

**Impact:** Low -- `backup_chk` is only meaningful for Gaussian jobs
(which do define `chkfile`), so this would only surface if a caller
passed `backup_chk=True` generically to a PyMOL job's `backup()`
call, e.g. via a shared CLI flag applied uniformly across job types.

**Suggested direction:** either drop the `backup_chk` parameter from
`PyMOLJob._backup_files` entirely (it doesn't apply to this job
type), or guard the block with `getattr(self, "chkfile", None)`.

## 49. `get_prepend_string_list_from_modred_free_format`'s single-list branch crashes with `UnboundLocalError` for an invalid `program`, instead of the clean `ValueError` the list-of-lists branch gives

**Location:** `chemsmart/utils/utils.py`, lines 1297-1308.

```python
elif isinstance(input_modred[0], int):
    # for a single list; e.g.: [2,3]
    prepend_string = get_prepend_string_for_modred(input_modred)
    if program == "gaussian" or program == "pymol":
        modred_string = convert_modred_list_to_string(input_modred)
    elif program == "orca":
        modred_string = convert_modred_list_to_string(
            [a - 1 for a in input_modred]
        )
    each_frozen_string = f"{prepend_string} {modred_string}"
    prepend_string_list.append(each_frozen_string)
```

Compare with the list-of-lists branch immediately above it (lines
1280-1296), which has a final `else: raise ValueError(...)` for an
unrecognized `program`. The single-list branch has no such `else` --
if `program` is neither `"gaussian"`/`"pymol"` nor `"orca"`,
`modred_string` is never assigned, and the very next line references
it, raising `UnboundLocalError: local variable 'modred_string'
referenced before assignment` instead of a clear, actionable
`ValueError`.

**Reproduce:**
```python
from chemsmart.utils.utils import get_prepend_string_list_from_modred_free_format

get_prepend_string_list_from_modred_free_format([1, 2], program="invalid")
# UnboundLocalError: local variable 'modred_string' referenced before assignment
```
See `tests/test_utils_extended.py::TestGetPrependStringListFromModredFreeFormat::test_single_list_invalid_program_crashes_with_unboundlocalerror`.

**Impact:** Low -- only triggered by an invalid `program` argument
combined with the single-list (not list-of-lists) input shape, which
requires a caller-side typo since valid callers only ever pass
`"gaussian"`, `"pymol"`, or `"orca"`. Still, the resulting crash is
far more confusing to debug than the list-of-lists branch's explicit
error.

**Suggested direction:** add the same `else: raise ValueError(...)`
to the single-list branch that the list-of-lists branch already has.

## 50. `sdf2molecule` crashes with `UnboundLocalError` for any input that's neither a `list` nor a `str`

**Location:** `chemsmart/utils/utils.py`, lines 1359-1364.

```python
if isinstance(sdf_lines, list):
    line_elements = sdf_lines
elif isinstance(sdf_lines, str):
    line_elements = sdf_lines.split("\n")

for line in line_elements:
    ...
```

Same shape as bug #49: there's no `else` branch, so if `sdf_lines` is
neither a `list` nor a `str` (e.g. `None`, or an int), `line_elements`
is never assigned, and the following `for line in line_elements:`
raises `UnboundLocalError: local variable 'line_elements' referenced
before assignment` instead of a clear `TypeError`/`ValueError`
describing what was actually expected.

**Reproduce:**
```python
from chemsmart.utils.utils import sdf2molecule

sdf2molecule(12345)
# UnboundLocalError: local variable 'line_elements' referenced before assignment
```
See `tests/test_utils.py::TestSdf2Molecule::test_invalid_type_crashes_with_unboundlocalerror`.

**Impact:** Low -- the function's own docstring already documents the
accepted types as `Union[list, str]`, so this only surfaces if a
caller passes something else by mistake; the confusing crash message
is the only issue.

**Suggested direction:** add an `else: raise TypeError(...)` (or
`ValueError`) describing the expected `Union[list, str]` input.

## 51. `chemsmart`'s top-level `--verbose` flag can never actually be turned off

**Location:** `chemsmart/cli/main.py`, lines 23-33.

```python
@click.option("--verbose", is_flag=True, default=True)
def entry_point(ctx, verbose):
    if verbose:
        debug = True
        stream = True
    else:
        debug = False
        stream = False
```

`--verbose` is declared as a plain `is_flag=True` option (not the
`--verbose/--no-verbose` dual-name form Click supports for toggle
flags), with `default=True`. A bare boolean flag can only ever be
*absent* (uses the default, `True`) or *present* (forces `True`) --
there is no way to pass `False` on the command line. `--no-verbose`
is not recognized (`Error: No such option: --no-verbose`). The
`else` branch handling `verbose=False` is therefore permanently dead
code when invoked through the real CLI; the option exists but cannot
do what its name implies.

**Reproduce:**
```
$ chemsmart --no-verbose sub gaussian ...
Error: No such option: --no-verbose Did you mean --verbose?
```
See `tests/test_main_cli.py::TestMainEntryPoint::test_verbose_false_branch_via_direct_callback_invocation`,
which reaches the `else` branch only by calling `entry_point.callback`
directly (bypassing Click's option parsing entirely) since there is
no real CLI invocation that produces `verbose=False`.

**Impact:** Low -- cosmetic/usability only. Logging always runs in
debug+stream mode regardless of user intent; nobody can quiet it via
this flag.

**Suggested direction:** declare the option as
`@click.option("--verbose/--no-verbose", default=True)` so the
negative form is actually reachable, matching what the flag's name
implies.

## 52. `chemsmart/cli/gaussian/qmmm.py`'s `_populate_charge_and_multiplicity_on_settings` is defined but never called

**Location:** `chemsmart/cli/gaussian/qmmm.py`, lines 385-413.

```python
def _populate_charge_and_multiplicity_on_settings(qs):
    charge = getattr(qs, "charge", None)
    mult = getattr(qs, "multiplicity", None)
    ...
```

This module-level helper is fully defined but never referenced
anywhere else in `qmmm.py`, nor anywhere else in the codebase apart
from `chemsmart/cli/orca/qmmm.py`, which has its own separate,
identically named (and near-identical) copy of the function that
*it* actually calls (`orca/qmmm.py` line 418). The Gaussian copy
appears to be an orphaned leftover -- likely copied from the ORCA
module during development but never wired into the Gaussian QMMM
settings-building flow (`qmmm()`'s charge/multiplicity handling at
lines 294-305 does its own inline `if charge_total is not None: ...`
assignment instead of delegating to this helper).

**Reproduce:** grep confirms zero call sites:
```
$ grep -rn "_populate_charge_and_multiplicity_on_settings" chemsmart/
chemsmart/cli/gaussian/qmmm.py:385:def _populate_charge_and_multiplicity_on_settings(qs):
chemsmart/cli/orca/qmmm.py:418:        _populate_charge_and_multiplicity_on_settings(qmmm_settings)
chemsmart/cli/orca/qmmm.py:460:def _populate_charge_and_multiplicity_on_settings(qs):
```
See `tests/test_gaussian_qmmm_cli.py::TestPopulateChargeAndMultiplicityOnSettings`,
which unit-tests the function directly since no code path reaches it.

**Impact:** Low -- dead code with no behavioral effect (the inline
logic in `qmmm()` already handles charge/multiplicity assignment).
Purely a maintenance/clarity issue: a future edit to the intended
(inline) logic could silently diverge from this unused duplicate
without anyone noticing.

**Suggested direction:** delete the function, or wire it in if it was
meant to replace the inline charge/multiplicity block.

## 53. `PDBFile._infer_element_from_atom_name`'s exception handling around `to_element` is unreachable

**Location:** `chemsmart/io/pdb/pdbfile.py`, lines 283-289 and
317-336 (the `except Exception:` clauses wrapping calls to
`PeriodicTable.to_element`).

```python
if len(cleaned) == 1:
    try:
        return p.to_element(cleaned[0].upper())
    except Exception:
        raise ValueError(
            f"Unable to infer element from atom name '{atom_name}'"
        )
...
if cleaned[:2].upper() not in ambiguous_biomolecular_names:
    try:
        return p.to_element(normalized_two_letter)
    except Exception:
        pass
...
try:
    return p.to_element(candidate)
except Exception:
    raise ValueError(
        f"Unable to infer element from atom name '{atom_name}'"
    )
```

`chemsmart/utils/periodictable.py`'s `PeriodicTable.to_element` only
raises `ValueError` when its input is empty after stripping
non-alphabetic characters (see its own `if not cleaned: raise
ValueError(...)` guards). For any non-empty alphabetic candidate --
which is guaranteed here, since `cleaned` was already validated
non-empty and is built purely from `[A-Za-z]` characters -- it always
returns *something*, falling back to a best-effort guess (e.g.
`to_element("Q")` returns `"Q"`, `to_element("Xy")` returns `"X"`)
rather than raising. This means none of the three `except Exception`
clauses above can ever trigger through genuine parsing of a PDB atom
name; they are only reachable by mocking `to_element` directly (see
`tests/test_converter.py::TestPDBFile::test_infer_element_single_letter_lookup_failure_raises`,
`test_infer_element_two_letter_failure_falls_back_to_single`, and
`test_infer_element_two_and_single_letter_failure_raises`, all of
which patch `pdbfile.p.to_element` with a raising side effect purely
to exercise these lines for coverage).

**Reproduce:**
```python
>>> from chemsmart.utils.periodictable import PeriodicTable
>>> PeriodicTable().to_element("Q")
'Q'
>>> PeriodicTable().to_element("Zz")
'Zz'
```
Neither call raises, even though `"Q"` and `"Zz"` are not real
element symbols -- so `_infer_element_from_atom_name` can silently
return a bogus symbol like `"Q"` for a malformed atom name with blank
element columns, instead of raising `ValueError` as its docstring
promises.

**Impact:** Low-to-moderate -- a PDB file with unrecognized atom
names and blank element columns (77-78) will get an invalid element
symbol silently assigned to the `Molecule`, rather than a clear
`ValueError` at parse time. Downstream code (element lookups, atomic
number/mass calculations) would then fail confusingly later instead
of failing fast here.

**Suggested direction:** either make `PeriodicTable.to_element`
validate its fallback guess against `self.PERIODIC_TABLE` and raise
when it isn't a real element, or have
`_infer_element_from_atom_name` itself validate the returned symbol
against the periodic table before returning it.

## 54. `cli/sub.py`'s `_replace_batch_table_tokens` guards against `label`/`scheme` being `None`, but its only two callers never pass `None`

**Location:** `chemsmart/cli/sub.py`, lines 210-213 and 258-269
(the `if batch_label is not None:` / `if batch_scheme is not None:`
guards inside `_replace_batch_table_tokens`, a closure nested in
`process_pipeline`).

```python
batch_label = batch_entry.get("label")
if batch_label is not None:
    option_map["--label"] = str(batch_label)
    option_map["-l"] = str(batch_label)
...
if batch_label is not None:
    _set_option(args, "--label", "-l", insert_before="pka")

batch_scheme = batch_entry.get("scheme")
if batch_scheme is not None:
    ...
```

`batch_entry` is only ever constructed in four places --
`_create_pka_jobs_from_table` and `_create_pka_jobs_from_molecules`
in both `chemsmart/cli/gaussian/pka.py` and
`chemsmart/cli/orca/pka.py` -- and all four unconditionally set both
`"label"` and `"scheme"` to real, non-`None` values (`label`/
`base_label`/`mol_label` and `row_shared["scheme"]`/`shared["scheme"]`
respectively). Only `"fragment_index"` is conditionally absent (the
table-based callers omit it entirely, unlike the CDXML-fragment
callers), so that particular `is not None` guard is genuinely
reachable, but the `label`/`scheme` ones are not reachable through any
current call path.

**Reproduce:** grep confirms all four `_batch_entry` construction
sites always populate `"label"` and `"scheme"`:
```
$ grep -n '"label"\|"scheme"' chemsmart/cli/gaussian/pka.py chemsmart/cli/orca/pka.py
```
See `tests/test_pka.py::TestSubProcessPipelineDirectInvocation::test_batch_scheme_and_label_none_skip_their_rewrite_blocks`,
which drives these branches directly via a hand-built `batch_entry`
since no real CLI invocation can produce one with `label`/`scheme` set
to `None`.

**Impact:** Low -- purely defensive dead code. If a future call site
ever added a `_batch_entry` without `"label"`/`"scheme"` keys, the
`.get(...)` fallback to `None` combined with these guards would
silently skip rewriting those options rather than raising, which is
probably the desired behavior anyway -- so this isn't a functional
risk, just unreachable code inflating the function's apparent branch
complexity.

**Suggested direction:** no action needed unless a future caller
intentionally omits `label`/`scheme`; if so, this guard already does
the right thing. Otherwise, the guards could be simplified to plain
assignments (dropping the `is not None` checks) once it's confirmed no
call site needs the optional behavior.

## 55. `cli/gaussian/gaussian.py`'s `click_gaussian_qmmm_options` is defined but never applied

**Location:** `chemsmart/cli/gaussian/gaussian.py`, lines 361-490.

```python
def click_gaussian_qmmm_options(f):
    """Common click options for QMMM jobs."""

    @click.option(
        "-hx",
        "--high-level-functional",
        ...
```

This decorator factory defines the full set of QMMM CLI options
(high/medium/low-level functional/basis/force-field, real/int/model
charge and multiplicity, atom-range options, bonded-atoms, and
scale-factors) but is never applied to any command. The actual `qmmm`
subcommand, built by `create_qmmm_subcommand` in
`chemsmart/cli/gaussian/qmmm.py`, defines its own separate,
near-identical set of `@click.option` decorators directly instead of
using this one -- an apparent copy that was never wired in, or a
leftover from before the options were inlined into `qmmm.py`.

**Reproduce:** grep confirms zero call sites:
```
$ grep -rn "click_gaussian_qmmm_options" chemsmart/
chemsmart/cli/gaussian/gaussian.py:361:def click_gaussian_qmmm_options(f):
```
See `tests/test_gaussian_qmmm_cli.py::TestClickGaussianQmmmOptionsIsUnusedDeadCode`,
which applies the decorator directly to a throwaway command purely to
exercise it, since no real command does.

**Impact:** Low -- dead code with no behavioral effect. A future edit
to QMMM's real CLI options (in `qmmm.py`) could drift from this unused
duplicate without anyone noticing, similar to bug #52.

**Suggested direction:** delete `click_gaussian_qmmm_options`, or
refactor `create_qmmm_subcommand` to use it instead of its own inline
copy.

## 56. `cli/gaussian/gaussian.py`'s QMMM-molecule early-conversion block never runs, and its own fallback can't actually succeed for real molecules

**Location:** `chemsmart/cli/gaussian/gaussian.py`, lines 819-858.

```python
    # If the user requested the qmmm subcommand, ensure molecules are
    # represented as QMMMMolecule so the subcommand sees QMMM-specific
    # attributes early (e.g., high_level_atoms, bonded_atoms).
    try:
        if ctx.invoked_subcommand == "qmmm":
            ...
```

Two separate defects in this block:

1. **The guard can never be true.** `qmmm` is always attached as a
   sub-subcommand of a jobtype group (`opt qmmm`, `ts qmmm`, `sp qmmm`,
   etc.) -- never as a direct child of the `gaussian` group itself. In
   Click, `ctx.invoked_subcommand` on a given group's own context only
   ever reflects the *immediate* next command from that context (e.g.
   `"opt"`), never a subcommand of a subcommand. Verified directly:
   ```python
   >>> # outer group -> inner group -> leaf command
   >>> CliRunner().invoke(outer, ["inner", "leaf"])
   outer invoked_subcommand: inner   # never "leaf"
   ```
   So `if ctx.invoked_subcommand == "qmmm":` is always `False` for any
   real `chemsmart run gaussian ... opt qmmm ...` invocation, and this
   whole block never executes. Impact is mitigated because
   `chemsmart/jobs/gaussian/writer.py:360` independently converts the
   job's molecule to `QMMMMolecule` at write time, and
   `chemsmart/cli/gaussian/qmmm.py`'s own callback only ever sets
   plain attributes on the molecule (which works identically on a
   plain `Molecule`), so no QMMM job is known to actually break from
   this -- the intended "early" conversion is simply a no-op.

2. **Even if reached, its own fallback rarely helps.** When
   `QMMMMolecule(molecule=m)` raises, the code retries with
   `QMMMMolecule(**getattr(m, "__dict__", {}))`. For a real `Molecule`
   instance, `__dict__` always contains private/derived keys (e.g.
   `_positions`, `_num_atoms`, `_energy`) that aren't valid
   `Molecule.__init__`/`QMMMMolecule.__init__` parameters, so this
   retry itself raises `TypeError` and falls through to the final
   `except Exception as exc2: ... converted.append(m)` (keeping the
   original molecule) rather than ever actually producing a converted
   `QMMMMolecule` via the fallback path:
   ```
   >>> QMMMMolecule(**vars(Molecule.from_filepath("some.xyz")))
   TypeError: Molecule.__init__() got an unexpected keyword argument '_positions'
   ```

**Reproduce:** see
`tests/test_gaussian_qmmm_cli.py::TestQmmmMoleculeConversionBlockInGaussianGroup`,
which forces `ctx.invoked_subcommand = "qmmm"` directly (a state Click
itself can never produce) to exercise the block's logic in isolation,
including a demonstration that the dict-based fallback needs a
custom-mocked `__init__` to succeed at all, since the real one can't
for a genuine `Molecule`.

**Impact:** Low -- the primary write-time conversion in
`writer.py` already ensures QMMM jobs get a proper `QMMMMolecule`
before input files are written, so no known job output is affected.
This is a latent correctness gap (the stated intent of "early"
QMMM-attribute visibility never happens) rather than an active bug.

**Suggested direction:** either move this conversion into
`create_qmmm_subcommand`'s own callback (where `ctx.invoked_subcommand`
would correctly resolve, or where the subcommand already has direct
access to its own molecules), or remove the block entirely and rely on
`writer.py`'s existing write-time conversion.

## 57. `cli/gaussian/gaussian.py` has several `is not None` guards that are unreachable given earlier validation

**Location:** `chemsmart/cli/gaussian/gaussian.py`, lines 774-790
(the `elif record_index is not None:` arms inside both the
`append_label` and default-label chemsmart-db suffix chains) and lines
811-814 (`if molecule_indices is not None and not isinstance(...)`).

The `structure_id`/`record_id`/`record_index` chain:
```python
if is_chemsmart_db:
    if structure_id is not None:
        label = f"{label}_SID-{structure_id}"
    elif record_id is not None:
        label = f"{label}_RID-{record_id}"
    elif record_index is not None:
        label = f"{label}_RI-{record_index}"
```
is only ever reached when `is_chemsmart_db` is `True`, but an earlier
guard in the same function (lines 552-558) already requires that
`sum([record_index is not None, record_id is not None]) +
(structure_id is not None) == 1` whenever `is_chemsmart_db` is `True`
-- i.e. *exactly one* of the three is guaranteed set, or the function
raises `click.UsageError` before ever reaching this chain. So the
"none of the three matched" fall-through (falling past the final
`elif` with no branch taken) can never happen through any genuine
invocation.

Similarly, `molecule_indices` at line 811 is only ever set to a
non-`None`, non-list value in lock-step with `molecules` becoming
non-list (both come from the same
`return_objects_and_indices_from_string_index` call at lines 803-807),
so whenever the outer `if molecules is not None and not
isinstance(molecules, list):` at line 809 is true via that path,
`molecule_indices` is simultaneously non-list too, making the `False`
arm of line 811's check effectively unreachable in the same way.

**Reproduce:** grep confirms the invariant-establishing guard:
```
$ grep -n "sum(record_selectors)" chemsmart/cli/gaussian/gaussian.py
chemsmart/cli/gaussian/gaussian.py:553:        if sum(record_selectors) + (structure_id is not None) != 1:
```

**Impact:** Low -- purely defensive dead code, consistent with
several similar findings this session (e.g. bugs #42-43, #52, #54).
No behavioral risk since the guarantee they depend on is enforced
immediately upstream in the same function.

**Suggested direction:** no action needed; documenting for awareness
only. These three branch arms remain uncovered by tests since
constructing a call that reaches them would require bypassing the
earlier validation that guarantees they can't occur.

## 58. `GaussianDIASJob.all_molecules_jobs` silently returns `None` for an invalid mode, unlike its sibling properties

**Location:** `chemsmart/jobs/gaussian/dias.py`, lines 312-338.

```python
@property
def all_molecules_jobs(self):
    ...
    if self.mode.lower() == "irc":
        ...
        return jobs
    elif self.mode.lower() == "ts":
        ...
        return [...]
    # no else/raise here
```

`fragment1_jobs` and `fragment2_jobs` both end their identical
`if/elif` chain with an `else: raise ValueError(f"Invalid mode: ...")`
for any mode other than `"irc"`/`"ts"`. `all_molecules_jobs` has the
same two branches but no `else` clause at all, so an invalid
`self.mode` makes it fall through and implicitly return `None`
instead of raising. Any caller that then does
`len(job.all_molecules_jobs)` or iterates over it gets a confusing
`TypeError: object of type 'NoneType' has no len()` /
`'NoneType' object is not iterable` instead of the clear "Invalid
mode" message the other two properties give for the exact same input.

**Reproduce:** see
`tests/test_gaussian_dias_job_unit.py::TestAllMoleculesJobsInvalidMode::test_invalid_mode_silently_returns_none`,
which constructs a `GaussianDIASJob` with `mode="bogus"` and confirms
`job.all_molecules_jobs is None` (whereas the equivalent
`fragment1_jobs`/`fragment2_jobs` calls raise `ValueError`).

**Impact:** Low -- `mode` is only ever set from a CLI option validated
upstream to `"irc"`/`"ts"`, so this isn't reachable through normal
usage. But it's an easy trap for anyone constructing `GaussianDIASJob`
programmatically (e.g. in a script or future refactor) with a typo'd
mode.

**Suggested direction:** add the same
`else: raise ValueError(f"Invalid mode: {self.mode}. Must be 'irc' or 'ts'.")`
clause to `all_molecules_jobs` for consistency with its two siblings.

## 59. `canonicalize_positions`'s equal-mass diatomic sign check can never flip

**Location:** `chemsmart/utils/geometry.py`, lines 352-354 (inside the
diatomic branch of `canonicalize_positions`).

```python
elif masses[0] == masses[1]:
    if rotated[0, 2] > 0:
        rotated[:, 2] *= -1
```

For a two-atom system with equal masses, the centre of mass is the
midpoint of the two atoms, so `shifted[0] = (r0 - r1) / 2 = -vec / 2`
where `vec = r1 - r0` and `z_hat = vec / norm(vec)`. Since
`rotated[0, 2] = shifted[0] . z_hat = -norm(vec) / 2`, which is always
negative (`norm(vec) > 0` for any two distinct points), the condition
`rotated[0, 2] > 0` can never be true for a genuine equal-mass
diatomic -- the flip on line 354 is dead code. (The unequal-mass case
right above it, `masses[0] > masses[1]`, does an *unconditional* flip
and is unaffected by this.)

**Reproduce:**
```python
>>> from chemsmart.utils.geometry import canonicalize_positions
>>> for coords in [[[0,0,0],[1,0,0]], [[1,1,1],[2,2,2]], [[3,1,2],[1,4,0]]]:
...     print(canonicalize_positions([14.0, 14.0], coords))
```
Every case places atom 0 at `z = -bond_length/2` without ever needing
the line 354 flip.

**Impact:** Low -- purely dead code; the function's actual output
(atom 0 always at negative z) is self-consistent and matches its
docstring's "deterministic sign convention" goal regardless. Covered
in `tests/test_geometry.py::TestCanonicalizePositions` only via the
surrounding branches (heavier-first/heavier-second-atom flip tests),
since the equal-mass flip itself is unreachable.

**Suggested direction:** no action needed; the surrounding branches
already guarantee deterministic output. Could be simplified by
removing the dead `if` check and just asserting the invariant, or left
as-is as a defensive check against a future refactor of the COM/z_hat
computation above it.

## 60. `DatabaseQuery.parse_query` has two guards that can never trigger through any real query string

**Location:** `chemsmart/database/query.py`, lines 237-241 (unsupported
operator) and 266-267 (empty clause parts).

```python
if operator not in SUPPORTED_OPERATORS:
    raise ValueError(
        f"Unsupported operator: '{operator}'. "
        f"Supported: {', '.join(sorted(SUPPORTED_OPERATORS))}"
    )
...
if not clause_parts:
    raise ValueError("Empty query string.")
```

Both guards are unreachable in practice:

1. **Unsupported operator.** The operator substring can only ever be
   whatever `query_condition_pattern` (in `chemsmart/utils/repattern.py`)
   itself matched: `==|<=|>=|!=|~|<|>|=` -- exactly the same seven
   operators as `SUPPORTED_OPERATORS`. Any other operator text fails
   the regex match earlier and raises "Invalid condition" instead of
   ever reaching this check.
2. **Empty clause parts.** `_LOGIC_SPLIT_RE.split(...)` on a non-empty
   string never returns an empty list (worst case `['']` for an
   all-whitespace input), and every token in the loop either matches
   `AND`/`OR` (appended), matches a condition (appended), or fails to
   match (raises) -- there is no path where a token is silently
   skipped. So `clause_parts` can never be empty when this check runs.

**Reproduce:** grep confirms the regex's operator alternation exactly
mirrors `SUPPORTED_OPERATORS`:
```
$ grep -n "query_condition_pattern\|SUPPORTED_OPERATORS" chemsmart/utils/repattern.py chemsmart/database/query.py
```
See `tests/test_database.py::TestDatabaseQuery::test_parse_query_unsupported_operator_is_unreachable_defensively`
and `test_parse_query_empty_clause_parts_is_unreachable_defensively`,
both of which monkeypatch the module's compiled regex objects to force
these states directly, since no real query string can produce them.

**Impact:** Low -- purely defensive dead code with helpful error
messages that happen to never fire. No functional risk.

**Suggested direction:** no action needed; these guards are harmless
and provide a clear error message if the regex/operator-set invariant
is ever broken by a future edit (e.g. adding an operator to one list
but not the other).

## 61. `GaussianFileMixin.jobtype`'s setter never actually persists a value

**Location:** `chemsmart/utils/mixins.py`, `GaussianFileMixin.route_object`
(lines ~704-719), `jobtype` getter/setter (lines ~734-758), and
`_get_modredundant_conditions` (lines 623/628, which rely on the
setter to record the detected job type).

```python
@property
def route_object(self):
    try:
        route_object = GaussianRoute(route_string=self.route_string)
        return route_object
    except TypeError as err:
        print(err)

@property
def jobtype(self):
    return self.route_object.jobtype

@jobtype.setter
def jobtype(self, value):
    self.route_object.jobtype = value
```

`route_object` is a plain `@property` (not `@cached_property`), so
every access -- including the one inside the setter -- constructs a
brand-new `GaussianRoute` from `route_string`. `jobtype = "scan"`
therefore sets an attribute on a throwaway `GaussianRoute` instance
that is discarded immediately; the very next `self.jobtype` read
builds yet another fresh `GaussianRoute` and re-derives the job type
purely from `route_string`, never seeing the earlier assignment. This
makes `_get_modredundant_conditions`'s `self.jobtype = "scan"` /
`self.jobtype = "modred"` calls (used to disambiguate a modredundant
block into a scan vs. a frozen-coordinate job) complete no-ops.

**Reproduce:**
```python
>>> from chemsmart.utils.mixins import GaussianFileMixin
>>> class F(GaussianFileMixin):
...     def __init__(self):
...         self.filename = "x"
...     route_string = "modred"
>>> f = F()
>>> f.jobtype = "scan"
>>> f.jobtype
'sp'  # not "scan" -- the assignment never took effect
```
See `tests/test_mixins.py::TestGaussianFileMixin::test_modred_scan_coords_skips_non_scan_lines`,
whose docstring documents this in place of asserting the (non-functional)
jobtype mutation.

**Impact:** Low-to-moderate -- `read_settings()` builds
`GaussianJobSettings(jobtype=self.jobtype, ..., modred=self.modred,
...)` in a single call, evaluating `self.jobtype` *before*
`self.modred` runs (kwarg values are evaluated left-to-right), so even
if the mutation worked it would already be too late for this call
site. In practice this means a `GaussianJobSettings` built via
`read_settings()` from a file with a modredundant scan block gets a
`jobtype` that reflects the route string's own (generic) job type
rather than `"scan"`/`"modred"`, while `modred` itself is still parsed
correctly.

**Suggested direction:** cache `route_object` (e.g. via
`@cached_property`, invalidated appropriately if `route_string` can
change), or have `jobtype`'s setter store the override on `self`
directly (e.g. `self._jobtype_override`) and have the getter check
that first, rather than routing through a freshly-constructed
`GaussianRoute` each time.

## 62. `ConnectivityGrouper._check_isomorphism` is defined but never called

**Location:** `chemsmart/jobs/grouper/connectivity.py`, lines 99-116.

```python
def _check_isomorphism(
    self, idx_pair: Tuple[int, int]
) -> Tuple[int, int, bool]:
    """
    Check graph isomorphism between two molecules for multiprocessing.
    ...
    """
    i, j = idx_pair
    return i, j, self._are_isomorphic(self.graphs[i], self.graphs[j])
```

Despite its docstring describing it as "multiprocessing-compatible" (implying it's meant to be dispatched via `joblib.Parallel`/`delayed` for the pairwise isomorphism checks), `group()`'s actual pairwise-check loop calls `self._are_isomorphic(...)` directly in a plain Python `for` loop (see lines 180-184), never `_check_isomorphism`. Only the graph-conversion step (`to_graph_wrapper` via `Parallel`) is actually parallelized; the isomorphism checks themselves run serially regardless of `num_procs`.

**Reproduce:** grep confirms zero call sites:
```
$ grep -rn "_check_isomorphism" chemsmart/
chemsmart/jobs/grouper/connectivity.py:99:    def _check_isomorphism(
```
See `tests/test_groupers.py::TestGroupers::test_check_isomorphism_direct`, which unit-tests the function directly since no code path reaches it.

**Impact:** Low-to-moderate -- purely a performance gap, not a correctness bug: pairwise isomorphism checks (O(n^2) for n molecules) never benefit from `num_procs > 1`, unlike the graph-conversion step. For large conformer sets this could be a meaningful missed optimization, but results are still correct.

**Suggested direction:** either wire `_check_isomorphism` into the pairwise-check loop via `Parallel(n_jobs=self.num_procs)(delayed(self._check_isomorphism)(pair) for pair in indices)` (mirroring the graph-conversion step), or delete it if serial checking is intentional (e.g. because Union-Find's sequential `union()` calls don't parallelize cleanly).

## 63. `cli/mol/align.py` has several validation checks that duplicate upstream guarantees and are unreachable through any input

**Location:** `chemsmart/cli/mol/align.py`:
- lines 68-77 (the `if index is None: if ... elif ... elif ...:` chain's third branch, `75->79`)
- lines 136-141 (a second "directory without filetype" `click.BadParameter`)
- lines 181-185 (a second "no input files" `click.BadParameter`)
- lines 187-189 (`if not isinstance(molecules, list): molecules = list(molecules) if molecules else []`)
- lines 202-208 (`base_label = "molecules"` fallback when `base_file_for_label` is falsy)

Each of these re-checks a condition already guaranteed impossible earlier in the *same function* or by `chemsmart/cli/mol/mol.py`'s own group callback:

- The lines 61-65 guard already raises if `directory` is truthy without `filetype`, so by the time execution reaches the `if directory: ... else: raise ...` at line 136-141, `filetype` is guaranteed truthy whenever `directory` is -- the `else` can never execute.
- The line 56-59 guard already raises if neither `filenames` nor `directory` is set, so the final `else` at line 181-185 (reached only when neither `if directory:` nor `elif filenames:` matched) can never execute either -- confirmed empirically: even calling the undecorated callback directly with a hand-built `ctx.obj` (bypassing Click and `mol.py`'s own group logic entirely) still hits the *first* guard before ever reaching this one.
- For the same reason, once execution passes line 56-59, at least one of `directory`/`filenames` is truthy, and combined with the 61-65 guard (directory implies filetype), the `if/elif/elif` index-defaulting chain's third branch's False arm (line 75, falling through to line 79) can never actually occur -- one of the three branches always matches.
- `molecules` is initialized as `molecules = []` at line 79 and only ever mutated via `.extend(...)` afterwards, so it is always already a `list` by line 188; the `isinstance` guard's True arm can't fire.
- `base_file_for_label` is unconditionally assigned inside both live branches of the `if directory: ... elif filenames: ...` chain (`matched_files[0]` and `filenames[0]` respectively) before the label-generation block runs, so the `"molecules"` fallback at line 208 can only fire from the (also dead) `else` branch at line 181.

Additionally, `mol.py`'s own group callback (`chemsmart/cli/mol/mol.py`) independently validates and populates `ctx.obj["filenames"]`/`ctx.obj["directory"]` *before* `align()` ever runs, so the CLI can never even reach `align()` in a state where these dead branches' guarded conditions hold true.

**Reproduce:** see `tests/test_mol_cli.py::TestMolCLIAlignDirectInvocation`, which calls `align`'s undecorated callback directly (via `inspect.unwrap`) with a hand-built `ctx.obj`, bypassing both Click's parsing and `mol.py`'s group logic, to reach the *reachable* dead branches (the second directory/filetype guard, the `isinstance(filenames, str)` normalization, and the "may not select enough structures" hint's False arm). Attempting to also reach the "no input files" second guard (line 183) via the same direct-invocation technique still hits the *first* guard (line 56-59) instead, proving it's unreachable via any input, not just via the real CLI.

**Impact:** Low -- entirely redundant defensive code with no behavioral effect; consistent with several other findings this session (e.g. bugs #42-43, #52, #60) where copy-pasted or belt-and-suspenders validation duplicates an invariant already enforced earlier in the same call path.

**Suggested direction:** no action needed; could be simplified by deleting the unreachable duplicate checks (139, 183, and the `isinstance`/`molecules` normalizations) now that they're confirmed dead, trusting the earlier guards instead.

---

## 64. `ThermochemistryJob.__init__`'s `molecule`/fallback label branches are unreachable dead code

**Location:** `chemsmart/jobs/thermochemistry/job.py:108-114`

```python
if label is None:
    if filename is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
    elif molecule is not None:
        label = molecule.get_chemical_formula(empirical=True)
    else:
        label = "thermochemistry_job"
```

`__init__` raises `ValueError: 'filename' must be provided.` a few lines earlier (line 72-73, see bug #1) whenever `filename is None`, so by the time this block runs `filename` is guaranteed truthy. The `elif molecule is not None:` and `else:` arms can never execute through any real construction path -- `if filename is not None:` always matches.

**Reproduce:** confirmed empirically -- there is no way to call `ThermochemistryJob(...)` with `molecule` set and `filename=None` without hitting the earlier `ValueError` first. See `tests/test_thermochemistry_job_unit.py::TestThermochemistryJobConstruction` for the reachable label-from-filename behavior.

**Impact:** None -- purely redundant defensive code; consistent with bug #1, which already documents that the `filename` guard was added later and makes this fallback chain moot.

**Suggested direction:** no action needed; could be simplified to unconditionally derive `label` from `filename` (dropping the `elif`/`else`), now that `molecule`-only construction is confirmed unreachable.

---

## 65. `cli/gaussian/irc.py`'s `irc()` command has five `is not None` update-guards that are unreachable through the real CLI

**Location:** `chemsmart/cli/gaussian/irc.py:62-71`

```python
if recalc_step is not None:
    irc_settings.recalc_step = recalc_step
if maxpoints is not None:
    irc_settings.maxpoints = maxpoints
if maxcycles is not None:
    irc_settings.maxcycles = maxcycles
if stepsize is not None:
    irc_settings.stepsize = stepsize
if flat_irc is not None:
    irc_settings.flat_irc = flat_irc
```

The corresponding Click options (`chemsmart/cli/gaussian/gaussian.py`'s `click_gaussian_irc_options`, lines ~125-178) all declare concrete, non-`None` defaults: `--recalc-step` defaults to `6`, `--maxpoints` to `512`, `--maxcycles` to `128`, `--stepsize` to `20`, and `--flat-irc/--no-flat-irc` to `False`. Click therefore always passes a concrete value for these five parameters -- never `None` -- so each `is not None` check is always `True` through any real CLI invocation; the implicit "skip" (False) arm can never execute. This differs from the two neighboring guards for `predictor` and `recorrect` (and `direction`), whose Click options do default to `None` and whose False arms are genuinely reachable (see `tests/test_gaussian_cli.py::TestGaussianCLIIrcCommand::test_basic_irc_job_creation`, which omits `-pt`/`-rc` and leaves `irc_settings.predictor`/`recorrect` at the project defaults).

**Reproduce:** confirmed via `coverage report -m` on `chemsmart/cli/gaussian/irc.py` after running the full `TestGaussianCLIIrcCommand` suite -- branches `62->64`, `64->66`, `66->68`, `68->70`, `70->72` (the False arm of each guard) remain unexecuted no matter which combination of real CLI flags is supplied, because Click never produces `None` for these five parameters.

**Impact:** None -- purely redundant defensive code that mirrors a pattern (guards duplicating an invariant already enforced upstream, here by Click's own default-value mechanism) seen elsewhere in this file, e.g. bugs #57 and #63.

**Suggested direction:** no action needed; the guards are harmless. Could be simplified by assigning these five attributes unconditionally (dropping the `is not None` checks) since Click guarantees a value, or by changing the five defaults to `None` if "unset by the user" needs to be distinguishable from "user explicitly passed the default" for merge-with-project-settings purposes -- but that would be a behavior change requiring product input, not a pure test-coverage fix.

---

## 66. `GenGenECPSection.from_bse_api`'s header-line falsy check is unreachable dead code

**Location:** `chemsmart/io/gaussian/gengenecp.py:312-316`

```python
header_block = heavy_atoms_gengenecp_basis_blocks[0]
for line in header_block:
    if line:
        genecp_string += line + "\n"
```

`heavy_atoms_gengenecp_basis_blocks` comes from `content_blocks_by_paragraph` (`chemsmart/utils/utils.py:230-247`), which groups a line list with `itertools.groupby(string_list, lambda x: x == "")` and explicitly discards groups where the key is `True` (`if not k`). By construction, no returned block can ever contain an empty-string element -- every line in `header_block` is truthy. The `if line:` check's False arm can therefore never execute.

**Reproduce:** confirmed by inspecting `content_blocks_by_paragraph`'s `groupby`/`if not k` filter -- any all-blank group is filtered out before blocks are returned, so a non-blank block (like `header_block`) cannot contain a blank entry. See `tests/test_GaussianGenECP.py::TestGenGenECPSectionEdgeCases::test_string_blocks_splits_on_blank_lines` for the general splitting behavior this relies on.

**Impact:** None -- purely redundant defensive code with no behavioral effect.

**Suggested direction:** no action needed; could be simplified to `genecp_string += "\n".join(header_block) + "\n"` (or similar) now that the per-line truthiness check is confirmed dead.

---

## 67. `GaussianRoute.get_additional_solvent_options`'s trailing `return None` is unreachable dead code

**Location:** `chemsmart/io/gaussian/route.py:529-571`

```python
def get_additional_solvent_options(self):
    if "scrf" in self.route_string:
        scrf_string = ""
        for each_input in self.route_inputs:
            if "scrf" in each_input:
                scrf_string = each_input
                ...
                if len(scrf_line_elements) <= 2:
                    return None
                ...
                return (
                    ",".join(filtered_elements)
                    if filtered_elements
                    else None
                )
    return None
```

`self.route_inputs` is `self.route_string.split()` (whitespace-separated tokens), so any contiguous 4-character substring like `"scrf"` found in `self.route_string` must lie entirely within a single token -- whitespace can't split it apart. That means whenever the outer `if "scrf" in self.route_string:` guard passes, the `for` loop's `if "scrf" in each_input:` is guaranteed to match on at least one token. Every code path inside that inner `if` block explicitly returns (either `return None` at the length check, or the final `return (...)` expression) -- there is no path where the loop body executes without returning. Consequently the `for` loop can never complete a full pass and fall through to the function-level `return None` at line 571 while the outer guard is true; the only way to reach line 571 is when `"scrf" not in self.route_string`, which is already handled by the identical early `return None` one line 571 shares indentation with (functionally redundant with the guard's implicit "no scrf" case).

**Reproduce:** confirmed by exhaustive reasoning over the tokenization invariant above (`str.split()` cannot separate a substring with no internal whitespace across two tokens); could not construct any route string where `"scrf" in self.route_string` is true but no single token in `self.route_inputs` contains `"scrf"`. See `tests/test_GaussianIO.py::TestRouteString::test_solvent_id_none_when_scrf_has_no_solvent_keyword` and `test_solvent_id_appends_read_when_read_precedes_solvent` for the reachable branches in the sibling `get_solvent_id` method, which has the identical structural pattern.

**Impact:** None -- purely redundant defensive code with no behavioral effect.

**Suggested direction:** no action needed; the trailing `return None` is harmless (and still correctly handles the "no scrf at all" case via the outer guard).

---

## 68. `remove_phantom_metal_carbons` in `utils/io.py` has ~46 lines of orphaned, unreachable code after its `return`

**Location:** `chemsmart/utils/io.py:810-924`

`remove_phantom_metal_carbons` returns unconditionally at line 877 (`return new_mol, new_metal_idxs`). Immediately after that `return`, still indented as part of the same function body, lines 879-924 contain an orphaned triple-quoted docstring followed by a second, complete block of logic that de-aromatizes 5-membered all-carbon rings into `[cH-]` Cp anions (setting a formal charge of -1 and an explicit H on one ring atom, then de-aromatizing the ring bonds) and returns `rw.GetMol()`. This second block can never execute -- Python returns from the function at line 877 before reaching it. It reads like a second function body (perhaps a `dearomatize_cp_ring`-style helper) that got merged into `remove_phantom_metal_carbons` during a refactor, losing its own `def` line in the process.

```python
    new_mol.UpdatePropertyCache(strict=False)
    return new_mol, new_metal_idxs                      # <-- line 877, always returns here

    """
    RDKit cannot sanitize a neutral aromatic 5-member carbon ring (c1cccc1).
    ...
    """
    rw = Chem.RWMol(mol)                                 # <-- unreachable from here on
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        ...
    return rw.GetMol()                                   # <-- line 924, dead
```

**Reproduce:** confirmed by code inspection (an unconditional `return` at the same indentation level unconditionally exits the function) and empirically by the coverage-improvement agent that found it while adding tests for this file -- no test input can reach lines 879-924 through the public `remove_phantom_metal_carbons` API. See `tests/test_io_utils_organometallic.py` (added while raising this file from 77% to 98% coverage) for the reachable behavior of this function.

**Impact:** None on current behavior (dead code, no effect), but the *intended* Cp-dearomatization logic itself never runs anywhere in the codebase -- if some other caller was meant to invoke this as a separate helper (e.g., to handle ChemDraw's neutral-aromatic-Cp drawing convention before `attach_eta_bonds_for_cp_rings` runs), that functionality is effectively missing/silently absent.

**Suggested direction:** needs a decision, not a pure test-coverage fix: either (a) delete the dead lines 879-924 if the logic is truly obsolete/superseded by `attach_eta_bonds_for_cp_rings`'s own dearomatization (lines 1054-1078 in the same file already do something similar inline), or (b) extract lines 879-924 into their own properly-named function (e.g. `_dearomatize_neutral_cp_rings`) and wire it into the actual CDX-parsing pipeline in `chemsmart/io/file.py` if the Cp-anion handling it implements is still needed for some ChemDraw drawing style not currently covered.

---

## 69. `read_molecular_job_yaml`'s qmmm-fallback `except` block is unreachable -- a bad key in `gas_config` always crashes instead of falling back

**Location:** `chemsmart/jobs/settings.py:310-331`

```python
for job in gas_phase_jobs:  # jobs using gas config
    all_project_configs[job] = default_config.copy()
    all_project_configs[job]["jobtype"] = job
    all_project_configs[job] = update_dict_with_existing_keys(
        all_project_configs[job], gas_config
    )                                              # <-- line 315-317, unguarded
    try:
        # Try updating with gas_config first
        all_project_configs[job] = update_dict_with_existing_keys(
            all_project_configs[job], gas_config
        )                                          # <-- line 320-322, identical, inside try
    except Exception as e:
        logger.warning(
            f"Updating job '{job}' with gas_config failed ({e}). "
            f"Falling back to qmmm_config."
        )
        all_project_configs[job] = update_dict_with_existing_keys(
            all_project_configs[job], qmmm_config
        )
```

There are two back-to-back, identical calls to `update_dict_with_existing_keys(all_project_configs[job], gas_config)` -- one unguarded (lines 315-317) immediately followed by the *same* call wrapped in a `try` (lines 320-322) whose `except` is meant to fall back to `qmmm_config` on failure. Since `update_dict_with_existing_keys` raises `ValueError` for any key in `gas_config` not already present in `all_project_configs[job]` (see `chemsmart/utils/utils.py:1191-1216`), the *first*, unguarded call always raises before the `try` block is ever entered whenever `gas_config` contains a bad key. The `except` branch's `logger.warning` and qmmm-fallback (lines 323-329) can therefore never execute -- the `ValueError` from the unguarded duplicate propagates straight out of `read_molecular_job_yaml` instead.

**Reproduce:** `tests/test_GaussianSettings.py::TestGaussianJobSettings::test_get_settings_from_yaml_gas_config_bad_key_raises` -- calling `read_molecular_job_yaml` with a `gas` section containing a key not in `defaults.yaml` (e.g. a QMMM-only key like `high_level_functional`) raises `ValueError: Keyword 'high_level_functional' is not in list of keywords ...` from line 315, even though a `qmmm` section with a valid fallback value is also present in the project YAML.

**Impact:** Medium -- this appears to be the intended mechanism for QMMM project YAMLs to let `gas`/`solv` sections use QMMM-specific keys with a fallback to a separate `qmmm` section (mirroring how `chemsmart/jobs/settings.py:355-368`'s explicit `qmmm` block handling works via direct key assignment, not `update_dict_with_existing_keys`). As written, any project YAML relying on this fallback instead crashes with a raw `ValueError` naming the bad key, rather than silently falling back. In practice this may be low-impact if no real project YAMLs currently rely on the fallback (the existing `qmmm.yaml` test fixture's `gas`/`solv` sections only use keys already in `defaults.yaml`, so they don't trigger this path at all).

**Suggested direction:** delete the unguarded duplicate call at lines 315-317 (it's a copy-paste duplicate of the one inside the `try`) so the `try`/`except` actually gets a chance to run and the qmmm fallback becomes reachable, if that fallback behavior is still wanted; otherwise remove the dead `try`/`except`/fallback entirely and let the `ValueError` propagate directly with a clearer message, if failing fast is the intended behavior.

---

## 70. `sort_structure_dicts_by_energy`'s `sorted_frames` falsy branch is unreachable dead code

**Location:** `chemsmart/database/utils.py:389-413`

```python
def sort_structure_dicts_by_energy(db_file, struct_dicts):
    if not struct_dicts:
        return struct_dicts
    frames = [
        {...}
        for s in struct_dicts
    ]
    sorted_frames = sort_frames_by_energy(frames)
    if sorted_frames:
        ...
    else:
        primary_mb = None
```

`frames` is built with a list comprehension over `struct_dicts` -- one entry per input dict, so `len(frames) == len(struct_dicts)`. `struct_dicts` is guaranteed non-empty by the early `if not struct_dicts: return struct_dicts` guard, so `frames` is always non-empty by the time `sort_frames_by_energy(frames)` is called. `sort_frames_by_energy` itself either returns its input unchanged (`return frames` early exit when no `(method, basis)` pairs are covered) or `sorted(frames, key=sort_key)`, both of which preserve length. So `sorted_frames` can never be empty/falsy at this point, and the `else: primary_mb = None` branch can never execute.

**Reproduce:** `tests/test_database.py::TestDatabaseUtilities::test_sort_structure_dicts_by_energy_empty_input` and `test_sort_structure_dicts_no_energy` cover the two adjacent reachable cases (empty `struct_dicts` hits the earlier guard; non-empty `struct_dicts` with zero energies still produces one frame per input dict, so `sorted_frames` stays non-empty and takes the `if` branch, not the `else`).

**Impact:** None -- purely redundant defensive code; the `if sorted_frames:` guard's `else` was presumably written defensively without noticing the length-preservation invariant of the two functions it wraps.

**Suggested direction:** no action needed; could be simplified by removing the `if sorted_frames: ... else: primary_mb = None` branch and unconditionally taking the `if` body, now that `sorted_frames` is confirmed always non-empty when reached.

---

## 71. `ORCAOutput.final_scf_energy`'s and `final_structure`'s single-point/abnormal-termination branches, and all of `_get_sp_scf_energy`, are unreachable dead code

**Location:** `chemsmart/io/orca/output.py:333-343, 973-987, 1104-1128`

```python
@cached_property
def optimized_output_lines(self):
    """... FOR SP CALCULATION, THIS WILL BE EMPTY! """
    optimized_output_lines = []
    for i, line_i in enumerate(self.contents):
        if "THE OPTIMIZATION HAS CONVERGED" in line_i:
            optimized_output_lines = list(self.contents[i:])
    return optimized_output_lines

@property
def _get_sp_scf_energy(self):
    if self.optimized_output_lines is None:
        for line in self.contents:
            ...
            return energy_in_hartree

@property
def final_scf_energy(self):
    if self.optimized_output_lines is not None:
        return self._get_optimized_scf_energy()
    return self._get_sp_scf_energy()

@property
def final_structure(self):
    if self.optimized_output_lines is not None:
        return self.optimized_structure
    try:
        return self.last_structure
    except (ValueError, IndexError):
        return self._get_molecule_from_sp_output_file()
```

`optimized_output_lines` always returns a `list` -- it's initialized to `[]` and only ever reassigned to another list via `list(...)`, never to `None`. So `self.optimized_output_lines is None` is always `False`, and `self.optimized_output_lines is not None` is always `True`. This affects two independent properties that both guard on it the same (broken) way:

1. `final_scf_energy`/`_get_sp_scf_energy`: `_get_sp_scf_energy`'s entire body is dead -- the docstring-implied "single point energy" fallback can never run, and the property implicitly returns `None` in every case. `final_scf_energy` always takes the `_get_optimized_scf_energy()` branch and never calls `_get_sp_scf_energy` (whose `else`-branch call site at line 1128 is itself unreachable). For a genuine single-point job, `optimized_output_lines` is `[]` (per its own docstring), so `_get_optimized_scf_energy` iterates zero times and `final_scf_energy` silently returns `None` instead of the SP energy the dead `_get_sp_scf_energy` was written to supply. As a bonus latent bug, even if line 1128 were somehow reached, `self._get_sp_scf_energy()` would fail: `_get_sp_scf_energy` is itself decorated `@property`, so `self._get_sp_scf_energy` already invokes it (returning `None`), and appending `()` would then try to call `None`, raising `TypeError`.

2. `final_structure`: always takes `return self.optimized_structure` and never reaches the `try`/`except` block, so the `except (ValueError, IndexError): return self._get_molecule_from_sp_output_file()` fallback (clearly intended for abnormally-terminated jobs) is dead code. For a job that crashed before producing any usable coordinate block, `self.optimized_structure` -> `self._get_optimized_final_structure()` can raise `ValueError` (via `CoordinateBlock([]).molecule` finding no symbols) that propagates straight out of `final_structure` uncaught, instead of falling back to `_get_molecule_from_sp_output_file()` as designed.

**Reproduce:** `tests/test_ORCAIO.py::TestORCAOutputDirectPropertyCoverage::test_final_scf_energy_and_single_point_energy_for_sp_job` -- for `water_dlpno_ccsdt_sp.out` (a pure single-point job), `oo.optimized_output_lines == []` and `oo.final_scf_energy is None`, even though the file contains a `Total Energy       :` line that `_get_sp_scf_energy` was clearly written to parse. `final_energy` happens to still work because it separately falls back to `self.single_point_energy` (parsed from `FINAL SINGLE POINT ENERGY`) whenever `final_scf_energy` is `None`. Separately, `test_abnormal_termination_all_structures_and_final_structure` shows `ORCAOutput(filename=".../GTOInt_error.out").optimized_structure` (and therefore `.final_structure`) raising an uncaught `ValueError` rather than falling back gracefully.

**Impact:** Low-to-moderate. `final_energy` papers over the `final_scf_energy` bug via its own fallback, but `final_scf_energy` itself is silently broken for every single-point ORCA job -- any caller relying on `final_scf_energy` directly (rather than `final_energy`) gets `None` instead of the SCF energy. `final_structure`'s dead fallback is more consequential: any abnormally-terminated job without a parseable coordinate block raises an uncaught `ValueError` from `final_structure`/`molecule` instead of the intended graceful fallback to the input-echo-derived structure.

As a further consequence, `_get_molecule_from_sp_output_file` (`chemsmart/io/orca/output.py:1011-1034`) and the helper it calls, `_get_input_structure_in_output` (`chemsmart/io/orca/output.py:1036-1069`), have no other call site in the codebase besides `final_structure`'s dead `except` block, so they are transitively unreachable through the only path that exists to them.

**Suggested direction:** either change `optimized_output_lines` to return `None` (not `[]`) when no `"THE OPTIMIZATION HAS CONVERGED"` line is found, or change all three `is None` / `is not None` checks (in `final_scf_energy`, `_get_sp_scf_energy`, and `final_structure`) to test truthiness (`if self.optimized_output_lines:` / `if not self.optimized_output_lines:`) instead of identity against `None`. Also fix the call site `self._get_sp_scf_energy()` to `self._get_sp_scf_energy` (drop the parens) once the branch becomes reachable, since it's a `@property`.

## 72. `ORCANEBOutput.product`'s except clause catches the wrong exception type and can never fire

**Location:** `chemsmart/io/orca/output.py:3571-3577`

```python
@property
def product(self):
    try:
        return self._get_geometries()[1]
    except (TypeError, ValueError):
        # product geometry may not be found for a free end NEB
        return None
```

`_get_geometries()` always returns a plain `list` of `Molecule` objects (possibly empty, possibly length 1 for a free-end NEB with only a `REACTANT (ANGSTROEM)` block and no `PRODUCT (ANGSTROEM)` block). Indexing a list with `[1]` when the list has fewer than 2 elements raises `IndexError` -- never `TypeError` or `ValueError`. So the `except (TypeError, ValueError):` clause, clearly intended to gracefully handle exactly this "free end NEB with no product geometry" case (per its own comment), can never actually catch anything.

**Reproduce:** `tests/test_ORCAIO.py::TestORCANEB::test_product_raises_indexerror_when_only_reactant_present` -- a synthetic output with only a `REACTANT (ANGSTROEM)` block makes `_get_geometries()` return a length-1 list, and `.product` raises an uncaught `IndexError` instead of returning `None`.

**Impact:** Real users hit this: any free-end NEB job (no product geometry) makes `.product` raise instead of returning `None` as documented by the surrounding comment.

**Suggested direction:** change the caught exception type to `IndexError` (or broaden to `(IndexError, TypeError, ValueError)` to stay defensive against other malformed-input failure modes too).

## 73. `ORCAOutput._get_all_structures` is defined but never called

**Location:** `chemsmart/io/orca/output.py:856-871`

```python
def _get_all_structures(self):
    """Extract all Cartesian coordinate blocks from the ORCA output.
    This does not however include energy and forces."""
    structures = []
    for i, line in enumerate(self.contents):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            ...
    return structures
```

A grep across `chemsmart/` and `tests/` finds no call site for `self._get_all_structures()` (or `ORCAOutput._get_all_structures`) anywhere outside its own definition. Its functionality is a strict subset of `_get_all_orientations` (used by the actively-maintained `all_structures` property), which does the same coordinate-block scan but additionally returns raw `np.array` coordinate data (rather than fully-built `Molecule` objects) for use with the energies/forces/PBC assembly pipeline in `all_structures`.

**Reproduce:** `grep -rn "_get_all_structures" chemsmart/ tests/` matches only the method's own `def` line.

**Impact:** None -- dead code with no runtime effect; presumably an earlier implementation of the same functionality later superseded by `_get_all_orientations` + `all_structures`, left behind after the refactor.

**Suggested direction:** remove `_get_all_structures`, or if some future caller needs a plain `Molecule`-list view without energies/forces, wire it in explicitly.

## 74. `ORCAOutput._get_input_structure_coordinates_block_in_output` is defined but never called

**Location:** `chemsmart/io/orca/output.py:170-202`

```python
def _get_input_structure_coordinates_block_in_output(self):
    """In ORCA output file, the input structure
    is rewritten and for single points,
    is same as the output structure.
    ...
    """
    coordinates_block_lines_list = []
    pattern = re.compile(orca_input_coordinate_in_output)
    for i, line in enumerate(self.contents):
        if "INPUT FILE" in line:
            ...
    cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
    return cb
```

`input_coordinates_block` (the only plausibly-related public property, right above this method) actually calls `self._get_first_structure_coordinates_block_in_output()` -- a differently-named, differently-implemented method a few lines below (`chemsmart/io/orca/output.py:204-233`) that scans for `"CARTESIAN COORDINATES (ANGSTROEM)"` rather than `"INPUT FILE"`. A grep across `chemsmart/` and `tests/` finds no call site for `_get_input_structure_coordinates_block_in_output` anywhere.

**Reproduce:** `grep -rn "_get_input_structure_coordinates_block_in_output" chemsmart/ tests/` matches only the method's own `def` line.

**Impact:** None -- dead code with no runtime effect. Likely an earlier implementation (parsing the "INPUT FILE" echo section) superseded by `_get_first_structure_coordinates_block_in_output` (parsing the "CARTESIAN COORDINATES (ANGSTROEM)" section directly), with the old version left behind.

**Suggested direction:** remove `_get_input_structure_coordinates_block_in_output`.

---

## 75. `XTBMainOut.only_rot_calc` searches for the wrong SETUP key, so it can never return anything but `None` for real xTB output

**Location:** `chemsmart/io/xtb/file.py:897-904`

```python
@property
def only_rot_calc(self):
    """compute only rotational contributions to Hessian,
    rather than the full vibrational analysis."""
    only_rot = self._get_setup_information("only rotational calc.")
    if only_rot:
        return only_rot.lower() == "true"
    return None
```

`_get_setup_information` does a literal substring search for `"only rotational calc."` inside the Hessian `SETUP` block. Real xTB output, however, prints this key as `only rotor calc.` (see e.g. `tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out`, which contains the line `:  only rotor calc.                    false      :` but never the string `"only rotational calc."`). Since the two strings never match, `_get_setup_information` always returns `None` for this keyword, so `only_rot` is always falsy and the property always takes the `return None` branch -- the `return only_rot.lower() == "true"` line is unreachable dead code for any real xTB output, and the property can never actually report whether a run was a rotation-only Hessian calculation.

**Reproduce:** `tests/test_XTBIO.py::TestXTBMainOutSyntheticEdgeCases::test_only_rot_calc_key_never_matches_real_xtb_output` confirms that `co2_ohess.out` contains `"only rotor calc."` but not `"only rotational calc."`, and that `only_rot_calc` returns `None` rather than `False` even though the run's SETUP block explicitly states `only rotor calc. false`.

**Impact:** Low but real -- any caller relying on `only_rot_calc` to distinguish a rotation-only Hessian run from a full vibrational analysis always gets `None`, never `True`/`False`, regardless of what the actual xTB output says.

**Suggested direction:** change the search keyword from `"only rotational calc."` to `"only rotor calc."` to match real xTB output.

---

## 76. `XTBMainOut.molecular_dipole_qonly`'s search loop can never advance past its first element

**Location:** `chemsmart/io/xtb/file.py:571-579`

```python
@property
def molecular_dipole_qonly(self):
    """Charge only dipole, computed only from atomic partial charges
    (electrostatic contribution)."""
    if self.molecular_dipole_lines is not None:
        for line in self.molecular_dipole_lines:
            if line.startswith("q only:"):
                return np.array([float(x) for x in line.split()[-3:]])
    return None
```

`molecular_dipole_lines` always builds its list as `self.contents[i + 2 : i + 4]`, i.e. exactly the two lines immediately following the `molecular dipole:` heading and its column-label row. Real xTB always prints the `q only:` row first and the `full:` row second within that block (see e.g. `tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out`), so whenever `molecular_dipole_lines` is non-`None` (i.e. non-empty), its very first element already starts with `"q only:"`. The `for` loop in `molecular_dipole_qonly` therefore always returns on its first iteration -- it can never advance to a second element, and the loop's own back-edge (continuing past a non-matching first line) is unreachable for any real xTB output. This differs from the sibling properties `molecular_dipole_full`/`total_molecular_dipole_moment`, which search for `"full:"` and *do* need to skip past the `q only:` row first, so their loops don't have this issue.

**Reproduce:** `tests/test_XTBIO.py::TestXTBMainOutSyntheticEdgeCases::test_dipole_full_and_total_none_when_full_line_truncated` shows `molecular_dipole_qonly` matching immediately, while `molecular_dipole_full`/`total_molecular_dipole_moment` are the ones that need a full 2-line window to find their marker.

**Impact:** None -- purely a coverage/dead-code observation; the loop still produces the correct result, it just never needs more than one iteration given the fixed, guaranteed line ordering.

**Suggested direction:** no action needed; if desired, the loop could be simplified to a direct index into `self.molecular_dipole_lines[0]` instead of a `for`/`startswith` scan, now that the ordering invariant is confirmed.

---

## 77. `IterateJobRunner.run_combinations`'s "no results at all" summary-skip branch is unreachable dead code

**Location:** `chemsmart/jobs/iterate/runner.py:290-292, 400-430`

```python
if not combinations:
    logger.warning("No combinations to process.")
    return []
...
# Check for missing results (crashes that didn't write to queue)
for comb in combinations:
    if comb.label not in results_dict:
        failed_labels.append(comb.label)
        results_dict[comb.label] = None
    elif (
        results_dict[comb.label] is None
        and comb.label not in timed_out_labels
    ):
        if comb.label not in failed_labels:
            failed_labels.append(comb.label)
...
if successful_labels or timed_out_labels or failed_labels:
    logger.info("=" * 40)
    logger.info("       SUMMARY OF RESULTS")
    ...
```

An early guard at the top of `run_combinations` already returns `[]` when `combinations` is empty, so by the time the summary-printing block runs, `combinations` is guaranteed non-empty -- at least one `IterateCombination` exists. For every combination, the bookkeeping loop guarantees its label ends up in exactly one of three buckets: `successful_labels` (built from any `results_dict` entry with a non-`None` molecule), `timed_out_labels` (populated directly by the watchdog when a worker is killed for exceeding its timeout), or `failed_labels` (populated either when the label never appears in `results_dict` at all -- a crashed worker -- or when it maps to an explicit `None` and wasn't a timeout). There is no path for a combination to avoid all three: a timeout always adds to `timed_out_labels`; a missing/crashed result always adds to `failed_labels`; an explicit `None` result either adds to `failed_labels` or (if a duplicate label already added it) leaves `failed_labels` non-empty regardless. So `successful_labels or timed_out_labels or failed_labels` is always `True` when this line is reached, and the branch where all three are falsy (skipping the summary block entirely) can never execute.

**Reproduce:** confirmed by reasoning through every code path that appends to `results_dict`/`successful_labels`/`timed_out_labels`/`failed_labels` -- combined with the `if not combinations: return []` early guard, no combination can reach the summary check without having landed in one of the three lists. See `tests/test_iterate_run_combinations_unit.py` for direct coverage of each of the three population paths (timeout, missing/crashed, explicit-None-not-timed-out) using fake `multiprocessing.Process`/`Manager` objects.

**Impact:** None -- purely redundant defensive code with no behavioral effect.

**Suggested direction:** no action needed; could be simplified by removing the `if successful_labels or timed_out_labels or failed_labels:` guard and printing the summary unconditionally, now that at least one of the three is confirmed always non-empty whenever this code runs.

---

## 78. `cli/orca/qmmm.py`'s `qmmm()` command has three branches unreachable through any real CLI invocation

**Location:** `chemsmart/cli/orca/qmmm.py:289-297, 317-318, 394-395`

```python
if job_settings is not None:
    try:
        qmmm_merged = ORCAQMMMJobSettings(
            **getattr(job_settings, "__dict__", job_settings)
        )
    except Exception:
        qmmm_merged = qmmm_settings
else:
    qmmm_merged = qmmm_settings          # <-- line 296-297

...
if parent_jobtype is not None:            # <-- line 317-318
    qmmm_settings.parent_jobtype = parent_jobtype
...
if parent_settings is not None:           # <-- line 394-395
    inherited_keywords = [...]
    try:
        qmmm_settings = qmmm_settings.merge(parent_settings, keywords=inherited_keywords)
    except Exception as exc:
        ...
```

Two independent invariants make these guards' False arms unreachable:

1. `job_settings = ctx.obj["job_settings"]` (line 273) is always populated by the parent `orca` group's callback (`chemsmart/cli/orca/orca.py`) before any subcommand -- including `qmmm` -- ever runs; every code path in that callback assigns a real `ORCAJobSettings` instance to `job_settings` (default, from-file, or from-database), never `None`. So the `else: qmmm_merged = qmmm_settings` branch (line 296-297, guarded by `job_settings is not None`'s False arm) can never execute.
2. `create_orca_qmmm_subcommand` is attached to exactly seven parent groups: `opt`, `ts`, `sp` (singlepoint), `scan`, `qrc`, `modred`, and `neb` (`grep -rn "create_orca_qmmm_subcommand" chemsmart/cli/orca/*.py`). Every one of these seven groups' own callback sets both `ctx.obj["parent_settings"]` and `ctx.obj["parent_jobtype"]` to a real (non-`None`) value immediately before Click dispatches to the `qmmm` subcommand (e.g. `chemsmart/cli/orca/opt.py:95-96`: `ctx.obj["parent_settings"] = opt_settings; ctx.obj["parent_jobtype"] = "opt"`). Since `qmmm` has no other attachment point, `parent_jobtype` and `parent_settings` are always non-`None` whenever the `qmmm()` callback runs through the real CLI, making both `if parent_jobtype is not None:` (317-318) and `if parent_settings is not None:` (394-395) always take their True arm.

**Reproduce:** `tests/test_orca_qmmm_cli.py::TestOrcaQmmmSubcommand` exercises the *reachable* neighboring branches (the `qmmm_merged`/`qmmm_settings` exception-fallback paths when `job_settings is not None`, via `test_qmmm_merge_failure_falls_back_to_settings_reconstruction` and related tests); confirmed unreachability for the three branches above by exhaustively checking all seven `create_orca_qmmm_subcommand` call sites and the `orca` group's `job_settings` assignment paths.

**Impact:** None -- purely redundant defensive code with no behavioral effect through any real CLI invocation.

**Suggested direction:** no action needed; could be simplified by removing the `job_settings is not None` / `parent_jobtype is not None` / `parent_settings is not None` guards (trusting the invariants above), or left as-is as cheap insurance against a future eighth attachment point that might not set `parent_settings`/`parent_jobtype`.

---

## 79. `cli/orca/qmmm.py`'s `-h`/`--high-level-h-bond-length` option can never successfully reach `ast.literal_eval` (extends bug #28)

**Location:** `chemsmart/cli/orca/qmmm.py:438-442`

Bug #28 already documents that Click's `type=dict` declaration for `-h` rejects any non-empty CLI value before the callback body ever runs. This entry adds the missing piece: the *one* value that does survive Click's `dict(value)` conversion -- the empty string `""`, which becomes `{}` -- still can't reach a successful `ast.literal_eval` call. `high_level_h_bond_length` ends up as the dict `{}` (not a string), and `ast.literal_eval({})` raises `ValueError: malformed node or string: {}` (confirmed via `tests/test_orca_qmmm_cli.py::TestOrcaQmmmSubcommand::test_high_level_h_bond_length_empty_string_crashes_differently`). So line 442 (`molecule.scale_factors = high_level_h_bond_length_dict`, reached only after a *successful* `ast.literal_eval`) is unreachable for literally any value `-h` can take through the real CLI -- not just the "typical" values bug #28 already covers, but even the one edge case that gets past Click's broken converter.

**Reproduce:** `tests/test_orca_qmmm_cli.py::TestOrcaQmmmSubcommand::test_high_level_h_bond_length_empty_string_crashes_differently` -- `-h ""` reaches `chemsmart/cli/orca/qmmm.py:439` (`ast.literal_eval(high_level_h_bond_length)`) and raises `ValueError: malformed node or string: {}` rather than completing.

**Impact:** Same as bug #28 -- the `-h` option is completely unusable for its documented purpose through any real CLI input, not just typical ones.

**Suggested direction:** same as bug #28: fix the option's `type=` declaration (e.g. `type=str`, parsed with `ast.literal_eval` in the callback as clearly intended) so a real dict-literal string like `"{1: 1.1}"` can actually reach line 439-442 successfully.

---

## 80. `ORCAInputWriter._write_scf_convergence`'s error message references the wrong object, crashing with `AttributeError` instead of raising the intended `ValueError`

**Location:** `chemsmart/jobs/orca/writer.py:232-244`

```python
def _write_scf_convergence(self, f):
    if self.settings.scf_convergence:
        from chemsmart.io.orca import ORCA_SCF_CONVERGENCE

        scf_conv = self.settings.scf_convergence.lower().strip()
        if scf_conv.endswith("scf"):
            scf_conv = scf_conv[:-3]
            if scf_conv not in ORCA_SCF_CONVERGENCE:
                raise ValueError(
                    f"Warning: SCF convergence {self.scf_convergence} is not supported by ORCA!\n"
                    f"Available SCF convergence options are: {ORCA_SCF_CONVERGENCE}"
                )
        f.write(f"  convergence {scf_conv}\n")
```

The f-string building the `ValueError` message references `self.scf_convergence` -- an attribute of the *writer* (`ORCAInputWriter`), which has no such attribute (the writer only exposes `self.settings.scf_convergence`, used correctly two lines above). Evaluating `self.scf_convergence` inside the f-string raises `AttributeError: 'ORCAInputWriter' object has no attribute 'scf_convergence'` before the `ValueError` can even be constructed, so callers get an unrelated, confusing crash instead of the intended "SCF convergence X is not supported by ORCA" message.

**Reproduce:** `tests/test_orca_writer_blocks_unit.py::TestScfBlock::test_invalid_convergence_raises` -- calling `_write_scf_convergence` with `scf_convergence="bogusSCF"` (a value with a valid `...scf` suffix but not in `ORCA_SCF_CONVERGENCE` once stripped) raises `AttributeError`, not `ValueError`.

**Impact:** Low-to-moderate -- functionally the input is still rejected (the job still fails to write), but the error message is actively misleading: a user who mistypes an SCF convergence keyword sees an unrelated `AttributeError` about a missing writer attribute instead of the helpful "not supported by ORCA" message with the list of valid options.

**Suggested direction:** change `self.scf_convergence` to `self.settings.scf_convergence` in the f-string.

---

## 81. Three `elif`/`else` branches in `jobs/orca/writer.py` are unreachable dead code given their guards' fixed value domains

**Location:** `chemsmart/jobs/orca/writer.py:434-457, 460-476, 602-613`

Three separate `if`/`elif` chains in this file check a value against a fixed, small set of options, where an earlier `assert` (or a helper function with a documented, exhaustive return domain) already guarantees the value can only be one of the exact strings each `elif` branch checks for -- making each chain's *final* `elif`'s "condition true" arc's implicit fallthrough (i.e. the case where none of the branches match) unreachable:

1. `_write_mdci_block`'s `mdci_cutoff` chain (`if ...== "loose": elif ...== "normal": elif ...== "tight":`) is preceded by `assert mdci_cutoff.lower() in ["loose", "normal", "tight"]`, so by the time the `elif` chain runs, the value is guaranteed to match exactly one of the three branches -- the "tight" branch is only ever reached when the value truly is "tight", so its own condition is always true when evaluated.
2. Similarly, `_write_mdci_block`'s nested `mdci_density` chain is preceded by `assert mdci_density.lower() in ["none", "unrelaxed", "relaxed"]`.
3. `_write_modred_if_dict`'s `prepend_string.lower().startswith(...)` chain (`"b"`/`"a"`/`"d"`, else `scan_var = "variable"`) processes `prepend_string` values produced exclusively by `get_prepend_string_for_modred` (`chemsmart/utils/utils.py:1221-1243`), which raises `ValueError` for any coordinate list not of length 2, 3, or 4 and otherwise returns exactly `"B"`, `"A"`, or `"D"` -- so `prepend_string.lower()` can only ever be `"b"`, `"a"`, or `"d"`, and the `else: scan_var = "variable"; scan_unit = "unit"` fallback can never execute.

**Reproduce:** confirmed by reading the guarding `assert`s (cases 1-2) and `get_prepend_string_for_modred`'s exhaustive 3-value return domain (case 3); `tests/test_orca_writer_blocks_unit.py::TestMdciBlock::test_cutoff_levels`/`test_density_levels` and `TestModredBlock` already cover the reachable branches of all three chains.

**Impact:** None -- purely redundant code paths with no behavioral effect (each `elif`'s condition is always true whenever control reaches it, and the `else` case in #3 can never be entered at all).

**Suggested direction:** no action needed; harmless as defensive coding, though the final `elif` in each chain could be simplified to a plain `else` now that exhaustiveness is confirmed.

---

## 82. `Thermochemistry.qrrho_vibrational_entropy` silently returns `0.0` instead of raising when `entropy_method` isn't `"grimme"` or `"truhlar"`

**Location:** `chemsmart/analysis/thermochemistry.py:901-940`

```python
@property
def qrrho_vibrational_entropy(self):
    ...
    if self.s_freq_cutoff is None or self.v is None:
        return None
    vib_entropy = []
    if self.entropy_method == "grimme":
        ...
        for j in range(0, len(self.v)):
            vib_entropy.append(...)
    elif self.entropy_method == "truhlar":
        ...
        for j in range(0, len(self.v)):
            vib_entropy.append(...)
    return sum(vib_entropy)
```

`entropy_method` defaults to `None` in both `Thermochemistry.__init__` (`chemsmart/analysis/thermochemistry.py:73`) and `ThermochemistryJobSettings.__init__` (`chemsmart/jobs/thermochemistry/settings.py:48`), and neither class validates that `entropy_method` is set (or is one of the two supported values) whenever `s_freq_cutoff` is provided. If a caller sets `s_freq_cutoff` without also setting `entropy_method` to `"grimme"` or `"truhlar"` -- e.g. by constructing `Thermochemistry` or `ThermochemistryJobSettings` directly rather than going through the CLI's `resolve_entropy_cutoff` helper (`chemsmart/cli/thermochemistry/thermochemistry.py:62-76`), which always pairs a supplied cutoff with the matching method -- then neither `if` nor `elif` branch is taken, `vib_entropy` stays `[]`, and the property returns `sum([])`, i.e. `0.0`. This silently propagates: `qrrho_total_entropy` becomes translational + rotational + electronic entropy only (vibrational contribution dropped to zero), and `qrrho_gibbs_free_energy*` variants become correspondingly wrong -- with no error, warning, or `None` to signal that anything is off. A typo in `entropy_method` (e.g. `"Grimme"` with capital G, since the comparison is case-sensitive) hits the same silent-zero path.

**Reproduce:** `tests/test_thermochemistry.py::TestThermochemistryRemainingBranches::test_qrrho_vibrational_entropy_unknown_method_returns_zero` constructs a mock with `s_freq_cutoff=100.0`, real `v`, and `entropy_method=None` and shows `qrrho_vibrational_entropy` returns `0` rather than raising.

**Impact:** Low-to-moderate -- unreachable through the documented CLI entry points (which always pair `s_freq_cutoff` with a valid `entropy_method`), but reachable via any direct/programmatic use of `Thermochemistry`, `ThermochemistryJobSettings`, or `BoltzmannAverageThermochemistry` (e.g. scripting against the library, or a future caller that forgets the pairing). The failure mode is silent and physically wrong (vibrational entropy dropped entirely) rather than a loud error, which is the worse kind of bug to have in a thermochemistry engine.

**Suggested direction:** raise a `ValueError` (e.g. in `__init__`, mirroring the existing `energy_type` validation in `BoltzmannAverageThermochemistry`) when `s_freq_cutoff` is set but `entropy_method` is not `"grimme"` or `"truhlar"`, instead of letting `qrrho_vibrational_entropy` fall through to an empty list.

---

## 83. `BoltzmannAverageThermochemistry.__init__`'s empty-file-list guard is unreachable -- an empty list raises `IndexError`, not the documented `ValueError`

**Location:** `chemsmart/analysis/thermochemistry.py:1535-1554`

```python
def __init__(self, files, energy_type="gibbs", **kwargs):
    super().__init__(
        filename=files[
            0
        ],  # No single file, we will take molecule from first filename
        **kwargs,
    )
    """
    Initialize with a list of Gaussian or ORCA output files.
    ...
    """
    if not files:
        raise ValueError("List of files cannot be empty.")
```

The `if not files: raise ValueError(...)` guard is meant to give a clear error when `files=[]` is passed. However, `files[0]` is evaluated as part of the `super().__init__(...)` call *before* that guard is ever reached -- so an empty `files` list raises a bare `IndexError: list index out of range` from the `files[0]` subscript, and the guard clause a few lines later is dead code that can never execute for this input. (The docstring placed between the `super().__init__()` call and the guard is itself also somewhat telling -- it reads like the guard was originally intended to run first, ahead of `super().__init__()`, and got left behind after a reorder.)

**Reproduce:** `tests/test_thermochemistry.py::TestBoltzmannAverageThermochemistryRemainingBranches::test_empty_files_raises_indexerror_before_reaching_guard` confirms `BoltzmannAverageThermochemistry(files=[], temperature=298.15)` raises `IndexError`, not `ValueError`. This is also the sole remaining coverage gap in `chemsmart/analysis/thermochemistry.py` (line 1554) after this effort's test additions -- the `ValueError` branch is provably unreachable as the code is currently ordered.

**Impact:** Low -- the call still fails loudly (just with the wrong exception type/message), so callers checking for empty input via a broad `except Exception` won't notice, but ones specifically catching `ValueError` (as the docstring/API implies they should be able to) will see an uncaught `IndexError` instead.

**Suggested direction:** move the `if not files: raise ValueError(...)` check to the very top of `__init__`, before the `super().__init__(filename=files[0], **kwargs)` call.

---

## 84. Two identical `if not route_string.startswith("!")` guards are dead code because `route_string` is always freshly initialized to `""`

**Location:** `chemsmart/jobs/orca/settings.py:526-528` (`ORCAJobSettings._get_route_string_from_jobtype`) and `chemsmart/jobs/orca/settings.py:2619-2621` (`ORCANEBJobSettings._get_neb_route_string`)

```python
route_string = ""
if not route_string.startswith("!"):
    route_string += "! "
```

In both methods, `route_string` is assigned the literal `""` on the line immediately before this guard, and no code path modifies it in between. An empty string never starts with `"!"`, so `not route_string.startswith("!")` is always `True` and the body always executes -- the `if` is unconditionally true every time either method runs. The guard can never take its "false" arc (i.e. `route_string` already starting with `"!"`) because there is no way to reach the check with a non-empty `route_string`.

**Reproduce:** Coverage of `tests/test_orca_settings_unit.py` and `tests/test_orca_qmmm_neb_settings_unit.py` (all `route_string`/`_get_neb_route_string` tests) leaves the `if`'s false branch (`527->532` and `2620->2625` in `coverage report -m`) permanently unreached no matter what settings are exercised, confirming the guard is tautological.

**Impact:** None -- purely dead/defensive code; every call always takes the same path, so behavior is unaffected. It looks like a copy-paste of the pattern used in `_get_route_string_from_user_input` (where the check *is* meaningful, since there `route_string = self.route_to_be_written` can legitimately already start with `"!"`), applied to a context where it can't ever be false.

**Suggested direction:** no action needed; could be simplified to unconditionally `route_string = "! "` in both places, removing the always-true `if`.

---

## 85. `ORCAJobSettings._get_level_of_theory`'s `elif self.basis is None` is a tautological branch that is always true when reached

**Location:** `chemsmart/jobs/orca/settings.py:665-671`

```python
if self.basis is not None:
    level_of_theory += f" {self.basis}"
elif self.basis is None:
    # allow missing basis for QMMM-type jobs where basis may be
    # provided per-layer (or omitted)
    if not is_qmmm:
        raise ValueError("Warning: basis is missing!")
```

The `elif` branch is only ever evaluated when the preceding `if self.basis is not None` was `False`, i.e. exactly when `self.basis is None`. At that point `elif self.basis is None` is guaranteed to be `True` -- it can never be reached and evaluate to `False`, since that would require `self.basis` to be simultaneously not-None (to fail the `if`) and None (to fail the `elif`), which is impossible. The `elif`'s condition is therefore redundant; it is functionally an `else`.

**Reproduce:** In `coverage report -m --include="*jobs/orca/settings.py"`, branch `667->673` (jumping from the `elif` line straight past its body to the `aux_basis` check that follows the whole `if/elif`) never appears as covered no matter what `basis`/`jobtype` combinations `tests/test_orca_settings_unit.py::TestGetLevelOfTheory` exercises, because that arc requires the impossible "basis is not None and also None" state.

**Impact:** None -- purely a redundant/tautological condition; behavior is identical to writing a plain `else:`.

**Suggested direction:** no action needed; could be simplified by replacing `elif self.basis is None:` with `else:` for clarity, since the condition adds no information.

---

## 86. `ORCAQMMMJobSettings._get_level_of_theory_string`'s `if self.low_level_of_theory is not None` check is always true because `validate_and_assign_level` unconditionally returns `"MM"` for the low-level layer

**Location:** `chemsmart/jobs/orca/settings.py:2015-2016, 2091-2094`

```python
# inside validate_and_assign_level(self, functional, basis, built_in_method, level_name):
if level_name == "low_level":
    level_of_theory = "MM"
return level_of_theory

# inside _get_level_of_theory_string:
self.low_level_of_theory = self.validate_and_assign_level(
    None, None, self.low_level_method, level_name="low_level"
)
if self.low_level_of_theory is not None:
    ...
```

`validate_and_assign_level` is always called with `level_name="low_level"` from this call site, and the `if level_name == "low_level": level_of_theory = "MM"` line unconditionally overwrites whatever the earlier `built_in_method`/`functional`/`basis` branches computed (including the `""` default), regardless of whether `self.low_level_method` was actually provided. So `self.low_level_of_theory` is always the string `"MM"`, never `None` or falsy, and the subsequent `if self.low_level_of_theory is not None:` check on line 2094 can never take its false branch.

**Reproduce:** `coverage report -m --include="*jobs/orca/settings.py"` shows arc `2094->2106` (skipping the entire QM/QMMM-labeling block) as unreachable across all of `tests/test_orca_qmmm_neb_settings_unit.py::TestQMMMLevelOfTheoryString`, including jobs with no `low_level_method` set at all (plain `QM/QM2` jobs).

**Impact:** None observed -- the *actual* decision about whether to append an `"/MM"` suffix is separately (and correctly) gated a few lines later by `if self.low_level_method is not None:` (line 2104), which checks the real attribute rather than the always-`"MM"` `low_level_of_theory`. So a plain `QM/QM2` job (no MM layer) still produces the correct `"QM/<intermediate>"` route string without an erroneous `"/MM"` suffix, because that second, correct guard saves it. The outer `if self.low_level_of_theory is not None:` is simply vestigial dead code that happens not to matter because of the redundant, correctly-guarded check inside it.

**Suggested direction:** no action needed functionally, but this is worth a closer look if `_get_level_of_theory_string` is ever refactored, since the outer condition reads as if it's meaningfully distinguishing "has an MM layer" from "doesn't" when it cannot -- a future edit that removes the inner `low_level_method` guard (trusting the outer one instead) would silently break additive-QMMM-vs-QM/QM2 formatting.

---

## 87. `ORCAQMMMJobSettings._write_qmmm_block`'s `if crystal_sub is not None` check is always true because `_write_crystal_qmmm_subblock` never returns `None`

**Location:** `chemsmart/jobs/orca/settings.py:2427-2431, 2445-2471`

```python
crystal_sub = self._write_crystal_qmmm_subblock()
if crystal_sub is not None:
    full_qm_block += crystal_sub

...

def _write_crystal_qmmm_subblock(self):
    crystal_qmmm_subblock = ""
    if not self.conv_charges:
        ...
    ...
    return crystal_qmmm_subblock
```

`_write_crystal_qmmm_subblock` initializes `crystal_qmmm_subblock = ""` and only ever appends to it; every code path returns this (possibly-still-empty) string, never `None`. So `crystal_sub` is always a `str` (falsy-but-not-None when no crystal fields are set), and `if crystal_sub is not None:` is always `True`. Appending an empty string to `full_qm_block` is a no-op anyway, so the check has no observable effect either way.

**Reproduce:** `coverage report -m --include="*jobs/orca/settings.py"` shows arc `2428->2431` (skipping the `full_qm_block += crystal_sub` append) as unreachable across all of `tests/test_orca_qmmm_neb_settings_unit.py::TestQmmmBlockGeneration`, including non-crystal `QMMM`/`QM/QM2` jobs where the subblock is empty.

**Impact:** None -- appending `""` is a no-op, so whether the `if` guard is "true but appends nothing" or hypothetically "false and skips the append" produces byte-identical output. Purely dead/defensive code.

**Suggested direction:** no action needed; could be simplified to an unconditional `full_qm_block += self._write_crystal_qmmm_subblock()`, since the `is not None` guard can never be the deciding factor.

---

## 88. `GaussianQRCJob._prepare_both_qrc_jobs`'s `elif direction == "r":` is a tautological branch that is always true when reached

**Location:** `chemsmart/jobs/gaussian/qrc.py:126-135`

```python
for direction in ["f", "r"]:
    ...
    if direction == "f":
        mol = self.qrcf_molecule
    elif direction == "r":
        mol = self.qrcr_molecule
```

`direction` only ever takes the two values in the loop's own fixed iterable, `["f", "r"]`. Whenever the `if direction == "f":` branch is `False`, `direction` is guaranteed to be `"r"` (the only other value the loop produces), so `elif direction == "r":` is always `True` when evaluated -- its "false" arm (which would leave `mol` as the unmodified `self.molecule` from line 131) can never execute.

**Reproduce:** `tests/test_GaussianJobs.py::TestGaussianQRCJobs` exercises both loop iterations (via `test_run_both_jobs_runs_forward_and_reverse_jobs` and the new label/jobtype tests); `coverage report -m --include="*jobs/gaussian/qrc.py"` shows arc `134->136` (the `elif`'s false arm, skipping to the `logger.debug` call after the chain) as permanently unreachable regardless of which branch combinations are exercised.

**Impact:** None -- purely redundant code with no behavioral effect, since the loop's iterable is a fixed two-element list.

**Suggested direction:** no action needed; could be simplified to a plain `else:` now that the two-value domain is confirmed.

**Suggested direction:** no action needed; if desired, the loop could be simplified to a direct index into `self.molecular_dipole_lines[0]` instead of a `for`/`startswith` scan, now that the ordering invariant is confirmed.

---

## 89. `Gaussian16Output.hirshfeld_charges`/`hirshfeld_spin_densities`/`hirshfeld_dipoles`/`hirshfeld_cm5_charges` crash with `IndexError` instead of returning `None` when the file has no Hirshfeld section

**Location:** `chemsmart/io/gaussian/output.py:2340-2383`

```python
def _get_hirshfeld_charges_spins_dipoles_cm5(self):
    all_hirshfeld_charges = []
    all_spin_densities = []
    all_dipoles = []
    all_cm5_charges = []
    for i, line_i in enumerate(self.contents):
        ...
        if (
            "Hirshfeld charges, spin densities, dipoles, and CM5 charges"
            in line_i
        ):
            ...
            all_hirshfeld_charges.append(hirshfeld_charges)
            all_spin_densities.append(spin_densities)
            all_dipoles.append(dipoles)
            all_cm5_charges.append(cm5_charges)
    return (
        all_hirshfeld_charges[-1],
        all_spin_densities[-1],
        all_dipoles[-1],
        all_cm5_charges[-1],
    )
```

Unlike its heavy-atom counterpart `_get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms` (a few lines below, at `output.py:2385-2451`), which explicitly checks `if all_hirshfeld_charges_heavy_atoms and ...` / `elif ... :` / `else: return None, None, None` before indexing, this plain (non-heavy-atom) variant has no such guard: it unconditionally indexes `all_hirshfeld_charges[-1]` etc. If the "Hirshfeld charges, spin densities, dipoles, and CM5 charges" header is never found in the file (e.g. a calculation that doesn't request Hirshfeld population analysis), all four accumulator lists stay empty, and `all_hirshfeld_charges[-1]` raises `IndexError: list index out of range` instead of gracefully falling back to `None`, `None`, `None`, `None` the way the heavy-atom sibling does for the exact same "not found" case.

**Reproduce:** `tests/test_GaussianIO.py::TestGaussian16OutputAdditionalCoverage::test_hirshfeld_charges_raises_indexerror_when_section_absent` -- for `nhc_neutral_singlet.log` (a normal opt+freq file with no Hirshfeld analysis requested), `g16.hirshfeld_charges` raises `IndexError`, while `g16.hirshfeld_charges_heavy_atoms` (and its spin/CM5 counterparts) correctly return `None` for the identical "section absent" situation.

**Impact:** Moderate. Any caller that accesses `hirshfeld_charges`, `hirshfeld_spin_densities`, `hirshfeld_dipoles`, or `hirshfeld_cm5_charges` on a Gaussian output file without a Hirshfeld analysis section (the common case, since Hirshfeld population analysis must be explicitly requested) gets an uncaught `IndexError` crash instead of a `None` they could check for, unlike every other "not found" property in this class (including the heavy-atom sibling of this exact data).

**Suggested direction:** guard the final `return` the same way the heavy-atom variant does, e.g. `return (all_hirshfeld_charges[-1] if all_hirshfeld_charges else None, ...)` for each of the four values, or restructure with an explicit `if not all_hirshfeld_charges: return None, None, None, None` before the indexed return.

---

## 90. `Gaussian16Output.hirshfeld_cm5_charges_heavy_atoms` returns a one-element list-of-dicts instead of a dict when spin densities are also present

**Location:** `chemsmart/io/gaussian/output.py:2432-2440`

```python
if (
    all_hirshfeld_charges_heavy_atoms
    and all_hirshfeld_spin_densities_heavy_atoms
):
    return (
        all_hirshfeld_charges_heavy_atoms[-1],
        all_hirshfeld_spin_densities_heavy_atoms[-1],
        all_cm5_charges_heavy_atoms,          # <-- missing [-1]
    )
elif (
    all_hirshfeld_charges_heavy_atoms
    and not all_hirshfeld_spin_densities_heavy_atoms
):
    return (
        all_hirshfeld_charges_heavy_atoms[-1],
        None,
        all_cm5_charges_heavy_atoms[-1],      # <-- correctly indexed here
    )
```

`_get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms` accumulates one dict per matched section into `all_cm5_charges_heavy_atoms` (a list), exactly like the other two accumulators. Both the charges and spin-densities return values are correctly indexed with `[-1]` (the last/most recent occurrence) in *both* branches of the `if`/`elif`. The CM5 return value, however, only gets `[-1]` in the `elif` (closed-shell, no spin densities) branch; in the `if` branch (open-shell, spin densities present), it returns the raw list `all_cm5_charges_heavy_atoms` unindexed. `hirshfeld_cm5_charges_heavy_atoms` therefore returns a `list` containing one `dict` for any open-shell Hirshfeld calculation, instead of the `dict` itself that every other branch (and every sibling property) returns.

**Reproduce:** `tests/test_GaussianIO.py::TestGaussian16OutputAdditionalCoverage::test_hirshfeld_cm5_charges_heavy_atoms_wrong_type_when_spin_present` -- for `oxetane_rc_hirshfeld_sp_smd_n_n-DiMethylFormamide.log` (an open-shell radical-cation Hirshfeld calculation, where `hirshfeld_spin_densities_heavy_atoms` is non-`None`), `g16.hirshfeld_cm5_charges_heavy_atoms` returns `[{'O1': 0.012879, ...}]` (a one-element `list`) rather than `{'O1': 0.012879, ...}` (a `dict`), even though the closed-shell fixture `oxetane_hirshfeld_sp_smd_n_n-DiMethylFormamide.log` correctly returns a plain `dict` for the same property.

**Impact:** Moderate. Any caller doing `hirshfeld_cm5_charges_heavy_atoms["O1"]` (the natural/documented usage, matching every other charges dict in this class) raises `TypeError: list indices must be integers or slices, not str` for any open-shell Hirshfeld calculation, while working fine for closed-shell ones -- a silent, data-dependent type inconsistency.

**Suggested direction:** change `all_cm5_charges_heavy_atoms` to `all_cm5_charges_heavy_atoms[-1]` in the `if` branch's return tuple, matching the `elif` branch immediately below it.

---

## 91. `Gaussian16Output.moments_of_inertia`/`moments_of_inertia_principal_axes` crash with `TypeError` instead of returning `None` when the "Principal axes" section is absent

**Location:** `chemsmart/io/gaussian/output.py:2475-2530`

```python
@cached_property
def moments_of_inertia(self):
    moments_of_inertia, _ = (
        self._get_moments_of_inertia_and_principal_axes()
    )
    return moments_of_inertia

def _get_moments_of_inertia_and_principal_axes(self):
    for i, line in enumerate(self.contents):
        if "Principal axes and moments of inertia" in line:
            ...
            return np.array(moments_of_inertia), np.array(
                moments_of_inertia_principal_axes
            )
    # falls off the end here -> implicitly returns a single `None`,
    # not a `(None, None)` tuple
```

`_get_moments_of_inertia_and_principal_axes` only ever `return`s explicitly inside the `if "Principal axes..." in line:` branch, as a 2-tuple. If that line is never found in the file, the `for` loop exhausts and the function falls off the end, implicitly returning a single `None` (Python's default), not `(None, None)`. Both `moments_of_inertia` and `moments_of_inertia_principal_axes` unconditionally unpack the helper's return value as a 2-tuple (`moments_of_inertia, _ = ...` / `_, principal_axes = ...`), so `None` cannot be unpacked and both properties raise `TypeError: cannot unpack non-iterable NoneType object` instead of gracefully returning `None`, unlike essentially every other "section not found" property in this class (e.g. `zero_point_energy`, `rotational_symmetry_number`, `pbc` in the PBC subclass, etc., all of which have an explicit trailing `return None`).

**Reproduce:** `tests/test_GaussianIO.py::TestGaussian16OutputAdditionalCoverage::test_moments_of_inertia_crashes_when_section_absent` -- for `oxygen_openshell_singlet_sp_link.log` (a real single-point link-job output, which never prints the "Principal axes and moments of inertia" banner), both `g16.moments_of_inertia` and `g16.moments_of_inertia_principal_axes` raise `TypeError` rather than returning `None`.

**Impact:** Moderate. Any single-point (or otherwise banner-less) Gaussian output crashes on `moments_of_inertia`/`moments_of_inertia_principal_axes` access instead of returning `None`, which is surprising given the consistent "return `None` when not found" convention used everywhere else in this class.

**Suggested direction:** add an explicit `return None, None` (or `return None`) after the `for` loop in `_get_moments_of_inertia_and_principal_axes`.

---

## 92. Three separate unreachable defensive guards in `Gaussian16Output`'s gen/genecp and ECP parsing helpers

**Location:** `chemsmart/io/gaussian/output.py:173-180, 303-306, 382-388`

```python
# (a) _genecp_info, line 173-180
if self.gen_genecp is None:
    return result
try:
    atom_symbols = self.symbols
except Exception:
    return result
if not atom_symbols:          # <-- unreachable
    return result

# (b) _parse_pseudopotential_section, line 303-306
def _parse_pseudopotential_section(self):
    try:
        atom_symbols = self.symbols     # <-- self.symbols already
    except Exception:                   #     succeeded and is cached
        return {}                       #     by the only caller

# (c) _parse_pseudopotential_section, line 382-388
elif (
    current is not None
    and not center_re.match(line)      # <-- redundant with `is_center`
    and not term_re.match(line)        #     and `is_term` above, which
):                                      #     already consume any line
    _flush_channel()                   #     these regexes would match
    current["_channel_name"] = line
```

Three independent, unreachable branches, all discovered while writing coverage for the gen/genecp and ECP (pseudopotential) parsing helpers:

(a) `self.symbols` (from `GaussianFileMixin`, via `CoordinateBlock.chemical_symbols` -> `_get_symbols`) can never return an empty/falsy list: `_get_symbols` explicitly `raise`s `ValueError("No symbols found in the coordinate block...")` whenever its accumulated `symbols` list would be empty. So `self.symbols` either returns a non-empty list or raises -- it can never reach `_genecp_info`'s `if not atom_symbols:` check with a falsy value; that check is dead.

(b) `_parse_pseudopotential_section` is only ever called from one call site, `_genecp_info` (`output.py:253`), and only *after* `_genecp_info` has already evaluated `atom_symbols = self.symbols` successfully at its own top (line 176) -- if that raised, `_genecp_info` would already have returned early via its own `except Exception: return result` at line 177-178, never reaching the call to `_parse_pseudopotential_section` at all. Since `symbols` is a `cached_property`, the second access inside `_parse_pseudopotential_section` is guaranteed to hit the cache and return the same already-successful value. Its own `try/except` around `self.symbols` can therefore never catch anything.

(c) Within `_parse_pseudopotential_section`'s per-line dispatch, `is_center` (`len(tokens) in (2, 3) and all(t.isdigit() for t in tokens)`) and `center_re` (`r"^\d+(?:\s+\d+){1,2}\s*$"`) match the same set of strings, as do `is_term` and `term_re`. Any line for which `center_re.match(line)` or `term_re.match(line)` would be `True` is therefore already consumed earlier by the `if is_center:` / `elif is_term and current is not None:` branches above, so by the time control reaches the final `elif ... not center_re.match(line) and not term_re.match(line):`, those two sub-conditions are always `True` given `current is not None`. The only way this `elif`'s condition could be `False` is `current is None`, but that can't happen either for genuine Gaussian ECP output, since the very first data line after the second `===` separator is always a center header (matching `is_center`), which sets `current` before any other line type is seen.

**Reproduce:** covered indirectly by `tests/test_GaussianIO.py::TestGaussian16OutputAdditionalCoverage::test_genecp_info_empty_symbols_raises_and_is_caught` (demonstrates (a)/(b): `self.symbols` raises rather than returning empty) and `test_parse_pseudopotential_section_only_reachable_via_genecp_info` (demonstrates (b): direct call succeeds using the already-cached `symbols`). Branch (c) remains uncovered in the coverage report at `chemsmart/io/gaussian/output.py:382->339`; no fixture could be constructed to reach it without fabricating non-Gaussian-shaped ECP table input.

**Impact:** None -- pure dead code / redundant defensive programming with no observed effect on behavior for any real Gaussian output.

**Suggested direction:** for (a), drop the `if not atom_symbols: return result` check (or change `_get_symbols` to return `[]` instead of raising, if that's ever desired, and keep the check). For (b), drop the redundant inner `try/except` (or, if defensive-programming-by-habit is preferred, leave as documentation that the method assumes a pre-validated `self.symbols`). For (c), simplify the final `elif` to just `elif current is not None:` since the regex re-checks add nothing given the earlier `is_center`/`is_term` dispatch already filters those line shapes out.

---

## 93. Six vibrational-property "stop at Thermochemistry" early-exit checks are unreachable dead code

**Location:** `chemsmart/io/gaussian/output.py:1068-1077, 1086-1095, 1104-1113, 1122-1131, 1140-1150, 1161-1184` (`vibrational_frequencies`, `reduced_masses`, `force_constants`, `ir_intensities`, `vibrational_mode_symmetries`, `vibrational_modes`)

```python
@cached_property
def vibrational_frequencies(self):
    frequencies = []
    for line in self.contents:
        if line.startswith("Frequencies --"):
            freq_string = line.split("--")[1].strip()
            for freq in freq_string.split():
                frequencies.append(float(freq))
        else:
            continue
        if "Thermochemistry" in line:
            break
    return frequencies
```

All six of these properties share the identical `if <prefix match>: ... else: continue` / `if "Thermochemistry" in line: break` idiom, apparently intended as an early-exit once the frequency block has been fully consumed and the "-------------------- Thermochemistry --------------------" banner that follows it is reached. But because the `else` branch does `continue` (skipping straight to the next loop iteration), the `if "Thermochemistry" in line:` check is *only* ever reached immediately after the `if` branch already ran -- i.e. only for a line that just matched `line.startswith("Frequencies --")` (or the analogous `"Red. masses --"`, `"Frc consts  --"`, `"IR Inten    --"` prefix, or was consumed as a normal-mode data row). No such line can simultaneously also contain the substring `"Thermochemistry"` in Gaussian's fixed output format -- the banner line is its own separate line, never combined with a frequency/mass/constant/intensity data line. So the `break` can never fire for real Gaussian output; the loop always runs to completion over the full file instead of stopping early once past the last frequency block.

**Reproduce:** confirmed via `coverage report` on `chemsmart/io/gaussian/output.py` -- the six `break` statements at lines 1076, 1094, 1112, 1130, 1149, and 1184 remain in the missing-lines list even after `tests/test_GaussianIO.py` exercises dozens of real frequency-containing fixtures (`gaussian_singlet_opt_outfile`, `gaussian_triplet_opt_outfile`, `gaussian_ts_genecp_outfile`, etc.) through these exact properties.

**Impact:** Low -- purely a missed micro-optimization (the loop scans the rest of the file instead of stopping early), not a correctness bug; the collected frequencies/masses/constants/intensities/symmetries/modes are unaffected either way since nothing after the frequency block would match the `startswith` prefixes anyway.

**Suggested direction:** if the early-exit is still wanted, move the `"Thermochemistry" in line` check so it triggers on any line (not gated behind the `else: continue`), e.g. check it unconditionally at the top of the loop body before the `if line.startswith(...)`. Otherwise, remove the dead check entirely.

---

## 94. `frozen_coordinate_indices` and `_get_frozen_and_free_atoms`'s `if len(line_i) == 0: break` checks are unreachable

**Location:** `chemsmart/io/gaussian/output.py:1239-1242, 1285-1288`

```python
for i, line_i in enumerate(self.contents):
    if "Symbolic Z-matrix:" in line_i:
        if len(line_i) == 0:
            break
        for j, line_j in enumerate(self.contents[i + 2 :]):
            ...
```

Both `frozen_coordinate_indices` and `_get_frozen_and_free_atoms` guard their inner parsing loop with `if "Symbolic Z-matrix:" in line_i:` and then immediately check `if len(line_i) == 0: break` before doing anything else. But the outer condition `"Symbolic Z-matrix:" in line_i` already guarantees `line_i` contains that 19-character substring, so `len(line_i)` is necessarily >= 19 whenever the inner check runs -- `len(line_i) == 0` can never be `True` at that point. The `break` is unreachable.

**Reproduce:** confirmed via `coverage report` -- lines 1242 and 1288 remain in the missing-lines list even with `tests/test_GaussianIO.py::TestGaussian16Output::test_read_frozen_opt_outputfile` (which exercises both properties on a real frozen-coordinate fixture) passing.

**Impact:** None -- dead code with no effect on behavior; likely copy-pasted from a similar (also arguably unnecessary) pattern elsewhere in the file.

**Suggested direction:** remove both unreachable `if len(line_i) == 0: break` checks.

---

## 95. `_get_route`'s final `else: route = None` and `energies`'s implicit fallthrough are both unreachable given their preceding exhaustive conditions

**Location:** `chemsmart/io/gaussian/output.py:922-941, 1380-1385`

```python
# _get_route
elif line.startswith("#"):
    if lines[i + 1].startswith("------"):
        route = line.lower()
    elif not lines[i + 1].startswith("------") and lines[i + 2].startswith("------"):
        route = line.lower()
        route += lines[i + 1].strip().lower()
    elif not lines[i + 1].startswith("------") and not lines[i + 2].startswith("------"):
        route = line.lower()
        route += lines[i + 1].lower()
        route += lines[i + 2].lower()
    else:
        route = None          # <-- unreachable
    return route

# energies
if len(self.mp2_energies) == 0 and len(self.oniom_energies) == 0:
    return self.scf_energies
elif len(self.mp2_energies) != 0:
    return self.mp2_energies
elif len(self.oniom_energies) != 0:   # <-- always True when reached
    return self.oniom_energies
# implicit `return None` here is unreachable
```

Two separate exhaustive-but-not-recognized-as-such `if`/`elif` chains, both with a residual branch that can never fire:

- `_get_route`'s three conditions test, respectively, "`lines[i+1]` starts with `------`" (call it `A`), "not `A` and `lines[i+2]` starts with `------`" (`B`), and "not `A` and not (`lines[i+2]` starts with `------`)" (exactly `not A and not B`'s complement-of-`B`-given-`not A`). Since `A`, `B`-given-`not A`, and `not B`-given-`not A` together cover every possibility, the trailing `else: route = None` can never execute.
- `energies`'s first condition is `mp2==0 and oniom==0`; if that's `False`, at least one of `mp2!=0`/`oniom!=0` holds. The second `elif` (`mp2!=0`) then filters out the `mp2!=0` case, so by the time the third `elif` (`oniom!=0`) is evaluated, the only way to have arrived there is `mp2==0` (else the second `elif` would already have returned) and `not (mp2==0 and oniom==0)` (else the first `if` would already have returned) -- which together force `oniom!=0`. So the third `elif`'s condition is always `True` when reached, and the function can never actually fall off the end.

**Reproduce:** confirmed via `coverage report` -- `chemsmart/io/gaussian/output.py:940` (the `else: route = None` line) and `:1384->exit` (the `energies` fallthrough arc) remain in the missing-lines list even with `tests/test_GaussianIO.py::TestGaussian16OutputAdditionalCoverage::test_route_string_spanning_two_lines`/`test_route_string_spanning_three_lines` (which exercise all three real `_get_route` branches) and the various `energies`-dependent tests (which exercise all three real `energies` branches) passing.

**Impact:** None -- dead code with no effect on behavior; both chains already correctly handle every real input via their preceding branches.

**Suggested direction:** for `_get_route`, drop the `else: route = None` (the final `elif` could become a plain `else`). For `energies`, drop the third condition's redundant re-check and make it a plain `else: return self.oniom_energies`.

---

## 96. Two `if not entries: raise ValueError(...)` empty-results guards in `utils/datasets.py` are unreachable dead code

**Location:** `chemsmart/utils/datasets.py:419-422` (`PKaTableEntry.parse_pka_table`) and `:896-899` (`PKaOutputTableEntry.parse_pka_output_table`)

```python
entries = canonical_dataset.to_entries(entry_cls=PKaTableEntry, row_offset=2)
...
if not entries:
    raise ValueError(f"No valid entries found in pKa table: {table_path}")
```

Both `entries` lists are built via `TabularDataset.to_entries()` (`chemsmart/io/datasets.py:82-86`), which returns exactly one entry per row of `self.dataframe` via `self.dataframe.iterrows()` -- so `len(entries) == len(dataframe)`. The `dataframe` in both cases comes from `TabularDataset.parse_table(...)` (or a same-row-count column-subset/rename of it), which already raises its own `ValueError: No valid entries found in table: ...` whenever `df.empty` (`chemsmart/io/datasets.py:60-61`) *before* either of these two functions ever reaches its own `to_entries()` call. So by the time `entries = ...to_entries(...)` executes, `dataframe` (and therefore `entries`) is already guaranteed non-empty, and the `if not entries:` guard's true arm can never fire.

**Reproduce:** `tests/test_utils.py::TestPKaTableParsing::test_parse_pka_table_empty_raises` confirms a comment-only file raises `ValueError` matching "No valid entries" -- but from `TabularDataset.parse_table`'s own check (`chemsmart/io/datasets.py:60-61`, message "No valid entries found in **table**: ...") rather than either of these two message variants ("... in **pKa** table: ..." / "... in pKa **output** table: ..."), which the test's partial regex match doesn't distinguish. `coverage report -m --include="*utils/datasets.py"` confirms lines 420 and 897 remain unreached across the full test suite.

**Impact:** None -- purely redundant defensive code with no behavioral effect; the earlier, differently-worded guard in `TabularDataset.parse_table` already covers this case.

**Suggested direction:** no action needed; could be removed (or the message text merged into `TabularDataset.parse_table`'s guard) now that the duplication is confirmed.

---

## 97. `PyMOLJobRunner._add_coordinates_labels`'s `elif prepend_string.startswith("D"):` is a tautological branch that is always true when reached

**Location:** `chemsmart/jobs/mol/runner.py:543-559`

```python
prepend_string_list = get_prepend_string_list_from_modred_free_format(
    input_modred=job.coordinates, program="pymol"
)
for prepend_string in prepend_string_list:
    if prepend_string.startswith("B"):
        distances.append(...)
    elif prepend_string.startswith("A"):
        angles.append(...)
    elif prepend_string.startswith("D"):
        dihedrals.append(...)
```

Same pattern as bug #81 (`jobs/orca/writer.py`'s analogous modred/scan-coordinate loop): `get_prepend_string_for_modred` (`chemsmart/utils/utils.py:1221-1243`), which `get_prepend_string_list_from_modred_free_format` calls for every coordinate group, raises `ValueError` for any coordinate list not of length 2, 3, or 4 and otherwise returns exactly `"B"`, `"A"`, or `"D"`. So every `prepend_string` reaching this loop starts with one of exactly those three letters -- by the time the `elif ... "D"` condition is evaluated (i.e. the string didn't start with `"B"` or `"A"`), it is guaranteed to start with `"D"`, and the condition's false arm (which would silently skip the entry, appending it to none of `distances`/`angles`/`dihedrals`) can never execute.

**Reproduce:** `tests/test_PyMOLJobs.py::TestPyMOLJobRunnerHelpers::test_add_coordinates_labels_handles_angles_and_dihedrals` exercises a mix of bond/angle/dihedral coordinates (including a dihedral followed by another entry, to force the loop to continue past a matched "D"); `coverage report -m --include="*jobs/mol/runner.py"` still shows arc `556->549` (the elif's false arm) as unreachable regardless of coordinate combinations, confirming the same fixed-domain guarantee as bug #81.

**Impact:** None -- purely redundant code with no behavioral effect, since `get_prepend_string_for_modred`'s return domain is fixed to exactly three values.

**Suggested direction:** no action needed; could be simplified to a plain `else:` now that the three-value domain is confirmed (consistent with the suggested fix for bug #81's identical pattern).

---

## 98. Five `ORCAOutput` energy-component `*_eV` properties (and `xc_energy` itself) crash instead of returning `None` when the underlying data is absent

**Location:** `chemsmart/io/orca/output.py:1322-1529` (`max_cosx_asymmetry_energy_eV`/`potential_energy_eV`/`kinetic_energy_eV`/`xc_energy_eV`/`dfet_embed_energy_eV` and their `_get_*` helpers)

```python
def _get_max_cosx_asymmetry_energy(self):
    max_cosx_asymmetry_energy_hartree = []
    for line in self.contents:
        if "Max COSX asymmetry :" in line:
            ...
            max_cosx_asymmetry_energy_hartree.append(energy_in_hartree)
    if len(max_cosx_asymmetry_energy_hartree) != 0:
        return max_cosx_asymmetry_energy_hartree
    # implicit `return None` when no matching lines found

def _get_max_cosx_asymmetry_energy_eV(self):
    max_cosx_asymmetry_energy_hartree = self._get_max_cosx_asymmetry_energy()
    if len(max_cosx_asymmetry_energy_hartree) != 0:   # <-- len(None) crashes
        ...
```

Four of these five property pairs (`max_cosx_asymmetry_energy`, `potential_energy`, `kinetic_energy`, `dfet_embed_energy`) follow the pattern above: the plain (Hartree) `_get_*` helper correctly returns `None` when no matching lines are found (via `if len(...) != 0: return ...`, with an implicit `return None` otherwise), and the plain property correctly guards with `if ... is not None:`. But each `_get_*_eV` sibling calls `len(...)` (three of them) or directly iterates (`_get_kinetic_energy_eV`, via a list comprehension with no guard at all) over that same helper's result *without* checking for `None` first -- so whenever the underlying data is absent, calling the `*_eV` property raises `TypeError: object of type 'NoneType' has no len()` (or `'NoneType' object is not iterable` for `kinetic_energy_eV`) instead of returning `None` like its Hartree sibling does.

The fifth, `xc_energy`/`xc_energy_eV`, has a related but distinct bug: `_get_xc_energy_hartree()` has *no* `if len(...) != 0:` guard at all -- it always returns a list, even an empty one, never `None`. So both `xc_energy` (the plain Hartree property) and `xc_energy_eV` see `is not None` as always `True` and unconditionally index `[-1]` into what can be an *empty* list, raising `IndexError: list index out of range` instead of returning `None`.

A further consequence: because `_get_*_eV()` can only ever return real data or raise (never `None`) for the four `TypeError` cases, the outer `*_eV` *property*'s own `is not None` guard can never see a `None` to act on -- its False arm is unreachable not because of a fixed value domain (like bugs #81/#97) but *because the bug itself forecloses the only input that would reach it*.

**Reproduce:** `tests/test_ORCAIO.py::TestORCAOutputDirectPropertyCoverage::test_ev_sibling_properties_crash_when_hartree_data_absent` -- a minimal `.out` file with none of the five marker strings present reproduces all five crashes (`TypeError` for the first four, `IndexError` for `xc_energy`/`xc_energy_eV`), while confirming `max_cosx_asymmetry_energy`/`potential_energy`/`kinetic_energy`/`dfet_embed_energy` correctly return `None` for the same input.

**Impact:** Medium -- any ORCA output that doesn't happen to print these specific energy-decomposition lines (e.g. a calculation without the relevant SCF print level, or a non-DFT method for `xc_energy`) makes these five properties entirely unusable, crashing instead of the `None`-on-absent behavior every other property in this file provides.

**Suggested direction:** for the four `TypeError` cases, add an `is not None` (or truthiness) check on the Hartree helper's result before calling `len()`/iterating in each `_get_*_eV()` method, mirroring the pattern already used correctly by the plain (non-eV) properties. For `xc_energy`/`xc_energy_eV`, add the same `if len(xc_energy_hartree) != 0: return xc_energy_hartree` guard to `_get_xc_energy_hartree()` that all four sibling Hartree helpers already have.

---

## 99. `ORCAOutput.num_forces`, `constrained_bond_lengths`/`constrained_bond_angles`/`constrained_dihedral_angles`, and `_get_input_structure_coordinates_block_in_output` -- three more instances of the "crashes/dead code when marker absent" pattern

**Location:** `chemsmart/io/orca/output.py:105-142` (`forces`/`num_forces`/`_get_forces_for_molecules`), `:249-330` (`constrained_*`/`_get_constraints`), `:170-202` (`_get_input_structure_coordinates_block_in_output`)

**(a) `num_forces` crashes instead of returning `0`:**

```python
def _get_forces_for_molecules(self):
    ...
    if len(list_of_all_forces) == 0:
        return None          # <-- forces is None, not [], when absent
    return list_of_all_forces

@cached_property
def num_forces(self):
    return len(self.forces)  # <-- len(None) crashes
```

Same shape as bug #98: `forces` correctly returns `None` when no `"CARTESIAN GRADIENT"` section exists, but `num_forces` calls `len()` on it unconditionally, so `oo.num_forces` raises `TypeError: object of type 'NoneType' has no len()` instead of `0` on any output file without a gradient calculation.

**(b) `constrained_bond_lengths`/`_angles`/`_dihedral_angles` crash instead of returning `{}`:**

```python
@cached_property
def _get_constraints(self):
    ...
    for i, line in enumerate(self.contents):
        if "Redundant Internal Coordinates" in line:
            ...
            return (constrained_bond_lengths, constrained_bond_angles, constrained_dihedral_angles)
    # implicit `return None` if the marker is never found

@property
def constrained_bond_lengths(self):
    constrained_bond_lengths, _, _ = self._get_constraints  # <-- unpacking None crashes
    return constrained_bond_lengths
```

`_get_constraints` has no fallback `return ({}, {}, {})` after its loop, so it implicitly returns `None` for any output file that never prints a `"Redundant Internal Coordinates"` block (e.g. a single-point job, or an optimization with no active constraints). All three public properties then crash with `TypeError: cannot unpack non-iterable NoneType object` instead of returning an empty dict, unlike almost every other "marker absent" property in this file.

Additionally, `_get_constraints`'s inner loop (`for j, line_j in enumerate(self.contents[i + 5:])`) has no fallback for running off the end of the file without hitting a blank-line terminator -- it still falls through to the same `return (...)` afterward, so that arc is harmless, just previously uncovered.

**(c) `_get_input_structure_coordinates_block_in_output` is unreachable dead code:**

```python
@cached_property
def input_coordinates_block(self):
    return self._get_first_structure_coordinates_block_in_output()

def _get_input_structure_coordinates_block_in_output(self):
    """In ORCA output file, the input structure is rewritten..."""
    ...
```

Despite its name closely matching `input_coordinates_block`, that property actually calls `_get_first_structure_coordinates_block_in_output` (a different, similarly-named method that scans for `"CARTESIAN COORDINATES (ANGSTROEM)"` instead of `"INPUT FILE"`). Nothing else in the codebase calls `_get_input_structure_coordinates_block_in_output` -- it is entirely unreferenced.

**Reproduce:** `tests/test_ORCAIO.py::TestORCAOutputDirectPropertyCoverage::test_num_forces_crashes_when_forces_absent`, `::test_get_constraints_absent_crashes_dependent_properties`, `::test_get_constraints_runs_to_natural_exhaustion`, and `::test_get_input_structure_coordinates_block_in_output_is_dead_code` (the last calls the method directly, since nothing else does).

**Impact:** Low-medium for (a)/(b) -- both are common cases (any output without a gradient print, or without active geometry constraints) that would currently crash callers relying on these properties as a lightweight "is this present" check. Low for (c) -- purely wasted code, but harmless since `input_coordinates_block` (the only plausibly-intended caller) already works via the other method.

**Suggested direction:** (a) guard `num_forces` with `len(self.forces) if self.forces is not None else 0`. (b) add a fallback `return ({}, {}, {})` after `_get_constraints`'s loop, matching the "always return a value, never implicitly `None`" convention used elsewhere in this file. (c) delete `_get_input_structure_coordinates_block_in_output` (and its docstring/pattern usage) as dead code, or rename it and wire it up if it was meant to be used instead of `_get_first_structure_coordinates_block_in_output`.

---

## 100. A whole family of population-analysis/dipole/rotational-constant properties crash with `IndexError` instead of returning `None`/`{}` when their section marker is absent

**Location:** `chemsmart/io/orca/output.py:1707-2262` -- `mulliken_atomic_charges`, `loewdin_atomic_charges`, `mayer_mulliken_gross_atomic_population`, `mayer_total_nuclear_charge`, `mayer_mulliken_gross_atomic_charge`, `mayer_total_valence`, `mayer_bonded_valence`, `mayer_free_valence`, `mayer_bond_orders_larger_than_zero_point_one`, `total_integrated_alpha_density`, `total_integrated_beta_density`, `_get_hirshfeld_charges_and_spins` (and its `hirshfeld_charges`/`hirshfeld_spin_densities` wrappers), `dipole_moment_electric_contribution`, `dipole_moment_nuclear_contribution`, `dipole_moment_in_au`, `dipole_moment_magnitude_in_au`, `dipole_moment_magnitude_in_debye`, `dipole_moment_along_axis_in_au`, `dipole_moment_along_axis_in_debye`, `rotational_constants_in_wavenumbers`, `rotational_constants_in_MHz` (17 distinct properties, all following the identical shape)

```python
@property
def mulliken_atomic_charges(self):
    all_mulliken_atomic_charges = []
    for i, line_i in enumerate(self.contents):
        if "MULLIKEN ATOMIC CHARGES" in line_i:
            ...
            all_mulliken_atomic_charges.append(mulliken_atomic_charges)
    return all_mulliken_atomic_charges[-1]   # <-- crashes if never appended to
```

Every one of these properties/helpers builds a list by appending once per occurrence of its section marker, then unconditionally returns `accumulator[-1]` with no `if accumulator:`/`if len(accumulator) != 0:` guard (unlike, e.g., `energies`, `_get_max_cosx_asymmetry_energy`, or `all_vibrational_frequencies`, which all correctly guard this same shape elsewhere in the file). Any ORCA output that doesn't happen to print the corresponding section -- e.g. a job run without `%output Print[P_Mayer] 1`, without a dipole calculation, or a single-point job with no `Rotational spectrum` block -- makes the property raise `IndexError: list index out of range` instead of the `None`/`{}` that every other "marker absent" property in this file returns.

This has a further knock-on effect: `rotational_constants_in_Hz` (`:2264-2275`) and `rotational_temperatures` (`:2277-2288`) both guard with `if self.rotational_constants_in_MHz is None: return None` -- but since `rotational_constants_in_MHz` can now only ever return real data or raise (never `None`), that guard's `None` branch is unreachable *because the bug forecloses the only input that would reach it*, exactly the same secondary effect documented for bug #98's `*_eV` properties. The `IndexError` simply propagates up through both.

**Reproduce:** `tests/test_ORCAIO.py::TestORCAOutputDirectPropertyCoverage::test_population_dipole_rotational_properties_crash_when_absent` -- a minimal `.out` file with none of the relevant section markers reproduces `IndexError` for all 21 attributes (17 directly-affected properties/helpers plus `rotational_constants_in_Hz`/`rotational_temperatures` inheriting the crash).

**Impact:** Medium-high -- this is the single largest source of crash-instead-of-`None` behavior in the file by property count. Any ORCA output missing one of these fairly common but non-default print sections (population analyses, dipole moments, rotational spectra) makes the corresponding property entirely unusable rather than gracefully reporting "not present."

**Suggested direction:** add `if accumulator: return accumulator[-1]` / `return accumulator[-1] if accumulator else None` (or the dict-equivalent `{}`/`None`) to each of the 17 directly-affected properties, mirroring the guard pattern already used correctly by `energies`, `_get_max_cosx_asymmetry_energy`, and `all_vibrational_frequencies` elsewhere in this same file. No change needed to `rotational_constants_in_Hz`/`rotational_temperatures` once their dependency is fixed.
