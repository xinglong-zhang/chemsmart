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
