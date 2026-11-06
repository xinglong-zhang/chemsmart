# Review cases

Regression examples for the review rules in SKILL.md. Check the current code before applying a case; the symbols shown are illustrations, not a current API specification. [review-evidence.md](review-evidence.md) records the source of each case.

## Dynamic `getattr` lookup

Source: [User feedback: dynamic getattr lookups](review-evidence.md#user-feedback-dynamic-getattr-lookups).

### Situation

`Molecule.__init__` sets every dataset attribute, and a schema supplies the field names at runtime:

```python
for key, unit in schema.property_units.items():
    value = getattr(molecule, key)
    if value is None:
        continue
    properties[key] = {"value": value, "unit": unit}
```

### Expected behaviour

Do not report this `getattr` as a defect or an OOP violation. The attribute name is runtime data, the constructor defines the attribute, and a wrong key raises `AttributeError`. Do not suggest replacing it with a hand-written dictionary or per-dataset attribute access only to remove the call. In particular, do not suggest `molecule.dataset_values.get(key)`: a field missing from that dictionary would be skipped silently and never stored.

Separate design concerns may still be reported as an OOP/API recommendation. Here, that means asking whether the separate attributes belong on `Molecule`, since `molecule.homo` is eV for QM9 and Hartree for revQM9. Report that concern for its own reason, not because `getattr` appears.

### Counterexample: this should still be flagged

The name is fixed, and the default hides whether the attribute was ever created:

```python
scratch_dir = getattr(job, "scratch_parent_folder", None)
```

Ask for the attribute to be initialised in the constructor, as in PR #687, so its availability does not depend on which code path has run. The same concern applies to `getattr(molecule, key, None)` when the attribute should always exist.

## Real test data

Source: [User feedback: test data and test support layout](review-evidence.md#user-feedback-test-data-and-test-support-layout).

### Situation

A PR adds a reader for a program's results file, and a test builds a failed run by hand because no committed output covers one:

```python
write_results(path, status={"normal_termination": False}, results={})
assert Output(path).normal_termination is False
```

### Expected behaviour

Ask whether the failed run can come from the program itself: a real calculation stopped by its iteration limit, committed as a small fixture, or run in the test if it is fast. A hand-written file only shows that the reader agrees with what the test's author assumed the program writes. Keep this as a recommendation, and do not report it as a defect in the reader.

### Counterexample: this should not be flagged

Some inputs cannot be real or cannot run in the test environment, so synthetic data or a narrow patch is appropriate:

```python
(tmp_path / "corrupt.h5").write_bytes(b"not an hdf5 file")
assert Output(tmp_path / "corrupt.out").normal_termination is False

with patch.object(Molecule, "from_pubchem", return_value=molecule):
    job = Job.from_pubchem("962", settings=settings, label="fetched")
```

## Test support layout

Source: [User feedback: test data and test support layout](review-evidence.md#user-feedback-test-data-and-test-support-layout).

### Situation

A new test file defines its own fixture and assertion helper, even though only that file uses them so far:

```python
@pytest.fixture()
def runner_options(pbs_server):
    return {"server": pbs_server, "scratch": False, "fake": True}


def _assert_runner_was_configured(job, options):
    assert job.jobrunner.server is options["server"]
```

### Expected behaviour

Check `tests/conftest.py` for the current layout, then ask for the fixture to move to the matching `tests/fixtures/<area>.py` module and the helper to `tests/helpers.py`. Report this as a consistency finding, not a correctness defect.

### Counterexample: this should not be flagged

A class or function that only makes sense inside one test, because it captures that test's local state, can stay there:

```python
def test_rejects_arbitrary_pickle_globals(tmp_path):
    marker = tmp_path / "unpickled"

    class Payload:
        def __reduce__(self):
            return os.mkdir, (str(marker),)
```

Also do not ask for a mass migration of private helpers in test files the PR does not change.

## Cached values on mutable objects

Source: [User feedback: cached_property and property](review-evidence.md#user-feedback-cached_property-and-property).

### Situation

A property is cached on a class whose inputs can be reassigned, and a caller copies the object and then changes it:

```python
class Molecule:
    @cached_property
    def distance_matrix(self):
        return cdist(self.positions, self.positions)

    def displaced(self, new_positions):
        molecule = self.copy()
        molecule.positions = new_positions
        return molecule
```

### Expected behaviour

Report a correctness finding. If `distance_matrix` was read before `displaced()` is called, the new molecule returns the old distances, because the copy keeps the cached value. Ask for `@property`, a regression test that compares the displaced molecule with a freshly built one, and a check of callers:

```python
distances = self.distance_matrix  # compute once, not per atom pair
for i in range(n_atoms):
    for j in range(i + 1, n_atoms):
        bond_length = distances[i, j]
```

### Counterexample: this should not be flagged

The value comes from a file that the object parses and nothing reassigns its inputs, so caching avoids repeated parsing and cannot go stale in normal use:

```python
class OutputFile:
    @cached_property
    def contents(self):
        with open(self.filename) as f:
            return f.readlines()
```

Also do not report a stale-value bug when no code path changes the input after the value is read. Raise that as an OOP/API recommendation if the class permits such a change.

## Duplicating an existing method

Source: [User feedback: external level of theory, PR #732](review-evidence.md#user-feedback-external-level-of-theory-pr-732).

### Situation

A runner gains a helper to copy results back from scratch, including folders:

```python
class JobRunner:
    @staticmethod
    def _copy_to_job_folder(path, job):
        if os.path.isdir(path):
            shutil.copytree(path, os.path.join(job.folder, os.path.basename(path)), dirs_exist_ok=True)
        else:
            shutil.copy(path, job.folder)
```

### Expected behaviour

Report two OOP/API findings. Ownership: the helper only uses `job.folder`, so it belongs to the job, not the runner that calls it. Reuse: `Job.backup_file(path, folder=job.folder)` already copies files and folders, and its rename-aside behaviour for an existing folder keeps two runs apart, where `dirs_exist_ok=True` would merge them. Ask the runners to call `job.backup_file(...)` and not to add `Job.copy_to_job_folder()` or `Job.backup_folder()`. Read the existing method before judging it: it is public, and its semantics are the better fit.

### Counterexample: this should not be flagged

A new method whose behaviour differs from the existing one in a way callers need, with the difference stated:

```python
def backup_file(self, file, folder=None, remove=False): ...    # copies or moves one path
def backup(self, **kwargs): ...                                 # backs up this job's known outputs into bk.<timestamp>/
```

Different responsibilities, not duplicates. Here, extend `backup()`'s host implementations to include new outputs instead of adding another method.

## Nested subcommands of reused Click commands

Source: [User feedback: external level of theory, PR #732](review-evidence.md#user-feedback-external-level-of-theory-pr-732).

### Situation

A new group reuses existing commands to apply a new setting:

```python
@gaussian.group("external", cls=MyGroup)
def external(ctx, backend, ...):
    ctx.obj["job_settings"].external = spec

for command in (opt, ts, modred, scan, irc, sp):
    external.add_command(command)
```

`opt` and `ts` already have a `qmmm` subcommand, which builds fresh settings from a fixed list of inherited fields that omits `external`.

### Expected behaviour

Report a correctness finding: `external -b xtb opt qmmm ...` is listed by `--help`, exits 0, and runs ONIOM at the project's level of theory under the label `*_opt_xtb_qmmm`. Ask for the unsupported combination to raise `click.UsageError`, with a `WIP` note if support is planned, and a CLI test for each host program.

### Counterexample: this should not be flagged

The reused commands are leaf commands, or every nested subcommand inherits the new setting through the same merge path, verified by a test.

## Parser patterns from guessed output

Source: [User feedback: external level of theory, PR #732](review-evidence.md#user-feedback-external-level-of-theory-pr-732).

### Situation

A parser gains a pattern for a line no committed output contains, tested against a fixture made by rewriting another program's log:

```python
external_recovered_energy_pattern = r"Recovered energy=\s*([-+]?\d*\.\d+)"
```

### Expected behaviour

Ask for a real output (or a short real run) before merging the pattern, and check how many times the line appears per step. In this case the program prints two energy lines per evaluation, and matching both would double the energies and shift their pairing with structures without any error. Keep it as a recommendation unless a concrete mismatch is shown.

### Counterexample: this should not be flagged

The pattern is tested against a committed real output that contains the line, and the test asserts that energy i pairs with structure i.

## Program name in prose

Source: [User feedback: external level of theory, PR #732](review-evidence.md#user-feedback-external-level-of-theory-pr-732) and [PR #612](https://github.com/xinglong-zhang/chemsmart/pull/612#discussion_r3563653821).

### Situation

```python
"""Job runners set it to ``<label>_external.log``. Chemsmart's launcher ..."""
sys.stderr.write("chemsmart external calculation failed:\n")
f"# chemsmart external log: host {host}"
```

### Expected behaviour

Report a consistency finding listing every occurrence: use `CHEMSMART` in docstrings, comments, messages and file headers. Find them by searching the added lines, not by reading.

### Counterexample: this should not be flagged

Identifiers stay lowercase:

```python
from chemsmart.io.external import build_engine
EXTERNAL_LOG_ENV = "CHEMSMART_EXTERNAL_LOG"
script = "chemsmart_external_xtb_1a2b3c4d.py"   # chemsmart run gaussian ... external ...
```

## Constants and regex defined locally

Source: [User feedback: external level of theory, PR #732](review-evidence.md#user-feedback-external-level-of-theory-pr-732) and [PR #310](https://github.com/xinglong-zhang/chemsmart/pull/310).

### Situation

```python
# chemsmart/io/external.py
BOHR_TO_ANG = 0.529177210903   # CODATA 2018
EV_PER_HA = 27.211386245988

# chemsmart/jobs/orca/settings.py
tokens = re.findall(r"[^\s!]+", route_string)
```

### Expected behaviour

Ask for the conversion factors to move to `chemsmart/utils/constants.py`, derived from `ase.units` like the rest of the repository, so values agree across modules, and for the pattern to move to `chemsmart/utils/repattern.py` with a descriptive name and comment.

### Counterexample: this should not be flagged

Feature configuration is not a physical constant and can stay with its feature:

```python
HESSIAN_MODES = ("auto", "analytic", "fd")
DEFAULT_FD_STEP = 0.005  # Angstrom
```

## Removing a correct test to pass new code

Source: [Extended historical sample](review-evidence.md#extended-historical-sample).

### Situation

A refactor of output parsing changes an existing expected value and deletes a test data file:

```diff
-    assert len(g16_output.all_structures) == 6
+    assert len(g16_output.all_structures) == 5
-tests/data/GaussianTests/outputs/CH3CH3.log
```

### Expected behaviour

Ask why each expected value changed and why the file was removed, and check the original value against the output file. If the original was correct, report a correctness finding: the change hides a regression.

### Counterexample: this should not be flagged

The PR shows the old expected value was wrong, with evidence from the file or independent software, and updates the assertion with that explanation.

## Keyword catch-alls

Source: [Technical review of the guidance](review-evidence.md#technical-review-of-the-guidance).

### Situation

A new runner option is added through a catch-all, and a command reuses a decorator whose options it ignores:

```python
class NewJobRunner(JobRunner):
    def __init__(self, server, scratch=None, **kwargs):
        super().__init__(server, scratch=scratch, **kwargs)
        self.keep_files = kwargs.get("keep_files", False)

@click_jobrunner_options        # adds --num-cores, --mem-gb, ...
@click.command()
def summary(**kwargs): ...      # never uses them
```

### Expected behaviour

Report an OOP/API finding: `NewJobRunner(server, keep_file=True)` runs with `keep_files=False` and no error, and `summary --num-cores 8` silently does nothing. Ask for `keep_files` to be an explicit parameter and for `summary` to use only the option groups it needs. Do not ask to remove existing catch-alls in `Job` or `JobRunner` in the same PR.

### Counterexample: this should not be flagged

`**kwargs` forwarded unchanged to a callee whose signature validates them:

```python
def create(cls, molecule, **kwargs):
    return cls(molecule=molecule, **kwargs)   # cls.__init__ rejects unknown keys
```

## Net spin versus radical character

Source: [Technical review of the guidance](review-evidence.md#technical-review-of-the-guidance).

### Situation

A PR adds a property for radical analysis on top of the existing definition:

```python
@property
def is_radical(self):
    return self.num_unpaired_electrons > 0   # num_unpaired_electrons = multiplicity - 1
```

### Expected behaviour

Report a correctness finding: a broken-symmetry open-shell singlet (multiplicity 1, `<S**2>` near 1, e.g. after `stable=opt`) is a diradical, but `is_radical` returns `False`; `tests/data/GaussianTests/outputs/link/oxygen_openshell_singlet_sp_link.log` (`<S**2>` = 1.005) is a real test case. Ask for the property to use the spin-density or orbital information it needs, or to be named and documented as a net-spin check, with an open-shell singlet test using a real output.

### Counterexample: this should not be flagged

A property documented as the net spin, used only for high-spin bookkeeping such as xTB's `--uhf`:

```python
def uhf(self):
    """xTB --uhf: N_alpha - N_beta = multiplicity - 1."""
    return self.multiplicity - 1
```

