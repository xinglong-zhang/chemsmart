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
