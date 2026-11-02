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
