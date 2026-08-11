"""A method or basis name means the same chemistry in any case.

This made Gaussian unusable to the agent, and the failure was invisible.

A live session wrote the standard ``6-31G(d)``. ChemSmart's own writer emitted
it into the route line; ChemSmart's own reader returned ``6-31g(d)``; and the
preview validator compared the two case-sensitively and refused the input:

    preview.semantic.mismatch | basis | expected '6-31G(d)' | observed '6-31g(d)'

No spelling could have passed, because the observed side always comes back
from the parser that lowercases it. The session recompiled Gaussian four times
and ORCA twice, guessing at ``functional`` vs ``ab_initio``, ``forces``,
``numfreq``, ``reference``, ``ri_approximation`` and ``6-31G*`` -- none of
which were wrong -- then abandoned both programs.

The rule already existed for ORCA alone, applied through the ``native_input``
branch, which is why only Gaussian's refusal looked mysterious.

Deliberately NOT fixed here: ``6-31G*`` and ``6-31G(d)`` denote the same Pople
basis, but treating them as equal is a chemistry-aliasing judgement, not a
round-trip defect, and it could mask a real difference. Case is different: the
round trip provably does not preserve it.
"""

from chemsmart.agent.program_verifiers import _settings_match


class _Parsed:
    """Stands in for settings read back from a generated input file."""

    def __init__(self, **fields):
        self.__dict__.update(fields)


def test_the_standard_spelling_is_accepted():
    parsed = _Parsed(
        ab_initio="hf", basis="6-31g(d)", freq=True, charge=0, multiplicity=1
    )
    expected = {
        "ab_initio": "hf",
        "basis": "6-31G(d)",
        "freq": True,
        "charge": 0,
        "multiplicity": 1,
    }
    assert _settings_match(parsed, expected) == []


def test_case_folding_applies_without_an_orca_native_input():
    """The defect was that this only ran for ORCA."""

    parsed = _Parsed(functional="B3LYP")
    assert _settings_match(parsed, {"functional": "b3lyp"}) == []


def test_a_different_basis_is_still_refused():
    parsed = _Parsed(basis="6-31g(d)")
    findings = _settings_match(parsed, {"basis": "cc-pVTZ"})
    assert [item.field for item in findings] == ["basis"]


def test_a_different_method_is_still_refused():
    parsed = _Parsed(ab_initio="hf")
    findings = _settings_match(parsed, {"ab_initio": "mp2"})
    assert [item.field for item in findings] == ["ab_initio"]


def test_a_non_string_difference_is_still_refused():
    """Case folding must not become a general leniency."""

    parsed = _Parsed(freq=False, charge=0)
    findings = _settings_match(parsed, {"freq": True, "charge": 1})
    assert sorted(item.field for item in findings) == ["charge", "freq"]


def test_a_field_absent_from_the_generated_input_is_reported():
    parsed = _Parsed(basis="6-31g(d)")
    findings = _settings_match(parsed, {"solvent_model": "smd"})
    assert findings and findings[0].field == "solvent_model"
