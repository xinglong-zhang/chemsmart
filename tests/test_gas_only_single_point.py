"""A gas-phase single point must be expressible in a gas-phase project.

``sp`` reads ``solv:`` because ChemSmart's canonical workflow optimises in gas
phase and takes the single point in solvent.  The fallback around that rule used
to be one-sided: a ``solv:``-only project fed every job type, while a
``gas:``-only project fed ``sp`` nothing at all.  The settings then reached the
ORCA writer carrying no level of theory, which raised

    ValueError: Warning: neither ab initio nor DFT is specified!

after having already created a zero-byte ``.inp``.

A live agent session hit exactly this.  Asked for a gas-phase B3LYP/aug-cc-pVTZ
single point it wrote ``gas: {functional: b3lyp, basis: aug-cc-pVTZ}`` -- the
only section name that means "gas phase" -- and got three identical opaque
failures, including one on a minimal functional+basis project.  It concluded
that ORCA could not materialise a single point in this environment and switched
program.  The project YAML was right and the hub refused it, which is the
inverse of what the project exists to do.

These tests pin the repaired rule for both phase-keyed programs:

* a ``gas:``-only project feeds ``sp`` its level of theory;
* borrowing ``gas:`` does not borrow its ``freq``, because ``gas:`` describes
  the optimisation and a single point is not an opt+freq;
* an explicit ``freq:`` under ``solv:`` still overrides sp's default, because
  ``solv:`` *is* sp's own section.
"""

import os
import tempfile

import pytest

from chemsmart.jobs.settings import read_molecular_job_yaml

PHASE_KEYED_PROGRAMS = ("orca", "gaussian")


def _configs(body, program):
    """Read ``body`` as a project YAML for ``program``."""
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "project.yaml")
        with open(path, "w") as handle:
            handle.write(body)
        return read_molecular_job_yaml(path, program=program)


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_gas_only_project_gives_sp_a_level_of_theory(program):
    """The regression: ``sp`` used to receive neither functional nor basis."""
    configs = _configs(
        "gas:\n  functional: b3lyp\n  basis: aug-cc-pVTZ\n", program
    )

    assert configs["sp"]["functional"] == "b3lyp"
    assert configs["sp"]["basis"] == "aug-cc-pVTZ"


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_gas_only_single_point_does_not_inherit_gas_freq(program):
    """Borrowing ``gas:`` for its method must not turn ``sp`` into opt+freq."""
    configs = _configs(
        "gas:\n  functional: b3lyp\n  basis: aug-cc-pVTZ\n  freq: true\n",
        program,
    )

    assert configs["sp"]["functional"] == "b3lyp"
    assert configs["sp"]["freq"] is False
    # The section that owns frequencies keeps them.
    assert configs["opt"]["freq"] is True


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_solv_still_wins_for_sp_when_both_sections_exist(program):
    """The documented split is untouched: ``solv:`` feeds ``sp``, ``gas:`` the rest."""
    configs = _configs(
        "gas:\n  functional: b3lyp\n  basis: 6-31G*\n"
        "solv:\n  functional: m06-2x\n  basis: aug-cc-pVTZ\n",
        program,
    )

    assert configs["sp"]["functional"] == "m06-2x"
    assert configs["sp"]["basis"] == "aug-cc-pVTZ"
    assert configs["opt"]["functional"] == "b3lyp"
    assert configs["opt"]["basis"] == "6-31G*"


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_explicit_freq_under_solv_still_overrides_the_sp_default(program):
    """``solv:`` is sp's own section, so an author writing ``freq:`` there means it."""
    configs = _configs(
        "gas:\n  functional: b3lyp\n  basis: 6-31G*\n"
        "solv:\n  functional: m06-2x\n  basis: aug-cc-pVTZ\n  freq: true\n",
        program,
    )

    assert configs["sp"]["freq"] is True


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_solv_only_project_is_unchanged(program):
    """The pre-existing one-phase branch keeps feeding every job type."""
    configs = _configs(
        "solv:\n  functional: b3lyp\n  basis: aug-cc-pVTZ\n"
        "  solvent_model: cpcm\n  solvent_id: water\n",
        program,
    )

    assert configs["sp"]["functional"] == "b3lyp"
    assert configs["opt"]["functional"] == "b3lyp"
    assert configs["sp"]["solvent_model"] == "cpcm"


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_an_explicit_sp_section_still_wins_over_the_inherited_gas(program):
    """A stage key is the most specific statement, so it overrides inheritance."""
    configs = _configs(
        "gas:\n  functional: b3lyp\n  basis: 6-31G*\n"
        "sp:\n  functional: m06-2x\n  basis: def2-QZVP\n",
        program,
    )

    assert configs["sp"]["functional"] == "m06-2x"
    assert configs["sp"]["basis"] == "def2-QZVP"


@pytest.mark.parametrize("program", PHASE_KEYED_PROGRAMS)
def test_a_project_with_neither_phase_section_says_so(program):
    """Naming the mistake in the file a chemist edits stays the failure mode."""
    with pytest.raises(ValueError) as excinfo:
        _configs("functional: b3lyp\nbasis: aug-cc-pVTZ\n", program)

    message = str(excinfo.value)
    assert "'gas'" in message and "'solv'" in message
    assert "basis, functional" in message
