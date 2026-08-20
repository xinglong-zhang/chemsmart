"""Conformer degeneracies are not optional, and there is one implementation.

`chemsmart/analysis/aggregation.py` states the rule plainly: degeneracies "are
not optional in practice… omitting it is a scientific error rather than an
approximation", with the worked case that counting n-butane's two gauche wells
once gives 82% anti where the correct treatment gives 70%.

The ensemble-averaging path hand-rolled its own partition function and had no
degeneracy parameter at all, so that error was unreachable to fix even by a
caller who knew the multiplicities. These tests pin the arithmetic against the
documented case rather than against the implementation, and check that the
default is unchanged so no existing report moves silently.
"""

from __future__ import annotations

import pytest

from chemsmart.analysis.aggregation import boltzmann_populations

# Relative conformer energies of n-butane in kcal/mol: one anti well and two
# symmetry-equivalent gauche wells.
_ANTI_GAUCHE = (0.0, 0.9, 0.9)


def test_the_documented_n_butane_case_reproduces():
    """Two gauche wells counted as two states, not one."""

    populations = boltzmann_populations(_ANTI_GAUCHE, temperature=298.15)

    anti = populations[0]
    gauche = populations[1] + populations[2]
    assert anti == pytest.approx(0.70, abs=0.02)
    assert gauche == pytest.approx(0.30, abs=0.02)
    assert sum(populations) == pytest.approx(1.0)


def test_collapsing_the_pair_gives_the_documented_wrong_answer():
    """The error the rule exists to prevent, reproduced deliberately."""

    populations = boltzmann_populations((0.0, 0.9), temperature=298.15)

    assert populations[0] == pytest.approx(0.82, abs=0.02)


def test_a_degeneracy_of_two_equals_two_explicit_states():
    """Supplying the multiplicity is the same physics as listing the wells."""

    explicit = boltzmann_populations(_ANTI_GAUCHE, temperature=298.15)
    weighted = boltzmann_populations(
        (0.0, 0.9), temperature=298.15, degeneracies=(1.0, 2.0)
    )

    assert weighted[0] == pytest.approx(explicit[0])
    assert weighted[1] == pytest.approx(explicit[1] + explicit[2])


def test_the_ensemble_path_now_takes_degeneracies():
    """The parameter has to exist, or correctness is unreachable."""

    import inspect

    from chemsmart.analysis.thermochemistry import (
        BoltzmannAverageThermochemistry,
    )

    signature = inspect.signature(BoltzmannAverageThermochemistry.__init__)
    assert "degeneracies" in signature.parameters


def test_the_ensemble_path_uses_the_shared_weighting():
    """One implementation of this physics, not two.

    The duplicate computed its own partition function, so a fix to the shared
    routine could not reach it -- which is how the two drifted apart in the
    first place.
    """

    import inspect

    from chemsmart.analysis import thermochemistry

    source = inspect.getsource(thermochemistry.BoltzmannAverageThermochemistry)
    assert "boltzmann_populations(" in source
    assert (
        "partition_function" not in source
    ), "the hand-rolled partition function must be gone, not shadowed"


def test_the_default_is_unchanged():
    """No existing report moves silently: absent degeneracies mean one each."""

    assert boltzmann_populations(
        _ANTI_GAUCHE, temperature=298.15
    ) == pytest.approx(
        boltzmann_populations(
            _ANTI_GAUCHE, temperature=298.15, degeneracies=(1.0, 1.0, 1.0)
        )
    )
