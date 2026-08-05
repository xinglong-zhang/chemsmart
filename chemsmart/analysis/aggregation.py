"""Post-processing operations that turn computed artifacts into an observable.

A multi-stage workflow rarely ends at a program invocation.  Complete-basis-set
extrapolation, a state-energy difference, and a Boltzmann average all combine
several finished calculations into the number the question actually asked for.
None of them is a program call, so none of them has a place in a DAG whose only
node kind is "run a program" -- which is how a requested observable ends up with
no executor and gets dropped during repair.

These operations are deterministic host arithmetic.  They read energies that a
program already produced and never invent one.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

#: Hartree to kcal/mol (CODATA-consistent with the rest of the package).
HARTREE_TO_KCAL_PER_MOL = 627.5094740631
#: Hartree to kJ/mol.
HARTREE_TO_KJ_PER_MOL = 2625.4996394799
#: Hartree to eV.
HARTREE_TO_EV = 27.211386245988
#: Gas constant in kcal/(mol K), for Boltzmann populations.
GAS_CONSTANT_KCAL = 1.987204258640832e-3

ENERGY_UNIT_FACTORS = {
    "hartree": 1.0,
    "kcal/mol": HARTREE_TO_KCAL_PER_MOL,
    "kj/mol": HARTREE_TO_KJ_PER_MOL,
    "ev": HARTREE_TO_EV,
}


class AggregationError(ValueError):
    """Raised when an aggregation request is not physically well posed."""


def _cardinal(value) -> int:
    """Return a validated basis-set cardinal number."""

    number = int(value)
    if number < 2:
        raise AggregationError(
            f"basis cardinal number must be at least 2, got {value!r}"
        )
    return number


def convert_energy(value: float, unit: str) -> float:
    """Convert a Hartree energy into ``unit``."""

    key = str(unit).strip().lower()
    if key not in ENERGY_UNIT_FACTORS:
        raise AggregationError(
            f"unsupported energy unit {unit!r}; "
            f"expected one of {sorted(ENERGY_UNIT_FACTORS)}"
        )
    return float(value) * ENERGY_UNIT_FACTORS[key]


@dataclass(frozen=True)
class CBSExtrapolationResultV1:
    """A two-point complete-basis-set extrapolation, reported by component."""

    scf_energy: float
    correlation_energy: float
    total_energy: float
    smaller_cardinal: int
    larger_cardinal: int
    scf_alpha: float
    unit: str = "hartree"


def extrapolate_correlation_helgaker(
    *,
    smaller_cardinal,
    larger_cardinal,
    smaller_correlation_energy: float,
    larger_correlation_energy: float,
) -> float:
    """Two-point ``X**-3`` correlation extrapolation (Helgaker et al.).

    The correlation energy converges as ``E(X) = E_inf + B * X**-3``, so two
    cardinals determine the limit:

        ``E_inf = (X**3 * E_X - Y**3 * E_Y) / (X**3 - Y**3)``
    """

    x = _cardinal(smaller_cardinal)
    y = _cardinal(larger_cardinal)
    if x >= y:
        raise AggregationError(
            "the larger cardinal must exceed the smaller one; "
            f"got {x} and {y}"
        )
    x3, y3 = x**3, y**3
    return (x3 * float(smaller_correlation_energy) - y3 * float(larger_correlation_energy)) / (
        x3 - y3
    )


def extrapolate_scf_exponential(
    *,
    smaller_cardinal,
    larger_cardinal,
    smaller_scf_energy: float,
    larger_scf_energy: float,
    alpha: float = 3.9,
) -> float:
    """Two-point exponential SCF extrapolation.

    The Hartree-Fock energy converges exponentially rather than as a power of
    the cardinal number, so it must be extrapolated separately from the
    correlation energy:

        ``E(L) = E_inf + A * exp(-alpha * sqrt(L))``

    Solving the two-point system gives ``A``, then ``E_inf``.
    """

    x = _cardinal(smaller_cardinal)
    y = _cardinal(larger_cardinal)
    if x >= y:
        raise AggregationError(
            "the larger cardinal must exceed the smaller one; "
            f"got {x} and {y}"
        )
    if float(alpha) <= 0:
        raise AggregationError("alpha must be positive")
    fx = math.exp(-float(alpha) * math.sqrt(x))
    fy = math.exp(-float(alpha) * math.sqrt(y))
    amplitude = (float(smaller_scf_energy) - float(larger_scf_energy)) / (fx - fy)
    return float(smaller_scf_energy) - amplitude * fx


def cbs_extrapolation(
    *,
    smaller_cardinal,
    larger_cardinal,
    smaller_scf_energy: float,
    larger_scf_energy: float,
    smaller_correlation_energy: float,
    larger_correlation_energy: float,
    scf_alpha: float = 3.9,
) -> CBSExtrapolationResultV1:
    """Extrapolate SCF and correlation separately and report both.

    The two components converge by different laws, so a single extrapolation of
    the total energy is not equivalent to this and the components are reported
    rather than folded away.
    """

    scf = extrapolate_scf_exponential(
        smaller_cardinal=smaller_cardinal,
        larger_cardinal=larger_cardinal,
        smaller_scf_energy=smaller_scf_energy,
        larger_scf_energy=larger_scf_energy,
        alpha=scf_alpha,
    )
    correlation = extrapolate_correlation_helgaker(
        smaller_cardinal=smaller_cardinal,
        larger_cardinal=larger_cardinal,
        smaller_correlation_energy=smaller_correlation_energy,
        larger_correlation_energy=larger_correlation_energy,
    )
    return CBSExtrapolationResultV1(
        scf_energy=scf,
        correlation_energy=correlation,
        total_energy=scf + correlation,
        smaller_cardinal=_cardinal(smaller_cardinal),
        larger_cardinal=_cardinal(larger_cardinal),
        scf_alpha=float(scf_alpha),
    )


@dataclass(frozen=True)
class StateEnergyDifferenceV1:
    """A signed state-to-state energy difference with its direction stated."""

    value: float
    unit: str
    upper_state_label: str
    lower_state_label: str
    convention: str


def state_energy_difference(
    *,
    upper_state_energy: float,
    lower_state_energy: float,
    upper_state_label: str = "upper",
    lower_state_label: str = "lower",
    unit: str = "kcal/mol",
) -> StateEnergyDifferenceV1:
    """Return ``E(upper) - E(lower)``, carrying its direction.

    A difference is meaningless without its direction, so the record names both
    states and states the convention explicitly rather than leaving the reader
    to infer the sign.
    """

    difference = float(upper_state_energy) - float(lower_state_energy)
    return StateEnergyDifferenceV1(
        value=convert_energy(difference, unit),
        unit=str(unit).strip().lower(),
        upper_state_label=str(upper_state_label),
        lower_state_label=str(lower_state_label),
        convention="E(upper) - E(lower)",
    )


def boltzmann_populations(
    energies: tuple[float, ...],
    *,
    temperature: float = 298.15,
    unit: str = "kcal/mol",
) -> tuple[float, ...]:
    """Return normalised Boltzmann populations for relative energies.

    ``energies`` may be absolute or relative; only differences matter, and the
    lowest is used as the reference to keep the exponentials finite.
    """

    if not energies:
        raise AggregationError("at least one energy is required")
    if float(temperature) <= 0:
        raise AggregationError("temperature must be positive")
    key = str(unit).strip().lower()
    if key == "hartree":
        values = [convert_energy(item, "kcal/mol") for item in energies]
    elif key == "kcal/mol":
        values = [float(item) for item in energies]
    elif key == "kj/mol":
        values = [float(item) / 4.184 for item in energies]
    else:
        raise AggregationError(
            f"unsupported energy unit {unit!r} for Boltzmann weighting"
        )
    reference = min(values)
    weights = [
        math.exp(-(item - reference) / (GAS_CONSTANT_KCAL * float(temperature)))
        for item in values
    ]
    total = sum(weights)
    return tuple(item / total for item in weights)


def boltzmann_average(
    values: tuple[float, ...],
    energies: tuple[float, ...],
    *,
    temperature: float = 298.15,
    unit: str = "kcal/mol",
) -> float:
    """Return the Boltzmann-weighted average of ``values``."""

    if len(values) != len(energies):
        raise AggregationError(
            "values and energies must have the same length"
        )
    populations = boltzmann_populations(
        energies, temperature=temperature, unit=unit
    )
    return sum(
        float(value) * weight for value, weight in zip(values, populations)
    )


#: Operations an ``aggregate`` workflow node may declare.  ChemSmart itself is
#: the engine for these; they never call an external program.
AGGREGATE_OPERATIONS = (
    "boltzmann_average",
    "cbs_extrapolation",
    "state_energy_difference",
)


__all__ = [
    "AGGREGATE_OPERATIONS",
    "AggregationError",
    "CBSExtrapolationResultV1",
    "HARTREE_TO_EV",
    "HARTREE_TO_KCAL_PER_MOL",
    "HARTREE_TO_KJ_PER_MOL",
    "StateEnergyDifferenceV1",
    "boltzmann_average",
    "boltzmann_populations",
    "cbs_extrapolation",
    "convert_energy",
    "extrapolate_correlation_helgaker",
    "extrapolate_scf_exponential",
    "state_energy_difference",
]


def extrapolate_exponential_three_point(
    energies: tuple[float, float, float],
) -> float:
    """Return the complete-basis-set limit of a three-point exponential series.

    The Hartree-Fock energy converges exponentially, ``E(n) = E_inf + A r**n``
    with ``r = exp(-alpha)``, so three energies at *equally spaced* cardinal
    numbers determine all three parameters without a nonlinear solver:

        ``E_inf = E3 - (E3 - E2)**2 / ((E3 - E2) - (E2 - E1))``

    This is the three-parameter form papers mean by "extrapolated with a
    three-parameter exponential formula".  The two-point
    :func:`extrapolate_scf_exponential` is the fixed-``alpha`` special case and
    remains available when only two basis sets exist.
    """

    if len(energies) != 3:
        raise AggregationError(
            "a three-parameter exponential extrapolation needs exactly three "
            f"energies at equally spaced cardinals, got {len(energies)}"
        )
    first, second, third = (float(item) for item in energies)
    near, far = second - first, third - second
    curvature = far - near
    if curvature == 0.0 or not math.isfinite(curvature):
        raise AggregationError(
            "the energies show no exponential curvature, so the limit is "
            "undetermined; check that they are ordered by increasing basis "
            "and come from equally spaced cardinals"
        )
    ratio = far / near if near else 0.0
    if not 0.0 < ratio < 1.0:
        raise AggregationError(
            "successive differences must shrink monotonically for an "
            f"exponential fit; observed ratio {ratio:.4f} is outside (0, 1)"
        )
    return third - (far * far) / curvature
