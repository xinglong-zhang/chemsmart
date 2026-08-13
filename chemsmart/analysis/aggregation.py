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

from ase import units as ase_units

from chemsmart.utils.constants import energy_conversion as _energy_conversion

#: Hartree to kcal/mol (CODATA-consistent with the rest of the package).
HARTREE_TO_KCAL_PER_MOL = _energy_conversion("hartree", "kcal/mol")
#: Hartree to kJ/mol.
HARTREE_TO_KJ_PER_MOL = _energy_conversion("hartree", "kJ/mol")
#: Hartree to eV.
HARTREE_TO_EV = _energy_conversion("hartree", "eV")
#: Gas constant in kcal/(mol K), for Boltzmann populations.
GAS_CONSTANT_KCAL = 1.987204258640832e-3
#: One wavenumber, as a molar energy: h * c * N_A, in kJ/mol per cm^-1.
#: This is the spectroscopic conversion, not a unit identity -- a wavenumber
#: and a molar energy have different dimensions, so the generic unit engine
#: refuses to relate them and is right to.
WAVENUMBER_TO_KJ_PER_MOL = 1.1962656362961932e-2
#: Transition-state crossover temperature per spectroscopic wavenumber.
#: ``T_c = h*c*|nu_tilde| / (2*pi*k_B)`` for ``nu_tilde`` in cm^-1.
TRANSITION_STATE_CROSSOVER_K_CM = (
    ase_units._hplanck
    * ase_units._c
    * 100.0
    / (2.0 * math.pi * ase_units._k)
)

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

    return extrapolate_correlation_inverse_power(
        smaller_cardinal=smaller_cardinal,
        larger_cardinal=larger_cardinal,
        smaller_correlation_energy=smaller_correlation_energy,
        larger_correlation_energy=larger_correlation_energy,
        exponent=3.0,
    )


def extrapolate_correlation_inverse_power(
    *,
    smaller_cardinal,
    larger_cardinal,
    smaller_correlation_energy: float,
    larger_correlation_energy: float,
    exponent: float,
) -> float:
    """Two-point inverse-power CBS extrapolation with an explicit exponent.

    This is the reusable form of the Helgaker correlation expression,
    ``E(X) = E_inf + B X**(-p)``.  The exponent is deliberately an explicit
    method/protocol input: the host must not silently assume that every
    correlated method uses the conventional ``p=3`` value.
    """

    x = _cardinal(smaller_cardinal)
    y = _cardinal(larger_cardinal)
    if x >= y:
        raise AggregationError(
            "the larger cardinal must exceed the smaller one; "
            f"got {x} and {y}"
        )
    power = float(exponent)
    if not math.isfinite(power) or power <= 0.0:
        raise AggregationError("inverse-power CBS exponent must be positive")
    xp, yp = x**power, y**power
    return (
        yp * float(larger_correlation_energy)
        - xp * float(smaller_correlation_energy)
    ) / (yp - xp)


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
    amplitude = (float(smaller_scf_energy) - float(larger_scf_energy)) / (
        fx - fy
    )
    return float(smaller_scf_energy) - amplitude * fx


def extrapolate_scf_inverse_power(
    *,
    smaller_cardinal,
    larger_cardinal,
    smaller_scf_energy: float,
    larger_scf_energy: float,
    exponent: float,
) -> float:
    """Two-point inverse-power SCF complete-basis-set limit.

    Some published protocols extrapolate the SCF component with the same
    algebraic form used for correlation, but with a protocol-specific exponent
    such as ``3.9``::

        E_inf = (X**p * E_X - Y**p * E_Y) / (X**p - Y**p)

    This is scientifically distinct from :func:`extrapolate_scf_exponential`.
    Keeping separate operations prevents a model or report from calling two
    different convergence laws equivalent merely because both produce a CBS
    energy.
    """

    return extrapolate_correlation_inverse_power(
        smaller_cardinal=smaller_cardinal,
        larger_cardinal=larger_cardinal,
        smaller_correlation_energy=smaller_scf_energy,
        larger_correlation_energy=larger_scf_energy,
        exponent=exponent,
    )


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
    degeneracies: tuple[float, ...] | None = None,
) -> tuple[float, ...]:
    """Return normalised Boltzmann populations for relative energies.

    ``energies`` may be absolute or relative; only differences matter, and the
    lowest is used as the reference to keep the exponentials finite.

    ``degeneracies`` are the multiplicities of each state, defaulting to one
    each.  They are not optional in practice: an enantiomeric pair, a rotamer
    related by an internal rotation, or a conformer reachable by symmetry all
    contribute their multiplicity, and omitting it is a scientific error rather
    than an approximation.  For n-butane at 298 K, weighting the two gauche
    forms as one state gives 82% anti where the correct treatment gives 70%.
    The values themselves belong to the protocol being reproduced; the
    weighting they enter belongs here.
    """

    if not energies:
        raise AggregationError("at least one energy is required")
    if float(temperature) <= 0:
        raise AggregationError("temperature must be positive")
    if degeneracies is None:
        multiplicities = [1.0] * len(tuple(energies))
    else:
        multiplicities = [float(item) for item in degeneracies]
        if len(multiplicities) != len(tuple(energies)):
            raise AggregationError(
                f"got {len(multiplicities)} degeneracies for "
                f"{len(tuple(energies))} states; supply one per state"
            )
        if any(item <= 0 for item in multiplicities):
            raise AggregationError("degeneracies must be positive")
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
        multiplicity
        * math.exp(
            -(item - reference) / (GAS_CONSTANT_KCAL * float(temperature))
        )
        for item, multiplicity in zip(values, multiplicities)
    ]
    total = sum(weights)
    return tuple(item / total for item in weights)


def boltzmann_average(
    values: tuple[float, ...],
    energies: tuple[float, ...],
    *,
    temperature: float = 298.15,
    unit: str = "kcal/mol",
    degeneracies: tuple[float, ...] | None = None,
) -> float:
    """Return the Boltzmann-weighted average of ``values``."""

    if len(values) != len(energies):
        raise AggregationError("values and energies must have the same length")
    populations = boltzmann_populations(
        energies,
        temperature=temperature,
        unit=unit,
        degeneracies=degeneracies,
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
    "TRANSITION_STATE_CROSSOVER_K_CM",
    "StateEnergyDifferenceV1",
    "boltzmann_average",
    "boltzmann_populations",
    "cbs_extrapolation",
    "convert_energy",
    "count_imaginary_modes",
    "extrapolate_correlation_helgaker",
    "extrapolate_correlation_inverse_power",
    "extrapolate_scf_exponential",
    "extrapolate_scf_inverse_power",
    "state_energy_difference",
    "transition_state_crossover_temperature",
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


def count_imaginary_modes(
    frequencies: tuple[float, ...],
    *,
    cutoff_cm1: float = 0.0,
) -> int:
    """Return how many harmonic modes are imaginary.

    A mode is imaginary when its frequency is negative; ChemSmart records those
    as negative wavenumbers.  ``cutoff_cm1`` ignores modes whose magnitude is
    at or below it, which is the standard treatment of numerically noisy
    near-zero modes and is what the PySCF stationary-point policy already
    applies through ``imaginary_mode_cutoff_cm1``.

    This is one of the most universal post-processing steps in the field: every
    optimization must show none, every transition state exactly one.  Left
    unexposed, a live session built a twenty-two node arithmetic contraption to
    count negative numbers in four frequency vectors, which put the definition
    of "imaginary" into the model's arithmetic instead of the toolkit's.
    """

    values = [float(item) for item in frequencies]
    if not values:
        raise AggregationError("at least one frequency is required")
    if any(not math.isfinite(item) for item in values):
        raise AggregationError("frequencies must be finite")
    cutoff = float(cutoff_cm1)
    if not math.isfinite(cutoff) or cutoff < 0:
        raise AggregationError("cutoff_cm1 must be finite and non-negative")
    return sum(1 for item in values if item < 0.0 and abs(item) > cutoff)


def transition_state_crossover_temperature(
    frequencies_cm1: tuple[float, ...],
) -> float:
    """Return the semiclassical crossover temperature of a first-order TS.

    For a full harmonic-frequency vector, exactly one negative mode is
    required and its magnitude is used.  A single scalar may instead be the
    signed imaginary frequency or an already-selected positive magnitude.
    ChemSmart owns the spectroscopic conversion and the ``2*pi`` factor, so a
    workflow need not introduce a model-authored physical constant.
    """

    values = [float(item) for item in frequencies_cm1]
    if not values:
        raise AggregationError("at least one transition frequency is required")
    if any(not math.isfinite(item) for item in values):
        raise AggregationError("transition frequencies must be finite")
    if len(values) == 1:
        magnitude = abs(values[0])
    else:
        imaginary = [item for item in values if item < 0.0]
        if len(imaginary) != 1:
            raise AggregationError(
                "a transition-state frequency vector must contain exactly "
                f"one imaginary mode; observed {len(imaginary)}"
            )
        magnitude = abs(imaginary[0])
    if magnitude == 0.0:
        raise AggregationError("the transition frequency must be non-zero")
    return magnitude * TRANSITION_STATE_CROSSOVER_K_CM


def harmonic_zero_point_energy(
    frequencies: tuple[float, ...],
    *,
    unit: str = "kj/mol",
) -> float:
    """Return the harmonic zero-point vibrational energy of a mode set.

    ``E_ZPVE = 1/2 * sum(h * c * nu_tilde) * N_A`` over the real modes, with
    frequencies given as wavenumbers in cm^-1.  ChemSmart already owns this as
    ``Thermochemistry.zero_point_energy``, which reaches it through vibrational
    temperatures; this is the same quantity taken straight from a frequency
    vector, which is the form a post-processing DAG holds.

    Exposing it removes two model-authored numbers from every ZPVE plan.  A
    live session asked for a scaled ZPVE built ``sum -> scale 0.5 -> convert ->
    scale`` and typed the 1/2 itself, then asked the generic unit engine to
    turn cm^-1 into kJ/mol.  That request is refused, correctly: a wavenumber
    and a molar energy are different dimensions, and relating them is a
    spectroscopic convention rather than a unit conversion.  The session
    predicted the refusal and proposed to hand-multiply by ``h*c*N_A`` instead,
    which would have put a physical constant into model arithmetic.

    Imaginary modes are excluded rather than added as negative contributions:
    an imaginary mode has no zero-point level.  Use
    :func:`count_imaginary_modes` to decide whether the structure should have
    been used at all.
    """

    values = [float(item) for item in frequencies]
    if not values:
        raise AggregationError("at least one frequency is required")
    if any(not math.isfinite(item) for item in values):
        raise AggregationError("frequencies must be finite")
    real_modes = [item for item in values if item > 0.0]
    if not real_modes:
        raise AggregationError(
            "a zero-point energy needs at least one real vibrational mode; "
            "every supplied frequency is imaginary or zero"
        )
    total_kj = 0.5 * sum(real_modes) * WAVENUMBER_TO_KJ_PER_MOL
    key = str(unit).strip().lower()
    if key not in ENERGY_UNIT_FACTORS:
        raise AggregationError(
            f"unsupported energy unit {unit!r}; "
            f"expected one of {sorted(ENERGY_UNIT_FACTORS)}"
        )
    return convert_energy(total_kj / HARTREE_TO_KJ_PER_MOL, key)
