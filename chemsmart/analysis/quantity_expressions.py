"""A bounded dimensional expression evaluator for parsed result quantities.

The evaluator intentionally accepts a typed directed acyclic graph rather than
Python, a formula string, or a shell fragment.  Nodes are evaluated once in
declared order, every number must be finite, and arithmetic is rejected when
dimensions are incompatible.  It is suitable for geometry measurements,
thermochemical identities, energy differences, and simple statistical
reductions without becoming a benchmark-specific calculator.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass, replace
from typing import Any, Iterable, Mapping

import numpy as np
from ase import units as ase_units

from chemsmart.analysis.literature_constants import (
    UnknownLiteratureConstantError,
    literature_constant,
)
from chemsmart.analysis.result_quantities import (
    ANGLE,
    CHARGE,
    DIMENSIONLESS,
    DIPOLE_MOMENT,
    ELECTRIC_POTENTIAL,
    ENERGY,
    ENTROPY,
    FREQUENCY,
    IR_INTENSITY,
    LENGTH,
    MASS,
    MOMENT_OF_INERTIA,
    PRESSURE,
    TEMPERATURE,
    Dimension,
    QuantityContractError,
    QuantityValueV1,
    canonical_quantity_sha256,
    make_quantity_value,
)
from chemsmart.utils.constants import au_to_debye, energy_conversion
from chemsmart.utils.geometry import (
    internal_angle,
    internal_dihedral,
    internal_distance,
)

MAX_EXPRESSION_NODES = 128
MAX_NODE_INPUTS = 64
_RECEIPT_REF = re.compile(r"(?:^|[;])receipt:([0-9a-f]{64})(?:$|[;])")
_QUANTITY_REF = re.compile(r"(?:^|[;])quantity:([^;]+)(?:$|[;])")
_SEMANTIC_ROLE_REF = re.compile(r"(?:^|[;])semantic-role:([^;]+)(?:$|[;])")

_OPERATIONS = frozenset(
    {
        "ref",
        "literal",
        "constant",
        "add",
        "subtract",
        "multiply",
        "divide",
        "scale",
        "abs",
        "sqrt",
        "power",
        "exp",
        "log",
        "sum",
        "mean",
        "min",
        "max",
        "coordinate_at_maximum",
        "coordinate_at_minimum",
        "distance",
        "angle",
        "dihedral",
        "convert",
        "linear_fit_slope",
        "linear_fit_intercept",
        "exponential_cbs_limit",
        "scf_exponential_cbs_limit",
        "scf_inverse_power_cbs_limit",
        "correlation_inverse_power_cbs_limit",
        "photon_wavelength",
        "gibbs_to_pka",
        "gibbs_to_redox_potential",
        "boltzmann_populations",
        "boltzmann_average",
        "imaginary_mode_count",
        "harmonic_zero_point_energy",
        "wavenumber_to_energy",
        "energy_to_wavenumber",
        "transition_state_crossover_temperature",
        "center_of_mass",
        "principal_moments_of_inertia",
        "linear_rotor_constant",
        "rigid_rotor_constants",
        "connectivity_difference_count",
    }
)


#: What each operation computes, and the input shape it expects.
#:
#: A bare enum of operation names asks a model to infer both.  Observed live:
#: a session that needed a three-point exponential basis-set limit rebuilt the
#: closed form from fifteen multiply/subtract/scale/divide nodes rather than
#: call the one operation that owns it, because nothing said the operation was
#: the right instrument or what it wanted.  Reconstructing a domain convention
#: by hand is the failure mode this vocabulary exists to prevent, so the
#: descriptions are part of the contract and are pinned to the operation set.
OPERATION_DESCRIPTIONS: Mapping[str, str] = {
    "ref": "name an earlier value or expression input; use indices to select",
    "literal": "a constant you supply; recorded as model-authored",
    "constant": (
        "a host-owned literature value selected by registered name via "
        "constant_name; no inputs; its unit and standard-state convention "
        "are owned by the registry, never restated by hand"
    ),
    "add": "sum of two values of one dimension",
    "subtract": "first input minus the second, one dimension",
    "multiply": "product of two values; dimensions multiply",
    "divide": "first input over the second; dimensions divide",
    "scale": "multiply one input by scale_factor; recorded as model-authored",
    "abs": "absolute value",
    "sqrt": "square root; the dimension must be an even power",
    "power": "raise one input to literal_value; recorded as model-authored",
    "exp": "exponential of a dimensionless input",
    "log": "natural logarithm of a positive dimensionless input",
    "sum": "sum over inputs or over one vector input",
    "mean": "arithmetic mean over inputs or over one vector input",
    "min": "smallest of the inputs",
    "max": "largest of the inputs",
    "coordinate_at_maximum": (
        "where along a series its largest value falls. Takes exactly two "
        "equal-length vectors in order: the values being extremised, then the "
        "coordinate they were measured at. Returns that coordinate, carrying "
        "the coordinate's own dimension. For a scanned surface this is the "
        "position of the barrier, which max alone cannot give"
    ),
    "coordinate_at_minimum": (
        "where along a series its smallest value falls; same two ordered "
        "vectors as coordinate_at_maximum. For a scanned surface this is the "
        "position of the well"
    ),
    "distance": "distance between two indexed coordinate vectors",
    "angle": "angle at the middle of three indexed coordinate vectors",
    "dihedral": (
        "signed torsion about the middle bond of four indexed coordinate "
        "vectors in bonded order a-b-c-d, in (-180, 180]. The third "
        "standard internal coordinate alongside distance and angle; do not "
        "rebuild it from cross products"
    ),
    "convert": "restate a value in target_unit; arithmetic stays canonical",
    "linear_fit_slope": "slope of a least-squares line through x and y",
    "linear_fit_intercept": "intercept of that same least-squares line",
    "exponential_cbs_limit": (
        "complete-basis-set limit of an exponentially convergent series such "
        "as Hartree-Fock. Takes the energies at three equally spaced cardinal "
        "numbers ordered by increasing basis, as three scalar inputs or one "
        "three-element input, and fits the decay from the data. Introduces no "
        "constant of your own; prefer it whenever a protocol says the energy "
        "was extrapolated exponentially and you have three points"
    ),
    "scf_exponential_cbs_limit": (
        "two-point SCF exponential limit from two energies at "
        "cardinal_numbers, using an extrapolation_exponent you supply"
    ),
    "scf_inverse_power_cbs_limit": (
        "two-point SCF complete-basis-set limit in inverse powers of the "
        "cardinal number, using an explicit extrapolation_exponent. This is "
        "a different convergence law from scf_exponential_cbs_limit"
    ),
    "correlation_inverse_power_cbs_limit": (
        "two-point correlation-energy limit in inverse powers of the cardinal "
        "number, using an extrapolation_exponent you supply. Takes correlation "
        "energies, not total energies"
    ),
    "photon_wavelength": "wavelength of a positive excitation energy",
    "gibbs_to_pka": (
        "pKa from a deprotonation free energy and a temperature, in that "
        "order; owns pKa = deltaG / (RT ln 10)"
    ),
    "gibbs_to_redox_potential": (
        "electrode potential from a half-reaction free energy and the number "
        "of electrons transferred, in that order; owns E = -deltaG / (n F) "
        "and with it the IUPAC sign, so a reduction free energy that is "
        "negative gives a positive potential. The result is an absolute "
        "potential on the same scale as the free energy it came from; "
        "subtract a reference-electrode constant separately, so which "
        "electrode you referenced stays visible in the expression"
    ),
    "boltzmann_populations": (
        "normalized Boltzmann populations of a set of states. Takes the state "
        "energies as one vector or as separate scalar inputs, followed by a "
        "temperature and optionally a dimensionless vector of per-state "
        "degeneracies, and owns the "
        "weighting, the gas constant and the unit handling. Supply the "
        "degeneracies whenever states are multiply realizable -- an "
        "enantiomeric pair counts twice. Do not rebuild any of this from exp, "
        "scale, divide and sum"
    ),
    "boltzmann_average": (
        "Boltzmann-weighted **linear** average of a property. Takes the "
        "per-state values, the state energies and a temperature, in that "
        "order, optionally followed by the per-state degeneracies. Supply the "
        "degeneracies whenever states are multiply realizable, exactly as for "
        "boltzmann_populations. The linear mean is the ensemble value only for "
        "a property that averages linearly. A quantity reported as the "
        "magnitude of a vector -- a dipole moment is the common case -- does "
        "not: its measured ensemble value is the root-mean-square average, so "
        "average the squares and take the square root instead of averaging the "
        "magnitudes"
    ),
    "imaginary_mode_count": (
        "how many harmonic modes are imaginary, from one frequency vector, "
        "optionally followed by a frequency cutoff below which a near-zero "
        "mode is ignored. Owns the sign convention: a minimum has zero, a "
        "transition state exactly one. Do not rebuild it from comparisons. "
        "The cutoff defaults to zero, so every negative counts; ChemSmart's "
        "thermochemistry instead treats a mode within 20 cm^-1 of zero as "
        "numerical noise in a translation or rotation rather than a reaction "
        "coordinate. A floppy torsion landing at, say, -12 cm^-1 therefore "
        "counts as imaginary here while the same result derives "
        "thermochemistry as a minimum. Pass 20 cm^-1 to ask the question the "
        "thermochemistry answers, or leave the default to ask the strict "
        "one -- but say which you asked, because the two disagree exactly "
        "where a structure is nearly, not quite, relaxed"
    ),
    "wavenumber_to_energy": (
        "restate a wavenumber (cm^-1) as a molar energy through h*c*N_A, "
        "scalar or vector; target_unit names the energy unit (default "
        "kJ/mol). Exchange couplings, zero-field splittings and spin-orbit "
        "gaps are quoted in cm^-1 and computed as energy differences; this "
        "operation owns the factor, so a declaration in cm^-1 is answered "
        "by a claim in cm^-1 without inventing a conversion"
    ),
    "energy_to_wavenumber": (
        "restate a molar energy as a wavenumber in cm^-1 through h*c*N_A, "
        "scalar or vector; the inverse of wavenumber_to_energy"
    ),
    "harmonic_zero_point_energy": (
        "harmonic zero-point vibrational energy from one frequency vector in "
        "cm^-1, returned as a molar energy. Owns both the factor of one half "
        "and the spectroscopic h*c*N_A conversion, and skips imaginary modes. "
        "Use this instead of sum -> scale 0.5 -> convert: a wavenumber and a "
        "molar energy are different dimensions, so convert refuses that step"
    ),
    "transition_state_crossover_temperature": (
        "semiclassical crossover temperature from one transition-state "
        "frequency vector, or from one already-selected imaginary-frequency "
        "magnitude. Owns selection of the unique imaginary mode and the "
        "h*c/(2*pi*k_B) spectroscopic conversion; do not supply that physical "
        "constant as a literal"
    ),
    "center_of_mass": (
        "mass-weighted Cartesian center from a coordinate matrix followed by "
        "an equal-length atomic-mass vector"
    ),
    "principal_moments_of_inertia": (
        "three ascending principal moments for a coordinate matrix and an "
        "equal-length atomic-mass vector; owns COM translation, construction "
        "of the full inertia tensor and symmetric diagonalization"
    ),
    "linear_rotor_constant": (
        "the single finite rotational constant B in cm^-1 from the three "
        "principal moments of a linear molecule. Requires one near-zero axial "
        "moment and two equal positive perpendicular moments, and owns "
        "h/(8*pi^2*c*I) and the spectroscopic unit conversion"
    ),
    "rigid_rotor_constants": (
        "A, B, C in descending frequency order from three positive principal "
        "moments of a nonlinear molecule; owns h/(8*pi^2*I) and the "
        "spectroscopic unit conversion"
    ),
    "connectivity_difference_count": (
        "number of undirected edges present in exactly one of two "
        "geometry-perceived connectivity matrices. Takes four inputs ordered "
        "as first connectivity, first symbols, second connectivity, second "
        "symbols; requires identical atom identity and order, square symmetric "
        "binary matrices, and zero diagonals. Connectivity is inferred from "
        "geometry and covalent radii, not an electronic bond-order claim"
    ),
}

if set(OPERATION_DESCRIPTIONS) != set(_OPERATIONS):  # pragma: no cover
    raise RuntimeError(
        "every expression operation must be described: "
        f"{sorted(set(_OPERATIONS) ^ set(OPERATION_DESCRIPTIONS))}"
    )

#: How many inputs each operation accepts, stated where the operations live.
#:
#: A frozenset lists the exact admissible counts; a (minimum, None) tuple
#: means that many or more.  The evaluator has always enforced these counts,
#: but only while evaluating -- after every engine had finished.  A live
#: six-species pKa workflow wrote an isodesmic exchange as one four-input
#: `subtract`, passed planning, preview, and approval, ran six solvated
#: opt+freq jobs, and lost the whole payload to "subtract requires two
#: inputs".  The count is a structural fact about the node, fully
#: determinable at planning, so it is declared here and checked at node
#: construction; the guard below pins the table to the operation set so the
#: two cannot drift apart.
OPERATION_INPUT_COUNTS: Mapping[str, frozenset[int] | tuple[int, None]] = {
    "ref": frozenset({0, 1}),
    "literal": frozenset({0}),
    "constant": frozenset({0}),
    "add": frozenset({2}),
    "subtract": frozenset({2}),
    "multiply": frozenset({2}),
    "divide": frozenset({2}),
    "scale": frozenset({1}),
    "abs": frozenset({1}),
    "sqrt": frozenset({1}),
    "power": frozenset({1}),
    "exp": frozenset({1}),
    "log": frozenset({1}),
    "sum": (1, None),
    "mean": (1, None),
    "min": (1, None),
    "max": (1, None),
    "coordinate_at_maximum": frozenset({2}),
    "coordinate_at_minimum": frozenset({2}),
    "distance": frozenset({2}),
    "angle": frozenset({3}),
    "dihedral": frozenset({4}),
    "convert": frozenset({1}),
    "linear_fit_slope": frozenset({2}),
    "linear_fit_intercept": frozenset({2}),
    "exponential_cbs_limit": frozenset({1, 3}),
    "scf_exponential_cbs_limit": frozenset({2}),
    "scf_inverse_power_cbs_limit": frozenset({2}),
    "correlation_inverse_power_cbs_limit": frozenset({2}),
    "photon_wavelength": frozenset({1}),
    "gibbs_to_pka": frozenset({2}),
    "gibbs_to_redox_potential": frozenset({2}),
    "boltzmann_populations": (2, None),
    "boltzmann_average": (3, None),
    "imaginary_mode_count": frozenset({1, 2}),
    "harmonic_zero_point_energy": frozenset({1}),
    "wavenumber_to_energy": frozenset({1}),
    "energy_to_wavenumber": frozenset({1}),
    "transition_state_crossover_temperature": frozenset({1}),
    "center_of_mass": frozenset({2}),
    "principal_moments_of_inertia": frozenset({2}),
    "linear_rotor_constant": frozenset({1}),
    "rigid_rotor_constants": frozenset({1}),
    "connectivity_difference_count": frozenset({4}),
}

if set(OPERATION_INPUT_COUNTS) != set(_OPERATIONS):  # pragma: no cover
    raise RuntimeError(
        "every expression operation must state its input count: "
        f"{sorted(set(_OPERATIONS) ^ set(OPERATION_INPUT_COUNTS))}"
    )


def _admissible_count_phrase(
    admissible: frozenset[int] | tuple[int, None],
) -> str:
    """State an operation's admissible input counts in one readable clause."""
    if isinstance(admissible, tuple):
        return f"{admissible[0]} or more inputs"
    counts = sorted(admissible)
    if counts == [0]:
        return "no inputs"
    words = " or ".join(str(item) for item in counts)
    return f"exactly {words} input(s)"


def require_expression_input_count(operation: str, count: int) -> None:
    """Refuse an input count the named operation can never evaluate.

    Raises the expression vocabulary's own error, carrying the
    operation's description: the refusal is the moment a session learns
    what the operation wanted, so the message states the accepted
    shape rather than only the count.
    """
    admissible = OPERATION_INPUT_COUNTS[operation]
    if isinstance(admissible, tuple):
        if count >= admissible[0]:
            return
    elif count in admissible:
        return
    noun = "input" if count == 1 else "inputs"
    raise QuantityExpressionError(
        f"{operation} accepts {_admissible_count_phrase(admissible)}; "
        f"got {count} {noun} -- {operation}: "
        f"{OPERATION_DESCRIPTIONS[operation]}"
    )


#: The family a chemist would look for each operation under.
#:
#: Until now only fifteen operations carried a family, and they carried it
#: indirectly: ``guides.LEAF_OPERATIONS`` said which guide hid them, so
#: "family" meant "what the stem does not show" and the other twenty-nine
#: had none.  That is a surface fact, not a scientific one.  This map is
#: the scientific one: it is what the searchable catalogue groups by and
#: what the capability ladder reports, and it exists so an operation can be
#: *found* by the quantity it computes rather than by the tool tree it
#: happened to sit in.  Every operation is a member of exactly one family;
#: ``test_every_operation_declares_one_family`` pins the map to the
#: operation set, so a new operation cannot arrive unfamilied.
OPERATION_FAMILIES: Mapping[str, str] = {
    # Algebra and plumbing: the dimension-aware arithmetic every expression
    # is built out of, which owns no chemistry convention of its own.
    "abs": "arithmetic",
    "add": "arithmetic",
    "convert": "arithmetic",
    "divide": "arithmetic",
    "exp": "arithmetic",
    "literal": "arithmetic",
    "log": "arithmetic",
    "max": "arithmetic",
    "mean": "arithmetic",
    "min": "arithmetic",
    "multiply": "arithmetic",
    "power": "arithmetic",
    "ref": "arithmetic",
    "scale": "arithmetic",
    "sqrt": "arithmetic",
    "subtract": "arithmetic",
    "sum": "arithmetic",
    # Molecular geometry read off a coordinate matrix.
    "angle": "geometry",
    "center_of_mass": "geometry",
    "connectivity_difference_count": "geometry",
    "dihedral": "geometry",
    "distance": "geometry",
    # Where a scanned or fitted series turns over, and the line through it.
    "coordinate_at_maximum": "series",
    "coordinate_at_minimum": "series",
    "linear_fit_intercept": "series",
    "linear_fit_slope": "series",
    # Moments of inertia and the rotational constants built from them.
    "linear_rotor_constant": "rotational",
    "principal_moments_of_inertia": "rotational",
    "rigid_rotor_constants": "rotational",
    # Harmonic frequencies read as chemistry.
    "harmonic_zero_point_energy": "vibrational",
    "imaginary_mode_count": "vibrational",
    "transition_state_crossover_temperature": "vibrational",
    # Energy and wavelength as one another.
    "energy_to_wavenumber": "spectroscopy",
    "photon_wavelength": "spectroscopy",
    "wavenumber_to_energy": "spectroscopy",
    # Populations over a set of states at a temperature.
    "boltzmann_average": "ensemble",
    "boltzmann_populations": "ensemble",
    # Complete-basis-set extrapolation.
    "correlation_inverse_power_cbs_limit": "cbs",
    "exponential_cbs_limit": "cbs",
    "scf_exponential_cbs_limit": "cbs",
    "scf_inverse_power_cbs_limit": "cbs",
    # A registered literature value, and the two conversions whose
    # constants and standard states the registry owns.
    "constant": "constants",
    "gibbs_to_pka": "constants",
    "gibbs_to_redox_potential": "constants",
}


#: Operations that carry a computational-chemistry convention ChemSmart owns.
#: Reaching a reported quantity through one of these means the convention came
#: from the toolkit.  Reaching it through arithmetic instead means the model
#: supplied the convention, which is the situation the project exists to avoid
#: and which no per-paper check would catch in general.
CONVENTION_OPERATIONS = frozenset(
    {
        "angle",
        "dihedral",
        "boltzmann_average",
        "boltzmann_populations",
        "correlation_inverse_power_cbs_limit",
        "harmonic_zero_point_energy",
        "wavenumber_to_energy",
        "energy_to_wavenumber",
        "imaginary_mode_count",
        "transition_state_crossover_temperature",
        "distance",
        "center_of_mass",
        "gibbs_to_pka",
        "gibbs_to_redox_potential",
        "principal_moments_of_inertia",
        "linear_rotor_constant",
        "rigid_rotor_constants",
        "connectivity_difference_count",
        "exponential_cbs_limit",
        "linear_fit_intercept",
        "linear_fit_slope",
        "photon_wavelength",
        "scf_exponential_cbs_limit",
        "scf_inverse_power_cbs_limit",
    }
)

#: Operations that move or restate a value without computing anything with it.
#: They are neither a convention nor arithmetic, so they are counted as
#: neither.
_PLUMBING_OPERATIONS = frozenset({"ref", "literal", "constant", "convert"})

#: Everything else: general arithmetic and reductions.  A reported quantity
#: reached through many of these, and through no convention operation, was
#: assembled rather than computed by the toolkit.
ARITHMETIC_OPERATIONS = (
    frozenset(_OPERATIONS) - CONVENTION_OPERATIONS - (_PLUMBING_OPERATIONS)
)


class QuantityExpressionError(ValueError):
    """Raised when an expression is unsafe, ill-typed, or non-finite."""


def _canonical_dimension(dimension: Dimension) -> Dimension:
    """Drop only cancelled appended bases, preserving the six-base legacy."""

    values = list(dimension)
    while len(values) > len(DIMENSIONLESS) and values[-1] == 0:
        values.pop()
    return tuple(values)


def _add_dimensions(left: Dimension, right: Dimension) -> Dimension:
    size = max(len(left), len(right))
    return _canonical_dimension(
        tuple(
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        )
    )


def _subtract_dimensions(left: Dimension, right: Dimension) -> Dimension:
    size = max(len(left), len(right))
    return _canonical_dimension(
        tuple(
            (left[index] if index < len(left) else 0)
            - (right[index] if index < len(right) else 0)
            for index in range(size)
        )
    )


def canonical_unit_for_dimension(dimension: Dimension) -> str:
    dimension = _canonical_dimension(dimension)
    known = {
        DIMENSIONLESS: "1",
        ENERGY: "hartree",
        LENGTH: "angstrom",
        TEMPERATURE: "K",
        ANGLE: "radian",
        FREQUENCY: "cm^-1",
        PRESSURE: "atm",
        ENTROPY: "hartree K^-1",
        DIPOLE_MOMENT: "debye",
        MASS: "u",
        CHARGE: "e",
        ELECTRIC_POTENTIAL: "hartree e^-1",
        IR_INTENSITY: "km/mol",
    }
    if dimension in known:
        return known[dimension]
    labels = (
        "hartree",
        "angstrom",
        "K",
        "radian",
        "cm^-1",
        "atm",
        "debye",
        "u",
        "e",
        "km/mol",
    )
    terms = []
    for label, exponent in zip(labels, dimension):
        if exponent == 0:
            continue
        terms.append(label if exponent == 1 else f"{label}^{exponent}")
    return " ".join(terms) if terms else "1"


def _normalized_unit_key(unit: str) -> str:
    normalized = (
        str(unit)
        .strip()
        .lower()
        .replace("·", " ")
        .replace("*", " ")
        .replace("å", "angstrom")
        .replace("−", "-")
        # xTB labels its native IR column ``km·mol⁻¹``.  This is the only
        # superscript spelling that reaches the public unit vocabulary here;
        # do not broaden normalization for unrelated Unicode numerals.
        .replace("⁻¹", "^-1")
    )
    normalized = " ".join(normalized.split())
    # A bracketed denominator is the SI-recommended spelling of a compound
    # unit, and it is how a chemist writes molar entropy: "J/(mol K)" is the
    # same unit as the already supported "J/mol/K".  Spacing around the
    # solidus and inside the bracket carries no meaning, so normalise it away
    # and distribute the solidus over the bracketed factors.  That keeps one
    # vocabulary for both spellings instead of duplicating every alias.
    normalized = re.sub(r"\s*/\s*", "/", normalized)
    normalized = re.sub(r"\(\s*", "(", normalized)
    normalized = re.sub(r"\s*\)", ")", normalized)
    return re.sub(
        r"/\(([^()]+)\)",
        lambda match: "/" + "/".join(match.group(1).split()),
        normalized,
    )


def _compound_unit_spec(
    key: str,
    aliases: Mapping[str, tuple[Dimension, str, float]],
) -> tuple[Dimension, str, float]:
    """Parse products of supported one-token units and integer powers.

    Arithmetic already emits canonical products such as ``angstrom^2 u``.
    Accepting the same vocabulary on input makes those values convertible and
    claimable without introducing a formula language.  Exact historical
    aliases (including slash and molar forms) are resolved before this helper.
    """

    token_aliases = {
        name: spec
        for name, spec in aliases.items()
        if name
        and " " not in name
        and "/" not in name
        and spec[0] != DIMENSIONLESS
    }
    bases = sorted(token_aliases, key=len, reverse=True)
    dimension = DIMENSIONLESS
    factor = 1.0
    for token in key.split():
        exponent = 1
        base = token if token in token_aliases else ""
        if not base:
            for candidate in bases:
                prefix = f"{candidate}^"
                if not token.startswith(prefix):
                    continue
                suffix = token[len(prefix) :]
                try:
                    exponent = int(suffix)
                except ValueError:
                    continue
                base = candidate
                break
        if not base:
            raise QuantityExpressionError(f"unsupported unit: {key!r}")
        base_dimension, _, base_factor = token_aliases[base]
        scaled_dimension = tuple(value * exponent for value in base_dimension)
        dimension = _add_dimensions(dimension, scaled_dimension)
        factor *= base_factor**exponent
    if not key:
        raise QuantityExpressionError(f"unsupported unit: {key!r}")
    return dimension, canonical_unit_for_dimension(dimension), factor


def unit_dimension(unit: str) -> Dimension:
    """Public dimension of a unit string, for plan-time checking."""

    dimension, _, _ = _unit_spec(unit)
    return dimension


#: Operations whose result carries their (single, agreeing) input
#: dimension. A disagreement among known input dimensions here is a
#: provable planning error, not an unknown.
_SAME_DIMENSION_OPERATIONS = frozenset(
    {
        "add",
        "subtract",
        "min",
        "max",
        "sum",
        "mean",
        "abs",
        "exponential_cbs_limit",
        "scf_exponential_cbs_limit",
        "scf_inverse_power_cbs_limit",
        "correlation_inverse_power_cbs_limit",
    }
)

#: Operations whose result dimension is fixed by the operation itself.
#: Stated via unit strings so the table cannot drift from the unit
#: registry. Deliberately incomplete: an operation absent here and
#: from every other rule yields "unknown", never a guess.
_FIXED_DIMENSION_OPERATION_UNITS = {
    "distance": "angstrom",
    "angle": "degree",
    "dihedral": "degree",
    "photon_wavelength": "angstrom",
    "gibbs_to_pka": "1",
    "boltzmann_populations": "1",
    "imaginary_mode_count": "1",
    "connectivity_difference_count": "1",
    "harmonic_zero_point_energy": "hartree",
    "wavenumber_to_energy": "hartree",
    "energy_to_wavenumber": "cm^-1",
    "transition_state_crossover_temperature": "K",
}


def operation_result_dimension(
    operation: str,
    input_dimensions: tuple[Dimension | None, ...],
    *,
    literal_unit: str = "",
    target_unit: str = "",
    constant_name: str = "",
) -> Dimension | None:
    """Dimension an expression node must produce, derivable at plan time.

    Returns ``None`` when the rule is unknown or an input dimension is
    unknown -- the caller skips rather than guesses. Raises
    QuantityExpressionError only on a provable inconsistency (a
    same-dimension operation fed unlike dimensions; a convert whose
    target dimension its input cannot reach). Motivated live by a
    claim render dying on a dimension mismatch after every engine had
    finished -- the dimensional twin of the four-input-subtract loss
    this file's arity table already refuses. The walk trusts declared
    producer units as seeds, so declared-versus-actual extraction lies
    are the extraction-side gate's job, not this rule's.
    """

    known = [item for item in input_dimensions if item is not None]
    if operation in {"ref", "scale"}:
        return input_dimensions[0] if input_dimensions else None
    if operation == "literal":
        return unit_dimension(literal_unit) if literal_unit else None
    if operation == "constant":
        if not constant_name:
            return None
        try:
            from chemsmart.analysis.literature_constants import (
                literature_constant,
            )

            return unit_dimension(literature_constant(constant_name).unit)
        except Exception:
            return None
    if operation in _SAME_DIMENSION_OPERATIONS:
        if len(known) != len(input_dimensions) or not known:
            return None
        first = _canonical_dimension(known[0])
        for item in known[1:]:
            if _canonical_dimension(item) != first:
                raise QuantityExpressionError(
                    f"{operation} requires inputs of one dimension; the "
                    "plan's declared units give "
                    f"{canonical_unit_for_dimension(first)!r} and "
                    f"{canonical_unit_for_dimension(item)!r}"
                )
        return first
    if operation in {"multiply", "divide"}:
        if len(known) != 2 or len(input_dimensions) != 2:
            return None
        first = _canonical_dimension(known[0])
        second = _canonical_dimension(known[1])
        if operation == "divide":
            second = tuple(-value for value in second)
        return _add_dimensions(first, second)
    if operation == "convert":
        if not input_dimensions or input_dimensions[0] is None:
            return None
        source = _canonical_dimension(input_dimensions[0])
        if target_unit:
            target = _canonical_dimension(unit_dimension(target_unit))
            if target != source:
                raise QuantityExpressionError(
                    f"convert to {target_unit!r} is unreachable from a "
                    "quantity of dimension "
                    f"{canonical_unit_for_dimension(source)!r}"
                )
            return target
        return source
    if operation in {"coordinate_at_minimum", "coordinate_at_maximum"}:
        # The answer is a coordinate: it carries the second input's
        # dimension, exactly as the evaluator computes it.
        if len(input_dimensions) == 2:
            return input_dimensions[1]
        return None
    if operation in _FIXED_DIMENSION_OPERATION_UNITS:
        return unit_dimension(_FIXED_DIMENSION_OPERATION_UNITS[operation])
    return None


def _unit_spec(unit: str) -> tuple[Dimension, str, float]:
    key = _normalized_unit_key(unit)
    aliases: dict[str, tuple[Dimension, str, float]] = {
        "": (DIMENSIONLESS, "1", 1.0),
        "1": (DIMENSIONLESS, "1", 1.0),
        "dimensionless": (DIMENSIONLESS, "1", 1.0),
        "hartree": (ENERGY, "hartree", 1.0),
        "ha": (ENERGY, "hartree", 1.0),
        "eh": (ENERGY, "hartree", 1.0),
        "ev": (ENERGY, "hartree", energy_conversion("eV", "hartree", 1.0)),
        "j/mol": (
            ENERGY,
            "hartree",
            energy_conversion("J/mol", "hartree", 1.0),
        ),
        "j mol^-1": (
            ENERGY,
            "hartree",
            energy_conversion("J/mol", "hartree", 1.0),
        ),
        "kj/mol": (
            ENERGY,
            "hartree",
            energy_conversion("kJ/mol", "hartree", 1.0),
        ),
        "kj mol^-1": (
            ENERGY,
            "hartree",
            energy_conversion("kJ/mol", "hartree", 1.0),
        ),
        "kcal/mol": (
            ENERGY,
            "hartree",
            energy_conversion("kcal/mol", "hartree", 1.0),
        ),
        "kcal mol^-1": (
            ENERGY,
            "hartree",
            energy_conversion("kcal/mol", "hartree", 1.0),
        ),
        "hartree/k": (ENTROPY, "hartree K^-1", 1.0),
        "hartree k^-1": (ENTROPY, "hartree K^-1", 1.0),
        "j/mol/k": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("J/mol", "hartree", 1.0),
        ),
        "j mol^-1 k^-1": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("J/mol", "hartree", 1.0),
        ),
        "kj/mol/k": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kJ/mol", "hartree", 1.0),
        ),
        "kj mol^-1 k^-1": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kJ/mol", "hartree", 1.0),
        ),
        "cal/mol/k": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kcal/mol", "hartree", 0.001),
        ),
        "cal mol^-1 k^-1": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kcal/mol", "hartree", 0.001),
        ),
        "kcal/mol/k": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kcal/mol", "hartree", 1.0),
        ),
        "kcal mol^-1 k^-1": (
            ENTROPY,
            "hartree K^-1",
            energy_conversion("kcal/mol", "hartree", 1.0),
        ),
        "angstrom": (LENGTH, "angstrom", 1.0),
        "ang": (LENGTH, "angstrom", 1.0),
        "bohr": (LENGTH, "angstrom", 0.529177210903),
        "a0": (LENGTH, "angstrom", 0.529177210903),
        "nm": (LENGTH, "angstrom", 10.0),
        "debye": (DIPOLE_MOMENT, "debye", 1.0),
        "d": (DIPOLE_MOMENT, "debye", 1.0),
        "e bohr": (DIPOLE_MOMENT, "debye", au_to_debye),
        "e a0": (DIPOLE_MOMENT, "debye", au_to_debye),
        "u": (MASS, "u", 1.0),
        "da": (MASS, "u", 1.0),
        "dalton": (MASS, "u", 1.0),
        "amu": (MASS, "u", 1.0),
        "k": (TEMPERATURE, "K", 1.0),
        "kelvin": (TEMPERATURE, "K", 1.0),
        "radian": (ANGLE, "radian", 1.0),
        "rad": (ANGLE, "radian", 1.0),
        "degree": (ANGLE, "radian", math.pi / 180.0),
        "degrees": (ANGLE, "radian", math.pi / 180.0),
        "deg": (ANGLE, "radian", math.pi / 180.0),
        "cm^-1": (FREQUENCY, "cm^-1", 1.0),
        "cm-1": (FREQUENCY, "cm^-1", 1.0),
        "hz": (FREQUENCY, "cm^-1", 1.0 / 29_979_245_800.0),
        "mhz": (FREQUENCY, "cm^-1", 1.0e6 / 29_979_245_800.0),
        "ghz": (FREQUENCY, "cm^-1", 1.0e9 / 29_979_245_800.0),
        "thz": (FREQUENCY, "cm^-1", 1.0e12 / 29_979_245_800.0),
        "atm": (PRESSURE, "atm", 1.0),
        "bar": (PRESSURE, "atm", 1.0 / 1.01325),
        "pa": (PRESSURE, "atm", 1.0 / 101_325.0),
        "e": (CHARGE, "e", 1.0),
        # One volt acting on one elementary charge is one electronvolt, by
        # definition of the electronvolt.  So the volt conversion factor is
        # exactly the eV-to-hartree factor already in this table, and no new
        # physical constant enters here.
        "v": (
            ELECTRIC_POTENTIAL,
            "hartree e^-1",
            energy_conversion("eV", "hartree", 1.0),
        ),
        "volt": (
            ELECTRIC_POTENTIAL,
            "hartree e^-1",
            energy_conversion("eV", "hartree", 1.0),
        ),
        "mv": (
            ELECTRIC_POTENTIAL,
            "hartree e^-1",
            energy_conversion("eV", "hartree", 1.0e-3),
        ),
        "km/mol": (IR_INTENSITY, "km/mol", 1.0),
        "km mol^-1": (IR_INTENSITY, "km/mol", 1.0),
        "km mol-1": (IR_INTENSITY, "km/mol", 1.0),
    }
    try:
        return aliases[key]
    except KeyError:
        try:
            return _compound_unit_spec(key, aliases)
        except QuantityExpressionError as exc:
            raise QuantityExpressionError(
                f"unsupported unit: {unit!r}"
            ) from exc


def normalize_numeric_value(
    value: Any, unit: str
) -> tuple[Any, str, Dimension]:
    """Convert a finite numeric value to the canonical unit for its dimension."""

    dimension, canonical_unit, factor = _unit_spec(unit)
    array = np.asarray(value, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise QuantityExpressionError(
            "numeric values must be finite and non-empty"
        )
    normalized = array * factor
    if normalized.ndim == 0:
        payload: Any = float(normalized)
    elif normalized.ndim <= 2:
        payload = tuple(
            (
                tuple(float(value) for value in row)
                if isinstance(row, np.ndarray)
                else float(row)
            )
            for row in normalized
        )
    else:
        raise QuantityExpressionError(
            "numeric arrays may have rank at most two"
        )
    return payload, canonical_unit, dimension


def convert_normalized_value(
    value: Any,
    dimension: Dimension,
    target_unit: str,
) -> Any:
    """Convert a canonical numeric value to a compatible display unit."""

    target_dimension, _, factor = _unit_spec(target_unit)
    if target_dimension != dimension:
        raise QuantityExpressionError(
            "target unit has an incompatible dimension"
        )
    array = np.asarray(value, dtype=float) / factor
    if not np.all(np.isfinite(array)):
        raise QuantityExpressionError(
            "unit conversion produced a non-finite value"
        )
    if array.ndim == 0:
        return float(array)
    if array.ndim <= 2:
        return tuple(
            (
                tuple(float(value) for value in row)
                if isinstance(row, np.ndarray)
                else float(row)
            )
            for row in array
        )
    raise QuantityExpressionError("numeric arrays may have rank at most two")


@dataclass(frozen=True)
class QuantityExpressionNodeV1:
    """One operation in a topologically ordered, bounded expression DAG."""

    node_id: str
    operation: str
    input_ids: tuple[str, ...] = ()
    reference: str = ""
    #: For ``ref``: one index per axis, outermost first. ``None`` keeps that
    #: whole axis, so ``(None, 3)`` is column 3 of a matrix -- one value per
    #: row -- and a selection by what a result measured can be planned
    #: before the result exists.
    indices: tuple[int | None, ...] = ()
    literal_value: Any = None
    literal_unit: str = "1"
    constant_name: str = ""
    scale_factor: float | None = None
    target_unit: str = ""
    cardinal_numbers: tuple[int, ...] = ()
    extrapolation_exponent: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "input_ids", tuple(self.input_ids))
        object.__setattr__(self, "indices", tuple(self.indices))
        object.__setattr__(
            self, "cardinal_numbers", tuple(self.cardinal_numbers)
        )
        if not self.node_id or len(self.node_id) > 128:
            raise QuantityContractError("expression node_id is invalid")
        if self.operation not in _OPERATIONS:
            raise QuantityContractError(
                f"unsupported expression operation: {self.operation!r}"
            )
        if len(self.input_ids) > MAX_NODE_INPUTS:
            raise QuantityContractError("expression node has too many inputs")
        require_expression_input_count(self.operation, len(self.input_ids))
        if self.operation == "constant":
            if not self.constant_name or len(self.constant_name) > 128:
                raise QuantityContractError(
                    "constant requires a registered constant_name"
                )
        elif self.constant_name:
            raise QuantityContractError(
                "constant_name applies only to the constant operation"
            )
        if any(index is not None and index < 0 for index in self.indices):
            raise QuantityContractError(
                "reference indices must be non-negative, or null to keep a "
                "whole axis"
            )
        if self.scale_factor is not None and not math.isfinite(
            self.scale_factor
        ):
            raise QuantityContractError("scale_factor must be finite")
        cbs_operations = {
            "scf_exponential_cbs_limit",
            "scf_inverse_power_cbs_limit",
            "correlation_inverse_power_cbs_limit",
        }
        if self.operation in cbs_operations:
            if (
                len(self.cardinal_numbers) != 2
                or any(
                    isinstance(value, bool)
                    or not isinstance(value, int)
                    or value < 2
                    for value in self.cardinal_numbers
                )
                or self.cardinal_numbers[0] >= self.cardinal_numbers[1]
            ):
                raise QuantityContractError(
                    "two-point CBS operations require increasing integer "
                    "cardinal_numbers >= 2"
                )
            if (
                self.extrapolation_exponent is None
                or not math.isfinite(self.extrapolation_exponent)
                or self.extrapolation_exponent <= 0.0
            ):
                raise QuantityContractError(
                    "two-point CBS operations require a positive explicit exponent"
                )
        elif self.cardinal_numbers or self.extrapolation_exponent is not None:
            raise QuantityContractError(
                "CBS cardinal numbers and exponent apply only to CBS operations"
            )


def expression_node_from_plan(
    item: Mapping[str, Any],
) -> QuantityExpressionNodeV1:
    """Coerce one planned expression dict into the typed node.

    This is the single dict-to-node reading, shared by the planning
    admission and the evaluation tool, so a node the plan admits is by
    construction a node the evaluator recognises: operation, input
    count, and per-operation field shape are all checked here, at the
    moment the plan is written, rather than after the engines have run.
    """
    return QuantityExpressionNodeV1(
        node_id=str(item.get("node_id", "")),
        operation=str(item.get("operation", "")),
        input_ids=tuple(str(value) for value in item.get("input_ids", ())),
        reference=str(item.get("reference", "") or ""),
        indices=tuple(
            None if value is None else int(value)
            for value in item.get("indices", ())
        ),
        literal_value=item.get("literal_value"),
        literal_unit=str(item.get("literal_unit", "1")),
        constant_name=str(item.get("constant_name", "")),
        scale_factor=(
            float(item["scale_factor"]) if "scale_factor" in item else None
        ),
        target_unit=str(item.get("target_unit", "")),
        cardinal_numbers=tuple(
            int(value) for value in item.get("cardinal_numbers", ())
        ),
        extrapolation_exponent=(
            float(item["extrapolation_exponent"])
            if "extrapolation_exponent" in item
            else None
        ),
    )


@dataclass(frozen=True)
class QuantityExpressionRequestV1:
    schema_version: str
    expression_id: str
    inputs: tuple[QuantityValueV1, ...]
    nodes: tuple[QuantityExpressionNodeV1, ...]
    output_node_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "inputs", tuple(self.inputs))
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(
            self, "output_node_ids", tuple(self.output_node_ids)
        )
        if self.schema_version != "chemsmart.quantity-expression-request.v1":
            raise QuantityContractError(
                "unsupported expression request schema"
            )
        if not self.expression_id or len(self.expression_id) > 128:
            raise QuantityContractError("expression_id is invalid")
        if not self.nodes or len(self.nodes) > MAX_EXPRESSION_NODES:
            raise QuantityContractError(
                f"expression requires 1..{MAX_EXPRESSION_NODES} nodes"
            )
        input_ids = [quantity.quantity_id for quantity in self.inputs]
        node_ids = [node.node_id for node in self.nodes]
        if len(input_ids) != len(set(input_ids)):
            raise QuantityContractError("input quantity IDs must be unique")
        if len(node_ids) != len(set(node_ids)):
            raise QuantityContractError("expression node IDs must be unique")
        if set(input_ids).intersection(node_ids):
            raise QuantityContractError("input and node IDs must not overlap")
        if not self.output_node_ids:
            raise QuantityContractError(
                "at least one expression output is required"
            )
        if len(self.output_node_ids) != len(set(self.output_node_ids)):
            raise QuantityContractError("expression output IDs must be unique")
        if not set(self.output_node_ids).issubset(node_ids):
            raise QuantityContractError(
                "every output must identify an expression node"
            )


#: The roles by which a number a model typed can enter an expression result.
#: Every other input to the arithmetic traces back to a measurement receipt.
MODEL_AUTHORED_CONSTANT_ROLES = (
    "extrapolation_exponent",
    "literal_value",
    "power_exponent",
    "scale_factor",
)


@dataclass(frozen=True)
class ModelAuthoredConstantV1:
    """One number that entered a result because a model wrote it down.

    The receipt closure already proves which measurements a value depends on.
    It could not, until now, distinguish a limit computed entirely from
    executed jobs from one where the model also supplied a decay exponent, a
    scale factor, or an outright literal -- and those constants move the answer
    by more than the calculations they are applied to.  Naming them makes the
    difference auditable instead of leaving it implicit in the request body.
    """

    node_id: str
    role: str
    value: str

    def __post_init__(self) -> None:
        if not self.node_id or len(self.node_id) > 128:
            raise QuantityContractError(
                "model-authored constant node ID is invalid"
            )
        if self.role not in MODEL_AUTHORED_CONSTANT_ROLES:
            raise QuantityContractError(
                f"unsupported model-authored constant role: {self.role!r}; "
                f"expected one of {list(MODEL_AUTHORED_CONSTANT_ROLES)}"
            )
        if not self.value or len(self.value) > 256:
            raise QuantityContractError(
                "model-authored constant value must be a short canonical text"
            )

    def sort_key(self) -> tuple[str, str, str]:
        return (self.node_id, self.role, self.value)


@dataclass(frozen=True)
class QuantityExpressionOutputDependencyV1:
    """Receipt closure for one expression output."""

    output_id: str
    source_receipt_sha256s: tuple[str, ...]
    model_authored_constants: tuple[ModelAuthoredConstantV1, ...] = ()
    #: Which conventions ChemSmart supplied on the way to this value.
    convention_operations: tuple[str, ...] = ()
    #: How many general-arithmetic nodes reach it.  Read together with the
    #: field above, this answers one paper-independent question: was this
    #: number produced by the toolkit's vocabulary, or assembled by the model?
    arithmetic_node_count: int = 0

    def __post_init__(self) -> None:
        if not self.output_id or len(self.output_id) > 128:
            raise QuantityContractError(
                "expression output dependency ID is invalid"
            )
        object.__setattr__(
            self,
            "model_authored_constants",
            tuple(self.model_authored_constants),
        )
        object.__setattr__(
            self, "convention_operations", tuple(self.convention_operations)
        )
        if self.convention_operations != tuple(
            sorted(set(self.convention_operations))
        ):
            raise QuantityContractError(
                "convention operations must be sorted and unique"
            )
        unknown = sorted(
            set(self.convention_operations) - CONVENTION_OPERATIONS
        )
        if unknown:
            raise QuantityContractError(
                f"unregistered convention operations: {unknown}"
            )
        if (
            isinstance(self.arithmetic_node_count, bool)
            or not isinstance(self.arithmetic_node_count, int)
            or self.arithmetic_node_count < 0
        ):
            raise QuantityContractError(
                "arithmetic_node_count must be a non-negative integer"
            )
        keys = tuple(item.sort_key() for item in self.model_authored_constants)
        if keys != tuple(sorted(set(keys))):
            raise QuantityContractError(
                "model-authored constants must be sorted and unique"
            )
        if self.source_receipt_sha256s != tuple(
            sorted(set(self.source_receipt_sha256s))
        ):
            raise QuantityContractError(
                "expression output receipt dependencies must be sorted and unique"
            )
        for digest in self.source_receipt_sha256s:
            if len(digest) != 64:
                raise QuantityContractError(
                    "expression dependency must be a SHA-256 digest"
                )
            try:
                int(digest, 16)
            except ValueError as exc:
                raise QuantityContractError(
                    "expression dependency must be a SHA-256 digest"
                ) from exc


@dataclass(frozen=True)
class QuantityExpressionReceiptV1:
    schema_version: str
    expression_id: str
    request_sha256: str
    semantic_signature_sha256: str
    node_values: tuple[QuantityValueV1, ...]
    outputs: tuple[QuantityValueV1, ...]
    output_dependencies: tuple[QuantityExpressionOutputDependencyV1, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "node_values", tuple(self.node_values))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(
            self, "output_dependencies", tuple(self.output_dependencies)
        )
        if self.schema_version != "chemsmart.quantity-expression-receipt.v1":
            raise QuantityContractError(
                "unsupported expression receipt schema"
            )
        if self.status != "derived":
            raise QuantityContractError("invalid expression receipt status")
        dependency_ids = tuple(
            dependency.output_id for dependency in self.output_dependencies
        )
        output_ids = tuple(output.quantity_id for output in self.outputs)
        if dependency_ids != output_ids:
            raise QuantityContractError(
                "expression output dependency order must match outputs"
            )
        for digest, label in (
            (self.request_sha256, "request_sha256"),
            (self.semantic_signature_sha256, "semantic_signature_sha256"),
        ):
            if len(digest) != 64:
                raise QuantityContractError(
                    f"{label} must be a SHA-256 digest"
                )
            try:
                int(digest, 16)
            except ValueError as exc:
                raise QuantityContractError(
                    f"{label} must be a SHA-256 digest"
                ) from exc
        body = {
            "schema_version": self.schema_version,
            "expression_id": self.expression_id,
            "request_sha256": self.request_sha256,
            "semantic_signature_sha256": self.semantic_signature_sha256,
            "node_values": self.node_values,
            "outputs": self.outputs,
            "output_dependencies": self.output_dependencies,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_quantity_sha256(body):
            raise QuantityContractError(
                "quantity expression receipt digest mismatch"
            )


def quantity_expression_semantic_signature(
    request: QuantityExpressionRequestV1,
) -> str:
    """Hash an identifier-independent symbolic form of an expression DAG.

    Model-chosen local aliases and intermediate node names are deliberately
    excluded.  Leaves are identified by an explicit semantic role when one was
    supplied, otherwise by the source quantity ID carried in deterministic
    evidence.  This lets equivalent formulae grade identically across papers
    and artifacts while still distinguishing, for example, ``reactant`` and
    ``product`` energies when those roles are declared.

    Numerical results, artifact IDs, and receipt hashes are never part of this
    signature.  Literal constants, units, indices, operation ordering, and the
    requested output roles remain authoritative.
    """

    values: dict[str, dict[str, Any]] = {}
    source_quantity_ids = {
        quantity.quantity_id: (
            matches[-1]
            if (matches := _QUANTITY_REF.findall(quantity.evidence_ref))
            else ""
        )
        for quantity in request.inputs
    }
    source_quantity_counts = Counter(source_quantity_ids.values())
    source_quantity_counts.pop("", None)
    role_owners: dict[str, str] = {}
    for quantity in request.inputs:
        role_matches = _SEMANTIC_ROLE_REF.findall(quantity.evidence_ref)
        source_quantity_id = source_quantity_ids[quantity.quantity_id]
        # A source quantity drawn once names its own role, which keeps the
        # signature readable.  Once the same source name arrives more than
        # once it can no longer identify anything, so a declared role decides,
        # and failing that the input's own id does: it is already required to
        # be unique within the expression, so it is a role by construction and
        # can never be ambiguous.
        #
        # Falling back to the shared source name instead, as this did, made
        # the collision automatic: comparing two structures repeats a quantity
        # name by definition, so any session that did not hand-author a role
        # for every occurrence was refused.  Three cycles of explaining that
        # requirement each halved the failure without clearing it.  The host
        # can derive what it was asking the model to supply.
        semantic_role = (
            source_quantity_id
            if source_quantity_id
            and source_quantity_counts[source_quantity_id] == 1
            else (role_matches[-1] if role_matches else quantity.quantity_id)
        )
        owner = role_owners.setdefault(semantic_role, quantity.quantity_id)
        if owner != quantity.quantity_id:
            raise QuantityContractError(
                "repeated source quantities require distinct semantic_role values"
            )
        values[quantity.quantity_id] = {
            "kind": "input",
            "semantic_role": semantic_role,
            "data_kind": quantity.data_kind,
            "unit": quantity.unit,
            "dimension": quantity.dimension,
        }

    node_by_id = {node.node_id: node for node in request.nodes}
    visiting: set[str] = set()

    def canonical_value(value_id: str) -> dict[str, Any]:
        if value_id in values:
            return values[value_id]
        try:
            node = node_by_id[value_id]
        except KeyError as exc:
            raise QuantityContractError(
                "semantic signature references an unknown expression value"
            ) from exc
        if value_id in visiting:
            raise QuantityContractError("expression graph contains a cycle")
        visiting.add(value_id)
        if node.operation == "literal":
            normalized, canonical_unit, dimension = normalize_numeric_value(
                node.literal_value, node.literal_unit
            )
            result: dict[str, Any] = {
                "operation": "literal",
                "value": normalized,
                "unit": canonical_unit,
                "dimension": dimension,
            }
        elif node.operation == "constant":
            try:
                entry = literature_constant(node.constant_name)
            except UnknownLiteratureConstantError as exc:
                raise QuantityContractError(str(exc)) from exc
            normalized, canonical_unit, dimension = normalize_numeric_value(
                entry.value, entry.unit
            )
            result = {
                "operation": "constant",
                "name": entry.name,
                "value": normalized,
                "unit": canonical_unit,
                "dimension": dimension,
            }
        elif node.operation == "ref":
            reference = node.reference
            if not reference and len(node.input_ids) == 1:
                reference = node.input_ids[0]
            if not reference:
                raise QuantityContractError(
                    "semantic ref requires exactly one prior value"
                )
            source = canonical_value(reference)
            result = (
                {
                    "operation": "ref",
                    "source": source,
                    "indices": node.indices,
                }
                if node.indices
                else source
            )
        else:
            inputs = [canonical_value(input_id) for input_id in node.input_ids]
            if node.operation == "convert" and len(inputs) == 1:
                # Arithmetic values are already normalized to canonical units;
                # display-unit choice is enforced by the claim contract.
                result = inputs[0]
                visiting.remove(value_id)
                values[value_id] = result
                return result
            if node.operation in {
                "add",
                "multiply",
                "sum",
                "mean",
                "min",
                "max",
            }:
                inputs = sorted(inputs, key=canonical_quantity_sha256)
            elif node.operation == "distance" and len(inputs) == 2:
                inputs = sorted(inputs, key=canonical_quantity_sha256)
            elif node.operation == "angle" and len(inputs) == 3:
                outer = sorted(
                    (inputs[0], inputs[2]), key=canonical_quantity_sha256
                )
                inputs = [outer[0], inputs[1], outer[1]]
            result = {
                "operation": node.operation,
                "inputs": tuple(inputs),
            }
            if node.scale_factor is not None:
                result["scale_factor"] = node.scale_factor
            if node.target_unit:
                target_dimension, target_unit, _ = _unit_spec(node.target_unit)
                result["target_unit"] = target_unit
                result["target_dimension"] = target_dimension
            if node.indices:
                result["indices"] = node.indices
            if node.cardinal_numbers:
                result["cardinal_numbers"] = node.cardinal_numbers
            if node.extrapolation_exponent is not None:
                result["extrapolation_exponent"] = node.extrapolation_exponent
        visiting.remove(value_id)
        values[value_id] = result
        return result

    body = {
        "schema_version": "chemsmart.quantity-expression-semantic-signature.v2",
        "expression_id": request.expression_id,
        "outputs": tuple(
            {
                "output_id": output_id,
                "expression": canonical_value(output_id),
            }
            for output_id in sorted(request.output_node_ids)
        ),
    }
    return canonical_quantity_sha256(body)


def _receipt_dependencies(evidence_ref: str) -> frozenset[str]:
    return frozenset(
        match.group(1) for match in _RECEIPT_REF.finditer(evidence_ref)
    )


def _numeric(quantity: QuantityValueV1) -> np.ndarray:
    if quantity.data_kind in {"text", "text_vector"}:
        raise QuantityExpressionError(
            f"quantity {quantity.quantity_id!r} is metadata, not numeric data"
        )
    array = np.asarray(quantity.value, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise QuantityExpressionError(
            "expression input is non-finite or empty"
        )
    return array


def _connectivity_observation(
    matrix_quantity: QuantityValueV1,
    symbols_quantity: QuantityValueV1,
    label: str,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Validate one geometry-perceived adjacency and its atom-order witness."""

    if matrix_quantity.dimension != DIMENSIONLESS:
        raise QuantityExpressionError(
            f"{label} connectivity must be dimensionless"
        )
    try:
        matrix = np.asarray(matrix_quantity.value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise QuantityExpressionError(
            f"{label} connectivity must be a numeric matrix"
        ) from exc
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise QuantityExpressionError(
            f"{label} connectivity must be a square matrix; got {matrix.shape}"
        )
    if matrix.shape[0] == 0 or not np.all(np.isfinite(matrix)):
        raise QuantityExpressionError(
            f"{label} connectivity is empty or non-finite"
        )
    if not np.all(np.isin(matrix, (0.0, 1.0))):
        raise QuantityExpressionError(
            f"{label} connectivity must contain only binary 0/1 values"
        )
    if not np.array_equal(matrix, matrix.T):
        raise QuantityExpressionError(
            f"{label} connectivity must be symmetric for undirected edges"
        )
    if np.any(np.diag(matrix) != 0.0):
        raise QuantityExpressionError(
            f"{label} connectivity must have a zero diagonal"
        )
    if symbols_quantity.data_kind != "text_vector":
        raise QuantityExpressionError(
            f"{label} symbols must be a text-vector atom-order witness"
        )
    symbols = tuple(str(item) for item in symbols_quantity.value)
    if len(symbols) != matrix.shape[0] or any(not item for item in symbols):
        raise QuantityExpressionError(
            f"{label} symbol count must equal connectivity size"
        )
    return matrix, symbols


def _payload(array: np.ndarray | float) -> Any:
    value = np.asarray(array, dtype=float)
    if value.size == 0 or not np.all(np.isfinite(value)):
        raise QuantityExpressionError("expression produced a non-finite value")
    if value.ndim == 0:
        return float(value)
    if value.ndim == 1:
        return tuple(float(item) for item in value)
    if value.ndim == 2:
        return tuple(tuple(float(item) for item in row) for row in value)
    raise QuantityExpressionError("expression output rank exceeds two")


def _node_value(
    *,
    node: QuantityExpressionNodeV1,
    values: Mapping[str, QuantityValueV1],
    evidence_ref: str,
) -> QuantityValueV1:
    operation = node.operation
    if operation == "ref":
        reference = node.reference
        if not reference and len(node.input_ids) == 1:
            reference = node.input_ids[0]
        elif node.input_ids:
            raise QuantityExpressionError(
                "ref accepts reference or one input_id, never both"
            )
        if not reference:
            raise QuantityExpressionError(
                "ref requires one prior value reference"
            )
        try:
            source = values[reference]
        except KeyError as exc:
            raise QuantityExpressionError(
                f"unknown quantity reference: {reference!r}"
            ) from exc
        if not node.indices:
            return make_quantity_value(
                quantity_id=node.node_id,
                source_value=source.source_value,
                source_unit=source.source_unit,
                value=source.value,
                unit=source.unit,
                dimension=source.dimension,
                evidence_ref=source.evidence_ref,
                data_kind=source.data_kind,
            )
        # One index per axis, taken together, so a null can keep its whole
        # axis: a mode's row was selectable and an atom's column was not, and
        # the mode in which given atoms move most could be chosen only by an
        # index read after the Hessian existed (R9 xtb g3).
        key = tuple(
            slice(None) if index is None else index for index in node.indices
        )
        try:
            selected = np.asarray(_numeric(source)[key])
        except IndexError as exc:
            raise QuantityExpressionError(
                "reference index is out of range"
            ) from exc
        payload = _payload(selected)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=source.unit,
            value=payload,
            unit=source.unit,
            dimension=source.dimension,
            evidence_ref=source.evidence_ref,
        )

    if operation == "literal":
        if node.input_ids or node.reference or node.literal_value is None:
            raise QuantityExpressionError(
                "literal requires literal_value and no references or inputs"
            )
        normalized, canonical_unit, dimension = normalize_numeric_value(
            node.literal_value, node.literal_unit
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=node.literal_value,
            source_unit=node.literal_unit,
            value=normalized,
            unit=canonical_unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "constant":
        if node.input_ids or node.reference or node.literal_value is not None:
            raise QuantityExpressionError(
                "constant takes a registered constant_name and no "
                "references, inputs, or literal value"
            )
        try:
            entry = literature_constant(node.constant_name)
        except UnknownLiteratureConstantError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        normalized, canonical_unit, dimension = normalize_numeric_value(
            entry.value, entry.unit
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=entry.value,
            source_unit=entry.unit,
            value=normalized,
            unit=canonical_unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    try:
        inputs = [values[input_id] for input_id in node.input_ids]
    except KeyError as exc:
        raise QuantityExpressionError(
            f"node references an unavailable prior value: {exc.args[0]!r}"
        ) from exc

    if operation == "connectivity_difference_count":
        if len(inputs) != 4:
            raise QuantityExpressionError(
                "connectivity_difference_count requires first connectivity, "
                "first symbols, second connectivity, and second symbols"
            )
        first, first_symbols = _connectivity_observation(
            inputs[0], inputs[1], "first"
        )
        second, second_symbols = _connectivity_observation(
            inputs[2], inputs[3], "second"
        )
        if first.shape != second.shape:
            raise QuantityExpressionError(
                "connectivity matrices have different atom counts"
            )
        if first_symbols != second_symbols:
            raise QuantityExpressionError(
                "connectivity inputs differ in atom identity or atom order"
            )
        changed_edges = int(np.count_nonzero(np.triu(first != second, k=1)))
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=changed_edges,
            source_unit="1",
            value=changed_edges,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )

    if operation in {"add", "subtract", "multiply", "divide"}:
        if len(inputs) != 2:
            raise QuantityExpressionError(f"{operation} requires two inputs")
        left, right = inputs
        left_value, right_value = _numeric(left), _numeric(right)
        if operation in {"add", "subtract"}:
            if left.dimension != right.dimension:
                raise QuantityExpressionError(
                    f"{operation} requires identical dimensions"
                )
            # A scalar broadcasts across a vector, as it already does for
            # multiply and divide: a relaxed scan's profile against its own
            # minimum is subtract(energies, min(energies)), and refusing it
            # here after every engine had run lost a delivery's whole
            # analysis chain (live, 2026-09-02). Two vectors still need
            # identical shapes.
            if (
                left_value.ndim
                and right_value.ndim
                and left_value.shape != right_value.shape
            ):
                raise QuantityExpressionError(
                    f"{operation} accepts scalar broadcasting or identical "
                    "shapes"
                )
            result = (
                left_value + right_value
                if operation == "add"
                else left_value - right_value
            )
            dimension = left.dimension
        elif operation == "multiply":
            if (
                left_value.ndim
                and right_value.ndim
                and left_value.shape != right_value.shape
            ):
                raise QuantityExpressionError(
                    "multiply accepts scalar broadcasting or identical shapes"
                )
            result = left_value * right_value
            dimension = _add_dimensions(left.dimension, right.dimension)
        else:
            if np.any(right_value == 0.0):
                raise QuantityExpressionError("division by zero is forbidden")
            if (
                left_value.ndim
                and right_value.ndim
                and left_value.shape != right_value.shape
            ):
                raise QuantityExpressionError(
                    "divide accepts scalar broadcasting or identical shapes"
                )
            result = left_value / right_value
            dimension = _subtract_dimensions(left.dimension, right.dimension)
        unit = canonical_unit_for_dimension(dimension)
        payload = _payload(result)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "scale":
        if len(inputs) != 1 or node.scale_factor is None:
            raise QuantityExpressionError(
                "scale requires one input and scale_factor"
            )
        source = inputs[0]
        payload = _payload(_numeric(source) * node.scale_factor)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=source.unit,
            value=payload,
            unit=source.unit,
            dimension=source.dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "abs":
        if len(inputs) != 1:
            raise QuantityExpressionError("abs requires one input")
        source = inputs[0]
        payload = _payload(np.abs(_numeric(source)))
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=source.unit,
            value=payload,
            unit=source.unit,
            dimension=source.dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "sqrt":
        if len(inputs) != 1:
            raise QuantityExpressionError("sqrt requires one input")
        source = inputs[0]
        if any(exponent % 2 for exponent in source.dimension):
            raise QuantityExpressionError(
                "sqrt requires even exponents in every physical dimension"
            )
        source_value = _numeric(source)
        if np.any(source_value < 0.0):
            raise QuantityExpressionError("sqrt requires non-negative values")
        dimension = tuple(exponent // 2 for exponent in source.dimension)
        payload = _payload(np.sqrt(source_value))
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "power":
        if len(inputs) != 1 or node.literal_value is None:
            raise QuantityExpressionError(
                "power requires one input and a literal exponent"
            )
        exponent_array = np.asarray(node.literal_value, dtype=float)
        if exponent_array.ndim != 0 or not np.isfinite(exponent_array):
            raise QuantityExpressionError(
                "power exponent must be a finite scalar"
            )
        exponent = float(exponent_array)
        source = inputs[0]
        if source.dimension != DIMENSIONLESS:
            rounded = round(exponent)
            if not math.isclose(
                exponent, rounded, rel_tol=0.0, abs_tol=1.0e-12
            ):
                raise QuantityExpressionError(
                    "dimensioned power requires an integer exponent"
                )
            if abs(rounded) > 12:
                raise QuantityExpressionError(
                    "dimensioned power exponent is outside the bounded range"
                )
            dimension = tuple(
                int(value * rounded) for value in source.dimension
            )
        else:
            dimension = DIMENSIONLESS
        source_value = _numeric(source)
        if np.any(source_value < 0.0) and not math.isclose(
            exponent, round(exponent), rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise QuantityExpressionError(
                "fractional power of a negative value is unsupported"
            )
        if exponent < 0.0 and np.any(source_value == 0.0):
            raise QuantityExpressionError(
                "negative power of zero is forbidden"
            )
        payload = _payload(np.power(source_value, exponent))
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {"exp", "log"}:
        if len(inputs) != 1 or inputs[0].dimension != DIMENSIONLESS:
            raise QuantityExpressionError(
                f"{operation} requires one dimensionless input"
            )
        source_value = _numeric(inputs[0])
        if operation == "log" and np.any(source_value <= 0.0):
            raise QuantityExpressionError(
                "log requires strictly positive values"
            )
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            result = (
                np.exp(source_value)
                if operation == "exp"
                else np.log(source_value)
            )
        payload = _payload(result)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="1",
            value=payload,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )

    if operation == "photon_wavelength":
        if len(inputs) != 1 or inputs[0].dimension != ENERGY:
            raise QuantityExpressionError(
                "photon_wavelength requires one positive energy input"
            )
        energy_hartree = _numeric(inputs[0])
        if np.any(energy_hartree <= 0.0):
            raise QuantityExpressionError(
                "photon_wavelength requires strictly positive energies"
            )
        energy_ev = energy_hartree * energy_conversion("hartree", "eV", 1.0)
        hc_ev_angstrom = (
            ase_units._hplanck * ase_units._c / ase_units._e * 1.0e10
        )
        payload = _payload(hc_ev_angstrom / energy_ev)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="angstrom",
            value=payload,
            unit="angstrom",
            dimension=LENGTH,
            evidence_ref=evidence_ref,
        )

    if operation == "gibbs_to_pka":
        from chemsmart.analysis.aggregation import (
            GAS_CONSTANT_KCAL,
            HARTREE_TO_KCAL_PER_MOL,
        )

        if (
            len(inputs) != 2
            or inputs[0].dimension != ENERGY
            or inputs[1].dimension != TEMPERATURE
        ):
            raise QuantityExpressionError(
                "gibbs_to_pka takes one deprotonation free energy and one "
                "temperature, in that order"
            )
        delta_g = _numeric(inputs[0])
        temperature = _numeric(inputs[1])
        if delta_g.size != 1 or temperature.size != 1:
            raise QuantityExpressionError("gibbs_to_pka takes scalar inputs")
        temperature_k = float(temperature.reshape(()))
        if temperature_k <= 0.0:
            raise QuantityExpressionError(
                "gibbs_to_pka requires a positive temperature"
            )
        delta_g_kcal = float(delta_g.reshape(())) * HARTREE_TO_KCAL_PER_MOL
        payload = _payload(
            delta_g_kcal / (GAS_CONSTANT_KCAL * temperature_k * math.log(10.0))
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="1",
            value=payload,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )

    if operation == "gibbs_to_redox_potential":
        # E = -deltaG / (n F).  Free energies here are per mole and F is
        # Avogadro's number times the elementary charge, so per electron the
        # arithmetic is just -deltaG / n; the elementary charge enters as the
        # unit of the result rather than as a factor.  That is why no Faraday
        # constant appears in this function or in the constants registry.
        if (
            len(inputs) != 2
            or inputs[0].dimension != ENERGY
            or inputs[1].dimension != DIMENSIONLESS
        ):
            raise QuantityExpressionError(
                "gibbs_to_redox_potential takes one half-reaction free "
                "energy and one electron count, in that order"
            )
        delta_g = _numeric(inputs[0])
        electrons = _numeric(inputs[1])
        if delta_g.size != 1 or electrons.size != 1:
            raise QuantityExpressionError(
                "gibbs_to_redox_potential takes scalar inputs"
            )
        electron_count = float(electrons.reshape(()))
        if electron_count <= 0.0:
            raise QuantityExpressionError(
                "gibbs_to_redox_potential requires a positive electron count"
            )
        payload = _payload(-float(delta_g.reshape(())) / electron_count)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="hartree e^-1",
            value=payload,
            unit="hartree e^-1",
            dimension=ELECTRIC_POTENTIAL,
            evidence_ref=evidence_ref,
        )

    if operation in {"center_of_mass", "principal_moments_of_inertia"}:
        if (
            len(inputs) != 2
            or inputs[0].dimension != LENGTH
            or inputs[1].dimension != MASS
        ):
            raise QuantityExpressionError(
                f"{operation} requires a coordinate matrix followed by an "
                "atomic-mass vector"
            )
        positions = _numeric(inputs[0])
        masses = _numeric(inputs[1]).reshape(-1)
        if (
            positions.ndim != 2
            or positions.shape[1] != 3
            or positions.shape[0] != masses.size
        ):
            raise QuantityExpressionError(
                f"{operation} requires positions shaped (N, 3) and N masses"
            )
        if np.any(masses <= 0.0):
            raise QuantityExpressionError("atomic masses must be positive")
        if operation == "center_of_mass":
            payload = _payload(np.average(positions, axis=0, weights=masses))
            unit = "angstrom"
            dimension = LENGTH
        else:
            from chemsmart.utils.geometry import calculate_moments_of_inertia

            _, moments, _ = calculate_moments_of_inertia(masses, positions)
            moments = np.asarray(moments, dtype=float)
            tolerance = max(1.0, float(np.max(np.abs(moments)))) * 1.0e-12
            moments[np.abs(moments) < tolerance] = 0.0
            if np.any(moments < 0.0):
                raise QuantityExpressionError(
                    "principal moments contain a negative eigenvalue"
                )
            payload = _payload(np.sort(moments))
            unit = "angstrom^2 u"
            dimension = MOMENT_OF_INERTIA
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {"linear_rotor_constant", "rigid_rotor_constants"}:
        if len(inputs) != 1 or inputs[0].dimension != MOMENT_OF_INERTIA:
            raise QuantityExpressionError(
                f"{operation} requires one principal-moment vector"
            )
        moments = np.sort(_numeric(inputs[0]).reshape(-1))
        if moments.size != 3:
            raise QuantityExpressionError(
                f"{operation} requires exactly three principal moments"
            )
        if operation == "linear_rotor_constant":
            perpendicular_moment = float(np.mean(moments[1:]))
            tolerance = max(1.0, perpendicular_moment) * 1.0e-6
            if (
                moments[0] < -tolerance
                or abs(float(moments[0])) > tolerance
                or np.any(moments[1:] <= 0.0)
                or abs(float(moments[2] - moments[1])) > tolerance
            ):
                raise QuantityExpressionError(
                    "linear_rotor_constant requires one near-zero axial "
                    "moment and two equal positive perpendicular moments"
                )
            moments = np.asarray((perpendicular_moment,), dtype=float)
        elif np.any(moments <= 0.0):
            raise QuantityExpressionError(
                "rigid_rotor_constants requires three positive moments for "
                "a nonlinear molecule"
            )
        moments_si = (
            moments * ase_units._amu * (ase_units.Ang / ase_units.m) ** 2
        )
        frequencies_hz = ase_units._hplanck / (8.0 * np.pi**2 * moments_si)
        rotational_constants = frequencies_hz / (ase_units._c * 100.0)
        payload = (
            float(rotational_constants[0])
            if operation == "linear_rotor_constant"
            else _payload(rotational_constants)
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="cm^-1",
            value=payload,
            unit="cm^-1",
            dimension=FREQUENCY,
            evidence_ref=evidence_ref,
        )

    if operation == "exponential_cbs_limit":
        # A Hartree-Fock basis-set series converges exponentially, so the
        # complete-basis limit is a three-parameter fit rather than a linear
        # one. Papers state it as "extrapolated with a three-parameter
        # exponential formula"; without this operation such a limit has no
        # expressible producer and the workflow stops one step short of the
        # number it was built to obtain.
        from chemsmart.analysis.aggregation import (
            AggregationError,
            extrapolate_exponential_three_point,
        )

        # Accept either shape the rest of the harness can actually produce.
        # Extraction yields one scalar energy per calculation, so demanding a
        # single three-element vector left a model with three separate results
        # no way to reach this operation -- and a live session responded by
        # rebuilding the closed form out of multiply, subtract and scale nodes,
        # reintroducing by hand the convention this operation exists to own.
        if len(inputs) == 3:
            if any(_numeric(item).ndim != 0 for item in inputs):
                raise QuantityExpressionError(
                    "exponential_cbs_limit takes three scalar energies, "
                    "ordered by increasing basis cardinal, or one input "
                    "holding all three"
                )
            if len({item.dimension for item in inputs}) != 1:
                raise QuantityExpressionError(
                    "exponential_cbs_limit needs three energies of one "
                    "dimension"
                )
            series = np.asarray(
                [float(_numeric(item)) for item in inputs], dtype=float
            )
        elif len(inputs) == 1:
            series = _numeric(inputs[0]).reshape(-1)
        else:
            raise QuantityExpressionError(
                "exponential_cbs_limit takes the energies at three equally "
                "spaced cardinal numbers, ordered by increasing basis, either "
                f"as three scalar inputs or as one three-element input; got "
                f"{len(inputs)} inputs"
            )
        if series.size != 3:
            raise QuantityExpressionError(
                "exponential_cbs_limit needs exactly three energies, got "
                f"{series.size}"
            )
        try:
            payload = extrapolate_exponential_three_point(
                (float(series[0]), float(series[1]), float(series[2]))
            )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        dimension = inputs[0].dimension
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {
        "scf_exponential_cbs_limit",
        "scf_inverse_power_cbs_limit",
        "correlation_inverse_power_cbs_limit",
    }:
        if (
            len(inputs) != 2
            or any(item.dimension != ENERGY for item in inputs)
            or any(_numeric(item).ndim != 0 for item in inputs)
        ):
            raise QuantityExpressionError(
                f"{operation} requires two scalar energy inputs ordered by "
                "increasing basis cardinal"
            )
        smaller_cardinal, larger_cardinal = node.cardinal_numbers
        smaller_energy = float(_numeric(inputs[0]))
        larger_energy = float(_numeric(inputs[1]))
        from chemsmart.analysis.aggregation import (
            AggregationError,
            extrapolate_correlation_inverse_power,
            extrapolate_scf_exponential,
            extrapolate_scf_inverse_power,
        )

        try:
            if operation == "scf_exponential_cbs_limit":
                payload = extrapolate_scf_exponential(
                    smaller_cardinal=smaller_cardinal,
                    larger_cardinal=larger_cardinal,
                    smaller_scf_energy=smaller_energy,
                    larger_scf_energy=larger_energy,
                    alpha=float(node.extrapolation_exponent),
                )
            elif operation == "scf_inverse_power_cbs_limit":
                payload = extrapolate_scf_inverse_power(
                    smaller_cardinal=smaller_cardinal,
                    larger_cardinal=larger_cardinal,
                    smaller_scf_energy=smaller_energy,
                    larger_scf_energy=larger_energy,
                    exponent=float(node.extrapolation_exponent),
                )
            else:
                payload = extrapolate_correlation_inverse_power(
                    smaller_cardinal=smaller_cardinal,
                    larger_cardinal=larger_cardinal,
                    smaller_correlation_energy=smaller_energy,
                    larger_correlation_energy=larger_energy,
                    exponent=float(node.extrapolation_exponent),
                )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="hartree",
            value=payload,
            unit="hartree",
            dimension=ENERGY,
            evidence_ref=evidence_ref,
        )

    if operation in {"sum", "mean", "min", "max"}:
        if not inputs:
            raise QuantityExpressionError(f"{operation} requires input values")
        dimension = inputs[0].dimension
        if any(item.dimension != dimension for item in inputs):
            raise QuantityExpressionError(
                f"{operation} requires quantities of one dimension"
            )
        arrays = [_numeric(item).reshape(-1) for item in inputs]
        joined = np.concatenate(arrays)
        reducer = {
            "sum": np.sum,
            "mean": np.mean,
            "min": np.min,
            "max": np.max,
        }[operation]
        payload = float(reducer(joined))
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {"coordinate_at_maximum", "coordinate_at_minimum"}:
        # max() gives a barrier's height; nothing gave its position, because
        # the reducers above discard the index they reduced over.  Locating an
        # extremum is therefore not reachable by composing them, which is why
        # this is an operation rather than a recipe.
        if len(inputs) != 2:
            raise QuantityExpressionError(
                f"{operation} requires exactly two inputs in order, the "
                "values being extremised and the coordinate they were "
                f"measured at; got {len(inputs)}"
            )
        values = _numeric(inputs[0]).reshape(-1)
        coordinates = _numeric(inputs[1]).reshape(-1)
        if values.size != coordinates.size:
            raise QuantityExpressionError(
                f"{operation} needs one coordinate per value; got "
                f"{values.size} values and {coordinates.size} coordinates"
            )
        if values.size < 2:
            raise QuantityExpressionError(
                f"{operation} needs a series; a single point has no extremum"
            )
        index = int(
            np.argmax(values)
            if operation == "coordinate_at_maximum"
            else np.argmin(values)
        )
        # The answer is a coordinate, so it carries the coordinate's dimension
        # and not the extremised quantity's.
        dimension = inputs[1].dimension
        payload = float(coordinates[index])
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "harmonic_zero_point_energy":
        from chemsmart.analysis.aggregation import (
            AggregationError,
            harmonic_zero_point_energy,
        )

        if len(inputs) != 1 or inputs[0].dimension != FREQUENCY:
            raise QuantityExpressionError(
                "harmonic_zero_point_energy takes exactly one frequency "
                f"vector in cm^-1; got {len(inputs)} inputs"
            )
        unit = str(node.target_unit or "kJ/mol")
        try:
            payload = _payload(
                float(
                    harmonic_zero_point_energy(
                        tuple(
                            float(item)
                            for item in _numeric(inputs[0]).reshape(-1)
                        ),
                        unit=unit,
                    )
                )
            )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        normalized, canonical_unit, dimension = normalize_numeric_value(
            payload, unit
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=normalized,
            unit=canonical_unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {"wavenumber_to_energy", "energy_to_wavenumber"}:
        # h*c*N_A, stated once: one hartree is 219474.6313705 cm^-1
        # (CODATA 2018). A session asked four times to convert an
        # exchange coupling from hartree to cm^-1 and was refused each
        # time because the two dimensions never met in the vocabulary
        # (NOVEL-3 ino2, 2026-09-05).
        hartree_in_cm1 = 219474.6313705
        if len(inputs) != 1:
            raise QuantityExpressionError(
                f"{operation} takes exactly one input; got {len(inputs)}"
            )
        source = inputs[0]
        array = _numeric(source)
        if operation == "wavenumber_to_energy":
            if source.dimension != FREQUENCY:
                raise QuantityExpressionError(
                    "wavenumber_to_energy takes a wavenumber (cm^-1); got "
                    f"dimension {canonical_unit_for_dimension(source.dimension)!r}"
                )
            unit = str(node.target_unit or "kJ/mol")
            target_dimension, _target_unit, target_scale = _unit_spec(unit)
            if target_dimension != ENERGY:
                raise QuantityExpressionError(
                    f"wavenumber_to_energy target_unit {unit!r} is not an "
                    "energy unit"
                )
            converted = (array / hartree_in_cm1) / target_scale
        else:
            if source.dimension != ENERGY:
                raise QuantityExpressionError(
                    "energy_to_wavenumber takes a molar energy; got "
                    f"dimension {canonical_unit_for_dimension(source.dimension)!r}"
                )
            unit = "cm^-1"
            converted = array * hartree_in_cm1
        payload = _payload(
            float(converted.reshape(-1)[0])
            if converted.size == 1 and source.data_kind == "scalar"
            else converted
        )
        normalized, canonical_unit, dimension = normalize_numeric_value(
            payload, unit
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=normalized,
            unit=canonical_unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
            data_kind=source.data_kind,
        )

    if operation == "imaginary_mode_count":
        from chemsmart.analysis.aggregation import (
            AggregationError,
            count_imaginary_modes,
        )

        if len(inputs) not in (1, 2) or inputs[0].dimension != FREQUENCY:
            raise QuantityExpressionError(
                "imaginary_mode_count takes one frequency vector, optionally "
                f"followed by a frequency cutoff; got {len(inputs)} inputs"
            )
        cutoff = 0.0
        if len(inputs) == 2:
            if inputs[1].dimension != FREQUENCY:
                raise QuantityExpressionError(
                    "the optional cutoff for imaginary_mode_count is a "
                    "frequency"
                )
            cutoff = abs(float(_numeric(inputs[1])))
        try:
            payload = _payload(
                float(
                    count_imaginary_modes(
                        tuple(
                            float(item)
                            for item in _numeric(inputs[0]).reshape(-1)
                        ),
                        cutoff_cm1=cutoff,
                    )
                )
            )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="1",
            value=payload,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=evidence_ref,
        )

    if operation == "transition_state_crossover_temperature":
        from chemsmart.analysis.aggregation import (
            AggregationError,
            transition_state_crossover_temperature,
        )

        if len(inputs) != 1 or inputs[0].dimension != FREQUENCY:
            raise QuantityExpressionError(
                "transition_state_crossover_temperature takes one frequency "
                "vector or one selected frequency in cm^-1"
            )
        try:
            payload = _payload(
                float(
                    transition_state_crossover_temperature(
                        tuple(
                            float(item)
                            for item in _numeric(inputs[0]).reshape(-1)
                        )
                    )
                )
            )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="K",
            value=payload,
            unit="K",
            dimension=TEMPERATURE,
            evidence_ref=evidence_ref,
        )

    if operation in {"boltzmann_populations", "boltzmann_average"}:
        # ChemSmart owns this weighting.  Left unexposed, a paper reporting
        # conformer populations forces a model to rebuild exp(-dG/RT)/sum from
        # exp, divide and sum, with R, T and the unit conversion entering as
        # constants it chose -- the same failure the basis-set limit showed.
        from chemsmart.analysis.aggregation import (
            AggregationError,
            boltzmann_average,
            boltzmann_populations,
        )

        wants_values = operation == "boltzmann_average"
        expected = 3 if wants_values else 2
        # An optional trailing dimensionless input carries the per-state
        # multiplicities.  Without it a model that knows two gauche forms are
        # enantiomers has to fold the factor into the energies or scale the
        # result by hand, and the factor stops being visible as a scientific
        # input.
        degeneracies = None
        if len(inputs) >= expected + 1 and (
            inputs[-1].dimension == DIMENSIONLESS
        ):
            candidate = inputs[-1]
            degeneracies = tuple(
                float(item) for item in _numeric(candidate).reshape(-1)
            )
            inputs = inputs[:-1]
        elif (
            len(inputs) >= expected + 1 and inputs[-2].dimension == TEMPERATURE
        ):
            raise QuantityExpressionError(
                "per-state degeneracies must be dimensionless; the trailing "
                f"input to {operation} carries unit '{inputs[-1].unit}'"
            )
        if wants_values and len(inputs) != expected:
            raise QuantityExpressionError(
                f"{operation} requires "
                + (
                    "the values, their state energies, and a temperature"
                    if wants_values
                    else "the state energies and a temperature"
                )
                + ", optionally followed by the per-state degeneracies; got "
                + f"{len(inputs)} inputs"
            )
        if not wants_values and len(inputs) < expected:
            raise QuantityExpressionError(
                "boltzmann_populations requires one energy vector or two or "
                "more scalar state energies followed by a temperature, "
                "optionally followed by per-state degeneracies"
            )
        temperature_value = inputs[-1]
        if temperature_value.dimension != TEMPERATURE:
            raise QuantityExpressionError(
                f"the last input to {operation} must be a temperature"
            )
        temperature = float(_numeric(temperature_value))
        if wants_values:
            energy_inputs = inputs[1:2]
        else:
            energy_inputs = inputs[:-1]
        if any(item.dimension != ENERGY for item in energy_inputs):
            raise QuantityExpressionError(
                f"{operation} weights states by energy"
            )
        if len(energy_inputs) == 1:
            energies = _numeric(energy_inputs[0]).reshape(-1)
        else:
            scalar_energies = []
            for item in energy_inputs:
                values = _numeric(item).reshape(-1)
                if values.size != 1:
                    raise QuantityExpressionError(
                        "separate boltzmann_populations energy inputs must be "
                        "scalars; pass multiple states in one vector instead"
                    )
                scalar_energies.append(float(values[0]))
            energies = np.asarray(scalar_energies, dtype=float)
        try:
            populations = boltzmann_populations(
                tuple(float(item) for item in energies),
                temperature=temperature,
                unit="hartree",
                degeneracies=degeneracies,
            )
            if not wants_values:
                payload = _payload(np.asarray(populations, dtype=float))
                return make_quantity_value(
                    quantity_id=node.node_id,
                    source_value=payload,
                    source_unit="1",
                    value=payload,
                    unit="1",
                    dimension=DIMENSIONLESS,
                    evidence_ref=evidence_ref,
                )
            values = _numeric(inputs[0]).reshape(-1)
            if values.size != energies.size:
                raise QuantityExpressionError(
                    f"{operation} needs one value per state; got "
                    f"{values.size} values and {energies.size} energies"
                )
            payload = _payload(
                boltzmann_average(
                    tuple(float(item) for item in values),
                    tuple(float(item) for item in energies),
                    temperature=temperature,
                    unit="hartree",
                    degeneracies=degeneracies,
                )
            )
        except AggregationError as exc:
            raise QuantityExpressionError(str(exc)) from exc
        dimension = inputs[0].dimension
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=canonical_unit_for_dimension(dimension),
            value=payload,
            unit=canonical_unit_for_dimension(dimension),
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation in {"linear_fit_slope", "linear_fit_intercept"}:
        if len(inputs) != 2:
            raise QuantityExpressionError(
                f"{operation} requires x and y inputs"
            )
        x, y = inputs
        x_values, y_values = _numeric(x), _numeric(y)
        if (
            x_values.ndim != 1
            or y_values.ndim != 1
            or x_values.shape != y_values.shape
            or x_values.size < 2
        ):
            raise QuantityExpressionError(
                "linear fit requires equal one-dimensional vectors of length >= 2"
            )
        x_centered = x_values - np.mean(x_values)
        denominator = float(np.dot(x_centered, x_centered))
        if denominator <= 0.0 or not math.isfinite(denominator):
            raise QuantityExpressionError(
                "linear fit requires non-constant finite x values"
            )
        slope = float(
            np.dot(x_centered, y_values - np.mean(y_values)) / denominator
        )
        intercept = float(np.mean(y_values) - slope * np.mean(x_values))
        if operation == "linear_fit_slope":
            payload = slope
            dimension = _subtract_dimensions(y.dimension, x.dimension)
        else:
            payload = intercept
            dimension = y.dimension
        unit = canonical_unit_for_dimension(dimension)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=dimension,
            evidence_ref=evidence_ref,
        )

    if operation == "distance":
        if len(inputs) != 2 or any(
            item.dimension != LENGTH for item in inputs
        ):
            raise QuantityExpressionError(
                "distance requires two length-coordinate vectors"
            )
        left, right = (_numeric(item) for item in inputs)
        if left.ndim != 1 or right.ndim != 1 or left.shape != right.shape:
            raise QuantityExpressionError(
                "distance requires equal one-dimensional coordinate vectors"
            )
        payload = internal_distance(left, right)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit="angstrom",
            value=payload,
            unit="angstrom",
            dimension=LENGTH,
            evidence_ref=evidence_ref,
        )

    if operation == "angle":
        if len(inputs) != 3 or any(
            item.dimension != LENGTH for item in inputs
        ):
            raise QuantityExpressionError(
                "angle requires three length-coordinate vectors"
            )
        first, center, last = (_numeric(item) for item in inputs)
        if any(vector.ndim != 1 for vector in (first, center, last)):
            raise QuantityExpressionError(
                "angle inputs must be coordinate vectors"
            )
        if first.shape != center.shape or center.shape != last.shape:
            raise QuantityExpressionError("angle coordinate shapes must match")
        left = first - center
        right = last - center
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        if denominator == 0.0:
            raise QuantityExpressionError(
                "angle is undefined for zero-length vectors"
            )
        radians = internal_angle(first, center, last)
        degrees = math.degrees(radians)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=degrees,
            source_unit="degree",
            value=radians,
            unit="radian",
            dimension=ANGLE,
            evidence_ref=evidence_ref,
        )

    if operation == "dihedral":
        # The third standard internal coordinate.  distance and angle were
        # owned and this was not, which leaves a torsion -- the coordinate a
        # rotational barrier is defined along -- to be rebuilt from cross
        # products and an atan2 the model would have to get the sign of right.
        if len(inputs) != 4 or any(
            item.dimension != LENGTH for item in inputs
        ):
            raise QuantityExpressionError(
                "dihedral requires four length-coordinate vectors, in bonded "
                "order a-b-c-d"
            )
        a, b, c, d = (_numeric(item) for item in inputs)
        if any(vector.ndim != 1 for vector in (a, b, c, d)):
            raise QuantityExpressionError(
                "dihedral inputs must be coordinate vectors"
            )
        if len({vector.shape for vector in (a, b, c, d)}) != 1:
            raise QuantityExpressionError(
                "dihedral coordinate shapes must match"
            )
        b1, b2, b3 = b - a, c - b, d - c
        norm = float(np.linalg.norm(b2))
        if norm == 0.0:
            raise QuantityExpressionError(
                "dihedral is undefined when the central atoms coincide"
            )
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        if (
            float(np.linalg.norm(n1)) == 0.0
            or float(np.linalg.norm(n2)) == 0.0
        ):
            raise QuantityExpressionError(
                "dihedral is undefined for three collinear atoms"
            )
        radians = internal_dihedral(a, b, c, d)
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=math.degrees(radians),
            source_unit="degree",
            value=radians,
            unit="radian",
            dimension=ANGLE,
            evidence_ref=evidence_ref,
        )

    if operation == "convert":
        if len(inputs) != 1 or not node.target_unit:
            raise QuantityExpressionError(
                "convert requires one input and target_unit"
            )
        source = inputs[0]
        display = convert_normalized_value(
            source.value, source.dimension, node.target_unit
        )
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=display,
            source_unit=node.target_unit,
            value=source.value,
            unit=source.unit,
            dimension=source.dimension,
            evidence_ref=evidence_ref,
        )

    raise QuantityExpressionError(
        f"operation is not implemented: {operation!r}"
    )


def _node_authored_constants(
    node: QuantityExpressionNodeV1,
) -> frozenset[ModelAuthoredConstantV1]:
    """Name every number this node contributes that no measurement produced."""

    found: list[ModelAuthoredConstantV1] = []
    if node.literal_value is not None:
        role = (
            "power_exponent" if node.operation == "power" else "literal_value"
        )
        found.append(
            ModelAuthoredConstantV1(
                node_id=node.node_id,
                role=role,
                value=_constant_text(node.literal_value),
            )
        )
    if node.scale_factor is not None:
        found.append(
            ModelAuthoredConstantV1(
                node_id=node.node_id,
                role="scale_factor",
                value=_constant_text(node.scale_factor),
            )
        )
    if node.extrapolation_exponent is not None:
        found.append(
            ModelAuthoredConstantV1(
                node_id=node.node_id,
                role="extrapolation_exponent",
                value=_constant_text(node.extrapolation_exponent),
            )
        )
    return frozenset(found)


def _constant_text(value: Any) -> str:
    """Render a model-supplied constant compactly and reproducibly."""

    if isinstance(value, (list, tuple)):
        rendered = ",".join(_constant_text(item) for item in value)
        text = f"[{rendered}]"
    else:
        text = repr(float(value))
    return text if len(text) <= 256 else text[:253] + "..."


def evaluate_quantity_expression(
    request: QuantityExpressionRequestV1,
) -> QuantityExpressionReceiptV1:
    """Evaluate a finite, topologically ordered expression DAG once."""

    request_sha256 = canonical_quantity_sha256(request)
    values: dict[str, QuantityValueV1] = {
        quantity.quantity_id: quantity for quantity in request.inputs
    }
    dependencies: dict[str, frozenset[str]] = {
        quantity.quantity_id: _receipt_dependencies(quantity.evidence_ref)
        for quantity in request.inputs
    }
    authored: dict[str, frozenset[ModelAuthoredConstantV1]] = {
        quantity.quantity_id: frozenset() for quantity in request.inputs
    }
    conventions: dict[str, frozenset[str]] = {
        quantity.quantity_id: frozenset() for quantity in request.inputs
    }
    arithmetic: dict[str, frozenset[str]] = {
        quantity.quantity_id: frozenset() for quantity in request.inputs
    }
    derived: list[QuantityValueV1] = []
    evidence_ref = f"expression:{request.expression_id}#{request_sha256}"
    for node in request.nodes:
        if node.node_id in values:
            raise QuantityExpressionError(
                f"expression node ID collides with an existing value: {node.node_id}"
            )
        value = _node_value(
            node=node, values=values, evidence_ref=evidence_ref
        )
        values[node.node_id] = value
        derived.append(value)
        if node.operation == "literal":
            sources: tuple[str, ...] = ()
            dependencies[node.node_id] = frozenset()
        elif node.operation == "ref":
            sources = (node.reference or node.input_ids[0],)
            dependencies[node.node_id] = dependencies[sources[0]]
        else:
            sources = tuple(node.input_ids)
            dependencies[node.node_id] = frozenset().union(
                *(dependencies[input_id] for input_id in sources)
            )
        authored[node.node_id] = frozenset().union(
            *(authored[source] for source in sources), frozenset()
        ) | _node_authored_constants(node)
        conventions[node.node_id] = frozenset().union(
            *(conventions[source] for source in sources), frozenset()
        ) | (
            {node.operation}
            if node.operation in CONVENTION_OPERATIONS
            else frozenset()
        )
        # Count nodes, not operation names: rebuilding a convention shows up
        # as many arithmetic nodes, and collapsing them by name would hide
        # exactly the thing worth seeing.
        arithmetic[node.node_id] = frozenset().union(
            *(arithmetic[source] for source in sources), frozenset()
        ) | (
            {node.node_id}
            if node.operation in ARITHMETIC_OPERATIONS
            else frozenset()
        )
    outputs = tuple(values[node_id] for node_id in request.output_node_ids)
    output_dependencies = tuple(
        QuantityExpressionOutputDependencyV1(
            output_id=node_id,
            source_receipt_sha256s=tuple(sorted(dependencies[node_id])),
            model_authored_constants=tuple(
                sorted(authored[node_id], key=lambda item: item.sort_key())
            ),
            convention_operations=tuple(sorted(conventions[node_id])),
            arithmetic_node_count=len(arithmetic[node_id]),
        )
        for node_id in request.output_node_ids
    )
    semantic_signature_sha256 = quantity_expression_semantic_signature(request)
    body = {
        "schema_version": "chemsmart.quantity-expression-receipt.v1",
        "expression_id": request.expression_id,
        "request_sha256": request_sha256,
        "semantic_signature_sha256": semantic_signature_sha256,
        "node_values": tuple(derived),
        "outputs": outputs,
        "output_dependencies": output_dependencies,
        "status": "derived",
    }
    return QuantityExpressionReceiptV1(
        **body, receipt_sha256=canonical_quantity_sha256(body)
    )


#: What a level says about the Hamiltonian a number was computed with.
#: Numerics -- grid, density fitting, convergence -- are not here: two
#: programs at one functional differ by ~1e-4 Eh in totals at their default
#: numerics and agree to ~2e-6 Eh at matched tight numerics (CUHK Slurm
#: 2149277/2149278), which is the same Hamiltonian computed twice.
LEVEL_IDENTITY_FIELDS = (
    "method",
    "basis",
    "basis_functions",
    "ecp_core_electrons",
    "dispersion",
    "solvation",
    "frozen_core",
    "response_method",
)
#: Level fields compared only between operands whose level states them: a
#: record minted before a reader stated the field, or a program with no
#: basis at all, says nothing about it, which is not a difference.
STATED_ONLY_LEVEL_FIELDS = ("basis_functions",)
#: Level fields that are facts about the molecule as well as the level,
#: and so are compared by what they mean rather than by equality: the
#: electrons each element's core potential replaced agree when every
#: element two operands share is treated alike (HI and H are one level
#: with 28 and 0 in total), and a frozen core agrees when the operands'
#: counts share one rule (HI freezes 4 orbitals, H none, under one rule).
MOLECULE_DEPENDENT_LEVEL_FIELDS = ("ecp_core_electrons", "frozen_core")
#: Level fields that describe how an excited root was computed, and so
#: are compared only between operands that are excited-root values: a
#: ground-state energy read from a TDA run and one read from a full TD-DFT
#: run are one SCF, while two excitations from those runs are two
#: approximations of one root (TDA lies above full TD-DFT, by 0.44 eV for
#: acrolein's bright pi->pi*, CUHK Slurm 2150076).  The manifold is not
#: here: a singlet and a triplet root of one reference subtract into a
#: gap by design, and which states were combined is an identity, not a
#: level.
EXCITED_ROOT_LEVEL_FIELDS = ("response_method",)


def expression_output_sources(
    request: QuantityExpressionRequestV1,
) -> dict[str, frozenset[tuple[str, str]]]:
    """The ``(extraction receipt, quantity id)`` pairs each output descends from.

    The same walk as the receipt-level dependencies, one grain finer: an
    extraction receipt can carry a reference value and an excited-root
    value side by side, and a level field that describes only one of them
    must be asked of the quantity actually consumed.
    """

    sources: dict[str, frozenset[tuple[str, str]]] = {}
    for quantity in request.inputs:
        receipts = _RECEIPT_REF.findall(quantity.evidence_ref)
        quantities = _QUANTITY_REF.findall(quantity.evidence_ref)
        sources[quantity.quantity_id] = frozenset(
            (receipt, name)
            for receipt in receipts[-1:]
            for name in quantities[-1:]
        )
    for node in request.nodes:
        if node.operation in {"literal", "constant"}:
            sources[node.node_id] = frozenset()
        elif node.operation == "ref":
            sources[node.node_id] = sources.get(
                node.reference or node.input_ids[0], frozenset()
            )
        else:
            sources[node.node_id] = frozenset().union(
                *(sources.get(item, frozenset()) for item in node.input_ids),
                frozenset(),
            )
    return {
        output_id: sources.get(output_id, frozenset())
        for output_id in request.output_node_ids
    }


def _level_identity(level: Mapping[str, Any]) -> dict[str, Any]:
    """The comparable identity of one producer's level record."""

    def _word(value: Any) -> Any:
        if value is None:
            return None
        return str(value).strip().lower() or None

    method = (
        level.get("functional")
        or level.get("ab_initio")
        or level.get("method")
    )
    basis = _word(level.get("basis"))
    model = _word(level.get("solvent_model"))
    cores = level.get("ecp_core_electrons")
    conventions = level.get("frozen_core_conventions")
    return {
        "method": _word(method),
        # def2-SVP and Gaussian's def2svp are one basis, and so are
        # 6-31G* and 6-31G(d): Pople's star is the parenthesised name.
        "basis": _canonical_basis_name(basis) if basis else None,
        "basis_functions": _word(level.get("basis_functions")),
        "ecp_core_electrons": (
            {str(key): int(value) for key, value in cores.items()}
            if isinstance(cores, Mapping)
            else None
        ),
        "dispersion": _word(level.get("dispersion")),
        "solvation": (
            f"{model}:{_word(level.get('solvent')) or ''}" if model else "gas"
        ),
        "frozen_core": (
            (level.get("frozen_core"), tuple(conventions))
            if conventions is not None
            else level.get("frozen_core")
        ),
        "response_method": _word(level.get("response_method")),
    }


#: The 6-31G and 6-311G names whose star notation is a synonym by
#: definition: one star is (d) on the heavy atoms, two add (p) on
#: hydrogen.  A session wrote Gaussian's 6-31g(d) and PySCF's 6-31g* for
#: one basis set in one goal (R10 q12 g2-nh3).
_POPLE_STAR = re.compile(r"(6311|631)(\+{0,2})g(\*{1,2})")


def _canonical_basis_name(basis: str) -> str:
    """One spelling for one basis name, as the level compares it."""

    word = basis.replace("-", "")
    match = _POPLE_STAR.fullmatch(word)
    if match is None:
        return word
    family, diffuse, stars = match.groups()
    return f"{family}{diffuse}g" + ("(d,p)" if stars == "**" else "(d)")


def _ecp_cores_differ(values: Mapping[str, Any]) -> bool:
    """Whether two operands treat one shared element with different cores."""

    per_element: dict[str, set[int]] = {}
    for mapping in values.values():
        if not isinstance(mapping, Mapping):
            continue
        for element, count in mapping.items():
            per_element.setdefault(element, set()).add(int(count))
    return any(len(counts) > 1 for counts in per_element.values())


def _frozen_cores_differ(values: Mapping[str, Any]) -> bool:
    """Whether operands' frozen cores follow no one rule.

    A value is ``(count, conventions)`` when its reader classified the
    count, and a bare count from a record minted before it did; bare
    counts compare as counts, among themselves.  A molecule with no core
    is consistent with every rule.
    """

    from chemsmart.analysis.result_readers import FROZEN_CORE_NO_CORE

    classified = [
        set(value[1])
        for value in values.values()
        if isinstance(value, tuple) and FROZEN_CORE_NO_CORE not in value[1]
    ]
    bare = {value for value in values.values() if not isinstance(value, tuple)}
    if len(bare) > 1:
        return True
    if classified and not set.intersection(*classified):
        return True
    return False


def expression_level_observations(
    receipt: QuantityExpressionReceiptV1,
    levels_by_receipt: Mapping[str, Mapping[str, Any] | None],
    *,
    request: QuantityExpressionRequestV1 | None = None,
    provenance_by_receipt: Mapping[str, Mapping[str, str]] | None = None,
    source_extractions: Mapping[str, Iterable[str]] | None = None,
) -> tuple[dict[str, Any], ...]:
    """Say when an output combines numbers computed at different levels.

    For each output of ``receipt``, the producer levels of the extraction
    receipts it descends from are compared field by field over
    ``LEVEL_IDENTITY_FIELDS``.  A field on which they differ is reported
    with each receipt's value; a receipt whose producer states no level is
    reported as unstated rather than assumed equal.  An observation, never
    a refusal: a composite method mixes levels on purpose, and a high-level
    single point on a low-level geometry is an ordinary protocol -- the
    number stands and the reader is told what it is made of.

    An operand that is an earlier expression's output stands for the
    extraction receipts that output descends from, which
    ``source_extractions`` names per expression receipt: a difference of
    two per-program bond energies combines every result both energies
    came from, not two unlevelled numbers (R10 q12: both live goals built
    their cross-program comparison that way, and nothing was compared).

    ``EXCITED_ROOT_LEVEL_FIELDS`` are compared only between the receipts
    whose consumed quantities are excited-root values, which needs the
    ``request`` (which quantity of which receipt feeds each output) and
    each receipt's electronic provenance per quantity; without them the
    response is not compared, rather than compared on the wrong numbers.
    """

    output_sources = (
        expression_output_sources(request) if request is not None else {}
    )
    provenance_by_receipt = provenance_by_receipt or {}
    source_extractions = source_extractions or {}
    observations: list[dict[str, Any]] = []
    for dependency in receipt.output_dependencies:
        sources = tuple(
            sorted(
                {
                    leaf
                    for digest in dependency.source_receipt_sha256s
                    for leaf in (source_extractions.get(digest) or (digest,))
                }
            )
        )
        if len(sources) < 2:
            continue
        stated = {
            digest: _level_identity(levels_by_receipt[digest])
            for digest in sources
            if levels_by_receipt.get(digest)
        }
        unstated = tuple(digest for digest in sources if digest not in stated)
        if len(stated) < 2:
            continue
        excited = {
            digest
            for digest, quantity_id in output_sources.get(
                dependency.output_id, frozenset()
            )
            if (provenance_by_receipt.get(digest) or {}).get(quantity_id)
            == "excited_root"
        }
        differing = {}
        for field in LEVEL_IDENTITY_FIELDS:
            compared = (
                {
                    digest: stated[digest]
                    for digest in excited
                    if digest in stated
                }
                if field in EXCITED_ROOT_LEVEL_FIELDS
                else stated
            )
            if field in STATED_ONLY_LEVEL_FIELDS:
                compared = {
                    digest: identity
                    for digest, identity in compared.items()
                    if identity[field] is not None
                }
            if len(compared) < 2:
                continue
            values = {
                digest: identity[field]
                for digest, identity in compared.items()
            }
            if field in MOLECULE_DEPENDENT_LEVEL_FIELDS:
                differs = {
                    "ecp_core_electrons": _ecp_cores_differ,
                    "frozen_core": _frozen_cores_differ,
                }[field](values)
            else:
                differs = len(set(values.values())) > 1
            if differs:
                differing[field] = {
                    digest[:12]: (
                        list(value) if isinstance(value, tuple) else value
                    )
                    for digest, value in sorted(values.items())
                }
        if not differing:
            continue
        observations.append(
            {
                "kind": "operands_at_different_levels",
                "output_id": dependency.output_id,
                "differing_fields": differing,
                "receipts_without_level": tuple(d[:12] for d in unstated),
            }
        )
    return tuple(observations)


#: What a free energy, an enthalpy or an entropy means beyond its unit.
#: Two Gibbs energies in hartree subtract without complaint whether one is
#: harmonic and the other Grimme's, one at 1 atm and the other at 1 mol/L,
#: or one printed by a program and the other derived by the host; each of
#: those moves a reaction free energy by up to kcal/mol.
THERMOCHEMICAL_CONVENTION_FIELDS = (
    "source",
    "treatment",
    "temperature_k",
    "standard_state",
    "frequency_scale_factor",
    "masses",
)

#: Receipt quantities whose value depends on the translational standard
#: state: those carrying an entropy.
_STANDARD_STATE_QUANTITIES = frozenset(
    {
        "entropy",
        "entropy_times_temperature",
        "gibbs_free_energy",
        "thermal_gibbs_correction",
        "quasi_harmonic_entropy",
        "quasi_harmonic_entropy_times_temperature",
        "quasi_harmonic_gibbs_free_energy",
        "quasi_harmonic_thermal_gibbs_correction",
    }
)
#: Receipt quantities no thermochemical convention changes.
_CONVENTION_FREE_QUANTITIES = frozenset(
    {"electronic_energy", "temperature", "pressure", "near_zero_mode_count"}
)
#: Selectors under which a program's own printed thermochemistry is read.
_PRINTED_THERMOCHEMISTRY_SELECTORS = frozenset(
    {"gibbs_free_energy", "entropy_times_temperature"}
)


def thermochemical_convention(
    receipt: Any, quantity_id: str
) -> dict[str, Any] | None:
    """The convention one receipt's quantity was computed under, or None.

    A host thermochemistry receipt states its conditions and treatment in
    typed fields; which of them a quantity carries depends on the quantity
    -- a Grimme receipt's ``gibbs_free_energy`` is harmonic and its
    ``quasi_harmonic_gibbs_free_energy`` is Grimme's, a zero-point energy
    has no temperature and an enthalpy no standard state.  A program's
    printed thermochemistry read through an extraction receipt is that
    program's own convention and is named as such.  Anything else --
    electronic energies, geometries, derived expressions -- is None.
    """

    quantity_id = str(quantity_id)
    if getattr(receipt, "engine_id", None) is not None:
        if quantity_id in _CONVENTION_FREE_QUANTITIES:
            return None
        method = str(getattr(receipt, "entropy_method", "rrho") or "rrho")
        entropy_cutoff = getattr(receipt, "entropy_cutoff_cm1", None)
        enthalpy_cutoff = getattr(receipt, "enthalpy_cutoff_cm1", None)
        alpha = getattr(receipt, "alpha", 4)
        if method == "grimme":
            entropy_treatment = (
                f"grimme entropy ({entropy_cutoff:g} cm-1, alpha {alpha})"
            )
        elif method == "truhlar":
            entropy_treatment = f"truhlar entropy ({entropy_cutoff:g} cm-1)"
        else:
            entropy_treatment = "harmonic entropy"
        enthalpy_treatment = (
            f"head-gordon enthalpy ({enthalpy_cutoff:g} cm-1, alpha {alpha})"
            if enthalpy_cutoff is not None
            else "harmonic enthalpy"
        )
        if quantity_id in {
            "quasi_harmonic_entropy",
            "quasi_harmonic_entropy_times_temperature",
        }:
            treatment = entropy_treatment
        elif quantity_id == "quasi_harmonic_enthalpy":
            treatment = enthalpy_treatment
        elif quantity_id.startswith("quasi_harmonic_"):
            treatment = f"{entropy_treatment}; {enthalpy_treatment}"
        else:
            treatment = "harmonic (RRHO)"
        concentration = getattr(receipt, "concentration_mol_l", None)
        standard_state = None
        if quantity_id in _STANDARD_STATE_QUANTITIES:
            standard_state = (
                f"{concentration:g} mol/L"
                if concentration is not None
                else f"ideal gas at {getattr(receipt, 'pressure_atm', 1.0):g} atm"
            )
        return {
            "source": "host derivation",
            "treatment": treatment,
            "temperature_k": (
                None
                if quantity_id == "zero_point_energy"
                else getattr(receipt, "temperature_k", None)
            ),
            "standard_state": standard_state,
            "frequency_scale_factor": getattr(
                receipt, "frequency_scale_factor", 1.0
            ),
            "masses": (
                "natural-abundance weighted"
                if getattr(receipt, "use_weighted_mass", False)
                else "most abundant"
            ),
        }
    bindings = dict(getattr(receipt, "selector_bindings", ()) or ())
    selector = bindings.get(quantity_id)
    if selector in _PRINTED_THERMOCHEMISTRY_SELECTORS:
        from chemsmart.analysis.result_readers import (
            PRINTED_THERMOCHEMISTRY_CONVENTIONS,
        )

        program = str(getattr(receipt, "program", "") or "program")
        return {
            "source": f"printed by {program}",
            "treatment": PRINTED_THERMOCHEMISTRY_CONVENTIONS.get(
                program, f"{program}'s own"
            ),
            "temperature_k": None,
            "standard_state": None,
            "frequency_scale_factor": None,
            "masses": None,
        }
    return None


def expression_thermochemical_convention_observations(
    request: QuantityExpressionRequestV1,
    conventions_by_input: Mapping[str, Mapping[str, Any] | None],
) -> tuple[dict[str, Any], ...]:
    """Say when an output combines thermochemistry under different conventions.

    Each output is followed back through the expression DAG to the inputs
    it reads, and the conventions of those that are thermochemical are
    compared field by field over ``THERMOCHEMICAL_CONVENTION_FIELDS``; a
    field one side does not carry (a zero-point energy's temperature) is
    not a difference.  An observation, never a refusal: the spread between
    a harmonic and a quasi-harmonic Gibbs energy is a measurement worth
    making, and a composite free energy is built from parts on purpose.
    """

    reads: dict[str, frozenset[str]] = {
        quantity.quantity_id: frozenset({quantity.quantity_id})
        for quantity in request.inputs
    }
    for node in request.nodes:
        if node.operation in {"literal", "constant"}:
            reads[node.node_id] = frozenset()
        elif node.operation == "ref":
            reads[node.node_id] = reads.get(
                node.reference or node.input_ids[0], frozenset()
            )
        else:
            reads[node.node_id] = frozenset().union(
                *(reads.get(name, frozenset()) for name in node.input_ids)
            )
    observations: list[dict[str, Any]] = []
    for output_id in request.output_node_ids:
        stated = {
            input_id: conventions_by_input[input_id]
            for input_id in sorted(reads.get(output_id, frozenset()))
            if conventions_by_input.get(input_id)
        }
        if len(stated) < 2:
            continue
        differing = {}
        for field in THERMOCHEMICAL_CONVENTION_FIELDS:
            values = {
                input_id: convention.get(field)
                for input_id, convention in stated.items()
                if convention.get(field) is not None
            }
            if len(set(values.values())) > 1:
                differing[field] = values
        if not differing:
            continue
        observations.append(
            {
                "kind": "operands_at_different_thermochemical_conventions",
                "output_id": output_id,
                "differing_fields": differing,
                "meaning": (
                    "this output combines thermochemical quantities "
                    "computed under different conventions ("
                    + ", ".join(sorted(differing))
                    + "); the number stands as computed -- say whether "
                    "the difference is the measurement or a mismatch"
                ),
            }
        )
    return tuple(observations)


@dataclass(frozen=True)
class ExpressionOperandV1:
    """What the host knows about one number an expression reads.

    ``name`` is the quantity's meaning -- the selector an extraction bound
    it to, or a thermochemistry receipt's quantity id -- which
    ``result_quantities.ENERGY_KINDS`` reads as a kind. ``species`` is
    ``(formula, charge, multiplicity)`` of the structure the result
    describes, read by its reader; ``structure`` the digest of that result;
    ``not_stationary`` the sentence ``structure_stationarity`` says when
    the result's structure is shown not to be a stationary point, and
    empty otherwise. ``distances`` is the geometry the number belongs to --
    the sorted interatomic distances of the structural state its selector
    declares (``result_quantities.geometry_of_selector``) -- and None where
    the host cannot say which structure that is.
    """

    name: str
    species: tuple[str, Any, Any] | None = None
    structure: str = ""
    not_stationary: str = ""
    distances: tuple[float, ...] | None = None


#: Two numbers describe one geometry when their structures' sorted
#: interatomic distances agree to this (Angstrom). Measured on archived
#: results: a converged ORCA saddle and the single points run on the
#: geometry saved from it agree to 1.8e-4 A (ax41 po3, esterc4 and esterc5);
#: one stationary point located by two methods differs by the methods'
#: difference (H2O2 held at 90 deg by Gaussian and ORCA B3LYP: 1.01e-3 A;
#: water at MP2 and at B3LYP: 0.014 A), and that is two geometries.
ONE_GEOMETRY_ANGSTROM = 1e-3


def interatomic_distances(
    positions: Any, limit: int = 300
) -> tuple[float, ...] | None:
    """The sorted interatomic distances of one structure, in Angstrom.

    They do not move with orientation or atom order, so two results
    describe one geometry exactly when these agree. An atom has none (and
    every atom is one geometry); None when the positions are unread or the
    structure has more than ``limit`` atoms.
    """

    if positions is None:
        return None
    try:
        array = np.asarray(positions, dtype=float)
    except (TypeError, ValueError):
        return None
    if array.ndim != 2 or array.shape[1] != 3 or not 0 < len(array) <= limit:
        return None
    upper = np.triu_indices(len(array), 1)
    separations = np.linalg.norm(
        array[:, None, :] - array[None, :, :], axis=-1
    )
    return tuple(float(value) for value in np.sort(separations[upper]))


def geometry_difference(
    first: tuple[float, ...] | None, second: tuple[float, ...] | None
) -> float | None:
    """The largest difference between two structures' sorted distances.

    None when either geometry is unknown or the two have different numbers
    of atoms.
    """

    if first is None or second is None or len(first) != len(second):
        return None
    if not first:
        return 0.0
    return float(max(abs(a - b) for a, b in zip(first, second)))


def hill_formula(symbols: Iterable[str]) -> str:
    """Carbon first, hydrogen second, the rest alphabetical (Hill)."""

    counts = Counter(str(symbol) for symbol in symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(item for item in counts if item not in {"C", "H"})
    else:
        order = sorted(counts)
    return "".join(
        element + (str(counts[element]) if counts[element] != 1 else "")
        for element in order
    )


def _formula_counts(formula: str) -> Counter:
    counts: Counter = Counter()
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", formula or ""):
        counts[element] += int(number) if number else 1
    return counts


def _scalar_value(value: Any) -> float | None:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return None
    return float(array) if array.ndim == 0 else None


def _expression_input_sources(
    request: QuantityExpressionRequestV1,
) -> dict[str, tuple[str, str]]:
    """input id -> (receipt digest, quantity id) read from its evidence."""

    sources = {}
    for quantity in request.inputs:
        receipts = _RECEIPT_REF.findall(quantity.evidence_ref)
        quantities = _QUANTITY_REF.findall(quantity.evidence_ref)
        if receipts and quantities:
            sources[quantity.quantity_id] = (receipts[-1], quantities[-1])
    return sources


def _expression_linear_terms(
    request: QuantityExpressionRequestV1,
    receipt: QuantityExpressionReceiptV1,
    operand_for: Any,
    depth: int = 0,
) -> dict[str, list[tuple[float, str, Any]] | None]:
    """Each output as a linear combination of the numbers it reads.

    A term is ``(coefficient, label, operand)`` where the operand is an
    ``ExpressionOperandV1`` or ``("literal"|"constant", text)`` for a
    number the session or the registry supplied. ``None`` where the
    output is not linear in its operands (a product, a logarithm, a
    vector reduction). An operand that is an earlier expression's output
    is expanded through that expression's own request, which
    ``operand_for`` returns as ``("expression", request, receipt)``.
    """

    from chemsmart.analysis.result_quantities import energy_kind

    node_values = {
        item.quantity_id: _scalar_value(item.value)
        for item in receipt.node_values
    }
    sources = _expression_input_sources(request)
    lin: dict[str, list[tuple[float, str, Any]] | None] = {}
    single: dict[str, Any] = {}
    for quantity in request.inputs:
        node_values[quantity.quantity_id] = _scalar_value(quantity.value)
        source = sources.get(quantity.quantity_id)
        facts = operand_for(*source) if source else None
        if isinstance(facts, ExpressionOperandV1):
            single[quantity.quantity_id] = facts
            lin[quantity.quantity_id] = (
                [(1.0, quantity.quantity_id, facts)]
                if node_values[quantity.quantity_id] is not None
                else None
            )
        elif (
            isinstance(facts, tuple)
            and len(facts) == 3
            and facts[0] == "expression"
            and depth < 8
        ):
            nested = _expression_linear_terms(
                facts[1], facts[2], operand_for, depth + 1
            ).get(source[1])
            lin[quantity.quantity_id] = (
                None
                if nested is None
                else [
                    (coefficient, f"{quantity.quantity_id}/{label}", operand)
                    for coefficient, label, operand in nested
                ]
            )
        else:
            lin[quantity.quantity_id] = None
    for node in request.nodes:
        name = node.node_id
        inputs = list(node.input_ids)
        operation = node.operation
        parts = [lin.get(item) for item in inputs]
        if operation == "ref":
            source = node.reference or (inputs[0] if inputs else "")
            single[name] = single.get(source)
            lin[name] = lin.get(source) if not node.indices else None
            # A null keeps a whole axis: the selection is a vector, not one
            # element, and a vector is no linear term.
            if node.indices and None not in node.indices and source in single:
                facts = single[source]
                kind = (
                    energy_kind(facts.name)
                    if isinstance(facts, ExpressionOperandV1)
                    else None
                )
                if kind is not None and kind.indexes_structures:
                    # One point of a scan, a path or an optimisation is its
                    # own structure, whose geometry the result's reached
                    # positions are not.
                    facts = replace(
                        facts,
                        structure=f"{facts.structure}{list(node.indices)}",
                        distances=None,
                    )
                lin[name] = [(1.0, f"{source}{list(node.indices)}", facts)]
        elif operation == "convert":
            single[name] = single.get(inputs[0]) if inputs else None
            lin[name] = parts[0] if parts else None
        elif operation in {"add", "subtract", "sum", "mean"}:
            if not parts or any(part is None for part in parts):
                lin[name] = None
                continue
            if operation == "subtract":
                weights = [1.0, -1.0]
            elif operation == "mean":
                weights = [1.0 / len(parts)] * len(parts)
            else:
                weights = [1.0] * len(parts)
            lin[name] = [
                (weight * coefficient, label, operand)
                for weight, part in zip(weights, parts)
                for coefficient, label, operand in part
            ]
        elif operation == "scale":
            lin[name] = (
                None
                if not parts or parts[0] is None or node.scale_factor is None
                else [
                    (float(node.scale_factor) * coefficient, label, operand)
                    for coefficient, label, operand in parts[0]
                ]
            )
        elif operation == "abs":
            source_value = node_values.get(inputs[0]) if inputs else None
            lin[name] = (
                None
                if not parts or parts[0] is None or source_value is None
                else [
                    (
                        (-1.0 if source_value < 0 else 1.0) * coefficient,
                        label,
                        operand,
                    )
                    for coefficient, label, operand in parts[0]
                ]
            )
        elif operation in {"min", "max"}:
            # The result is one of its operands; which one is read from
            # the values the evaluation recorded.
            chosen = [
                item
                for item in inputs
                if node_values.get(item) is not None
                and node_values.get(item) == node_values.get(name)
            ]
            lin[name] = lin.get(chosen[0]) if chosen else None
        elif operation == "harmonic_zero_point_energy":
            facts = single.get(inputs[0]) if inputs else None
            lin[name] = (
                [
                    (
                        1.0,
                        name,
                        replace(facts, name="zero_point_energy"),
                    )
                ]
                if isinstance(facts, ExpressionOperandV1)
                else None
            )
        elif operation in {"literal", "constant"}:
            text = (
                str(node.constant_name)
                if operation == "constant"
                else f"{node.literal_value} {node.literal_unit}"
            )
            lin[name] = [(1.0, name, (operation, text))]
        else:
            lin[name] = None
    return {output: lin.get(output) for output in request.output_node_ids}


def _expression_reached_operands(
    request: QuantityExpressionRequestV1, operand_for: Any
) -> dict[str, dict[str, ExpressionOperandV1]]:
    """output id -> the host-known operands its value was computed from."""

    sources = _expression_input_sources(request)
    reads: dict[str, dict[str, ExpressionOperandV1]] = {}
    for quantity in request.inputs:
        source = sources.get(quantity.quantity_id)
        facts = operand_for(*source) if source else None
        reads[quantity.quantity_id] = (
            {quantity.quantity_id: facts}
            if isinstance(facts, ExpressionOperandV1)
            else {}
        )
    for node in request.nodes:
        if node.operation in {"literal", "constant"}:
            reads[node.node_id] = {}
        elif node.operation == "ref":
            reads[node.node_id] = dict(
                reads.get(node.reference or node.input_ids[0], {})
            )
        else:
            merged: dict[str, ExpressionOperandV1] = {}
            for item in node.input_ids:
                merged.update(reads.get(item, {}))
            reads[node.node_id] = merged
    return {
        output: reads.get(output, {}) for output in request.output_node_ids
    }


def _species_text(species: tuple[str, Any, Any]) -> str:
    formula, charge, _multiplicity = species
    if not isinstance(charge, (int, float)) or int(charge) == 0:
        return formula
    magnitude = abs(int(charge))
    return (
        formula
        + (str(magnitude) if magnitude > 1 else "")
        + ("+" if charge > 0 else "-")
    )


def _structures_of_one_composition(
    output_id: str,
    by_structure: Mapping[tuple, list[float]],
    structure_labels: Mapping[tuple, set[str]],
    structure_facts: Mapping[tuple, tuple[Any, str]],
    energy_layers: tuple[str, ...],
    layer_blocks: Mapping[tuple[int, ...], str],
) -> list[dict[str, Any]]:
    """Name what each structure of one composition contributes to an output.

    A conformer difference, a barrier within one formula, or a free energy
    assembled from one structure's electronic energy and another's thermal
    correction all cancel at the level of the formula, so the reaction
    statement says nothing about them. Per structure they do not. Numbers
    are grouped by the geometry they belong to (``distances``, read by the
    structural state their selectors declare), so a single point run on a
    frequency job's structure is one structure -- a composite, and nothing
    to say -- while a thermal part that enters beside another geometry's
    electronic energy is named, with how far apart the geometries are and
    whether the one whose energy it joins is a stationary point at all
    (R10 Q21 g1-hooh's first plan built "G at 90 deg" as G(cis saddle) +
    E(held 90) - E(cis saddle)). A number whose geometry the host cannot
    read is its own structure, and the sentence says so.
    """

    by_species: dict[tuple, list[tuple]] = {}
    for key in by_structure:
        by_species.setdefault(key[0], []).append(key)
    observations = []
    for species, keys in sorted(by_species.items(), key=str):
        clusters: list[dict[str, Any]] = []
        for key in sorted(
            keys, key=lambda item: sorted(structure_labels[item])
        ):
            distances, not_stationary = structure_facts.get(key, (None, ""))
            home = None
            for cluster in clusters if distances is not None else ():
                difference = geometry_difference(
                    cluster["distances"], distances
                )
                if difference is not None and (
                    difference <= ONE_GEOMETRY_ANGSTROM
                ):
                    home = cluster
                    break
            if home is None:
                home = {
                    "distances": distances,
                    "vector": [0.0] * len(energy_layers),
                    "labels": set(),
                    "not_stationary": "",
                }
                clusters.append(home)
            home["vector"] = [
                total + value
                for total, value in zip(home["vector"], by_structure[key])
            ]
            home["labels"] |= set(structure_labels[key])
            home["not_stationary"] = home["not_stationary"] or not_stationary
        entering = []
        for cluster in clusters:
            nonzero = [
                (index, value)
                for index, value in enumerate(cluster["vector"])
                if abs(value) > 1e-9
            ]
            if not nonzero:
                continue
            entering.append(
                {
                    "name": "[" + ", ".join(sorted(cluster["labels"])) + "]",
                    "coefficient": nonzero[0][1],
                    "block": tuple(index for index, _value in nonzero),
                    "values": [
                        (energy_layers[index], value)
                        for index, value in nonzero
                    ],
                    "uniform": len({round(value, 9) for _i, value in nonzero})
                    == 1,
                    "distances": cluster["distances"],
                    "not_stationary": cluster["not_stationary"],
                }
            )
        if len(entering) < 2:
            continue

        def _block_text(item: Mapping[str, Any]) -> str:
            if not item["uniform"]:
                return ", ".join(
                    f"{value:+g} x {layer}" for layer, value in item["values"]
                )
            named = layer_blocks.get(item["block"])
            if named:
                return f"{item['coefficient']:+g} x {named}"
            layers_text = "+".join(
                energy_layers[index] for index in item["block"]
            )
            if 0 in item["block"]:
                return f"{item['coefficient']:+g} x {layers_text}"
            return (
                f"{item['coefficient']:+g} x thermal part only ({layers_text})"
            )

        known = [item for item in entering if item["distances"] is not None]
        unread = [item for item in entering if item["distances"] is None]
        differences = [
            geometry_difference(first["distances"], second["distances"])
            for index, first in enumerate(known)
            for second in known[index + 1 :]
        ]
        largest = max(
            (value for value in differences if value is not None),
            default=None,
        )
        meaning = (
            f"this output combines {len(entering)} structures of one "
            f"composition ({_species_text(species)}): "
            + "; ".join(
                f"{item['name']} {_block_text(item)}" for item in entering
            )
        )
        if largest is not None:
            meaning += (
                f"; their geometries differ by up to {largest:.3g} A in "
                "interatomic distance"
            )
        if unread:
            meaning += (
                "; the host could not read which geometry "
                + ", ".join(item["name"] for item in unread)
                + (" describes" if len(unread) == 1 else " describe")
                + ", so each is counted as its own structure"
            )
        thermal_only = [item for item in entering if 0 not in item["block"]]
        # An energy entering bare: a structure that brings no thermal part
        # of its own, so the one beside it is borrowed.
        bare = [item for item in entering if item["block"] == (0,)]
        borrowed = bool(thermal_only and bare)
        judgments = []
        if borrowed:
            judgment = (
                "the thermal part of "
                + ", ".join(item["name"] for item in thermal_only)
                + " stands beside the electronic energy of "
                + ", ".join(item["name"] for item in bare)
            )
            apart = [
                geometry_difference(first["distances"], second["distances"])
                for first in thermal_only
                for second in bare
            ]
            if apart and all(value is not None for value in apart):
                judgment += (
                    f", another geometry ({max(apart):.3g} A apart): a "
                    "composite across levels when both are one stationary "
                    "point located by two methods, a correction borrowed "
                    "from another structure when they are not"
                )
            else:
                judgment += (
                    ": a composite when both describe one structure, a "
                    "correction borrowed from another structure when they "
                    "do not"
                )
            judgments.append(judgment)
            for item in bare:
                if item["not_stationary"]:
                    judgments.append(
                        f"{item['name']} is {item['not_stationary']}; a "
                        "thermal part beside its energy is therefore the "
                        "free energy of no state -- at most the "
                        "approximation that another structure's thermal "
                        "correction holds at this one"
                    )
        elif len({item["block"] for item in entering}) > 1:
            judgments.append(
                "the structures enter at different layers, so the output "
                "is a change at no single layer"
            )
        for item in entering:
            if not item["uniform"]:
                judgments.append(
                    f"{item['name']}'s own layers enter with different "
                    "coefficients, so it contributes no state function of "
                    "that structure"
                )
        if judgments:
            meaning += " -- " + "; ".join(judgments)
        observations.append(
            {
                "kind": "structures_of_one_composition",
                "output_id": output_id,
                "species": _species_text(species),
                "structures": [
                    {
                        "operands": item["name"][1:-1],
                        "coefficient": item["coefficient"],
                        "layers": [layer for layer, _value in item["values"]],
                        "geometry_read": item["distances"] is not None,
                        "not_stationary": item["not_stationary"],
                    }
                    for item in entering
                ],
                "largest_geometry_difference_angstrom": largest,
                "borrowed_thermal_part": borrowed,
                "meaning": meaning,
            }
        )
    return observations


def expression_kind_observations(
    request: QuantityExpressionRequestV1,
    receipt: QuantityExpressionReceiptV1,
    operand_for: Any,
) -> tuple[dict[str, Any], ...]:
    """Say what each output is, where its operands' kinds decide it.

    Units say two numbers can be added; they cannot say whether the sum
    means anything. For each output this reads the kind of every energy
    it was computed from (``result_quantities.ENERGY_KINDS``) and states:

    - ``output_is_an_orbital_rotation_curvature``: an output built from
      SCF stability eigenvalues is a curvature, not an energy between
      states, and eigenvalues of different matrices compare in sign only
      (R10 Q13 dans: ``min`` over PySCF's internal, external and
      real -> complex roots was claimed as "the lowest eigenvalue" and a
      spread of two was read as "~3.2 kcal/mol");
    - ``orbital_energy_combined_with_a_state_energy``: a one-electron
      eigenvalue added to or subtracted from an energy of a state;
    - ``one_species_at_two_coefficients``: where the output is linear in
      its energies, one species whose own layers (E, ZPE, thermal, pV,
      -TS) enter with different coefficients -- the change of no reaction
      (R10 Q14 G2: ``[E(OH) - E(O) - E(H)] - ZPE(OH)`` delivered as a D0 of
      -105.90 kcal/mol);
    - ``reaction_the_output_measures``: otherwise, the reaction the
      coefficients describe, reactants to products, each species' layer
      block, and whether the atoms balance -- the sign convention of the
      number, stated by the host rather than left to be inferred;
    - ``structures_of_one_composition``: where two or more geometries of
      one formula enter (conformers, saddles, a held structure, one
      minimum located by two methods), each one's coefficient and layer
      block, how far apart they are, a thermal part that enters beside
      another geometry's electronic energy, and an energy joined by a
      thermal part at a structure shown not to be stationary. A single
      point on a frequency run's own structure is one geometry: a
      composite, stated by the reaction it enters, not here.

    Observations, never refusals: a curvature compared with zero is a
    legitimate question, and a number stands as computed. What the host
    adds is what the number is.
    """

    from chemsmart.analysis.result_quantities import (
        ENERGY_LAYER_BLOCKS,
        ENERGY_LAYERS,
        energy_kind,
    )

    reached = _expression_reached_operands(request, operand_for)
    terms = _expression_linear_terms(request, receipt, operand_for)
    observations: list[dict[str, Any]] = []
    for output_id in request.output_node_ids:
        operands = reached.get(output_id, {})
        kinds = {
            input_id: energy_kind(facts.name)
            for input_id, facts in operands.items()
            if energy_kind(facts.name) is not None
        }
        curvatures = {
            input_id: kind
            for input_id, kind in kinds.items()
            if kind.kind == "orbital_rotation_curvature"
        }
        if curvatures:
            normalisations = sorted(
                {kind.normalisation for kind in curvatures.values()}
            )
            others = sorted(
                {
                    kind.kind
                    for kind in kinds.values()
                    if kind.kind != "orbital_rotation_curvature"
                }
            )
            meaning = (
                "this output is built from eigenvalues of SCF stability "
                "matrices: curvatures of the energy along rotations of the "
                "orbitals, not energy differences between states -- a "
                "hartree here is no gap of any kind, and it is never a "
                "kcal/mol of anything"
            )
            if len(normalisations) > 1:
                meaning += (
                    f"; its operands come from {len(normalisations)} "
                    "different matrices, whose signs compare (each says "
                    "whether its own rotation space holds a lower "
                    "solution) and whose magnitudes do not"
                )
            if others:
                meaning += (
                    "; it also combines them with "
                    + ", ".join(item.replace("_", " ") for item in others)
                    + ", which are not curvatures"
                )
            observations.append(
                {
                    "kind": "output_is_an_orbital_rotation_curvature",
                    "output_id": output_id,
                    "operands": {
                        input_id: operands[input_id].name
                        for input_id in sorted(curvatures)
                    },
                    "normalisations": normalisations,
                    "magnitudes_comparable": len(normalisations) == 1
                    and not others,
                    "meaning": meaning,
                }
            )
        linear = terms.get(output_id)
        if not linear:
            continue
        linear_kinds = {
            label: energy_kind(operand.name)
            for _coefficient, label, operand in linear
            if isinstance(operand, ExpressionOperandV1)
            and energy_kind(operand.name) is not None
        }
        orbital = sorted(
            label
            for label, kind in linear_kinds.items()
            if kind.kind == "orbital_energy"
        )
        stated = sorted(
            label
            for label, kind in linear_kinds.items()
            if kind.kind in {"state_energy", "correction", "excitation_energy"}
        )
        if orbital and stated:
            observations.append(
                {
                    "kind": "orbital_energy_combined_with_a_state_energy",
                    "output_id": output_id,
                    "orbital_operands": orbital,
                    "state_operands": stated,
                    "meaning": (
                        "an orbital eigenvalue is a one-electron energy; "
                        "adding it to or subtracting it from the energy of "
                        "a state is no energy difference between states"
                    ),
                }
            )
        # A zero-point or thermal part built from the modes of a structure
        # shown not to be stationary is the curvature there, which the
        # thermochemistry derivation refuses to call a free energy; the
        # same modes reached through arithmetic are named, not refused.
        unstationary = sorted(
            {
                (label, operand.not_stationary)
                for _coefficient, label, operand in linear
                if isinstance(operand, ExpressionOperandV1)
                and operand.not_stationary
                and energy_kind(operand.name) is not None
                and any(
                    layer > 0 for layer in energy_kind(operand.name).layers
                )
            }
        )
        if unstationary:
            observations.append(
                {
                    "kind": "vibrational_energy_of_a_structure_not_stationary",
                    "output_id": output_id,
                    "operands": dict(unstationary),
                    "meaning": (
                        "the zero-point or thermal part of this output comes "
                        "from harmonic modes of a structure that is not a "
                        "stationary point ("
                        + "; ".join(reason for _label, reason in unstationary)
                        + "): it is the curvature there, not the energy of "
                        "any state, and a free energy or zero-point-corrected "
                        "energy built on it describes no state"
                    ),
                }
            )
        layers: dict[tuple, list[float]] = {}
        placed: dict[tuple, list[tuple[float, str, str]]] = {}
        # The same accounting per structure, within each species: structures
        # of one composition (conformers, saddles, a held structure) all
        # share a formula, and a species-level sum cannot say whose
        # thermal part stands beside whose electronic energy.
        by_structure: dict[tuple, list[float]] = {}
        structure_labels: dict[tuple, set[str]] = {}
        structure_facts: dict[tuple, tuple[Any, str]] = {}
        authored = []
        for coefficient, label, operand in linear:
            if not isinstance(operand, ExpressionOperandV1):
                authored.append(f"{coefficient:+g} x {operand[1]} ({label})")
                continue
            kind = energy_kind(operand.name)
            if kind is None or not kind.layers or operand.species is None:
                continue
            vector = layers.setdefault(
                operand.species, [0.0] * len(ENERGY_LAYERS)
            )
            key = (operand.species, operand.structure, operand.distances)
            per_structure = by_structure.setdefault(
                key, [0.0] * len(ENERGY_LAYERS)
            )
            for layer in kind.layers:
                vector[layer] += coefficient * kind.sign
                per_structure[layer] += coefficient * kind.sign
            structure_labels.setdefault(key, set()).add(label)
            structure_facts[key] = (operand.distances, operand.not_stationary)
            placed.setdefault(operand.species, []).append(
                (coefficient, label, operand.name)
            )
        if not layers:
            continue
        observations.extend(
            _structures_of_one_composition(
                output_id,
                by_structure,
                structure_labels,
                structure_facts,
                ENERGY_LAYERS,
                ENERGY_LAYER_BLOCKS,
            )
        )
        inconsistent = {
            species: vector
            for species, vector in layers.items()
            if len({round(value, 9) for value in vector if abs(value) > 1e-9})
            > 1
        }
        if inconsistent:
            for species, vector in sorted(inconsistent.items(), key=str):
                entering = ", ".join(
                    f"{value:+g} x {ENERGY_LAYERS[index]}"
                    for index, value in enumerate(vector)
                    if abs(value) > 1e-9
                )
                observations.append(
                    {
                        "kind": "one_species_at_two_coefficients",
                        "output_id": output_id,
                        "species": _species_text(species),
                        "layers": {
                            ENERGY_LAYERS[index]: value
                            for index, value in enumerate(vector)
                            if abs(value) > 1e-9
                        },
                        "operands": [
                            f"{coefficient:+g} x {name} ({label})"
                            for coefficient, label, name in placed[species]
                        ],
                        "meaning": (
                            f"{_species_text(species)} enters this output "
                            f"as {entering}: the layers of one species' "
                            "energy carry different coefficients, so the "
                            "output is the energy change of no reaction and "
                            "no state function of that species -- a sign "
                            "was flipped on part of it, or a part belongs to "
                            "another step"
                        ),
                    }
                )
            continue
        sides: dict[str, list[str]] = {"reactants": [], "products": []}
        blocks = {}
        net_atoms: Counter = Counter()
        net_charge = 0.0
        electronic = False
        for species, vector in sorted(layers.items(), key=str):
            nonzero = [
                (index, value)
                for index, value in enumerate(vector)
                if abs(value) > 1e-9
            ]
            if not nonzero:
                continue
            coefficient = nonzero[0][1]
            block = tuple(index for index, _value in nonzero)
            blocks[_species_text(species)] = ENERGY_LAYER_BLOCKS.get(
                block, "+".join(ENERGY_LAYERS[index] for index in block)
            )
            text = (
                f"{abs(coefficient):g} "
                if abs(abs(coefficient) - 1) > 1e-9
                else ""
            ) + _species_text(species)
            sides["products" if coefficient > 0 else "reactants"].append(text)
            if 0 in block:
                electronic = True
                for element, count in _formula_counts(species[0]).items():
                    net_atoms[element] += coefficient * count
                if isinstance(species[1], (int, float)):
                    net_charge += coefficient * float(species[1])
        if not sides["reactants"] or not sides["products"]:
            continue
        imbalance = {
            element: round(value, 6)
            for element, value in net_atoms.items()
            if abs(value) > 1e-9
        }
        if not electronic:
            balance = (
                "no electronic energy enters: a difference of corrections "
                "between the species named"
            )
        elif not imbalance:
            balance = "the atoms balance"
        elif (
            set(imbalance) == {"H"} and abs(imbalance["H"] - net_charge) < 1e-9
        ):
            balance = (
                f"the atoms balance with {abs(imbalance['H']):g} bare "
                "proton(s), which carry no electronic energy (their "
                "thermochemistry does not vanish)"
            )
        elif authored:
            balance = (
                f"the atoms do not balance (net {imbalance}) except "
                "through number(s) the expression itself supplies: "
                + "; ".join(authored)
            )
        else:
            balance = (
                f"the atoms do not balance (net {imbalance}): this is not "
                "the energy change of any reaction"
            )
        reaction = (
            " + ".join(sides["reactants"])
            + " -> "
            + " + ".join(sides["products"])
        )
        treatments = sorted(set(blocks.values()))
        observations.append(
            {
                "kind": "reaction_the_output_measures",
                "output_id": output_id,
                "reaction": reaction,
                "treatments": blocks,
                "atoms_balance": electronic and not imbalance,
                "meaning": (
                    f"this output is the change of {reaction} (products "
                    "minus reactants), "
                    + (
                        f"in {treatments[0]}"
                        if len(treatments) == 1
                        else "with species at different treatments "
                        + ", ".join(f"{k} as {v}" for k, v in blocks.items())
                    )
                    + f"; {balance}"
                ),
            }
        )
    return tuple(observations)


def _quantity_text(value: float, unit: str) -> str:
    text = f"{float(value):.6g}"
    return text if unit in {"", "1"} else f"{text} {unit}"


def expression_extremum_observations(
    request: QuantityExpressionRequestV1,
    receipt: QuantityExpressionReceiptV1,
) -> tuple[dict[str, Any], ...]:
    """Name the runner-up beside every coordinate_at_maximum and _minimum.

    Where an extremum falls says nothing about how far ahead of the rest it
    was. Summed C=O participation picked the carbonyl stretches of both R9
    xtb g3 Hessians at shares of 0.96-0.98 while the next modes in rank --
    low-frequency motions loading the same two atoms -- held 0.36-0.63 (R11
    probe M, CUHK 2157069). Each selection's runner-up is read from the
    arrays the evaluation read, for every such node whether or not it is an
    output: a fact beside the receipt, never a verdict, and never inside
    it, so the receipt's digest stays what the arithmetic makes it.
    """

    known = {quantity.quantity_id: quantity for quantity in request.inputs}
    known.update({item.quantity_id: item for item in receipt.node_values})
    observations: list[dict[str, Any]] = []
    for node in request.nodes:
        if node.operation not in {
            "coordinate_at_maximum",
            "coordinate_at_minimum",
        }:
            continue
        operands = [known.get(item) for item in node.input_ids]
        if len(operands) != 2 or None in operands:
            continue
        series, measured_at = operands
        values = np.asarray(series.value, dtype=float).reshape(-1)
        coordinates = np.asarray(measured_at.value, dtype=float).reshape(-1)
        if values.size != coordinates.size or values.size < 2:
            continue
        largest = node.operation == "coordinate_at_maximum"
        pick = np.argmax if largest else np.argmin
        chosen = int(pick(values))
        rest = np.delete(np.arange(values.size), chosen)
        runner_up = int(rest[pick(values[rest])])
        separation = abs(float(values[chosen]) - float(values[runner_up]))

        def point(index: int) -> dict[str, Any]:
            return {
                "index": index,
                "coordinate": float(coordinates[index]),
                "value": float(values[index]),
            }

        rank, side = (
            ("largest", "lower") if largest else ("smallest", "higher")
        )
        observations.append(
            {
                "kind": "extremum_runner_up",
                "node_id": node.node_id,
                "operation": node.operation,
                "points": int(values.size),
                "selected": point(chosen),
                "runner_up": point(runner_up),
                "separation": separation,
                "coordinate_unit": measured_at.unit,
                "value_unit": series.unit,
                "meaning": (
                    f"{node.operation} returned "
                    f"{_quantity_text(coordinates[chosen], measured_at.unit)}"
                    f", where the {rank} of its {values.size} values falls "
                    f"({_quantity_text(values[chosen], series.unit)}); the "
                    f"next {rank} "
                    f"({_quantity_text(values[runner_up], series.unit)}) "
                    "falls at "
                    f"{_quantity_text(coordinates[runner_up], measured_at.unit)}"
                    f", {_quantity_text(separation, series.unit)} {side}"
                ),
            }
        )
    return tuple(observations)


def quantity_expression_receipt_from_record(
    record: Mapping[str, Any], *, receipt_sha256: str
) -> QuantityExpressionReceiptV1:
    """Rehydrate an expression receipt persisted by Runtime V2."""

    from chemsmart.analysis.result_quantities import quantity_value_from_record

    values = dict(record)
    values["node_values"] = tuple(
        quantity_value_from_record(item)
        for item in values.get("node_values") or ()
    )
    values["outputs"] = tuple(
        quantity_value_from_record(item)
        for item in values.get("outputs") or ()
    )
    values["output_dependencies"] = tuple(
        QuantityExpressionOutputDependencyV1(
            output_id=str(item["output_id"]),
            source_receipt_sha256s=tuple(
                item.get("source_receipt_sha256s") or ()
            ),
            model_authored_constants=tuple(
                ModelAuthoredConstantV1(
                    node_id=str(entry["node_id"]),
                    role=str(entry["role"]),
                    value=str(entry["value"]),
                )
                for entry in item.get("model_authored_constants") or ()
            ),
            convention_operations=tuple(
                item.get("convention_operations") or ()
            ),
            arithmetic_node_count=int(item.get("arithmetic_node_count") or 0),
        )
        for item in values.get("output_dependencies") or ()
    )
    return QuantityExpressionReceiptV1(**values, receipt_sha256=receipt_sha256)


__all__ = [
    "MAX_EXPRESSION_NODES",
    "MAX_NODE_INPUTS",
    "ARITHMETIC_OPERATIONS",
    "CONVENTION_OPERATIONS",
    "OPERATION_FAMILIES",
    "MODEL_AUTHORED_CONSTANT_ROLES",
    "OPERATION_DESCRIPTIONS",
    "ModelAuthoredConstantV1",
    "QuantityExpressionError",
    "QuantityExpressionNodeV1",
    "QuantityExpressionOutputDependencyV1",
    "QuantityExpressionReceiptV1",
    "QuantityExpressionRequestV1",
    "canonical_unit_for_dimension",
    "convert_normalized_value",
    "evaluate_quantity_expression",
    "ExpressionOperandV1",
    "expression_extremum_observations",
    "expression_kind_observations",
    "expression_level_observations",
    "expression_output_sources",
    "expression_thermochemical_convention_observations",
    "geometry_difference",
    "hill_formula",
    "interatomic_distances",
    "ONE_GEOMETRY_ANGSTROM",
    "EXCITED_ROOT_LEVEL_FIELDS",
    "LEVEL_IDENTITY_FIELDS",
    "THERMOCHEMICAL_CONVENTION_FIELDS",
    "thermochemical_convention",
    "normalize_numeric_value",
    "quantity_expression_semantic_signature",
    "quantity_expression_receipt_from_record",
]
