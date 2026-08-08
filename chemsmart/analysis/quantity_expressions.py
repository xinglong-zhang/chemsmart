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
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
from ase import units as ase_units

from chemsmart.analysis.result_quantities import (
    ANGLE,
    DIMENSIONLESS,
    ENERGY,
    ENTROPY,
    FREQUENCY,
    LENGTH,
    PRESSURE,
    TEMPERATURE,
    Dimension,
    QuantityContractError,
    QuantityValueV1,
    canonical_quantity_sha256,
    make_quantity_value,
)
from chemsmart.utils.constants import energy_conversion

MAX_EXPRESSION_NODES = 128
MAX_NODE_INPUTS = 64
_RECEIPT_REF = re.compile(r"(?:^|[;])receipt:([0-9a-f]{64})(?:$|[;])")
_QUANTITY_REF = re.compile(r"(?:^|[;])quantity:([^;]+)(?:$|[;])")
_SEMANTIC_ROLE_REF = re.compile(r"(?:^|[;])semantic-role:([^;]+)(?:$|[;])")

_OPERATIONS = frozenset(
    {
        "ref",
        "literal",
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
        "boltzmann_populations",
        "boltzmann_average",
        "imaginary_mode_count",
        "harmonic_zero_point_energy",
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
    "distance": "distance between two indexed coordinate vectors",
    "angle": "angle at the middle of three indexed coordinate vectors",
    "dihedral": (
        "signed torsion about the middle bond of four indexed coordinate "
        "vectors, in (-180, 180]. The third standard internal coordinate "
        "alongside distance and angle; do not rebuild it from cross products"
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
    "boltzmann_populations": (
        "normalized Boltzmann populations of a set of states. Takes the state "
        "energies and a temperature, in that order, optionally followed by a "
        "dimensionless vector of per-state degeneracies, and owns the "
        "weighting, the gas constant and the unit handling. Supply the "
        "degeneracies whenever states are multiply realizable -- an "
        "enantiomeric pair counts twice. Do not rebuild any of this from exp, "
        "scale, divide and sum"
    ),
    "boltzmann_average": (
        "Boltzmann-weighted average of a property. Takes the per-state values, "
        "the state energies and a temperature, in that order, optionally "
        "followed by the per-state degeneracies"
    ),
    "imaginary_mode_count": (
        "how many harmonic modes are imaginary, from one frequency vector, "
        "optionally followed by a frequency cutoff below which a near-zero "
        "mode is ignored. Owns the sign convention: a minimum has zero, a "
        "transition state exactly one. Do not rebuild it from comparisons"
    ),
    "harmonic_zero_point_energy": (
        "harmonic zero-point vibrational energy from one frequency vector in "
        "cm^-1, returned as a molar energy. Owns both the factor of one half "
        "and the spectroscopic h*c*N_A conversion, and skips imaginary modes. "
        "Use this instead of sum -> scale 0.5 -> convert: a wavenumber and a "
        "molar energy are different dimensions, so convert refuses that step"
    ),
}

if set(OPERATION_DESCRIPTIONS) != set(_OPERATIONS):  # pragma: no cover
    raise RuntimeError(
        "every expression operation must be described: "
        f"{sorted(set(_OPERATIONS) ^ set(OPERATION_DESCRIPTIONS))}"
    )

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
        "imaginary_mode_count",
        "distance",
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
_PLUMBING_OPERATIONS = frozenset({"ref", "literal", "convert"})

#: Everything else: general arithmetic and reductions.  A reported quantity
#: reached through many of these, and through no convention operation, was
#: assembled rather than computed by the toolkit.
ARITHMETIC_OPERATIONS = (
    frozenset(_OPERATIONS) - CONVENTION_OPERATIONS - (_PLUMBING_OPERATIONS)
)


class QuantityExpressionError(ValueError):
    """Raised when an expression is unsafe, ill-typed, or non-finite."""


def _add_dimensions(left: Dimension, right: Dimension) -> Dimension:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _subtract_dimensions(left: Dimension, right: Dimension) -> Dimension:
    return tuple(a - b for a, b in zip(left, right))  # type: ignore[return-value]


def canonical_unit_for_dimension(dimension: Dimension) -> str:
    known = {
        DIMENSIONLESS: "1",
        ENERGY: "hartree",
        LENGTH: "angstrom",
        TEMPERATURE: "K",
        ANGLE: "radian",
        FREQUENCY: "cm^-1",
        PRESSURE: "atm",
        ENTROPY: "hartree K^-1",
    }
    if dimension in known:
        return known[dimension]
    labels = ("hartree", "angstrom", "K", "radian", "cm^-1", "atm")
    terms = []
    for label, exponent in zip(labels, dimension):
        if exponent == 0:
            continue
        terms.append(label if exponent == 1 else f"{label}^{exponent}")
    return " ".join(terms) if terms else "1"


def _normalized_unit_key(unit: str) -> str:
    return (
        str(unit)
        .strip()
        .lower()
        .replace("·", " ")
        .replace("å", "angstrom")
        .replace("−", "-")
        .replace("⁻", "-")
        .replace("  ", " ")
    )


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
        "thz": (FREQUENCY, "cm^-1", 1.0e12 / 29_979_245_800.0),
        "atm": (PRESSURE, "atm", 1.0),
        "bar": (PRESSURE, "atm", 1.0 / 1.01325),
        "pa": (PRESSURE, "atm", 1.0 / 101_325.0),
    }
    try:
        return aliases[key]
    except KeyError as exc:
        raise QuantityExpressionError(f"unsupported unit: {unit!r}") from exc


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
    indices: tuple[int, ...] = ()
    literal_value: Any = None
    literal_unit: str = "1"
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
        if any(index < 0 for index in self.indices):
            raise QuantityContractError(
                "reference indices must be non-negative"
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
        semantic_role = (
            source_quantity_id
            if source_quantity_id
            and source_quantity_counts[source_quantity_id] == 1
            else (
                role_matches[-1]
                if role_matches
                else (
                    source_quantity_id
                    if source_quantity_id
                    else quantity.quantity_id
                )
            )
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
        selected = _numeric(source)
        try:
            for index in node.indices:
                selected = np.asarray(selected[index])
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

    try:
        inputs = [values[input_id] for input_id in node.input_ids]
    except KeyError as exc:
        raise QuantityExpressionError(
            f"node references an unavailable prior value: {exc.args[0]!r}"
        ) from exc

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
            if left_value.shape != right_value.shape:
                raise QuantityExpressionError(
                    f"{operation} requires identical numeric shapes"
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
        return make_quantity_value(
            quantity_id=node.node_id,
            source_value=payload,
            source_unit=unit,
            value=payload,
            unit=unit,
            dimension=ENERGY,
            evidence_ref=evidence_ref,
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
        # input.  For n-butane it moves the answer by twelve points.
        degeneracies = None
        if len(inputs) == expected + 1:
            candidate = inputs[-1]
            if candidate.dimension != DIMENSIONLESS:
                raise QuantityExpressionError(
                    f"the optional last input to {operation} is the "
                    "per-state degeneracies and must be dimensionless"
                )
            degeneracies = tuple(
                float(item) for item in _numeric(candidate).reshape(-1)
            )
            inputs = inputs[:-1]
        if len(inputs) != expected:
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
        temperature_value = inputs[-1]
        if temperature_value.dimension != TEMPERATURE:
            raise QuantityExpressionError(
                f"the last input to {operation} must be a temperature"
            )
        temperature = float(_numeric(temperature_value))
        energies = _numeric(inputs[expected - 2]).reshape(-1)
        if inputs[expected - 2].dimension != ENERGY:
            raise QuantityExpressionError(
                f"{operation} weights states by energy"
            )
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
        payload = float(np.linalg.norm(left - right))
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
        cosine = float(np.dot(left, right) / denominator)
        radians = math.acos(max(-1.0, min(1.0, cosine)))
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
        radians = math.atan2(
            float(np.dot(np.cross(n1, n2), b2 / norm)), float(np.dot(n1, n2))
        )
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
    "normalize_numeric_value",
    "quantity_expression_semantic_signature",
    "quantity_expression_receipt_from_record",
]
