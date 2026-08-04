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
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

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
        "convert",
        "linear_fit_slope",
        "linear_fit_intercept",
    }
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
        "j/mol": (ENERGY, "hartree", energy_conversion("J/mol", "hartree", 1.0)),
        "j mol^-1": (ENERGY, "hartree", energy_conversion("J/mol", "hartree", 1.0)),
        "kj/mol": (ENERGY, "hartree", energy_conversion("kJ/mol", "hartree", 1.0)),
        "kj mol^-1": (ENERGY, "hartree", energy_conversion("kJ/mol", "hartree", 1.0)),
        "kcal/mol": (ENERGY, "hartree", energy_conversion("kcal/mol", "hartree", 1.0)),
        "kcal mol^-1": (ENERGY, "hartree", energy_conversion("kcal/mol", "hartree", 1.0)),
        "hartree/k": (ENTROPY, "hartree K^-1", 1.0),
        "hartree k^-1": (ENTROPY, "hartree K^-1", 1.0),
        "j/mol/k": (ENTROPY, "hartree K^-1", energy_conversion("J/mol", "hartree", 1.0)),
        "j mol^-1 k^-1": (ENTROPY, "hartree K^-1", energy_conversion("J/mol", "hartree", 1.0)),
        "kj/mol/k": (ENTROPY, "hartree K^-1", energy_conversion("kJ/mol", "hartree", 1.0)),
        "kj mol^-1 k^-1": (ENTROPY, "hartree K^-1", energy_conversion("kJ/mol", "hartree", 1.0)),
        "cal/mol/k": (ENTROPY, "hartree K^-1", energy_conversion("kcal/mol", "hartree", 0.001)),
        "cal mol^-1 k^-1": (ENTROPY, "hartree K^-1", energy_conversion("kcal/mol", "hartree", 0.001)),
        "kcal/mol/k": (ENTROPY, "hartree K^-1", energy_conversion("kcal/mol", "hartree", 1.0)),
        "kcal mol^-1 k^-1": (ENTROPY, "hartree K^-1", energy_conversion("kcal/mol", "hartree", 1.0)),
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


def normalize_numeric_value(value: Any, unit: str) -> tuple[Any, str, Dimension]:
    """Convert a finite numeric value to the canonical unit for its dimension."""

    dimension, canonical_unit, factor = _unit_spec(unit)
    array = np.asarray(value, dtype=float)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise QuantityExpressionError("numeric values must be finite and non-empty")
    normalized = array * factor
    if normalized.ndim == 0:
        payload: Any = float(normalized)
    elif normalized.ndim <= 2:
        payload = tuple(
            tuple(float(value) for value in row)
            if isinstance(row, np.ndarray)
            else float(row)
            for row in normalized
        )
    else:
        raise QuantityExpressionError("numeric arrays may have rank at most two")
    return payload, canonical_unit, dimension


def convert_normalized_value(
    value: Any,
    dimension: Dimension,
    target_unit: str,
) -> Any:
    """Convert a canonical numeric value to a compatible display unit."""

    target_dimension, _, factor = _unit_spec(target_unit)
    if target_dimension != dimension:
        raise QuantityExpressionError("target unit has an incompatible dimension")
    array = np.asarray(value, dtype=float) / factor
    if not np.all(np.isfinite(array)):
        raise QuantityExpressionError("unit conversion produced a non-finite value")
    if array.ndim == 0:
        return float(array)
    if array.ndim <= 2:
        return tuple(
            tuple(float(value) for value in row)
            if isinstance(row, np.ndarray)
            else float(row)
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

    def __post_init__(self) -> None:
        object.__setattr__(self, "input_ids", tuple(self.input_ids))
        object.__setattr__(self, "indices", tuple(self.indices))
        if not self.node_id or len(self.node_id) > 128:
            raise QuantityContractError("expression node_id is invalid")
        if self.operation not in _OPERATIONS:
            raise QuantityContractError(
                f"unsupported expression operation: {self.operation!r}"
            )
        if len(self.input_ids) > MAX_NODE_INPUTS:
            raise QuantityContractError("expression node has too many inputs")
        if any(index < 0 for index in self.indices):
            raise QuantityContractError("reference indices must be non-negative")
        if self.scale_factor is not None and not math.isfinite(self.scale_factor):
            raise QuantityContractError("scale_factor must be finite")


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
        object.__setattr__(self, "output_node_ids", tuple(self.output_node_ids))
        if self.schema_version != "chemsmart.quantity-expression-request.v1":
            raise QuantityContractError("unsupported expression request schema")
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
            raise QuantityContractError("at least one expression output is required")
        if len(self.output_node_ids) != len(set(self.output_node_ids)):
            raise QuantityContractError("expression output IDs must be unique")
        if not set(self.output_node_ids).issubset(node_ids):
            raise QuantityContractError("every output must identify an expression node")


@dataclass(frozen=True)
class QuantityExpressionOutputDependencyV1:
    """Receipt closure for one expression output."""

    output_id: str
    source_receipt_sha256s: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.output_id or len(self.output_id) > 128:
            raise QuantityContractError("expression output dependency ID is invalid")
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
            raise QuantityContractError("unsupported expression receipt schema")
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
                raise QuantityContractError(f"{label} must be a SHA-256 digest")
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
            raise QuantityContractError("quantity expression receipt digest mismatch")


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
    role_owners: dict[str, str] = {}
    for quantity in request.inputs:
        role_matches = _SEMANTIC_ROLE_REF.findall(quantity.evidence_ref)
        quantity_matches = _QUANTITY_REF.findall(quantity.evidence_ref)
        semantic_role = (
            role_matches[-1]
            if role_matches
            else quantity_matches[-1]
            if quantity_matches
            else quantity.quantity_id
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
            if node.operation in {"add", "multiply", "sum", "mean", "min", "max"}:
                inputs = sorted(inputs, key=canonical_quantity_sha256)
            elif node.operation == "distance" and len(inputs) == 2:
                inputs = sorted(inputs, key=canonical_quantity_sha256)
            elif node.operation == "angle" and len(inputs) == 3:
                outer = sorted((inputs[0], inputs[2]), key=canonical_quantity_sha256)
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
        raise QuantityExpressionError("expression input is non-finite or empty")
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
            raise QuantityExpressionError("ref requires one prior value reference")
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
            raise QuantityExpressionError("reference index is out of range") from exc
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
            if left_value.ndim and right_value.ndim and left_value.shape != right_value.shape:
                raise QuantityExpressionError(
                    "multiply accepts scalar broadcasting or identical shapes"
                )
            result = left_value * right_value
            dimension = _add_dimensions(left.dimension, right.dimension)
        else:
            if np.any(right_value == 0.0):
                raise QuantityExpressionError("division by zero is forbidden")
            if left_value.ndim and right_value.ndim and left_value.shape != right_value.shape:
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
            raise QuantityExpressionError("scale requires one input and scale_factor")
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
        dimension = tuple(
            exponent // 2 for exponent in source.dimension
        )
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
            raise QuantityExpressionError("power exponent must be a finite scalar")
        exponent = float(exponent_array)
        source = inputs[0]
        if source.dimension != DIMENSIONLESS:
            rounded = round(exponent)
            if not math.isclose(exponent, rounded, rel_tol=0.0, abs_tol=1.0e-12):
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
            raise QuantityExpressionError("negative power of zero is forbidden")
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
            raise QuantityExpressionError("log requires strictly positive values")
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
        if len(inputs) != 2 or any(item.dimension != LENGTH for item in inputs):
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
        if len(inputs) != 3 or any(item.dimension != LENGTH for item in inputs):
            raise QuantityExpressionError(
                "angle requires three length-coordinate vectors"
            )
        first, center, last = (_numeric(item) for item in inputs)
        if any(vector.ndim != 1 for vector in (first, center, last)):
            raise QuantityExpressionError("angle inputs must be coordinate vectors")
        if first.shape != center.shape or center.shape != last.shape:
            raise QuantityExpressionError("angle coordinate shapes must match")
        left = first - center
        right = last - center
        denominator = float(np.linalg.norm(left) * np.linalg.norm(right))
        if denominator == 0.0:
            raise QuantityExpressionError("angle is undefined for zero-length vectors")
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

    if operation == "convert":
        if len(inputs) != 1 or not node.target_unit:
            raise QuantityExpressionError("convert requires one input and target_unit")
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

    raise QuantityExpressionError(f"operation is not implemented: {operation!r}")


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
    derived: list[QuantityValueV1] = []
    evidence_ref = f"expression:{request.expression_id}#{request_sha256}"
    for node in request.nodes:
        if node.node_id in values:
            raise QuantityExpressionError(
                f"expression node ID collides with an existing value: {node.node_id}"
            )
        value = _node_value(node=node, values=values, evidence_ref=evidence_ref)
        values[node.node_id] = value
        derived.append(value)
        if node.operation == "literal":
            dependencies[node.node_id] = frozenset()
        elif node.operation == "ref":
            reference = node.reference or node.input_ids[0]
            dependencies[node.node_id] = dependencies[reference]
        else:
            dependencies[node.node_id] = frozenset().union(
                *(dependencies[input_id] for input_id in node.input_ids)
            )
    outputs = tuple(values[node_id] for node_id in request.output_node_ids)
    output_dependencies = tuple(
        QuantityExpressionOutputDependencyV1(
            output_id=node_id,
            source_receipt_sha256s=tuple(sorted(dependencies[node_id])),
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
        )
        for item in values.get("output_dependencies") or ()
    )
    return QuantityExpressionReceiptV1(
        **values, receipt_sha256=receipt_sha256
    )


__all__ = [
    "MAX_EXPRESSION_NODES",
    "MAX_NODE_INPUTS",
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
