"""Deterministic evaluation of planned scientific-analysis predicates.

The model may choose the checks while it declares a scientific toolchain, but
it cannot restate or weaken them at execution time.  This module evaluates the
rules already sealed into one ``scientific_validation`` node against exact
typed quantities produced by upstream ChemSmart analysis receipts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.scientific_toolchain import AnalysisNodeIntentV1
from chemsmart.analysis.quantity_expressions import convert_normalized_value
from chemsmart.analysis.result_quantities import (
    DIMENSIONLESS,
    QuantityValueV1,
    make_quantity_value,
    quantity_value_from_record,
)


@dataclass(frozen=True)
class ScientificValidationInputBindingV1:
    """Exact typed quantity bound to one declared validation input."""

    input_id: str
    source_receipt_sha256: str
    quantity_id: str
    quantity_value_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.input_id, "validation input_id")
        require_sha256(
            self.source_receipt_sha256,
            "validation source_receipt_sha256",
        )
        require_identifier(self.quantity_id, "validation quantity_id")
        require_sha256(
            self.quantity_value_sha256,
            "validation quantity_value_sha256",
        )


@dataclass(frozen=True)
class ScientificValidationRuleResultV1:
    """Host-observed outcome of one rule from the planned analysis node."""

    rule_id: str
    predicate: str
    input_ids: tuple[str, ...]
    passed: bool
    observed_value: float | int
    threshold: float | None
    expected_count: int | None
    unit: str

    def __post_init__(self) -> None:
        require_identifier(self.rule_id, "validation rule_id")
        object.__setattr__(self, "input_ids", tuple(self.input_ids))
        if not self.input_ids or self.input_ids != tuple(
            sorted(set(self.input_ids))
        ):
            raise ContractError(
                "validation rule input IDs must be sorted and unique"
            )
        if isinstance(self.observed_value, bool) or not isinstance(
            self.observed_value, (int, float)
        ):
            raise ContractError("validation observed value must be numerical")
        if not math.isfinite(float(self.observed_value)):
            raise ContractError("validation observed value must be finite")


@dataclass(frozen=True)
class ScientificValidationReceiptV1:
    """Typed receipt proving that a declared validation node was evaluated."""

    schema_version: str
    workflow_id: str
    plan_sha256: str
    node_id: str
    input_bindings: tuple[ScientificValidationInputBindingV1, ...]
    source_receipt_sha256s: tuple[str, ...]
    rule_results: tuple[ScientificValidationRuleResultV1, ...]
    outputs: tuple[QuantityValueV1, ...]
    all_rules_passed: bool
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-validation-receipt.v1":
            raise ContractError("unsupported scientific validation receipt")
        require_identifier(self.workflow_id, "validation workflow_id")
        require_sha256(self.plan_sha256, "validation plan_sha256")
        require_identifier(self.node_id, "validation node_id")
        object.__setattr__(self, "input_bindings", tuple(self.input_bindings))
        object.__setattr__(
            self, "source_receipt_sha256s", tuple(self.source_receipt_sha256s)
        )
        object.__setattr__(self, "rule_results", tuple(self.rule_results))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        input_ids = tuple(item.input_id for item in self.input_bindings)
        if not input_ids or input_ids != tuple(sorted(set(input_ids))):
            raise ContractError(
                "validation input bindings must be sorted and unique"
            )
        if self.source_receipt_sha256s != tuple(
            sorted(set(self.source_receipt_sha256s))
        ):
            raise ContractError(
                "validation source receipts must be sorted and unique"
            )
        for digest in self.source_receipt_sha256s:
            require_sha256(digest, "validation source receipt")
        rule_ids = tuple(item.rule_id for item in self.rule_results)
        if not rule_ids or rule_ids != tuple(sorted(set(rule_ids))):
            raise ContractError(
                "validation rule results must be sorted and unique"
            )
        output_ids = tuple(item.quantity_id for item in self.outputs)
        if len(output_ids) != 1 or output_ids != tuple(
            sorted(set(output_ids))
        ):
            raise ContractError(
                "scientific validation requires one typed verdict output"
            )
        if (
            self.outputs[0].dimension != DIMENSIONLESS
            or self.outputs[0].data_kind != "integer"
            or self.outputs[0].value not in {0, 1}
        ):
            raise ContractError(
                "scientific validation verdict must be dimensionless 0 or 1"
            )
        if self.all_rules_passed is not all(
            item.passed for item in self.rule_results
        ):
            raise ContractError("validation aggregate verdict is inconsistent")
        if self.status != "evaluated":
            raise ContractError(
                "scientific validation status must be evaluated"
            )
        if self.receipt_sha256 != canonical_sha256(
            scientific_validation_receipt_body(self)
        ):
            raise ContractError(
                "scientific validation receipt digest mismatch"
            )


def scientific_validation_receipt_body(
    receipt: ScientificValidationReceiptV1,
) -> dict[str, Any]:
    return {
        "schema_version": receipt.schema_version,
        "workflow_id": receipt.workflow_id,
        "plan_sha256": receipt.plan_sha256,
        "node_id": receipt.node_id,
        "input_bindings": receipt.input_bindings,
        "source_receipt_sha256s": receipt.source_receipt_sha256s,
        "rule_results": receipt.rule_results,
        "outputs": receipt.outputs,
        "all_rules_passed": receipt.all_rules_passed,
        "status": receipt.status,
    }


def scientific_validation_receipt_from_record(
    record: Mapping[str, Any], *, receipt_sha256: str
) -> ScientificValidationReceiptV1:
    """Rehydrate a validation receipt from one canonical Runtime event."""

    values = dict(record)
    values["input_bindings"] = tuple(
        ScientificValidationInputBindingV1(**dict(item))
        for item in values.get("input_bindings") or ()
    )
    values["source_receipt_sha256s"] = tuple(
        values.get("source_receipt_sha256s") or ()
    )
    values["rule_results"] = tuple(
        ScientificValidationRuleResultV1(
            **{
                **dict(item),
                "input_ids": tuple(dict(item).get("input_ids") or ()),
            }
        )
        for item in values.get("rule_results") or ()
    )
    values["outputs"] = tuple(
        quantity_value_from_record(item)
        for item in values.get("outputs") or ()
    )
    return ScientificValidationReceiptV1(
        **values, receipt_sha256=receipt_sha256
    )


def evaluate_planned_scientific_validation(
    *,
    workflow_id: str,
    plan_sha256: str,
    node: AnalysisNodeIntentV1,
    inputs: Mapping[str, tuple[str, QuantityValueV1]],
) -> ScientificValidationReceiptV1:
    """Evaluate every rule sealed into one planned validation node."""

    if node.analysis_kind != "scientific_validation":
        raise ContractError("requested node is not scientific_validation")
    if node.support_state != "planned":
        raise ContractError("blocked validation nodes cannot be evaluated")
    if not node.validation_rules:
        raise ContractError("scientific validation node declares no rules")
    expected_inputs = tuple(sorted(item.input_id for item in node.inputs))
    if tuple(sorted(inputs)) != expected_inputs:
        raise ContractError(
            "scientific validation inputs differ from the planned node"
        )

    bindings = tuple(
        ScientificValidationInputBindingV1(
            input_id=input_id,
            source_receipt_sha256=inputs[input_id][0],
            quantity_id=inputs[input_id][1].quantity_id,
            quantity_value_sha256=inputs[input_id][1].value_sha256,
        )
        for input_id in expected_inputs
    )
    rule_results = tuple(
        _evaluate_rule(rule=rule, inputs=inputs)
        for rule in node.validation_rules
    )
    all_rules_passed = all(item.passed for item in rule_results)
    output = node.outputs[0]
    verdict = make_quantity_value(
        quantity_id=output.output_id,
        source_value=int(all_rules_passed),
        source_unit="1",
        value=int(all_rules_passed),
        unit="1",
        dimension=DIMENSIONLESS,
        evidence_ref=(
            f"scientific-validation:{plan_sha256}:{node.node_id};"
            + ";".join(item.rule_id for item in rule_results)
        ),
        data_kind="integer",
    )
    body = {
        "schema_version": "chemsmart.scientific-validation-receipt.v1",
        "workflow_id": workflow_id,
        "plan_sha256": plan_sha256,
        "node_id": node.node_id,
        "input_bindings": bindings,
        "source_receipt_sha256s": tuple(
            sorted({receipt for receipt, _ in inputs.values()})
        ),
        "rule_results": rule_results,
        "outputs": (verdict,),
        "all_rules_passed": all_rules_passed,
        # Evaluated is operationally green even when the scientific predicate
        # is false: a failed minimum test is still a completed scientific
        # determination and must remain reportable as such.
        "status": "evaluated",
    }
    return ScientificValidationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _evaluate_rule(
    *, rule: Any, inputs: Mapping[str, tuple[str, QuantityValueV1]]
):
    quantities = tuple(inputs[input_id][1] for input_id in rule.input_ids)
    predicate = rule.predicate

    if predicate == "count_equals":
        if len(quantities) != 1:
            raise ContractError("count_equals requires one scalar count")
        values = _numeric_values(quantities[0], target_unit="1")
        if len(values) != 1 or not float(values[0]).is_integer():
            raise ContractError("count_equals input must be an integer scalar")
        observed: float | int = int(values[0])
        passed = observed == rule.expected_count
    elif predicate == "all_finite":
        observed = 1
        for quantity in quantities:
            _numeric_values(quantity, target_unit="")
        passed = True
    elif predicate == "all_equal":
        if len({tuple(quantity.dimension) for quantity in quantities}) != 1:
            raise ContractError(
                "all_equal inputs have incompatible dimensions"
            )
        flattened = tuple(
            value
            for quantity in quantities
            for value in _numeric_values(quantity, target_unit="")
        )
        reference = flattened[0]
        observed = max(abs(value - reference) for value in flattened)
        passed = observed == 0.0
    elif predicate == "minimum_greater_equal":
        flattened = tuple(
            value
            for quantity in quantities
            for value in _numeric_values(quantity, target_unit=rule.unit)
        )
        observed = min(flattened)
        passed = observed >= float(rule.threshold)
    elif predicate == "maximum_absolute_less_equal":
        flattened = tuple(
            value
            for quantity in quantities
            for value in _numeric_values(quantity, target_unit=rule.unit)
        )
        observed = max(abs(value) for value in flattened)
        passed = observed <= float(rule.threshold)
    elif predicate == "symmetric_within":
        if len(quantities) != 1:
            raise ContractError("symmetric_within requires one matrix")
        converted = convert_normalized_value(
            quantities[0].value,
            quantities[0].dimension,
            rule.unit,
        )
        matrix = np.asarray(converted, dtype=float)
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
            raise ContractError(
                "symmetric_within input must be a square matrix"
            )
        observed = float(np.max(np.abs(matrix - matrix.T)))
        passed = observed <= float(rule.threshold)
    else:  # pragma: no cover - the planning contract owns this vocabulary
        raise ContractError("unsupported scientific validation predicate")

    return ScientificValidationRuleResultV1(
        rule_id=rule.rule_id,
        predicate=predicate,
        input_ids=rule.input_ids,
        passed=bool(passed),
        observed_value=observed,
        threshold=rule.threshold,
        expected_count=rule.expected_count,
        unit=rule.unit,
    )


def _numeric_values(
    quantity: QuantityValueV1, *, target_unit: str
) -> tuple[float, ...]:
    if quantity.data_kind in {"text", "text_vector"}:
        raise ContractError("scientific validation requires numerical inputs")
    try:
        value = (
            convert_normalized_value(
                quantity.value, quantity.dimension, target_unit
            )
            if target_unit
            else quantity.value
        )
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ContractError(
            "scientific validation input is not numerical"
        ) from exc
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ContractError(
            "scientific validation input must be finite and non-empty"
        )
    return tuple(float(value) for value in array.reshape(-1))


__all__ = [
    "ScientificValidationInputBindingV1",
    "ScientificValidationReceiptV1",
    "ScientificValidationRuleResultV1",
    "evaluate_planned_scientific_validation",
    "scientific_validation_receipt_body",
    "scientific_validation_receipt_from_record",
]
