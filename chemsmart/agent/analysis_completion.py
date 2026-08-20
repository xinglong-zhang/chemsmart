"""Task-owned completion policy for structured scientific analysis.

The policy is supplied by the user or calling application, bound to the exact
task digest by the host, and evaluated only against deterministic receipts.
It describes generic analysis stages and requested quantity classes without
embedding a task answer or paper-specific formula.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)

ANALYSIS_COMPLETION_STAGES = frozenset(
    {
        "quantity_extraction",
        "thermochemistry",
        "quantity_expression",
        "analysis_claims",
        "scientific_decision",
    }
)


class AnalysisIncompleteError(ContractError):
    """Expected absence of a task-owned analysis requirement."""


@dataclass(frozen=True)
class AnalysisConditionV1:
    """Task-owned numerical condition with an explicit evidence origin."""

    condition_id: str
    value: float
    unit: str
    origin: str
    evidence_ref: str

    def __post_init__(self) -> None:
        require_identifier(self.condition_id, "analysis_condition_id")
        if (
            isinstance(self.value, bool)
            or not isinstance(self.value, (int, float))
            or not math.isfinite(float(self.value))
        ):
            raise ContractError("analysis condition value must be finite")
        if not self.unit.strip() or not self.evidence_ref.strip():
            raise ContractError(
                "analysis condition requires a unit and evidence reference"
            )
        if self.origin not in {
            "paper",
            "supplementary_information",
            "benchmark",
            "user",
            "workspace",
            "computed",
        }:
            raise ContractError("unsupported analysis condition origin")


@dataclass(frozen=True)
class AnalysisExpressionSourceRequirementV1:
    """Upstream deterministic evidence required by one expression.

    Empty ``artifact_sha256s`` means every target artifact named by the
    enclosing policy.  A non-empty tuple permits an expression to depend on a
    scientifically meaningful subset, such as the reactant side of a larger
    reaction bundle, without encoding a predetermined answer.
    """

    stage: str
    artifact_sha256s: tuple[str, ...] = ()
    output_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.stage not in {"quantity_extraction", "thermochemistry"}:
            raise ContractError(
                "expression source stage must be quantity_extraction or "
                "thermochemistry"
            )
        if self.artifact_sha256s != tuple(sorted(set(self.artifact_sha256s))):
            raise ContractError(
                "expression source artifact hashes must be sorted and unique"
            )
        for digest in self.artifact_sha256s:
            require_sha256(digest, "expression_source_artifact_sha256")
        if self.output_ids != tuple(sorted(set(self.output_ids))):
            raise ContractError(
                "expression source output IDs must be sorted and unique"
            )
        for output_id in self.output_ids:
            require_identifier(output_id, "expression_source_output_id")


@dataclass(frozen=True)
class AnalysisExpressionRequirementV1:
    """Named dimensional-AST result required by a task policy."""

    expression_id: str
    required_output_ids: tuple[str, ...]
    required_sources: tuple[AnalysisExpressionSourceRequirementV1, ...] = ()
    semantic_signature_sha256: str = ""

    def __post_init__(self) -> None:
        require_identifier(self.expression_id, "expression_id")
        if not self.required_output_ids:
            raise ContractError("analysis expression requires output IDs")
        if self.required_output_ids != tuple(
            sorted(set(self.required_output_ids))
        ):
            raise ContractError(
                "analysis expression output IDs must be sorted and unique"
            )
        for output_id in self.required_output_ids:
            require_identifier(output_id, "required_output_id")
        for source in self.required_sources:
            if not set(source.output_ids).issubset(self.required_output_ids):
                raise ContractError(
                    "expression source output IDs must be required outputs"
                )
        if self.semantic_signature_sha256:
            require_sha256(
                self.semantic_signature_sha256,
                "expression_semantic_signature_sha256",
            )
        source_keys = tuple(
            (source.stage, source.artifact_sha256s, source.output_ids)
            for source in self.required_sources
        )
        if source_keys != tuple(sorted(set(source_keys))):
            raise ContractError(
                "analysis expression source requirements must be sorted and unique"
            )


@dataclass(frozen=True)
class AnalysisClaimRequirementV1:
    """Task-owned identity and display unit for one reported value."""

    claim_id: str
    source_kind: str
    quantity_id: str
    display_unit: str
    source_artifact_sha256s: tuple[str, ...] = ()
    source_expression_id: str = ""
    source_selector: str = ""

    def __post_init__(self) -> None:
        require_identifier(self.claim_id, "claim_id")
        if self.source_kind not in {
            "quantity_extraction",
            "thermochemistry",
            "quantity_expression",
        }:
            raise ContractError("unsupported analysis claim source kind")
        require_identifier(self.quantity_id, "claim_quantity_id")
        if not self.display_unit:
            raise ContractError(
                "analysis claim display unit must not be empty"
            )
        if self.source_artifact_sha256s != tuple(
            sorted(set(self.source_artifact_sha256s))
        ):
            raise ContractError(
                "analysis claim artifact hashes must be sorted and unique"
            )
        for digest in self.source_artifact_sha256s:
            require_sha256(digest, "analysis_claim_artifact_sha256")
        if len(self.source_artifact_sha256s) > 1:
            raise ContractError(
                "one direct claim must bind one artifact; use an expression for aggregates"
            )
        if self.source_expression_id:
            require_identifier(
                self.source_expression_id, "analysis_claim_expression_id"
            )
        if self.source_kind == "quantity_expression":
            if self.source_artifact_sha256s or self.source_selector:
                raise ContractError(
                    "expression claims cannot declare artifact or selector bindings"
                )
        elif self.source_expression_id:
            raise ContractError(
                "non-expression claims cannot declare an expression binding"
            )
        if self.source_selector and self.source_kind != "quantity_extraction":
            raise ContractError(
                "only extraction claims can declare a source selector"
            )


@dataclass(frozen=True)
class AnalysisCompletionPolicyV1:
    """Host-bound requirements for declaring an analysis task complete."""

    schema_version: str
    policy_id: str
    task_spec_sha256: str
    target_artifact_sha256s: tuple[str, ...]
    required_stages: tuple[str, ...]
    required_extraction_selectors: tuple[str, ...]
    required_thermochemistry_quantity_ids: tuple[str, ...]
    required_expressions: tuple[AnalysisExpressionRequirementV1, ...]
    required_claims: tuple[AnalysisClaimRequirementV1, ...]
    required_conditions: tuple[AnalysisConditionV1, ...]
    temperature_k: float | None
    pressure_atm: float | None
    minimum_expression_receipts: int
    require_decision_evidence_binding: bool
    policy_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.analysis-completion-policy.v1":
            raise ContractError(
                "unsupported analysis completion policy schema"
            )
        require_identifier(self.policy_id, "policy_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if self.target_artifact_sha256s != tuple(
            sorted(set(self.target_artifact_sha256s))
        ):
            raise ContractError(
                "target artifact hashes must be sorted and unique"
            )
        if not self.target_artifact_sha256s and {
            "quantity_extraction",
            "thermochemistry",
        }.intersection(self.required_stages):
            raise ContractError(
                "artifact analysis stages require target artifact hashes"
            )
        for digest in self.target_artifact_sha256s:
            require_sha256(digest, "target_artifact_sha256")
        if not self.required_stages:
            raise ContractError(
                "analysis completion requires at least one stage"
            )
        if self.required_stages != tuple(sorted(set(self.required_stages))):
            raise ContractError(
                "required analysis stages must be sorted and unique"
            )
        unknown = set(self.required_stages).difference(
            ANALYSIS_COMPLETION_STAGES
        )
        if unknown:
            raise ContractError(
                "unsupported analysis completion stage: "
                + ", ".join(sorted(unknown))
            )
        if self.required_thermochemistry_quantity_ids != tuple(
            sorted(set(self.required_thermochemistry_quantity_ids))
        ):
            raise ContractError(
                "required thermochemistry quantity IDs must be sorted and unique"
            )
        if self.required_extraction_selectors != tuple(
            sorted(set(self.required_extraction_selectors))
        ):
            raise ContractError(
                "required extraction selectors must be sorted and unique"
            )
        if self.required_extraction_selectors and (
            "quantity_extraction" not in self.required_stages
        ):
            raise ContractError(
                "extraction selectors require the quantity_extraction stage"
            )
        expression_ids = tuple(
            requirement.expression_id
            for requirement in self.required_expressions
        )
        if expression_ids != tuple(sorted(set(expression_ids))):
            raise ContractError(
                "analysis expression requirements must be sorted and unique"
            )
        if self.required_expressions and (
            "quantity_expression" not in self.required_stages
        ):
            raise ContractError(
                "expression requirements require the quantity_expression stage"
            )
        claim_ids = tuple(
            requirement.claim_id for requirement in self.required_claims
        )
        if claim_ids != tuple(sorted(set(claim_ids))):
            raise ContractError(
                "analysis claim requirements must be sorted and unique"
            )
        if (
            self.required_claims
            and "analysis_claims" not in self.required_stages
        ):
            raise ContractError(
                "claim requirements require the analysis_claims stage"
            )
        if (
            "analysis_claims" in self.required_stages
            and not self.required_claims
        ):
            raise ContractError(
                "analysis_claims requires named claim requirements"
            )
        for requirement in self.required_claims:
            source_stage = {
                "quantity_extraction": "quantity_extraction",
                "thermochemistry": "thermochemistry",
                "quantity_expression": "quantity_expression",
            }[requirement.source_kind]
            if source_stage not in self.required_stages:
                raise ContractError(
                    "analysis claim source stage must be enabled by the policy"
                )
            unknown_artifacts = set(
                requirement.source_artifact_sha256s
            ).difference(self.target_artifact_sha256s)
            if unknown_artifacts:
                raise ContractError(
                    "analysis claim source artifacts must be policy targets"
                )
            if requirement.source_expression_id and (
                requirement.source_expression_id not in expression_ids
            ):
                raise ContractError(
                    "analysis claim source expression must be policy-required"
                )
        condition_ids = tuple(
            condition.condition_id for condition in self.required_conditions
        )
        if condition_ids != tuple(sorted(set(condition_ids))):
            raise ContractError(
                "analysis condition IDs must be sorted and unique"
            )
        for quantity_id in self.required_thermochemistry_quantity_ids:
            require_identifier(
                quantity_id, "required_thermochemistry_quantity_id"
            )
        if self.required_thermochemistry_quantity_ids and (
            "thermochemistry" not in self.required_stages
        ):
            raise ContractError(
                "thermochemistry quantities require the thermochemistry stage"
            )
        for requirement in self.required_expressions:
            for source in requirement.required_sources:
                if source.stage not in self.required_stages:
                    raise ContractError(
                        "expression source stage must be enabled by the policy"
                    )
                unknown_artifacts = set(source.artifact_sha256s).difference(
                    self.target_artifact_sha256s
                )
                if unknown_artifacts:
                    raise ContractError(
                        "expression source artifacts must be policy targets"
                    )
        for value, label in (
            (self.temperature_k, "temperature_k"),
            (self.pressure_atm, "pressure_atm"),
        ):
            if value is not None and (
                not isinstance(value, (int, float))
                or isinstance(value, bool)
                or not math.isfinite(float(value))
                or float(value) <= 0.0
            ):
                raise ContractError(f"{label} must be finite and positive")
        if (
            self.temperature_k is not None or self.pressure_atm is not None
        ) and ("thermochemistry" not in self.required_stages):
            raise ContractError(
                "thermochemistry conditions require the thermochemistry stage"
            )
        if (
            not isinstance(self.minimum_expression_receipts, int)
            or isinstance(self.minimum_expression_receipts, bool)
            or self.minimum_expression_receipts < 0
        ):
            raise ContractError(
                "minimum_expression_receipts must be a non-negative integer"
            )
        if self.minimum_expression_receipts and (
            "quantity_expression" not in self.required_stages
        ):
            raise ContractError(
                "expression receipt minimum requires the quantity_expression stage"
            )
        if (
            "quantity_expression" in self.required_stages
            and self.minimum_expression_receipts < 1
        ):
            raise ContractError(
                "quantity_expression requires at least one expression receipt"
            )
        if self.minimum_expression_receipts != len(self.required_expressions):
            raise ContractError(
                "every required expression receipt must have a named contract"
            )
        if self.require_decision_evidence_binding and (
            "scientific_decision" not in self.required_stages
        ):
            raise ContractError(
                "decision evidence binding requires the scientific_decision stage"
            )
        body = _policy_body(self)
        if self.policy_sha256 != canonical_sha256(body):
            raise ContractError("analysis completion policy digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return {**_policy_body(self), "policy_sha256": self.policy_sha256}


@dataclass(frozen=True)
class AnalysisCompletionReceiptV1:
    """Aggregate gate over the deterministic receipts required by a policy."""

    schema_version: str
    policy_sha256: str
    task_spec_sha256: str
    source_receipt_sha256s: tuple[str, ...]
    status: str
    findings: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.analysis-completion-receipt.v1":
            raise ContractError(
                "unsupported analysis completion receipt schema"
            )
        require_sha256(self.policy_sha256, "policy_sha256")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if not self.source_receipt_sha256s:
            raise ContractError(
                "analysis completion receipt requires source receipts"
            )
        if self.source_receipt_sha256s != tuple(
            sorted(set(self.source_receipt_sha256s))
        ):
            raise ContractError(
                "analysis completion source receipts must be sorted and unique"
            )
        for digest in self.source_receipt_sha256s:
            require_sha256(digest, "source_receipt_sha256")
        if self.status != "passed" or self.findings:
            raise ContractError(
                "analysis completion receipt must represent a green evaluation"
            )
        body = _completion_receipt_body(self)
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("analysis completion receipt digest mismatch")


def _completion_receipt_body(
    receipt: AnalysisCompletionReceiptV1,
) -> dict[str, Any]:
    return {
        "schema_version": receipt.schema_version,
        "policy_sha256": receipt.policy_sha256,
        "task_spec_sha256": receipt.task_spec_sha256,
        "source_receipt_sha256s": receipt.source_receipt_sha256s,
        "status": receipt.status,
        "findings": receipt.findings,
    }


def build_analysis_completion_receipt(
    *,
    policy: AnalysisCompletionPolicyV1,
    source_receipt_sha256s: tuple[str, ...],
) -> AnalysisCompletionReceiptV1:
    body = {
        "schema_version": "chemsmart.analysis-completion-receipt.v1",
        "policy_sha256": policy.policy_sha256,
        "task_spec_sha256": policy.task_spec_sha256,
        "source_receipt_sha256s": tuple(sorted(set(source_receipt_sha256s))),
        "status": "passed",
        "findings": (),
    }
    return AnalysisCompletionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _policy_body(policy: AnalysisCompletionPolicyV1) -> dict[str, Any]:
    return {
        "schema_version": policy.schema_version,
        "policy_id": policy.policy_id,
        "task_spec_sha256": policy.task_spec_sha256,
        "target_artifact_sha256s": policy.target_artifact_sha256s,
        "required_stages": policy.required_stages,
        "required_extraction_selectors": policy.required_extraction_selectors,
        "required_thermochemistry_quantity_ids": (
            policy.required_thermochemistry_quantity_ids
        ),
        "required_expressions": policy.required_expressions,
        "required_claims": policy.required_claims,
        "required_conditions": policy.required_conditions,
        "temperature_k": policy.temperature_k,
        "pressure_atm": policy.pressure_atm,
        "minimum_expression_receipts": policy.minimum_expression_receipts,
        "require_decision_evidence_binding": (
            policy.require_decision_evidence_binding
        ),
    }


def load_analysis_completion_policy(
    path: str | Path,
    *,
    task_spec_sha256: str,
) -> AnalysisCompletionPolicyV1:
    """Load a data-only policy and bind it to the exact active task."""

    source = Path(path)
    try:
        values = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(
            "analysis completion policy is not valid JSON"
        ) from exc
    if not isinstance(values, Mapping):
        raise ContractError("analysis completion policy must be a JSON object")
    allowed = {
        "schema_version",
        "policy_id",
        "target_artifact_sha256s",
        "required_stages",
        "required_extraction_selectors",
        "required_thermochemistry_quantity_ids",
        "required_expressions",
        "required_claims",
        "required_conditions",
        "temperature_k",
        "pressure_atm",
        "minimum_expression_receipts",
        "require_decision_evidence_binding",
    }
    unknown = sorted(set(values).difference(allowed))
    if unknown:
        raise ContractError(
            "analysis completion policy fields not allowed: "
            + ", ".join(unknown)
        )
    if values.get("schema_version") != (
        "chemsmart.analysis-completion-policy-spec.v1"
    ):
        raise ContractError(
            "unsupported analysis completion policy input schema"
        )
    raw_stages = values.get("required_stages", ())
    raw_quantity_ids = values.get("required_thermochemistry_quantity_ids", ())
    raw_selectors = values.get("required_extraction_selectors", ())
    raw_expressions = values.get("required_expressions", ())
    raw_claims = values.get("required_claims", ())
    raw_conditions = values.get("required_conditions", ())
    if not isinstance(raw_stages, list) or not all(
        isinstance(item, str) for item in raw_stages
    ):
        raise ContractError("required_stages must be an array of strings")
    if not isinstance(raw_quantity_ids, list) or not all(
        isinstance(item, str) for item in raw_quantity_ids
    ):
        raise ContractError(
            "required_thermochemistry_quantity_ids must be an array of strings"
        )
    if not isinstance(raw_selectors, list) or not all(
        isinstance(item, str) for item in raw_selectors
    ):
        raise ContractError(
            "required_extraction_selectors must be an array of strings"
        )
    if not isinstance(raw_expressions, list):
        raise ContractError("required_expressions must be an array")
    if not isinstance(raw_claims, list):
        raise ContractError("required_claims must be an array")
    if not isinstance(raw_conditions, list):
        raise ContractError("required_conditions must be an array")
    expression_requirements = []
    for item in raw_expressions:
        if (
            not isinstance(item, Mapping)
            or not {
                "expression_id",
                "required_output_ids",
            }.issubset(item)
            or set(item).difference(
                {
                    "expression_id",
                    "required_output_ids",
                    "required_sources",
                    "semantic_signature_sha256",
                }
            )
        ):
            raise ContractError(
                "each required expression needs expression_id and "
                "required_output_ids, with optional required_sources"
            )
        output_ids = item["required_output_ids"]
        if not isinstance(output_ids, list) or not all(
            isinstance(value, str) for value in output_ids
        ):
            raise ContractError(
                "required expression output IDs must be an array of strings"
            )
        raw_sources = item.get("required_sources", [])
        if not isinstance(raw_sources, list):
            raise ContractError("required_sources must be an array")
        source_requirements = []
        for source in raw_sources:
            if (
                not isinstance(source, Mapping)
                or not {"stage", "artifact_sha256s"}.issubset(source)
                or set(source).difference(
                    {"stage", "artifact_sha256s", "output_ids"}
                )
            ):
                raise ContractError(
                    "each expression source needs stage and artifact_sha256s, "
                    "with optional output_ids"
                )
            source_artifacts = source["artifact_sha256s"]
            if not isinstance(source_artifacts, list) or not all(
                isinstance(value, str) for value in source_artifacts
            ):
                raise ContractError(
                    "expression source artifact_sha256s must be an array of strings"
                )
            source_outputs = source.get("output_ids", [])
            if not isinstance(source_outputs, list) or not all(
                isinstance(value, str) for value in source_outputs
            ):
                raise ContractError(
                    "expression source output_ids must be an array of strings"
                )
            source_requirements.append(
                AnalysisExpressionSourceRequirementV1(
                    stage=str(source["stage"]).strip().lower(),
                    artifact_sha256s=tuple(
                        sorted(
                            {
                                str(value).strip().lower()
                                for value in source_artifacts
                            }
                        )
                    ),
                    output_ids=tuple(
                        sorted(
                            {
                                str(value).strip().lower()
                                for value in source_outputs
                            }
                        )
                    ),
                )
            )
        expression_requirements.append(
            AnalysisExpressionRequirementV1(
                expression_id=str(item["expression_id"]).strip().lower(),
                required_output_ids=tuple(
                    sorted(
                        {str(value).strip().lower() for value in output_ids}
                    )
                ),
                required_sources=tuple(
                    sorted(
                        source_requirements,
                        key=lambda source: (
                            source.stage,
                            source.artifact_sha256s,
                            source.output_ids,
                        ),
                    )
                ),
                semantic_signature_sha256=str(
                    item.get("semantic_signature_sha256") or ""
                )
                .strip()
                .lower(),
            )
        )
    claim_requirements = []
    for item in raw_claims:
        required_claim_fields = {
            "claim_id",
            "source_kind",
            "quantity_id",
            "display_unit",
        }
        optional_claim_fields = {
            "source_artifact_sha256s",
            "source_expression_id",
            "source_selector",
        }
        if (
            not isinstance(item, Mapping)
            or not required_claim_fields.issubset(item)
            or set(item).difference(
                required_claim_fields | optional_claim_fields
            )
        ):
            raise ContractError(
                "each required claim needs claim_id, source_kind, quantity_id, "
                "and display_unit, with optional deterministic source bindings"
            )
        raw_claim_artifacts = item.get("source_artifact_sha256s", [])
        if not isinstance(raw_claim_artifacts, list) or not all(
            isinstance(value, str) for value in raw_claim_artifacts
        ):
            raise ContractError(
                "claim source_artifact_sha256s must be an array of strings"
            )
        claim_requirements.append(
            AnalysisClaimRequirementV1(
                claim_id=str(item["claim_id"]).strip().lower(),
                source_kind=str(item["source_kind"]).strip().lower(),
                quantity_id=str(item["quantity_id"]).strip().lower(),
                display_unit=str(item["display_unit"]).strip(),
                source_artifact_sha256s=tuple(
                    sorted(
                        {
                            str(value).strip().lower()
                            for value in raw_claim_artifacts
                        }
                    )
                ),
                source_expression_id=str(
                    item.get("source_expression_id") or ""
                )
                .strip()
                .lower(),
                source_selector=str(item.get("source_selector") or "")
                .strip()
                .lower(),
            )
        )
    condition_requirements = []
    for item in raw_conditions:
        if not isinstance(item, Mapping) or set(item) != {
            "condition_id",
            "value",
            "unit",
            "origin",
            "evidence_ref",
        }:
            raise ContractError(
                "each analysis condition needs ID, value, unit, origin, and evidence"
            )
        condition_requirements.append(
            AnalysisConditionV1(
                condition_id=str(item["condition_id"]).strip().lower(),
                value=float(item["value"]),
                unit=str(item["unit"]).strip(),
                origin=str(item["origin"]).strip().lower(),
                evidence_ref=str(item["evidence_ref"]).strip(),
            )
        )
    raw_artifacts = values.get("target_artifact_sha256s", ())
    if not isinstance(raw_artifacts, list) or not all(
        isinstance(item, str) for item in raw_artifacts
    ):
        raise ContractError(
            "target_artifact_sha256s must be an array of strings"
        )
    raw_minimum = values.get("minimum_expression_receipts", 0)
    if not isinstance(raw_minimum, int) or isinstance(raw_minimum, bool):
        raise ContractError("minimum_expression_receipts must be an integer")
    raw_binding = values.get("require_decision_evidence_binding", False)
    if not isinstance(raw_binding, bool):
        raise ContractError(
            "require_decision_evidence_binding must be a boolean"
        )
    body = {
        "schema_version": "chemsmart.analysis-completion-policy.v1",
        "policy_id": str(values.get("policy_id") or ""),
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "target_artifact_sha256s": tuple(
            sorted({str(item).strip().lower() for item in raw_artifacts})
        ),
        "required_stages": tuple(
            sorted({str(item).strip().lower() for item in raw_stages})
        ),
        "required_thermochemistry_quantity_ids": tuple(
            sorted({str(item).strip().lower() for item in raw_quantity_ids})
        ),
        "required_extraction_selectors": tuple(
            sorted({str(item).strip().lower() for item in raw_selectors})
        ),
        "required_expressions": tuple(
            sorted(
                expression_requirements, key=lambda item: item.expression_id
            )
        ),
        "required_claims": tuple(
            sorted(claim_requirements, key=lambda item: item.claim_id)
        ),
        "required_conditions": tuple(
            sorted(condition_requirements, key=lambda item: item.condition_id)
        ),
        "temperature_k": values.get("temperature_k"),
        "pressure_atm": values.get("pressure_atm"),
        "minimum_expression_receipts": raw_minimum,
        "require_decision_evidence_binding": raw_binding,
    }
    return AnalysisCompletionPolicyV1(
        **body, policy_sha256=canonical_sha256(body)
    )


__all__ = [
    "AnalysisConditionV1",
    "ANALYSIS_COMPLETION_STAGES",
    "AnalysisIncompleteError",
    "AnalysisExpressionSourceRequirementV1",
    "AnalysisExpressionRequirementV1",
    "AnalysisClaimRequirementV1",
    "AnalysisCompletionReceiptV1",
    "AnalysisCompletionPolicyV1",
    "build_analysis_completion_receipt",
    "load_analysis_completion_policy",
]
