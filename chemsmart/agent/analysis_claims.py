"""Host-rendered numerical claims over deterministic analysis receipts."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


@dataclass(frozen=True)
class AnalysisReportedQuantityV1:
    """One claim whose value is copied by the host from a typed quantity."""

    claim_id: str
    source_kind: str
    source_receipt_sha256: str
    quantity_id: str
    quantity_value_sha256: str
    display_value: Any
    display_unit: str
    canonical_value: Any
    canonical_unit: str
    dimension: tuple[int, ...]
    data_kind: str

    def __post_init__(self) -> None:
        require_identifier(self.claim_id, "claim_id")
        if self.source_kind not in {
            "quantity_extraction",
            "thermochemistry",
            "quantity_expression",
            "scientific_validation",
        }:
            raise ContractError("unsupported analysis claim source kind")
        require_sha256(self.source_receipt_sha256, "source_receipt_sha256")
        require_identifier(self.quantity_id, "quantity_id")
        require_sha256(self.quantity_value_sha256, "quantity_value_sha256")
        if len(self.dimension) not in {6, 7, 8} or not all(
            isinstance(value, int) for value in self.dimension
        ):
            raise ContractError(
                "analysis claim dimension must contain six legacy, seven "
                "dipole-extended, or eight mass-extended integers"
            )
        _require_finite_payload(self.display_value, "display_value")
        _require_finite_payload(self.canonical_value, "canonical_value")
        if not self.display_unit or not self.canonical_unit:
            raise ContractError("analysis claim units must not be empty")


@dataclass(frozen=True)
class AnalysisClaimRecordV1:
    """Task-bound set of exact numerical claims for reporting."""

    schema_version: str
    task_spec_sha256: str
    claims: tuple[AnalysisReportedQuantityV1, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.analysis-claim-record.v1":
            raise ContractError("unsupported analysis claim record schema")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if not self.claims:
            raise ContractError("analysis claim record requires claims")
        claim_ids = tuple(claim.claim_id for claim in self.claims)
        if claim_ids != tuple(sorted(set(claim_ids))):
            raise ContractError("analysis claims must be sorted and unique")
        if self.status != "recorded":
            raise ContractError("analysis claim record status must be recorded")
        body = analysis_claim_record_body(self)
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("analysis claim record digest mismatch")


def analysis_claim_record_body(record: AnalysisClaimRecordV1) -> dict[str, Any]:
    return {
        "schema_version": record.schema_version,
        "task_spec_sha256": record.task_spec_sha256,
        "claims": record.claims,
        "status": record.status,
    }


def build_analysis_claim_record(
    *,
    task_spec_sha256: str,
    claims: tuple[AnalysisReportedQuantityV1, ...],
) -> AnalysisClaimRecordV1:
    body = {
        "schema_version": "chemsmart.analysis-claim-record.v1",
        "task_spec_sha256": task_spec_sha256,
        "claims": tuple(sorted(claims, key=lambda claim: claim.claim_id)),
        "status": "recorded",
    }
    return AnalysisClaimRecordV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def analysis_claim_record_from_record(
    record: dict[str, Any], *, receipt_sha256: str
) -> AnalysisClaimRecordV1:
    """Rehydrate a host-rendered claim record from a Runtime V2 event."""

    values = dict(record)
    claims = []
    for item in values.get("claims") or ():
        claim = dict(item)
        claim["dimension"] = tuple(claim.get("dimension") or ())
        claims.append(AnalysisReportedQuantityV1(**claim))
    values["claims"] = tuple(claims)
    return AnalysisClaimRecordV1(
        **values, receipt_sha256=receipt_sha256
    )


def _require_finite_payload(value: Any, field: str) -> None:
    if isinstance(value, bool) or value is None or isinstance(value, str):
        raise ContractError(f"{field} must be numerical")
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise ContractError(f"{field} must be finite")
        return
    if isinstance(value, (tuple, list)) and value:
        for item in value:
            _require_finite_payload(item, field)
        return
    raise ContractError(f"{field} must be a finite numerical value or array")


__all__ = [
    "AnalysisClaimRecordV1",
    "AnalysisReportedQuantityV1",
    "analysis_claim_record_body",
    "analysis_claim_record_from_record",
    "build_analysis_claim_record",
]
