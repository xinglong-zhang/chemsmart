"""Host-rendered numerical claims over deterministic analysis receipts."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

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
    #: What the session attributes as this number's uncertainty, in the
    #: display unit, and whether it was measured here, inferred from a
    #: cited receipt, or asserted. Part of the claim rather than beside
    #: it: a number and the uncertainty its author gives it are one
    #: statement, so a corrected uncertainty is a different claim with a
    #: different receipt. While it lived outside the record, re-claiming
    #: with a corrected uncertainty produced the same digest, the same
    #: idempotency key and a different payload, and the event store
    #: refused it -- which made "claim it again with the uncertainty you
    #: measured", the route the sufficiency wake offers, unwalkable.
    uncertainty: float | None = None
    uncertainty_basis: str = ""
    #: Where the magnitude came from, and the terms it was built from.
    #: The host resolved the reference against its own receipt
    #: registries before admitting the word `measured` or `inferred`,
    #: and then dropped it: the only durable trace was one boolean, so
    #: nothing afterwards -- a human, a later cycle, an auditor --
    #: could ask which receipt backed the number that discharged the
    #: contract. They ride inside the digest for the same reason the
    #: uncertainty does: re-claiming the same value with a corrected
    #: budget must be a different receipt, not an idempotency
    #: collision.
    uncertainty_reference: str = ""
    uncertainty_components: tuple[Any, ...] = ()
    #: What the host observed about the cited magnitude while resolving
    #: it, and did not rule on. Three of these were refusals until the
    #: boundary was drawn (owner ruling, 2026-09-10): a chain naming a
    #: number the session supplied, a spread over a single receipt, and
    #: a magnitude of exactly zero. Each rejected legitimate science --
    #: a variance over many samples in one receipt, an equality symmetry
    #: enforces, a coefficient a definition fixes -- and none prevented
    #: the substitution it was aimed at, because an equivalent spelling
    #: walks past all three. So the host reports what it saw and the
    #: session owns what it means. They ride inside the digest with the
    #: reference and the components: what the host observed about a
    #: number is part of the statement that number appears in.
    uncertainty_observations: tuple[str, ...] = ()
    #: How the session combined its components into the total, what its
    #: inputs mean, what coverage the total claims, and what it assumes
    #: about dependence between the terms. The host copies it, renders
    #: it beside the verdict, and grades none of it: how the terms add
    #: is the science.
    #:
    #: Three windows in a row had the combination rule, not the
    #: evidence, decide the word: identical components gave a
    #: root-sum-square inside the tolerance and a linear sum outside it,
    #: and nothing in the assessment, the record or the settlement named
    #: which had produced the number a human read. Detecting a `sqrt`
    #: node would be the syntactic version -- an equivalent spelling
    #: walks past it -- and any host arithmetic relating a total to its
    #: components is forbidden, because magnitudes are signed and
    #: correlated terms may legitimately cancel. So the session declares
    #: it and the host carries it (SUFFICIENCY-5, 2026-09-10).
    uncertainty_combination: Mapping[str, Any] | None = None
    #: When this number approximates the declared quantity rather than
    #: being it: which declaration, by what relationship, on what
    #: basis. The relationship is the session's own words -- the host
    #: ships no closed vocabulary of chemistry kinds and never infers
    #: one from an identifier.
    #:
    #: A live claim delivered `quartet-minus-doublet-u0k-kjmol` against
    #: a declaration of `delta-g-quartet-minus-doublet`: a zero-point
    #: corrected electronic difference, evaluated on a structure the
    #: host had itself typed a first-order saddle, answering a
    #: declaration about a Gibbs difference between minima. Same id,
    #: same dimension, different quantity -- and the expectation row
    #: printed `agreed`. The number stays delivered, because a refusal
    #: that buries a finding is a defect in the refusal; what stops is
    #: the host asserting an unqualified agreement over an
    #: approximation (SUFFICIENCY-5, 2026-09-10).
    approximates: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        require_identifier(self.claim_id, "claim_id")
        if self.uncertainty_reference and self.uncertainty is None:
            raise ContractError(
                "an uncertainty_reference belongs to an uncertainty"
            )
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
        if self.uncertainty is not None:
            if float(self.uncertainty) < 0.0:
                raise ContractError(
                    "an uncertainty is a magnitude in the claim's display "
                    "unit and cannot be negative"
                )
            if self.uncertainty_basis not in {
                "measured",
                "inferred",
                "asserted",
            }:
                raise ContractError(
                    "an uncertainty needs uncertainty_basis: measured, "
                    "inferred, or asserted"
                )
        if len(self.dimension) not in {6, 7, 8, 9} or not all(
            isinstance(value, int) for value in self.dimension
        ):
            raise ContractError(
                "analysis claim dimension must contain six legacy, seven "
                "dipole-extended, eight mass-extended, or nine "
                "charge-extended integers"
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
            raise ContractError(
                "analysis claim record status must be recorded"
            )
        body = analysis_claim_record_body(self)
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("analysis claim record digest mismatch")


def analysis_claim_record_body(
    record: AnalysisClaimRecordV1,
) -> dict[str, Any]:
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
    return AnalysisClaimRecordV1(**body, receipt_sha256=canonical_sha256(body))


def analysis_claim_record_from_record(
    record: dict[str, Any], *, receipt_sha256: str
) -> AnalysisClaimRecordV1:
    """Rehydrate a host-rendered claim record from a Runtime V2 event."""

    values = dict(record)
    claims = []
    for item in values.get("claims") or ():
        claim = dict(item)
        claim["dimension"] = tuple(claim.get("dimension") or ())
        claim["uncertainty_components"] = tuple(
            claim.get("uncertainty_components") or ()
        )
        claims.append(AnalysisReportedQuantityV1(**claim))
    values["claims"] = tuple(claims)
    return AnalysisClaimRecordV1(**values, receipt_sha256=receipt_sha256)


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
