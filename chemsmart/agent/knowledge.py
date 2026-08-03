"""Evidence-bound candidate retrieval and explicit substitution decisions."""

from __future__ import annotations

from dataclasses import dataclass
from difflib import SequenceMatcher
from enum import Enum
from typing import Iterable, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.capabilities import (
    CapabilityQueryReceiptV1,
    CapabilityQueryStatus,
    ProgramCandidateProposalV1,
)


class SubstitutionDecision(str, Enum):
    EXACT = "exact"
    APPROVED = "approved"
    BLOCKED_AMBIGUOUS = "blocked_ambiguous"
    BLOCKED_NO_MATCH = "blocked_no_match"
    REJECTED_UNAPPROVED = "rejected_unapproved"


@dataclass(frozen=True)
class KnowledgeCandidateV1:
    candidate_id: str
    canonical_value: str
    normalized_value: str
    similarity: float
    source_id: str
    source_sha256: str
    exact: bool

    def __post_init__(self) -> None:
        if not self.candidate_id or not self.canonical_value:
            raise ContractError("candidate identity and value must not be empty")
        if not 0.0 <= self.similarity <= 1.0:
            raise ContractError("candidate similarity must be within [0, 1]")
        if self.exact != (self.similarity == 1.0):
            raise ContractError("candidate exact flag differs from similarity")
        require_sha256(self.source_sha256, "source_sha256")


@dataclass(frozen=True)
class CandidateQueryReceiptV1:
    schema_version: str
    query_id: str
    program: str
    setting_kind: str
    raw_query_sha256: str
    normalized_query: str
    candidate_source_sha256: str
    candidates: tuple[KnowledgeCandidateV1, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.candidate-query-receipt.v1":
            raise ContractError("unsupported candidate query receipt schema")
        require_identifier(self.program, "program")
        require_identifier(self.setting_kind, "setting_kind")
        order = tuple(
            sorted(
                self.candidates,
                key=lambda item: (-item.similarity, item.canonical_value),
            )
        )
        if order != self.candidates:
            raise ContractError("candidates must be deterministically ranked")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "query_id": self.query_id,
                "program": self.program,
                "setting_kind": self.setting_kind,
                "raw_query_sha256": self.raw_query_sha256,
                "normalized_query": self.normalized_query,
                "candidate_source_sha256": self.candidate_source_sha256,
                "candidates": self.candidates,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("candidate query receipt digest mismatch")


@dataclass(frozen=True)
class SubstitutionReceiptV1:
    schema_version: str
    candidate_receipt_sha256: str
    selected_candidate_id: str
    selected_value: str
    decision: SubstitutionDecision
    approval_ref: str
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.substitution-receipt.v1":
            raise ContractError("unsupported substitution receipt schema")
        if self.decision is SubstitutionDecision.APPROVED and not self.approval_ref:
            raise ContractError("fuzzy substitution requires an approval reference")
        expected = canonical_sha256(
            {
                "schema_version": self.schema_version,
                "candidate_receipt_sha256": self.candidate_receipt_sha256,
                "selected_candidate_id": self.selected_candidate_id,
                "selected_value": self.selected_value,
                "decision": self.decision,
                "approval_ref": self.approval_ref,
                "rule_ids": self.rule_ids,
            }
        )
        if self.receipt_sha256 != expected:
            raise ContractError("substitution receipt digest mismatch")


def normalize_setting_value(value: str) -> str:
    """Normalize only spelling noise; never change scientific meaning."""

    return "".join(character.lower() for character in str(value) if character.isalnum())


def rank_setting_candidates(
    *,
    query_id: str,
    program: str,
    setting_kind: str,
    raw_query: str,
    candidates: Mapping[str, tuple[str, str]],
    limit: int = 5,
) -> CandidateQueryReceiptV1:
    """Rank a trusted canonical vocabulary without selecting a replacement.

    ``candidates`` maps canonical values to ``(source_id, source_sha256)``.
    The raw query is represented by a digest in the receipt, while its
    normalized spelling is retained for deterministic rerendering.
    """

    if limit < 1:
        raise ContractError("candidate limit must be positive")
    normalized_query = normalize_setting_value(raw_query)
    if not normalized_query:
        raise ContractError("candidate query must contain alphanumeric text")
    source_rows = tuple(
        sorted(
            (value, source_id, source_sha256)
            for value, (source_id, source_sha256) in candidates.items()
        )
    )
    source_digest = canonical_sha256(source_rows)
    ranked = []
    for canonical_value, source_id, source_sha256 in source_rows:
        require_sha256(source_sha256, "candidate source_sha256")
        normalized = normalize_setting_value(canonical_value)
        similarity = SequenceMatcher(
            None, normalized_query, normalized, autojunk=False
        ).ratio()
        ranked.append(
            KnowledgeCandidateV1(
                candidate_id=canonical_sha256(
                    {
                        "program": program,
                        "setting_kind": setting_kind,
                        "canonical_value": canonical_value,
                        "source_sha256": source_sha256,
                    }
                )[:24],
                canonical_value=canonical_value,
                normalized_value=normalized,
                similarity=similarity,
                source_id=source_id,
                source_sha256=source_sha256,
                exact=similarity == 1.0,
            )
        )
    ranked.sort(key=lambda item: (-item.similarity, item.canonical_value))
    body = {
        "schema_version": "chemsmart.candidate-query-receipt.v1",
        "query_id": str(query_id),
        "program": require_identifier(program, "program"),
        "setting_kind": require_identifier(setting_kind, "setting_kind"),
        "raw_query_sha256": canonical_sha256({"raw_query": raw_query}),
        "normalized_query": normalized_query,
        "candidate_source_sha256": source_digest,
        "candidates": tuple(ranked[:limit]),
    }
    return CandidateQueryReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def decide_substitution(
    receipt: CandidateQueryReceiptV1,
    *,
    selected_candidate_id: str = "",
    approval_ref: str = "",
) -> SubstitutionReceiptV1:
    """Resolve an exact match or record an explicit fuzzy-choice boundary."""

    exact = [item for item in receipt.candidates if item.exact]
    selected = next(
        (
            item
            for item in receipt.candidates
            if item.candidate_id == selected_candidate_id
        ),
        None,
    )
    if len(exact) == 1 and (
        not selected_candidate_id
        or selected_candidate_id == exact[0].candidate_id
    ):
        selected = exact[0]
        decision = SubstitutionDecision.EXACT
        rules = ("knowledge.candidate.exact",)
    elif selected is not None and approval_ref:
        decision = SubstitutionDecision.APPROVED
        rules = ("knowledge.substitution.explicitly_approved",)
    elif selected is not None:
        decision = SubstitutionDecision.REJECTED_UNAPPROVED
        rules = ("knowledge.substitution.approval_required",)
    elif not receipt.candidates:
        decision = SubstitutionDecision.BLOCKED_NO_MATCH
        rules = ("knowledge.candidate.none",)
    else:
        decision = SubstitutionDecision.BLOCKED_AMBIGUOUS
        rules = ("knowledge.candidate.ambiguous",)
    body = {
        "schema_version": "chemsmart.substitution-receipt.v1",
        "candidate_receipt_sha256": receipt.receipt_sha256,
        "selected_candidate_id": selected.candidate_id if selected else "",
        "selected_value": selected.canonical_value if selected else "",
        "decision": decision,
        "approval_ref": approval_ref if decision is SubstitutionDecision.APPROVED else "",
        "rule_ids": rules,
    }
    return SubstitutionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__ = [
    "CandidateQueryReceiptV1",
    "KnowledgeCandidateV1",
    "SubstitutionDecision",
    "SubstitutionReceiptV1",
    "ProgramSubstitutionReceiptV1",
    "decide_substitution",
    "normalize_setting_value",
    "rank_setting_candidates",
]


@dataclass(frozen=True)
class ProgramSubstitutionReceiptV1:
    schema_version: str
    proposal_sha256: str
    substitution_request_sha256: str
    capability_receipt_sha256: str
    requested_program: str
    selected_program: str
    requested_engine: str
    selected_engine: str
    decision: str
    approval_ref: str
    readiness_authority: bool
    rule_ids: tuple[str, ...]
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-substitution-receipt.v1":
            raise ContractError("unsupported program substitution schema")
        if self.readiness_authority:
            raise ContractError("substitution assessment cannot establish readiness")
        if self.decision not in {
            "exact",
            "approved",
            "blocked",
            "rejected",
        }:
            raise ContractError("invalid program substitution decision")
        if self.decision == "approved" and not self.approval_ref:
            raise ContractError("program substitution requires explicit approval")
        body = {
            "schema_version": self.schema_version,
            "proposal_sha256": self.proposal_sha256,
            "substitution_request_sha256": (
                self.substitution_request_sha256
            ),
            "capability_receipt_sha256": self.capability_receipt_sha256,
            "requested_program": self.requested_program,
            "selected_program": self.selected_program,
            "requested_engine": self.requested_engine,
            "selected_engine": self.selected_engine,
            "decision": self.decision,
            "approval_ref": self.approval_ref,
            "readiness_authority": self.readiness_authority,
            "rule_ids": self.rule_ids,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("program substitution receipt digest mismatch")


def assess_program_substitution(
    proposal: ProgramCandidateProposalV1,
    capability: CapabilityQueryReceiptV1,
    *,
    approval_ref: str = "",
) -> ProgramSubstitutionReceiptV1:
    matches = (
        proposal.candidate_program == capability.query.program
        and (
            not proposal.candidate_engine
            or proposal.candidate_engine == capability.query.engine
        )
    )
    supported = capability.status in {
        CapabilityQueryStatus.SUPPORTED,
        CapabilityQueryStatus.PREVIEW_ONLY,
    }
    if not matches or not supported:
        decision = "blocked"
        rules = ("program.substitution.capability_red",)
    elif proposal.requested_program == proposal.candidate_program:
        decision = "exact"
        rules = ("program.substitution.identity",)
    else:
        decision = "blocked"
        rules = ("program.substitution.typed_eligibility_required",)
    body = {
        "schema_version": "chemsmart.program-substitution-receipt.v1",
        "proposal_sha256": proposal.proposal_sha256,
        "substitution_request_sha256": "",
        "capability_receipt_sha256": capability.receipt_sha256,
        "requested_program": proposal.requested_program,
        "selected_program": proposal.candidate_program if matches else "",
        "requested_engine": "",
        "selected_engine": proposal.candidate_engine if matches else "",
        "decision": decision,
        "approval_ref": approval_ref if decision == "approved" else "",
        "readiness_authority": False,
        "rule_ids": rules,
    }
    return ProgramSubstitutionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__.append("assess_program_substitution")


_PYSCF_TRANSFER_JOB_FAMILIES = frozenset(
    {"sp", "opt", "hess", "freq", "opt_freq"}
)
_PYSCF_FORBIDDEN_JOB_FAMILIES = frozenset(
    {"ts", "irc", "td", "scan", "qmmm", "neb"}
)


@dataclass(frozen=True)
class ScientificClaimEvidenceV1:
    """Host-verified literature claim; a bare digest is never evidence."""

    schema_version: str
    claim_sha256: str
    source_artifact_sha256: str
    source_locator: str
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-claim-evidence.v1":
            raise ContractError("unsupported scientific claim evidence schema")
        require_sha256(self.claim_sha256, "claim_sha256")
        require_sha256(self.source_artifact_sha256, "source_artifact_sha256")
        if not self.source_locator:
            raise ContractError("scientific claim requires a source locator")
        if self.status not in {"verified", "rejected"}:
            raise ContractError("invalid scientific claim evidence status")
        body = {
            "schema_version": self.schema_version,
            "claim_sha256": self.claim_sha256,
            "source_artifact_sha256": self.source_artifact_sha256,
            "source_locator": self.source_locator,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("scientific claim evidence digest mismatch")


@dataclass(frozen=True)
class FunctionalEquivalenceReceiptV1:
    """Independent host verdict for a cross-program method mapping."""

    schema_version: str
    requested_program: str
    selected_program: str
    method_name: str
    claim_evidence_receipt_sha256s: tuple[str, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.functional-equivalence-receipt.v1":
            raise ContractError("unsupported functional equivalence schema")
        require_identifier(self.requested_program, "requested_program")
        require_identifier(self.selected_program, "selected_program")
        if not self.method_name:
            raise ContractError("functional equivalence requires a method")
        if self.claim_evidence_receipt_sha256s != tuple(
            sorted(set(self.claim_evidence_receipt_sha256s))
        ):
            raise ContractError("equivalence evidence must be sorted and unique")
        for digest in self.claim_evidence_receipt_sha256s:
            require_sha256(digest, "claim evidence receipt")
        if self.status not in {"verified", "rejected"}:
            raise ContractError("invalid functional equivalence status")
        body = {
            "schema_version": self.schema_version,
            "requested_program": self.requested_program,
            "selected_program": self.selected_program,
            "method_name": self.method_name,
            "claim_evidence_receipt_sha256s": (
                self.claim_evidence_receipt_sha256s
            ),
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("functional equivalence receipt digest mismatch")


@dataclass(frozen=True)
class ProgramSubstitutionApprovalV1:
    """Non-model approval bound to one exact substitution request."""

    schema_version: str
    substitution_request_sha256: str
    actor: str
    decision: str
    approval_id: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-substitution-approval.v1":
            raise ContractError("unsupported substitution approval schema")
        require_sha256(
            self.substitution_request_sha256,
            "substitution_request_sha256",
        )
        require_identifier(self.approval_id, "approval_id")
        if not self.actor or self.decision not in {"approved", "rejected"}:
            raise ContractError("invalid substitution approval")
        body = {
            "schema_version": self.schema_version,
            "substitution_request_sha256": self.substitution_request_sha256,
            "actor": self.actor,
            "decision": self.decision,
            "approval_id": self.approval_id,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("substitution approval digest mismatch")


@dataclass(frozen=True)
class ProgramSubstitutionRequestV1:
    """Typed scientific eligibility request for changing program authority."""

    schema_version: str
    request_id: str
    requested_program: str
    selected_program: str
    requested_engine: str
    selected_engine: str
    job_families: tuple[str, ...]
    method_family: str
    method_name: str
    basis_mode: str
    constraint_kinds: tuple[str, ...]
    requires_post_hf: bool
    requires_double_hybrid: bool
    functional_semantics_confirmed: bool
    source_claim_sha256s: tuple[str, ...]
    request_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.program-substitution-request.v1":
            raise ContractError("unsupported program substitution request")
        require_identifier(self.request_id, "request_id")
        require_identifier(self.requested_program, "requested_program")
        require_identifier(self.selected_program, "selected_program")
        require_identifier(self.requested_engine, "requested_engine")
        require_identifier(self.selected_engine, "selected_engine")
        if self.job_families != tuple(sorted(set(self.job_families))):
            raise ContractError("substitution job families must be sorted")
        if not self.job_families:
            raise ContractError("substitution requires a job family")
        if self.basis_mode not in {"uniform", "mixed", "ecp", "mixed_ecp"}:
            raise ContractError("unsupported substitution basis mode")
        if self.constraint_kinds != tuple(sorted(set(self.constraint_kinds))):
            raise ContractError("constraint kinds must be sorted and unique")
        for digest in self.source_claim_sha256s:
            require_sha256(digest, "source claim digest")
        body = {
            "schema_version": self.schema_version,
            "request_id": self.request_id,
            "requested_program": self.requested_program,
            "selected_program": self.selected_program,
            "requested_engine": self.requested_engine,
            "selected_engine": self.selected_engine,
            "job_families": self.job_families,
            "method_family": self.method_family,
            "method_name": self.method_name,
            "basis_mode": self.basis_mode,
            "constraint_kinds": self.constraint_kinds,
            "requires_post_hf": self.requires_post_hf,
            "requires_double_hybrid": self.requires_double_hybrid,
            "functional_semantics_confirmed": (
                self.functional_semantics_confirmed
            ),
            "source_claim_sha256s": self.source_claim_sha256s,
        }
        if self.request_sha256 != canonical_sha256(body):
            raise ContractError("program substitution request digest mismatch")


def build_program_substitution_request(
    *,
    request_id: str,
    requested_program: str,
    selected_program: str,
    requested_engine: str,
    selected_engine: str,
    job_families: Iterable[str],
    method_family: str,
    method_name: str,
    basis_mode: str = "uniform",
    constraint_kinds: Iterable[str] = (),
    requires_post_hf: bool = False,
    requires_double_hybrid: bool = False,
    functional_semantics_confirmed: bool = False,
    source_claim_sha256s: Iterable[str] = (),
) -> ProgramSubstitutionRequestV1:
    body = {
        "schema_version": "chemsmart.program-substitution-request.v1",
        "request_id": require_identifier(request_id, "request_id"),
        "requested_program": require_identifier(
            requested_program, "requested_program"
        ),
        "selected_program": require_identifier(
            selected_program, "selected_program"
        ),
        "requested_engine": require_identifier(
            requested_engine, "requested_engine"
        ),
        "selected_engine": require_identifier(
            selected_engine, "selected_engine"
        ),
        "job_families": tuple(sorted(set(str(item) for item in job_families))),
        "method_family": require_identifier(method_family, "method_family"),
        "method_name": str(method_name).strip(),
        "basis_mode": str(basis_mode).strip().lower(),
        "constraint_kinds": tuple(
            sorted(set(str(item).strip().lower() for item in constraint_kinds))
        ),
        "requires_post_hf": bool(requires_post_hf),
        "requires_double_hybrid": bool(requires_double_hybrid),
        "functional_semantics_confirmed": bool(
            functional_semantics_confirmed
        ),
        "source_claim_sha256s": tuple(sorted(set(source_claim_sha256s))),
    }
    return ProgramSubstitutionRequestV1(
        **body, request_sha256=canonical_sha256(body)
    )


def assess_typed_program_substitution(
    request: ProgramSubstitutionRequestV1,
    capability: CapabilityQueryReceiptV1,
    *,
    approval_ref: str = "",
) -> ProgramSubstitutionReceiptV1:
    """Apply the closed Gaussian-to-PySCF transfer matrix.

    The assessment is deliberately advisory.  Project-loader, environment,
    preflight, preview, and result-verifier receipts remain mandatory.
    """

    failures = []
    if request.selected_program != capability.query.program:
        failures.append("program.substitution.selected_program_mismatch")
    if request.selected_engine != capability.query.engine:
        failures.append("program.substitution.selected_engine_mismatch")
    if capability.status not in {
        CapabilityQueryStatus.SUPPORTED,
        CapabilityQueryStatus.PREVIEW_ONLY,
    }:
        failures.append("program.substitution.capability_red")
    identity = request.requested_program == request.selected_program
    if not identity:
        if (request.requested_program, request.selected_program) != (
            "gaussian",
            "pyscf",
        ):
            failures.append("program.substitution.matrix_pair_unsupported")
        if not set(request.job_families).issubset(
            _PYSCF_TRANSFER_JOB_FAMILIES
        ):
            failures.append("program.substitution.job_family_unsupported")
        if set(request.job_families) & _PYSCF_FORBIDDEN_JOB_FAMILIES:
            failures.append("program.substitution.job_family_forbidden")
        if request.method_family not in {"hf", "dft"}:
            failures.append("program.substitution.method_family_unsupported")
        if request.requires_post_hf:
            failures.append("program.substitution.post_hf_unsupported")
        if request.requires_double_hybrid:
            failures.append("program.substitution.double_hybrid_unsupported")
        if request.basis_mode != "uniform":
            failures.append("program.substitution.mixed_basis_ecp_unsupported")
        if request.constraint_kinds:
            failures.append("program.substitution.constraints_unsupported")
        if (
            request.method_family == "dft"
            and not request.functional_semantics_confirmed
        ):
            failures.append("program.substitution.functional_semantics_unknown")
        if not request.source_claim_sha256s:
            failures.append("program.substitution.source_evidence_missing")

    if failures:
        decision = "blocked"
        rules = tuple(sorted(set(failures)))
    elif identity:
        decision = "exact"
        rules = ("program.substitution.identity",)
    elif not approval_ref:
        decision = "rejected"
        rules = ("program.substitution.approval_required",)
    else:
        decision = "approved"
        rules = (
            "program.substitution.gaussian_to_pyscf_matrix",
            "program.substitution.explicitly_approved",
        )
    body = {
        "schema_version": "chemsmart.program-substitution-receipt.v1",
        "proposal_sha256": "",
        "substitution_request_sha256": request.request_sha256,
        "capability_receipt_sha256": capability.receipt_sha256,
        "requested_program": request.requested_program,
        "selected_program": request.selected_program,
        "requested_engine": request.requested_engine,
        "selected_engine": request.selected_engine,
        "decision": decision,
        "approval_ref": approval_ref if decision == "approved" else "",
        "readiness_authority": False,
        "rule_ids": rules,
    }
    return ProgramSubstitutionReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__.extend(
    [
        "FunctionalEquivalenceReceiptV1",
        "ProgramSubstitutionRequestV1",
        "ProgramSubstitutionApprovalV1",
        "ScientificClaimEvidenceV1",
        "assess_typed_program_substitution",
        "build_program_substitution_request",
    ]
)
