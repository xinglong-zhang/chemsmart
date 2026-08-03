"""Versioned advisory knowledge packs with no readiness authority."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import ContractError, canonical_sha256


@dataclass(frozen=True)
class AdvisoryProgramKnowledgePackV1:
    schema_version: str
    pack_id: str
    pack_version: str
    target_program: str
    target_engine: str
    activation_terms: tuple[str, ...]
    exclusions: tuple[str, ...]
    advisory_topics: tuple[str, ...]
    source_ids: tuple[str, ...]
    reference_only: bool
    readiness_authority: bool
    pack_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.advisory-program-pack.v1":
            raise ContractError("unsupported advisory knowledge pack schema")
        if self.readiness_authority:
            raise ContractError("knowledge packs cannot establish readiness")
        for values, name in (
            (self.activation_terms, "activation_terms"),
            (self.exclusions, "exclusions"),
            (self.advisory_topics, "advisory_topics"),
            (self.source_ids, "source_ids"),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} must be sorted and unique")
        body = {
            "schema_version": self.schema_version,
            "pack_id": self.pack_id,
            "pack_version": self.pack_version,
            "target_program": self.target_program,
            "target_engine": self.target_engine,
            "activation_terms": self.activation_terms,
            "exclusions": self.exclusions,
            "advisory_topics": self.advisory_topics,
            "source_ids": self.source_ids,
            "reference_only": self.reference_only,
            "readiness_authority": self.readiness_authority,
        }
        if self.pack_sha256 != canonical_sha256(body):
            raise ContractError("advisory knowledge pack digest mismatch")


@dataclass(frozen=True)
class KnowledgePackActivationReceiptV1:
    schema_version: str
    request_sha256: str
    program: str
    engine: str
    activated_pack_sha256s: tuple[str, ...]
    excluded_pack_sha256s: tuple[str, ...]
    advisory_only: bool
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.knowledge-pack-activation.v1":
            raise ContractError("unsupported knowledge activation schema")
        if not self.advisory_only:
            raise ContractError("knowledge activation must remain advisory")
        body = {
            "schema_version": self.schema_version,
            "request_sha256": self.request_sha256,
            "program": self.program,
            "engine": self.engine,
            "activated_pack_sha256s": self.activated_pack_sha256s,
            "excluded_pack_sha256s": self.excluded_pack_sha256s,
            "advisory_only": self.advisory_only,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("knowledge activation receipt digest mismatch")


def build_pack(**values: object) -> AdvisoryProgramKnowledgePackV1:
    body = {
        "schema_version": "chemsmart.advisory-program-pack.v1",
        "pack_id": str(values["pack_id"]),
        "pack_version": str(values["pack_version"]),
        "target_program": str(values["target_program"]),
        "target_engine": str(values.get("target_engine", "")),
        "activation_terms": tuple(
            sorted(set(str(item).lower() for item in values["activation_terms"]))
        ),
        "exclusions": tuple(
            sorted(set(str(item).lower() for item in values["exclusions"]))
        ),
        "advisory_topics": tuple(
            sorted(set(str(item) for item in values["advisory_topics"]))
        ),
        "source_ids": tuple(
            sorted(set(str(item) for item in values["source_ids"]))
        ),
        "reference_only": bool(values.get("reference_only", False)),
        "readiness_authority": False,
    }
    return AdvisoryProgramKnowledgePackV1(
        **body, pack_sha256=canonical_sha256(body)
    )


__all__ = [
    "AdvisoryProgramKnowledgePackV1",
    "KnowledgePackActivationReceiptV1",
    "build_pack",
]
