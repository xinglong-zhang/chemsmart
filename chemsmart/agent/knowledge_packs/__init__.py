"""Advisory program knowledge activated independently from readiness."""

from __future__ import annotations

from typing import Iterable

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.knowledge_packs.builtin import BUILTIN_PROGRAM_PACKS
from chemsmart.agent.knowledge_packs.contracts import (
    AdvisoryProgramKnowledgePackV1,
    KnowledgePackActivationReceiptV1,
)


def activate_program_knowledge(
    *,
    request: str,
    program: str,
    engine: str = "",
    packs: Iterable[AdvisoryProgramKnowledgePackV1] = BUILTIN_PROGRAM_PACKS,
) -> KnowledgePackActivationReceiptV1:
    text = str(request).lower()
    activated = []
    excluded = []
    for pack in packs:
        if pack.target_program != str(program).lower():
            continue
        if pack.target_engine and pack.target_engine != str(engine).lower():
            continue
        if any(term in text for term in pack.exclusions):
            excluded.append(pack.pack_sha256)
            continue
        if any(term in text for term in pack.activation_terms):
            activated.append(pack.pack_sha256)
    body = {
        "schema_version": "chemsmart.knowledge-pack-activation.v1",
        "request_sha256": canonical_sha256({"request": request}),
        "program": str(program).lower(),
        "engine": str(engine).lower(),
        "activated_pack_sha256s": tuple(sorted(activated)),
        "excluded_pack_sha256s": tuple(sorted(excluded)),
        "advisory_only": True,
    }
    return KnowledgePackActivationReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


__all__ = [
    "AdvisoryProgramKnowledgePackV1",
    "BUILTIN_PROGRAM_PACKS",
    "KnowledgePackActivationReceiptV1",
    "activate_program_knowledge",
]
