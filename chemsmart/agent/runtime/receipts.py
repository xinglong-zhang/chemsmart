"""Deterministic artifact receipts and lineage extraction."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

from chemsmart.agent.runtime.contracts import ArtifactRef

_PATH_KEYS = {
    "artifact",
    "generated_input",
    "input_file",
    "path",
    "submit_artifact",
    "written_path",
}


def artifact_ref(
    path: str | Path,
    *,
    kind: str,
    producer_tool: str = "",
    metadata: dict[str, Any] | None = None,
) -> ArtifactRef | None:
    candidate = Path(path).expanduser()
    if not candidate.is_file():
        return None
    content = candidate.read_bytes()
    return ArtifactRef(
        artifact_id=f"artifact_{uuid4().hex[:12]}",
        kind=kind,
        path=str(candidate.resolve()),
        sha256=hashlib.sha256(content).hexdigest(),
        size_bytes=len(content),
        producer_tool=producer_tool,
        metadata=dict(metadata or {}),
    )


def collect_artifact_refs(
    value: Any,
    *,
    producer_tool: str,
) -> tuple[ArtifactRef, ...]:
    found: list[ArtifactRef] = []
    seen: set[str] = set()

    def visit(node: Any, key: str = "") -> None:
        if isinstance(node, dict):
            for child_key, child in node.items():
                visit(child, str(child_key))
            return
        if isinstance(node, (list, tuple)):
            for child in node:
                visit(child, key)
            return
        if key not in _PATH_KEYS or not isinstance(node, str):
            return
        receipt = artifact_ref(
            node,
            kind=key,
            producer_tool=producer_tool,
        )
        if receipt is None or receipt.sha256 in seen:
            return
        seen.add(receipt.sha256)
        found.append(receipt)

    visit(value)
    return tuple(found)


def receipt_hashes(receipts: Iterable[ArtifactRef]) -> tuple[str, ...]:
    return tuple(receipt.sha256 for receipt in receipts)


__all__ = ["artifact_ref", "collect_artifact_refs", "receipt_hashes"]
