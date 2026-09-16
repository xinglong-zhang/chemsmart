"""A cohort: the calculations one wave runs, named once and not again.

A *wave* is what the Agent chose to see before it reasons again -- a set
of currently-ready, scientifically independent approved calculations. The
host executes them concurrently and wakes the Agent once, when every
member has reached a terminal state, however it reached one.

The manifest is the only thing that maps a scheduler array index to an
approved calculation. An array index is scheduler representation and
never scientific identity: nothing durable keys a calculation by it, so a
reader asking what element 3 computed is answered with a node id.

Membership is fixed at dispatch. A set that could grow after submission
would be a different experiment than the one the barrier is waiting for.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)

#: Written in the run directory when a cycle dispatches a wave.
COHORT_MANIFEST_FILE = "cohort.json"


@dataclass(frozen=True)
class CohortManifestV1:
    """One wave's members, in the order their elements are numbered."""

    schema_version: str
    goal_id: str
    cycle: int
    bundle_sha256: str
    node_ids: tuple[str, ...]
    max_concurrent_tasks: int
    created_at: str
    cohort_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.cohort-manifest.v1":
            raise ContractError("unsupported cohort manifest schema")
        require_identifier(self.goal_id, "goal_id")
        require_sha256(self.bundle_sha256, "bundle_sha256")
        if int(self.cycle) < 1:
            raise ContractError("a cohort belongs to a cycle")
        if not self.node_ids:
            raise ContractError("a cohort has at least one calculation")
        for node_id in self.node_ids:
            require_identifier(node_id, "node_id")
        if len(set(self.node_ids)) != len(self.node_ids):
            raise ContractError(
                "a calculation appears in a cohort once: an element is a "
                "position and a node is an identity"
            )
        if int(self.max_concurrent_tasks) < 1:
            raise ContractError("concurrency is at least one")
        digest = canonical_sha256(self._body())
        if self.cohort_sha256 != digest:
            raise ContractError("cohort manifest digest mismatch")

    def _body(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "goal_id": self.goal_id,
            "cycle": int(self.cycle),
            "bundle_sha256": self.bundle_sha256,
            "node_ids": list(self.node_ids),
            "max_concurrent_tasks": int(self.max_concurrent_tasks),
            "created_at": self.created_at,
        }

    @property
    def cohort_id(self) -> str:
        return self.cohort_sha256

    @property
    def size(self) -> int:
        return len(self.node_ids)

    def node_for_element(self, element: int) -> str:
        """The approved calculation this array element runs.

        Args:
            element: The scheduler's task index, from zero.

        Returns:
            str: The node id.

        Raises:
            ContractError: If the index is outside this cohort.
        """

        try:
            index = int(element)
        except (TypeError, ValueError) as exc:
            raise ContractError(
                f"element {element!r} is not an index"
            ) from exc
        if index < 0 or index >= len(self.node_ids):
            raise ContractError(
                f"element {index} is outside this cohort of "
                f"{len(self.node_ids)}"
            )
        return self.node_ids[index]

    def element_for_node(self, node_id: str) -> int:
        try:
            return self.node_ids.index(str(node_id))
        except ValueError as exc:
            raise ContractError(
                f"{node_id!r} is not a member of this cohort"
            ) from exc

    def public_record(self) -> dict[str, Any]:
        return {**self._body(), "cohort_sha256": self.cohort_sha256}

    def write(self, run_directory: str | Path) -> Path:
        path = Path(run_directory) / COHORT_MANIFEST_FILE
        path.write_text(
            json.dumps(self.public_record(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return path


def cohort_frontier(
    ready: tuple[str, ...],
    cohort_node_ids: tuple[str, ...] | None,
) -> tuple[str, ...]:
    """The ready nodes this wave may run: its own members and no others.

    The executor's walk maximises throughput -- run whatever is ready --
    and that is the wrong shape here. The Agent asked for these
    calculations as one scientific experiment, and their *collective*
    evidence is what triggers the next reasoning turn. A node whose
    dependency cleared mid-wave is a decision the Agent has not made yet,
    so it waits for the wake even though the scheduler could start it.

    ``None`` means no cohort: a single-job dispatch and every run
    recorded before waves existed keep the walk they had.

    Args:
        ready: What the DAG authority says is runnable, in its order.
        cohort_node_ids: This wave's members, or None.

    Returns:
        tuple[str, ...]: The subset to run, in the order offered.
    """

    if cohort_node_ids is None:
        return tuple(ready)
    members = set(cohort_node_ids)
    return tuple(node_id for node_id in ready if node_id in members)


def build_cohort_manifest(
    *,
    goal_id: str,
    cycle: int,
    bundle_sha256: str,
    node_ids: tuple[str, ...],
    max_concurrent_tasks: int,
    created_at: str,
) -> CohortManifestV1:
    body = {
        "schema_version": "chemsmart.cohort-manifest.v1",
        "goal_id": str(goal_id),
        "cycle": int(cycle),
        "bundle_sha256": str(bundle_sha256),
        "node_ids": [str(item) for item in node_ids],
        "max_concurrent_tasks": int(max_concurrent_tasks),
        "created_at": str(created_at),
    }
    return CohortManifestV1(
        schema_version=body["schema_version"],
        goal_id=body["goal_id"],
        cycle=body["cycle"],
        bundle_sha256=body["bundle_sha256"],
        node_ids=tuple(body["node_ids"]),
        max_concurrent_tasks=body["max_concurrent_tasks"],
        created_at=body["created_at"],
        cohort_sha256=canonical_sha256(body),
    )


def read_cohort_manifest(
    run_directory: str | Path,
) -> CohortManifestV1 | None:
    """The cohort this run directory dispatched, or None.

    A single-job dispatch writes no manifest, and that is absence rather
    than failure. A manifest that does not reconstruct its own digest is
    a failure: membership is what the barrier waits on.
    """

    path = Path(run_directory) / COHORT_MANIFEST_FILE
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(raw, dict):
        return None
    return CohortManifestV1(
        schema_version=str(raw.get("schema_version") or ""),
        goal_id=str(raw.get("goal_id") or ""),
        cycle=int(raw.get("cycle") or 0),
        bundle_sha256=str(raw.get("bundle_sha256") or ""),
        node_ids=tuple(str(item) for item in raw.get("node_ids") or ()),
        max_concurrent_tasks=int(raw.get("max_concurrent_tasks") or 0),
        created_at=str(raw.get("created_at") or ""),
        cohort_sha256=str(raw.get("cohort_sha256") or ""),
    )


__all__ = [
    "COHORT_MANIFEST_FILE",
    "CohortManifestV1",
    "build_cohort_manifest",
    "cohort_frontier",
    "read_cohort_manifest",
]
