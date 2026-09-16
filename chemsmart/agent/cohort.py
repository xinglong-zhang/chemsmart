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
        """Write this cohort once.

        Membership is fixed at dispatch. An unconditional overwrite let a
        second correctly-digested manifest replace the first through the
        public writer -- a different experiment than the one the barrier
        is waiting for, with the digest check none the wiser.
        """

        path = Path(run_directory) / COHORT_MANIFEST_FILE
        if path.exists():
            existing = read_cohort_manifest(run_directory)
            if existing is not None and existing.cohort_id == self.cohort_id:
                return path
            raise ContractError(
                f"a cohort is already dispatched in {run_directory}: "
                "membership is fixed when the wave is submitted"
            )
        path.write_text(
            json.dumps(self.public_record(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return path


#: What the host can say about one proposed member of a wave.
COHORT_MEMBER_STATUSES = ("ready", "not_ready", "depends_on", "not_in_plan")


@dataclass(frozen=True)
class CohortMemberVerdictV1:
    """One proposed calculation, and what the host sees about it."""

    node_id: str
    status: str
    detail: str = ""

    def public_record(self) -> dict[str, Any]:
        return {
            "node_id": self.node_id,
            "status": self.status,
            "detail": self.detail,
        }


@dataclass(frozen=True)
class CohortValidityV1:
    """Whether a proposed wave can be dispatched, member by member.

    Typed evidence, never an exception. A wave the host cannot dispatch
    is something the Agent reads and selects again from; raising would
    teach a session to carry workarounds for a question the host owns.
    """

    members: tuple[str, ...]
    rows: tuple[CohortMemberVerdictV1, ...]
    summary: str

    @property
    def dispatchable(self) -> bool:
        return bool(self.members) and all(
            row.status == "ready" for row in self.rows
        )

    def public_record(self) -> dict[str, Any]:
        return {
            "dispatchable": self.dispatchable,
            "members": list(self.members),
            "rows": [row.public_record() for row in self.rows],
            "summary": self.summary,
        }


def validate_wave(
    *,
    proposed: tuple[str, ...],
    ready: tuple[str, ...],
    edges: tuple[tuple[str, str], ...],
    planned: tuple[str, ...] = (),
) -> CohortValidityV1:
    """Judge a proposed wave against the host's own ready frontier.

    Readiness is a host/DAG fact with a single authority, and this does
    not compute it -- it is handed the frontier that authority produced.
    Which ready calculations to run together is the Agent's scientific
    strategy, so the host says only what it sees: ready, not yet ready
    and what it waits on, or ordered against a sibling in the same wave.

    Two members joined by an edge are one experiment rather than two,
    however ready they both are: running them together would mean the
    consumer starting before the Agent had seen the producer's evidence,
    which is the barrier this exists to keep.

    Nothing raises. Every shape the model could get wrong -- an empty
    wave, a repeated member, a node outside the plan -- comes back as a
    verdict.
    """

    members: list[str] = []
    repeated: list[str] = []
    for node_id in proposed:
        name = str(node_id).strip()
        if not name:
            continue
        if name in members:
            if name not in repeated:
                repeated.append(name)
            continue
        members.append(name)
    member_set = set(members)
    known = {str(item) for item in planned}
    ready_set = {str(item) for item in ready}
    producers: dict[str, list[str]] = {}
    for source, target in edges:
        producers.setdefault(str(target), []).append(str(source))

    rows = []
    for node_id in members:
        blocking = sorted(
            producer
            for producer in producers.get(node_id, ())
            if producer in member_set
        )
        if blocking:
            rows.append(
                CohortMemberVerdictV1(
                    node_id=node_id,
                    status="depends_on",
                    detail=(
                        "consumes "
                        + ", ".join(blocking)
                        + " in this same wave, so the two are one "
                        "experiment: run the producer first and choose "
                        "the consumer after reading its evidence"
                    ),
                )
            )
            continue
        if node_id in ready_set:
            rows.append(CohortMemberVerdictV1(node_id=node_id, status="ready"))
            continue
        if known and node_id not in known:
            # A mistyped id read as `not_ready` with "is not in this
            # plan's ready frontier" -- word for word what a node waiting
            # on a producer says -- so the Agent was invited to wait
            # forever for a calculation that does not exist.
            rows.append(
                CohortMemberVerdictV1(
                    node_id=node_id,
                    status="not_in_plan",
                    detail=(
                        "is not a calculation in this workflow; check the "
                        "node id against the plan you made"
                    ),
                )
            )
            continue
        waiting = sorted(producers.get(node_id, ()))
        rows.append(
            CohortMemberVerdictV1(
                node_id=node_id,
                status="not_ready",
                detail=(
                    "waits on " + ", ".join(waiting)
                    if waiting
                    else "is not in this plan's ready frontier"
                ),
            )
        )

    note = (
        " (" + ", ".join(repeated) + " was named more than once and "
        "counts once)"
        if repeated
        else ""
    )
    if not members:
        summary = "the proposed wave is empty"
    elif all(row.status == "ready" for row in rows):
        summary = (
            f"{len(members)} calculations ready to run together" + note
        )
    else:
        summary = (
            "; ".join(
                f"{row.node_id}: {row.status}"
                for row in rows
                if row.status != "ready"
            )
            + note
        )
    return CohortValidityV1(
        members=tuple(members), rows=tuple(rows), summary=summary
    )


def authorise_cohort_element(
    run_directory: str | Path,
    *,
    element: int | None,
    bundle_sha256: str,
) -> str | None:
    """The calculation this element is authorised to run, or None.

    The approval bundle is one-shot, which is right: a human approved one
    execution of one plan. "One-shot" was implemented as "one process",
    and a wave is N processes running N members of that same approved
    plan -- so the second element was either admitted through the
    continuation path, whose own contract says it authorises nothing, or
    refused outright as "a second independent execution of a consumed
    bundle".

    Membership replaces the second claim. The manifest already names
    exactly which calculations this approval covers and which element
    runs which, it is digest-bound, and it was written before any element
    started. An element is authorised by being in it. The approval is
    still consumed once, by the dispatcher.

    Authority is checked against *this* approval: a manifest written for
    another bundle authorises nothing here, so a stale run directory
    cannot lend its membership to a different decision.

    Returns:
        str | None: The node id, or None when no cohort governs this run.

    Raises:
        ContractError: If the element is outside the cohort, or the
            cohort belongs to another approval.
    """

    manifest = read_cohort_manifest(run_directory)
    if manifest is None or element is None:
        return None
    if manifest.bundle_sha256 != str(bundle_sha256):
        raise ContractError(
            "this cohort was dispatched for another approval "
            f"({manifest.bundle_sha256[:12]}...), so it authorises nothing "
            "for this one"
        )
    return manifest.node_for_element(element)


def execution_result_file(
    run_directory: str | Path, *, element: int | None = None
) -> Path:
    """Where one element writes what it did.

    A cohort's elements each redirect to their own path. They used to
    share one: the last writer defined the cycle, every element reported
    `partial` by construction because none walks every approved node, and
    a shorter record written after a longer one left the first one's tail
    behind it.

    ``element=None`` is the single-job path, unchanged.
    """

    base = Path(run_directory)
    if element is None:
        return base / "execution-result.json"
    return base / f"execution-result.{int(element)}.json"


def cohort_completion(
    run_directory: str | Path,
) -> tuple[bool | None, tuple[str, ...]]:
    """Whether every member of this wave has reached a terminal state.

    Asked of the durable stream, never of a file's existence. The job
    script creates its redirect target before the engine starts, so
    ``is_file()`` was true from second zero -- for three hours of a
    running calculation, and forever after a job killed before it wrote
    anything.

    Terminality, not success: a member that failed or was cancelled has
    ended, and its outcome is evidence the Agent must see. A member still
    inside its launch lease has not ended, whatever else the stream says.

    Returns:
        tuple: ``(complete, pending)``. ``complete`` is ``None`` when no
        cohort was dispatched, which is the single-job path and not a
        failure.
    """

    from types import SimpleNamespace

    from chemsmart.agent.terminal_states import run_live_leases

    manifest = read_cohort_manifest(run_directory)
    if manifest is None:
        return None, ()
    # Read the stream as lines rather than as reconstructed events: this
    # question needs only the kind and the payload, and a stream one
    # element is still appending to must not fail the barrier.
    try:
        lines = (
            (Path(run_directory) / "events.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
        )
    except OSError:
        return False, manifest.node_ids
    events = []
    for line in lines:
        try:
            raw = json.loads(line)
        except json.JSONDecodeError:
            continue
        events.append(
            SimpleNamespace(
                kind=str(raw.get("kind") or ""),
                payload=raw.get("payload") or {},
            )
        )

    # Terminality, not an execution receipt. A member cancelled before
    # launch, or refused admission, reaches a terminal state without ever
    # running an engine; asking for a receipt made the wave wait on it
    # forever, so one cancelled calculation could never wake the Agent.
    from chemsmart.agent.execution import TERMINAL_NODE_RUN_STATES

    finished: set[str] = set()
    for event in events:
        kind = getattr(event, "kind", "")
        payload = getattr(event, "payload", None) or {}
        record = payload.get("record") or {}
        node_id = str(
            payload.get("node_id") or (record or {}).get("node_id") or ""
        )
        if not node_id:
            continue
        if kind == "program_execution_observed":
            finished.add(node_id)
        elif kind == "workflow_node_state_changed":
            # `node_state` is this node's word. `record` is the *run's*
            # state record, and `record["state"]` is the run's summary --
            # ambiguous, running, validated, failed, blocked, cancelled.
            # Five of those six are also node-terminal words, so reading
            # it was a type confusion that type-checks, and it failed in
            # both directions: a member cancelled while a sibling still
            # ran was written under the run's "running" and never counted
            # finished (the wave then waits forever, which is the case
            # this function's own docstring says it fixed), while once any
            # member ended "ambiguous" the run's word outranked the rest
            # and the next member's ordinary pending -> running
            # transition marked it finished with its engine starting.
            state = str(payload.get("node_state") or "")
            if state in TERMINAL_NODE_RUN_STATES:
                finished.add(node_id)
    live = set(run_live_leases(events))
    pending = tuple(
        node_id
        for node_id in manifest.node_ids
        if node_id not in finished or node_id in live
    )
    return (not pending), pending


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
        text = path.read_text(encoding="utf-8")
    except OSError:
        # No file at all: a single-job dispatch, which is absence.
        return None
    # A file that exists and cannot be read is **not** absence. `None`
    # means "no cohort", and `cohort_frontier(ready, None)` admits every
    # ready node -- so a truncated manifest would silently turn a bounded
    # wave back into the flowing walk it exists to prevent, executing
    # work the Agent did not ask for in this wave.
    try:
        raw = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ContractError(
            f"the cohort manifest at {path} is unreadable, and a wave "
            "whose membership cannot be read is not a wave without one"
        ) from exc
    if not isinstance(raw, dict):
        raise ContractError(
            f"the cohort manifest at {path} is damaged: a manifest is an "
            "object naming this wave's members"
        )
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
    "authorise_cohort_element",
    "build_cohort_manifest",
    "cohort_completion",
    "cohort_frontier",
    "execution_result_file",
    "CohortMemberVerdictV1",
    "CohortValidityV1",
    "COHORT_MEMBER_STATUSES",
    "validate_wave",
    "read_cohort_manifest",
]
