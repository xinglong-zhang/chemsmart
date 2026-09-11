"""The goal: one human decision that a bounded recovery loop consumes.

A frozen graph per approval means an agent cannot recover through
failure: the correct next move after a typed terminal state -- read it,
revise the route, continue -- was structurally unreachable, because the
approval was spent and every revision was a returning human action. The
goal moves the grain: the human approves the observables, the identity
bindings, the physical conditions, the envelope with its engine-call,
wall-clock, and revision budgets, and the complete initial plan; the
host then admits a revised workflow only when it cites the typed
terminal evidence it answers and preserves every invariant the human
approved. The model never approves. The goal approval is the sole
authority a revision consumes, and every admission names it beside the
human who granted it.

What lives here is deliberately deterministic: records, a ledger, the
condition extractor, and the admission checks. No provider call, no
retry policy, no judgement of whether a revision is scientifically
wise -- grading a route is what execution and validation are for.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

import yaml

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
    require_identifier,
)

GOAL_SCHEMA_VERSION = "chemsmart.goal.v1"

#: Terminal goal settlements. ``unreachable_from_evidence`` is the typed
#: refusal the loop exists to make as deliverable as an answer: settling
#: there requires receipts, never prose alone.
GOAL_SETTLEMENTS = (
    "achieved",
    "achieved_with_observations",
    "exhausted",
    "unreachable_from_evidence",
    "returned_to_human",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def extract_plan_conditions(
    *,
    project_settings_texts: tuple[str, ...],
    thermochemistry_controls: tuple[
        tuple[float | None, float | None, float | None], ...
    ],
) -> dict[str, Any]:
    """The physical conditions a plan runs under, as comparable data.

    Solvent names are conditions -- water versus toluene changes the
    chemistry -- while the continuum model that implements them is
    method, free to move under ruling 2. So the extractor collects
    every ``solvent`` value from the effective project settings and the
    (temperature, pressure, concentration) triples from planned
    thermochemistry, and nothing else. A revision must reproduce the
    sets exactly: adding a solvent, dropping solvation, or moving a
    temperature is a new goal and a returning human decision.
    """

    solvents: set[str] = set()
    for text in project_settings_texts:
        try:
            loaded = yaml.safe_load(text) or {}
        except yaml.YAMLError:
            continue

        def _walk(value: Any) -> None:
            if isinstance(value, Mapping):
                for key, item in value.items():
                    if str(key).strip().lower() == "solvent" and isinstance(
                        item, str
                    ):
                        cleaned = item.strip().lower()
                        if cleaned:
                            solvents.add(cleaned)
                    else:
                        _walk(item)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    _walk(item)

        _walk(loaded)
    thermo = sorted(
        (
            (
                None if t is None else round(float(t), 6),
                None if p is None else round(float(p), 6),
                None if c is None else round(float(c), 6),
            )
            for t, p, c in thermochemistry_controls
        ),
    )
    return {
        "solvents": tuple(sorted(solvents)),
        "thermochemistry": tuple(thermo),
    }


def conditions_from_review(review: Mapping[str, Any]) -> dict[str, Any]:
    """Extract the conditions a stored execution review displays."""

    node_reviews = review.get("node_reviews") or ()
    texts = tuple(
        str(item.get("project_settings_text") or "")
        for item in node_reviews
        if isinstance(item, Mapping)
    )
    toolchain = review.get("scientific_toolchain_plan") or {}
    controls: list[tuple[float | None, float | None, float | None]] = []
    for node in toolchain.get("analysis_nodes") or ():
        if not isinstance(node, Mapping):
            continue
        if str(node.get("analysis_kind") or "") != "thermochemistry":
            continue
        controls.append(
            (
                node.get("temperature_k"),
                node.get("pressure_atm"),
                node.get("concentration_mol_l"),
            )
        )
    return extract_plan_conditions(
        project_settings_texts=texts,
        thermochemistry_controls=tuple(controls),
    )


@dataclass(frozen=True)
class GoalRecordV1:
    """What the one human decision covers, digest-bound."""

    schema_version: str
    goal_id: str
    task_spec_sha256: str
    scientific_identity_sha256: str
    conditions: Mapping[str, Any]
    envelope: Mapping[str, Any]
    max_revisions: int
    granted_by: str
    initial_review_sha256: str
    created_at: str
    goal_sha256: str = ""

    def __post_init__(self) -> None:
        if self.schema_version != GOAL_SCHEMA_VERSION:
            raise ContractError("unsupported goal schema")
        require_identifier(self.goal_id, "goal_id")
        if not self.granted_by.strip():
            raise ContractError("a goal names the human who granted it")
        if (
            isinstance(self.max_revisions, bool)
            or not isinstance(self.max_revisions, int)
            or self.max_revisions < 0
        ):
            raise ContractError("max_revisions must be a non-negative integer")
        body = canonical_data(
            {
                "schema_version": self.schema_version,
                "goal_id": self.goal_id,
                "task_spec_sha256": self.task_spec_sha256,
                "scientific_identity_sha256": (
                    self.scientific_identity_sha256
                ),
                "conditions": canonical_data(self.conditions),
                "envelope": canonical_data(self.envelope),
                "max_revisions": self.max_revisions,
                "granted_by": self.granted_by,
                "initial_review_sha256": self.initial_review_sha256,
                "created_at": self.created_at,
            }
        )
        digest = canonical_sha256(body)
        if self.goal_sha256:
            if self.goal_sha256 != digest:
                raise ContractError("goal record digest mismatch")
        else:
            object.__setattr__(self, "goal_sha256", digest)

    @property
    def actor(self) -> str:
        """The composite authority a host-admitted revision records."""

        return f"goal-approval:{self.goal_id}"


@dataclass(frozen=True)
class GoalBudgetsV1:
    """What remains of the human's grant, from durable evidence only."""

    engine_calls_remaining: int
    wall_seconds_remaining: float
    revisions_remaining: int
    #: The excursion grant's own line; it never lends to, or borrows
    #: from, the engine-call budget.
    excursion_calls_remaining: int = 0


class GoalLedger:
    """Append-only account of one goal's cycles under its grant."""

    def __init__(self, directory: str | Path) -> None:
        self.directory = Path(directory)
        self.goal_path = self.directory / "goal.json"
        self.ledger_path = self.directory / "ledger.jsonl"

    def create(self, record: GoalRecordV1) -> None:
        self.directory.mkdir(parents=True, exist_ok=True)
        if self.goal_path.exists():
            raise ContractError("goal already exists; a goal is created once")
        body = canonical_data(
            {
                "schema_version": record.schema_version,
                "goal_id": record.goal_id,
                "task_spec_sha256": record.task_spec_sha256,
                "scientific_identity_sha256": (
                    record.scientific_identity_sha256
                ),
                "conditions": canonical_data(record.conditions),
                "envelope": canonical_data(record.envelope),
                "max_revisions": record.max_revisions,
                "granted_by": record.granted_by,
                "initial_review_sha256": record.initial_review_sha256,
                "created_at": record.created_at,
                "goal_sha256": record.goal_sha256,
            }
        )
        self.goal_path.write_text(
            json.dumps(body, indent=1, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        self.append("goal_created", {"goal_sha256": record.goal_sha256})

    def load(self) -> GoalRecordV1:
        raw = json.loads(self.goal_path.read_text(encoding="utf-8"))
        return GoalRecordV1(**raw)

    def append(self, kind: str, payload: Mapping[str, Any]) -> None:
        entry = {
            "kind": str(kind),
            "at": _utc_now(),
            "payload": canonical_data(payload),
        }
        with self.ledger_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(entry, sort_keys=True) + "\n")

    def entries(self) -> tuple[dict[str, Any], ...]:
        if not self.ledger_path.exists():
            return ()
        rows = []
        for line in self.ledger_path.read_text(encoding="utf-8").splitlines():
            text = line.strip()
            if text:
                rows.append(json.loads(text))
        return tuple(rows)

    def budgets(self, record: GoalRecordV1) -> GoalBudgetsV1:
        """Remaining grant, decremented by the ledger's recorded runs.

        The wall clock is durable here by construction: each recorded
        run carries the engine seconds its receipts state, so a goal
        resumed in a new process still knows what it has spent -- the
        per-invocation monotonic clock never was durable and is not
        consulted.
        """

        envelope = dict(record.envelope)
        calls = int(envelope.get("max_engine_calls") or 0)
        excursions = int(envelope.get("max_excursion_calls") or 0)
        wall = float(envelope.get("episode_wall_time_seconds") or 0.0)
        revisions = record.max_revisions
        for entry in self.entries():
            if entry["kind"] == "run_recorded":
                calls -= int(
                    entry["payload"].get("engine_calls_consumed") or 0
                )
                excursions -= int(
                    entry["payload"].get("excursion_calls_consumed") or 0
                )
                wall -= float(
                    entry["payload"].get("engine_wall_seconds") or 0.0
                )
            elif entry["kind"] in {"revision_admitted", "rewake_opened"}:
                # A re-wake after a cycle that delivered nothing is
                # charged as the revision it consumes (owner ruling,
                # 2026-09-05).
                revisions -= 1
        return GoalBudgetsV1(
            engine_calls_remaining=max(calls, 0),
            wall_seconds_remaining=max(wall, 0.0),
            revisions_remaining=max(revisions, 0),
            excursion_calls_remaining=max(excursions, 0),
        )

    def settle(
        self,
        state: str,
        *,
        reasons: tuple[str, ...],
        evidence: Mapping[str, Any] | None = None,
    ) -> None:
        if state not in GOAL_SETTLEMENTS:
            raise ContractError(f"unsupported goal settlement: {state!r}")
        if state in {
            "unreachable_from_evidence",
            "achieved_with_observations",
        } and not (evidence or {}):
            raise ContractError(
                f"{state} settles on receipts, never prose alone"
            )
        self.append(
            "goal_settled",
            {
                "state": state,
                "reasons": tuple(reasons),
                "evidence": canonical_data(evidence or {}),
            },
        )


def session_read_run_outcome(
    session_events_path: str | Path, run: str
) -> bool:
    """Whether this session pulled the typed outcome of one named run.

    The evidence gate does not ask the model to copy digests: when a
    session reads a run's outcome, the host itself records a
    ``run_outcome_inspected`` event naming the run and the exact
    stream bytes served. The first live goal round showed why the
    weaker check -- "the inspection tool succeeded at least once" --
    is unsound: a bare listing over an empty root succeeded, and the
    gate credited a session that had read nothing. What this proves
    is that the named run's typed outcome entered the session's
    context; comprehension is graded by the human reading, as always.
    """

    reference = str(run or "").strip()
    if not reference:
        return False
    try:
        lines = (
            Path(session_events_path).read_text(encoding="utf-8").splitlines()
        )
    except OSError:
        return False
    for line in lines:
        text = line.strip()
        if not text:
            continue
        try:
            event = json.loads(text)
        except json.JSONDecodeError:
            continue
        if str(event.get("kind") or "") != "run_outcome_inspected":
            continue
        payload = event.get("payload") or {}
        if str(payload.get("run") or "") == reference:
            return True
    return False


@dataclass(frozen=True)
class RevisionAdmissionV1:
    """One deterministic admission verdict, with every check named."""

    admitted: bool
    checks: Mapping[str, bool]
    reasons: tuple[str, ...] = ()
    cited_evidence_event_hashes: tuple[str, ...] = ()
    #: The scope this admission established, when the goal had none to
    #: preserve. Empty on every later revision, which compares instead.
    bound_scientific_identity_sha256: str = ""
    bound_conditions: Mapping[str, Any] | None = None


def goal_scope_is_unbound(goal: GoalRecordV1) -> bool:
    """True when no executable review has ever fixed this goal's scope.

    A goal created from an analysis-only first cycle carries
    ``initial_review_sha256 == ""`` and, with it, an empty identity and
    empty conditions -- not because the human approved "no molecule and
    no solvent", but because no plan with a molecule in it had been
    displayed yet. Comparing a later executable plan against those empty
    values refuses it for changing a scope that was never set, so a goal
    whose first cycle read registered results could never launch its
    first calculation (OPEN-2 ino3-qwen would have hit this the moment
    it was woken; found in review 2026-09-09, before it cost a window).

    The emptiness of ``conditions`` cannot carry this by itself, because
    a genuine gas-phase plan is empty too. The absence of an initial
    review digest can: it says no executable partition was ever
    displayed under this goal.
    """

    return not str(goal.initial_review_sha256 or "").strip()


def admit_revision(
    *,
    goal: GoalRecordV1,
    budgets: GoalBudgetsV1,
    revision_review: Mapping[str, Any],
    revision_scientific_identity_sha256: str,
    session_events_path: str | Path,
    prior_outcome_evidence_hashes: tuple[str, ...],
    previous_run_reference: str,
    wake_embedded_run: str = "",
    bound_scientific_identity_sha256: str | None = None,
    bound_conditions: Mapping[str, Any] | None = None,
) -> RevisionAdmissionV1:
    """Admit or return one revision, deterministically.

    Every check is a comparison the human's one decision already
    covered. A failed check returns the revision to the human with the
    named reason -- it never silently narrows, never retries, and never
    grades the science: whether the revised route is *wise* is what
    execution, validation, and the reading scientist are for.

    ``bound_*`` carry the scope a previous cycle established when the
    goal record itself had none; the driver reads them from the ledger.
    Absent them, an unbound goal's first executable revision *sets* the
    identity and conditions instead of comparing against nothing, and
    every revision after it compares against what was set.
    """

    checks: dict[str, bool] = {}
    reasons: list[str] = []

    effective_identity = (
        goal.scientific_identity_sha256
        if bound_scientific_identity_sha256 is None
        else str(bound_scientific_identity_sha256)
    )
    effective_conditions = (
        goal.conditions if bound_conditions is None else bound_conditions
    )
    binds_scope = (
        goal_scope_is_unbound(goal)
        and not str(effective_identity or "").strip()
    )

    revision_conditions = canonical_data(
        conditions_from_review(revision_review)
    )

    if binds_scope:
        # Nothing to preserve: this revision is the first executable
        # partition this goal has ever displayed, so it fixes the scope
        # that every later revision is held to.
        checks["identity_preserved"] = True
        checks["conditions_preserved"] = True
        checks["scope_bound_here"] = True
    else:
        checks["identity_preserved"] = (
            revision_scientific_identity_sha256 == effective_identity
        )
        if not checks["identity_preserved"]:
            reasons.append(
                "the revision binds different molecular identities or "
                "electronic states than the goal approved"
            )
        checks["conditions_preserved"] = revision_conditions == canonical_data(
            effective_conditions
        )
        if not checks["conditions_preserved"]:
            reasons.append(
                "the revision changes physical conditions (solvent or "
                "thermochemical state) the goal approved"
            )

    goal_envelope = dict(goal.envelope)
    review_envelope = dict(revision_review.get("execution_envelope") or {})
    goal_pairs = {
        str(program): tuple(engines)
        for program, engines in (
            goal_envelope.get("allowed_program_engines") or ()
        )
    }
    review_pairs = {
        str(program): tuple(engines)
        for program, engines in (
            review_envelope.get("allowed_program_engines") or ()
        )
    }
    checks["programs_within_envelope"] = all(
        program in goal_pairs and set(engines) <= set(goal_pairs[program])
        for program, engines in review_pairs.items()
    )
    if not checks["programs_within_envelope"]:
        reasons.append(
            "the revision names a program or engine outside the goal's "
            "envelope"
        )

    # An excursion node is charged to its own line: the plan record
    # names which nodes carry the tag, so the two counts never mix.
    excursion_ids = {
        str(node.get("node_id") or "")
        for node in (
            (revision_review.get("scientific_plan") or {}).get("nodes") or ()
        )
        if isinstance(node, Mapping) and node.get("excursion")
    }
    reviewed = tuple(
        item
        for item in (revision_review.get("node_reviews") or ())
        if isinstance(item, Mapping)
    )
    executable = tuple(
        item
        for item in reviewed
        if str(item.get("node_id") or "") not in excursion_ids
    )
    excursions = tuple(
        item
        for item in reviewed
        if str(item.get("node_id") or "") in excursion_ids
    )
    checks["engine_budget_remains"] = (
        len(executable) <= budgets.engine_calls_remaining
    )
    if not checks["engine_budget_remains"]:
        reasons.append(
            f"the revision plans {len(executable)} engine calls with "
            f"{budgets.engine_calls_remaining} remaining in the grant"
        )
    checks["excursion_budget_remains"] = (
        len(excursions) <= budgets.excursion_calls_remaining
    )
    if not checks["excursion_budget_remains"]:
        reasons.append(
            f"the revision plans {len(excursions)} excursion calls with "
            f"{budgets.excursion_calls_remaining} remaining in the grant"
        )

    checks["revision_budget_remains"] = budgets.revisions_remaining > 0
    if not checks["revision_budget_remains"]:
        reasons.append("the goal's revision budget is exhausted")

    # The revision must have been planned by a session that provably
    # HELD the typed outcome of the run it revises. Two host-attested
    # routes prove it: the wake context embedded that run's outcome
    # (recorded in the goal ledger when the host composed it), or the
    # session read it through inspect_run_outcome, which records a
    # run-bound event. Neither route grades comprehension.
    revised_run = str(previous_run_reference or "").strip()
    checks["evidence_read"] = bool(revised_run) and (
        str(wake_embedded_run or "").strip() == revised_run
        or session_read_run_outcome(session_events_path, revised_run)
    )
    if not checks["evidence_read"]:
        reasons.append(
            "the revising session never held the typed outcome of the "
            f"run it revises ({revised_run or 'unrecorded'}): the wake "
            "context did not embed it and no run-bound read appears in "
            "the session's stream; a plan that answers a failure it "
            "never read is answering a guess"
        )

    admitted = all(checks.values())
    return RevisionAdmissionV1(
        admitted=admitted,
        checks=checks,
        reasons=tuple(reasons),
        cited_evidence_event_hashes=(
            tuple(prior_outcome_evidence_hashes) if admitted else ()
        ),
        bound_scientific_identity_sha256=(
            str(revision_scientific_identity_sha256)
            if admitted and binds_scope
            else ""
        ),
        bound_conditions=(
            conditions_from_review(revision_review)
            if admitted and binds_scope
            else None
        ),
    )


__all__ = [
    "goal_scope_is_unbound",
    "GOAL_SCHEMA_VERSION",
    "GOAL_SETTLEMENTS",
    "GoalBudgetsV1",
    "GoalLedger",
    "GoalRecordV1",
    "RevisionAdmissionV1",
    "admit_revision",
    "conditions_from_review",
    "extract_plan_conditions",
    "session_read_run_outcome",
]
