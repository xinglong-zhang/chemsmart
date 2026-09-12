"""The goal driver: plan, one human decision, execute, read, revise.

Every entry point that drives a goal -- ``chemsmart agent goal``, the
terminal interface, a scheduler job re-entering a parked run -- is a
view of this one step machine. Each phase boundary is a durable ledger
entry, so a process may stop after any step and a later process may
resume from the ledger: that is what lets an approved partition be
submitted to a scheduler, the driver exit, and the job's own tail wake
the goal when the engines are done.

The driver cycles the existing authority chain -- ``run_live_agent_session``
for planning, the stored-review resolution for the decision, the
provider-free executor for engines -- and adds nothing to any of them.
What is new is only the connective tissue the chain lacked: a durable
goal ledger, a typed wake context built from the same terminal
derivation a session's own tool reads, and the deterministic revision
admission from :mod:`chemsmart.agent.goal`. Detach, typed completion,
re-invoke: the session never blocks on an engine, the executor never
sees a provider, and every wake reads sealed evidence rather than a
retelling.

Dependencies are injected so the loop's bookkeeping is testable without
a provider or an engine; production callers pass nothing and get the
real chain.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Mapping, Sequence

from chemsmart.agent._contracts import ContractError, canonical_data
from chemsmart.agent.delivery import (
    current_assessments,
    unresolved_requirement_ids,
)
from chemsmart.agent.execution import anomaly_standing
from chemsmart.agent.goal import (
    GOAL_SCHEMA_VERSION,
    GoalLedger,
    GoalRecordV1,
    admit_revision,
    conditions_from_review,
    goal_scope_is_unbound,
)
from chemsmart.agent.rules import rules_by_id
from chemsmart.agent.terminal_states import (
    PRIOR_ANOMALIES_FILE,
    REPAIRABLE_NODE_STATES,
    is_provider_transport_terminal,
)
from chemsmart.agent.workspace_record import (
    printed_modes,
    read_workspace_record,
    record_run,
    render_workspace_record,
    uncharacterised_artifacts,
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class GoalLoopResultV1:
    """How one goal ended, for the caller and the report."""

    goal_id: str
    settlement: str
    cycles: int
    revisions_admitted: int
    reasons: tuple[str, ...] = ()


def _default_plan_session(**kwargs: Any) -> Any:
    from chemsmart.agent.live_session import run_live_agent_session

    return run_live_agent_session(**kwargs)


def _default_resolve(
    *,
    review_file: Path,
    workspace: Path,
    decision: str,
    actor: str,
    approval_id: str,
) -> tuple[str, Path]:
    from chemsmart.agent.live_session import (
        inspect_workflow_execution_replay,
        resolve_workflow_execution_review,
    )

    report = inspect_workflow_execution_replay(
        review_file=review_file,
        workspace=workspace,
        task_spec_sha256="",
    )
    scope = (
        Path(workspace).resolve()
        / ".chemsmart-agent"
        / "replays"
        / approval_id
    )
    scope.mkdir(parents=True, exist_ok=True, mode=0o700)
    resolve_workflow_execution_review(
        review_file=review_file,
        reviewed_sha256=report["review_sha256"],
        decision=decision,
        actor=actor,
        output_file=scope / "bundle.json" if decision == "approve" else None,
        decision_log=scope / "decisions.jsonl",
        approval_id=approval_id,
    )
    return str(report["review_sha256"]), scope / "bundle.json"


def _default_execute(
    *,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    stop_file: Path | None = None,
) -> Any:
    from chemsmart.agent.executor import execute_approved_workflow

    # The charter gives the human a withdrawal at any node boundary, and
    # this hook dropped the stop file, so the executor's should_stop was
    # None and only the cycle transition read it. Twice observed: REACH-1
    # po3 ran 2.5 h of nodes after its STOP, and OPEN-1 po3 launched a
    # transition-state search fourteen minutes after mine (2026-09-07).
    return execute_approved_workflow(
        approval_file=approval_file,
        workspace=workspace,
        run_directory=run_directory,
        task_spec_sha256="",
        stop_file=stop_file,
    )


def _execute_hook_takes_stop_file(hook: Any) -> bool:
    """Whether this execute hook accepts the human's withdrawal.

    The hook is replaceable (tests and the terminal interface pass their
    own), so the withdrawal is offered rather than forced: a hook that
    does not take it is called exactly as before.
    """

    import inspect

    try:
        parameters = inspect.signature(hook).parameters
    except (TypeError, ValueError):
        return False
    return "stop_file" in parameters or any(
        item.kind is inspect.Parameter.VAR_KEYWORD
        for item in parameters.values()
    )


def _anomaly_evidence(
    evidence: Mapping[str, Any] | None,
    ledger_anomalies: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """The receipts a settlement with observations stands on.

    The anomaly receipts are the observation's own evidence: decisions
    live in a session's stream and the executor's run stream carries
    none, so a word that requires receipts must bring its own.
    """

    digests = {
        str(item.get("receipt_sha256") or "")
        for item in ledger_anomalies
        if str(item.get("receipt_sha256") or "")
    }
    merged = dict(evidence or {})
    receipts = set(merged.get("receipt_sha256s") or ()) | digests
    if receipts:
        merged["receipt_sha256s"] = tuple(sorted(receipts))
    return merged


def _achieved_word(
    delivery: "_AnalysisDelivery",
    ledger_anomalies: Sequence[Mapping[str, Any]] = (),
) -> tuple[str, tuple[str, ...]]:
    """The settlement word for a certified delivery.

    Plain ``achieved`` means the host saw nothing it could not explain;
    a delivery that carries host-recorded observations settles with the
    word that says so, because the one word a human reads first must
    not hide what the run found (owner ruling, 2026-09-03).
    """

    observed = tuple(
        sorted(
            set(delivery.anomaly_output_ids)
            | set(_anomaly_output_ids_from_records(ledger_anomalies))
        )
    )
    # Where a delivered number came from is true of either word, so it is
    # said under both: the anomaly-carrying word must not be the only one
    # that discloses provenance.
    checked = set(delivery.characterised_source_quantity_ids)
    provenance: tuple[str, ...] = ()
    if delivery.characterised_source_quantity_ids:
        provenance = provenance + (
            "delivered from a characterised result: "
            + ", ".join(delivery.characterised_source_quantity_ids),
        )
    unchecked = tuple(
        item
        for item in delivery.failed_source_quantity_ids
        if item not in checked
    )
    if unchecked:
        # The number may be exactly the finding; what a reader needs is
        # that it came from a run that did not meet the promise it was
        # launched under, and nobody had the host check the structure.
        provenance = provenance + (
            "delivered from a node that did not meet its promise, "
            "uncharacterised: " + ", ".join(unchecked),
        )
    if delivery.uncharacterised_source_quantity_ids:
        # The structure was never asked what it is: no frequencies were
        # printed, so the promise of a minimum or a saddle was never
        # checked, and the number is delivered saying so.
        provenance = provenance + (
            "delivered from a result whose stationary point is "
            "uncharacterised (no frequencies printed): "
            + ", ".join(delivery.uncharacterised_source_quantity_ids),
        )
    if delivery.delivered_in_earlier_cycles:
        provenance = provenance + (
            "delivered in an earlier cycle: "
            + ", ".join(delivery.delivered_in_earlier_cycles),
        )
    post_hoc = tuple(
        str(row.get("observable_id") or "")
        for row in delivery.prediction_rows
        if row.get("declared_after_evidence")
    )
    if post_hoc:
        # An expectation written once the numbers existed is a
        # restatement, not a pre-registration, and the settlement says
        # which rows those were (OPEN-1 ino3, 2026-09-07).
        provenance = provenance + (
            "declared after the evidence existed: " + ", ".join(post_hoc),
        )
    # A goal whose headline carried a precision requirement settles
    # without the settlement ever saying so: three reasons named the
    # anomaly, the earlier cycle and the decision's own uncertainties,
    # and not one named the tolerance, the number that answered it, or
    # what the host resolved to accept that number. No gate can decide
    # whether a stated magnitude is really the uncertainty of a claim --
    # that is chemistry. What the host owes instead is that its word is
    # auditable, so the word carries what it rests on (owner ruling,
    # 2026-09-09).
    # "the precision the task asked for was answered" claimed more than
    # the host checked, and could name the wrong requester. What the
    # host establishes is arithmetic on a number it resolved, against a
    # tolerance somebody declared; whether that number estimates the
    # relevant scientific error is the session's argument and stays
    # attributed. Two independent audits of the 2026-09-10 boundary
    # agreed on this wording, and the same run declared a tolerance of
    # its own for a decision the task set no precision for -- so the
    # sentence now says whose precision it was (SUFFICIENCY-5).
    discharged = tuple(
        (
            "the host resolved a stated uncertainty inside the tolerance "
            + {
                "task": "the task states",
                "session": "the session declared for its own decision "
                "question (the task fixes none here)",
            }.get(
                str(row.get("tolerance_origin") or ""),
                "on this contract (origin unstated)",
            )
            + ": "
        )
        + f"{row.get('observable_id')} within "
        f"{row.get('required_tolerance')} {row.get('unit')} at "
        f"{row.get('uncertainty')} {row.get('unit')}, "
        + (
            f"read by the host from {str(row.get('uncertainty_reference'))[:8]}"
            if row.get("uncertainty_reference")
            else f"on the session's own word ({row.get('uncertainty_basis')})"
        )
        # And what the host saw about that number without ruling on it.
        # A `met` standing on a spread of zero, on a single receipt, or
        # on a chain naming a coefficient the session supplied is
        # admitted now where it was refused before, so the word must
        # carry the observation to the human who reads it: the boundary
        # traded three refusals for one report, and the report is only
        # worth the trade where a reader meets it (owner ruling,
        # 2026-09-10).
        + (
            "; the host observed "
            + ", ".join(
                str(item) for item in row.get("uncertainty_observations") or ()
            )
            if row.get("uncertainty_observations")
            else ""
        )
        + (
            "; combined by "
            + str(
                (row.get("uncertainty_combination") or {}).get("rule")
                or "a rule the session did not state"
            )
            + (
                ", coverage "
                + str(
                    (row.get("uncertainty_combination") or {}).get("coverage")
                )
                if (row.get("uncertainty_combination") or {}).get("coverage")
                else ""
            )
            + (
                ", assuming "
                + str(
                    (row.get("uncertainty_combination") or {}).get(
                        "dependence"
                    )
                )
                if (row.get("uncertainty_combination") or {}).get("dependence")
                else ""
            )
        )
        + ". Adequacy of that estimate is the session's claim, not a "
        "host finding"
        for row in current_assessments(delivery.sufficiency).values()
        if str(row.get("state") or "") == "met"
    )
    provenance = provenance + discharged
    if delivery.decision_uncertainties:
        provenance = provenance + (
            "the recorded decision states its uncertainties: "
            + " | ".join(delivery.decision_uncertainties),
        )
    if delivery.flagged_quantity_ids:
        # Naming the anomaly is not naming the number. Five deliveries
        # across two windows reported a structure the host had flagged
        # as the answer, and the word said only that something somewhere
        # was observed (E4', 2026-09-03).
        provenance = provenance + (
            "delivered from the flagged result: "
            + ", ".join(delivery.flagged_quantity_ids),
        )
    if observed:
        return (
            "achieved_with_observations",
            (
                "the host completion gate certified the delivery; the "
                "host also recorded observations nobody asked for: "
                + ", ".join(observed),
            )
            + provenance,
        )
    return (
        "achieved",
        ("the host completion gate certified the delivery",) + provenance,
    )


def _open_requirement_reasons(delivery: Any) -> tuple[str, ...]:
    """What an open requirement stands at, for the human handed it back.

    The assessment reached the settlement only through the `met` branch
    of `_achieved_word`, so on returned_to_human, exhausted and
    unreachable_from_evidence it lived only in workspace-record.jsonl --
    and the live run that motivated this settled returned_to_human,
    which is how that leg went unexercised. A first attempt at this
    repair appended these rows to `_achieved_word`'s provenance, where
    a returned goal never reads them: the same computed-and-unconsumed
    pattern the repair exists to close, caught by its own witness
    (SUFFICIENCY-5, 2026-09-10).
    """

    rows = current_assessments(delivery.sufficiency).values()
    return tuple(
        f"{row.get('observable_id')} stands {row.get('state')} at "
        f"{row.get('uncertainty')} {row.get('unit')} against "
        f"{row.get('required_tolerance')} {row.get('unit')} "
        f"({row.get('tolerance_origin') or 'unstated'} tolerance)"
        + (
            "; the host observed "
            + ", ".join(
                str(item) for item in row.get("uncertainty_observations") or ()
            )
            if row.get("uncertainty_observations")
            else ""
        )
        for row in rows
        if str(row.get("state") or "") in {"short", "attested", "unstated"}
    )


def _settle_from_delivery(
    ledger: GoalLedger,
    *,
    goal_id: str,
    cycles: int,
    revisions_admitted: int,
    events_path: Path,
    terminal: str,
    workspace: Path | None = None,
) -> GoalLoopResultV1:
    """Settle a goal from one stream's typed delivery facts.

    Serves two shapes: a planning session that ended without an
    executable partition (its own stream carries the delivery), and an
    admitted revision whose approved bundle launched no engine -- the
    executor walks the analysis chain into the run directory's stream
    and no workflow run is recorded, which is a legitimate cycle
    ending, not a defect.
    """

    delivery = _analysis_delivery(
        events_path,
        flagged_artifact_sha256s=_flagged_artifact_sha256s(
            _goal_anomalies(ledger)
        ),
        uncharacterised_artifact_sha256s=(
            uncharacterised_artifacts(workspace) if workspace else ()
        ),
        goal_delivered_ids=_goal_delivered_ids(workspace, goal_id),
        declared_observables=_first_declarations(ledger),
    )
    evidence = _settlement_evidence(delivery)
    # What the goal declared and has not delivered under its id in any
    # cycle, read from the ledger's first declarations and the record --
    # a refusal made from the in-session route has no completion and so
    # no limitation ids to key on (NOVEL-3 po3, 2026-09-05).
    declared_ids = _required_declared_ids(ledger)
    delivered_ids = set(_goal_delivered_ids(workspace, goal_id))
    delivered_ids.update(delivery.delivered_quantity_ids)
    open_ids = tuple(
        dict.fromkeys(
            tuple(
                observable_id
                for observable_id in declared_ids
                if observable_id and observable_id not in delivered_ids
            )
            + tuple(delivery.undelivered_declared_ids)
        )
    )
    # An undelivered observable the session refused, and an observable
    # it delivered whose *precision* it refused, settle by one rule:
    # both are obligations the host verified as unreachable.
    refused_ids = tuple(
        dict.fromkeys(
            tuple(open_ids) + tuple(delivery.refused_requirement_ids)
        )
    )
    # The completion receipt is the host's certification; the session's
    # terminal word is its posture at exit. A live session delivered its
    # claims, passed the gate, then left an analysis-only draft and
    # ended "planned", and the goal returned to the human saying the
    # gate had not passed (E2 window, 2026-09-03). The receipt decides.
    certified = delivery.completion_status == "passed"
    if delivery.unresolved_requirement_ids:
        # The task said how good the answer had to be, the delivery says
        # it is not that good, and nothing computed, argued or refused
        # the difference. That is the human's to read, exactly as a
        # claim the session's own decision doubts is: the numbers stay
        # delivered and the word stops saying the goal was achieved.
        # No fifth settlement word (owner ruling 2026-09-09): a request
        # that is physically unreachable is a real ending, and the
        # session reaches it by refusing with receipts, which settles
        # unreachable_from_evidence -- a deliverable.
        settled = "returned_to_human"
        # "the precision the task asked for" is true of a restated
        # tolerance and false of one the session declared for a decision
        # question of its own, and both reach this branch. The host says
        # what it can check -- a requirement on this contract stands
        # unresolved -- and the assessment rows beneath name whose it
        # was (SUFFICIENCY-5, 2026-09-10).
        reasons = (
            "a delivered number does not answer a precision requirement "
            "on this contract, and the shortfall was neither computed "
            "away, re-claimed with the uncertainty it measured, nor "
            "refused with receipts: "
            + ", ".join(delivery.unresolved_requirement_ids),
        ) + _open_requirement_reasons(delivery)
    elif delivery.doubted_quantity_ids:
        # The session claimed from a receipt its own recorded decision
        # doubts. Whatever the completion word says -- partial by the
        # gate, or passed because the doubt came after -- a doubted
        # number is the human's to read.
        settled = "returned_to_human"
        reasons = (
            "a rendered claim stands on a receipt the session's "
            "own recorded decision doubts: "
            + ", ".join(delivery.doubted_quantity_ids),
        )
    elif (
        certified
        and delivery.blocked_output_ids
        and delivery.decisions
        and evidence
    ):
        # The plan's own completion receipt says a required output's
        # producer was declared blocked: the requested observable was
        # not delivered, and the recorded decision with its receipts is
        # the typed refusal.
        settled = "unreachable_from_evidence"
        reasons = (
            "the completion receipt names required outputs "
            "delivered without: " + ", ".join(delivery.blocked_output_ids),
        )
    elif (
        delivery.decisions
        and refused_ids
        and set(refused_ids) <= set(delivery.verified_unreachable_ids)
        and evidence
    ):
        # The typed refusal, reachable from either delivery route: the
        # session named each undelivered observable unreachable with the
        # producer it would need and receipts, and the host verified
        # each producer absent. Two sessions wrote this refusal in words
        # and neither reached the word (NOVEL-3 po3, ino3).
        # A refused *precision* arrives here too. Its observable was
        # claimed, so it is in no undelivered list, and the projections
        # simply removed it: a goal whose stated accuracy the host
        # certified as unreachable settled achieved, and the tool's own
        # reply had promised this word for it.
        settled = "unreachable_from_evidence"
        reasons = (
            "the recorded decision names these declared observables "
            "unreachable from the admissible evidence, and the host "
            "verified each: "
            + "; ".join(
                f"{observable_id} -- "
                f"{delivery.unreachable_bases.get(observable_id, '')}"
                for observable_id in refused_ids
            ),
        )
    elif (
        delivery.decisions
        and open_ids
        and set(open_ids)
        <= set(delivery.verified_unreachable_ids)
        | set(delivery.unverified_unreachable_ids)
    ):
        # Named unreachable, and the host could not verify at least one:
        # a human reads it, with the session's statement and the host's
        # basis side by side.
        settled = "returned_to_human"
        reasons = (
            "the recorded decision names declared observables unreachable "
            "that the host could not verify: "
            + "; ".join(
                f"{observable_id} -- "
                f"{delivery.unreachable_bases.get(observable_id, '')}"
                for observable_id in open_ids
                if observable_id in delivery.unverified_unreachable_ids
            ),
        )
    elif certified and delivery.undelivered_declared_ids:
        # The chain's kernels ran clean, and the headline the goal
        # declared was never claimed by its id in any cycle. That is not
        # a delivery a scientist would sign, whatever the completion word
        # says.
        settled = "returned_to_human"
        reasons = (
            "the completion certified the chain, but these declared "
            "observables have no claim carrying their id in any cycle: "
            + ", ".join(delivery.undelivered_declared_ids),
        ) + delivery.open_declared_misses
        if delivery.delivered_in_earlier_cycles:
            reasons = reasons + (
                "delivered in an earlier cycle: "
                + ", ".join(delivery.delivered_in_earlier_cycles),
            )
    elif certified and delivery.unanswered_verdicts:
        # The host itself rendered a verdict saying a delivered structure
        # is not what the task required, and no recorded decision cites
        # it. The required outputs are all present, so the completion
        # gate is green and the goal is not "achieved" in any sense a
        # scientist would sign: a human reads it.
        settled = "returned_to_human"
        reasons = (
            "a validation verdict failed and no recorded decision cites "
            "it: " + ", ".join(delivery.unanswered_verdicts),
        )
    elif certified and delivery.claims:
        goal_anomalies = _goal_anomalies(ledger)
        settled, reasons = _achieved_word(delivery, goal_anomalies)
        if settled == "achieved_with_observations":
            evidence = _anomaly_evidence(evidence, goal_anomalies)
        if terminal != "complete":
            reasons = reasons + (
                f"the session ended {terminal!r} after the completion "
                "passed; what it left unmaterialised is not part of the "
                "delivery",
            )
    elif delivery.stopped_by:
        # The stream says what stopped the run before a delivery could
        # exist -- a launch the executor refused, a review the host
        # refused to build. Quote it: a returned goal that names nothing
        # has not settled (two live goals, 2026-09-02).
        settled = "returned_to_human"
        reasons = delivery.stopped_by
    elif delivery.claims or delivery.decisions:
        # Something was recorded, but the host never certified
        # completion -- a human reads it, whatever the session's
        # terminal word was.
        settled = "returned_to_human"
        reasons = (
            f"the session ended {terminal!r} ({delivery.ending}); it "
            "recorded analysis but the host completion gate did not pass",
        )
    else:
        settled = "returned_to_human"
        reasons = (f"the session ended {terminal!r}: {delivery.ending}",)
    ledger.settle(settled, reasons=reasons, evidence=evidence)
    if workspace is not None and settled in {
        "achieved",
        "achieved_with_observations",
    }:
        _record_goal_qualification(
            ledger,
            workspace=workspace,
            goal_id=goal_id,
            current_run=f"goals/{goal_id}/runs/cycle-{cycles}",
            current_outcome=None,
        )
    return GoalLoopResultV1(
        goal_id=goal_id,
        settlement=settled,
        cycles=cycles,
        revisions_admitted=revisions_admitted,
        reasons=reasons,
    )


def _typed_error_settlement(
    ledger: GoalLedger,
    *,
    goal_id: str,
    cycles: int,
    revisions_admitted: int,
    stage: str,
    error: ContractError,
) -> GoalLoopResultV1:
    """Settle a typed error instead of letting it escape unsettled.

    Observed live (C5): a session attempted terminal completion against
    a red gate, the event store's ContractError propagated through the
    loop, the process died, and the goal's durable story ended at
    wake_composed with no settlement -- violating "a goal settles into
    one typed state". A typed contract error is an outcome the human
    reads, not a crash; a non-contract exception stays a crash, because
    a genuine defect must not be laundered into a settlement.
    """

    reason = f"cycle {cycles}, {stage}: {error}"
    ledger.settle("returned_to_human", reasons=(reason,))
    return GoalLoopResultV1(
        goal_id=goal_id,
        settlement="returned_to_human",
        cycles=cycles,
        revisions_admitted=revisions_admitted,
        reasons=(reason,),
    )


def _review_record(review_file: Path) -> Mapping[str, Any]:
    return json.loads(Path(review_file).read_text(encoding="utf-8"))


def _plan_identity_sha256(review: Mapping[str, Any]) -> str:
    plan = review.get("scientific_plan") or {}
    return str(plan.get("scientific_identity_sha256") or "")


def _session_events_path(session_result: Any, workspace: Path) -> Path:
    run_id = str(getattr(session_result, "run_id", "") or "")
    candidate = (
        Path(workspace) / ".chemsmart-agent" / "runs" / run_id / "events.jsonl"
    )
    if candidate.is_file():
        return candidate
    runs = sorted(
        (Path(workspace) / ".chemsmart-agent" / "runs").glob(
            "live-*/events.jsonl"
        )
    )
    if not runs:
        raise ContractError("the planning session left no event stream")
    return runs[-1]


def _goal_bound_scope(
    ledger: GoalLedger,
) -> tuple[str | None, Mapping[str, Any] | None]:
    """The scope a previous cycle established, when the goal had none.

    ``(None, None)`` when nothing bound one, which is every goal whose
    first cycle displayed an executable partition: those carry their
    scope in the digest-bound goal record itself.
    """

    for entry in reversed(list(ledger.entries())):
        if entry["kind"] != "goal_scope_bound":
            continue
        payload = entry["payload"]
        return (
            str(payload.get("scientific_identity_sha256") or ""),
            payload.get("conditions"),
        )
    return (None, None)


def _previous_run_reference(ledger: GoalLedger) -> str:
    """The evidence a revision answers -- the last run, or the last
    analysis-only cycle's own stream.

    Reading ``run_recorded`` alone left a goal whose cycles produced no
    engine run with no reference at all, so revision admission refused
    its first executable plan for never having read an outcome that was
    "unrecorded". A cycle that delivered claims and a decision over
    registered results produced real, durable, typed evidence; it simply
    launched no engine. The invariant ``evidence_read`` protects -- that
    a plan answering a failure it never read is answering a guess --
    is served by naming that evidence, not by inventing an engine run
    and not by waiving the check.
    """

    reference = ""
    for entry in ledger.entries():
        if entry["kind"] == "run_recorded":
            reference = str(entry["payload"].get("run") or "")
        elif entry["kind"] == "analysis_evidence_recorded":
            reference = str(entry["payload"].get("evidence") or "")
    return reference


#: The refusal affordance, present from cycle 1: the first live goal
#: round's honest refusal was invisible to the settlement layer partly
#: because no session was ever told the typed route exists.
#: The closing act is adversarial and points at a tool call, never at
#: re-reading one's own prose: intrinsic self-review without external
#: feedback degrades, while one further typed observation can refute.
#: The wake's rules render from the registry (chemsmart.agent.rules), so
#: they have ids and provenance there and a test can find them.
_WAKE_RULES = rules_by_id()
_GOAL_AUTHORITY = _WAKE_RULES["wake.goal_authority"].text
_OBSERVABLE_RESTATEMENT_ASK = _WAKE_RULES["wake.restate_observable"].text + " "
_ADVERSARIAL_CLOSE = _WAKE_RULES["wake.adversarial_close"].text + " "
_RECOVERY_ROUTE = " " + _WAKE_RULES["wake.recovery_route"].text
_DISPOSITION_BRANCH = " " + _WAKE_RULES["wake.disposition_branch"].text + " "
_CLAIM_BY_ID = " " + _WAKE_RULES["wake.claim_by_id_costs_no_engine_call"].text
_APPROACHES_TRIED = " " + _WAKE_RULES["wake.approaches_already_tried"].text
_MENU_DISPOSITIONS = (
    " " + _WAKE_RULES["wake.menu_dispositions_are_recorded"].text
)
_REFUSAL_AFFORDANCE = _WAKE_RULES["wake.refusal_is_a_deliverable"].text
_WORKSPACE_RECORD_RULE = " " + _WAKE_RULES["wake.workspace_record"].text
_EXCURSION_REPLICATION = (
    " " + _WAKE_RULES["wake.excursion_buys_replication"].text
)

#: The three ways a requirement short of its tolerance can be resolved.
#: Every one is a route the host can actually walk, because a named
#: route is a capability claim: three of seven repair-menu entries once
#: named routes the host could not walk, and the model asked five times
#: for the one that was promised.
#: An uncertainty is evaluated after the numbers exist, never before.
#: A planned claim node is written while the calculation is still a
#: plan, so it carries no uncertainty and cannot: the assessment is a
#: judgement about a result, and a plan has none yet. An executed
#: delivery therefore lands ``unstated`` by construction, and the honest
#: route is to look at what came back and say what it is worth. That
#: costs one further cycle and no engine call, and the cost is the
#: design rather than a defect in it (owner ruling, 2026-09-09).
SUFFICIENCY_UNSTATED_ROUTE = (
    "assess it: the number is delivered and no uncertainty stands beside "
    "it. Read your own result and claim the observable again with the "
    "uncertainty you attribute to it and its basis -- measured here, "
    "inferred from a receipt you cite, or asserted. An uncertainty is a "
    "judgement about a result, so it is made now and not at plan time; "
    "this costs no engine call."
)

#: A number within its tolerance on the session's own word, or with a
#: term of its own budget left unquantified. The number stands; the
#: requirement does not close, because a self-report cannot discharge a
#: contract that judges the self-report.
SUFFICIENCY_ATTESTED_ROUTE = (
    "back it or bound it: your uncertainty is within the tolerance and "
    "rests on your own word, or leaves a term unquantified, so the "
    "requirement stands open. Cite the receipt the magnitude comes from "
    "and claim again with basis measured or inferred; or measure the "
    "term you could not -- a conformer search, one method axis changed "
    "-- and cite that; or refuse it with the receipts that show the gap. "
    "Citing a receipt you already hold costs no engine call and so does "
    "refusing; only measuring a term you have not measured spends one. "
    "Naming a term you cannot quantify is the honest act and it is why "
    "this is open: saying so is a delivery, and dropping the term to "
    "close the word is not."
)

#: A number whose stated uncertainty does not reach the tolerance its
#: declaration asked for. Three routes, two of them free.
SUFFICIENCY_SHORT_ROUTE = (
    "compute it: plan the calculation that narrows the term you named -- "
    "the same identity, state and conditions, one method axis changed -- "
    "and claim the observable again with the uncertainty you measured "
    "and the receipt it came from; or ask a question this number can "
    "answer: declare the margin the decision actually turns on as its "
    "own observable, with the tolerance that decision needs, and deliver "
    "it like any other quantity; or refuse it: "
    "record_scientific_decision's unreachable_observable_ids naming the "
    "producer no program in this envelope provides and the receipts that "
    "show the gap. A tolerance nothing in the envelope can reach is a "
    "finding, not a failure, and saying so is a delivery."
)


def sufficiency_menu(states: Sequence[str]) -> str:
    """The routes that answer the states actually open.

    Offering "compute the calculation that narrows the term you named"
    to a session that has named no term is a route to nowhere; offering
    only "assess it" to a session whose assessment already misses the
    tolerance is no route at all. The host names what fits and never
    chooses among them.
    """

    by_state = {
        "unstated": SUFFICIENCY_UNSTATED_ROUTE,
        "attested": SUFFICIENCY_ATTESTED_ROUTE,
        "short": SUFFICIENCY_SHORT_ROUTE,
    }
    # Validated against the routes this function actually has, not
    # against the wider vocabulary: `met` is in the vocabulary and has
    # no route, so it passed the guard and reached a KeyError, and the
    # empty string passed it and returned a silently empty route.
    unknown = sorted(state for state in set(states) if state not in by_state)
    if unknown:
        # A state with no route is a state the host cannot ask anything
        # about, and falling back to "compute the calculation that
        # narrows the term you named" is the route this function exists
        # to stop offering blindly.
        raise ContractError(
            "no route answers requirement state(s) "
            + ", ".join(repr(state) for state in unknown)
            + "; a state without a route cannot be asked about"
        )
    routes = [state for state in by_state if state in set(states)]
    routes = [by_state[state] for state in routes]
    return " Or, for the rest: ".join(routes)


#: Node endings a revision can answer with ordinary work, and what that
#: work is. The host names the route; the physics decides whether it was
#: right, after the session acts. Endings absent here -- a launch that
#: never happened, an admission refusal, a cancellation, any ambiguous
#: termination -- are not evidence a revision can stand on and return to
#: the human.
REPAIR_MENU: Mapping[str, str] = {
    "failed_wrong_stationary_point": (
        "The search converged onto a stationary point of the wrong order. "
        "Read which atoms move in the offending mode "
        "(vibrational_mode_atom_participation), step the structure along "
        "it with displace_along_vibrational_mode and optimise again, or "
        "change the internal coordinate that mode moves with "
        "edit_molecular_geometry; for a transition state, seed the search "
        "from a validated frequency-bearing producer's Hessian. Or the "
        "saddle is the finding: its energy, its mode and the heavy atoms "
        "that carry it are on the node's anomalies -- name what it is "
        "before you decide whether to leave it. Its numbers were always "
        "readable as they stand; characterise_stationary_point has the "
        "host check the order you name against the printed modes, so a "
        "barrier or a free energy you deliver from that result says "
        "which structure it belongs to."
    ),
    "failed_nonconverged_scf": (
        "An SCF that will not converge is usually a state problem before "
        "it is a solver problem: check the multiplicity and charge you "
        "bound against the chemistry, then change the initial guess or "
        "the convergence route within the approved conditions. Or the "
        "failure is the finding: an SCF oscillating between two solutions "
        "or unstable under analysis is telling you about the state."
    ),
    "failed_nonconverged_geometry": (
        "Restart from the last geometry the run reached rather than the "
        "original coordinates -- a fresh start repeats the same path, "
        "and on this host one did, to the last decimal. "
        "bind_reached_geometry carries that structure forward; bind its "
        "charge and multiplicity, then optimise it in a new workflow. "
        "The optimiser settings the project exposes are the other lever: "
        "for ORCA those are the geometry maximum-iteration and "
        "convergence fields. Or the walk is the finding: a geometry that "
        "keeps moving may have left one basin for another, and the "
        "reached geometry says which."
    ),
    "failed_nonconverged_scan_step": (
        "A scan step failed to converge: the surface reached so far is "
        "readable as it stands, and bind_scan_point_geometry carries a "
        "converged point of a completed scan forward. Where the scan "
        "itself did not complete, bind_reached_geometry carries the "
        "structure the run reached. Or loosen the step's own "
        "optimisation controls."
    ),
    "timeout_terminated": (
        "The engine ran out of the time the envelope granted. Restart "
        "from the geometry the run reached -- bind_reached_geometry "
        "carries it forward -- inside the remaining budget, or reduce "
        "the method's cost within the approved conditions; conditions "
        "themselves may not move."
    ),
    "memory_limit_terminated": (
        "The engine exceeded its memory. Reduce what it holds -- basis, "
        "auxiliary basis, integral storage -- within the approved "
        "conditions, or split the calculation; the resources are the "
        "envelope's and cannot be raised here."
    ),
    "failed_native": (
        "The program stopped on its own error. Read the native findings "
        "and the engine's last lines on the run outcome; a program error "
        "that names a setting is repaired in project YAML, one that "
        "names the molecule is a new decision for the human. An input-check "
        "abort reached no chemistry: nothing was reached to restart from, "
        "and the repair is the field the engine named."
    ),
}

#: Node endings a revision can answer (the repair menu's keys), which
#: are the one set terminal_states owns; the menu may not drift from it.
REPAIRABLE_TERMINAL_STATES = REPAIRABLE_NODE_STATES
if frozenset(REPAIR_MENU) != REPAIRABLE_TERMINAL_STATES:  # pragma: no cover
    raise ImportError("the repair menu and the repairable endings differ")


def _goal_terms_context(
    *,
    goal_id: str,
    granted_by: str,
    envelope_record: Mapping[str, Any],
    max_revisions: int,
    workspace: Path | None = None,
) -> dict[str, Any]:
    """Cycle 1's context: the goal's terms before any run exists.

    The first live goal round handed cycle 1 nothing -- an
    analysis-only goal is single-cycle by construction, so its session
    never saw the budgets (a zero engine-call grant would have said
    "analysis only" before the session drafted engine work), the
    authority in force, or the refusal affordance. The terms are all
    in hand when the loop starts; only the trajectory is empty.
    """

    return {
        "schema_version": "chemsmart.goal-wake-context.v1",
        "goal_id": goal_id,
        "granted_by": granted_by,
        "conditions": {},
        "budgets": {
            "binding_line": (
                "nearest exhausted: engine calls "
                f"{int(envelope_record.get('max_engine_calls') or 0)} of "
                f"{int(envelope_record.get('max_engine_calls') or 0)} "
                "remaining (100%)"
            ),
            "engine_calls_remaining": int(
                envelope_record.get("max_engine_calls") or 0
            ),
            "wall_seconds_remaining": float(
                envelope_record.get("episode_wall_time_seconds") or 0.0
            ),
            "revisions_remaining": int(max_revisions),
            # The line was absent here, so the host read None and every
            # first cycle displayed 0 excursion calls while the envelope
            # granted 2 (REACH-1, both goals).
            "excursion_calls_remaining": int(
                envelope_record.get("max_excursion_calls") or 0
            ),
        },
        # Cycle 1 has delivered nothing yet, but it carries the same
        # keys a wake does: one shape across every cycle is what lets a
        # session read the record the same way each time.
        "deliverables": _deliverables_record(
            _AnalysisDelivery("", (), 0, 0, ())
        ),
        "trajectory": (),
        "previous_run": "",
        "previous_run_outcome": {},
        "declared_observables": (),
        # What earlier goals in this workspace computed and claimed,
        # host-written from their receipts; empty in a fresh workspace.
        "workspace_record": (
            render_workspace_record(workspace) if workspace is not None else {}
        ),
        "authority": (
            _GOAL_AUTHORITY
            + " This session plans cycle 1; the budgets above are the "
            "whole grant. "
            + _OBSERVABLE_RESTATEMENT_ASK
            + _ADVERSARIAL_CLOSE
            + _REFUSAL_AFFORDANCE
            + _WORKSPACE_RECORD_RULE
        ),
    }


def _goal_envelope_record(shown: Mapping[str, Any]) -> dict[str, Any]:
    """The budget lines a goal is granted, from the envelope it was shown.

    Two hand-listed copies of this record dropped the excursion line, so
    every woken cycle of every excursion arm in two sealed windows was
    told it had zero excursions remaining and the plan gate refused the
    line the human had granted (E4, E4', 2026-09-03). One record, every
    line the ledger's budgets read.
    """

    return {
        "allowed_program_engines": canonical_data(
            shown.get("allowed_program_engines") or ()
        ),
        "max_engine_calls": int(shown.get("max_engine_calls") or 0),
        "episode_wall_time_seconds": float(
            shown.get("episode_wall_time_seconds") or 0.0
        ),
        "max_excursion_calls": int(shown.get("max_excursion_calls") or 0),
    }


def _deliverables_record(delivery: _AnalysisDelivery) -> dict[str, Any]:
    """What the previous run's own stream says stands delivered.

    Names quantities and stated limitations, never values: the goal's
    demand is in the task, and this record lets a wake session see what
    it has already delivered, what the chain declared it could not, and
    what its own decisions doubt -- so the next action can follow the
    gap rather than the tool list.
    """

    return {
        "delivered_quantity_ids": delivery.delivered_quantity_ids,
        "limitation_output_ids": delivery.limitation_output_ids,
        "doubted_quantity_ids": delivery.doubted_quantity_ids,
        "unanswered_failed_verdicts": delivery.unanswered_verdicts,
        "stale_quantity_ids": delivery.stale_quantity_ids,
        "unclaimed_output_ids": delivery.unclaimed_output_ids,
        "undelivered_declared_observable_ids": (
            delivery.undelivered_declared_ids
        ),
        "unresolved_requirement_ids": delivery.unresolved_requirement_ids,
        "sufficiency": delivery.sufficiency,
        "flagged_quantity_ids": delivery.flagged_quantity_ids,
        "failed_source_quantity_ids": delivery.failed_source_quantity_ids,
        "characterised_source_quantity_ids": (
            delivery.characterised_source_quantity_ids
        ),
        "uncharacterised_source_quantity_ids": (
            delivery.uncharacterised_source_quantity_ids
        ),
        "expectation_rows": delivery.prediction_rows,
    }


def _goal_anomalies(ledger: GoalLedger) -> tuple[dict[str, Any], ...]:
    """Every anomaly any cycle of the goal recorded, one per receipt."""

    seen: dict[str, dict[str, Any]] = {}
    for entry in ledger.entries():
        if entry["kind"] != "anomalies_observed":
            continue
        for record in entry["payload"].get("anomalies") or ():
            digest = str(record.get("receipt_sha256") or "")
            if digest and digest not in seen:
                seen[digest] = dict(record)
    return tuple(seen.values())


def _flagged_artifact_sha256s(
    records: Sequence[Mapping[str, Any]],
) -> tuple[str, ...]:
    """Artifacts of every node an anomaly still standing has flagged."""

    return tuple(
        sorted(
            {
                str(digest)
                for record in anomaly_standing(records)
                for digest in record.get("flagged_artifact_sha256s") or ()
                if digest
            }
        )
    )


def _anomaly_output_ids_from_records(
    records: Sequence[Mapping[str, Any]],
) -> tuple[str, ...]:
    # Standing is the head of each chain: a replicated or refuted
    # receipt speaks for the anomaly it supersedes.
    return tuple(
        sorted(
            {
                "anomaly:"
                f"{record.get('signal_id')}:{record.get('status')}:"
                f"{str(record.get('receipt_sha256') or '')[:8]}"
                for record in anomaly_standing(records)
                if record.get("signal_id") and record.get("receipt_sha256")
            }
        )
    )


def _first_declarations(ledger: GoalLedger) -> tuple[dict[str, Any], ...]:
    """The goal's first declaration of each observable, in ledger order.

    A woken session re-declared its expectations with a flipped sign
    convention and wider bands, and the completion row printed agreed
    over a falsified first prior (live, 2026-09-02). An expectation is
    a prediction only if it predates the physics, so the earliest
    declaration of an id is the one every later session's host is
    seeded with.
    """

    first: dict[str, dict[str, Any]] = {}
    for entry in ledger.entries():
        if entry["kind"] != "observables_declared":
            continue
        for record in entry["payload"].get("observables") or ():
            observable_id = str(record.get("observable_id") or "")
            if observable_id and observable_id not in first:
                first[observable_id] = dict(record)
    return tuple(first.values())


def _goal_delivered_ids(
    workspace: Path | None, goal_id: str
) -> dict[str, dict[str, Any]]:
    """Declared ids this goal has delivered under their id, in any cycle.

    Read from the workspace record, which every recorded run and every
    in-session delivery appends to: the claim's claim_id or quantity_id
    is the key, the latest cycle wins, and the row says where the number
    came from. A settlement that read one stream called two observables
    undelivered that cycle 2 had delivered under their exact ids (NOVEL-3
    ino3, 2026-09-05); the record held both rows in the same directory.
    """

    if workspace is None:
        return {}
    delivered: dict[str, dict[str, Any]] = {}
    for entry in read_workspace_record(workspace):
        if entry.get("kind") != "claim" or entry.get("goal_id") != goal_id:
            continue
        for id_field in ("claim_id", "quantity_id"):
            key = str(entry.get(id_field) or "")
            if not key:
                continue
            previous = delivered.get(key)
            if previous is None or int(entry.get("cycle") or 0) >= int(
                previous.get("cycle") or 0
            ):
                delivered[key] = {
                    "cycle": int(entry.get("cycle") or 0),
                    "run": entry.get("run"),
                    "value": entry.get("value"),
                    "unit": entry.get("unit"),
                    "dimension": tuple(entry.get("dimension") or ()),
                    "claim_id": entry.get("claim_id"),
                    # The assessment travels with the number, so a later
                    # cycle can read what an earlier one established
                    # about its precision.
                    "sufficiency": entry.get("sufficiency"),
                }
    return delivered


def _required_declared_ids(ledger: GoalLedger) -> tuple[str, ...]:
    """The declared ids the goal must deliver: every first declaration
    except the diagnostics.

    A diagnostic is the session's own prediction about the route, given
    standing so it is scored (owner ruling R3, 2026-09-06). It is never
    a deliverable: it does not hold a settlement open, does not earn a
    re-wake, and never prints as an undelivered id.
    """

    from chemsmart.agent.delivery import superseded_observable_ids

    declarations = _first_declarations(ledger)
    retired = superseded_observable_ids(declarations)
    return tuple(
        str(record.get("observable_id") or "")
        for record in declarations
        if str(record.get("role") or "requested") != "diagnostic"
        and str(record.get("observable_id") or "") not in retired
    )


def _session_dispositions(events_path: Path | None) -> tuple[dict, ...]:
    """Every repair-menu disposition a session's decisions recorded."""

    if events_path is None:
        return ()
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    dispositions: list[dict] = []
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "scientific_decision_recorded":
            continue
        payload = event.get("payload") or {}
        for item in payload.get("menu_route_dispositions") or ():
            if isinstance(item, Mapping):
                dispositions.append(dict(item))
    return tuple(dispositions)


def _session_input_checks(events_path: Path | None) -> dict[str, Any]:
    """How the session's input-check probes concluded, and why.

    The counts alone were the whole row, and that is half a repair:
    843e04f2 restored the *number* to the ledger and left the *reason*
    behind, so a wake could say `aborted: 2` without saying what the
    program objected to. The engine's own lines are what a next cycle
    can act on -- ORCA names the three legal keywords in its abort --
    and re-learning them from a dead run costs an engine call per node.
    """

    counts: dict[str, Any] = {"passed": 0, "aborted": 0, "not_run": 0}
    reasons: dict[str, list[str]] = {}
    if events_path is None:
        return counts
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return counts
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "input_check_probed":
            continue
        payload = event.get("payload") or {}
        status = str(payload.get("status") or "")
        if status in counts:
            counts[status] += 1
        if status == "aborted":
            node = str(payload.get("node_id") or "")
            lines = [str(line) for line in (payload.get("engine_lines") or ())]
            if lines:
                reasons[node] = lines[:6]
    if reasons:
        # Keyed by node, because the repair is per node and the model
        # needs to know which input the program refused.
        counts["aborted_engine_lines"] = {
            node: tuple(lines) for node, lines in sorted(reasons.items())
        }
    return counts


#: The typed reads a cycle can perform without launching an engine.
#: Any one of them is the host reading a real quantity out of real
#: bytes and minting a receipt for it, which is evidence a later plan
#: can answer -- the thing the revision gate actually asks about.
_TYPED_READ_KINDS = (
    "result_quantities_extracted",
    "thermochemistry_derived",
    "quantity_expression_evaluated",
    "stationary_point_characterised",
)


def _session_typed_reads(events_path: Path | None) -> int:
    """How many typed reads this session's own stream recorded.

    po3-r17 cycle 1 (2026-09-11) extracted eight quantities and recorded
    a decision, rendered no claim, and was therefore named as evidence by
    nothing -- so its re-wake's first executable plan was refused for
    answering an outcome that was "unrecorded". Reading a quantity is not
    claiming one, and the gate's own invariant is about having *read*.
    """

    if events_path is None:
        return 0
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return 0
    total = 0
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") in _TYPED_READ_KINDS:
            total += 1
    return total


def _recorded_input_checks(ledger: GoalLedger) -> dict[str, Any]:
    """Every input-check probe the goal's sessions ran, summed, never
    charged: engine calls derive from execution receipts alone."""

    total = {"passed": 0, "aborted": 0, "not_run": 0}
    #: Why each abort happened, in the program's own words, carried into
    #: the wake beside the count. A wake that says `aborted: 2` and
    #: nothing else tells the next cycle that something is wrong and not
    #: what -- so the diagnosis was re-bought from the dead run at one
    #: engine call per node, twice, in two windows.
    reasons: dict[str, tuple[str, ...]] = {}
    for entry in ledger.entries():
        if entry["kind"] != "input_checks_probed":
            continue
        payload = entry["payload"]
        for key in total:
            total[key] += int(payload.get(key) or 0)
        for node, lines in (payload.get("aborted_engine_lines") or {}).items():
            reasons[str(node)] = tuple(str(line) for line in lines)
    summary: dict[str, Any] = {**total, "charged": 0}
    if reasons:
        summary["aborted_engine_lines"] = reasons
    return summary


def _session_approaches(events_path: Path | None) -> tuple[dict, ...]:
    """What this session tried and rejected, in its own words.

    A rejected alternative carries the mechanism the session reasoned
    with, which is exactly what the next cycle needs and exactly what
    the wake did not carry: REACH-1 po3 rejected restarting from a
    reached geometry because the structure was a separated pair, then
    the next cycle re-seeded a structurally identical dual contact
    (2026-09-06).
    """

    if events_path is None:
        return ()
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    approaches: list[dict] = []
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "scientific_decision_recorded":
            continue
        record = (event.get("payload") or {}).get("record") or {}
        for text in record.get("alternatives") or ():
            statement = str(text).strip()
            if statement:
                approaches.append(
                    {
                        "approach": statement[:400],
                        "outcome": "rejected_by_the_session",
                    }
                )
    return tuple(approaches)


def _run_approaches(run_events_path: Path | None) -> tuple[dict, ...]:
    """How each node that did not deliver actually ended, with numbers.

    The repair menu names routes for a terminal state; this names the
    node, the state and the sensor numbers under it, so a later cycle
    can see that this approach has already been paid for.
    """

    if run_events_path is None or not run_events_path.is_file():
        return ()
    states: dict[str, str] = {}
    signals: dict[str, list[str]] = {}
    for line in run_events_path.read_text(encoding="utf-8").splitlines():
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        payload = event.get("payload") or {}
        node_id = str(payload.get("node_id") or "")
        if not node_id:
            continue
        if event.get("kind") == "program_result_verified":
            states[node_id] = str(payload.get("status") or "")
        elif event.get("kind") == "anomaly_observed":
            record = payload.get("record") or {}
            values = record.get("values") or {}
            numbers = ", ".join(
                f"{key}={value}"
                for key, value in sorted(values.items())
                if isinstance(value, (int, float))
            )
            signal = str(record.get("signal_id") or "")
            if signal:
                signals.setdefault(node_id, []).append(
                    f"{signal}({numbers})" if numbers else signal
                )
    approaches = []
    for node_id, status in sorted(states.items()):
        if status == "valid":
            continue
        approaches.append(
            {
                "approach": node_id,
                "outcome": status or "invalid",
                "mechanism": "; ".join(signals.get(node_id, ())) or "",
            }
        )
    return tuple(approaches)


def _recorded_approaches(
    ledger: GoalLedger, limit: int = 12
) -> tuple[dict[str, Any], ...]:
    """Everything this goal has already tried, most recent first."""

    entries = []
    for entry in ledger.entries():
        if entry["kind"] != "approaches_recorded":
            continue
        cycle = entry["payload"].get("cycle")
        for item in entry["payload"].get("approaches") or ():
            entries.append({**dict(item), "cycle": cycle})
    return tuple(reversed(entries))[:limit]


def _recorded_dispositions(ledger: GoalLedger) -> tuple[dict[str, Any], ...]:
    """The previous cycle's dispositions of the menu it was offered,
    each carrying the cycle that made it."""

    latest: tuple[dict[str, Any], ...] = ()
    for entry in ledger.entries():
        if entry["kind"] != "repair_menu_dispositions":
            continue
        cycle = entry["payload"].get("cycle")
        latest = tuple(
            {**dict(item), "cycle": cycle}
            for item in entry["payload"].get("dispositions") or ()
        )
    return latest


def _session_declarations(events_path: Path | None) -> tuple[dict, ...]:
    """Every observable a session stream declared, in stream order."""

    if events_path is None:
        return ()
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    declared: list[dict] = []
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "requested_observable_declared":
            continue
        for record in (event.get("payload") or {}).get("observables") or ():
            declared.append(dict(record))
    return tuple(declared)


def _host_seconds_spent(ledger: GoalLedger) -> float:
    """The host's own post-processing time this goal has paid, summed from
    its recorded runs. Shown beside the engine lines and charged to none
    of them: the goal clock counts engine seconds (owner's ruling)."""

    total = 0.0
    for entry in ledger.entries():
        if entry["kind"] == "run_recorded":
            total += float(entry["payload"].get("host_seconds") or 0.0)
    return float(f"{total:.1f}")


def _binding_budget_line(
    goal: GoalRecordV1, budgets: Any, ledger: GoalLedger
) -> str:
    """The budget line nearest exhaustion, and whether the slowest node
    this goal ran would fit in the engine wall that remains.

    Four remaining numbers with no grant beside them read as plenty: a
    session with 12 of 30 engine calls left and 764 of 21 600 wall
    seconds planned three optimisations of a metal complex whose last
    one had taken 4 000 s (NOVEL-3 ino3, 2026-09-05). The wall, not the
    calls, was the binding line, and nothing said so.
    """

    envelope = dict(goal.envelope)
    lines = [
        (
            "engine calls",
            float(budgets.engine_calls_remaining),
            float(envelope.get("max_engine_calls") or 0),
            "",
        ),
        (
            "engine wall",
            float(budgets.wall_seconds_remaining),
            float(envelope.get("episode_wall_time_seconds") or 0.0),
            " s",
        ),
        (
            "revisions",
            float(budgets.revisions_remaining),
            float(goal.max_revisions),
            "",
        ),
    ]
    excursions = float(envelope.get("max_excursion_calls") or 0)
    if excursions:
        lines.append(
            (
                "excursion calls",
                float(budgets.excursion_calls_remaining),
                excursions,
                "",
            )
        )
    name, remaining, grant, unit = min(
        lines, key=lambda item: (item[1] / item[2]) if item[2] else 1.0
    )
    share = int(round(100.0 * remaining / grant)) if grant else 100
    line = (
        f"nearest exhausted: {name} {remaining:g} of {grant:g}{unit} "
        f"remaining ({share}%)"
    )
    slowest_seconds = 0.0
    slowest_node = ""
    for entry in ledger.entries():
        if entry["kind"] != "run_recorded":
            continue
        seconds = float(entry["payload"].get("slowest_node_seconds") or 0.0)
        if seconds > slowest_seconds:
            slowest_seconds = seconds
            slowest_node = str(entry["payload"].get("slowest_node_id") or "")
    if slowest_seconds > 0:
        fits = int(budgets.wall_seconds_remaining // slowest_seconds)
        line += (
            f"; the slowest node this goal ran took {slowest_seconds:.0f} s"
            f"{f' ({slowest_node})' if slowest_node else ''}, and the "
            f"engine wall remaining fits {fits} such node"
            f"{'' if fits == 1 else 's'}"
        )
    return line


def _wake_context(
    goal: GoalRecordV1,
    ledger: GoalLedger,
    outcome: Any,
    *,
    workspace: Path | None = None,
    failure_report: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    budgets = ledger.budgets(goal)
    trajectory = tuple(
        {
            "kind": entry["kind"],
            "at": entry["at"],
            "payload": {
                key: value
                for key, value in entry["payload"].items()
                if key
                in {
                    "cycle",
                    "engine_calls_consumed",
                    "excursion_calls_consumed",
                    "engine_wall_seconds",
                    "workflow_state",
                    "run",
                    "reasons",
                }
            },
        }
        for entry in ledger.entries()
        if entry["kind"]
        in {"run_recorded", "revision_admitted", "revision_returned"}
    )
    # The evidence a revision answers, whether an engine produced it or
    # an analysis-only cycle did. Gating this on `outcome is not None`
    # meant an analysis-only cycle embedded nothing, so admission
    # compared a real reference against an empty one and refused the
    # goal's first calculation -- A10 named the evidence and the wake
    # then declined to carry it.
    previous_run = _previous_run_reference(ledger)
    deliverables: dict[str, Any] = {
        "delivered_quantity_ids": (),
        "limitation_output_ids": (),
        "doubted_quantity_ids": (),
    }
    if previous_run and workspace is not None:
        deliverables = _deliverables_record(
            _analysis_delivery(
                Path(workspace)
                / ".chemsmart-agent"
                / Path(*previous_run.split("/"))
                / "events.jsonl",
                uncharacterised_artifact_sha256s=uncharacterised_artifacts(
                    workspace
                ),
                declared_observables=_first_declarations(ledger),
                goal_delivered_ids=_goal_delivered_ids(
                    workspace, goal.goal_id
                ),
            )
        )
    failure_report = dict(failure_report or {})
    approaches_tried = _recorded_approaches(ledger)
    repair_menu = {
        state: REPAIR_MENU[state]
        for state in sorted(
            {
                str(getattr(node, "state", "") or "")
                for node in (outcome.nodes if outcome is not None else ())
            }
        )
        if state in REPAIR_MENU
    }
    # A requirement left unresolved by an executed run re-opens the goal
    # through recovery_opened, which is keyed on how nodes ended and so
    # carried no route for it at all: the session was woken and told
    # nothing. The route rides the menu the session already reads, keyed
    # by the state of the requirement rather than of a node, so the
    # dispositions machinery records what was done with it as it does
    # for any other offered route.
    for state in sorted(
        {
            str(row.get("state") or "")
            for row in deliverables.get("sufficiency") or ()
            if str(row.get("observable_id") or "")
            in set(deliverables.get("unresolved_requirement_ids") or ())
        }
    ):
        if state == "unstated":
            repair_menu["requirement_unstated"] = SUFFICIENCY_UNSTATED_ROUTE
        elif state == "attested":
            repair_menu["requirement_attested"] = SUFFICIENCY_ATTESTED_ROUTE
        elif state == "short":
            repair_menu["requirement_short"] = SUFFICIENCY_SHORT_ROUTE
    repair_sentence = (
        " repair_menu names, for each way a node of the previous run "
        "ended, the ordinary route that answers it; the host does not "
        "say which route is right, the next run does. "
        if repair_menu
        else ""
    )
    return {
        "schema_version": "chemsmart.goal-wake-context.v1",
        "goal_id": goal.goal_id,
        "granted_by": goal.granted_by,
        # The conditions this goal is actually held to. A goal whose
        # first cycle read registered results carries empty conditions
        # in its digest-bound record and binds the real ones on its
        # first executable revision, so rendering the record showed a
        # woken session an empty scope while admission held it to the
        # bound one.
        "conditions": canonical_data(
            _goal_bound_scope(ledger)[1] or goal.conditions
        ),
        "budgets": {
            "binding_line": _binding_budget_line(goal, budgets, ledger),
            "host_seconds_spent": _host_seconds_spent(ledger),
            "engine_calls_remaining": budgets.engine_calls_remaining,
            "wall_seconds_remaining": budgets.wall_seconds_remaining,
            "revisions_remaining": budgets.revisions_remaining,
            "excursion_calls_remaining": budgets.excursion_calls_remaining,
            "input_check_probes": _recorded_input_checks(ledger),
        },
        "deliverables": deliverables,
        "trajectory": trajectory,
        "previous_run": previous_run,
        "previous_run_outcome": (
            outcome.public_record() if outcome is not None else {}
        ),
        "repair_menu": repair_menu,
        "repair_menu_dispositions": _recorded_dispositions(ledger),
        "approaches_tried": approaches_tried,
        "declared_observables": _first_declarations(ledger),
        "diagnostic_observables": tuple(
            record
            for record in _first_declarations(ledger)
            if str(record.get("role") or "requested") == "diagnostic"
        ),
        "anomalies": _goal_anomalies(ledger),
        "workspace_record": (
            render_workspace_record(workspace) if workspace is not None else {}
        ),
        # Present only on a re-wake: the previous cycle of this same
        # goal ended without delivering, refusing, or planning, and this
        # is the one further wake the budget grants it.
        **({"failure_report": failure_report} if failure_report else {}),
        "authority": (
            _GOAL_AUTHORITY
            + " The previous run's typed outcome is embedded above as "
            "previous_run_outcome; inspect_run re-reads it by the "
            "previous_run reference or by its run_id, or any earlier "
            "run. "
            + _RECOVERY_ROUTE
            + repair_sentence
            + (_APPROACHES_TRIED + " " if approaches_tried else "")
            + (_MENU_DISPOSITIONS + " " if repair_menu else "")
            + _DISPOSITION_BRANCH
            + _CLAIM_BY_ID
            + _OBSERVABLE_RESTATEMENT_ASK
            + _ADVERSARIAL_CLOSE
            + _REFUSAL_AFFORDANCE
            + _WORKSPACE_RECORD_RULE
            + (
                _EXCURSION_REPLICATION
                if budgets.excursion_calls_remaining
                else ""
            )
        ),
    }


def _achieved(execute_result: Any) -> bool:
    status = str(getattr(execute_result, "status", "") or "")
    analysis = str(getattr(execute_result, "analysis_status", "") or "")
    return status == "completed" and analysis in {"completed", ""}


@dataclass(frozen=True)
class _AnalysisDelivery:
    """What one session's durable stream says it delivered."""

    completion_status: str
    limitation_output_ids: tuple[str, ...]
    claims: int
    decisions: int
    receipt_sha256s: tuple[str, ...]
    #: Claim quantities whose supporting receipt a recorded decision
    #: doubts (``doubt:{receipt}`` evidence references intersected with
    #: the rendered claims' source receipts) -- computed here from the
    #: stream so a doubt recorded after the completion event still
    #: reaches the settlement.
    doubted_quantity_ids: tuple[str, ...] = ()
    #: Every quantity a rendered claim delivered, by its own id.
    delivered_quantity_ids: tuple[str, ...] = ()
    #: Delivered quantities that descend from a node whose run did not
    #: meet the promise it was launched under. Never a verdict: the
    #: number may be exactly the finding worth reporting, and it is
    #: readable by design. It is the join a reader needs, because the
    #: settlement never said which delivered value came from there.
    failed_source_quantity_ids: tuple[str, ...] = ()
    #: The subset of those whose source the session had the host check,
    #: so the delivered number names the structure it belongs to.
    characterised_source_quantity_ids: tuple[str, ...] = ()
    #: Quantities standing on an optimisation or transition-state search
    #: whose result printed no frequencies: its stationary point is
    #: uncharacterised, and the number says so.
    uncharacterised_source_quantity_ids: tuple[str, ...] = ()
    #: The recorded decision's own stated uncertainties, verbatim.
    decision_uncertainties: tuple[str, ...] = ()
    #: Result artifacts a session characterised, with the host checking
    #: the stationary-point order against the printed modes. A delivery
    #: standing on one of these stands on a checked statement.
    characterised_artifact_sha256s: tuple[str, ...] = ()
    #: Delivered quantities that descend from a node the host flagged
    #: with an anomaly. Never a verdict: the flagged run may be exactly
    #: right. It is the join a reader needs, because the settlement word
    #: named the anomaly and never said the number came from it.
    flagged_quantity_ids: tuple[str, ...] = ()
    #: Host-rendered verdicts that failed and that no recorded decision
    #: has cited. A failed verdict is the host saying the delivered
    #: structure is not what the task required -- a minimum that is a
    #: saddle, a transition state with the wrong imaginary-mode count.
    #: It never made a goal partial and it never opened a cycle, so a
    #: run could deliver every required output, display "failed" in its
    #: own report, and settle achieved with budget in hand. A decision
    #: that cites the validation receipt has answered it; one that does
    #: not has left it open.
    unanswered_verdicts: tuple[str, ...] = ()
    #: Delivered quantities whose own receipt lineage traces back to a
    #: result a failed verdict rejected. A recovery cycle that replaces
    #: the structure does not replace the numbers computed from the old
    #: one, and the wake record used to name those numbers as delivered
    #: in the same breath as the verdict invalidating them -- so a live
    #: session cleared the verdict, recomputed the value into an
    #: expression receipt, never rendered it as a claim, and settled
    #: with the superseded number standing as the answer. Stale is not
    #: wrong: the arithmetic held, the structure beneath it did not.
    stale_quantity_ids: tuple[str, ...] = ()
    #: Artifact digests a verdict rejected. A rejection is a fact about
    #: bytes and does not expire, so a goal carries these across cycles:
    #: a later cycle may render a claim from a result an earlier cycle
    #: already rejected.
    rejected_artifact_sha256s: tuple[str, ...] = ()
    #: Whether this run rendered any claim at all. A run that rendered
    #: none did not replace the standing delivery -- which is exactly how
    #: a recovery that fixed the structure and claimed nothing left the
    #: superseded number standing as the goal's answer.
    claims_rendered: bool = False
    #: Quantities an expression exported and no claim ever rendered. The
    #: host computed them from real program output and no reader of the
    #: delivery can see them: across the recorded campaign 89 of 242
    #: exported quantities were never claimed, and five goals exported
    #: quantities and claimed none at all -- four of those settled
    #: achieved with engine calls still unspent. An expression that is
    #: evaluated and never claimed is not delivered.
    unclaimed_output_ids: tuple[str, ...] = ()
    #: What the host observed that nobody asked for, from the completion
    #: receipt; a certified delivery carrying any settles with the word.
    anomaly_output_ids: tuple[str, ...] = ()
    #: The completion's own words for each declared-observable miss; the
    #: settlement quoted only the ids and dropped "of matching dimension"
    #: (NOVEL-3 ino2, 2026-09-05).
    declared_observable_misses: tuple[str, ...] = ()
    #: Declared observables the recorded decision named unreachable, split
    #: by whether the host could verify the named producer's absence
    #: (owner ruling 2026-09-06: only a verified refusal settles the word).
    verified_unreachable_ids: tuple[str, ...] = ()
    unverified_unreachable_ids: tuple[str, ...] = ()
    #: The completion's expectation rows, requested and diagnostic: what
    #: each pre-registered prediction said beside what was delivered,
    #: so the next cycle reads its own score rather than re-declaring.
    prediction_rows: tuple[dict[str, Any], ...] = ()
    #: What the session did with each route the wake's repair menu
    #: offered, host-verified against that menu.
    route_dispositions: tuple[dict[str, Any], ...] = ()
    unreachable_bases: Mapping[str, str] = field(default_factory=dict)
    #: What each claim of this stream carries, under both the names it
    #: answers to, so a current-cycle claim is judged in the dimension
    #: its declaration asked for exactly as a record row is.
    claim_rows: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    #: The runtime's own terminal reason, kept typed rather than read
    #: back out of the rendered ending sentence: a cycle the provider's
    #: transport ended is a different fact from a cycle the science
    #: ended, and only one of them should spend a scientific wake.
    terminal_reason: str = ""
    #: How the session ended, in the stream's own facts: workflows
    #: planned, nodes previewed, a review the host refused, the last
    #: refused plan call, the runtime's terminal reason. A settlement or a
    #: re-wake names this, never a bare terminal word (REACH-1 ino3 was
    #: returned as "runtime terminal state is absorbing").
    ending: str = ""
    #: Declared ids this goal delivered under their id in an earlier cycle
    #: (id -> {cycle, run, value, unit}); a settlement that read one
    #: stream called two of them undelivered (NOVEL-3 ino3).
    #: The goal's first declaration of each observable, so this
    #: projection judges a record row in the dimension the declaration
    #: asked for -- the same comparison the completion gate makes.
    declared_observables: tuple[Mapping[str, Any], ...] = ()
    goal_delivered: Mapping[str, Mapping[str, Any]] = field(
        default_factory=dict
    )
    #: How each delivered number stands against the precision its own
    #: declaration asked for: met, attested, short, or unstated. The
    #: task's own tolerance had no home in the
    #: goal at all, so a session could deliver a number its own
    #: limitations said missed the requester's tolerance and the
    #: settlement quoted that sentence into its reasons and said
    #: achieved (OPEN-2 ino3-qwen, 2026-09-07).
    sufficiency: tuple[Mapping[str, Any], ...] = ()

    @property
    def retired_observable_ids(self) -> frozenset[str]:
        """Ids a later declaration explicitly replaced."""

        from chemsmart.agent.delivery import superseded_observable_ids

        return superseded_observable_ids(self.declared_observables)

    @property
    def unresolved_requirement_ids(self) -> tuple[str, ...]:
        """Requirements a delivered number has not answered yet.

        Open when the stated uncertainty exceeds the tolerance, open
        when it is within the tolerance on the session's own word
        alone, and open when a claim on a tolerance-bearing observable
        states no uncertainty at all -- silence is not sufficiency.
        """

        # A retired id stops being owed, so an assessment of the
        # observable it replaced never keeps the goal open -- and
        # neither does one the host itself verified as unreachable.
        # Refusing with receipts is route three of the menu the wake
        # offers and the charter calls it a deliverable, but only the
        # executed settlement subtracted refusals, with its own ad-hoc
        # filter; the analysis-only settlement checked unresolved first
        # and short-circuited, so on that path the route could not
        # change the word at all. One subtraction, here, for every
        # reader.
        closed = self.retired_observable_ids | frozenset(
            self.verified_unreachable_ids
        )
        return tuple(
            observable_id
            for observable_id in unresolved_requirement_ids(self.sufficiency)
            if observable_id not in closed
        )

    @property
    def refused_requirement_ids(self) -> tuple[str, ...]:
        """Open requirements the session refused and the host verified.

        These are the ids ``unresolved_requirement_ids`` subtracts. The
        subtraction is right -- a refused requirement is not owed --
        but nothing read what it removed, so a goal whose stated
        precision was refused as unreachable settled ``achieved``: the
        refusal bought the best word in the vocabulary rather than the
        deliverable the charter promises it.
        """

        from chemsmart.agent.delivery import refused_requirement_ids

        return tuple(
            observable_id
            for observable_id in refused_requirement_ids(
                self.sufficiency, self.verified_unreachable_ids
            )
            if observable_id not in self.retired_observable_ids
        )

    @property
    def undelivered_declared_ids(self) -> tuple[str, ...]:
        """Declared observables no claim carried by id.

        They ride the completion's limitation list under their own
        prefix. They are not a declared-blocked producer (a typed
        refusal) and they are not a delivered headline: two live goals
        settled achieved over them (E4 window, 2026-09-03).
        """

        prefix = "declared_observable:"
        return tuple(
            item[len(prefix) :]
            for item in self.limitation_output_ids
            if item.startswith(prefix)
            and not self._delivered_by_the_goal(item[len(prefix) :])
        )

    def answers_declaration(self, observable_id: str) -> bool:
        """Whether anything this goal holds answers the declaration.

        A record row from an earlier cycle or a claim of this stream,
        both judged in the dimension the declaration asked for. The
        re-wake used to dimension-check the first and take the second on
        its id alone, so a claim in the wrong dimension counted as
        delivered here and undelivered at settlement -- and the cycle
        that could have repaired it was never opened.
        """

        if self._delivered_by_the_goal(observable_id):
            return True
        from chemsmart.agent.delivery import (
            declarations_by_id,
            observable_is_delivered,
        )

        row = self.claim_rows.get(observable_id)
        if row is None:
            return False
        declaration = declarations_by_id(self.declared_observables).get(
            observable_id
        )
        if declaration is None:
            # Nothing declared it, so nothing constrains its dimension.
            return True
        return observable_is_delivered(declaration, row)

    def _delivered_by_the_goal(self, observable_id: str) -> bool:
        """Whether a record row answers this declaration, in dimension.

        Reading every cycle rescued claims a single stream forgot; the
        id-only union then let a dimension mismatch through, and a goal
        settled achieved over six observables its own completion gate
        had called undelivered (OPEN-1 ino3, 2026-09-07). Both readers
        now call one predicate.
        """

        from chemsmart.agent.delivery import (
            declarations_by_id,
            observable_is_delivered,
        )

        row = self.goal_delivered.get(observable_id)
        if row is None:
            return False
        declaration = declarations_by_id(self.declared_observables).get(
            observable_id
        )
        if declaration is None:
            # No declaration reached this projection: the id join is all
            # this reader can honestly make, and it says so by keeping
            # the older behaviour rather than inventing a dimension.
            return True
        return observable_is_delivered(declaration, row)

    @property
    def delivered_in_earlier_cycles(self) -> tuple[str, ...]:
        """Declared ids this stream missed that an earlier cycle delivered."""

        prefix = "declared_observable:"
        return tuple(
            f"{item[len(prefix):]} (delivered in cycle "
            f"{self.goal_delivered[item[len(prefix):]].get('cycle')})"
            for item in self.limitation_output_ids
            if item.startswith(prefix)
            and self._delivered_by_the_goal(item[len(prefix) :])
        )

    @property
    def open_declared_misses(self) -> tuple[str, ...]:
        """The completion's miss text for each id still undelivered."""

        open_ids = set(self.undelivered_declared_ids)
        return tuple(
            text
            for text in self.declared_observable_misses
            if any(f"'{observable_id}'" in text for observable_id in open_ids)
        )

    @property
    def blocked_output_ids(self) -> tuple[str, ...]:
        """Limitations that are declared-blocked producers, not id misses."""

        return tuple(
            item
            for item in self.limitation_output_ids
            if not item.startswith("declared_observable:")
        )

    #: What stopped the run before a delivery could exist: a node the
    #: executor refused to launch, or an execution review the session's
    #: host refused to build. Each is one sentence the settlement quotes,
    #: because a returned goal that names nothing has not settled.
    stopped_by: tuple[str, ...] = ()


def _stale_quantity_ids(
    *,
    claim_pairs: Sequence[tuple[str, str]],
    rejected_bindings: Sequence[tuple[str, str]],
    expression_outputs: Sequence[tuple[str, str, tuple[str, ...]]],
    artifact_by_receipt: Mapping[str, str],
    inherited_rejected_artifacts: Sequence[str] = (),
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Delivered quantities standing on a result its own verdict rejected.

    What a verdict rejects is a **result**, and a result is not a
    receipt: one finished calculation is routinely read by several
    extraction calls, each with its own receipt.  Keying the rejection
    on the receipt the failed rule happened to read lets every other
    read of the same result escape -- one live run extracted
    frequencies and coordinates from one saddle in two calls, and the
    torsion of that rejected structure, which was the task's own
    requested observable, was reported as delivered beside the
    zero-point energy that was correctly withheld.  So the seed
    resolves to artifact digests and taints every receipt that read
    them.

    A rule that reads an expression output rather than an extraction
    resolves backwards through that output's own sources, because the
    verdict is about the structure underneath, not about the
    arithmetic: one run computed a complexation energy and a barrier in
    a single expression, and only the barrier stood on the rejected
    result.

    Returns the stale quantity ids and the rejected artifact digests,
    the latter so a goal can carry them across cycles -- a rejection is
    a fact about bytes and does not expire.
    """

    sources_by_output = {
        (receipt, output_id): sources
        for receipt, output_id, sources in expression_outputs
    }
    rejected_artifacts = {
        str(item) for item in inherited_rejected_artifacts if item
    }

    def _resolve(receipt: str, quantity_id: str, depth: int = 0) -> None:
        """Name the results a rejected binding ultimately rests on."""

        if depth > 8:
            return
        artifact = artifact_by_receipt.get(receipt)
        if artifact:
            rejected_artifacts.add(artifact)
            return
        for source in sources_by_output.get((receipt, quantity_id), ()):
            if not source:
                continue
            producer = artifact_by_receipt.get(source)
            if producer:
                rejected_artifacts.add(producer)
                continue
            for other, output_id in sources_by_output:
                if other == source:
                    _resolve(other, output_id, depth + 1)

    for receipt, quantity_id in rejected_bindings:
        if receipt:
            _resolve(receipt, quantity_id)

    # Every read of a rejected result, not only the one the rule saw.
    tainted_receipts = {
        receipt
        for receipt, artifact in artifact_by_receipt.items()
        if artifact in rejected_artifacts
    }
    seeded = set(tainted_receipts)
    tainted_outputs: set[tuple[str, str]] = set()
    while True:
        grew = False
        for receipt, output_id, sources in expression_outputs:
            if (receipt, output_id) in tainted_outputs:
                continue
            if not tainted_receipts.intersection(sources):
                continue
            tainted_outputs.add((receipt, output_id))
            tainted_receipts.add(receipt)
            grew = True
        if not grew:
            break
    stale = tuple(
        sorted(
            {
                quantity_id
                for receipt, quantity_id in claim_pairs
                if quantity_id
                and (
                    receipt in seeded
                    or (receipt, quantity_id) in tainted_outputs
                )
            }
        )
    )
    return stale, tuple(sorted(rejected_artifacts))


#: Job types whose result promises a stationary point; a run of one that
#: printed no vibrational modes never checked that promise.
_STATIONARY_POINT_JOBTYPES = frozenset({"opt", "ts"})


def _printed_no_modes(record: Mapping[str, Any]) -> bool:
    """Whether a verified opt/ts result carries no frequency block.

    Read from the validator's own observations: every program's block
    records ``vibrational_mode_count``, and a count of zero on a job
    type that promises a stationary point means the promise was never
    checked, for a valid result as much as a failed one.
    """

    observations = record.get("observations") or {}
    jobtype = str(observations.get("jobtype") or record.get("jobtype") or "")
    if jobtype not in _STATIONARY_POINT_JOBTYPES:
        return False
    for key, value in observations.items():
        if isinstance(value, Mapping) and "vibrational_mode_count" in value:
            try:
                return int(value.get("vibrational_mode_count") or 0) == 0
            except (TypeError, ValueError):
                return False
    return False


def _analysis_delivery(
    events_path: Path,
    *,
    goal_delivered_ids: Mapping[str, Mapping[str, Any]] | None = None,
    declared_observables: Sequence[Mapping[str, Any]] = (),
    uncharacterised_artifact_sha256s: tuple[str, ...] = (),
    inherited_rejected_artifacts: Sequence[str] = (),
    flagged_artifact_sha256s: Sequence[str] = (),
    failed_artifact_sha256s: Sequence[str] = (),
    inherited_unreachable: Mapping[str, str] = {},
) -> _AnalysisDelivery:
    """Read the delivery facts a settlement stands on.

    Every field is a typed record the host itself wrote: the
    completion receipt with its stated limitations, the claim and
    decision records, and the receipt digests a settlement cites. The
    first live goal round's classifier read none of these -- it
    watched validation rules, the one place an honest refusal leaves
    no trace -- and settled a receipts-backed refusal as achieved.
    """

    completion_status = ""
    prediction_rows: tuple[dict[str, Any], ...] = ()
    route_dispositions: tuple[dict[str, Any], ...] = ()

    declared_misses: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    failed_verdicts: list[tuple[str, str, str]] = []
    decision_refs: set[str] = set()
    claims = 0
    decisions = 0
    decision_uncertainties: list[str] = []
    sufficiency_rows: list[dict[str, Any]] = []
    verified_unreachable: set[str] = set()
    unverified_unreachable: set[str] = set()
    # A refusal the host verified is recorded by the planning session
    # that wrote it, and an executed run's own stream never carries a
    # scientific decision -- the provider-free walker records
    # extraction, thermochemistry, expressions, verdicts and claims and
    # nothing else. So the assessment the run wrote and the refusal that
    # answers it live in two streams, and neither delivery could see the
    # other: a requirement the session refused with receipts re-opened
    # the goal and spent a recovery cycle on an obligation the host had
    # already certified as refused. One delivery inherits the other's
    # verified refusals explicitly, because two projections cannot each
    # subtract the other's.
    #
    # Only the verified ones, and the verification travels with the id
    # rather than beside it. The first version of this inheritance
    # passed the *bases* mapping, which explains every refusal the
    # session wrote including the ones the host declined to verify, and
    # the receiver promoted every key it received: an unverified refusal
    # arrived carrying authority it had been denied, removed the run's
    # open requirement, and left a settlement that could reach achieved
    # -- the previous audit's own pattern, inside its repair. An
    # explanation string never confers authority.
    unreachable_bases: dict[str, str] = {
        observable_id: str(basis)
        for observable_id, basis in inherited_unreachable.items()
    }
    verified_unreachable |= set(unreachable_bases)
    receipts: list[str] = []
    doubt_refs: set[str] = set()
    claim_pairs: list[tuple[str, str]] = []
    claim_rows: dict[str, dict[str, Any]] = {}
    rejected_bindings: list[tuple[str, str]] = []
    expression_outputs: list[tuple[str, str, tuple[str, ...]]] = []
    artifact_by_receipt: dict[str, str] = {}
    failed_artifacts: set[str] = set()
    uncharacterised_artifacts: set[str] = set(uncharacterised_artifact_sha256s)
    characterised: set[str] = set()
    # The handoff edges and what each node left, so a validated Hessian
    # characterises the optimisation whose reached geometry it consumed.
    handoffs: dict[str, str] = {}
    node_outputs: dict[str, tuple[str, ...]] = {}
    characterising: set[str] = set()
    stopped_by: list[str] = []
    workflows_planned = 0
    nodes_previewed = 0
    last_plan_refusal = ""
    terminal_reason = ""
    anomaly_ids: tuple[str, ...] = ()
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        lines = []
    for line in lines:
        text = line.strip()
        if not text:
            continue
        try:
            event = json.loads(text)
        except json.JSONDecodeError:
            continue
        kind = str(event.get("kind") or "")
        payload = event.get("payload") or {}
        digest = str(payload.get("receipt_sha256") or "")
        if kind == "scientific_decision_recorded":
            decisions += 1
            record = payload.get("record") or {}
            route_dispositions += tuple(
                dict(item)
                for item in payload.get("menu_route_dispositions") or ()
                if isinstance(item, Mapping)
            )
            for item in payload.get("unreachable_observables") or ():
                observable_id = str(item.get("observable_id") or "")
                if not observable_id:
                    continue
                unreachable_bases[observable_id] = (
                    f"{item.get('statement') or ''} [{item.get('basis') or ''}]"
                )
                if bool(item.get("verified")):
                    verified_unreachable.add(observable_id)
                else:
                    unverified_unreachable.add(observable_id)
                # The receipts a refusal cites are the settlement's
                # evidence when nothing else was claimed.
                receipts.extend(
                    str(digest)
                    for digest in item.get("receipt_sha256s") or ()
                    if digest
                )
            # What the session said it was unsure of, carried where the
            # human reads first. A woken session wrote that the author's
            # 2 kJ/mol design threshold sits inside the method's 1-2
            # kJ/mol error -- the most useful sentence either window
            # produced -- and it reached no receipt, no claim and no
            # settlement (NOVEL-2 po2, 2026-09-04).
            decision_uncertainties.extend(
                str(item)
                for item in record.get("uncertainties") or ()
                if str(item).strip()
            )
            for reference in record.get("evidence_refs") or ():
                text_ref = str(reference)
                if text_ref.startswith("doubt:"):
                    doubt_refs.add(text_ref[len("doubt:") :])
                # A decision that cites the validation receipt has looked
                # at the failed verdict and stood by its delivery. That is
                # the scientist's call to make, and citing the receipt is
                # how it is made without the host grading prose.
                decision_refs.add(text_ref.split(":")[-1])
        elif kind == "analysis_claims_recorded":
            claims += 1
            if digest:
                receipts.append(digest)
            # How each delivered number stands against the precision its
            # own declaration asked for, judged by the one function the
            # gate and the settlement share.
            for row in payload.get("sufficiency") or ():
                if isinstance(row, Mapping) and row.get("observable_id"):
                    sufficiency_rows.append(dict(row))
            record = payload.get("record") or {}
            for claim in record.get("claims") or ():
                receipt_digest = str(claim.get("source_receipt_sha256") or "")
                claim_pairs.append(
                    (receipt_digest, str(claim.get("quantity_id") or ""))
                )
                # What a stream claim carries, under both the names it
                # answers to: a claim of this cycle used to be counted
                # delivered on its id alone while the settlement checked
                # the dimension, so a wrong-dimension claim suppressed
                # the very wake that could have repaired it. The row is
                # kept whole -- a claim written before dimensions
                # travelled carries only its display unit, and the
                # shared predicate resolves that through the same unit
                # table the analysis layer uses.
                row = {
                    "dimension": tuple(
                        int(item) for item in claim.get("dimension") or ()
                    ),
                    "display_unit": str(claim.get("display_unit") or ""),
                }
                for name in (
                    str(claim.get("claim_id") or ""),
                    str(claim.get("quantity_id") or ""),
                ):
                    if name:
                        # Last-wins, with the gate and the record: a
                        # corrected claim is the one that answers.
                        claim_rows[name] = row
                # The id a claim carries is delivered under both names
                # it holds. A wake listed six ids as delivered (by
                # quantity_id) and the same six as undelivered (by
                # claim_id) in one record (NOVEL-2 po2, 2026-09-04); the
                # gate now joins on either field, and so does this.
                claim_id = str(claim.get("claim_id") or "")
                if claim_id and claim_id != claim.get("quantity_id"):
                    claim_pairs.append((receipt_digest, claim_id))
        elif kind == "workflow_node_launch_refused":
            stopped_by.append(
                f"node {payload.get('node_id')} never launched: "
                f"{payload.get('reason')}"
            )
        elif kind == "execution_review_refused":
            stopped_by.append(
                f"execution review refused: {payload.get('reason')}"
            )
        elif kind == "command_workflow_planned":
            workflows_planned += 1
        elif kind == "program_node_preflighted":
            nodes_previewed += 1
        elif kind == "tool_failed" and (
            str(payload.get("tool") or "") == "plan_scientific_workflow"
        ):
            report = payload.get("failure_report") or {}
            last_plan_refusal = str(
                report.get("gate")
                or (payload.get("canonical_result") or {}).get("message")
                or payload.get("message")
                or "refused"
            )[:240]
        elif kind == "runtime_terminated":
            terminal_reason = str(
                payload.get("reason")
                or payload.get("terminal_reason")
                or payload.get("terminal_state")
                or ""
            )[:240]
        elif kind == "program_result_verified":
            record = payload.get("record") or {}
            node_name = str(
                record.get("node_id") or payload.get("node_id") or ""
            )
            if node_name:
                node_outputs[node_name] = tuple(
                    str(item.get("sha256") or "")
                    for item in record.get("output_artifacts") or ()
                    if item.get("sha256")
                )
                if str(record.get("state") or "") == "valid" and (
                    printed_modes(record)
                ):
                    characterising.add(node_name)
            if str(record.get("state") or "") != "valid":
                failed_artifacts.update(
                    str(item.get("sha256") or "")
                    for item in record.get("output_artifacts") or ()
                    if item.get("sha256")
                )
            # An optimisation that printed no frequencies made no claim
            # about its stationary point, and a number standing on it
            # says so. NOVEL-1's iron runs were opt-only; the trans
            # triplet delivered as an ordered gap was, under NOVEL-2's
            # Freq on the bit-identical geometry, a second-order saddle
            # (owner ruling 2026-09-05: a provenance word, never a
            # refusal).
            if _printed_no_modes(record):
                uncharacterised_artifacts.update(
                    str(item.get("sha256") or "")
                    for item in record.get("output_artifacts") or ()
                    if item.get("sha256")
                )
        elif kind == "stationary_point_characterised":
            record = payload.get("record") or {}
            digest = str(record.get("result_artifact_sha256") or "")
            if digest:
                characterised.add(digest)
        elif kind == "optimized_geometry_handed_off":
            if str(payload.get("status") or "") == "validated_handoff":
                handoffs[str(payload.get("consumer_node_id") or "")] = str(
                    payload.get("producer_node_id") or ""
                )
        elif kind == "analysis_completion_evaluated":
            completion_status = str(payload.get("status") or "")
            prediction_rows = tuple(
                dict(row)
                for row in (
                    payload.get("declared_observable_predictions") or ()
                )
                if isinstance(row, Mapping)
            )
            limitations = tuple(
                str(item)
                for item in (payload.get("limitation_output_ids") or ())
            )
            declared_misses = tuple(
                str(item)
                for item in (payload.get("declared_observable_misses") or ())
            )
            anomaly_ids = tuple(
                str(item) for item in (payload.get("anomaly_output_ids") or ())
            )
            if digest:
                receipts.append(digest)
        elif kind == "scientific_validation_evaluated":
            if digest:
                receipts.append(digest)
            if not bool(payload.get("all_rules_passed", True)):
                node_id = str(payload.get("node_id") or "")
                record = payload.get("record") or {}
                bindings = {
                    str(binding.get("input_id") or ""): (
                        str(binding.get("source_receipt_sha256") or ""),
                        str(binding.get("quantity_id") or ""),
                    )
                    for binding in record.get("input_bindings") or ()
                }
                for rule in record.get("rule_results") or ():
                    if not bool(rule.get("passed", True)):
                        failed_verdicts.append(
                            (
                                node_id,
                                str(rule.get("rule_id") or ""),
                                digest,
                            )
                        )
                        # The rule read these receipts and rejected what
                        # it found in them; everything else computed from
                        # the same receipts describes the same rejected
                        # structure.
                        for input_id in rule.get("input_ids") or ():
                            binding = bindings.get(str(input_id))
                            if binding and binding[0]:
                                rejected_bindings.append(binding)
        elif kind in {
            "result_quantities_extracted",
            "quantity_expression_evaluated",
            "thermochemistry_derived",
        }:
            if digest:
                receipts.append(digest)
            if kind == "result_quantities_extracted" and digest:
                # The result this receipt read. A verdict rejects the
                # result, and one result is routinely read by several
                # extraction calls.
                record = payload.get("record") or {}
                artifact = str(
                    payload.get("artifact_sha256")
                    or record.get("artifact_sha256")
                    or ""
                )
                if artifact:
                    artifact_by_receipt[digest] = artifact
            if kind == "quantity_expression_evaluated" and digest:
                record = payload.get("record") or {}
                for dependency in record.get("output_dependencies") or ():
                    expression_outputs.append(
                        (
                            digest,
                            str(dependency.get("output_id") or ""),
                            tuple(
                                str(item)
                                for item in (
                                    dependency.get("source_receipt_sha256s")
                                    or ()
                                )
                            ),
                        )
                    )
    unanswered = tuple(
        f"{node_id}/{rule_id}"
        for node_id, rule_id, digest in failed_verdicts
        if digest not in decision_refs
    )
    # Quantities standing on a node that did not meet the promise it was
    # launched under, split by whether the session had the host check what
    # that structure is. Both are statements, never refusals: the number
    # may be exactly the finding worth reporting.
    # A two-node minimum: the Hessian that consumed an optimisation's
    # reached geometry through a validated handoff characterised that
    # geometry, byte for byte, and the number read from the optimisation
    # is not "uncharacterised (no frequencies printed)" -- it was worded
    # so on two live PySCF goals while the Hessian beside it validated
    # (PySCF round, 2026-09-12).
    for consumer, producer in handoffs.items():
        if consumer in characterising:
            characterised.update(node_outputs.get(producer, ()))
    failed_seed = set(failed_artifacts) | set(failed_artifact_sha256s)
    failed_quantities, _failed_walk = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=(),
        expression_outputs=expression_outputs,
        artifact_by_receipt=artifact_by_receipt,
        inherited_rejected_artifacts=tuple(sorted(failed_seed)),
    )
    characterised_quantities, _checked_walk = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=(),
        expression_outputs=expression_outputs,
        artifact_by_receipt=artifact_by_receipt,
        inherited_rejected_artifacts=tuple(
            sorted(failed_seed & characterised)
        ),
    )
    uncharacterised_quantities, _unchecked_walk = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=(),
        expression_outputs=expression_outputs,
        artifact_by_receipt=artifact_by_receipt,
        inherited_rejected_artifacts=tuple(
            sorted(uncharacterised_artifacts - characterised)
        ),
    )
    stale, rejected_artifacts = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=rejected_bindings,
        expression_outputs=expression_outputs,
        artifact_by_receipt=artifact_by_receipt,
        inherited_rejected_artifacts=inherited_rejected_artifacts,
    )
    # The same walk, seeded with the artifacts of every node an anomaly
    # flagged: which delivered numbers stand on a flagged result. Not a
    # verdict and not a rejection -- the flagged run may be right.
    flagged, _flagged_artifacts = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=(),
        expression_outputs=expression_outputs,
        artifact_by_receipt=artifact_by_receipt,
        inherited_rejected_artifacts=flagged_artifact_sha256s,
    )
    # An expression's exported outputs, as opposed to the intermediate
    # node_values it computed on the way: the receipt contract pins
    # output_dependencies' ids to outputs' quantity_ids, in order, so the
    # lineage already collected names exactly what was exported.
    claimed_ids = {
        quantity_id for _receipt, quantity_id in claim_pairs if quantity_id
    }
    exported_output_ids = {
        output_id
        for _digest, output_id, _sources in expression_outputs
        if output_id
    }
    if stopped_by:
        ending = "; ".join(stopped_by)
    elif workflows_planned == 0:
        ending = "no workflow was planned" + (
            f"; the last plan call was refused ({last_plan_refusal})"
            if last_plan_refusal
            else ""
        )
    else:
        ending = (
            f"{workflows_planned} workflow(s) planned, {nodes_previewed} "
            "node preview(s), no execution review built"
        )
    if terminal_reason:
        ending += f"; the session ended: {terminal_reason}"
    # An assessment recorded in an earlier cycle reaches this one
    # through the workspace record, and this stream's own rows come
    # after it so the current assessment wins.
    carried = tuple(
        row["sufficiency"]
        for row in (goal_delivered_ids or {}).values()
        if isinstance(row, Mapping)
        and isinstance(row.get("sufficiency"), Mapping)
    )
    return _AnalysisDelivery(
        ending=ending,
        terminal_reason=terminal_reason,
        sufficiency=carried + tuple(sufficiency_rows),
        unclaimed_output_ids=tuple(sorted(exported_output_ids - claimed_ids)),
        stopped_by=tuple(stopped_by),
        anomaly_output_ids=anomaly_ids,
        unanswered_verdicts=unanswered,
        completion_status=completion_status,
        limitation_output_ids=limitations,
        claims=claims,
        decisions=decisions,
        receipt_sha256s=tuple(receipts),
        doubted_quantity_ids=tuple(
            sorted(
                {
                    quantity_id
                    for receipt, quantity_id in claim_pairs
                    if quantity_id and receipt in doubt_refs
                }
            )
        ),
        claim_rows=dict(claim_rows),
        delivered_quantity_ids=tuple(
            sorted(
                {
                    quantity_id
                    for _receipt, quantity_id in claim_pairs
                    if quantity_id and quantity_id not in stale
                }
            )
        ),
        stale_quantity_ids=stale,
        flagged_quantity_ids=flagged,
        rejected_artifact_sha256s=rejected_artifacts,
        claims_rendered=bool(claims),
        failed_source_quantity_ids=failed_quantities,
        characterised_source_quantity_ids=characterised_quantities,
        characterised_artifact_sha256s=tuple(sorted(characterised)),
        uncharacterised_source_quantity_ids=uncharacterised_quantities,
        decision_uncertainties=tuple(decision_uncertainties),
        declared_observable_misses=declared_misses,
        goal_delivered=dict(goal_delivered_ids or {}),
        declared_observables=tuple(
            dict(item) for item in declared_observables
        ),
        verified_unreachable_ids=tuple(sorted(verified_unreachable)),
        prediction_rows=prediction_rows,
        route_dispositions=route_dispositions,
        unverified_unreachable_ids=tuple(
            sorted(unverified_unreachable - verified_unreachable)
        ),
        unreachable_bases=dict(unreachable_bases),
    )


def _settlement_evidence(delivery: _AnalysisDelivery) -> dict[str, Any]:
    """Receipts a settlement cites, from the session's own stream."""

    if delivery.decisions and delivery.receipt_sha256s:
        evidence: dict[str, Any] = {
            "scientific_decisions": delivery.decisions,
            "receipt_sha256s": delivery.receipt_sha256s,
        }
        if delivery.decision_uncertainties:
            evidence["decision_uncertainties"] = (
                delivery.decision_uncertainties
            )
        return evidence
    return {}


#: The goal's task text, kept beside its ledger so a resumed driver plans
#: the next cycle against exactly what the human asked.
TASK_FILE = "task.md"

#: How the driver was constructed -- envelope, grant, provider, dispatch
#: -- so a process that wakes the goal needs only the workspace and id.
DRIVER_FILE = "driver.json"

#: Where an approved run is executed: in this process, or handed to the
#: server profile's scheduler with the job waking the goal when done.
DISPATCH_MODES = ("local", "scheduler")

#: The phases one goal cycle passes through, in order. ``parked`` and
#: ``settled`` are where a process may stop: a parked goal has a run in
#: a scheduler's hands and resumes at ``outcome``; a settled goal is
#: finished.
GOAL_PHASES = (
    "plan",
    "decide",
    "execute",
    "outcome",
    "settle",
    "parked",
    "settled",
)


@dataclass(frozen=True)
class GoalStepV1:
    """One phase the driver just performed."""

    phase: str
    next_phase: str
    cycle: int
    result: GoalLoopResultV1 | None = None


class GoalDriver:
    """Drive one goal to settlement under one human decision, one step at
    a time.

    ``initial_decision`` is the human's: "approve" records it under
    ``granted_by`` exactly as ``agent review --decision approve`` would,
    and anything else settles the goal ``returned_to_human`` before any
    engine runs. Every later cycle's decision is host-admitted under the
    goal and recorded with the composite actor.

    ``dispatch_run``, when given, replaces the in-process engine walk:
    it hands the approved bundle and run directory to a scheduler and
    returns a receipt naming the job; the driver records
    ``run_dispatched`` and parks. :meth:`resume` rebuilds a driver from
    the ledger at the ``outcome`` phase once the job has finished.
    """

    def __init__(
        self,
        *,
        task: str,
        workspace: str | Path,
        execution_envelope_file: str | Path | None,
        goal_id: str,
        granted_by: str,
        max_revisions: int = 5,
        provider: str | None = None,
        provider_config_file: str | Path | None = None,
        analysis_completion_file: str | Path | None = None,
        plan_session: Callable[..., Any] = _default_plan_session,
        resolve_review: Callable[..., tuple[str, Path]] = _default_resolve,
        execute_bundle: Callable[..., Any] = _default_execute,
        dispatch_run: Callable[..., Any] | None = None,
        dispatch: str = "local",
        server: str | None = None,
        initial_decision: str = "approve",
        stop_file: str | Path | None = None,
        session_kwargs: Mapping[str, Any] | None = None,
        _resuming: bool = False,
    ) -> None:
        if dispatch not in DISPATCH_MODES:
            raise ContractError(
                f"unsupported dispatch mode {dispatch!r}; "
                f"choose one of {DISPATCH_MODES}"
            )
        self.task = task
        self.workspace = Path(workspace).resolve()
        self.goal_id = goal_id
        self.granted_by = granted_by
        self.max_revisions = max_revisions
        self.provider = provider
        self.provider_config_file = provider_config_file
        self.execution_envelope_file = execution_envelope_file
        self.analysis_completion_file = analysis_completion_file
        self.plan_session = plan_session
        self.resolve_review = resolve_review
        self.execute_bundle = execute_bundle
        self.dispatch = dispatch
        self.server = server
        if dispatch_run is None and dispatch == "scheduler":
            from functools import partial

            from chemsmart.agent.dispatch import dispatch_run_to_scheduler

            dispatch_run = partial(dispatch_run_to_scheduler, server=server)
        self.dispatch_run = dispatch_run
        self.initial_decision = initial_decision
        self.stop_file = stop_file
        self.session_kwargs = dict(session_kwargs or {})

        self.goal_dir = self.workspace / ".chemsmart-agent" / "goals" / goal_id
        self.ledger = GoalLedger(self.goal_dir)
        self.ledger.directory.mkdir(parents=True, exist_ok=True)
        if not _resuming and self.ledger.goal_path.exists():
            raise ContractError(
                "this goal already exists; a goal is one human decision, "
                "not a resumable queue. A parked goal is resumed with "
                "GoalDriver.resume (chemsmart agent wake)."
            )
        if not _resuming:
            # The task text is what every later cycle plans against; a
            # process that resumes a parked goal reads it from here, and
            # the construction record beside it says how to drive on.
            (self.goal_dir / TASK_FILE).write_text(str(task), encoding="utf-8")
            (self.goal_dir / DRIVER_FILE).write_text(
                json.dumps(
                    {
                        "schema_version": "chemsmart.goal-driver.v1",
                        "execution_envelope_file": _resolved_or_none(
                            execution_envelope_file
                        ),
                        "granted_by": granted_by,
                        "max_revisions": int(max_revisions),
                        "provider": provider,
                        "provider_config_file": _resolved_or_none(
                            provider_config_file
                        ),
                        "analysis_completion_file": _resolved_or_none(
                            analysis_completion_file
                        ),
                        "dispatch": dispatch,
                        "server": server,
                        "stop_file": _resolved_or_none(stop_file),
                    },
                    indent=2,
                    sort_keys=True,
                ),
                encoding="utf-8",
            )

        self.envelope = None
        self.envelope_record: dict[str, Any] = {}
        if execution_envelope_file is not None:
            from chemsmart.agent.execution_envelope import (
                load_bounded_execution_envelope,
            )

            self.envelope = load_bounded_execution_envelope(
                execution_envelope_file
            )
            self.envelope_record = _goal_envelope_record(
                {
                    "allowed_program_engines": (
                        self.envelope.allowed_program_engines
                    ),
                    "max_engine_calls": self.envelope.max_engine_calls,
                    "episode_wall_time_seconds": (
                        self.envelope.episode_wall_time_seconds
                    ),
                    "max_excursion_calls": self.envelope.max_excursion_calls,
                }
            )

        self.goal: GoalRecordV1 | None = None
        self.outcome: Any = None
        #: The failure report a re-woken cycle carries: what the previous
        #: cycle left undone, said as gate, diagnosis, route and cost.
        self.failure_report: dict[str, Any] | None = None
        self.cycles = 0
        self.revisions_admitted = 0
        # A rejection is a fact about bytes and does not expire, so the
        # rejected results accumulate; the *standing* delivery is
        # whatever the most recent claim-rendering cycle said, because
        # that is what the goal currently answers with. Keying this on
        # quantity ids instead would ask a later cycle to reuse an
        # earlier cycle's names: one live recovery re-derived a torsion
        # correctly under a new id and would have been held open forever
        # over a number it had already replaced.
        self.rejected_artifacts: set[str] = set()
        #: (cycle, stream) pairs already projected into the workspace
        #: record, because record_run appends and never deduplicates.
        self._recorded_streams: set[tuple[int, str]] = set()
        self.standing_stale: tuple[str, ...] = ()

        self.phase = "plan"
        self.result: GoalLoopResultV1 | None = None
        self.wake: dict[str, Any] | None = None
        self.session: Any = None
        self.events_path: Path | None = None
        self.review_file: Path | None = None
        self.bundle_file: Path | None = None
        self.run_directory: Path | None = None
        self.execute_result: Any = None
        self._pending_declarations: tuple[dict[str, Any], ...] = ()
        #: Ledger rows a first cycle produced before its goal record
        #: existed, in the order they were produced. Cycle 1 runs every
        #: cycle-ending recorder while ``self.goal`` is still None, and
        #: four of the five used to return there instead of deferring --
        #: only declarations had learned this. Round 14 paid for it
        #: twice: ino3-r17 lost six aborted ORCA input checks, so the
        #: next wake read ``aborted: 0`` and re-bought the diagnosis for
        #: six of forty engine calls; po3-r17 lost the evidence
        #: reference and its re-wake could never be admitted, ending
        #: with no calculation at all.
        self._pending_ledger_rows: tuple[tuple[str, dict[str, Any]], ...] = ()
        self.dispatch_receipt: Any = None

    # -- public surface ---------------------------------------------------

    @classmethod
    def resume(
        cls, *, workspace: str | Path, goal_id: str, **kwargs: Any
    ) -> "GoalDriver":
        """Rebuild a parked goal's driver from its ledger, at ``outcome``.

        Admitted only when the ledger's last run was dispatched and never
        recorded, and the goal is not settled: the same one human
        decision continues in its own run directory, and nothing here
        creates a second one.
        """

        goal_dir = (
            Path(workspace).resolve() / ".chemsmart-agent" / "goals" / goal_id
        )
        try:
            recorded = json.loads(
                (goal_dir / DRIVER_FILE).read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError):
            recorded = {}
        for key in (
            "execution_envelope_file",
            "granted_by",
            "max_revisions",
            "provider",
            "provider_config_file",
            "analysis_completion_file",
            "dispatch",
            "server",
            "stop_file",
        ):
            if key not in kwargs and recorded.get(key) is not None:
                kwargs[key] = recorded[key]
        if "granted_by" not in kwargs:
            raise ContractError(
                f"goal {goal_id!r} has no recorded grant to resume under"
            )
        driver = cls(
            task=str(kwargs.pop("task", "")),
            workspace=workspace,
            goal_id=goal_id,
            _resuming=True,
            **kwargs,
        )
        if not driver.ledger.goal_path.exists():
            raise ContractError(f"goal {goal_id!r} has no ledger to resume")
        entries = driver.ledger.entries()
        if any(entry["kind"] == "goal_settled" for entry in entries):
            raise ContractError(f"goal {goal_id!r} is settled")
        dispatched = [
            entry for entry in entries if entry["kind"] == "run_dispatched"
        ]
        recorded_cycles = {
            int(entry["payload"].get("cycle", 0))
            for entry in entries
            if entry["kind"] == "run_recorded"
        }
        parked = [
            entry
            for entry in dispatched
            if int(entry["payload"].get("cycle", 0)) not in recorded_cycles
        ]
        dispatched_cycles = {
            int(entry["payload"].get("cycle", 0)) for entry in dispatched
        }
        interrupted = [
            entry
            for entry in entries
            if entry["kind"] == "run_started"
            and int(entry["payload"].get("cycle", 0)) not in recorded_cycles
            and int(entry["payload"].get("cycle", 0)) not in dispatched_cycles
        ]
        if not parked and not interrupted:
            raise ContractError(
                f"goal {goal_id!r} has no parked or interrupted run to resume"
            )
        entry = (parked or interrupted)[-1]
        driver.goal = driver.ledger.load()
        if not driver.task:
            try:
                driver.task = (driver.goal_dir / TASK_FILE).read_text(
                    encoding="utf-8"
                )
            except OSError as exc:
                raise ContractError(
                    f"goal {goal_id!r} has no recorded task text to resume"
                ) from exc
        driver.cycles = int(entry["payload"]["cycle"])
        driver.revisions_admitted = sum(
            1 for item in entries if item["kind"] == "revision_admitted"
        )
        # Replay the standing delivery from every earlier recorded run, in
        # order, so a rejection made in cycle 1 still taints cycle 3.
        for item in entries:
            if item["kind"] != "run_recorded":
                continue
            delivery = _analysis_delivery(
                driver.workspace
                / ".chemsmart-agent"
                / Path(*str(item["payload"].get("run") or "").split("/"))
                / "events.jsonl",
                inherited_rejected_artifacts=tuple(
                    sorted(driver.rejected_artifacts)
                ),
            )
            driver.rejected_artifacts.update(
                delivery.rejected_artifact_sha256s
            )
            if delivery.claims_rendered:
                driver.standing_stale = delivery.stale_quantity_ids
        if entry["kind"] == "run_started":
            # An interrupted local run: re-enter the execute phase with the
            # same bundle and run directory; the executor's continuation
            # replays finished nodes and names the interrupted one.
            driver.bundle_file = Path(str(entry["payload"]["approval_file"]))
            driver.phase = "execute"
            return driver
        driver.run_directory = (
            driver.goal_dir / "runs" / f"cycle-{driver.cycles}"
        )
        driver.dispatch_receipt = dict(entry["payload"])
        driver.execute_result = _recorded_execution_result(
            driver.run_directory
        )
        driver.phase = "outcome"
        return driver

    @classmethod
    def from_review(
        cls,
        *,
        review_file: str | Path,
        task_spec_sha256: str = "",
        **kwargs: Any,
    ) -> "GoalDriver":
        """A driver at ``decide`` over a review a session already wrote.

        This is how a stored review is re-presented for one fresh human
        decision -- the terminal interface's resume, `agent review` --
        without planning again: the review file stands in for cycle 1's
        session, and everything from the decision on is the same
        machine.
        """

        driver = cls(**kwargs)
        driver.cycles = 1
        driver.review_file = Path(review_file)
        driver.session = SimpleNamespace(
            terminal_state="waiting_for_approval",
            task_spec_sha256=str(task_spec_sha256 or ""),
        )
        driver.phase = "decide"
        return driver

    def run(self) -> GoalLoopResultV1:
        """Step until the goal settles or parks."""

        while self.phase not in {"settled", "parked"}:
            self.step()
        assert self.result is not None
        return self.result

    def step(self) -> GoalStepV1:
        """Perform exactly one phase and say which comes next."""

        phase = self.phase
        if phase in {"settled", "parked"}:
            raise ContractError(f"goal {self.goal_id!r} is {phase}")
        handler = {
            "plan": self._plan,
            "decide": self._decide,
            "execute": self._execute,
            "outcome": self._outcome,
            "settle": self._settle,
        }[phase]
        handler()
        return GoalStepV1(
            phase=phase,
            next_phase=self.phase,
            cycle=self.cycles,
            result=self.result,
        )

    # -- phases -------------------------------------------------------------

    def _stopped(self) -> bool:
        return bool(self.stop_file and Path(self.stop_file).exists())

    def _settled(
        self, settlement: str, reasons: tuple[str, ...]
    ) -> GoalLoopResultV1:
        self.result = GoalLoopResultV1(
            goal_id=self.goal_id,
            settlement=settlement,
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            reasons=reasons,
        )
        self.phase = "settled"
        return self.result

    def _defer_or_append(self, kind: str, payload: dict[str, Any]) -> None:
        """Enter a cycle-ending row, or hold it until the goal exists.

        A first cycle produces these before ``goal_created``, and a row
        that is dropped there is gone: the ledger is what every later
        process reads. So the row waits rather than returning, and the
        flush that follows ``ledger.create`` enters it in the order it
        was produced. Nothing new is asked of the model.
        """

        if self.goal is None:
            self._pending_ledger_rows += ((kind, payload),)
            return
        self.ledger.append(kind, payload)

    def _record_declarations(self) -> None:
        """Enter this cycle's new observable declarations in the ledger."""

        known = {
            str(record.get("observable_id") or "")
            for record in _first_declarations(self.ledger)
        }
        fresh = tuple(
            record
            for record in _session_declarations(self.events_path)
            if str(record.get("observable_id") or "") not in known
        )
        if self.goal is None:
            # Cycle 1 declares before the goal record exists; the entry
            # follows goal_created.
            self._pending_declarations = fresh
            return
        if fresh:
            self.ledger.append(
                "observables_declared",
                {"cycle": self.cycles, "observables": fresh},
            )

    def _record_approaches(self) -> None:
        """Enter what this cycle tried, and how it ended, in the ledger.

        The wake carried budgets, anomalies and a repair menu, and not
        what had already been attempted: po3 re-seeded an approach its
        own previous cycle had diagnosed as the failure (REACH-1,
        2026-09-06). Nothing new is asked of the model; this is the
        stream's own words, kept where the next cycle reads.
        """

        approaches = _session_approaches(self.events_path) + _run_approaches(
            (self.run_directory / "events.jsonl")
            if self.run_directory is not None
            else None
        )
        if approaches:
            self._defer_or_append(
                "approaches_recorded",
                {"cycle": self.cycles, "approaches": approaches},
            )

    def _record_dispositions(self) -> None:
        """Enter this cycle's repair-menu dispositions in the ledger, so
        the next wake shows what was done with the menu it inherits."""

        dispositions = _session_dispositions(self.events_path)
        if dispositions:
            self._defer_or_append(
                "repair_menu_dispositions",
                {"cycle": self.cycles, "dispositions": dispositions},
            )

    def _record_input_checks(self) -> None:
        """Enter this cycle's input-check probe counts in the ledger."""

        counts = _session_input_checks(self.events_path)
        if any(counts.values()):
            self._defer_or_append(
                "input_checks_probed", {"cycle": self.cycles, **counts}
            )

    def _flush_declarations(self) -> None:
        """Enter every row a first cycle produced before its goal record.

        Declarations first, because the goal's contract is unreadable
        without them and every later join is by declared id; then the
        rest in the order they were produced.
        """

        fresh, self._pending_declarations = self._pending_declarations, ()
        if fresh:
            self.ledger.append(
                "observables_declared",
                {"cycle": self.cycles, "observables": fresh},
            )
        deferred, self._pending_ledger_rows = self._pending_ledger_rows, ()
        for kind, payload in deferred:
            self.ledger.append(kind, payload)

    def _declared_non_executable_ids(self) -> tuple[str, ...]:
        """Node ids the cycle's displayed review retained as intent only."""

        review_file = self.goal_dir / "reviews" / f"cycle-{self.cycles}.json"
        try:
            review = _review_record(review_file)
        except (OSError, json.JSONDecodeError):
            return ()
        packet = review.get("workflow_execution_review") or review
        return tuple(
            str(item) for item in packet.get("non_executable_node_ids") or ()
        )

    def _record_workspace(self, events_path: Path, run_reference: str) -> None:
        """Append what this run proved to the workspace's own record.

        Host-written from receipts; never ends a goal -- a record that
        cannot be written is a ledger line, not an exception.
        """

        # Once per (goal, cycle, stream). record_run appends and does
        # not deduplicate -- two calls on one stream wrote 57 rows
        # twice -- and this projection now has more than one caller,
        # because every route that ends a cycle must project the stream
        # it ended on.
        marker = (self.cycles, str(events_path))
        if marker in self._recorded_streams:
            return
        self._recorded_streams.add(marker)
        try:
            appended = record_run(
                self.workspace,
                goal_id=self.goal_id,
                cycle=self.cycles,
                run_events_path=events_path,
                run=run_reference,
                review_file=self.review_file,
            )
        except (
            Exception
        ) as exc:  # noqa: BLE001 -- the record never ends a goal
            self.ledger.append(
                "workspace_record_skipped",
                {
                    "cycle": self.cycles,
                    "reason": f"{type(exc).__name__}: {exc}",
                },
            )
            return
        if appended:
            self.ledger.append(
                "workspace_recorded",
                {"cycle": self.cycles, "rows": int(appended)},
            )

    def _record_analysis_evidence(self) -> None:
        """Name this cycle's own stream as the evidence a revision answers.

        A cycle that claimed from registered results and recorded a
        decision left durable typed evidence and no engine run, and the
        revision gate reads only runs -- so the next cycle's first
        executable plan was refused for answering an "unrecorded"
        outcome. The reference is the session's own run directory,
        relative to the workspace, and the wake embeds it exactly as it
        embeds a run, so no new act is asked of the model.
        """

        if self.events_path is None:
            return
        # Only a cycle that actually read something: a stream that ended
        # on a transport loss or in silence produced no evidence a
        # revision could answer, and naming it would let an empty cycle
        # satisfy the gate that exists to stop exactly that.
        #
        # What counts as having read is the question, and requiring a
        # *claim* answered it wrongly. po3-r17 cycle 1 (2026-09-11)
        # extracted eight quantities through the typed layer, measured
        # that two task-supplied structures carried swapped labels, and
        # recorded a decision -- then rendered no claim, because its
        # review was refused before any physics existed and there was
        # nothing yet to claim. Declining to claim was the honest act.
        # It left the cycle named by nothing, so the re-wake the charter
        # grants could never be admitted and the window ended with no
        # calculation. A typed read is the host reading a real quantity
        # out of real bytes and minting a receipt; that is evidence a
        # later plan can answer, whether or not a claim stood on it.
        delivery = _analysis_delivery(self.events_path)
        read_something = bool(
            delivery.claims or _session_typed_reads(self.events_path)
        )
        if not (read_something and delivery.decisions):
            return
        # Relative to .chemsmart-agent, the base run_recorded uses and
        # every consumer assumes: the wake resolves
        # workspace/.chemsmart-agent/<reference>/events.jsonl. Written
        # relative to the workspace it doubled the directory and
        # resolved to nothing.
        try:
            evidence = str(
                self.events_path.parent.relative_to(
                    self.workspace / ".chemsmart-agent"
                )
            )
        except (ValueError, AttributeError, TypeError):
            return
        self._defer_or_append(
            "analysis_evidence_recorded",
            {"cycle": self.cycles, "evidence": evidence},
        )

    def _record_analysis_only_revision(self) -> None:
        """An analysis-only plan the host walked at wake is a revision.

        It changes no identity, state or condition and launches no
        engine, so the admission checks hold by construction; the
        owner's ruling (2026-09-05) admits and executes it with no
        displayed review. The ledger records it as the revision it is,
        so the budget charges it and the wake's trajectory shows it.
        """

        if self.goal is None or self.events_path is None:
            return
        try:
            lines = self.events_path.read_text(encoding="utf-8").splitlines()
        except OSError:
            return
        for line in lines:
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get("kind") != "analysis_only_plan_executed":
                continue
            payload = event.get("payload") or {}
            self.revisions_admitted += 1
            self.ledger.append(
                "revision_admitted",
                {
                    "cycle": self.cycles,
                    "review_sha256": "",
                    "actor": self.goal.actor,
                    "granted_by": self.goal.granted_by,
                    "analysis_only": True,
                    "toolchain_plan_sha256": str(
                        payload.get("toolchain_plan_sha256") or ""
                    ),
                    "analysis_status": str(
                        payload.get("analysis_status") or ""
                    ),
                    "checks": {
                        "identity_preserved": True,
                        "conditions_preserved": True,
                        "programs_within_envelope": True,
                        "engine_calls_within_budget": True,
                    },
                    "cited_evidence_event_hashes": (),
                },
            )

    def _settle_delivery(self, terminal: str) -> GoalLoopResultV1:
        if self.events_path is None:
            # Nothing durable to read a delivery from: the goal still ends
            # in a typed state, and the reason says why a human reads it.
            reason = (
                f"session terminal state: {terminal}; the planning session "
                "left no event stream"
            )
            self.ledger.settle("returned_to_human", reasons=(reason,))
            return self._settled("returned_to_human", (reason,))
        self.result = _settle_from_delivery(
            self.ledger,
            goal_id=self.goal_id,
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            events_path=self.events_path,
            terminal=terminal,
            workspace=self.workspace,
        )
        self.phase = "settled"
        return self.result

    def _project_before_settling(self) -> None:
        """Make this cycle's delivery readable before the goal ends.

        A session that raised a typed error still produced everything it
        produced: SUFFICIENCY-2 recorded 57 claims, 11 declarations, an
        assessment and a scientific decision, and settled with a ledger
        holding one line and a workspace record holding none, because
        the projection runs after the planning session returns and the
        error returned first. The evidence was never destroyed -- it was
        made unreachable, which is the same thing to every later reader.
        Surviving an error and preserving what the error interrupted are
        two different properties, and only the first had been repaired.
        """

        events_path = self.events_path
        if events_path is None:
            # The session raised before the driver resolved its stream,
            # which is exactly the case that lost SUFFICIENCY-2's
            # delivery. The resolver's own fallback -- the newest
            # live-* stream in this workspace -- is what production
            # already trusts when a session returns no run id, and the
            # session that just raised is its newest writer.
            try:
                events_path = _session_events_path(None, self.workspace)
            except ContractError:
                return
            self.events_path = events_path
        if self.goal is None:
            # The goal record is created after the planning session
            # returns, so an error raised inside it left the ledger with
            # a settlement and no goal at all -- a malformed story, and
            # nothing for the record's rows to name. Create it from what
            # is known, exactly as the analysis-only path does: no
            # identity, no conditions, no approved review, which is
            # already the shape of a goal that has displayed no
            # executable partition. A settled goal never resumes, so
            # this grants nothing.
            self.goal = self._goal_record(
                identity="", conditions={}, review_sha256=""
            )
            self.ledger.create(self.goal)
        # What the session declared is part of what it delivered: the
        # goal's contract is unreadable without it, and the settlement
        # joins claims to declarations by id.
        self._record_declarations()
        self._flush_declarations()
        self._record_workspace(events_path, "")
        self._record_analysis_evidence()
        # And the run's own stream, when one exists: an executor error
        # raised after nodes completed returns before _outcome, which
        # owns the only other projection, so a cycle could lose real
        # engine wall time's evidence the same way.
        if self.run_directory is not None:
            run_events = self.run_directory / "events.jsonl"
            if run_events.is_file():
                self._record_workspace(
                    run_events,
                    f"goals/{self.goal_id}/runs/cycle-{self.cycles}",
                )

    def _typed_error(self, stage: str, error: ContractError) -> None:
        self._project_before_settling()
        self.result = _typed_error_settlement(
            self.ledger,
            goal_id=self.goal_id,
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            stage=stage,
            error=error,
        )
        self.phase = "settled"

    def _goal_record(
        self,
        *,
        identity: str,
        conditions: Mapping[str, Any],
        review_sha256: str,
    ) -> GoalRecordV1:
        return GoalRecordV1(
            schema_version=GOAL_SCHEMA_VERSION,
            goal_id=self.goal_id,
            task_spec_sha256=str(
                getattr(self.session, "task_spec_sha256", "") or ""
            ),
            scientific_identity_sha256=identity,
            conditions=conditions,
            envelope=self.envelope_record,
            max_revisions=self.max_revisions,
            granted_by=self.granted_by,
            initial_review_sha256=review_sha256,
            created_at=_utc_now(),
        )

    def _plan(self) -> None:
        self.cycles += 1
        if self._stopped():
            self.ledger.append(
                "cancelled_by_human",
                {"cycle": self.cycles, "at": _utc_now()},
            )
            self.ledger.settle(
                "returned_to_human",
                reasons=("cancelled by the human's stop file",),
            )
            self._settled("returned_to_human", ("cancelled",))
            return
        self.review_file = (
            self.goal_dir / "reviews" / f"cycle-{self.cycles}.json"
        )
        self.review_file.parent.mkdir(parents=True, exist_ok=True)
        self.wake = (
            _wake_context(
                self.goal,
                self.ledger,
                self.outcome,
                workspace=self.workspace,
                failure_report=self.failure_report,
            )
            if self.goal is not None
            else _goal_terms_context(
                goal_id=self.goal_id,
                granted_by=self.granted_by,
                envelope_record=self.envelope_record,
                max_revisions=self.max_revisions,
                workspace=self.workspace,
            )
        )
        if self.wake is not None and self.wake.get("previous_run"):
            # Host attestation for the evidence gate: this cycle's
            # session was handed the named run's typed outcome.
            self.ledger.append(
                "wake_composed",
                {"cycle": self.cycles, "run": self.wake["previous_run"]},
            )
        try:
            self.session = self.plan_session(
                task=self.task,
                provider=self.provider,
                provider_config_file=self.provider_config_file,
                workspace=self.workspace,
                execution_enabled=False,
                approval_file=None,
                execution_envelope_file=self.execution_envelope_file,
                analysis_completion_file=self.analysis_completion_file,
                review_file=self.review_file,
                goal_context=self.wake,
                **self.session_kwargs,
            )
        except ContractError as exc:
            self._typed_error("planning session", exc)
            return
        terminal = str(getattr(self.session, "terminal_state", "") or "")
        try:
            self.events_path = _session_events_path(
                self.session, self.workspace
            )
        except ContractError:
            # A session that left no stream can still be decided on if it
            # prepared a review; it cannot settle a delivery or carry the
            # evidence a revision must cite, and those paths say so.
            self.events_path = None
        self._record_declarations()
        self._record_dispositions()
        self._record_approaches()
        self._record_input_checks()
        if terminal != "waiting_for_approval":
            # No executable partition was planned. Either the session
            # delivered over registered results, refused with receipts,
            # or stopped; each settles the goal from durable evidence.
            if self.goal is None:
                self.goal = self._goal_record(
                    identity="",
                    conditions={"solvents": (), "thermochemistry": ()},
                    review_sha256="",
                )
                self.ledger.create(self.goal)
                self._flush_declarations()
            self._record_analysis_only_revision()
            self._record_analysis_evidence()
            try:
                reopened = self._rewake(terminal)
            except ContractError as exc:
                # The re-wake refuses a requirement state it has no
                # route for, deliberately -- but it is called outside
                # every handler, so the refusal would have ended the
                # goal unsettled. No node may end the goal.
                self._typed_error("goal re-wake", exc)
                return
            # Project first, then decide whether to re-open. A re-woken
            # cycle used to return here without recording anything, so
            # the one cycle the whole sufficiency mechanism exists to
            # produce was the one cycle whose delivery never reached the
            # record: SUFFICIENCY-3 lost 26 rows including its own
            # `attested` assessment, and had its next cycle delivered
            # anything else the goal-grain join would have fallen back
            # to a staler, worse row from two cycles earlier. The
            # projection is idempotent at its own layer, so calling it
            # on both paths is safe.
            if self.events_path is not None:
                self._record_workspace(self.events_path, "")
            if reopened:
                return
            self._settle_delivery(terminal)
            return
        self.phase = "decide"

    def _rewake(self, terminal: str) -> bool:
        """One further wake for a woken cycle that left the declared
        observables undelivered.

        NOVEL-2's po2 (2026-09-04) ended its woken cycle with a correct
        plan nothing could run and no claim, and the goal settled with
        sixteen engine calls, eight revisions and three and a half hours
        unspent. The owner ruled (2026-09-05) that such a cycle gets one
        more wake, carrying a failure report, charged as one revision.

        The first version keyed "delivered nothing" on the absence of
        claims. NOVEL-3 supplied the A/B in one window: po3 claimed four
        wall observations, none the declared observable, and settled with
        28 calls and 7 revisions unspent; ino3 claimed nothing, was
        re-woken, and delivered. A proxy for the invariant fails exactly
        where it diverges from it, and it rewarded claiming anything. The
        condition is now the invariant the settlement itself checks: a
        declared observable still undelivered by its id in any cycle,
        with budget in hand. A host-verified refusal, a session the host
        stopped, and a second silence settle as before.
        """

        # A cycle-1 session has no previous outcome and a session the
        # review builder refused has a stopped_by: both left the declared
        # observables undelivered with the grant in hand, and both are
        # what the one further wake is for (REACH-1 ino3 returned with 40
        # of 40 calls unspent on exactly this shape). Only a typed refusal
        # of the observables themselves, or a second silence, settles.
        if self.goal is None:
            return False
        if self.events_path is None:
            return False
        # `self.failure_report` used to be read here too, as a second
        # way of asking "has this goal already had its wake". It holds
        # the report composed for *this* cycle's wake, is never
        # cleared, and answered the question by a different means than
        # the ledger scan below -- which A12 taught to exempt a
        # transport continuation while this reader, two lines above it,
        # blocked unconditionally. Observed live: SUFFICIENCY-4 arm A
        # lost cycle 1 to four inter-event timeouts, delivered an
        # `attested` requirement at cycle 2, and settled with forty
        # engine calls and every revision unspent, so the arm that was
        # to receive the sufficiency consequence never did and the
        # window was void. Two organs answering one question call one
        # function; the ledger is the one that can see a transport
        # continuation.
        # One further wake per goal -- but a cycle the provider's
        # transport ended produced no scientific evidence and must not
        # spend it. SUFFICIENCY-1 lost cycle 1 to three inter-event
        # timeouts, the wake that answered them consumed the goal's only
        # opportunity, and the requirement wake could never fire
        # (2026-09-09). Transport continuation and scientific revision
        # are different budgets; the revision and wall budgets bound
        # both.
        if any(
            entry["kind"] == "rewake_opened"
            and not entry["payload"].get("transport_continuation")
            for entry in self.ledger.entries()
        ):
            return False
        # This reader used to build its own id-only union while the
        # settlement built a dimension-checked one from the same stream,
        # which is the defect b46290f1 was written for, one organ over:
        # a claim in the wrong dimension looked delivered here and
        # undelivered there, so the cycle that could have repaired it
        # was never opened and the goal settled naming it. Two organs
        # that answer one question call one function.
        delivery = _analysis_delivery(
            self.events_path,
            goal_delivered_ids=_goal_delivered_ids(
                self.workspace, self.goal_id
            ),
            declared_observables=_first_declarations(self.ledger),
        )
        if delivery.blocked_output_ids:
            return False
        declared = _required_declared_ids(self.ledger)
        delivered = set(
            observable_id
            for observable_id in set(
                _goal_delivered_ids(self.workspace, self.goal_id)
            )
            | set(delivery.delivered_quantity_ids)
            if delivery.answers_declaration(observable_id)
        )
        delivered.update(delivery.verified_unreachable_ids)
        undelivered = tuple(
            observable_id
            for observable_id in declared
            if observable_id and observable_id not in delivered
        )
        # A number can be delivered under its id and still not answer the
        # precision the task asked for. That gap had no home: OPEN-2's
        # ino3-qwen delivered eleven of eleven, wrote that its method
        # carries 0.2-0.4 V against the requester's +/-0.2 V, priced the
        # experiment that would narrow it at two of forty engine calls,
        # and the goal settled achieved with all forty unspent. The
        # obligation is to *resolve* the requirement -- compute it, show
        # re-claim it with a measured uncertainty, or refuse it -- and
        # never to meet it, because a tolerance no method in the
        # envelope can reach is a real answer and a better one than a
        # number.
        unresolved = delivery.unresolved_requirement_ids
        if not undelivered and not unresolved:
            return False
        budgets = self.ledger.budgets(self.goal)
        if budgets.revisions_remaining <= 0 or (
            budgets.wall_seconds_remaining <= 0
        ):
            return False
        rendered = tuple(
            item
            for item in delivery.delivered_quantity_ids
            if item not in declared
        )
        if not undelivered:
            self._open_requirement_rewake(terminal, delivery, unresolved)
            return True
        diagnosis = (
            f"the previous cycle ended {terminal!r} ({delivery.ending}) "
            "with these declared observables still undelivered by their "
            "id in any cycle: "
            + ", ".join(undelivered)
            + (
                "; it rendered claims under other names instead: "
                + ", ".join(rendered)
                if rendered
                else "; it rendered no claim"
            )
            + (
                "; no completion was certified"
                if delivery.completion_status != "passed"
                else "; the completion certified the chain without them"
            )
            + "."
        )
        self.failure_report = {
            "gate": "goal.cycle_delivers_or_returns",
            "invariant": (
                "a woken cycle ends by delivering every declared observable "
                "under its id, by a typed refusal the host can verify, or by "
                "an executable plan for review."
            ),
            "diagnosis": diagnosis,
            "route": (
                "claim each undelivered id from receipts in hand -- "
                "extract_result_quantities, derive_thermochemistry, "
                "evaluate_quantity_expression, record_analysis_claims with "
                "claim_id set to the declared id -- or "
                "plan_scientific_workflow with no calculation_nodes, which "
                "the host executes when planned; or record_scientific_"
                "decision naming each id that cannot be delivered, with its "
                "required producer and the receipts that show it"
            ),
            "cost": (
                "no engine call; this re-wake is charged one revision and "
                "is the last for this goal"
            ),
        }
        self.ledger.append(
            "rewake_opened",
            {
                "cycle": self.cycles,
                "terminal_state": terminal,
                "transport_continuation": is_provider_transport_terminal(
                    delivery.terminal_reason
                ),
                "undelivered_declared_observable_ids": list(undelivered),
                "failure_report": dict(self.failure_report),
            },
        )
        self.phase = "plan"
        return True

    def _open_requirement_rewake(
        self,
        terminal: str,
        delivery: "_AnalysisDelivery",
        unresolved: tuple[str, ...],
    ) -> None:
        """Wake a cycle whose numbers arrived short of what was asked.

        The host names the routes and never chooses one; the physics
        decides after the session acts. Two of the three cost no engine
        call, and refusing is a deliverable -- "no method in this
        envelope reaches that tolerance" is an answer, and for a
        requester it is often the more useful one.
        """

        rows = {
            str(row.get("observable_id") or ""): row
            for row in delivery.sufficiency
        }
        lines = []
        states = []
        for observable_id in unresolved:
            row = rows.get(observable_id, {})
            unit = str(row.get("unit") or "")
            state = str(row.get("state") or "")
            states.append(state)
            if state == "unstated":
                lines.append(
                    f"{observable_id}: delivered, and no uncertainty stands "
                    f"beside it against a required tolerance of "
                    f"{row.get('required_tolerance')} {unit}".rstrip()
                )
            elif state == "attested":
                unquantified = tuple(row.get("unquantified_components") or ())
                lines.append(
                    f"{observable_id}: uncertainty "
                    f"{row.get('uncertainty')} {unit} is within the "
                    f"required {row.get('required_tolerance')} {unit} and "
                    + (
                        "leaves these terms unquantified: "
                        + "; ".join(unquantified)
                        if unquantified
                        else "rests on your own word"
                    )
                )
            else:
                lines.append(
                    f"{observable_id}: uncertainty "
                    f"{row.get('uncertainty')} {unit} against a required "
                    f"{row.get('required_tolerance')} {unit}".rstrip()
                )
        self.failure_report = {
            "gate": "goal.requirement_is_resolved",
            "invariant": (
                "a requested observable is resolved when its stated "
                "uncertainty is within the tolerance the task asked for "
                "and rests on evidence the host can resolve, or when the "
                "session refuses it with receipts. "
                "A number delivered with no uncertainty beside it, or "
                "short of what was asked, with budget in hand, is none of "
                "those. An uncertainty is a judgement about a result, so "
                "a delivery that came from an executed plan carries none "
                "until a session looks at what came back."
            ),
            "diagnosis": (
                f"the previous cycle ended {terminal!r} ({delivery.ending}) "
                "having delivered every declared observable, with these "
                "requirements unresolved: " + "; ".join(lines)
            ),
            "route": sufficiency_menu(states),
            "cost": (
                "two of the three routes cost no engine call; this "
                "re-wake is charged one revision and is the last for this "
                "goal"
            ),
        }
        self.ledger.append(
            "rewake_opened",
            {
                "cycle": self.cycles,
                "terminal_state": terminal,
                "transport_continuation": is_provider_transport_terminal(
                    delivery.terminal_reason
                ),
                "unresolved_requirement_ids": list(unresolved),
                "sufficiency": [dict(row) for row in delivery.sufficiency],
                "failure_report": dict(self.failure_report),
            },
        )
        self.phase = "plan"

    def _decide(self) -> None:
        assert self.review_file is not None
        review = _review_record(self.review_file)
        if self.goal is None and not self.envelope_record:
            # No envelope file was given (the terminal interface
            # over a stored review): the envelope the human is
            # deciding on is the one the review itself displays.
            shown = dict(review.get("execution_envelope") or {})
            self.envelope_record = _goal_envelope_record(shown)
        # A goal record's existence is not evidence that an execution
        # grant exists. An analysis-only first cycle creates the goal
        # in _plan with an empty initial review, so a later cycle's
        # first executable review found `self.goal` already set, took
        # the revision path, and resolved its own review with
        # decision="approve" -- and `--initial-decision deny` was never
        # consulted. Observed live: deny, one engine partition
        # launched, settled achieved. The grant is what the human gave,
        # so the gate reads the grant.
        granted = self.goal is not None and not goal_scope_is_unbound(
            self.goal
        )
        if not granted:
            if self.initial_decision != "approve":
                if self.goal is None:
                    self.ledger.create(
                        self._goal_record(
                            identity=_plan_identity_sha256(review),
                            conditions=conditions_from_review(review),
                            review_sha256="",
                        )
                    )
                self.ledger.settle(
                    "returned_to_human",
                    reasons=("the human declined the initial review",),
                )
                self._settled(
                    "returned_to_human", ("initial review declined",)
                )
                return
        if self.goal is None:
            review_sha256, self.bundle_file = self.resolve_review(
                review_file=self.review_file,
                workspace=self.workspace,
                decision="approve",
                actor=self.granted_by,
                approval_id=f"goal-{self.goal_id}-cycle-{self.cycles}",
            )
            self.goal = self._goal_record(
                identity=_plan_identity_sha256(review),
                conditions=conditions_from_review(review),
                review_sha256=review_sha256,
            )
            self.ledger.create(self.goal)
            self._flush_declarations()
            self.phase = "execute"
            return
        if self.events_path is None:
            raise ContractError(
                "a revision must come from a session with an event stream"
            )
        budgets = self.ledger.budgets(self.goal)
        bound_identity, bound_conditions = _goal_bound_scope(self.ledger)
        verdict = admit_revision(
            goal=self.goal,
            budgets=budgets,
            revision_review=review,
            revision_scientific_identity_sha256=_plan_identity_sha256(review),
            session_events_path=self.events_path,
            prior_outcome_evidence_hashes=tuple(
                digest
                for node in (self.outcome.nodes if self.outcome else ())
                for digest in node.evidence_event_hashes
            ),
            previous_run_reference=_previous_run_reference(self.ledger),
            wake_embedded_run=str((self.wake or {}).get("previous_run") or ""),
            bound_scientific_identity_sha256=bound_identity,
            bound_conditions=bound_conditions,
        )
        if verdict.admitted and verdict.bound_scientific_identity_sha256:
            # The goal record is digest-bound and never rewritten, so the
            # scope this revision established lives on the ledger, where
            # every later cycle reads it and is held to it.
            self.ledger.append(
                "goal_scope_bound",
                {
                    "cycle": self.cycles,
                    "scientific_identity_sha256": (
                        verdict.bound_scientific_identity_sha256
                    ),
                    "conditions": canonical_data(
                        verdict.bound_conditions or {}
                    ),
                },
            )
        if not verdict.admitted:
            self.ledger.append(
                "revision_returned",
                {"cycle": self.cycles, "reasons": verdict.reasons},
            )
            self.ledger.settle("returned_to_human", reasons=verdict.reasons)
            self._settled("returned_to_human", verdict.reasons)
            return
        review_sha256, self.bundle_file = self.resolve_review(
            review_file=self.review_file,
            workspace=self.workspace,
            decision="approve",
            actor=self.goal.actor,
            approval_id=f"goal-{self.goal_id}-cycle-{self.cycles}",
        )
        self.revisions_admitted += 1
        self.ledger.append(
            "revision_admitted",
            {
                "cycle": self.cycles,
                "review_sha256": review_sha256,
                "actor": self.goal.actor,
                "granted_by": self.goal.granted_by,
                "checks": dict(verdict.checks),
                "cited_evidence_event_hashes": (
                    verdict.cited_evidence_event_hashes
                ),
            },
        )
        self.phase = "execute"

    def _execute(self) -> None:
        assert self.bundle_file is not None
        self.run_directory = self.goal_dir / "runs" / f"cycle-{self.cycles}"
        self.run_directory.mkdir(parents=True, exist_ok=True)
        run_reference = f"goals/{self.goal_id}/runs/cycle-{self.cycles}"
        if self.dispatch_run is not None:
            try:
                receipt = self.dispatch_run(
                    approval_file=self.bundle_file,
                    workspace=self.workspace,
                    run_directory=self.run_directory,
                    goal_id=self.goal_id,
                    cycle=self.cycles,
                )
            except ContractError as exc:
                self._typed_error("scheduler dispatch", exc)
                return
            self.dispatch_receipt = receipt
            payload = {
                "cycle": self.cycles,
                "run": run_reference,
                **{
                    key: value
                    for key, value in _record_of(receipt).items()
                    if key
                    in {
                        "scheduler",
                        "job_id",
                        "submitted_at",
                        "submit_script",
                        "wake_job_id",
                    }
                },
            }
            self.ledger.append("run_dispatched", payload)
            self.result = GoalLoopResultV1(
                goal_id=self.goal_id,
                settlement="parked",
                cycles=self.cycles,
                revisions_admitted=self.revisions_admitted,
                reasons=(
                    f"cycle {self.cycles} submitted as "
                    f"{payload.get('scheduler', 'scheduler')} job "
                    f"{payload.get('job_id', '?')}; resume with "
                    f"chemsmart agent wake --goal {self.goal_id}",
                ),
            )
            self.phase = "parked"
            return
        # The execute boundary is a ledger entry like every other, so a
        # process killed mid-engine leaves a run that agent wake can
        # re-enter: the run directory and the one-shot bundle are named
        # here, and the executor's own continuation replays what
        # finished (live, 2026-09-03: a launcher timeout killed a goal
        # three relaxations into its second cycle and nothing could
        # resume it).
        prior = _goal_anomalies(self.ledger)
        if prior:
            (self.run_directory / PRIOR_ANOMALIES_FILE).write_text(
                json.dumps(canonical_data(prior), sort_keys=True) + "\n",
                encoding="utf-8",
            )
        self.ledger.append(
            "run_started",
            {
                "cycle": self.cycles,
                "run": run_reference,
                "approval_file": str(self.bundle_file),
            },
        )
        try:
            self.execute_result = self.execute_bundle(
                approval_file=self.bundle_file,
                workspace=self.workspace,
                run_directory=self.run_directory,
                **(
                    {"stop_file": Path(self.stop_file)}
                    if self.stop_file is not None
                    and _execute_hook_takes_stop_file(self.execute_bundle)
                    else {}
                ),
            )
        except ContractError as exc:
            self._typed_error("approved execution", exc)
            return
        self.phase = "outcome"

    def _outcome(self) -> None:
        assert self.run_directory is not None
        from chemsmart.agent.terminal_states import (
            derive_run_outcome,
            read_run_events,
        )

        run_reference = f"goals/{self.goal_id}/runs/cycle-{self.cycles}"
        events_path = self.run_directory / "events.jsonl"
        try:
            self.outcome = derive_run_outcome(read_run_events(events_path))
        except (ValueError, OSError) as exc:
            if isinstance(exc, ValueError) and "found 0" not in str(exc):
                raise
            # An admitted revision whose approved bundle launched no
            # engine: the executor walked the analysis chain into the
            # run stream and recorded no workflow run. Observed live
            # (C8): the cycle is a legitimate delivery-or-return, and
            # letting the derivation's own contract error escape left
            # the goal unsettled. Settle from the run stream's typed
            # delivery, exactly as a no-partition planning cycle does.
            self.ledger.append(
                "run_recorded",
                {
                    "cycle": self.cycles,
                    "run": run_reference,
                    "workflow_state": "analysis_only",
                    "engine_calls_consumed": 0,
                    "engine_wall_seconds": 0.0,
                    "stopped_by": list(
                        _analysis_delivery(events_path).stopped_by
                    ),
                },
            )
            self.events_path = events_path
            self._record_workspace(events_path, run_reference)
            self._settle_delivery("complete")
            return
        payload = {
            "cycle": self.cycles,
            "run": run_reference,
            "workflow_state": self.outcome.workflow_state,
            "engine_calls_consumed": self.outcome.engine_calls_consumed,
            "engine_wall_seconds": self.outcome.engine_wall_seconds,
            "host_seconds": float(
                getattr(self.outcome, "host_seconds", 0.0) or 0.0
            ),
            "excursion_calls_consumed": (
                self.outcome.excursion_calls_consumed
            ),
        }
        timed = [
            (float(node.wall_seconds), node.node_id)
            for node in self.outcome.nodes
            if getattr(node, "wall_seconds", None) is not None
        ]
        if timed:
            slowest_seconds, slowest_node = max(timed)
            payload["slowest_node_seconds"] = slowest_seconds
            payload["slowest_node_id"] = slowest_node
        if self.dispatch_receipt is not None:
            payload["queue_wait_seconds"] = _queue_wait_seconds(
                self.dispatch_receipt, events_path
            )
        self.ledger.append("run_recorded", payload)
        self._record_workspace(events_path, run_reference)
        # What the host detected belongs to the goal, not to the host
        # that detected it: two live goals settled plain "achieved" over
        # recorded anomalies because the receipt died with its cycle
        # (E2 window, 2026-09-03). The ledger carries every anomaly and
        # every later host is seeded from it.
        observed = tuple(
            {**dict(item), "node_id": node.node_id}
            for node in self.outcome.nodes
            for item in node.anomalies
        )
        if observed:
            self.ledger.append(
                "anomalies_observed",
                {
                    "cycle": self.cycles,
                    "run": run_reference,
                    "anomalies": observed,
                },
            )
        if self.execute_result is None:
            self.execute_result = _recorded_execution_result(
                self.run_directory
            ) or SimpleNamespace(
                status=(
                    "completed"
                    if self.outcome.workflow_state == "completed"
                    else self.outcome.workflow_state
                ),
                analysis_status="",
            )
        self.phase = "settle"

    def _record_qualification(self) -> None:
        _record_goal_qualification(
            self.ledger,
            workspace=self.workspace,
            goal_id=self.goal_id,
            current_run=f"goals/{self.goal_id}/runs/cycle-{self.cycles}",
            current_outcome=self.outcome,
        )

    def _settle(self) -> None:
        assert self.run_directory is not None and self.goal is not None
        # Built first, because the run delivery inherits its verified
        # refusals: the refusal lives in this stream and the assessment
        # it answers lives in the run's.
        session_delivery = (
            _analysis_delivery(
                self.events_path,
                # The one reader of six still built without them, so a
                # session that retired a wrong-unit observable and
                # declared its replacement kept the retired id open.
                declared_observables=_first_declarations(self.ledger),
            )
            if self.events_path is not None
            else None
        )
        run_delivery = _analysis_delivery(
            self.run_directory / "events.jsonl",
            # The verified refusals only: the bases mapping explains
            # every refusal the session wrote, verified or not, and
            # passing it handed authority to refusals the host had
            # declined to verify.
            inherited_unreachable=(
                {
                    observable_id: session_delivery.unreachable_bases.get(
                        observable_id, ""
                    )
                    for observable_id in (
                        session_delivery.verified_unreachable_ids
                    )
                }
                if session_delivery is not None
                else {}
            ),
            # The settlement was the one reader of five that omitted
            # these, so a retired observable's assessment kept an
            # executed goal open for ever after a legitimate unit
            # correction, and the reason printed an id the record says
            # was retired. I dropped it here to keep an earlier commit
            # scoped and never put it back.
            declared_observables=_first_declarations(self.ledger),
            inherited_rejected_artifacts=tuple(
                sorted(self.rejected_artifacts)
            ),
            uncharacterised_artifact_sha256s=uncharacterised_artifacts(
                self.workspace
            ),
            flagged_artifact_sha256s=_flagged_artifact_sha256s(
                _goal_anomalies(self.ledger)
            ),
            goal_delivered_ids=_goal_delivered_ids(
                self.workspace, self.goal_id
            ),
        )
        self.rejected_artifacts.update(run_delivery.rejected_artifact_sha256s)
        if run_delivery.claims_rendered:
            self.standing_stale = run_delivery.stale_quantity_ids
        unrefreshed = self.standing_stale
        budgets = self.ledger.budgets(self.goal)
        # A verdict or a stale number needs an engine to answer it; a
        # claim that was never rendered, or a declared observable no
        # claim carries by id, is answered from receipts already in
        # hand and costs no engine call. Two goals settled exhausted
        # with every receipt they needed on disk (E4 window).
        engine_needed = bool(run_delivery.unanswered_verdicts or unrefreshed)
        recovery_affordable = (
            budgets.revisions_remaining > 0
            and budgets.wall_seconds_remaining > 0
            and (budgets.engine_calls_remaining > 0 or not engine_needed)
        )
        # A typed refusal the host verified, recorded by the planning
        # session of this cycle, closes the ids it names.
        refused = set(
            session_delivery.verified_unreachable_ids
            if session_delivery is not None
            else ()
        )
        open_declared = tuple(
            observable_id
            for observable_id in run_delivery.undelivered_declared_ids
            if observable_id not in refused
        )
        # A requirement is open when the number that answers it does not
        # answer the precision it was asked for. The claim may have been
        # rendered by the executor's own analysis walk or by this
        # cycle's session, so both streams are read; a host-verified
        # refusal closes the id here as it does for an undelivered one.
        # This settlement consulted neither, so a goal that actually ran
        # engines reached achieved with its requirement unresolved --
        # the one path the whole mechanism was built for
        # (SUFFICIENCY-1, 2026-09-09).
        # The predicate subtracts refusals and retirements itself, for
        # every reader; this settlement used to do it here and the
        # analysis-only one did not do it at all.
        open_requirements = tuple(
            dict.fromkeys(
                run_delivery.unresolved_requirement_ids
                + (
                    session_delivery.unresolved_requirement_ids
                    if session_delivery is not None
                    else ()
                )
            )
        )
        open_delivery = bool(
            run_delivery.unanswered_verdicts
            or unrefreshed
            or run_delivery.unclaimed_output_ids
            or open_declared
            or open_requirements
        )
        achieved = _achieved(self.execute_result)
        # An observable the session refused, and an observable it
        # delivered whose precision it refused, settle by one rule: the
        # host verified both as unreachable, and the charter calls that
        # ending a deliverable. Only the first was ever read, so a goal
        # that ran engines and had its stated accuracy certified
        # unreachable settled achieved.
        refused_ids = tuple(
            dict.fromkeys(
                tuple(run_delivery.undelivered_declared_ids)
                + tuple(run_delivery.refused_requirement_ids)
            )
        )
        if (
            achieved
            and not open_delivery
            and refused
            and refused_ids
            and set(refused_ids) <= refused
            and session_delivery is not None
            and session_delivery.decisions
        ):
            reason = (
                f"cycle {self.cycles}: the recorded decision names these "
                "declared observables unreachable from the admissible "
                "evidence, and the host verified each: "
                + "; ".join(
                    f"{observable_id} -- "
                    f"{session_delivery.unreachable_bases.get(observable_id, '')}"
                    for observable_id in refused_ids
                )
            )
            self.ledger.settle(
                "unreachable_from_evidence",
                reasons=(reason,),
                evidence=_settlement_evidence(session_delivery),
            )
            self._settled("unreachable_from_evidence", (reason,))
            return
        if achieved and open_delivery and recovery_affordable:
            # Every required output arrived and a host-rendered verdict
            # says one of the structures behind them is not what the
            # task required. Two live cases settled achieved here in one
            # cycle with four engine calls unspent, so no session was
            # ever asked whether it wanted to answer the failure it had
            # just been shown. Withhold the word and wake another cycle
            # while a recovery is still affordable; the session may
            # recover, or cite the validation receipt in its decision
            # and stand by the delivery. What it may not do is have the
            # choice made for it by a settlement.
            self.ledger.append(
                "recovery_opened",
                {
                    "cycle": self.cycles,
                    "verdicts": list(run_delivery.unanswered_verdicts),
                    "stale_quantity_ids": list(unrefreshed),
                    "unclaimed_output_ids": list(
                        run_delivery.unclaimed_output_ids
                    ),
                    "undelivered_declared_observable_ids": list(open_declared),
                    "unresolved_requirement_ids": list(open_requirements),
                    "engine_calls_remaining": budgets.engine_calls_remaining,
                },
            )
            self.phase = "plan"
            return
        if achieved and open_delivery:
            # Same failure, nothing left to answer it with.
            if run_delivery.unanswered_verdicts:
                reason = (
                    f"cycle {self.cycles}: a validation verdict failed and "
                    "no budget remains to answer it: "
                    + ", ".join(run_delivery.unanswered_verdicts)
                )
                open_items = run_delivery.unanswered_verdicts
            elif unrefreshed:
                reason = (
                    f"cycle {self.cycles}: a verdict rejected the result "
                    "these quantities were computed from and no budget "
                    "remains to re-derive them: " + ", ".join(unrefreshed)
                )
                open_items = unrefreshed
            elif run_delivery.unclaimed_output_ids:
                # The host computed these from real program output and
                # no reader of the delivery can see them.
                reason = (
                    f"cycle {self.cycles}: these quantities were computed "
                    "and never rendered as a claim, and no budget remains "
                    "to deliver them: "
                    + ", ".join(run_delivery.unclaimed_output_ids)
                )
                open_items = run_delivery.unclaimed_output_ids
            elif open_requirements:
                reason = (
                    f"cycle {self.cycles}: these delivered numbers do not "
                    "answer the precision their declaration asked for, and "
                    "no revision remains to resolve them: "
                    + ", ".join(open_requirements)
                )
                open_items = open_requirements
            else:
                reason = (
                    f"cycle {self.cycles}: these declared observables have "
                    "no claim carrying their id in any cycle, and no "
                    "revision remains to claim them: "
                    + ", ".join(run_delivery.undelivered_declared_ids)
                    + (
                        "; " + " | ".join(run_delivery.open_declared_misses)
                        if run_delivery.open_declared_misses
                        else ""
                    )
                )
                open_items = run_delivery.undelivered_declared_ids
            self.ledger.settle("returned_to_human", reasons=(reason,))
            self._settled("returned_to_human", open_items)
            return
        if achieved:
            goal_anomalies = _goal_anomalies(self.ledger)
            word, why = _achieved_word(run_delivery, goal_anomalies)
            evidence = _settlement_evidence(run_delivery)
            if word == "achieved_with_observations":
                evidence = _anomaly_evidence(evidence, goal_anomalies)
            self.ledger.settle(
                word,
                reasons=(
                    f"cycle {self.cycles}: workflow completed with its "
                    "analysis chain; " + why[0],
                    # Every reason the word carries, not only the first:
                    # the second one names which delivered number stands
                    # on a flagged result.
                    *why[1:],
                ),
                evidence=evidence,
            )
            self._record_qualification()
            self._settled(word, why)
            return
        if (
            budgets.engine_calls_remaining <= 0
            and budgets.revisions_remaining > 0
            and budgets.wall_seconds_remaining > 0
            and run_delivery.undelivered_declared_ids
        ):
            # The engine line is spent and the run did not complete,
            # yet declared observables sit unclaimed with receipts in
            # hand: an analysis-only cycle can still claim them by id,
            # and the plan-time budget gate refuses any engine node.
            self.ledger.append(
                "recovery_opened",
                {
                    "cycle": self.cycles,
                    "undelivered_declared_observable_ids": list(
                        run_delivery.undelivered_declared_ids
                    ),
                    "engine_calls_remaining": 0,
                    "analysis_only": True,
                },
            )
            self.phase = "plan"
            return
        if (
            budgets.revisions_remaining <= 0
            or budgets.engine_calls_remaining <= 0
            or budgets.wall_seconds_remaining <= 0
        ):
            self.ledger.settle(
                "exhausted",
                reasons=(
                    f"engine calls remaining "
                    f"{budgets.engine_calls_remaining}, wall seconds "
                    f"remaining {budgets.wall_seconds_remaining:.0f}, "
                    f"revisions remaining "
                    f"{budgets.revisions_remaining}",
                ),
            )
            self._settled("exhausted", ("budgets exhausted",))
            return
        # A run that did not complete, with budget in hand. Which way
        # its nodes ended decides whether a revision can answer it:
        # a wrong stationary point, a convergence failure, a timeout are
        # ordinary work; a launch that never happened, an admission
        # refusal, a cancellation or an ambiguous termination are not
        # evidence a revision can stand on, and the human reads them.
        # A stage the plan declared non-executable was displayed with the
        # review, never approved and never launched: it has no ending to
        # answer. Counting its not_launched as one returned a goal whose
        # every executable node had validated (live, 2026-09-03).
        retained = set(self._declared_non_executable_ids())
        terminal_states = {
            str(node.node_id): str(node.state)
            for node in (self.outcome.nodes if self.outcome else ())
            if str(node.state) != "validated"
            and not (
                str(node.state) == "not_launched"
                and str(node.node_id) in retained
            )
        }
        repairable = {
            node_id: state
            for node_id, state in terminal_states.items()
            if state in REPAIRABLE_TERMINAL_STATES
        }
        if terminal_states and not repairable:
            reason = (
                f"cycle {self.cycles}: the run ended in a state no revision "
                "can answer: "
                + ", ".join(
                    f"{node_id}={state}"
                    for node_id, state in sorted(terminal_states.items())
                )
            )
            self.ledger.settle("returned_to_human", reasons=(reason,))
            self._settled("returned_to_human", (reason,))
            return
        self.ledger.append(
            "recovery_opened",
            {
                "cycle": self.cycles,
                "terminal_states": dict(sorted(repairable.items())),
                # Observed live (W1): every node validated and the
                # analysis chain was partial, so the entry named no
                # node; say what opened the cycle.
                "analysis_status": str(
                    getattr(self.execute_result, "analysis_status", "") or ""
                ),
                "verdicts": [],
                "stale_quantity_ids": list(unrefreshed),
                "unclaimed_output_ids": list(
                    run_delivery.unclaimed_output_ids
                ),
                "engine_calls_remaining": budgets.engine_calls_remaining,
            },
        )
        self.phase = "plan"


def _resolved_or_none(path: str | Path | None) -> str | None:
    return str(Path(path).resolve()) if path is not None else None


def _record_goal_qualification(
    ledger: GoalLedger,
    *,
    workspace: Path,
    goal_id: str,
    current_run: str,
    current_outcome: Any,
) -> tuple[dict[str, Any], ...]:
    """An achieved goal is live evidence: write it to the goal ledger and
    the host's qualification store, so "qualified" is a fact the
    capability registry can read rather than a claim.

    Called from every path that settles ``achieved``: the executed run's
    and the analysis-only cycle's. The second had none -- g5 settled
    ``achieved_with_observations`` from a cycle that launched nothing,
    with two validated PySCF nodes in its cycle 1, and qualified nothing
    while the repaired reader sat one call away (PySCF round,
    2026-09-12): right where computed, unconnected where consumed.
    """

    from chemsmart.agent.capability_registry import (
        record_host_qualification,
    )
    from chemsmart.agent.terminal_states import (
        derive_run_outcome,
        read_run_events,
    )

    agent_root = Path(workspace) / ".chemsmart-agent"

    def _read_outcome(run: str) -> Any:
        return derive_run_outcome(
            read_run_events(agent_root / run / "events.jsonl")
        )

    entries = _qualification_entries_for_goal(
        ledger_entries=ledger.entries(),
        goal_id=goal_id,
        current_run=current_run,
        current_outcome=current_outcome,
        read_outcome=_read_outcome,
    )
    for entry in entries:
        ledger.append("qualified", entry)
    try:
        record_host_qualification(entries)
    except OSError:
        # The host store is a convenience mirror; the ledger is the
        # durable record.
        pass
    return entries


def _qualification_entries_for_goal(
    *,
    ledger_entries: Sequence[Mapping[str, Any]],
    goal_id: str,
    current_run: str,
    current_outcome: Any,
    read_outcome: Callable[[str], Any],
) -> tuple[dict[str, Any], ...]:
    """Every validated node of every run the goal recorded.

    The rows were read from the settling cycle's outcome alone, so a goal
    that executed in cycle one and settled ``achieved`` in an analysis-
    only cycle two wrote none for the nodes it ran (PySCF round, g2,
    2026-09-12). Each recorded run is read once; one it cannot read
    contributes nothing.
    """

    outcomes: dict[str, Any] = {}
    for entry in ledger_entries:
        if str(entry.get("kind") or "") != "run_recorded":
            continue
        run = str((entry.get("payload") or {}).get("run") or "")
        if not run or run in outcomes or run == current_run:
            continue
        try:
            outcomes[run] = read_outcome(run)
        except (KeyError, ValueError, OSError):
            continue
    outcomes[current_run] = current_outcome
    entries: list[dict[str, Any]] = []
    for run, outcome in outcomes.items():
        entries.extend(
            _qualification_entries(outcome, goal_id=goal_id, run=run)
        )
    return tuple(entries)


def _qualification_entries(
    outcome: Any, *, goal_id: str, run: str
) -> tuple[dict[str, Any], ...]:
    """What an achieved run qualifies: every validated node's program,
    engine and jobtype, with the run and receipts that say so."""

    entries: list[dict[str, Any]] = []
    for node in getattr(outcome, "nodes", ()):
        if str(getattr(node, "state", "")) != "validated":
            continue
        program = str(getattr(node, "program", "") or "")
        jobtype = str(getattr(node, "jobtype", "") or "")
        if not program or not jobtype:
            continue
        entries.append(
            {
                "kind": "program_jobtype",
                "id": f"{program}:cpu:{jobtype}",
                "goal": goal_id,
                "run": run,
                "node": str(getattr(node, "node_id", "")),
                "evidence_event_hashes": list(
                    getattr(node, "evidence_event_hashes", ())
                ),
                "date": _utc_now(),
            }
        )
    return tuple(entries)


def _record_of(receipt: Any) -> dict[str, Any]:
    """A receipt's public fields, whatever kind of object carries them."""

    if isinstance(receipt, Mapping):
        return dict(receipt)
    if hasattr(receipt, "public_record"):
        return dict(receipt.public_record())
    if hasattr(receipt, "__dataclass_fields__"):
        from dataclasses import asdict

        return asdict(receipt)
    return dict(vars(receipt))


#: Where a dispatched run's job writes the executor's own result record,
#: so the driver that resumes it reads the same status words a local run
#: returns in memory.
EXECUTION_RESULT_FILE = "execution-result.json"


def _recorded_execution_result(run_directory: Path) -> Any:
    path = Path(run_directory) / EXECUTION_RESULT_FILE
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(record, Mapping):
        return None
    return SimpleNamespace(
        status=str(record.get("status") or ""),
        analysis_status=str(record.get("analysis_status") or ""),
    )


def _queue_wait_seconds(receipt: Any, events_path: Path) -> float | None:
    """Seconds between submission and the run's first recorded event."""

    submitted = str(_record_of(receipt).get("submitted_at") or "")
    if not submitted:
        return None
    try:
        first = json.loads(
            events_path.read_text(encoding="utf-8").splitlines()[0]
        )
        # The runtime stream stamps each event as ``timestamp``; the
        # goal ledger stamps its own entries as ``at``. Observed live:
        # reading the ledger's word here left every queue wait None.
        started = str(first.get("timestamp") or first.get("at") or "")
        start = datetime.fromisoformat(started)
        submit = datetime.fromisoformat(submitted)
    except (OSError, IndexError, ValueError, json.JSONDecodeError):
        return None
    return max(0.0, (start - submit).total_seconds())


def run_goal_loop(**kwargs: Any) -> GoalLoopResultV1:
    """Drive one goal to settlement or parking in this process."""

    return GoalDriver(**kwargs).run()


__all__ = [
    "REPAIRABLE_TERMINAL_STATES",
    "REPAIR_MENU",
    "DISPATCH_MODES",
    "DRIVER_FILE",
    "EXECUTION_RESULT_FILE",
    "TASK_FILE",
    "GOAL_PHASES",
    "GoalDriver",
    "GoalLoopResultV1",
    "GoalStepV1",
    "run_goal_loop",
]
