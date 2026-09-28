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
import time
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable, Collection, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    require_identifier,
)
from chemsmart.agent.delivery import (
    current_assessments,
    superseded_observable_ids,
    unresolved_requirement_ids,
)
from chemsmart.agent.execution import STRUCTURE_MOVING_STAGES, anomaly_standing
from chemsmart.agent.goal import (
    GOAL_SCHEMA_VERSION,
    FailedCriterionV1,
    GoalLedger,
    GoalRecordV1,
    admit_revision,
    cited_receipts,
    conditions_from_review,
    failed_criteria,
    goal_scope_is_unbound,
    results_read,
)
from chemsmart.agent.rules import rules_by_id
from chemsmart.agent.terminal_states import (
    GEOMETRY_SEARCH_JOBTYPES,
    PRIOR_ANOMALIES_FILE,
    REPAIRABLE_NODE_STATES,
    is_provider_transport_terminal,
)
from chemsmart.agent.workspace_record import (
    failed_artifacts,
    printed_modes,
    read_workspace_record,
    record_run,
    recorded_surface,
    render_workspace_record,
    uncharacterised_artifacts,
)
from chemsmart.analysis.result_readers import surfaces_agree


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


#: How a settlement quotes verified refusals: the session's statement, and
#: in brackets what the host checked. "the host verified each: <statement>"
#: read as the host vouching for the session's reason -- a PubChem outage,
#: a writer that dropped the IRC block (r9 xtb g1, r9 gaussian g2) -- when
#: it had checked only that the plan retained a node as blocked.
_VERIFIED_REFUSAL_LEAD = (
    "the recorded decision names these declared observables unreachable "
    "from the admissible evidence; the host verified each on the basis in "
    "brackets, and the text before the brackets is the session's: "
)


def _anomaly_evidence(
    evidence: Mapping[str, Any] | None,
    ledger_anomalies: Sequence[Mapping[str, Any]],
    delivery: "_AnalysisDelivery | None" = None,
    superseded: Collection[str] = (),
) -> dict[str, Any]:
    """The receipts a settlement with observations stands on.

    The anomaly receipts are the observation's own evidence: decisions
    live in a session's stream and the executor's run stream carries
    none, so a word that requires receipts must bring its own. So does an
    expectation an earlier cycle's completion scored: the word names it
    from that completion, and cites it. So does a failed criterion a
    decision in another stream answered: the word cites the receipts that
    state the verdict, the one the decision cites among them.
    """

    digests = (
        {
            str(item.get("receipt_sha256") or "")
            for item in ledger_anomalies
            if str(item.get("receipt_sha256") or "")
        }
        | {
            str(row.get("completion_receipt_sha256") or "")
            for row in (delivery.carried_expectations if delivery else ())
            if str(row.get("completion_receipt_sha256") or "")
            and str(row.get("observable_id") or "") not in set(superseded)
        }
        | {
            str(receipt)
            for verdict, _standing in (
                delivery.answered_criteria if delivery else ()
            )
            for receipt in (*verdict.receipt_sha256s, *verdict.answered_by)
            if receipt
        }
    )
    merged = dict(evidence or {})
    receipts = set(merged.get("receipt_sha256s") or ()) | digests
    if receipts:
        merged["receipt_sha256s"] = tuple(sorted(receipts))
    return merged


def _finding_reasons(
    findings: Sequence[Mapping[str, Any]],
) -> tuple[str, ...]:
    """One settlement reason per standing finding, in the session's words.

    The session's conclusions in its own words, each on relations the
    host checked and on nothing else the host vouches for; a finding
    standing on a result a sensor had already flagged says which,
    because repeating a sensor is not a discovery. A declared category's
    answer is the word the host read, stated first; the sentence beside
    it is the session's interpretation.
    """

    return tuple(
        (
            f"{row.get('answers_observable_id')} = "
            + ", ".join(
                f"{word.get('word')!r} (read by the host: "
                f"{word.get('selector') or 'selector unrecorded'} on "
                f"{str(word.get('source_receipt_sha256') or '')[:8]})"
                for word in row.get("answer") or ()
            )
            + f"; the session's finding {row.get('finding_id')}, "
            "its interpretation: "
            if row.get("answers_observable_id") and row.get("answer")
            else f"the session's finding {row.get('finding_id')}"
            + (
                f" (names {row.get('answers_observable_id')} and rests "
                "on no word the host read, so it answers nothing)"
                if row.get("answers_observable_id")
                else (
                    " (not asked for)"
                    if row.get("standing") == "unrequested"
                    else " (on the requested answer)"
                )
            )
            + ", its words, on relations the host checked: "
        )
        + str(row.get("statement"))
        + (
            "; host anomalies already under its evidence: "
            + ", ".join(row.get("host_signals") or ())
            if row.get("host_signals")
            else ""
        )
        for row in findings
    )


def _claimed_answer_reasons(
    claimed_answers: Mapping[str, Sequence[Mapping[str, Any]]],
) -> tuple[str, ...]:
    """One reason per declared category a claim answered under its id.

    The word the host read, stated first and in the same form a finding's
    answer is, with what read it: a category is delivered as a number
    is, by the claim that carries its id.
    """

    return tuple(
        f"{observable_id} = "
        + ", ".join(
            f"{word.get('word')!r} (read by the host: "
            f"{word.get('selector') or 'selector unrecorded'} on "
            f"{str(word.get('source_receipt_sha256') or '')[:8]})"
            for word in words
        )
        + ", claimed under the declared id"
        for observable_id, words in sorted(claimed_answers.items())
    )


def _answered_criterion_reasons(
    answered: Sequence[tuple[FailedCriterionV1, tuple[str, ...]]],
) -> tuple[str, ...]:
    """One reason per failed criterion a recorded decision answered.

    The verdict with the number that failed it, the receipt the decision
    cites, and the delivered numbers standing on it: the delivery carries
    the failed expectation and the session's reading of it, and the word
    says which finding it rests on.
    """

    return tuple(
        "the plan's own acceptance criterion did not hold and the recorded "
        f"decision cites it ({verdict.answered_by[-1][:8]}): "
        f"{verdict.statement()}"
        + (
            "; delivered standing on it: " + ", ".join(standing)
            if standing
            else ""
        )
        for verdict, standing in answered
    )


def _unanswered_verdicts_named(delivery: "_AnalysisDelivery") -> str:
    """Each unanswered verdict as the host read it: rule, number, receipt."""

    return "; ".join(
        f"{verdict.statement()} (receipt {verdict.receipt_sha256s[-1][:8]})"
        for verdict in delivery.unanswered_criteria
    ) or ", ".join(delivery.unanswered_verdicts)


def _undelivered_declared_named(
    delivery: "_AnalysisDelivery", observable_ids: Sequence[str]
) -> str:
    """Which undelivered declared observables a claim carries, and which none.

    "No claim carrying their id in any cycle" was written over ids a claim
    carried: G-h2's three category questions were claimed under their own
    ids as the verdict numbers 1 and 0 in two cycles (R10 Q22, CUHK
    2153627), and ax41's ino3-cont and ino3-r13a carried sixteen declared
    ids as quantity ids in another unit -- 4 of the 11 such statements the
    archives let one check. A claim of this stream or a row of the goal's
    record under the id carries it; what it lacks is on the completion
    receipt's own miss text, which the callers carry beside this.
    """

    carried = tuple(
        observable_id
        for observable_id in observable_ids
        if observable_id in delivery.claim_rows
        or observable_id in delivery.goal_delivered
    )
    uncarried = tuple(
        observable_id
        for observable_id in observable_ids
        if observable_id not in carried
    )
    return "; ".join(
        part
        for part in (
            (
                "no claim in any cycle carries these declared observables: "
                + ", ".join(uncarried)
                if uncarried
                else ""
            ),
            (
                "a claim carries each of these declared observables without "
                "answering its declaration: " + ", ".join(carried)
                if carried
                else ""
            ),
        )
        if part
    )


def _inherited_verdict_reason(delivery: "_AnalysisDelivery") -> str:
    """Why numbers standing on another cycle's unanswered verdict wait."""

    return (
        "these delivered quantities stand on results the goal's own "
        "acceptance criterion rejected in another cycle, and no recorded "
        "decision cites that verdict: "
        + "; ".join(
            f"{', '.join(standing)} on {verdict.statement()} (receipt "
            f"{verdict.receipt_sha256s[-1][:8]})"
            for verdict, standing in delivery.inherited_unanswered
        )
    )


def _holding_verdicts(delivery: "_AnalysisDelivery") -> tuple[str, ...]:
    """Every unanswered verdict that holds this delivery open: those its
    own stream typed, then those another stream of the goal typed that
    its numbers stand on."""

    return tuple(
        dict.fromkeys(
            delivery.unanswered_verdicts
            + tuple(
                verdict.label for verdict, _ in delivery.inherited_unanswered
            )
        )
    )


def _held_quantity_ids(delivery: "_AnalysisDelivery") -> tuple[str, ...]:
    """The numbers a run's delivery holds on a criterion nobody answered:
    those standing on a verdict its own stream typed, and those standing
    on a verdict another stream of the goal typed."""

    return tuple(
        dict.fromkeys(
            delivery.stale_quantity_ids
            + tuple(
                quantity_id
                for _verdict, standing in delivery.inherited_unanswered
                for quantity_id in standing
            )
        )
    )


def _claims_a_later_refusal_supersedes(
    refusing: "_AnalysisDelivery",
    claimed_here: Collection[str],
    goal_delivered: Mapping[str, Mapping[str, Any]],
    declarations: Sequence[Mapping[str, Any]],
    verified: Collection[str] | None = None,
) -> dict[str, Mapping[str, Any]]:
    """Declared ids an earlier cycle claimed that this cycle refused, with
    the host's verification: the latest typed word about an id governs it.

    A claim made in an earlier cycle delivers an id the later cycles did
    not claim again (the goal-grain rule), and it used to deliver one a
    later cycle refused as well. R10 Q24's live goal g2r (CUHK 2153691)
    claimed dg-torsion-90deg in cycle 1 from a saddle search seeded at 90
    degrees; cycle 2 read that the search had reached the cis saddle, H-O-O-H
    0.047 degrees, renamed the number, and refused the observable through a
    blocked node the host verified, and its passed completion listed the id
    as delivered without. The goal settled achieved_with_observations
    saying "delivered in an earlier cycle: dg-torsion-90deg" -- the cis
    barrier, 7.95 kcal/mol, presented as the 90-degree free energy.

    Not asked: an id this cycle claimed itself, and an id declared with a
    required tolerance. A refusal of a claimed requirement is the
    sufficiency menu's route for a precision no admissible calculation
    reaches, and the claimed number stands (ax41 po3-r19 refused its
    0.5 kcal/mol, through a blocked node, over the value it delivered).
    An id declared with no tolerance has no precision to refuse, so its
    refusal can only be of the number itself.
    """

    allowed = set(
        refusing.verified_unreachable_ids if verified is None else verified
    )
    here = {str(item) for item in claimed_here}
    requirements = {
        str(row.get("observable_id") or "")
        for row in declarations
        if row.get("required_tolerance") is not None
    }
    superseded: dict[str, Mapping[str, Any]] = {}
    for observable_id in sorted(allowed):
        selector, _jobtype, blocked = refusing.unreachable_producers.get(
            observable_id, ("", "", "")
        )
        if (
            not (selector or blocked)
            or observable_id in here
            or observable_id in requirements
        ):
            continue
        row = goal_delivered.get(observable_id)
        if row is not None:
            superseded[observable_id] = row
    return superseded


def _superseded_claim_named(row: Mapping[str, Any]) -> str:
    """The clause a refusal adds about the earlier claim it supersedes."""

    value = row.get("value")
    number = (
        f" ({value} {row.get('unit') or ''})".replace(" )", ")")
        if value is not None
        else ""
    )
    return (
        "; this refusal supersedes the claim cycle "
        f"{row.get('cycle')} rendered under this id{number}"
    )


def _earlier_deliveries(
    delivery: "_AnalysisDelivery", superseded: Collection[str] = ()
) -> tuple[str, ...]:
    """The ids an earlier cycle delivered that this stream missed, less
    the ones this cycle's verified refusal superseded."""

    return tuple(
        item
        for item in delivery.delivered_in_earlier_cycles
        if item.split(" (delivered in cycle", 1)[0] not in set(superseded)
    )


def _what_the_delivery_carries(
    delivery: "_AnalysisDelivery",
    ledger_anomalies: Sequence[Mapping[str, Any]],
    evidence: Mapping[str, Any],
    superseded: Collection[str] = (),
) -> tuple[tuple[str, ...], dict[str, Any]]:
    """The lines and receipts a delivery carries beside its word.

    `_achieved_word` names them for an achieved goal -- the certificate
    it stands on, the observations the host recorded, the provenance of
    each number, the session's findings. A goal that settles on a typed
    refusal delivered the rest of its answer as well, and its word used
    to carry the refusal alone: 10 of 17 archived
    unreachable_from_evidence settlements left out anomaly receipts,
    falsified expectations or findings their own delivery held (R10
    Q24 census), which is what the owner's ruling that the first word
    must not hide what the run found was made against (2026-09-03).
    """

    word, lines = _achieved_word(
        delivery, ledger_anomalies, superseded=superseded
    )
    if word == "achieved_with_observations":
        evidence = _anomaly_evidence(
            evidence, ledger_anomalies, delivery, superseded
        )
    return lines, dict(evidence)


def _achieved_word(
    delivery: "_AnalysisDelivery",
    ledger_anomalies: Sequence[Mapping[str, Any]] = (),
    superseded: Collection[str] = (),
) -> tuple[str, tuple[str, ...]]:
    """The settlement word for a certified delivery.

    Plain ``achieved`` means the host saw nothing it could not explain;
    a delivery that carries host-recorded observations settles with the
    word that says so, because the one word a human reads first must
    not hide what the run found (owner ruling, 2026-09-03).
    """

    # A failed acceptance criterion the session answered is a finding the
    # delivery carries, and the word says so. The charter's words for
    # achieved_with_observations name "a pre-registered expectation the
    # physics left"; an answered criterion rides beside it whether the
    # session stated it before the engine ran or after it had read the
    # number -- the records read here do not say which. The settlement
    # names each from the records as they stand now; a completion's
    # listing, minted before a decision that answered it, describes that
    # earlier moment.
    # So is a pre-registered expectation the physics left: read as the
    # goal's latest score of each number it delivers, whichever cycle's
    # completion scored it, except a claim a later verified refusal
    # superseded -- the number that expectation judged is not delivered.
    carried = tuple(
        row
        for row in delivery.carried_expectations
        if str(row.get("observable_id") or "") not in set(superseded)
    )
    observed = tuple(
        sorted(
            {
                item
                for item in delivery.anomaly_output_ids
                if not str(item).startswith("failed_criterion:")
            }
            | set(_anomaly_output_ids_from_records(ledger_anomalies))
            | {
                verdict.observation_id
                for verdict, _standing in delivery.answered_criteria
            }
            | {
                f"falsified_expectation:{row.get('observable_id')}"
                for row in carried
            }
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
    if delivery.surface_mismatched_characterisations:
        provenance = provenance + (
            "characterised on another surface, so the geometry's own "
            "stationary point is still unchecked: "
            + ", ".join(
                f"{producer} by {consumer}"
                for producer, consumer in (
                    delivery.surface_mismatched_characterisations
                )
            ),
        )
    if delivery.surface_uncompared_characterisations:
        provenance = provenance + (
            "characterised by a Hessian whose surface the host could not "
            "compare with the geometry's: "
            + ", ".join(
                f"{producer} by {consumer}"
                for producer, consumer in (
                    delivery.surface_uncompared_characterisations
                )
            ),
        )
    earlier = _earlier_deliveries(delivery, superseded)
    if earlier:
        provenance = provenance + (
            "delivered in an earlier cycle: " + ", ".join(earlier),
        )
    if delivery.answered_criteria:
        provenance = provenance + _answered_criterion_reasons(
            delivery.answered_criteria
        )
    if delivery.claimed_answers:
        provenance = provenance + _claimed_answer_reasons(
            delivery.claimed_answers
        )
    if delivery.findings:
        provenance = provenance + _finding_reasons(delivery.findings)
    post_hoc = tuple(
        dict.fromkeys(
            str(row.get("observable_id") or "")
            for row in (*delivery.prediction_rows, *carried)
            if row.get("declared_after_evidence")
        )
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
    # "Certified" names a completion receipt that passed. A workflow that
    # declared nothing and carried no analysis chain has none, and the
    # sentence was written over it all the same.
    certified = (
        "the host completion gate certified the delivery"
        if _delivery_certified(delivery)
        else "no completion gate certified this delivery"
    )
    if observed:
        # What the host detected unasked and what the session itself stated
        # about its results are not one kind of thing: a criterion of the
        # session's own plan is not an observation "nobody asked for". Nor
        # is it always an expectation stated before the physics: L-S2 (R10
        # Q19, CUHK Slurm 2153514) planned its real->complex rule after it
        # had extracted -0.0288 Eh, and read its failed external rule as
        # the answer it had meant to test for. When a criterion was stated
        # is not read here, so the sentence says what was read: it did not
        # hold.
        stated = tuple(
            item
            for item in observed
            if item.startswith(("falsified_expectation:", "failed_criterion:"))
        )
        unasked = tuple(item for item in observed if item not in stated)
        lead = certified
        if unasked:
            lead += (
                "; the host also recorded observations nobody asked for: "
                + ", ".join(unasked)
            )
        if stated:
            lead += (
                "; criteria and predictions the session itself stated that "
                "did not hold: " + ", ".join(stated)
            )
        return ("achieved_with_observations", (lead,) + provenance)
    return ("achieved", (certified,) + provenance)


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

    settled, reasons, evidence = _delivery_settlement(
        ledger,
        goal_id=goal_id,
        events_path=events_path,
        terminal=terminal,
        workspace=workspace,
    )
    return _write_delivery_settlement(
        ledger,
        goal_id=goal_id,
        cycles=cycles,
        revisions_admitted=revisions_admitted,
        settled=settled,
        reasons=reasons,
        evidence=evidence,
        workspace=workspace,
    )


def _write_delivery_settlement(
    ledger: GoalLedger,
    *,
    goal_id: str,
    cycles: int,
    revisions_admitted: int,
    settled: str,
    reasons: tuple[str, ...],
    evidence: Mapping[str, Any],
    workspace: Path | None = None,
) -> GoalLoopResultV1:
    """Write a delivery's settlement and qualify what an achieved goal ran."""

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


def _delivery_settlement(
    ledger: GoalLedger,
    *,
    goal_id: str,
    events_path: Path,
    terminal: str,
    workspace: Path | None = None,
) -> tuple[str, tuple[str, ...], dict[str, Any]]:
    """The word, reasons and evidence one stream's delivery settles with.

    Computed and not written, so a turn that reads the delivery before
    the goal settles is handed the word the host would write without it
    and cannot change it.
    """

    delivery = _analysis_delivery(
        events_path,
        flagged_artifact_sha256s=_flagged_artifact_sha256s(
            _goal_anomalies(ledger)
        ),
        # The results an earlier run typed failed, so a number this
        # cycle claimed on one of them is worded as standing on it --
        # characterised or not -- exactly as it would be in that run's
        # own stream.
        failed_artifact_sha256s=(
            failed_artifacts(workspace) if workspace else ()
        ),
        uncharacterised_artifact_sha256s=(
            uncharacterised_artifacts(workspace) if workspace else ()
        ),
        goal_delivered_ids=_goal_delivered_ids(workspace, goal_id),
        declared_observables=_first_declarations(ledger),
        goal_findings=_goal_findings(workspace, goal_id),
        goal_streams=_goal_streams(ledger, workspace, goal_id),
        expectation_streams=_goal_streams(ledger, workspace, goal_id),
    )
    evidence = _settlement_evidence(delivery)
    # An earlier cycle's claim this cycle's verified refusal superseded:
    # the id is open again, and the refusal answers it.
    superseded = _claims_a_later_refusal_supersedes(
        delivery,
        delivery.claim_rows,
        _goal_delivered_ids(workspace, goal_id),
        _first_declarations(ledger),
    )
    # What the goal declared and has not delivered under its id in any
    # cycle, read from the ledger's first declarations and the record --
    # a refusal made from the in-session route has no completion and so
    # no limitation ids to key on (NOVEL-3 po3, 2026-09-05).
    open_ids = tuple(
        dict.fromkeys(
            _goal_open_declared_ids(
                ledger,
                workspace,
                goal_id,
                delivered_here=delivery.delivered_quantity_ids,
            )
            + tuple(delivery.undelivered_declared_ids)
            + tuple(superseded)
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
    certified = _delivery_certified(delivery)
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
        carried, evidence = _what_the_delivery_carries(
            delivery, _goal_anomalies(ledger), evidence, superseded
        )
        reasons = (
            "the completion receipt names required outputs "
            "delivered without: " + ", ".join(delivery.blocked_output_ids),
        ) + carried
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
        carried, evidence = _what_the_delivery_carries(
            delivery, _goal_anomalies(ledger), evidence, superseded
        )
        reasons = (
            _VERIFIED_REFUSAL_LEAD
            + "; ".join(
                f"{observable_id} -- "
                f"{delivery.unreachable_bases.get(observable_id, '')}"
                + (
                    _superseded_claim_named(superseded[observable_id])
                    if observable_id in superseded
                    else ""
                )
                for observable_id in refused_ids
            ),
        ) + carried
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
        # declared was never delivered by its id in any cycle. That is not
        # a delivery a scientist would sign, whatever the completion word
        # says.
        settled = "returned_to_human"
        reasons = (
            "the completion certified the chain, but "
            + _undelivered_declared_named(
                delivery, delivery.undelivered_declared_ids
            ),
        ) + delivery.open_declared_misses
        earlier = _earlier_deliveries(delivery, superseded)
        if earlier:
            reasons = reasons + (
                "delivered in an earlier cycle: " + ", ".join(earlier),
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
            "it: " + _unanswered_verdicts_named(delivery),
        )
    elif certified and delivery.inherited_unanswered:
        # The delivery stands on results the goal's own criterion rejected
        # in another cycle, and no recorded decision answered that
        # verdict: certifying it would say achieved over a finding nobody
        # read.
        settled = "returned_to_human"
        reasons = (_inherited_verdict_reason(delivery),)
    elif certified and delivery.claims:
        goal_anomalies = _goal_anomalies(ledger)
        settled, reasons = _achieved_word(
            delivery, goal_anomalies, superseded=superseded
        )
        if settled == "achieved_with_observations":
            evidence = _anomaly_evidence(
                evidence, goal_anomalies, delivery, superseded
            )
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
    elif (delivery.claims or delivery.decisions) and (
        delivery.unanswered_verdicts or delivery.inherited_unanswered
    ):
        # The gate did not certify because the plan's own criterion failed
        # and nobody answered it. Say which finding the goal is waiting
        # on: "the gate did not pass" named nothing a human could act on
        # (L1, R10 Q16).
        settled = "returned_to_human"
        reasons = (
            f"the session ended {terminal!r} ({delivery.ending}); the "
            "completion is not certified",
        )
        if delivery.unanswered_verdicts:
            reasons = reasons + (
                "a validation verdict failed and no recorded decision "
                "cites it: " + _unanswered_verdicts_named(delivery),
            )
        if delivery.inherited_unanswered:
            reasons = reasons + (_inherited_verdict_reason(delivery),)
    elif delivery.claims or delivery.decisions:
        # Something was recorded, but the host never certified
        # completion -- a human reads it, whatever the session's
        # terminal word was. The word says what it read: "the host
        # completion gate did not pass" was written over 25 archived
        # streams still settled this way, and 24 of them held no
        # completion receipt at all -- no gate had run; the 25th (L1,
        # R10 Q16) held a partial one whose findings the word dropped.
        settled = "returned_to_human"
        reasons = (
            f"the session ended {terminal!r} ({delivery.ending}); it "
            "recorded analysis and "
            + (
                "this stream holds no completion receipt, so nothing "
                "certifies it"
                if not delivery.completion_receipt_sha256
                else f"its completion receipt "
                f"{delivery.completion_receipt_sha256[:8]} is "
                f"{delivery.completion_status or 'unstated'}"
                + (
                    ", naming " + ", ".join(delivery.completion_findings)
                    if delivery.completion_findings
                    else ""
                )
                + (
                    "; limitations: "
                    + ", ".join(delivery.limitation_output_ids)
                    if delivery.limitation_output_ids
                    else ""
                )
            ),
        )
    else:
        settled = "returned_to_human"
        reasons = (f"the session ended {terminal!r}: {delivery.ending}",)
    return settled, tuple(reasons), dict(evidence)


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


def _session_run_id(session_result: Any) -> str:
    """The run directory a planning session wrote its stream into.

    A live session names it by ``session_id`` and carries no ``run_id``;
    reading only ``run_id`` -- the attribute the tests supplied -- left
    every live stream unnamed on the goal's spine, so all four goals of
    one campaign recorded no ``session_stream_recorded`` row and every
    wake resolved its stream by the newest-first glob that row exists to
    replace (pak campaign, 2026-09-19).
    """

    return str(
        getattr(session_result, "run_id", "")
        or getattr(session_result, "session_id", "")
        or ""
    )


def _session_events_path(session_result: Any, workspace: Path) -> Path:
    run_id = _session_run_id(session_result)
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
    streams: dict[int, str] = {}
    for entry in ledger.entries():
        if entry["kind"] == "run_recorded":
            reference = str(entry["payload"].get("run") or "")
        elif entry["kind"] == "analysis_evidence_recorded":
            reference = str(entry["payload"].get("evidence") or "")
        elif entry["kind"] == "session_stream_recorded":
            payload = entry["payload"]
            run_id = str(payload.get("run_id") or "")
            if run_id:
                streams[int(payload.get("cycle") or 0)] = run_id
        elif entry["kind"] == "rewake_opened":
            # A cycle that planned a workflow nobody could approve launched
            # no engine and, if it read nothing, recorded no analysis
            # evidence -- yet the host composed a typed report of how it
            # ended, and the one further wake embeds it. That report is
            # what the next plan answers, and the stream it describes is
            # the cycle's own: named here as analysis evidence is named, so
            # the wake and admission compare one reference the host wrote.
            # Trans-glyoxal (2026-09-19) was returned to the human for
            # never having read a run that did not exist. A transport
            # continuation carries no such report, and still names nothing.
            payload = entry["payload"]
            report = payload.get("failure_report") or {}
            stream = streams.get(int(payload.get("cycle") or 0), "")
            if (
                stream
                and not payload.get("transport_continuation")
                and isinstance(report, Mapping)
                and str(report.get("diagnosis") or "").strip()
            ):
                reference = f"runs/{stream}"
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
_FAILED_VALIDATION_CITATION = (
    " " + _WAKE_RULES["wake.failed_validation_receipt_answers_verdict"].text
)
_REFUSAL_AFFORDANCE = _WAKE_RULES["wake.refusal_is_a_deliverable"].text
_WORKSPACE_RECORD_RULE = " " + _WAKE_RULES["wake.workspace_record"].text
_EXCURSION_REPLICATION = (
    " " + _WAKE_RULES["wake.excursion_buys_replication"].text
)
_COHORT_EVIDENCE = " " + _WAKE_RULES["wake.cohort_evidence"].text
_READING_TURN = _WAKE_RULES["wake.reading_turn"].text + " "

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
        "A scan step failed to converge. Every step before it that "
        "converged is a point of the surface -- a constrained minimum at "
        "its held value, never a saddle: inspect_run on the result lists "
        "them with their energies, and bind_scan_point_geometry carries any "
        "of them forward. The step that failed is not a point, even where "
        "the program wrote a file for it, and bind_reached_geometry does "
        "not read a scan. Or loosen the step's own optimisation controls."
    ),
    "failed_nonconverged_excited_state": (
        "The response solver left a root unconverged, or the followed root "
        "fell below PySCF's positive-eigenvalue filter and vanished from "
        "the spectrum; the run outcome names which roots and the per-root "
        "flags stay inspectable on the result. Raise td_max_cycle in the "
        "project section (PySCF's default is 100), request fewer or more "
        "roots so the one you follow is well separated, or restart an "
        "excited-root optimisation from the reached geometry with "
        "bind_reached_geometry. Or the collapse is the finding: a root "
        "that meets the ground state or its neighbour is what a "
        "single-reference response cannot describe, and the recorded gap "
        "says so."
    ),
    "failed_nonconverged_correlation": (
        "The coupled-cluster amplitudes (or lambda equations) did not "
        "converge within cc_max_cycle; the SCF beneath them did. Raise "
        "cc_max_cycle in the project section (PySCF's default is 50), or "
        "check the reference: an amplitude set that will not settle often "
        "sits on a multireference case, where the spin diagnostic on the "
        "result is the finding rather than the iteration cap."
    ),
    "timeout_terminated": (
        "The engine ran out of the time the envelope granted. An "
        "optimisation moved and was cut off: restart from the geometry it "
        "reached -- bind_reached_geometry carries it forward -- inside the "
        "remaining budget. A relaxed scan keeps every step that converged "
        "before the clock, each a constrained minimum at its held value "
        "and none a saddle: inspect_run on the result lists them with "
        "their energies, and bind_scan_point_geometry carries any of them "
        "forward (bind_reached_geometry does not read a scan). Or reduce "
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
        "and the repair is the field the engine named. A relaxed scan keeps "
        "every step that converged before the error: inspect_run on the "
        "result lists them, and bind_scan_point_geometry carries any of "
        "them forward."
    ),
    "failed_result_validation": (
        "The program finished normally and the host's check of its result "
        "refused it; the native findings name the rule (a requested charge "
        "or method the output does not state, an output the reader cannot "
        "count as one result). A finding that names a setting is repaired "
        "in project YAML. Running the same input again gives the same "
        "output and the same refusal: when the output shows the program "
        "did what was asked, the refusal is the host's to answer, and the "
        "decision says so."
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
            # Cycle 1's plan is the goal's initial decision, not a
            # revision, so every revision is still a wake it can open.
            "wakes_after_this_cycle": int(max_revisions),
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
            + _FAILED_VALIDATION_CITATION
            + _OBSERVABLE_RESTATEMENT_ASK
            + _ADVERSARIAL_CLOSE
            + _REFUSAL_AFFORDANCE
            + _WORKSPACE_RECORD_RULE
        ),
    }


def _dispatch_ledger_keys() -> frozenset[str]:
    """Dispatch-receipt fields the goal ledger keeps, named once.

    ``scheduler_request`` carries what the scheduler was asked for beside
    what the envelope requested and the ceiling that bounded it: a clamp
    living only in a sidecar file is a clamp the goal record cannot be
    audited for, and the goal record is what a later process reads.

    ``wake_command`` is the only durable record of how this goal is meant
    to be resumed, and was dropped.

    ``wake_job_id`` names the job that will re-enter the goal when this
    cycle's elements are done -- a cohort's own elements deliberately
    wake nothing, because N tails would wake the model N times. It had no
    producer when this was written and the comment said so; the
    dispatcher populates it now (live: goal `butane-wave-2`, arrays
    2135242/2135266/2135285 each with their own wake job), so a reader of
    a parked goal can say which job is going to wake it.
    """

    return frozenset(
        {
            "scheduler",
            "job_id",
            "submitted_at",
            "submit_script",
            "wake_command",
            "wake_job_id",
            "scheduler_request",
        }
    )


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
        # The cores and memory the human approved. They were dropped from
        # this record, so a goal whose scheduler request disagreed with
        # its own review could not be audited from the ledger afterwards
        # -- which is how a live goal displayed 4 cores / 8 GB and was
        # allocated 32 / 160 with nothing on disk saying so.
        "cores": int(shown.get("cores") or 0),
        "memory_gb": float(shown.get("memory_gb") or 0.0),
    }


def _failed_criterion_record(
    verdict: FailedCriterionV1, minted_by: str = ""
) -> dict[str, Any]:
    """One failed criterion as a woken session is told it.

    The rule, the number it read against what it was held to, and the
    receipts that state the verdict, with the stream that minted them:
    a decision answers a verdict only by citing one, and three of three
    live sessions woken with the rule's name alone re-planned and
    re-evaluated their criteria to mint one to cite (o2r, L1, L-S2; the
    host had named the same receipt in its settlement all along).
    """

    return {
        "verdict": verdict.label,
        "statement": verdict.statement(),
        "receipt_sha256s": tuple(verdict.receipt_sha256s),
        **({"minted_by": minted_by} if minted_by else {}),
    }


def _deliverables_record(
    delivery: _AnalysisDelivery, previous_run: str = ""
) -> dict[str, Any]:
    """What the previous run's own stream says stands delivered.

    Names quantities and stated limitations, never values: the goal's
    demand is in the task, and this record lets a wake session see what
    it has already delivered, what the chain declared it could not, and
    what its own decisions doubt -- so the next action can follow the
    gap rather than the tool list. A failed criterion is the exception:
    it is named with the number it read and the receipts that state it,
    because answering it means citing one of them.
    """

    return {
        "delivered_quantity_ids": delivery.delivered_quantity_ids,
        "limitation_output_ids": delivery.limitation_output_ids,
        "doubted_quantity_ids": delivery.doubted_quantity_ids,
        "unanswered_failed_verdicts": tuple(
            _failed_criterion_record(verdict, previous_run)
            for verdict in delivery.unanswered_criteria
        ),
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


def _approved_bundle_digest_or_empty(approval_file: Any) -> str:
    """The bundle's declared digest, or empty when it has none.

    A local run's manifest binds to the same digest the scheduler path
    binds to, so an element and a local walk are admitted by one rule.
    A bundle written before that field existed simply binds to nothing,
    which is what it always did.
    """

    from chemsmart.agent.dispatch import _approved_bundle_digest

    try:
        return _approved_bundle_digest(approval_file)
    except ContractError:
        return ""


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
        if entry.get("goal_id") != goal_id:
            continue
        if entry.get("kind") == "answer":
            # A declared category a claim answered under its own id, as the
            # completion gate certified it; the shared predicate delivers
            # only a category by it, as it does a finding's answer.
            key = str(entry.get("claim_id") or "")
            previous = delivered.get(key)
            if key and (
                previous is None
                or int(entry.get("cycle") or 0)
                >= int(previous.get("cycle") or 0)
            ):
                delivered[key] = {
                    "cycle": int(entry.get("cycle") or 0),
                    "run": entry.get("run"),
                    "claim_id": key,
                    "claim_receipt_sha256": entry.get("claim_receipt_sha256"),
                    "answer": tuple(
                        dict(item) for item in entry.get("answer") or ()
                    ),
                }
            continue
        if entry.get("kind") == "finding":
            # A declared question the session answered with a finding;
            # the shared predicate delivers only a category by it.
            key = str(entry.get("claim_id") or "")
            previous = delivered.get(key)
            if key and (
                previous is None
                or int(entry.get("cycle") or 0)
                >= int(previous.get("cycle") or 0)
            ):
                delivered[key] = {
                    "cycle": int(entry.get("cycle") or 0),
                    "run": entry.get("run"),
                    "claim_id": key,
                    "finding_id": entry.get("finding_id"),
                    "statement": entry.get("statement"),
                    "finding_receipt_sha256": entry.get(
                        "finding_receipt_sha256"
                    ),
                    # The words the host read; the predicate delivers a
                    # category by them and by nothing else.
                    "answer": tuple(
                        dict(item) for item in entry.get("answer") or ()
                    ),
                }
            continue
        if entry.get("kind") != "claim":
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


def _goal_findings(
    workspace: Path | None, goal_id: str
) -> tuple[dict[str, Any], ...]:
    """Every finding this goal's sessions recorded, from the record.

    A finding written in one cycle is settled on in another: a session
    that reads a run, records what it found and plans the next stage
    hands the settlement a run stream that carries no decision. The
    record keeps each finding; the latest statement of an id stands and
    a supersession retires the id it names.
    """

    if workspace is None:
        return ()
    by_id: dict[str, dict[str, Any]] = {}
    for entry in read_workspace_record(workspace):
        if entry.get("kind") != "finding" or entry.get("goal_id") != goal_id:
            continue
        finding_id = str(entry.get("finding_id") or "")
        receipt = str(entry.get("finding_receipt_sha256") or "")
        if not finding_id or not receipt:
            continue
        by_id[finding_id] = {
            "finding_id": finding_id,
            "receipt_sha256": receipt,
            "statement": str(entry.get("statement") or ""),
            "answers_observable_id": str(entry.get("claim_id") or ""),
            "standing": str(entry.get("standing") or ""),
            "host_signals": tuple(entry.get("host_signals") or ()),
            "supersedes_finding_id": str(
                entry.get("supersedes_finding_id") or ""
            ),
            "answer": tuple(dict(item) for item in entry.get("answer") or ()),
            "cycle": int(entry.get("cycle") or 0),
        }
    return tuple(by_id.values())


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


def _goal_open_declared_ids(
    ledger: GoalLedger,
    workspace: Path | None,
    goal_id: str,
    delivered_here: Sequence[str] = (),
) -> tuple[str, ...]:
    """Declared ids the goal must deliver that no cycle delivered by id.

    The goal-grain join: the first declarations against every claim and
    answered category the record holds for this goal, plus what the
    stream being settled delivered. Both settlements ask this question;
    only the analysis-only one used to, so a run whose own stream listed
    nothing settled achieved over headlines no cycle had claimed.
    """

    delivered = set(_goal_delivered_ids(workspace, goal_id))
    delivered.update(str(item) for item in delivered_here)
    return tuple(
        observable_id
        for observable_id in _required_declared_ids(ledger)
        if observable_id and observable_id not in delivered
    )


def _stream_completion_status(events_path: Path) -> str:
    """The last completion status a stream recorded, or ``""``."""

    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ""
    status = ""
    for line in lines:
        if "analysis_completion_evaluated" not in line:
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") == "analysis_completion_evaluated":
            status = str((event.get("payload") or {}).get("status") or "")
    return status


def _session_wave_selections(
    events_path: Path | None,
) -> tuple[dict[str, Any], ...]:
    """Every wave a session selected, in stream order, as the host replied.

    The reply is the host's own record of the selection: its status, the
    workflow it names and the members it will submit.
    """

    if events_path is None:
        return ()
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    selections: list[dict[str, Any]] = []
    for line in lines:
        if "select_execution_wave" not in line:
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        payload = event.get("payload") or {}
        if (
            event.get("kind") != "tool_succeeded"
            or payload.get("tool") != "select_execution_wave"
        ):
            continue
        result = (payload.get("canonical_result") or {}).get("result") or {}
        selections.append(
            {
                "status": str(result.get("status") or ""),
                "workflow_id": str(result.get("workflow_id") or ""),
                "node_ids": tuple(
                    str(item) for item in result.get("node_ids") or ()
                ),
            }
        )
    return tuple(selections)


def _pending_decision_answered_in_text(events_path: Path | None) -> bool:
    """Whether the session was told once that a decision was pending, made
    no decision call after it, and ended on a turn with no tool call.

    Read from the session's own stream: the notice event the loop wrote,
    the calls that succeeded after it, and the last provider turn. It says
    how the session ended and nothing about what its text said -- the
    host reads no decision from text, and says so where it parks.
    """

    from chemsmart.agent.exposure import EXECUTION_DECISION_TOOLS

    if events_path is None:
        return False
    try:
        lines = events_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return False
    told = False
    decided = False
    last_turn_called: bool | None = None
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        kind = event.get("kind")
        payload = event.get("payload") or {}
        if kind == "execution_wave_decision_pending":
            told, decided = True, False
        elif kind == "tool_succeeded" and told:
            decided = (
                decided or payload.get("tool") in EXECUTION_DECISION_TOOLS
            )
        elif kind == "provider_turn_observed" and (
            "tool_calls_present" in payload
        ):
            last_turn_called = bool(payload.get("tool_calls_present"))
    return told and not decided and last_turn_called is False


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
    previous_run: str | None = None,
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
                    # Why a recovery opened when the previous run's own
                    # stream cannot say: a run with no analysis chain lists
                    # nothing undelivered, and a refusal re-read at
                    # settlement lives only on this row.
                    "uncertified",
                    "refusals_reread",
                    "undelivered_declared_observable_ids",
                    # What opened each earlier recovery. Without these a
                    # recovery row read {"cycle": 1} in the woken cycle's
                    # own trajectory (L-S2, R10 Q19), and a later cycle
                    # could not see which criterion or ending it was.
                    "verdicts",
                    "terminal_states",
                    "analysis_status",
                }
            },
        }
        for entry in ledger.entries()
        if entry["kind"]
        in {
            "run_recorded",
            "revision_admitted",
            "revision_returned",
            "recovery_opened",
        }
    )
    # The evidence a revision answers, whether an engine produced it or
    # an analysis-only cycle did. Gating this on `outcome is not None`
    # meant an analysis-only cycle embedded nothing, so admission
    # compared a real reference against an empty one and refused the
    # goal's first calculation -- A10 named the evidence and the wake
    # then declined to carry it. A reading turn names the stream it
    # reads, because a delivery a session made without recording a
    # decision leaves the ledger naming no evidence at all.
    if previous_run is None:
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
                # Whether a verdict is answered is a question about the
                # goal, and the settlement asks it of every stream the goal
                # holds. Read from the previous stream alone, G-h2's cycle-3
                # wake named val-real-stab/real-stable unanswered with the
                # receipt cycle 2 had minted by judging it again, while
                # cycle 2's decision had cited the run's receipt of the
                # same verdict and the settlement called it answered; the
                # woken session judged it twice more and cited that (R10
                # Q22, CUHK 2153627). Two organs, one function.
                goal_streams=_goal_streams(ledger, workspace, goal.goal_id),
            ),
            previous_run,
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
            # A woken cycle's plan is admitted as a revision, and a wake
            # after its run needs one more. g3-ethane and g5-methoxy
            # (2026-09-20) each selected a wave in their last cycle, were
            # told they would be woken when it ended, deferred their
            # claims to that wake, and were settled without it.
            "wakes_after_this_cycle": max(
                0, int(budgets.revisions_remaining) - 1
            ),
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
            + _FAILED_VALIDATION_CITATION
            + _CLAIM_BY_ID
            + _OBSERVABLE_RESTATEMENT_ASK
            + _ADVERSARIAL_CLOSE
            + _REFUSAL_AFFORDANCE
            + _WORKSPACE_RECORD_RULE
            + _COHORT_EVIDENCE
            + (
                _EXCURSION_REPLICATION
                if budgets.excursion_calls_remaining
                else ""
            )
        ),
    }


#: The two words a certified delivery settles with; only these are held
#: for a reading turn, because every other word already put a session
#: over the evidence or handed it to the human.
_CERTIFIED_WORDS = frozenset({"achieved", "achieved_with_observations"})


def _reading_context(
    goal: GoalRecordV1,
    ledger: GoalLedger,
    outcome: Any,
    *,
    workspace: Path | None,
    opened: Mapping[str, Any],
) -> dict[str, Any]:
    """The context of the one session that reads a certified delivery.

    The wake's own composition -- the typed outcome, the deliverables,
    the declarations, the anomalies, the workspace record -- over the
    stream the delivery came from, because a reading is judged against
    the same record a revision would be. What differs is what the turn
    may do: its budgets are zero, so a plan with a calculation node is
    refused where it is made, and it is told the word the host already
    computed, which nothing it records can change.
    """

    context = _wake_context(
        goal,
        ledger,
        outcome,
        workspace=workspace,
        previous_run=str(opened.get("run") or ""),
    )
    # Nothing ended in a way a revision answers, and no revision follows.
    for key in ("repair_menu", "repair_menu_dispositions", "failure_report"):
        context.pop(key, None)
    context["schema_version"] = "chemsmart.goal-reading-context.v1"
    context["budgets"] = {
        "binding_line": (
            "this reading turn spends nothing: it launches no engine and "
            "opens no revision"
        ),
        "engine_calls_remaining": 0,
        "excursion_calls_remaining": 0,
        "wall_seconds_remaining": 0.0,
        "revisions_remaining": 0,
        "wakes_after_this_cycle": 0,
    }
    context["reading"] = {
        "settlement_before_reading": {
            "state": str(opened.get("state") or ""),
            "reasons": tuple(opened.get("reasons") or ()),
        },
    }
    context["authority"] = _READING_TURN + _ADVERSARIAL_CLOSE
    return context


def _session_provider_cost(events_path: Path | None) -> dict[str, Any]:
    """What one session cost the provider, from its own stream.

    Every request the transport made is an ``api_attempt_observed`` row
    carrying the tokens the provider reported; the stream's first and
    last timestamps bound the session's wall time.
    """

    cost: dict[str, Any] = {
        "provider_requests": 0,
        "provider_turns": 0,
        "input_tokens": 0,
        "output_tokens": 0,
        "reasoning_tokens": 0,
        "stream_wall_seconds": 0.0,
    }
    if events_path is None:
        return cost
    try:
        lines = Path(events_path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return cost
    stamps: list[datetime] = []
    for line in lines:
        text = line.strip()
        if not text:
            continue
        try:
            event = json.loads(text)
        except json.JSONDecodeError:
            continue
        try:
            stamps.append(datetime.fromisoformat(str(event.get("timestamp"))))
        except (TypeError, ValueError):
            pass
        kind = str(event.get("kind") or "")
        payload = event.get("payload") or {}
        if kind == "api_attempt_observed":
            cost["provider_requests"] += 1
            for key in ("input_tokens", "output_tokens", "reasoning_tokens"):
                try:
                    cost[key] += int(payload.get(key) or 0)
                except (TypeError, ValueError):
                    pass
        elif kind == "provider_turn_observed":
            cost["provider_turns"] += 1
    if len(stamps) >= 2:
        cost["stream_wall_seconds"] = round(
            max(0.0, (max(stamps) - min(stamps)).total_seconds()), 3
        )
    return cost


def _achieved(execute_result: Any) -> bool:
    status = str(getattr(execute_result, "status", "") or "")
    analysis = str(getattr(execute_result, "analysis_status", "") or "")
    return status == "completed" and analysis in {"completed", ""}


def _delivery_certified(delivery: "_AnalysisDelivery") -> bool:
    """Did a completion gate certify this delivery? Its receipt says.

    The one answer every settlement path reads: the planning path, the
    run path and the achieved word's own sentence asked it three ways,
    and the run path never asked at all -- it signed achieved beside the
    sentence "no completion gate certified this delivery" (g2-hooh).
    """

    return delivery.completion_status == "passed"


def _analysis_nodes_run(events_path: Path) -> int:
    """How many analysis nodes a run's stream records as having run.

    A node the walk settled ``blocked_unsupported`` was declared
    non-executable intent and never ran; every other settlement --
    executed, failed, skipped -- is a node the chain reached.
    """

    count = 0
    try:
        lines = Path(events_path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return 0
    for line in lines:
        if '"workflow_analysis_node_settled"' not in line:
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue
        if event.get("kind") != "workflow_analysis_node_settled":
            continue
        state = str((event.get("payload") or {}).get("state") or "")
        if state and state != "blocked_unsupported":
            count += 1
    return count


@dataclass(frozen=True)
class _AnalysisDelivery:
    """What one session's durable stream says it delivered."""

    completion_status: str
    limitation_output_ids: tuple[str, ...]
    claims: int
    decisions: int
    receipt_sha256s: tuple[str, ...]
    #: The latest completion receipt of the stream and the findings it
    #: named, so a word about the gate says which receipt it read and
    #: what that receipt held -- or that the stream holds none.
    completion_receipt_sha256: str = ""
    completion_findings: tuple[str, ...] = ()
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
    #: Producer/consumer node pairs whose Hessian was computed on a
    #: different electronic surface from the geometry it consumed. The
    #: Hessian is a real number about a real structure; it just does not
    #: characterise the surface the producer walked on, and a reader who
    #: is not told that reads a ground-state spectrum as an excited
    #: minimum's (PySCF round 2 E2, 2026-09-13).
    surface_mismatched_characterisations: tuple[tuple[str, str], ...] = ()
    #: Producer/consumer node pairs credited although the two surfaces
    #: could not be compared -- a reader that records no surface, or a
    #: field it cannot determine. The credit stands; the comparison that
    #: was not made is said, because None is not agreement.
    surface_uncompared_characterisations: tuple[tuple[str, str], ...] = ()
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
    #: Verdicts that failed and that no recorded decision has cited, by
    #: ``node/rule``. A failed verdict is one of the plan's own acceptance
    #: criteria that did not hold on a result the host read -- a minimum
    #: that is a saddle, a reference that is unstable, a spread outside
    #: the tolerance the session stated. It never made a goal partial and
    #: it never opened a cycle, so a run could deliver every required
    #: output, display "failed" in its own report, and settle achieved
    #: with budget in hand. A decision that cites the validation receipt
    #: has answered it; one that does not has left it open.
    unanswered_verdicts: tuple[str, ...] = ()
    #: The same verdicts as the host read them -- the rule, the number it
    #: judged, what it was held against -- each carrying only the
    #: receipts of *this* stream that state it. The labels above are how
    #: the host counts; these are what a woken session is told, because
    #: a decision answers a verdict only by citing a receipt, and a wake
    #: that named the rule alone left every woken session re-minting one
    #: (R10 Q22: o2r, L1 and L-S2, 3 of 3).
    unanswered_criteria: tuple[FailedCriterionV1, ...] = ()
    #: The plan's own acceptance criteria that failed and that a recorded
    #: decision answered, each with the delivered quantities standing on
    #: the results it judged. The session read the finding and stands by
    #: its numbers; the delivery carries both, and the word says so.
    answered_criteria: tuple[
        tuple[FailedCriterionV1, tuple[str, ...]], ...
    ] = ()
    #: Failed criteria another stream of the goal typed and no recorded
    #: decision answered, each with the quantities this stream delivers
    #: from the results it rejected. Nothing standing on them is
    #: certified.
    inherited_unanswered: tuple[
        tuple[FailedCriterionV1, tuple[str, ...]], ...
    ] = ()
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
    #: Expectations another completion of the goal scored diverged, for a
    #: number the goal still delivers, which this stream's latest
    #: completion did not score again. Each row as it was scored, with
    #: the ``completion_receipt_sha256`` that scored it.
    carried_expectations: tuple[dict[str, Any], ...] = ()
    #: What the session did with each route the wake's repair menu
    #: offered, host-verified against that menu.
    route_dispositions: tuple[dict[str, Any], ...] = ()
    unreachable_bases: Mapping[str, str] = field(default_factory=dict)
    #: The (selector, jobtype, blocked node) each refusal named as the
    #: producer it needs.
    unreachable_producers: Mapping[str, tuple[str, str, str]] = field(
        default_factory=dict
    )
    #: What each claim of this stream carries, under both the names it
    #: answers to, so a current-cycle claim is judged in the dimension
    #: its declaration asked for exactly as a record row is.
    claim_rows: Mapping[str, Mapping[str, Any]] = field(default_factory=dict)
    #: Each claim of this stream as the pair of names it holds (claim_id,
    #: quantity_id), so a reader can tell the second name of a claim that
    #: carries a declared id from a claim under another name.
    claim_names: tuple[tuple[str, str], ...] = ()
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
    #: The session's standing findings in this stream: its sentence, the
    #: question it answers if any, the receipt the host minted over the
    #: relations it checked, and the host anomalies already under its
    #: evidence. A conclusion that reached no reader was prose.
    findings: tuple[Mapping[str, Any], ...] = ()
    #: Declared categories a claim answered under their own id, as the
    #: completion gate certified them: id -> the words the host read,
    #: each with its selector and receipt. A finding answering the same
    #: id speaks for it in ``findings`` instead.
    claimed_answers: Mapping[str, tuple[Mapping[str, Any], ...]] = field(
        default_factory=dict
    )

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

    for receipt, quantity_id in rejected_bindings:
        if receipt:
            # The results a rejected binding ultimately rests on, read by
            # the same resolver that gives a failed criterion its
            # identity; a source the records cannot resolve rejects
            # nothing.
            rejected_artifacts.update(
                result
                for result in results_read(
                    receipt,
                    quantity_id,
                    result_artifacts=artifact_by_receipt,
                    expression_sources=sources_by_output,
                )
                if not result.startswith("receipt:")
            )

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


def _printed_no_modes(record: Mapping[str, Any]) -> bool:
    """Whether a verified opt/ts result carries no frequency block.

    Read from the validator's own observations: every program's block
    records ``vibrational_mode_count``, and a count of zero on a job
    type that promises a stationary point means the promise was never
    checked, for a valid result as much as a failed one.
    """

    observations = record.get("observations") or {}
    jobtype = str(observations.get("jobtype") or record.get("jobtype") or "")
    if jobtype not in GEOMETRY_SEARCH_JOBTYPES:
        return False
    for key, value in observations.items():
        if isinstance(value, Mapping) and "vibrational_mode_count" in value:
            try:
                return int(value.get("vibrational_mode_count") or 0) == 0
            except (TypeError, ValueError):
                return False
    return False


@dataclass(frozen=True)
class _VerdictRecords:
    """What a stream holds that judges a failed acceptance criterion.

    The validation receipts, the digests recorded decisions cite, and the
    two maps that resolve a rule's inputs to the results it read. The
    settlement and the tool host build these from one representation --
    the receipt, carried whole in the stream's event -- and hand them to
    one function (``goal.failed_criteria``).
    """

    validations: tuple[Mapping[str, Any], ...] = ()
    cited: frozenset[str] = frozenset()
    result_artifacts: Mapping[str, str] = field(default_factory=dict)
    expression_sources: Mapping[tuple[str, str], tuple[str, ...]] = field(
        default_factory=dict
    )


def _stream_lines(path: Path) -> list[str]:
    try:
        return Path(path).read_text(encoding="utf-8").splitlines()
    except OSError:
        return []


def _carried_expectations(
    streams: Sequence[Path],
    events_path: Path,
    *,
    scored_here: Collection[str],
    retired: Collection[str],
) -> tuple[dict[str, Any], ...]:
    """Diverged expectations the goal scored outside this stream's latest
    completion, for numbers it still delivers.

    A pre-registered expectation the physics left is a result of the goal,
    and the settlement read it from one completion: the one its delivery
    stood on. A later cycle that re-rendered nothing scored nothing, so
    the earlier score of the number the goal still delivers never reached
    the word. R11 truth's census: CUHK r10/q7 g2-scan-modred settled plain
    achieved over cis-barrier 8.30 kcal/mol (declared 2-8) and
    oo160-torsion 141.5 deg (declared 100-135), which its cycle-1 run had
    scored diverged.

    Every completion of the goal's streams is read in the order the goal
    wrote them, this stream last; the latest row that scored a delivered
    claim is each id's score, so a later claim that agrees replaces an
    earlier divergence. Only an id the settling completion did not score
    (``scored_here``: the ids its own rows scored against a delivered
    claim) is carried -- what that completion scored reaches the word
    through its own receipt, as it always did -- and an id a later
    declaration retired is not scored at all.
    """

    if not streams:
        return ()
    here = Path(events_path).resolve()
    ordered = [path for path in streams if Path(path).resolve() != here] + [
        Path(events_path)
    ]
    latest: dict[str, dict[str, Any]] = {}
    for path in ordered:
        for line in _stream_lines(path):
            if '"analysis_completion_evaluated"' not in line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get("kind") != "analysis_completion_evaluated":
                continue
            payload = event.get("payload") or {}
            receipt = str(payload.get("receipt_sha256") or "")
            for row in payload.get("declared_observable_predictions") or ():
                if not isinstance(row, Mapping):
                    continue
                observable_id = str(row.get("observable_id") or "")
                if observable_id and row.get("delivered_claim_id"):
                    latest[observable_id] = {
                        **dict(row),
                        "completion_receipt_sha256": receipt,
                    }
    return tuple(
        row
        for observable_id, row in sorted(latest.items())
        if row.get("agreement") == "diverged"
        and observable_id not in set(retired)
        and observable_id not in set(scored_here)
    )


def _verdict_records(lines: Sequence[str]) -> _VerdictRecords:
    validations: list[Mapping[str, Any]] = []
    cited: set[str] = set()
    result_artifacts: dict[str, str] = {}
    expression_sources: dict[tuple[str, str], tuple[str, ...]] = {}
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
        record = payload.get("record") or {}
        if kind == "scientific_validation_evaluated":
            validations.append(payload)
        elif kind == "scientific_decision_recorded":
            # A decision that cites the validation receipt has looked at
            # the failed verdict and stood by its delivery. That is the
            # scientist's call to make, and citing the receipt is how it
            # is made without the host grading prose.
            cited |= cited_receipts(record.get("evidence_refs") or ())
        elif (
            kind in {"result_quantities_extracted", "thermochemistry_derived"}
            and digest
        ):
            # A thermochemistry derivation reads its result as surely as
            # an extraction does: a zero-point energy of a rejected saddle
            # stands on the verdict that rejected it.
            artifact = str(
                payload.get("artifact_sha256")
                or record.get("artifact_sha256")
                or ""
            )
            if artifact:
                result_artifacts[digest] = artifact
        elif kind == "quantity_expression_evaluated" and digest:
            for dependency in record.get("output_dependencies") or ():
                expression_sources[
                    (digest, str(dependency.get("output_id") or ""))
                ] = tuple(
                    str(item)
                    for item in dependency.get("source_receipt_sha256s") or ()
                )
    return _VerdictRecords(
        validations=tuple(validations),
        cited=frozenset(cited),
        result_artifacts=result_artifacts,
        expression_sources=expression_sources,
    )


def _merge_verdict_records(*parts: _VerdictRecords) -> _VerdictRecords:
    """The goal's records, in stream order, each receipt once."""

    validations: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    cited: set[str] = set()
    result_artifacts: dict[str, str] = {}
    expression_sources: dict[tuple[str, str], tuple[str, ...]] = {}
    for part in parts:
        for item in part.validations:
            digest = str(item.get("receipt_sha256") or "")
            if digest and digest in seen:
                continue
            seen.add(digest)
            validations.append(item)
        cited |= part.cited
        result_artifacts.update(part.result_artifacts)
        expression_sources.update(part.expression_sources)
    return _VerdictRecords(
        validations=tuple(validations),
        cited=frozenset(cited),
        result_artifacts=result_artifacts,
        expression_sources=expression_sources,
    )


def _goal_streams(
    ledger: GoalLedger, workspace: Path | None, goal_id: str
) -> tuple[Path, ...]:
    """Every stream the goal's own spine names: planning sessions and runs.

    A spine names a session's stream by its session row, and by the
    analysis evidence a cycle recorded when it read results and decided
    (``.chemsmart-agent``-relative, as the wake resolves it). Every goal
    written before the session row existed (2026-09-17) names its sessions
    only the second way, and so does a session that raised. Read by the
    first way alone, ax41 goal-ino3-r17's word left out three expectations
    its cycle-3 session had scored diverged, for numbers it still delivered
    (R11 truth-4). The first mention of a stream keeps its place.
    """

    if workspace is None:
        return ()
    agent = Path(workspace) / ".chemsmart-agent"
    streams: list[Path] = []
    for entry in ledger.entries():
        payload = entry.get("payload") or {}
        if entry["kind"] == "session_stream_recorded":
            run_id = str(payload.get("run_id") or "")
            if run_id:
                streams.append(agent / "runs" / run_id / "events.jsonl")
        elif entry["kind"] == "run_recorded":
            run = str(payload.get("run") or "")
            if run:
                streams.append(agent / Path(*run.split("/")) / "events.jsonl")
        elif entry["kind"] == "analysis_evidence_recorded":
            evidence = str(payload.get("evidence") or "")
            if evidence:
                streams.append(
                    agent / Path(*evidence.split("/")) / "events.jsonl"
                )
    return tuple(dict.fromkeys(path for path in streams if path.is_file()))


def goal_verdict_records(
    goal_directory: str | Path, *, excluding: str | Path | None = None
) -> _VerdictRecords:
    """The verdict records of every stream a goal's ledger names.

    A completion certificate asks of its goal what the settlement asks of
    it: which failed acceptance criteria stand, and whether a recorded
    decision answered each. The settlement reads every stream this ledger
    names; a host that mints a certificate reads the same streams through
    this, less its own (``excluding``), whose records it holds itself. A
    certificate that read its own host alone was partial over a verdict a
    decision in another stream had answered, and the goal returned to the
    human (R11 truth-3, item 1).
    """

    directory = Path(goal_directory)
    workspace = directory.parent.parent.parent
    skip = Path(excluding).resolve() if excluding is not None else None
    return _merge_verdict_records(
        *(
            _verdict_records(_stream_lines(path))
            for path in _goal_streams(
                GoalLedger(directory), workspace, directory.name
            )
            if skip is None or path.resolve() != skip
        )
    )


def _analysis_delivery(
    events_path: Path,
    *,
    goal_delivered_ids: Mapping[str, Mapping[str, Any]] | None = None,
    declared_observables: Sequence[Mapping[str, Any]] = (),
    uncharacterised_artifact_sha256s: tuple[str, ...] = (),
    flagged_artifact_sha256s: Sequence[str] = (),
    failed_artifact_sha256s: Sequence[str] = (),
    inherited_unreachable: Mapping[str, str] = {},
    goal_findings: Sequence[Mapping[str, Any]] = (),
    goal_streams: Sequence[Path] = (),
    expectation_streams: Sequence[Path] = (),
) -> _AnalysisDelivery:
    """Read the delivery facts a settlement stands on.

    ``goal_streams`` are the goal's other streams -- its earlier cycles'
    planning sessions and runs. A settlement reads the plan's own failed
    acceptance criteria across all of them: a verdict one cycle's run was
    typed with is answered by the decision a later session records, and a
    number this stream delivers from a result an earlier verdict rejected
    stands on that verdict whichever stream typed it.

    ``expectation_streams`` are the goal's streams whose completions
    scored its pre-registered expectations; a settlement passes them so
    the word carries the goal's latest score of every number it delivers
    (``_carried_expectations``).

    Every field is a typed record the host itself wrote: the
    completion receipt with its stated limitations, the claim and
    decision records, and the receipt digests a settlement cites. The
    first live goal round's classifier read none of these -- it
    watched validation rules, the one place an honest refusal leaves
    no trace -- and settled a receipts-backed refusal as achieved.
    """

    completion_status = ""
    completion_receipt = ""
    completion_findings: tuple[str, ...] = ()
    categorical_answers: dict[str, tuple[dict[str, Any], ...]] = {}
    prediction_rows: tuple[dict[str, Any], ...] = ()
    route_dispositions: tuple[dict[str, Any], ...] = ()

    declared_misses: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
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
    # The producer each refusal named, so a settlement can read the
    # results again for it once a run has written them.
    unreachable_producers: dict[str, tuple[str, str, str]] = {}
    receipts: list[str] = []
    doubt_refs: set[str] = set()
    claim_pairs: list[tuple[str, str]] = []
    # "<receipt>:<quantity_id>" of every quantity a claim carries as its
    # uncertainty: rendered on that claim, where every reader of it sees it.
    uncertainty_references: set[str] = set()
    claim_rows: dict[str, dict[str, Any]] = {}
    claim_names: list[tuple[str, str]] = []
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
    node_surfaces: dict[str, Mapping[str, Any]] = {}
    stopped_by: list[str] = []
    workflows_planned = 0
    nodes_previewed = 0
    last_plan_refusal = ""
    terminal_reason = ""
    anomaly_ids: tuple[str, ...] = ()
    # The goal's earlier findings first, from the record; this stream's
    # own statement of an id comes after and stands.
    findings_by_id: dict[str, dict[str, Any]] = {
        str(row.get("finding_id") or ""): {
            "finding_id": str(row.get("finding_id") or ""),
            "receipt_sha256": str(row.get("receipt_sha256") or ""),
            "statement": str(row.get("statement") or ""),
            "answers_observable_id": str(
                row.get("answers_observable_id") or ""
            ),
            "standing": str(row.get("standing") or ""),
            "host_signals": tuple(row.get("host_signals") or ()),
            "supersedes_finding_id": str(
                row.get("supersedes_finding_id") or ""
            ),
            "answer": tuple(dict(item) for item in row.get("answer") or ()),
        }
        for row in goal_findings
        if row.get("finding_id") and row.get("receipt_sha256")
    }
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
            # The session's findings, each a receipt the host minted over
            # relations it checked. The latest statement of an id is the
            # standing one; a supersession retires the id it names.
            for item in payload.get("findings") or ():
                if not isinstance(item, Mapping):
                    continue
                finding_id = str(item.get("finding_id") or "")
                digest_of = str(item.get("receipt_sha256") or "")
                if not finding_id or not digest_of:
                    continue
                findings_by_id[finding_id] = {
                    "finding_id": finding_id,
                    "receipt_sha256": digest_of,
                    "statement": str(item.get("statement") or ""),
                    "answers_observable_id": str(
                        item.get("answers_observable_id") or ""
                    ),
                    "standing": str(item.get("standing") or ""),
                    "host_signals": tuple(
                        str(signal)
                        for signal in item.get("host_signals") or ()
                    ),
                    "supersedes_finding_id": str(
                        item.get("supersedes_finding_id") or ""
                    ),
                    "answer": tuple(
                        dict(word)
                        for word in item.get("answer") or ()
                        if isinstance(word, Mapping)
                    ),
                }
            for item in payload.get("unreachable_observables") or ():
                observable_id = str(item.get("observable_id") or "")
                if not observable_id:
                    continue
                unreachable_bases[observable_id] = (
                    f"{item.get('statement') or ''} [{item.get('basis') or ''}]"
                )
                unreachable_producers[observable_id] = (
                    str(item.get("selector") or ""),
                    str(item.get("jobtype") or ""),
                    str(item.get("blocked_node_id") or ""),
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
                reference = str(claim.get("uncertainty_reference") or "")
                if reference:
                    uncertainty_references.add(reference)
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
                claim_names.append(
                    (claim_id, str(claim.get("quantity_id") or ""))
                )
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
                surface = recorded_surface(record)
                if surface:
                    node_surfaces[node_name] = surface
                # A result characterises the structure it was handed only
                # if it kept that structure: a re-optimisation prints the
                # modes of the structure it reached. 12 archived credits
                # went to producers whose consumer had moved (an ORCA
                # scan's refined well, an opt re-optimised).
                if (
                    str(record.get("state") or "") == "valid"
                    and printed_modes(record)
                    and str(record.get("jobtype") or "")
                    not in STRUCTURE_MOVING_STAGES
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
            completion_receipt = digest
            # What the gate certified for each declared category, the last
            # completion's word: a claim under the category's id answers
            # it exactly as a finding does (R10 Q23).
            categorical_answers = {
                str(observable_id): tuple(
                    dict(word)
                    for word in words or ()
                    if isinstance(word, Mapping)
                )
                for observable_id, words in (
                    payload.get("declared_categorical_answers") or {}
                ).items()
            }
            completion_findings = tuple(
                str(item)
                for item in (
                    (payload.get("record") or {}).get("findings") or ()
                )
            )
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
    # The plan's own acceptance criteria, judged by the one function the
    # session's completion and the executor's completion also call. A
    # verdict is answered when a recorded decision cites a receipt that
    # states it -- in this stream or in any other stream of the goal, since
    # a woken session answers the verdict its previous run was typed with.
    here = _verdict_records(lines)
    goal = _merge_verdict_records(
        here, *(_verdict_records(_stream_lines(path)) for path in goal_streams)
    )
    verdicts = failed_criteria(
        goal.validations,
        cited=goal.cited,
        result_artifacts=goal.result_artifacts,
        expression_sources=goal.expression_sources,
    )
    here_receipts = {
        str(item.get("receipt_sha256") or "") for item in here.validations
    }
    unanswered_criteria = tuple(
        replace(
            verdict,
            receipt_sha256s=tuple(
                receipt
                for receipt in verdict.receipt_sha256s
                if receipt in here_receipts
            ),
        )
        for verdict in verdicts
        if not verdict.answered
        and here_receipts.intersection(verdict.receipt_sha256s)
    )
    unanswered = tuple(verdict.label for verdict in unanswered_criteria)
    # The rule read these receipts and rejected what it found in them;
    # everything else computed from the same receipts describes the same
    # rejected structure -- unless a recorded decision answered the
    # verdict, in which case the session has read the finding and stands
    # by what it computed, and the numbers are its delivery, carrying it.
    rejected_bindings.extend(
        binding
        for verdict in verdicts
        if not verdict.answered
        and here_receipts.intersection(verdict.receipt_sha256s)
        for binding in verdict.bindings
        if binding[0]
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
    # ...and only when the Hessian was computed on the producer's own
    # surface. A ground-state Hessian on the geometry an excited-root
    # optimisation reached describes a different potential energy
    # surface: it is a real number about a real structure and it says
    # nothing about whether that structure is a minimum of the surface
    # the optimisation walked on. Where the two cannot be compared --
    # a program that records no surface, a field its reader cannot
    # determine -- the Hessian still characterises, as it did before,
    # and the run story says the comparison was not available.
    surface_mismatches: list[tuple[str, str]] = []
    # ``surfaces_agree`` answers None when it cannot compare, and its own
    # contract says a caller must not read None as agreement. The credit
    # stands, as it did before, and the comparison that was not made is
    # recorded beside it rather than living only in this comment: 36 of
    # 37 archived credits were given on None (every ORCA, Gaussian and
    # xTB pair; any PySCF result written before surfaces were recorded).
    surface_uncompared: list[tuple[str, str]] = []
    for consumer, producer in handoffs.items():
        if consumer not in characterising:
            continue
        verdict = surfaces_agree(
            node_surfaces.get(producer), node_surfaces.get(consumer)
        )
        if verdict is False:
            surface_mismatches.append((producer, consumer))
            continue
        if verdict is None:
            surface_uncompared.append((producer, consumer))
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
    stale, _rejected_artifacts = _stale_quantity_ids(
        claim_pairs=claim_pairs,
        rejected_bindings=rejected_bindings,
        expression_outputs=expression_outputs,
        # The verdict join reads every receipt that read a result,
        # thermochemistry derivations included.
        artifact_by_receipt=here.result_artifacts,
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
    # The same walk once per failed criterion, over the whole goal's
    # records: which numbers this stream delivers stand on the results
    # that criterion judged. An answered verdict rides the word with the
    # numbers it carries. An unanswered one typed in another stream holds
    # every number of this stream standing on it: the planning path read
    # only its own stream, so a woken session could deliver from results
    # an earlier run's criterion had rejected, rewrite the criterion
    # without the rules that failed, and settle achieved (h1b, ax41
    # general round, 2026-09-02) -- while the run path, which carries
    # rejections across cycles, would have held the same numbers stale.
    goal_expressions = tuple(
        (receipt, output_id, tuple(sources))
        for (receipt, output_id), sources in goal.expression_sources.items()
    )
    answered_criteria: list[tuple[FailedCriterionV1, tuple[str, ...]]] = []
    inherited_unanswered: list[tuple[FailedCriterionV1, tuple[str, ...]]] = []
    for verdict in verdicts:
        standing, _judged = _stale_quantity_ids(
            claim_pairs=claim_pairs,
            rejected_bindings=verdict.bindings,
            expression_outputs=goal_expressions,
            artifact_by_receipt=goal.result_artifacts,
        )
        typed_here = bool(here_receipts.intersection(verdict.receipt_sha256s))
        if verdict.answered and (standing or typed_here):
            answered_criteria.append((verdict, standing))
        elif not verdict.answered and standing and not typed_here:
            inherited_unanswered.append((verdict, standing))
    # An expression's exported outputs, as opposed to the intermediate
    # node_values it computed on the way: the receipt contract pins
    # output_dependencies' ids to outputs' quantity_ids, in order, so the
    # lineage already collected names exactly what was exported.
    claimed_ids = {
        quantity_id for _receipt, quantity_id in claim_pairs if quantity_id
    }
    # An output a delivered claim carries as its uncertainty was rendered:
    # it is on the claim. Counting it "computed and never rendered" held
    # two goals whose completions had passed -- r10/q3 g2 (the D0 spread,
    # 0.50 kJ/mol and 41.8 cm-1, on both headline claims) and r10/q9 g1
    # (the BDE's 6.0 kJ/mol) -- and opened a recovery, a revision spent on
    # nothing, in eight more.
    exported_output_ids = {
        output_id
        for digest, output_id, _sources in expression_outputs
        if output_id and f"{digest}:{output_id}" not in uncertainty_references
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
    # A finding answers the declared question it names, and one the task
    # did not ask for is an observation the word carries under its own
    # prefix -- the session's, on relations the host checked. Its receipt
    # is the settlement's evidence like any other the stream minted.
    retired_findings = {
        row["supersedes_finding_id"]
        for row in findings_by_id.values()
        if row["supersedes_finding_id"]
    }
    standing_findings = tuple(
        row
        for finding_id, row in findings_by_id.items()
        if finding_id not in retired_findings
    )
    answered_ids: set[str] = set()
    # A declared category a claim answered under its own id, as the gate
    # certified it; a finding answering the same id below speaks for it.
    claimed_answers = {
        observable_id: words
        for observable_id, words in categorical_answers.items()
        if words and not any(word.get("finding_id") for word in words)
    }
    for observable_id, words in claimed_answers.items():
        answered_ids.add(observable_id)
        claim_rows[observable_id] = {
            "answer": words,
            "claim_record_sha256": str(
                words[0].get("claim_record_sha256") or ""
            ),
        }
    for row in standing_findings:
        receipts.append(row["receipt_sha256"])
        if row["answers_observable_id"] and row["answer"]:
            # Delivered only by the words the host read; a finding recorded
            # without them (before the answer was bound to a word) answers
            # nothing.
            answered_ids.add(row["answers_observable_id"])
            claim_rows[row["answers_observable_id"]] = {
                "finding_receipt_sha256": row["receipt_sha256"],
                "finding_id": row["finding_id"],
                "answer": row["answer"],
            }
        # A finding never joins the observations the word names. The word
        # is the host's: what its sensors detected and what the physics
        # made of a prediction written before it. Four development
        # sessions (2026-09-24) typed a process remark -- "the same pair
        # reads 3.296 A in the other isomer, so the observable
        # distinguishes them" -- as a finding nobody asked for, in both
        # arms of a matched pair, and the word said the run had seen
        # something. The session's findings ride the reasons and the
        # evidence under every word instead, as its own.
    return _AnalysisDelivery(
        findings=standing_findings,
        claimed_answers={
            observable_id: words
            for observable_id, words in claimed_answers.items()
            if not any(
                row["answers_observable_id"] == observable_id and row["answer"]
                for row in standing_findings
            )
        },
        ending=ending,
        terminal_reason=terminal_reason,
        sufficiency=carried + tuple(sufficiency_rows),
        unclaimed_output_ids=tuple(sorted(exported_output_ids - claimed_ids)),
        stopped_by=tuple(stopped_by),
        anomaly_output_ids=anomaly_ids,
        unanswered_verdicts=unanswered,
        unanswered_criteria=unanswered_criteria,
        answered_criteria=tuple(answered_criteria),
        inherited_unanswered=tuple(inherited_unanswered),
        completion_status=completion_status,
        completion_receipt_sha256=completion_receipt,
        completion_findings=completion_findings,
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
        claim_names=tuple(dict.fromkeys(claim_names)),
        delivered_quantity_ids=tuple(
            sorted(
                {
                    quantity_id
                    for _receipt, quantity_id in claim_pairs
                    if quantity_id and quantity_id not in stale
                }
                | answered_ids
            )
        ),
        stale_quantity_ids=stale,
        flagged_quantity_ids=flagged,
        claims_rendered=bool(claims),
        failed_source_quantity_ids=failed_quantities,
        characterised_source_quantity_ids=characterised_quantities,
        characterised_artifact_sha256s=tuple(sorted(characterised)),
        uncharacterised_source_quantity_ids=uncharacterised_quantities,
        surface_mismatched_characterisations=tuple(surface_mismatches),
        surface_uncompared_characterisations=tuple(surface_uncompared),
        decision_uncertainties=tuple(decision_uncertainties),
        declared_observable_misses=declared_misses,
        goal_delivered=dict(goal_delivered_ids or {}),
        declared_observables=tuple(
            dict(item) for item in declared_observables
        ),
        verified_unreachable_ids=tuple(sorted(verified_unreachable)),
        prediction_rows=prediction_rows,
        carried_expectations=_carried_expectations(
            expectation_streams,
            events_path,
            scored_here={
                str(row.get("observable_id") or "")
                for row in prediction_rows
                if row.get("delivered_claim_id")
            },
            retired=superseded_observable_ids(declared_observables),
        ),
        route_dispositions=route_dispositions,
        unverified_unreachable_ids=tuple(
            sorted(unverified_unreachable - verified_unreachable)
        ),
        unreachable_bases=dict(unreachable_bases),
        unreachable_producers=dict(unreachable_producers),
    )


def _settlement_evidence(delivery: _AnalysisDelivery) -> dict[str, Any]:
    """Receipts a settlement cites, from the session's own stream."""

    evidence: dict[str, Any] = {}
    if delivery.decisions and delivery.receipt_sha256s:
        evidence = {
            "scientific_decisions": delivery.decisions,
            "receipt_sha256s": delivery.receipt_sha256s,
        }
        if delivery.decision_uncertainties:
            evidence["decision_uncertainties"] = (
                delivery.decision_uncertainties
            )
    elif delivery.anomaly_output_ids and delivery.receipt_sha256s:
        # An observation brings its own receipts. The executor's stream
        # never carries a decision, and a diverged pre-registration --
        # the session's expectation, scored by the completion receipt
        # in this stream -- raised the word that settles on receipts
        # and handed it none: r9 g5 and r8 goal-irc2 delivered every
        # declared observable and returned to the human with a
        # contract error in place of the delivery.
        evidence = {"receipt_sha256s": delivery.receipt_sha256s}
    if delivery.findings:
        # Under every word, not only the one that names them: a
        # conclusion the session bound to receipts is part of what the
        # goal delivered whatever else the settlement says -- including
        # a settlement read off an executor's stream, which holds none.
        evidence["findings"] = tuple(dict(row) for row in delivery.findings)
    return evidence


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
    "read",
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
        sealed: bool = True,
        initial_decision: str = "approve",
        stop_file: str | Path | None = None,
        session_kwargs: Mapping[str, Any] | None = None,
        manual_execution_intent: bool = False,
        reading_turn: bool = False,
        _resuming: bool = False,
    ) -> None:
        if dispatch not in DISPATCH_MODES:
            raise ContractError(
                f"unsupported dispatch mode {dispatch!r}; "
                f"choose one of {DISPATCH_MODES}"
            )
        self.task = task
        self.workspace = Path(workspace).resolve()
        # A goal id is a public identifier like every id the host mints
        # from it; normalising it here, through the one rule, is what keeps
        # the approval the driver names and the approval the bundle checks
        # one string (PySCF round 2 E1, 2026-09-13).
        self.goal_id = require_identifier(goal_id, "goal_id")
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
        # Sealed is the default for a scheduler dispatch: the memory
        # ceiling is the profile's maximum less the declared headroom.
        # An unsealed run holds the request to the profile alone, which
        # is the escape for a profile too small to leave anything above
        # the headroom.
        self.sealed = bool(sealed)
        if dispatch_run is None and dispatch == "scheduler":
            from functools import partial

            from chemsmart.agent.dispatch import dispatch_run_to_scheduler

            dispatch_run = partial(
                dispatch_run_to_scheduler, server=server, sealed=sealed
            )
        self.dispatch_run = dispatch_run
        self.initial_decision = initial_decision
        self.stop_file = stop_file
        self.session_kwargs = dict(session_kwargs or {})
        # `from_review` is the one human-direct review surface. It names a
        # human execution decision, not an Agent's omitted wave, and retains
        # its deliberately serial/manual execution behavior below.
        self.manual_execution_intent = bool(manual_execution_intent)
        #: Host policy: whether a certified delivery is read by one
        #: session before the goal settles. The reading launches nothing,
        #: admits no revision and cannot change the word the host had
        #: already computed; it adds what the session found in the
        #: results to the settlement beside that word. Nothing before
        #: the settlement reads this, so a goal's cycles are the same
        #: with it on or off up to the word the host would have written.
        self.reading_turn = bool(reading_turn)

        self.goal_dir = (
            self.workspace / ".chemsmart-agent" / "goals" / self.goal_id
        )
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
                        "sealed": bool(sealed),
                        "server": server,
                        "stop_file": _resolved_or_none(stop_file),
                        "reading_turn": bool(reading_turn),
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
                    "cores": self.envelope.resources.cores,
                    "memory_gb": self.envelope.resources.memory_gb,
                }
            )

        self.goal: GoalRecordV1 | None = None
        self.outcome: Any = None
        #: The failure report a re-woken cycle carries: what the previous
        #: cycle left undone, said as gate, diagnosis, route and cost.
        self.failure_report: dict[str, Any] | None = None
        self.cycles = 0
        self.revisions_admitted = 0
        # A rejection is a verdict about bytes and does not expire: the
        # settlement reads every verdict of the goal's streams again when
        # it signs, with every decision that answered one, so a result an
        # earlier cycle's criterion rejected holds a later claim on it
        # until a decision cites that verdict (R11 truth). The *standing*
        # delivery is whatever the most recent claim-rendering cycle
        # held, because that is what the goal currently answers with.
        # Keying this on quantity ids instead would ask a later cycle to
        # reuse an earlier cycle's names: one live recovery re-derived a
        # torsion correctly under a new id and would have been held open
        # forever over a number it had already replaced.
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
        self._resuming_started_run = False
        #: A reading a previous process opened and never recorded: the
        #: resumed driver settles the held word without reading again.
        self._reading_interrupted = False

    # -- public surface ---------------------------------------------------

    @classmethod
    def resume(
        cls, *, workspace: str | Path, goal_id: str, **kwargs: Any
    ) -> "GoalDriver":
        """Rebuild a parked goal's driver from its ledger, at ``outcome``.

        Admitted only when the ledger's last run was dispatched and never
        recorded, and the goal is not settled: the same one human
        decision continues in its own run directory, and nothing here
        creates a second one. A reading turn a previous process opened
        and never recorded resumes at ``read`` and settles the word it
        held without reading again.
        """

        goal_id = require_identifier(goal_id, "goal_id")
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
            "sealed",
            "server",
            "stop_file",
            "reading_turn",
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
        pending_decisions = [
            entry
            for entry in entries
            if entry["kind"] == "execution_wave_decision_pending"
            and int(entry["payload"].get("cycle", 0)) not in dispatched_cycles
            and int(entry["payload"].get("cycle", 0))
            not in {
                int(item["payload"].get("cycle", 0)) for item in interrupted
            }
        ]
        # A reading turn holds a settlement the host had already computed;
        # a process that died inside it left that word unwritten. The
        # delivery is not the reading's to lose, so the resumed driver
        # writes the held word and reads nothing.
        read_cycles = {
            int(item["payload"].get("cycle", 0))
            for item in entries
            if item["kind"] == "reading_recorded"
        }
        held_readings = [
            item
            for item in entries
            if item["kind"] == "reading_opened"
            and int(item["payload"].get("cycle", 0)) not in read_cycles
        ]
        if (
            held_readings
            and not parked
            and not interrupted
            and not pending_decisions
        ):
            driver.goal = driver.ledger.load()
            driver.cycles = int(held_readings[-1]["payload"]["cycle"])
            driver.revisions_admitted = sum(
                1 for item in entries if item["kind"] == "revision_admitted"
            )
            driver._reading_interrupted = True
            driver.phase = "read"
            return driver
        if not parked and not interrupted and not pending_decisions:
            raise ContractError(
                f"goal {goal_id!r} has no parked, interrupted, or pending "
                "execution decision to resume"
            )
        entry = (parked or interrupted or pending_decisions)[-1]
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
        driver._restore_standing_delivery(entries)
        if entry["kind"] == "run_started":
            # An interrupted local run: re-enter the execute phase with the
            # same bundle and run directory; the executor's continuation
            # replays finished nodes and names the interrupted one.
            driver.bundle_file = Path(str(entry["payload"]["approval_file"]))
            driver._resuming_started_run = True
            driver.phase = "execute"
            return driver
        if entry["kind"] == "execution_wave_decision_pending":
            from chemsmart.agent.cohort import (
                execution_wave_decision_from_record,
            )

            decision = execution_wave_decision_from_record(
                dict(entry["payload"].get("execution_wave_decision") or {})
            )
            driver.bundle_file = Path(
                str(entry["payload"].get("approval_file") or "")
            )
            review_file = str(entry["payload"].get("review_file") or "")
            driver.review_file = Path(review_file) if review_file else None
            driver.session = SimpleNamespace(
                execution_wave_decision=decision,
                selected_execution_wave=tuple(decision.node_ids),
            )
            driver.phase = "execution_decision_pending"
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
        driver.manual_execution_intent = True
        driver.phase = "decide"
        return driver

    def _restore_standing_delivery(
        self, entries: Sequence[Mapping[str, Any]]
    ) -> None:
        """The numbers the recorded runs among ``entries`` left held.

        Replayed in order and read as the settlement reads them -- against
        every verdict and every decision of the goal's streams -- so a
        rejection made in cycle 1 still holds a number cycle 3 claimed
        from the same result, unless a decision has answered it. One
        function, because a resumed goal and a replayed settlement must
        rebuild the same state: a census that restated this loop broke on
        the first tree that changed it (R11 truth).
        """

        streams = _goal_streams(self.ledger, self.workspace, self.goal_id)
        for item in entries:
            if item["kind"] != "run_recorded":
                continue
            delivery = _analysis_delivery(
                self.workspace
                / ".chemsmart-agent"
                / Path(*str(item["payload"].get("run") or "").split("/"))
                / "events.jsonl",
                goal_streams=streams,
            )
            if delivery.claims_rendered:
                self.standing_stale = _held_quantity_ids(delivery)

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
            "execution_decision_pending": (
                self._execution_decision_pending_resume
            ),
            "outcome": self._outcome,
            "settle": self._settle,
            "read": self._read,
        }[phase]
        try:
            handler()
        except ContractError as exc:
            # A typed refusal raised inside a phase but outside that
            # handler's own net ended the goal unsettled: E1 (PySCF round
            # 2, 2026-09-13) died in decide before the ledger existed --
            # the bundle refused the driver's own approval id -- and the
            # process exited with no goal_created and no settlement.
            # Every ending is a settlement. A ContractError is an outcome
            # the human reads; a non-contract exception stays a crash,
            # because a genuine defect must not be laundered into one.
            if self.phase in {"settled", "parked"}:
                raise
            if self.goal is None:
                self.goal = self._goal_record(
                    identity="",
                    conditions={"solvents": (), "thermochemistry": ()},
                    review_sha256="",
                )
                self.ledger.create(self.goal)
                self._flush_declarations()
            self._typed_error(f"{phase} phase", exc)
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

    def _record_session_stream(self, session: Any) -> None:
        """Name the stream this cycle planned in, on the goal's own spine.

        The settling path resolves a missing `events_path` by taking the
        newest `live-*` stream in the workspace. That is right for the
        case it was written for -- a session that raised seconds ago in
        this very process, where the guess cannot be wrong -- and wrong
        for a woken cycle, where the window is minutes to hours *by
        design*, because the point of `--dispatch scheduler` is that the
        scientist does something else while the array runs. A session id
        carries a task-spec digest and not a goal id, so two goals in one
        workspace are indistinguishable by name and the glob has no goal
        filter: goal A's wake would project goal B's declarations and
        claims into goal A's ledger.

        One workspace with two goals is the ordinary way a scientist
        works on one system, so the guess is recorded away rather than
        improved.
        """

        run_id = _session_run_id(session)
        if not run_id:
            return
        self._defer_or_append(
            "session_stream_recorded",
            {"cycle": self.cycles, "run_id": run_id},
        )

    def _named_its_own_stream(self) -> bool:
        """Whether this cycle recorded which stream it planned in.

        The distinction `_planned_events_path`'s `None` cannot carry:
        "no row for this cycle", which is every goal written before the
        row existed and must keep the old fallback, against "a row whose
        stream has since gone", which must not fall back to guessing.
        A goal that named its stream never reaches the workspace-wide
        glob -- that is the whole content of naming it.
        """

        return any(
            entry.get("kind") == "session_stream_recorded"
            and int((entry.get("payload") or {}).get("cycle") or 0)
            == self.cycles
            for entry in self.ledger.entries()
        )

    def _planned_events_path(self) -> Path | None:
        """The stream this cycle planned in, from the goal's own record.

        ``None`` when this cycle recorded none -- every goal written
        before this did -- or when the stream it named is gone, which is
        absence rather than licence to substitute another goal's.
        """

        run_id = ""
        for entry in reversed(self.ledger.entries()):
            if entry.get("kind") != "session_stream_recorded":
                continue
            payload = entry.get("payload") or {}
            if int(payload.get("cycle") or 0) != self.cycles:
                continue
            run_id = str(payload.get("run_id") or "")
            break
        if not run_id:
            return None
        candidate = (
            self.workspace
            / ".chemsmart-agent"
            / "runs"
            / run_id
            / "events.jsonl"
        )
        return candidate if candidate.is_file() else None

    def _planned_streams(self) -> tuple[Path, ...]:
        """Every stream this goal planned in, from its own record.

        A reached geometry bound in one cycle's session can be consumed
        by any later cycle's plan, and the lineage receipt lives only in
        the stream that bound it.
        """

        found = []
        for entry in self.ledger.entries():
            if entry.get("kind") != "session_stream_recorded":
                continue
            run_id = str((entry.get("payload") or {}).get("run_id") or "")
            candidate = (
                self.workspace
                / ".chemsmart-agent"
                / "runs"
                / run_id
                / "events.jsonl"
            )
            if run_id and candidate.is_file() and candidate not in found:
                found.append(candidate)
        return tuple(found)

    def _review_file_for_cycle(self) -> Path | None:
        """This cycle's displayed review, wherever the driver came from.

        `_plan` sets `self.review_file`, and a parked cycle's workspace
        record is written by the *wake* -- which `resume` rebuilds at the
        outcome phase without planning. So every scheduler-dispatched
        cycle recorded its results with `review_file=None`, and
        `record_run` reads the level of every result from exactly that
        file: seven live result rows went to disk with an empty
        `level_sha256` while the review beside them carried the settings
        text and its digest for every node. The per-claim level
        attribution this round built had nothing to attribute, and
        nothing failed.

        The path is deterministic; `_plan` computes the same one.
        """

        if self.review_file is not None and Path(self.review_file).is_file():
            return Path(self.review_file)
        candidate = self.goal_dir / "reviews" / f"cycle-{self.cycles}.json"
        return candidate if candidate.is_file() else None

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

    def _reviewed_node_ids(self) -> tuple[str, ...]:
        """Every node the cycle's displayed review carried.

        Empty for a review packet that names none, which is absence and
        not a reason to drop a wave.
        """

        review_file = self.goal_dir / "reviews" / f"cycle-{self.cycles}.json"
        try:
            review = _review_record(review_file)
        except (OSError, json.JSONDecodeError):
            return ()
        packet = review.get("workflow_execution_review") or review
        return tuple(
            str(row.get("node_id") or "")
            for row in packet.get("node_reviews") or ()
            if isinstance(row, Mapping) and row.get("node_id")
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
                review_file=self._review_file_for_cycle(),
                lineage_events_paths=self._planned_streams(),
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

    def _settle_delivery(
        self, terminal: str, extra_reasons: tuple[str, ...] = ()
    ) -> GoalLoopResultV1 | None:
        if self.events_path is None:
            # Nothing durable to read a delivery from: the goal still ends
            # in a typed state, and the reason says why a human reads it.
            #
            # This arm is unreachable from the two callers that exist.
            # `_plan` sets `events_path` a few lines before it calls here,
            # and `_outcome`'s analysis-only branch assigns it
            # immediately before -- so a woken cycle, which `resume`
            # rebuilds without planning, never arrives with None. It is
            # *absent* from that path rather than rescued on it: nothing
            # resolves a missing stream here the way
            # `_project_before_settling` does for `_typed_error`, and
            # that was deliberate, because a lookup for a state no caller
            # can produce is speculative code that would also mask the
            # signal a third caller ought to trip.
            #
            # A third caller is therefore what this arm is for, and it is
            # the worst kind of dead branch to leave unannotated: it does
            # not crash, it *ends the goal* -- returned_to_human, with a
            # reason naming a cause that would be wrong. Anyone adding
            # one should give this arm a real answer first.
            reason = (
                f"session terminal state: {terminal}; the planning session "
                "left no event stream"
            )
            self.ledger.settle("returned_to_human", reasons=(reason,))
            return self._settled("returned_to_human", (reason,))
        settled, reasons, evidence = _delivery_settlement(
            self.ledger,
            goal_id=self.goal_id,
            events_path=self.events_path,
            terminal=terminal,
            workspace=self.workspace,
        )
        reasons = tuple(reasons) + tuple(extra_reasons)
        try:
            stream = str(
                self.events_path.parent.relative_to(
                    self.workspace / ".chemsmart-agent"
                )
            )
        except ValueError:
            stream = ""
        if self._open_reading(
            path="delivery",
            run=stream,
            state=settled,
            reasons=reasons,
            result_reasons=reasons,
            evidence=evidence,
        ):
            return None
        self.result = _write_delivery_settlement(
            self.ledger,
            goal_id=self.goal_id,
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            settled=settled,
            reasons=reasons,
            evidence=evidence,
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
            events_path = self._planned_events_path()
            if events_path is not None:
                self.events_path = events_path
            elif self._named_its_own_stream():
                # Named, and gone. Absence, not licence to substitute
                # another goal's -- which is what the glob below does.
                return
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
        # engine wall time's evidence the same way. Only this cycle's
        # run: an error raised while planning finds the previous cycle's
        # run directory here, and projecting it re-recorded that run
        # under the label of a run this cycle never made (o2r).
        if (
            self.run_directory is not None
            and self.run_directory
            == self.goal_dir / "runs" / f"cycle-{self.cycles}"
        ):
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
        # The stream this cycle plans in is not known until its session
        # returns. Until then the previous cycle's stream is the previous
        # cycle's, and a session that raises is projected from whatever
        # this names: o2r's cycle 2 (R10 Q13, CUHK 2152079) raised after
        # claiming its whole answer, and the typed-error projection read
        # cycle 1's stream, so none of cycle 2's claims or findings
        # reached the workspace record.
        self.events_path = None
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
        # Named on the goal's own spine, so the wake that settles this
        # cycle hours later resolves the stream it planned in rather than
        # whatever is newest in a workspace it may be sharing.
        self._record_session_stream(self.session)
        self._record_declarations()
        self._record_dispositions()
        self._record_approaches()
        self._record_input_checks()
        if terminal == "waiting_for_approval" and self.events_path is not None:
            # A session that read results, claimed, found, and then
            # planned the next stage delivered those rows as surely as one
            # that stopped. Only a stopping session's stream was projected,
            # so ino3-r12's cycle-2 claims (26 of them) reached no record
            # and a finding written before a further run could not reach
            # the settlement that run ends in.
            self._record_workspace(self.events_path, "")
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
            # And across the goal's streams, as the settlement reads it: a
            # number standing on a verdict another stream's decision
            # answered is delivered there and was stale here, so this
            # reader re-woke a goal for an observable the settlement
            # calls delivered (R10 Q22).
            goal_streams=_goal_streams(
                self.ledger, self.workspace, self.goal_id
            ),
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
        # A claim that carries a declared id is not a claim "under another
        # name", and neither is the second name it holds. G-h2's re-wake
        # listed verdict-bs, verdict-complex and verdict-real as claims
        # under other names: they were the quantity ids of the three claims
        # that carried its three declared category questions (R10 Q22).
        carrying = {
            name
            for pair in delivery.claim_names
            if set(pair).intersection(declared)
            for name in pair
        }
        rendered = tuple(
            item
            for item in delivery.delivered_quantity_ids
            if item not in declared and item not in carrying
        )
        if not undelivered:
            self._open_requirement_rewake(terminal, delivery, unresolved)
            return True
        # What the completion receipt says each one lacks, when it says
        # it: the host held the sentence that answers the diagnosis ("a
        # category is answered by a word the host read, bound by a
        # finding") and the report withheld it.
        misses = tuple(
            text
            for text in delivery.declared_observable_misses
            if any(f"'{item}'" in text for item in undelivered)
        )
        diagnosis = (
            f"the previous cycle ended {terminal!r} ({delivery.ending}); "
            + _undelivered_declared_named(delivery, undelivered)
            + (
                "; its completion receipt says: " + " | ".join(misses)
                if misses
                else ""
            )
            + (
                "; it rendered claims under other names: "
                + ", ".join(rendered)
                if rendered
                else (
                    ""
                    if delivery.claims_rendered
                    else "; it rendered no claim"
                )
            )
            + (
                "; no completion was certified"
                if not _delivery_certified(delivery)
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
            # Every ending the invariant admits, each with what it costs.
            # The route named only the first two and the cost said "no
            # engine call", while R10 Q15 g1 had been woken by a review the
            # host refused over one unpreviewed node, with its whole grant
            # in hand: the ending that answers such a diagnosis is the plan
            # itself, repaired.
            "route": (
                "claim each undelivered id from receipts in hand -- "
                "extract_result_quantities, derive_thermochemistry, "
                "evaluate_quantity_expression, record_analysis_claims with "
                "claim_id set to the declared id, in the declared unit; a "
                "category is answered the same way, by the word (or the "
                "integer) the host read from the program's output claimed "
                "with claim_id set to the declared id -- or "
                "plan_scientific_workflow with no calculation_nodes, which "
                "the host executes when planned; or record_scientific_"
                "decision naming each id that cannot be delivered, with its "
                "required producer and the receipts that show it; or, where "
                "the calculation that produces them has not run, end with "
                "an executable plan for review: repair what the diagnosis "
                "names and give every initial node a green preview with "
                "compile_command, and the host builds the review under the "
                "goal's standing decision"
            ),
            "cost": (
                "claiming, an analysis-only plan and a refusal cost no "
                "engine call; an executable plan spends engine calls from "
                "what remains. This re-wake is charged one revision and is "
                "the only re-wake this goal is granted"
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
                "re-wake is charged one revision and is the only re-wake "
                "this goal is granted"
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
        decision = self._execution_wave_decision()
        if (
            decision is None
            and not self.manual_execution_intent
            and not self._resuming_started_run
        ):
            from chemsmart.agent.cohort import build_execution_wave_decision

            decision = build_execution_wave_decision()
        if decision is not None and decision.state != "selected":
            if (
                str(getattr(decision, "workflow_id", "") or "")
                and not tuple(getattr(decision, "ready_node_ids", ()) or ())
                and self.events_path is not None
                and not _session_wave_selections(self.events_path)
            ):
                # A park waits for a decision, and none can be made here:
                # the workflow the session planned last holds nothing to
                # select, the session selected no wave on any workflow,
                # and a resumed pending decision stays pending. R10 Q15 g2
                # (CUHK 2153334) parked this way for good -- a diagnostic
                # probe had been reviewed, the session then planned an
                # analysis-only refusal and recorded a refusal the host
                # verified, and nothing ever read it. The goal settles on
                # what the session's stream holds, and says which reviewed
                # workflow did not run and why. A session that did select
                # a wave a later plan replaced still parks naming it: that
                # is the wave's boundary, not a missing decision.
                self._settle_delivery(
                    str(getattr(self.session, "terminal_state", "") or ""),
                    extra_reasons=(
                        "the reviewed workflow was not run: the session "
                        "selected no wave, and the workflow it planned last "
                        f"({decision.workflow_id or '(unnamed)'}) holds no "
                        "calculation a wave could select",
                    ),
                )
                return
            self._park_for_execution_wave_decision(
                decision,
                reason=(
                    "the Agent explicitly continued scientific reasoning"
                    if decision.state == "continue_reasoning"
                    else self._undecided_boundary_reason(decision)
                ),
            )
            return
        wave = self._dispatchable_wave()
        if decision is not None and not wave:
            # A selected member can fall outside the final approval after a
            # re-plan or a non-executable review finding.  That is evidence
            # to return to the Agent, never permission to reinterpret an
            # explicit cohort as the old unbounded serial path.
            self._park_for_execution_wave_decision(
                decision,
                reason=(
                    "the selected wave is not covered by this approval; "
                    "the host did not substitute a serial execution"
                ),
            )
            return
        self.run_directory = self.goal_dir / "runs" / f"cycle-{self.cycles}"
        self.run_directory.mkdir(parents=True, exist_ok=True)
        run_reference = f"goals/{self.goal_id}/runs/cycle-{self.cycles}"
        if self.dispatch_run is not None:
            # The claim precedes the submission, because the submission
            # is the irreversible act. `run_dispatched` is keyed too, but
            # its key was taken *after* sbatch: two attempts on one cycle
            # both reached the scheduler, got different job ids, and only
            # then did the conflicting payload raise -- with a second job
            # already queued, two jobs writing one run directory, and the
            # ledger and the receipt sidecar naming different jobs.
            claimed = self.ledger.append(
                "run_dispatch_claimed",
                {"cycle": self.cycles, "run": run_reference},
                idempotency_key=(
                    f"run-dispatch-claimed:{self.goal_id}:{self.cycles}"
                ),
            )
            if not claimed:
                existing = next(
                    (
                        entry
                        for entry in reversed(self.ledger.entries())
                        if entry.get("kind") == "run_dispatched"
                        and int((entry.get("payload") or {}).get("cycle") or 0)
                        == self.cycles
                    ),
                    None,
                )
                abandoned = any(
                    entry.get("kind") == "run_dispatch_abandoned"
                    and int((entry.get("payload") or {}).get("cycle") or 0)
                    == self.cycles
                    for entry in self.ledger.entries()
                )
                if existing is None and not abandoned:
                    # Claimed and never dispatched: a controller killed
                    # inside the submission window. Whether a job exists
                    # is not decidable from here, and submitting again is
                    # how one grant buys two jobs. This is the state
                    # CHEMSMART already has a word for.
                    self._typed_error(
                        "scheduler dispatch",
                        ContractError(
                            f"cycle {self.cycles} was claimed for dispatch "
                            "and no job was recorded: the submission is "
                            "ambiguous and pending human reconciliation. "
                            "Check the scheduler for a job naming "
                            f"{run_reference} before resubmitting."
                        ),
                    )
                    return
                if existing is None:
                    # Claimed, abandoned, and recorded as abandoned:
                    # this cycle may be dispatched again, because
                    # nothing reached a scheduler.
                    self.ledger.append(
                        "run_dispatch_reclaimed",
                        {"cycle": self.cycles, "run": run_reference},
                    )
                else:
                    payload = dict(existing.get("payload") or {})
                    self.result = GoalLoopResultV1(
                        goal_id=self.goal_id,
                        settlement="parked",
                        cycles=self.cycles,
                        revisions_admitted=self.revisions_admitted,
                        reasons=(
                            f"cycle {self.cycles} is already submitted "
                            f"as {payload.get('scheduler', 'scheduler')} "
                            f"job {payload.get('job_id', '?')}; resume "
                            "with chemsmart agent wake --goal "
                            f"{self.goal_id}",
                        ),
                    )
                    self.phase = "parked"
                    return
            try:
                receipt = self.dispatch_run(
                    approval_file=self.bundle_file,
                    workspace=self.workspace,
                    run_directory=self.run_directory,
                    goal_id=self.goal_id,
                    cycle=self.cycles,
                    # The driver has held the approved resources all
                    # along; the scheduler request is now made from them
                    # rather than from the operator's profile alone.
                    resources=getattr(
                        getattr(self, "envelope", None), "resources", None
                    ),
                    envelope=getattr(self, "envelope", None),
                    # The explicitly selected wave, in the order the Agent
                    # selected it.  `_execute` above rejects an undecided
                    # boundary before this irreversible scheduler call.
                    cohort_node_ids=wave,
                )
            except (ContractError, ValueError, OSError) as exc:
                # A dispatch that could not happen is not an ambiguous
                # submission. The claim is written before `sbatch` so two
                # attempts cannot both reach the scheduler, and a claimed
                # cycle with no job is reported as pending human
                # reconciliation -- right for a controller killed inside
                # the submission window, and wrong for a submission that
                # provably never left this process. `_require_array_support`
                # raises `ValueError` and `parse_submission` raises
                # `ProbeUnitError(ValueError)`, and catching only
                # `ContractError` left the goal unsettled with a claim on
                # the ledger, so the next invocation sent a human looking
                # for a job that was never submitted.
                # Which side of the irreversible act this failed on is
                # a fact on disk, not an inference from the exception
                # type: the dispatcher writes the receipt before the
                # first `sbatch` with no job id and rewrites it the
                # moment the scheduler names one. Everything after that
                # -- the wake script, the wake submission, the final
                # receipt -- can fail with an array already queued, and
                # recording *that* as abandoned would tell a reader
                # nothing reached the scheduler while the allocation
                # burned, and would disable the one branch that sends
                # them to look for it.
                from chemsmart.agent.dispatch import read_dispatch_receipt

                partial = read_dispatch_receipt(self.run_directory)
                submitted = str(getattr(partial, "job_id", "") or "")
                self.ledger.append(
                    (
                        "run_dispatch_submitted"
                        if submitted
                        else "run_dispatch_abandoned"
                    ),
                    {
                        "cycle": self.cycles,
                        "run": run_reference,
                        "reason": str(exc),
                        **({"job_id": submitted} if submitted else {}),
                    },
                )
                self._typed_error("scheduler dispatch", exc)
                return
            self.dispatch_receipt = receipt
            payload = {
                "cycle": self.cycles,
                "run": run_reference,
                **{
                    key: value
                    for key, value in _record_of(receipt).items()
                    if key in _dispatch_ledger_keys()
                },
            }
            # One cycle is dispatched once. A retried dispatch that
            # reached the scheduler twice would park the goal on the
            # second job and orphan the first.
            self.ledger.append(
                "run_dispatched",
                payload,
                idempotency_key=(
                    f"run-dispatched:{self.goal_id}:{self.cycles}"
                ),
            )
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
        # A wave the session selected bounds a local run too. The stem
        # promises every session, without qualification, that a wave is
        # what runs and then it reasons again; on the local path -- the
        # default, and what the terminal interface uses -- no manifest
        # was written, so the walk was unbounded and the executor ran
        # straight past the barrier. Concurrency was never the promise
        # and one process cannot offer it; the barrier is the promise,
        # and one process keeps it exactly. Written before the executor
        # is handed the directory, for the same reason the scheduler
        # path writes it before submitting.
        local_wave = wave
        if local_wave:
            from chemsmart.agent.cohort import build_cohort_manifest

            build_cohort_manifest(
                goal_id=self.goal_id,
                cycle=self.cycles,
                bundle_sha256=_approved_bundle_digest_or_empty(
                    self.bundle_file
                ),
                node_ids=local_wave,
                max_concurrent_tasks=1,
                created_at=_utc_now(),
            ).write(self.run_directory)
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
            recorded = self.ledger.append(
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
                # The same cycle by another route: recorded once, however
                # it is recorded.
                idempotency_key=(f"run-recorded:{self.goal_id}:{self.cycles}"),
            )
            if not recorded:
                self._concede(run_reference)
                return
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
        # The charge is keyed on the goal, the cycle and the bundle that
        # ran, so a second wake on one parked cycle -- a duplicate
        # scheduler notification, a human running `agent wake` twice, a
        # dependent wake job that fires beside a tail -- records the same
        # fact once instead of subtracting its engine calls from the
        # grant again. `resume` filtering already-recorded cycles is a
        # read-then-act; this is the guarantee.
        if not self.ledger.append(
            "run_recorded",
            payload,
            idempotency_key=(f"run-recorded:{self.goal_id}:{self.cycles}"),
        ):
            # Another process recorded this cycle. The key deduplicated
            # the row, and for years that was read as "nothing to do" and
            # returned silently -- so the loser carried straight on into
            # workspace recording, settlement and, for a repairable
            # ending, a second planning turn over the same evidence. One
            # cohort is one wake is one turn: the row that charges the
            # cycle is also what hands over the continuation, and this
            # process did not win it.
            self._concede(run_reference)
            return
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

    def _concede(self, run_reference: str) -> None:
        """Stand down: another process owns this cycle's continuation.

        Exactly-once is about the *turn*, not the row. Two wakes on one
        parked cycle -- a duplicate scheduler notification, a dependent
        wake job firing beside a tail, a human running the public
        command twice -- both pass ``resume`` before either writes, and
        only the keyed append can tell them apart. The loser stops here
        with the evidence intact; the winner records the workspace,
        settles, and takes the one model turn the cohort earned.

        It projects the workspace record first. Standing down closed the
        duplicated turn and opened a smaller hole: if the *winner* dies
        between writing the row and writing the record, nothing projects
        that cycle's receipts and `resume` computes the parked set as
        dispatched-minus-recorded, which is now empty, so the goal is
        unresumable. `record_run` is durably idempotent, so doing it here
        costs the winner nothing and the charter's sentence holds --
        evidence that survives on disk and reaches no projection is
        unreachable to every later reader, which is the same as lost.
        """

        if self.run_directory is not None:
            events_path = (
                self.events_path
                if getattr(self, "events_path", None) is not None
                else self.run_directory / "events.jsonl"
            )
            self._record_workspace(events_path, run_reference)
        self.result = GoalLoopResultV1(
            goal_id=self.goal_id,
            settlement="parked",
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            reasons=(
                f"cycle {self.cycles} was recorded by another process; "
                f"its continuation owns {run_reference}",
            ),
        )
        self.phase = "parked"

    def _dispatchable_wave(self) -> tuple[str, ...]:
        """The selected wave, less anything this approval cannot launch.

        The model selects against the *planned* draft's ready frontier,
        which is what exists at selection time -- before any review is
        resolved. A stage the review then retained as non-executable was
        displayed as intent and never approved, so it can appear in a
        wave and can never run.

        The failure was silent and total: every element's cohort scope
        would be a node nothing can launch, so each runs nothing, the run
        outcome derivation finds no run, and the analysis-only branch
        settles a partition that never ran as a complete delivery.

        This is the host checking two things it owns against each other,
        not a refusal shown to the Agent; what it drops it records, and
        an empty result is a pending execution decision rather than an array
        of nothing or an unbounded serial execution.
        """

        decision = self._execution_wave_decision()
        selected = (
            tuple(str(item) for item in decision.node_ids)
            if decision is not None and decision.state == "selected"
            else ()
        )
        if not selected:
            return ()
        # Subtracting the retained set is not intersecting with the
        # approved one: a member in neither -- an id from a workflow the
        # session re-planned away from, and the selection carries no
        # workflow id -- went straight through.
        retained = set(self._declared_non_executable_ids())
        reviewed = set(self._reviewed_node_ids())
        approved = (reviewed - retained) if reviewed else None

        def _covered(node_id: str) -> bool:
            if node_id in retained:
                return False
            return approved is None or node_id in approved

        wave = tuple(item for item in selected if _covered(item))
        # Two causes, and one fixed sentence told a reader the wrong one
        # for half of them. A retained stage was displayed and refused
        # execution; an id the review never carried is a selection made
        # against a workflow this approval is not for, which is the
        # stale-selection case a re-plan produces -- and a re-plan is how
        # the live sequence's second wave was made.
        dropped = tuple(
            {
                "node_id": item,
                "reason": (
                    "the displayed review retained this as "
                    "non-executable intent, so this approval cannot "
                    "launch it"
                    if item in retained
                    else "this approval's review does not carry this "
                    "calculation, so the selection was made against "
                    "another workflow"
                ),
            }
            for item in selected
            if not _covered(item)
        )
        if dropped:
            self.ledger.append(
                "wave_members_dropped",
                {
                    "cycle": self.cycles,
                    "selected": list(selected),
                    "dropped": list(dropped),
                },
            )
        if wave:
            # The goal's own spine has to answer "which calculations ran
            # together", because the barrier is the whole design and the
            # dispatch receipt is overwritten every cycle. The order is
            # the Agent's and is recorded unsorted.
            self.ledger.append(
                "wave_selected",
                {
                    "cycle": self.cycles,
                    "node_ids": list(wave),
                },
            )
        return wave

    def _execution_wave_decision(self) -> Any | None:
        """Return the planning session's explicit execution boundary.

        Every current live Agent session carries this field. A missing field
        is therefore undecided at the production planning boundary; only the
        separate human-direct review surface and an already-started resumed
        run can proceed without it.
        """

        decision = getattr(self.session, "execution_wave_decision", None)
        return decision

    def _park_for_execution_wave_decision(
        self, decision: Any, *, reason: str
    ) -> None:
        """Persist an unexecuted scientific boundary without inventing one."""

        payload = {
            "cycle": self.cycles,
            "approval_file": str(self.bundle_file),
            "review_file": str(self.review_file) if self.review_file else "",
            "execution_wave_decision": decision.public_record(),
            "reason": reason,
        }
        self.ledger.append(
            "execution_wave_decision_pending",
            payload,
            idempotency_key=(
                f"execution-wave-decision-pending:{self.goal_id}:{self.cycles}"
            ),
        )
        ready = tuple(getattr(decision, "ready_node_ids", ()) or ())
        self.result = GoalLoopResultV1(
            goal_id=self.goal_id,
            settlement="execution_wave_decision_pending",
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            reasons=(
                reason,
                self._nothing_launched_this_cycle()
                + (
                    "; the Agent may make an explicit execution decision"
                    if ready or decision.state != "undecided"
                    else "; no calculation of workflow "
                    f"{decision.workflow_id or '(unnamed)'} is ready to run, "
                    "so no wave can be selected on it"
                ),
            ),
        )
        self.phase = "parked"

    def _undecided_boundary_reason(self, decision: Any) -> str:
        """Say what the Agent decided, when its last plan left none.

        A plan resets the execution boundary to undecided for the
        workflow it plans, which is right -- a wave names one workflow --
        and the park then said "the Agent made no execution-boundary
        decision". losartan-micropka-r2 cycle 4 (CUHK, 2026-09-18)
        selected c-neutral-opt on losartan-micropka-rev4a, was told
        "this wave is what will be submitted", then planned
        losartan-micropka-r4-settlement, which has no calculation ready
        to run; the park denied the selection it had made. The session's
        own stream holds every selection; the reason names the last one
        a later plan replaced.
        """

        replaced = [
            item
            for item in _session_wave_selections(self.events_path)
            if item["status"] == "ready"
            and item["workflow_id"] != str(decision.workflow_id or "")
        ]
        planned = str(decision.workflow_id or "") or "(unnamed)"
        if not replaced:
            # "No decision" was the whole word even when the session had
            # been told, once, that a decision was pending and then ended
            # on text -- R10 Q20 G1's cycle 4 wrote its wave that way (CUHK
            # 2153658). The park says how the session ended; it never says
            # what the text meant, because no text is read as a decision.
            told = (
                "; the session was told once that the decision was pending "
                "and ended on text, calling neither select_execution_wave "
                "nor continue_execution_reasoning, and the host reads no "
                "decision from text"
                if _pending_decision_answered_in_text(self.events_path)
                else ""
            )
            return (
                "the Agent made no execution-boundary decision on workflow "
                + planned
                + told
            )
        last = replaced[-1]
        return (
            f"the Agent selected {', '.join(last['node_ids'])} on workflow "
            f"{last['workflow_id']}; its later plan of workflow {planned} "
            "replaced that workflow, and it selected no wave on "
            f"{planned}; the host holds the boundary of the last planned "
            f"workflow only, so the selection on {last['workflow_id']} is "
            "not submitted"
        )

    def _nothing_launched_this_cycle(self) -> str:
        """Say what did not launch without denying what already ran.

        A goal parks on the wave it has not decided, which may be its
        second: the ledger that records the pending decision then also
        records the runs before it. "No engine launch occurred" was true
        of the parked cycle and read as a statement about the goal -- a
        live xTB goal (CUHK 2142404) ran two optimisations in cycle 1 and
        settled with a sentence from which a reader concludes nothing was
        computed. The sentence is scoped to its cycle and names what the
        goal already holds.
        """

        sentence = (
            "no scheduler submission or engine launch occurred in cycle "
            f"{self.cycles}"
        )
        earlier = [
            entry["payload"]
            for entry in self.ledger.entries()
            if entry["kind"] == "run_recorded"
            and int(entry["payload"].get("engine_calls_consumed", 0) or 0) > 0
        ]
        if earlier:
            calls = sum(
                int(item.get("engine_calls_consumed", 0) or 0)
                for item in earlier
            )
            cycles = ", ".join(str(item.get("cycle")) for item in earlier)
            sentence += (
                f"; this goal already records {calls} engine call(s) in "
                f"cycle(s) {cycles}, and that evidence stands"
            )
        return sentence

    def _execution_decision_pending_resume(self) -> None:
        """Keep a resumed pending decision pending; never replay it serially."""

        decision = self._execution_wave_decision()
        assert decision is not None
        self.result = GoalLoopResultV1(
            goal_id=self.goal_id,
            settlement="execution_wave_decision_pending",
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            reasons=(
                "the durable execution boundary remains " + decision.state,
                self._nothing_launched_this_cycle(),
            ),
        )
        self.phase = "parked"

    def _unanswerable_terminal_states(self) -> dict[str, str]:
        """How each node ended, for the endings a human has to read.

        Two exclusions, and they are the same exclusion twice. A stage
        the plan declared non-executable was displayed with the review,
        never approved and never launched, so it has no ending to answer
        -- counting its `not_launched` returned a goal whose every
        executable node had validated (live, 2026-09-03).

        A calculation outside *this wave* is the same: the Agent chose to
        see this wave's evidence before deciding the next, which is what
        the barrier is for, so every approved node outside the cohort is
        `not_launched` by construction and by design. Without this the
        first barrier of every multi-wave goal settled
        `returned_to_human` -- and the only cohort that escaped was one
        holding the whole approved partition, where the barrier does
        nothing. A member of *this* wave that never launched is still
        read: that is a launch that should have happened.

        The cohort manifest is the authority, because it is digest-bound
        to this approval and was written before any element started. No
        manifest means no cohort, and then every unlaunched node counts
        exactly as it always did.
        """

        from chemsmart.agent.cohort import read_cohort_manifest

        retained = set(self._declared_non_executable_ids())
        deferred: set[str] = set()
        if self.run_directory is not None:
            manifest = read_cohort_manifest(self.run_directory)
            if manifest is not None:
                members = set(manifest.node_ids)
                deferred = {
                    str(node.node_id)
                    for node in (self.outcome.nodes if self.outcome else ())
                    if str(node.node_id) not in members
                }
        return {
            str(node.node_id): str(node.state)
            for node in (self.outcome.nodes if self.outcome else ())
            if str(node.state) != "validated"
            and not (
                str(node.state) == "not_launched"
                and (
                    str(node.node_id) in retained
                    or str(node.node_id) in deferred
                )
            )
        }

    def _refusals_the_results_now_answer(
        self, session_delivery: "_AnalysisDelivery | None"
    ) -> dict[str, str]:
        """Verified refusals the registered results, read now, contradict.

        The session verifies a refusal against the results registered when
        it writes it, and a planning session writes it before the run it
        plans. A refusal of a quantity no reader serves, made before a
        Gaussian run whose log then prints it, was true when written and
        would settle unreachable_from_evidence over the printed value. The
        same reading is made here, over what the workspace now holds.
        """

        if session_delivery is None:
            return {}
        named = {
            observable_id: session_delivery.unreachable_producers.get(
                observable_id, ("", "", "")
            )
            for observable_id in session_delivery.verified_unreachable_ids
        }
        # A refusal of presence -- a named selector or a blocked node; a
        # refused precision stands on the delivered number instead.
        named = {
            key: value for key, value in named.items() if value[0] or value[2]
        }
        if not named:
            return {}
        from chemsmart.agent.live_session import (
            discover_registered_result_artifacts,
        )
        from chemsmart.agent.tool_runtime import (
            refusal_read_against_results,
            selector_declared_by,
        )
        from chemsmart.analysis.result_readers import (
            registered_reader_programs,
        )

        # The envelope says which producers exist; every registered result
        # is read, whichever program wrote it.
        programs = tuple(
            str(program)
            for program, _engines in getattr(
                getattr(self, "envelope", None),
                "allowed_program_engines",
                (),
            )
            or ()
        ) or tuple(registered_reader_programs())
        try:
            artifacts = {
                artifact.artifact_id: artifact
                for artifact in discover_registered_result_artifacts(
                    self.workspace
                )
            }
        except Exception:  # noqa: BLE001 - nothing registered reads as none
            artifacts = {}
        reread: dict[str, str] = {}
        for observable_id, (selector, jobtype, _node) in sorted(named.items()):
            still, basis = refusal_read_against_results(
                artifacts=artifacts,
                observable_id=observable_id,
                selector=selector,
                jobtype=jobtype,
                programs=tuple(registered_reader_programs()),
                selector_declared=bool(
                    selector_declared_by(programs, selector, jobtype)
                ),
                is_verified=True,
                basis="",
            )
            if not still:
                reread[observable_id] = basis.lstrip("; ")
        return reread

    def _latest_completion_stream(self) -> tuple[Path, str] | None:
        """The newest stream of this goal that holds a completion receipt.

        Streams are ordered by cycle, a cycle's planning session before
        its run; the run being settled is skipped, because it is the one
        that holds none. Returns the stream and how a reader names it.
        """

        agent_dir = self.workspace / ".chemsmart-agent"
        candidates: list[tuple[int, int, Path, str]] = []
        for entry in self.ledger.entries():
            if entry.get("kind") != "session_stream_recorded":
                continue
            payload = entry.get("payload") or {}
            run_id = str(payload.get("run_id") or "")
            if not run_id:
                continue
            cycle = int(payload.get("cycle") or 0)
            candidates.append(
                (
                    cycle,
                    0,
                    agent_dir / "runs" / run_id / "events.jsonl",
                    f"cycle {cycle}'s planning session",
                )
            )
        for stream in self.goal_dir.glob("runs/cycle-*/events.jsonl"):
            try:
                cycle = int(stream.parent.name.split("-", 1)[1])
            except (IndexError, ValueError):
                continue
            candidates.append((cycle, 1, stream, f"cycle {cycle}'s run"))
        settling = (
            (self.run_directory / "events.jsonl").resolve()
            if self.run_directory is not None
            else None
        )
        for _cycle, _order, stream, label in sorted(
            candidates, key=lambda item: (item[0], item[1]), reverse=True
        ):
            if settling is not None and stream.resolve() == settling:
                continue
            if _stream_completion_status(stream):
                return stream, label
        return None

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
        # A refusal is verified against the evidence of the moment it is
        # written, and a planning session writes it before the run it
        # plans: "no registered result exists it could be read from" was
        # true then and is read again now, against what the run wrote.
        reread = self._refusals_the_results_now_answer(session_delivery)
        delivery_kwargs: dict[str, Any] = dict(
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
                    if observable_id not in reread
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
            # The plan's failed acceptance criteria, read at the goal's
            # grain as the planning path reads them: every verdict of the
            # goal's streams, answered by a decision in any of them. A
            # run's own stream never holds a decision, and this path read
            # only it -- with earlier cycles' rejections carried as a set
            # of results no later decision could answer -- so a verdict a
            # woken session had answered by citing the failed receipt, as
            # the wake prescribes, was signed "no budget remains to answer
            # it", and a number standing on it was held stale (R11 truth,
            # test_a_goal_word_reads_every_criterion_the_goal_holds).
            goal_streams=_goal_streams(
                self.ledger, self.workspace, self.goal_id
            ),
            failed_artifact_sha256s=failed_artifacts(self.workspace),
            uncharacterised_artifact_sha256s=uncharacterised_artifacts(
                self.workspace
            ),
            flagged_artifact_sha256s=_flagged_artifact_sha256s(
                _goal_anomalies(self.ledger)
            ),
            goal_delivered_ids=_goal_delivered_ids(
                self.workspace, self.goal_id
            ),
            # The executor's stream holds no decision, so a finding the
            # session recorded before this run reaches the word's
            # reasons only through the record.
            goal_findings=_goal_findings(self.workspace, self.goal_id),
            # The word carries the goal's latest score of every number it
            # delivers, not only the one this run's completion scored.
            expectation_streams=_goal_streams(
                self.ledger, self.workspace, self.goal_id
            ),
        )
        run_delivery = _analysis_delivery(
            self.run_directory / "events.jsonl", **delivery_kwargs
        )
        # What the executor refused to launch, in this run's own stream,
        # kept before a chainless run is re-read from another stream.
        launch_refusals = run_delivery.stopped_by
        # A run whose stream holds no completion receipt carried no
        # analysis chain: nothing it computed was read, so it delivered
        # nothing and certified nothing. The settlement read that empty
        # stream, found no limitation list, took "nothing listed" for
        # "nothing undelivered", and wrote "workflow completed with its
        # analysis chain; the host completion gate certified the
        # delivery" over declared observables no cycle had claimed: four
        # of six in r9/gaussian g1 and g3, two of three in the R9 merged
        # smoke goal; r10/q2 g1-hono's only completion was partial and
        # its two falsified expectations never reached the word. The
        # goal's delivery is the one its latest completion receipt holds.
        run_completion = run_delivery.completion_status
        # Read from the run's own records, never from the executor's word.
        # This used to need two witnesses -- the executor's empty analysis
        # status and no receipt -- and the executor's word for an approved
        # toolchain with no analysis node is "completed" (all() over
        # nothing), so the first witness never agreed: every archived
        # chainless run was such an empty chain, the six repaired above
        # among them, and each still signed achieved. R10 Q21's g2-hooh
        # (CUHK 2153668) was the seventh, over the words "no completion
        # gate certified this delivery". A run whose stream holds no
        # completion receipt and in which no analysis node ran -- none
        # planned, or every one declared non-executable -- certified
        # nothing, whatever its executor reported.
        chainless = not run_delivery.completion_status and not (
            _analysis_nodes_run(self.run_directory / "events.jsonl")
        )
        stands_on = ""
        if chainless:
            latest = self._latest_completion_stream()
            if latest is not None:
                latest_stream, stands_on = latest
                run_delivery = _analysis_delivery(
                    latest_stream, **delivery_kwargs
                )
        goal_undelivered = (
            _goal_open_declared_ids(self.ledger, self.workspace, self.goal_id)
            if chainless
            else ()
        )
        # Something was owed and the completion the delivery stands on did
        # not pass: a partial chain, or none at all over declared
        # observables. Asked of every run, as the planning path asks it
        # (`_delivery_certified`), so an achieved word is never signed over
        # an uncertified delivery by a path that forgot to look.
        uncertified = not _delivery_certified(run_delivery) and bool(
            run_delivery.completion_status
            or _required_declared_ids(self.ledger)
        )
        if run_delivery.claims_rendered:
            self.standing_stale = _held_quantity_ids(run_delivery)
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
        # session of this cycle, closes the ids it names -- unless the
        # run's results, read now, hold what it refused.
        refused = set(
            session_delivery.verified_unreachable_ids
            if session_delivery is not None
            else ()
        ) - set(reread)
        open_declared = tuple(
            observable_id
            for observable_id in dict.fromkeys(
                run_delivery.undelivered_declared_ids + goal_undelivered
            )
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
            or uncertified
        )
        achieved = _achieved(self.execute_result)
        # What a run without a chain leaves the settlement to say, in
        # front of whatever it says: no quantity of it was extracted or
        # claimed (its validators still read its outputs), and which
        # completion the goal's delivery stands on.
        chainless_prefix = (
            f"cycle {self.cycles}: the workflow ran without an analysis "
            "chain, so no quantity it computed was extracted or claimed; "
            + (
                f"the goal's delivery stands on {stands_on}, whose "
                f"completion is {run_delivery.completion_status}"
                if stands_on
                else "no completion receipt in any cycle of this goal "
                "certified a delivery"
            )
            if chainless
            else ""
        )
        # The uncertified word for a run that carried a chain names the
        # receipt its own chain minted.
        uncertified_reason = chainless_prefix or (
            f"cycle {self.cycles}: this run's completion receipt "
            f"{run_delivery.completion_receipt_sha256[:8]} is "
            f"{run_delivery.completion_status or 'absent'}, so no "
            "completion gate certified the delivery"
        )
        # An observable the session refused, and an observable it
        # delivered whose precision it refused, settle by one rule: the
        # host verified both as unreachable, and the charter calls that
        # ending a deliverable. Only the first was ever read, so a goal
        # that ran engines and had its stated accuracy certified
        # unreachable settled achieved.
        # An earlier cycle's claim this cycle's verified refusal superseded
        # (`_claims_a_later_refusal_supersedes`): the refusal answers it.
        superseded = (
            _claims_a_later_refusal_supersedes(
                session_delivery,
                set(run_delivery.claim_rows)
                | set(session_delivery.claim_rows),
                _goal_delivered_ids(self.workspace, self.goal_id),
                _first_declarations(self.ledger),
                verified=refused,
            )
            if session_delivery is not None
            else {}
        )
        refused_ids = tuple(
            dict.fromkeys(
                tuple(run_delivery.undelivered_declared_ids)
                + tuple(run_delivery.refused_requirement_ids)
                + tuple(superseded)
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
                f"cycle {self.cycles}: "
                + _VERIFIED_REFUSAL_LEAD
                + "; ".join(
                    f"{observable_id} -- "
                    f"{session_delivery.unreachable_bases.get(observable_id, '')}"
                    + (
                        _superseded_claim_named(superseded[observable_id])
                        if observable_id in superseded
                        else ""
                    )
                    for observable_id in refused_ids
                )
            )
            carried, evidence = _what_the_delivery_carries(
                run_delivery,
                _goal_anomalies(self.ledger),
                _settlement_evidence(session_delivery),
                superseded,
            )
            self.ledger.settle(
                "unreachable_from_evidence",
                reasons=(reason, *carried),
                evidence=evidence,
            )
            self._settled("unreachable_from_evidence", (reason, *carried))
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
                    "verdicts": list(_holding_verdicts(run_delivery)),
                    "stale_quantity_ids": list(unrefreshed),
                    "unclaimed_output_ids": list(
                        run_delivery.unclaimed_output_ids
                    ),
                    "undelivered_declared_observable_ids": list(open_declared),
                    "unresolved_requirement_ids": list(open_requirements),
                    "engine_calls_remaining": budgets.engine_calls_remaining,
                    **(
                        {"uncertified": uncertified_reason}
                        if uncertified
                        and (
                            chainless
                            or not (
                                run_delivery.unanswered_verdicts
                                or unrefreshed
                                or run_delivery.unclaimed_output_ids
                                or open_declared
                                or open_requirements
                            )
                        )
                        else {}
                    ),
                    **({"refusals_reread": dict(reread)} if reread else {}),
                },
            )
            self.phase = "plan"
            return
        reread_reasons = (
            (
                f"cycle {self.cycles}: refusals verified before this run's "
                "results existed were read again against them and no longer "
                "hold: "
                + "; ".join(
                    f"{observable_id} -- {basis}"
                    for observable_id, basis in reread.items()
                ),
            )
            if reread
            else ()
        )
        if achieved and open_delivery:
            # Same failure, nothing left to answer it with.
            if uncertified and not (
                run_delivery.unanswered_verdicts
                or unrefreshed
                or run_delivery.unclaimed_output_ids
                or open_requirements
            ):
                reason = (
                    uncertified_reason
                    + ", and no revision remains to certify it"
                    + (
                        "; "
                        + _undelivered_declared_named(
                            run_delivery, open_declared
                        )
                        if open_declared
                        else ""
                    )
                )
                self.ledger.settle(
                    "returned_to_human", reasons=(reason, *reread_reasons)
                )
                self._settled(
                    "returned_to_human", open_declared or (uncertified_reason,)
                )
                return
            if run_delivery.unanswered_verdicts:
                reason = (
                    f"cycle {self.cycles}: a validation verdict failed and "
                    "no budget remains to answer it: "
                    + _unanswered_verdicts_named(run_delivery)
                )
                open_items = run_delivery.unanswered_verdicts
            elif run_delivery.inherited_unanswered:
                # The planning path's own sentence: the verdict, its
                # number and receipt, and the numbers standing on it.
                reason = (
                    f"cycle {self.cycles}: "
                    + _inherited_verdict_reason(run_delivery)
                    + "; no budget remains to answer it"
                )
                open_items = tuple(
                    dict.fromkeys(
                        quantity_id
                        for _verdict, standing in (
                            run_delivery.inherited_unanswered
                        )
                        for quantity_id in standing
                    )
                )
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
                    f"cycle {self.cycles}: "
                    + _undelivered_declared_named(
                        run_delivery, run_delivery.undelivered_declared_ids
                    )
                    + "; no revision remains to deliver them"
                    + (
                        "; " + " | ".join(run_delivery.open_declared_misses)
                        if run_delivery.open_declared_misses
                        else ""
                    )
                )
                open_items = run_delivery.undelivered_declared_ids
            self.ledger.settle(
                "returned_to_human",
                reasons=(
                    *((chainless_prefix,) if chainless else ()),
                    *reread_reasons,
                    reason,
                ),
            )
            self._settled("returned_to_human", open_items)
            return
        if achieved:
            goal_anomalies = _goal_anomalies(self.ledger)
            word, why = _achieved_word(
                run_delivery, goal_anomalies, superseded=superseded
            )
            evidence = _settlement_evidence(run_delivery)
            if word == "achieved_with_observations":
                evidence = _anomaly_evidence(
                    evidence, goal_anomalies, run_delivery, superseded
                )
            reasons = (
                (
                    chainless_prefix + "; " + why[0]
                    if chainless
                    else f"cycle {self.cycles}: workflow completed"
                    + (" with its analysis chain" if run_completion else "")
                    + "; "
                    + why[0]
                ),
                # Every reason the word carries, not only the first:
                # the second one names which delivered number stands
                # on a flagged result.
                *why[1:],
            )
            if self._open_reading(
                path="run",
                run=f"goals/{self.goal_id}/runs/cycle-{self.cycles}",
                state=word,
                reasons=reasons,
                result_reasons=why,
                evidence=evidence,
            ):
                return
            self.ledger.settle(word, reasons=reasons, evidence=evidence)
            self._record_qualification()
            self._settled(word, why)
            return
        if (
            budgets.engine_calls_remaining <= 0
            and budgets.revisions_remaining > 0
            and budgets.wall_seconds_remaining > 0
        ):
            # The engine line is spent, but a provider-free cycle can
            # still read the terminal outcome, extract the evidence that
            # exists, and decide what a scientific failure means.  This
            # used to require an analysis-completion receipt with an
            # undelivered declared id; a failed Hessian has native bytes
            # and a host-recorded anomaly before any such receipt exists,
            # so it settled exhausted before the Agent could interpret it
            # (CUHK acetamide r9, 2026-09-18).  The next plan is explicitly
            # analysis-only: the ordinary plan-time budget gate still
            # refuses every engine node, and the host neither selects a
            # recovery calculation nor decides what the finding means.
            terminal_states = self._unanswerable_terminal_states()
            self.ledger.append(
                "recovery_opened",
                {
                    "cycle": self.cycles,
                    "terminal_states": dict(sorted(terminal_states.items())),
                    "undelivered_declared_observable_ids": list(
                        run_delivery.undelivered_declared_ids
                    ),
                    "unclaimed_output_ids": list(
                        run_delivery.unclaimed_output_ids
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
        terminal_states = self._unanswerable_terminal_states()
        repairable = {
            node_id: state
            for node_id, state in terminal_states.items()
            if state in REPAIRABLE_TERMINAL_STATES
        }
        if terminal_states and not repairable:
            # A launch the executor refused is an event in the run's own
            # stream, and the settlement quotes it: R10 Q15 g1 returned
            # naming "pbnz-opt=not_launched, ts-search=not_launched" while
            # the stream held why -- a stale input check of another
            # program's bytes, false, which only the quote would have let a
            # reader see.
            reason = (
                f"cycle {self.cycles}: the run ended in a state no revision "
                "can answer: "
                + ", ".join(
                    f"{node_id}={state}"
                    for node_id, state in sorted(terminal_states.items())
                )
                + (
                    "; " + "; ".join(launch_refusals)
                    if launch_refusals
                    else ""
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
                # What the run left open, read here exactly as the branch
                # above reads it. Since the executor reads criteria (R10
                # Q19) a claim standing on a failed criterion makes the
                # chain partial, so every such run arrives here -- and
                # this row wrote "verdicts": [] over L-S2's failed
                # external-stability criterion (CUHK Slurm 2153514).
                "verdicts": list(_holding_verdicts(run_delivery)),
                "stale_quantity_ids": list(unrefreshed),
                "unclaimed_output_ids": list(
                    run_delivery.unclaimed_output_ids
                ),
                "undelivered_declared_observable_ids": list(open_declared),
                "unresolved_requirement_ids": list(open_requirements),
                "engine_calls_remaining": budgets.engine_calls_remaining,
                # A partial chain's row already says analysis_status
                # partial; only a chainless run needs the sentence.
                **(
                    {"uncertified": chainless_prefix}
                    if uncertified and chainless
                    else {}
                ),
                **({"refusals_reread": dict(reread)} if reread else {}),
            },
        )
        self.phase = "plan"

    # -- the reading turn ---------------------------------------------------

    def _reading_opened_for_cycle(self) -> Mapping[str, Any] | None:
        """The held settlement this cycle's reading turn reads under."""

        for entry in reversed(self.ledger.entries()):
            if entry["kind"] != "reading_opened":
                continue
            payload = entry["payload"]
            if int(payload.get("cycle") or 0) == self.cycles:
                return payload
        return None

    def _open_reading(
        self,
        *,
        path: str,
        run: str,
        state: str,
        reasons: tuple[str, ...],
        result_reasons: tuple[str, ...],
        evidence: Mapping[str, Any],
    ) -> bool:
        """Hold a certified delivery's settlement for one reading turn.

        A complete delivery used to settle with no session reading what
        the run computed: in 10 of 23 archived successful engine goals
        nothing ever looked at the results of the last run, so a
        phenomenon present only in computed results reached nobody. The
        word is computed first and recorded here, exactly as it would
        have been written, so the reading can add to the settlement and
        can never change it -- and the ledger keeps what the goal would
        have settled with had nobody read.
        """

        if not self.reading_turn or self.goal is None:
            return False
        if state not in _CERTIFIED_WORDS:
            return False
        if self._reading_opened_for_cycle() is not None:
            return False
        opened = self.ledger.append(
            "reading_opened",
            {
                "cycle": self.cycles,
                "path": path,
                "run": run,
                "state": state,
                "reasons": tuple(reasons),
                "result_reasons": tuple(result_reasons),
                "evidence": dict(evidence),
            },
            idempotency_key=f"reading-opened:{self.goal_id}:{self.cycles}",
        )
        if not opened:
            return False
        self.phase = "read"
        return True

    def _read(self) -> None:
        """One session reads the certified delivery; then the goal settles.

        It launches nothing and admits nothing: its context grants zero
        engine calls, so a calculation plan is refused where it is made,
        and whatever it plans is never decided. What it records -- claims
        and findings -- reaches the settlement beside the word the host
        had already computed. A reading that fails, or that a previous
        process began and never recorded, costs the delivery nothing.
        """

        opened = self._reading_opened_for_cycle()
        if opened is None:  # pragma: no cover - the phase is set with it
            raise ContractError("a reading turn opened no held settlement")
        if self._reading_interrupted:
            summary = {
                "cycle": self.cycles,
                "run_id": "",
                "terminal_state": "",
                "error": (
                    "the reading turn a previous process opened never "
                    "recorded; the goal settles on the held word"
                ),
                "findings": (),
            }
            self.ledger.append(
                "reading_recorded",
                summary,
                idempotency_key=(
                    f"reading-recorded:{self.goal_id}:{self.cycles}"
                ),
            )
            self._settle_after_reading(opened, summary)
            return
        context = _reading_context(
            self.goal,
            self.ledger,
            self.outcome,
            workspace=self.workspace,
            opened=opened,
        )
        started = time.monotonic()
        session: Any = None
        error = ""
        try:
            session = self.plan_session(
                task=self.task,
                provider=self.provider,
                provider_config_file=self.provider_config_file,
                workspace=self.workspace,
                execution_enabled=False,
                approval_file=None,
                execution_envelope_file=self.execution_envelope_file,
                analysis_completion_file=self.analysis_completion_file,
                review_file=None,
                goal_context=context,
                **self.session_kwargs,
            )
        except ContractError as exc:
            error = f"{type(exc).__name__}: {exc}"
        wall_seconds = time.monotonic() - started
        run_id = _session_run_id(session) if session is not None else ""
        events_path: Path | None = None
        if run_id:
            candidate = (
                self.workspace
                / ".chemsmart-agent"
                / "runs"
                / run_id
                / "events.jsonl"
            )
            # Only the stream the session named: the newest-first glob
            # would hand this goal another session's reading.
            events_path = candidate if candidate.is_file() else None
        findings: tuple[Mapping[str, Any], ...] = ()
        uncertainties: tuple[str, ...] = ()
        claims = decisions = typed_reads = 0
        if events_path is not None:
            self._record_workspace(events_path, "")
            delivery = _analysis_delivery(events_path)
            findings = delivery.findings
            uncertainties = tuple(delivery.decision_uncertainties)
            claims = delivery.claims
            decisions = delivery.decisions
            typed_reads = _session_typed_reads(events_path)
        summary = {
            "cycle": self.cycles,
            "run_id": run_id,
            "terminal_state": str(
                getattr(session, "terminal_state", "") or ""
            ),
            "error": error,
            "claims": claims,
            "decisions": decisions,
            "typed_reads": typed_reads,
            "findings": tuple(dict(row) for row in findings),
            "decision_uncertainties": uncertainties,
            "cost": {
                **_session_provider_cost(events_path),
                "driver_wall_seconds": round(wall_seconds, 3),
            },
        }
        self.ledger.append(
            "reading_recorded",
            summary,
            idempotency_key=f"reading-recorded:{self.goal_id}:{self.cycles}",
        )
        self._settle_after_reading(opened, summary)

    def _settle_after_reading(
        self, opened: Mapping[str, Any], summary: Mapping[str, Any]
    ) -> None:
        """Settle on the held word, with the reading beside it.

        The word is the one recorded before the reading began; the
        reading adds reasons and evidence and nothing else, so a goal's
        settlement with the policy on differs from the one it would have
        had with it off only by what a session found in its results.
        """

        state = str(opened.get("state") or "")
        findings = tuple(
            row
            for row in summary.get("findings") or ()
            if isinstance(row, Mapping)
        )
        run_id = str(summary.get("run_id") or "")
        if summary.get("error") or not run_id:
            line = "the reading turn recorded nothing: " + str(
                summary.get("error") or "the session named no stream"
            )
        else:
            # What the reading did, from its own stream's counts. The line
            # used to quote the session's terminal word and assert that it
            # "read the delivered results": a reading that recorded only a
            # decision ends 'blocked' -- the planning word for stopping
            # before a workflow, which is how a reading stops -- and two
            # sealed Q6 goals settled saying "ended blocked ... read the
            # delivered results" of a session that had read nothing. The
            # terminal word stays in the evidence, and is quoted here only
            # when it says the session failed.
            reads = int(summary.get("typed_reads") or 0)
            claims = int(summary.get("claims") or 0)
            decisions = int(summary.get("decisions") or 0)
            terminal = str(summary.get("terminal_state") or "")
            line = (
                f"the reading turn ({run_id}) made "
                + (
                    f"{reads} typed read(s) of the delivered results"
                    if reads
                    else "no typed read of the delivered results"
                )
                + " and recorded "
                + (
                    f"{len(findings)} finding(s), in its own words beneath"
                    if findings
                    else "no finding"
                )
                + (
                    f"; it recorded {claims} claim(s) and "
                    f"{decisions} decision(s)"
                    if claims or decisions
                    else ""
                )
                + (
                    f"; its session ended {terminal}"
                    if terminal in {"failed", "cancelled"}
                    else ""
                )
            )
        # The rule sends a check that the delivery holds, and "nothing
        # else bears on it", into the decision's words; the planning
        # session's recorded uncertainties reach the settlement, and the
        # reading's were never read (R10 Q29 census: all 14 archived
        # readings recorded some; 48 of their 60 appear nowhere in the
        # settlement, the rest only where a planning session had stated
        # the same words).
        uncertainties = tuple(
            str(item) for item in summary.get("decision_uncertainties") or ()
        )
        added = (
            (line,)
            + (
                (
                    "the reading's recorded decision states its "
                    "uncertainties: " + " | ".join(uncertainties),
                )
                if uncertainties
                else ()
            )
            + _finding_reasons(findings)
        )
        evidence = dict(opened.get("evidence") or {})
        evidence["reading"] = {
            key: summary.get(key)
            for key in (
                "run_id",
                "terminal_state",
                "error",
                "claims",
                "decisions",
                "decision_uncertainties",
                "typed_reads",
                "cost",
            )
        }
        if findings:
            known = {
                str(row.get("receipt_sha256") or "")
                for row in evidence.get("findings") or ()
                if isinstance(row, Mapping)
            }
            evidence["findings"] = tuple(
                evidence.get("findings") or ()
            ) + tuple(
                dict(row)
                for row in findings
                if str(row.get("receipt_sha256") or "") not in known
            )
            evidence["receipt_sha256s"] = tuple(
                sorted(
                    set(evidence.get("receipt_sha256s") or ())
                    | {
                        str(row.get("receipt_sha256") or "")
                        for row in findings
                        if row.get("receipt_sha256")
                    }
                )
            )
        reasons = tuple(opened.get("reasons") or ()) + added
        if str(opened.get("path") or "") == "run":
            self.ledger.settle(state, reasons=reasons, evidence=evidence)
            self._record_qualification()
            self._settled(
                state, tuple(opened.get("result_reasons") or ()) + added
            )
            return
        self.result = _write_delivery_settlement(
            self.ledger,
            goal_id=self.goal_id,
            cycles=self.cycles,
            revisions_admitted=self.revisions_admitted,
            settled=state,
            reasons=reasons,
            evidence=evidence,
            workspace=self.workspace,
        )
        self.phase = "settled"


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
