"""The loop's bookkeeping, driven with stubbed session/resolve/execute
so no provider and no engine runs: cycle accounting, the wake context,
admission wiring, and every settlement path. The live qualification
drives the real chain; these tests pin that the connective tissue does
exactly what the charter text says and nothing more.
"""

import json
import shutil
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.driver import GoalDriver, run_goal_loop
from chemsmart.agent.execution import build_program_execution_receipt
from chemsmart.agent.goal import GoalLedger
from chemsmart.agent.runtime.event_store import RuntimeEventStore

from .test_runtime_v2_launch_fence import _reserve


def _envelope_file(tmp_path, calls=6, excursions=0):
    target = tmp_path / "execution-envelope.yaml"
    target.write_text(
        "\n".join(
            (
                *(
                    (f"max_excursion_calls: {excursions}",)
                    if excursions
                    else ()
                ),
                "schema_version: chemsmart.bounded-execution-envelope.v1",
                "mode: bounded-local",
                "allowed_program_engines:",
                "  orca:",
                "  - cpu",
                "  pyscf:",
                "  - cpu",
                "resources:",
                "  execution_target: run",
                "  cores: 4",
                "  memory_gb: 16",
                "  gpu_count: 0",
                "  scratch_policy: server",
                "  node_timeout_seconds: 600",
                "episode_wall_time_seconds: 7200",
                "postprocess_reserve_seconds: 600",
                f"max_engine_calls: {calls}",
                f"scratch_root: {tmp_path / 'scratch'}",
            )
        )
        + "\n"
    )
    return target


def _review_payload(identity="b" * 64):
    return {
        "review_sha256": "d" * 64,
        "scientific_plan": {"scientific_identity_sha256": identity},
        "execution_envelope": {
            "allowed_program_engines": (("orca", ("cpu",)),),
        },
        "node_reviews": ({"project_settings_text": "orca:\n  gas: {}\n"},),
        "scientific_toolchain_plan": {"analysis_nodes": ()},
    }


def _write_session_stream(workspace, name, rows):
    run_dir = workspace / ".chemsmart-agent" / "runs" / name
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "events.jsonl").write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n",
        encoding="utf-8",
    )


def _engine_stream(
    tmp_path, target, *, failed, findings=("execution.process.timeout",)
):
    build_dir = tmp_path / f"build-{target.name}"
    store = RuntimeEventStore(
        build_dir / "events.jsonl", session_id="water-session"
    )
    _, plan, _m, _a, invocation = _reserve(store, build_dir)
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="failed" if failed else "engine_complete",
        exit_status=1 if failed else 0,
        child_exit_status=1 if failed else 0,
        engine_complete=not failed,
        validated=False,
        findings=tuple(findings) if failed else (),
        started_at="2026-08-04T00:00:00+00:00",
        finished_at="2026-08-04T00:00:05+00:00",
    )
    store.record_program_execution_receipt(
        turn_id="turn-1",
        workflow_id=plan.workflow_id,
        run_id="run.water-approval",
        receipt=receipt,
    )
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy(build_dir / "events.jsonl", target / "events.jsonl")


def _loop(
    tmp_path,
    *,
    sessions,
    executes,
    calls=6,
    max_revisions=5,
    excursions=0,
    goal_id="goal-t1",
    resolve=None,
):
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    session_iter = iter(sessions)
    execute_iter = iter(executes)

    def plan_session(**kwargs):
        step = next(session_iter)
        return step(workspace, kwargs)

    def resolve_review(**kwargs):
        if resolve is not None:
            return resolve(**kwargs)
        return ("d" * 64, tmp_path / "bundle.json")

    def execute_bundle(*, approval_file, workspace, run_directory):
        step = next(execute_iter)
        return step(run_directory)

    return run_goal_loop(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(
            tmp_path, calls, excursions=excursions
        ),
        goal_id=goal_id,
        granted_by="claude-owner-delegated-reviewer",
        max_revisions=max_revisions,
        plan_session=plan_session,
        resolve_review=resolve_review,
        execute_bundle=execute_bundle,
    )


def _planning_session(
    name, *, review=None, terminal="waiting_for_approval", wake_rows=()
):
    def step(workspace, kwargs):
        rows = list(wake_rows) or [{"kind": "session_started", "payload": {}}]
        _write_session_stream(workspace, name, rows)
        if review is not None:
            review_file = kwargs["review_file"]
            review_file.parent.mkdir(parents=True, exist_ok=True)
            review_file.write_text(json.dumps(review), encoding="utf-8")
        return SimpleNamespace(
            terminal_state=terminal, task_spec_sha256="a" * 64
        )

    return step


def _execute(tmp_path, *, failed, status, analysis="completed"):
    def step(run_directory):
        _engine_stream(tmp_path, run_directory, failed=failed)
        return SimpleNamespace(status=status, analysis_status=analysis)

    return step


_READ_OUTCOME_ROWS = (
    {
        "kind": "tool_started",
        "payload": {"request_id": "r1", "tool": "inspect_run_outcome"},
    },
    {"kind": "tool_succeeded", "payload": {"request_id": "r1"}},
    {
        "kind": "run_outcome_inspected",
        "payload": {
            "run": "goals/goal-t1/runs/cycle-1",
            "stream_sha256": "a" * 64,
        },
    },
)


def test_every_cycle_sees_the_goal_terms(tmp_path):
    """Cycle 1 used to receive nothing -- an analysis-only goal is
    single-cycle by construction, so its session never learned the
    budgets, the authority, or that a typed refusal is a deliverable.
    The terms are in hand when the loop starts; every cycle gets them,
    and cycles with a previous run get its typed outcome embedded."""

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    _loop(
        tmp_path,
        sessions=[
            capture(_planning_session("live-1", review=_review_payload())),
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=1,
    )
    first, second = contexts
    assert first["schema_version"] == "chemsmart.goal-wake-context.v1"
    assert first["budgets"] == {
        "binding_line": (
            "nearest exhausted: engine calls 6 of 6 remaining (100%)"
        ),
        "engine_calls_remaining": 6,
        "wall_seconds_remaining": 7200.0,
        "revisions_remaining": 1,
        "excursion_calls_remaining": 0,
    }
    assert first["previous_run"] == ""
    assert first["trajectory"] == ()
    assert "unreachable_observable_ids" in first["authority"]
    assert "attempt to refute" in first["authority"]
    assert first["deliverables"] == {
        "delivered_quantity_ids": (),
        "limitation_output_ids": (),
        "doubted_quantity_ids": (),
        "unanswered_failed_verdicts": (),
        "stale_quantity_ids": (),
        "unclaimed_output_ids": (),
        "undelivered_declared_observable_ids": (),
        # The ninth gap (2026-09-09): a number delivered under its id
        # can still miss the precision the task asked for, and the wake
        # names that too, so the next action follows the gap.
        "unresolved_requirement_ids": (),
        "sufficiency": (),
        "flagged_quantity_ids": (),
        "failed_source_quantity_ids": (),
        "characterised_source_quantity_ids": (),
        "uncharacterised_source_quantity_ids": (),
        "expectation_rows": (),
    }
    assert second["previous_run"] == "goals/goal-t1/runs/cycle-1"
    assert second["previous_run_outcome"]
    assert "unreachable_observable_ids" in second["authority"]
    assert "attempt to refute" in second["authority"]
    # The wake states what the previous run's own stream delivered --
    # here an engine failure with no claims, so every list is empty but
    # the record is present for the session to read.
    assert set(second["deliverables"]) == {
        "delivered_quantity_ids",
        "limitation_output_ids",
        "doubted_quantity_ids",
        # A host-rendered verdict saying a delivered structure is not
        # what the task required now reaches the wake beside what was
        # delivered: it is exactly the gap the next action should follow.
        "unanswered_failed_verdicts",
        # And beside it, the numbers that verdict took down with it --
        # rendered from the rejected result, sound arithmetic, no
        # structure left underneath.
        "stale_quantity_ids",
        # And what the host computed and no claim ever showed a reader.
        "unclaimed_output_ids",
        "undelivered_declared_observable_ids",
        # And which requirement a delivered number has not answered
        # yet, with the arithmetic behind that word.
        "unresolved_requirement_ids",
        "sufficiency",
        # And which delivered number stands on a result the host flagged.
        "flagged_quantity_ids",
        # And which delivered number stands on a run that did not meet
        # the promise it was launched under, checked or not.
        "failed_source_quantity_ids",
        "characterised_source_quantity_ids",
        "uncharacterised_source_quantity_ids",
        "expectation_rows",
    }


def test_cycle_one_approves_runs_and_settles_achieved(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", review=_review_payload()),
        ],
        executes=[
            _execute(tmp_path, failed=False, status="completed"),
        ],
    )
    assert result.settlement == "achieved"
    assert result.cycles == 1
    assert result.revisions_admitted == 0
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    kinds = [entry["kind"] for entry in ledger.entries()]
    assert kinds == [
        "goal_created",
        "run_started",
        "run_recorded",
        "goal_settled",
    ]


def test_a_failed_run_wakes_a_revision_that_recovers(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", review=_review_payload()),
            _planning_session(
                "live-2",
                review=_review_payload(),
                wake_rows=_READ_OUTCOME_ROWS,
            ),
        ],
        executes=[
            _execute(
                tmp_path,
                failed=True,
                status="partial",
                analysis="partial",
            ),
            _execute(tmp_path, failed=False, status="completed"),
        ],
    )
    assert result.settlement == "achieved"
    assert result.cycles == 2
    assert result.revisions_admitted == 1
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    admitted = [
        entry
        for entry in ledger.entries()
        if entry["kind"] == "revision_admitted"
    ]
    assert len(admitted) == 1
    payload = admitted[0]["payload"]
    assert payload["actor"] == "goal-approval:goal-t1"
    assert payload["granted_by"] == "claude-owner-delegated-reviewer"
    assert payload[
        "cited_evidence_event_hashes"
    ], "an admitted revision cites the terminal evidence it answered"
    assert payload["checks"]["evidence_read"] is True


def test_a_wake_handed_outcome_admits_without_ritual(tmp_path):
    """The host composed the wake context, embedded the previous run's
    typed outcome, and recorded that act in the ledger. A revision of
    that run is evidence-read by construction -- the first live goal
    round's gate demanded a re-read of what the host itself handed
    over, and blocked a scientifically sound revision on a wiring
    defect. The attestation, not the ritual, is the evidence."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", review=_review_payload()),
            _planning_session("live-2", review=_review_payload()),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=1,
    )
    assert result.settlement == "exhausted"
    assert result.revisions_admitted == 1
    entries = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    ).entries()
    by_kind = {}
    for entry in entries:
        by_kind.setdefault(entry["kind"], []).append(entry["payload"])
    (wake,) = by_kind["wake_composed"]
    assert wake == {"cycle": 2, "run": "goals/goal-t1/runs/cycle-1"}
    (admitted,) = by_kind["revision_admitted"]
    assert admitted["checks"]["evidence_read"] is True


def test_an_identity_change_returns_to_the_human(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", review=_review_payload()),
            _planning_session(
                "live-2",
                review=_review_payload(identity="9" * 64),
                wake_rows=_READ_OUTCOME_ROWS,
            ),
        ],
        executes=[
            _execute(
                tmp_path,
                failed=True,
                status="partial",
                analysis="partial",
            ),
        ],
    )
    assert result.settlement == "returned_to_human"
    assert any("identities" in r for r in result.reasons)


def _delivery_rows(
    *,
    completion="passed",
    limitations=(),
    claims=True,
    decision=True,
    failed_rule=False,
    claim_source="",
    doubt_ref="",
):
    """Stream shapes drawn from the live goal round's three sessions."""

    rows = [
        {
            "kind": "result_quantities_extracted",
            "payload": {"receipt_sha256": "e" * 64},
        },
    ]
    if failed_rule:
        rows.append(
            {
                "kind": "scientific_validation_evaluated",
                "payload": {
                    "receipt_sha256": "f" * 64,
                    "record": {
                        "rule_results": ({"rule_id": "same", "passed": False},)
                    },
                },
            }
        )
    if claims:
        claims_payload = {"receipt_sha256": "a1" + "a" * 62}
        if claim_source:
            claims_payload["record"] = {
                "claims": (
                    {
                        "source_receipt_sha256": claim_source,
                        "quantity_id": "dg_solv",
                    },
                )
            }
        rows.append(
            {
                "kind": "analysis_claims_recorded",
                "payload": claims_payload,
            }
        )
    if decision:
        decision_payload = {}
        if doubt_ref:
            decision_payload["record"] = {
                "evidence_refs": ("doubt:" + doubt_ref,)
            }
        rows.append(
            {
                "kind": "scientific_decision_recorded",
                "payload": decision_payload,
            }
        )
    if completion:
        rows.append(
            {
                "kind": "analysis_completion_evaluated",
                "payload": {
                    "receipt_sha256": "c1" + "c" * 62,
                    "status": completion,
                    "limitation_output_ids": list(limitations),
                },
            }
        )
    return tuple(rows)


def test_a_stated_limitation_settles_as_the_typed_refusal(tmp_path):
    """The live r3 shape: completion passed with a blocked required
    output, substitute claims delivered, decision recorded, every
    validation rule green. The old classifier watched failed rules --
    the one place an honest refusal leaves no trace -- and settled
    this achieved."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="complete",
                wake_rows=_delivery_rows(limitations=("dg",)),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "unreachable_from_evidence"
    assert any("dg" in reason for reason in result.reasons)
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    settled = ledger.entries()[-1]
    assert settled["payload"]["evidence"]["receipt_sha256s"]


def test_a_certified_delivery_with_observations_carries_the_word(tmp_path):
    """The one word a human reads first must not hide what the run
    found: a certified delivery whose completion receipt names host
    anomaly observations settles achieved_with_observations, on
    receipts, at the delivery path."""

    rows = _delivery_rows()
    for row in rows:
        if row["kind"] == "analysis_completion_evaluated":
            row["payload"]["anomaly_output_ids"] = (
                "anomaly:stationary_point.unexpected_order:unreplicated:"
                "0a1b2c3d",
            )
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="complete", wake_rows=rows),
        ],
        executes=[],
    )
    assert result.settlement == "achieved_with_observations"
    assert "stationary_point.unexpected_order" in result.reasons[0]
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    settled = ledger.entries()[-1]
    assert settled["payload"]["state"] == "achieved_with_observations"
    assert settled["payload"]["evidence"]["receipt_sha256s"]


def test_a_certified_clean_delivery_settles_achieved(tmp_path):
    """The live r1 shape: completion passed with no limitations."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1", terminal="complete", wake_rows=_delivery_rows()
            ),
        ],
        executes=[],
    )
    assert result.settlement == "achieved"


def test_an_uncertified_delivery_returns_naming_the_gate(tmp_path):
    """The live r2 shape: claims and a decision recorded, no
    completion event, terminal 'planned'. The word was right before,
    for a reason that would also fire on a clean delivery; the reason
    now names what is actually missing."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="planned",
                wake_rows=_delivery_rows(completion=""),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "returned_to_human"
    assert any("completion gate" in reason for reason in result.reasons)


def test_a_failed_rule_is_not_a_refusal(tmp_path):
    """The hazard the seal named and the round left unexercised: a
    delivering session whose stream holds one honestly failed
    validation rule must not settle as a refusal."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="complete",
                wake_rows=_delivery_rows(failed_rule=True),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "achieved"


def test_a_goal_is_not_a_resumable_queue(tmp_path):
    _loop(
        tmp_path,
        sessions=[_planning_session("live-1", review=_review_payload())],
        executes=[_execute(tmp_path, failed=False, status="completed")],
    )
    with pytest.raises(ContractError, match="one human decision"):
        _loop(
            tmp_path,
            sessions=[_planning_session("live-9", review=_review_payload())],
            executes=[],
        )


def test_the_stop_file_cancels_at_the_cycle_boundary(tmp_path):
    stop = tmp_path / "stop"
    stop.write_text("cancel\n")
    result = run_goal_loop(
        task="the goal task",
        workspace=(tmp_path / "ws2"),
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-t2",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kwargs: (_ for _ in ()).throw(
            AssertionError("no session may start after cancel")
        ),
        resolve_review=lambda **kwargs: ("d" * 64, tmp_path / "b.json"),
        execute_bundle=lambda **kwargs: None,
        stop_file=stop,
    )
    assert result.settlement == "returned_to_human"
    assert result.reasons == ("cancelled",)


def test_the_executor_accepts_a_stop_file_and_marks_cancelled():
    """The walk marks an unlaunched node cancelled with its typed id.

    The full executor path needs an approved bundle and engines; the
    live goal qualification exercises it. Here the contract is pinned
    at the seam: the entry accepts stop_file, and the cancelled node
    record carries the typed rule id beside the human sentence.
    """

    import inspect

    from chemsmart.agent.executor import (
        ExecutedNodeV1,
        execute_approved_workflow,
    )

    assert (
        "stop_file" in inspect.signature(execute_approved_workflow).parameters
    )
    node = ExecutedNodeV1(
        node_id="n",
        program="orca",
        jobtype="opt",
        state="cancelled",
        invocation_identity_sha256="",
        execution_receipt_sha256="",
        rule_ids=("execution.cancelled.human",),
        failure="cancelled by the human at a node boundary",
    )
    assert node.state == "cancelled"


def test_a_claim_from_a_doubted_receipt_returns_to_the_human(tmp_path):
    """The r3 residue: a session wrote the correct doubt and claimed
    the doubted number anyway. Typed as doubt:{receipt}, the settlement
    returns the number to the human even when the completion certified
    before the decision landed -- the intersection is computed from the
    stream, not from the completion word."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="complete",
                wake_rows=_delivery_rows(
                    claim_source="e" * 64, doubt_ref="e" * 64
                ),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "returned_to_human"
    assert any("dg_solv" in reason for reason in result.reasons)


def test_a_doubt_about_an_unclaimed_receipt_changes_nothing(tmp_path):
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="complete",
                wake_rows=_delivery_rows(
                    claim_source="e" * 64, doubt_ref="9" * 64
                ),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "achieved"


def test_the_wake_deliverables_come_from_the_previous_runs_stream(tmp_path):
    """Names quantities and stated limitations, never values: the wake
    session sees what already stands delivered, what the chain declared
    it could not produce, and what its own decisions doubt."""

    from chemsmart.agent.driver import (
        _analysis_delivery,
        _deliverables_record,
    )

    workspace = tmp_path / "ws"
    _write_session_stream(
        workspace,
        "cycle-1",
        _delivery_rows(
            limitations=("dg",),
            claim_source="e" * 64,
            doubt_ref="e" * 64,
        ),
    )
    record = _deliverables_record(
        _analysis_delivery(
            workspace
            / ".chemsmart-agent"
            / "runs"
            / "cycle-1"
            / "events.jsonl"
        )
    )

    assert record == {
        "delivered_quantity_ids": ("dg_solv",),
        "limitation_output_ids": ("dg",),
        "doubted_quantity_ids": ("dg_solv",),
        "unanswered_failed_verdicts": (),
        "stale_quantity_ids": (),
        "unclaimed_output_ids": (),
        "undelivered_declared_observable_ids": (),
        # The ninth gap (2026-09-09): a number delivered under its id
        # can still miss the precision the task asked for, and the wake
        # names that too, so the next action follows the gap.
        "unresolved_requirement_ids": (),
        "sufficiency": (),
        "flagged_quantity_ids": (),
        "failed_source_quantity_ids": (),
        "characterised_source_quantity_ids": (),
        "uncharacterised_source_quantity_ids": (),
        "expectation_rows": (),
    }


def test_goal_terms_are_restated_in_the_recency_slot(tmp_path):
    """Alphabetical canonical JSON lands the goal block mid-context --
    the attention trough. A goal session's coordinator message now ends
    with the goal terms restated verbatim; a plain session's message
    shape is untouched."""

    from chemsmart.agent.live_session import _coordinator_base_messages

    goal = {"goal_id": "goal-t1", "budgets": {"engine_calls_remaining": 2}}
    with_goal = _coordinator_base_messages(
        context={"task": "t", "goal": goal},
        approved_workflow=None,
    )
    without = _coordinator_base_messages(
        context={"task": "t"},
        approved_workflow=None,
    )

    assert len(without) == 2
    assert len(with_goal) == 3
    tail = with_goal[-1]
    assert tail["role"] == "user"
    assert "restated for recency" in tail["content"]
    assert '"goal_id":"goal-t1"' in tail["content"].replace(" ", "")


def test_a_typed_error_settles_instead_of_escaping(tmp_path):
    """C5's live crash: a session attempted terminal completion against
    a red gate, the ContractError escaped the loop, and the goal's
    durable story ended at wake_composed with no settlement. A typed
    contract error now settles returned_to_human naming the stage."""

    def raising_session(workspace, kwargs):
        raise ContractError("a required completion gate is red")

    result = _loop(
        tmp_path,
        sessions=[raising_session],
        executes=[],
    )

    assert result.settlement == "returned_to_human"
    assert any(
        "a required completion gate is red" in reason
        for reason in result.reasons
    )
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    assert ledger.entries()[-1]["kind"] == "goal_settled"


def test_a_typed_engine_error_settles_the_same_way(tmp_path):
    def raising_execute(run_directory):
        raise ContractError("approval bundle names a missing artifact")

    result = _loop(
        tmp_path,
        sessions=[_planning_session("live-1", review=_review_payload())],
        executes=[raising_execute],
    )

    assert result.settlement == "returned_to_human"
    assert any("approved execution" in reason for reason in result.reasons)


def test_an_engineless_cycle_settles_from_its_delivery(tmp_path):
    """C8's live crash: an admitted revision launched no engine, the
    run directory recorded no workflow run, and derive_run_outcome's
    ValueError escaped the loop unsettled. The cycle now settles from
    the run stream's typed delivery, exactly as a no-partition
    planning cycle does."""

    def analysis_only_execute(run_directory):
        store = RuntimeEventStore(
            run_directory / "events.jsonl", session_id="exec-1"
        )
        store.append(
            turn_id="t1",
            kind="session_started",
            payload={"phase": "route", "task_id": "t"},
        )
        return SimpleNamespace(status="partial", analysis_status="partial")

    result = _loop(
        tmp_path,
        sessions=[_planning_session("live-1", review=_review_payload())],
        executes=[analysis_only_execute],
    )

    assert result.settlement == "returned_to_human"
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    kinds = [entry["kind"] for entry in ledger.entries()]
    assert kinds[-1] == "goal_settled"
    run_rows = [
        entry for entry in ledger.entries() if entry["kind"] == "run_recorded"
    ]
    assert run_rows[-1]["payload"]["workflow_state"] == "analysis_only"


def test_a_refused_launch_is_named_by_the_settlement(tmp_path):
    """A revision reused a failed node's id, the node-workspace guard
    refused the launch, the executor kept the refusal in memory, and
    the goal settled naming only the session's terminal word. The
    refusal is now a run-stream event and the settlement quotes it."""

    def refused_execute(run_directory):
        store = RuntimeEventStore(
            run_directory / "events.jsonl", session_id="exec-1"
        )
        store.append(
            turn_id="t1",
            kind="session_started",
            payload={"phase": "route", "task_id": "t"},
        )
        store.append(
            turn_id="exec-refused-ts_b",
            kind="workflow_node_launch_refused",
            payload={
                "node_id": "ts_b",
                "program": "orca",
                "jobtype": "ts",
                "reason": "execution workspace already contains outputs",
            },
        )
        return SimpleNamespace(status="partial", analysis_status="partial")

    result = _loop(
        tmp_path,
        sessions=[_planning_session("live-1", review=_review_payload())],
        executes=[refused_execute],
    )

    assert result.settlement == "returned_to_human"
    assert result.reasons == (
        "node ts_b never launched: execution workspace already contains "
        "outputs",
    )
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    run_rows = [
        entry for entry in ledger.entries() if entry["kind"] == "run_recorded"
    ]
    assert run_rows[-1]["payload"]["stopped_by"] == list(result.reasons)


def test_a_refused_review_is_named_by_the_settlement(tmp_path):
    """A woken session planned twelve engine nodes against five
    remaining calls; the frontier called it approvable, the review
    builder refused it at session end into a log line, and the goal
    returned naming nothing. The refusal is now a session event and
    the settlement quotes it."""

    rows = [
        {"kind": "session_started", "payload": {}},
        {
            "kind": "execution_review_refused",
            "payload": {
                "workflow_id": "w",
                "reason": (
                    "scientific workflow exceeds bounded engine-call "
                    "budget: 12 nodes for 5 calls"
                ),
            },
        },
    ]
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="complete", wake_rows=rows)
        ],
        executes=[],
    )

    assert result.settlement == "returned_to_human"
    assert result.reasons == (
        "execution review refused: scientific workflow exceeds bounded "
        "engine-call budget: 12 nodes for 5 calls",
    )


def test_the_goals_first_declarations_ride_the_wake(tmp_path):
    """Session 1's expectations are the goal's; a woken session is
    seeded with them (the host keeps the first) and the ledger carries
    them from the cycle that declared them."""

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    declared = {
        "kind": "requested_observable_declared",
        "payload": {
            "observables": [
                {
                    "observable_id": "cis-barrier",
                    "unit": "kcal/mol",
                    "dimension": [1, 0, 0, 0, 0, 0],
                    "meaning": "syn barrier above anti",
                    "expectation_basis": "torsional barriers",
                    "expected_sign": "positive",
                    "expected_low": 3.0,
                    "expected_high": 8.0,
                }
            ],
            "declared_total": 1,
        },
    }
    _loop(
        tmp_path,
        sessions=[
            capture(
                _planning_session(
                    "live-1",
                    review=_review_payload(),
                    wake_rows=[
                        {"kind": "session_started", "payload": {}},
                        declared,
                    ],
                )
            ),
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
            _execute(tmp_path, failed=False, status="completed"),
        ],
    )

    assert contexts[0]["declared_observables"] == ()
    assert [
        item["observable_id"] for item in contexts[1]["declared_observables"]
    ] == ["cis-barrier"]
    assert contexts[1]["declared_observables"][0]["expected_high"] == 8.0
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    kinds = [entry["kind"] for entry in ledger.entries()]
    assert kinds.index("goal_created") < kinds.index("observables_declared")


def test_an_interrupted_local_run_resumes_through_wake(tmp_path):
    """A launcher timeout killed a goal mid-engine in its second cycle
    and nothing could resume it: wake knew only parked (dispatched)
    runs. The execute boundary is now a ledger entry naming the run and
    its bundle, and resume re-enters the execute phase with them."""

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    calls = []

    def killed(*, approval_file, workspace, run_directory):
        calls.append("killed")
        _engine_stream(tmp_path, run_directory, failed=False)
        raise KeyboardInterrupt

    def continued(*, approval_file, workspace, run_directory):
        calls.append(("continued", approval_file.name, run_directory.name))
        return SimpleNamespace(status="completed", analysis_status="")

    common = dict(
        execution_envelope_file=_envelope_file(tmp_path, 6),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: _planning_session(
            "live-1", review=_review_payload()
        )(workspace, kw),
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
    )
    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execute_bundle=killed,
        **common,
    )
    with pytest.raises(KeyboardInterrupt):
        driver.run()
    ledger = GoalLedger(workspace / ".chemsmart-agent" / "goals" / "goal-t1")
    kinds = [entry["kind"] for entry in ledger.entries()]
    assert "run_started" in kinds and "run_recorded" not in kinds

    resumed = GoalDriver.resume(
        workspace=workspace,
        goal_id="goal-t1",
        execute_bundle=continued,
        plan_session=common["plan_session"],
        resolve_review=common["resolve_review"],
    )
    assert resumed.phase == "execute"
    result = resumed.run()
    assert calls == ["killed", ("continued", "bundle.json", "cycle-1")]
    assert result.settlement == "achieved"
    with pytest.raises(ContractError, match="is settled"):
        GoalDriver.resume(workspace=workspace, goal_id="goal-t1")


def test_retained_intent_is_not_an_ending(tmp_path):
    """Four hess nodes a session retained as declared non-executable
    intent derived not_launched, and the settlement called them a state
    no revision can answer while every executable node had validated.
    The cycle's displayed review says which ids were retained."""

    workspace = tmp_path / "ws"
    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path, 6),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: None,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **_kw: None,
    )
    driver.cycles = 1
    assert driver._declared_non_executable_ids() == ()
    review = _review_payload()
    review["non_executable_node_ids"] = ["hess-a", "hess-b"]
    reviews = driver.goal_dir / "reviews"
    reviews.mkdir(parents=True)
    (reviews / "cycle-1.json").write_text(json.dumps(review))
    assert driver._declared_non_executable_ids() == ("hess-a", "hess-b")


def test_an_anomaly_belongs_to_the_goal_not_to_the_cycle(tmp_path):
    """Two live goals settled plain achieved over anomalies recorded a
    cycle earlier: the receipt died with its cycle. The ledger carries
    every anomaly, the wake re-seeds later hosts, the run directory
    hands them to the executor, and the word consults the ledger."""

    from chemsmart.agent.terminal_states import (
        PRIOR_ANOMALIES_FILE,
        derive_run_outcome,
        read_run_events,
    )

    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    def failed_with_anomaly(run_directory):
        _engine_stream(tmp_path, run_directory, failed=True)
        events = run_directory / "events.jsonl"
        node = derive_run_outcome(read_run_events(events)).nodes[0]
        session_id = json.loads(events.read_text().splitlines()[0])[
            "session_id"
        ]
        RuntimeEventStore(events, session_id=session_id).append(
            turn_id="exec-anomaly",
            kind="anomaly_observed",
            payload={
                "receipt_sha256": "c" * 64,
                "status": "unreplicated",
                "node_id": node.node_id,
                "signal_id": "stationary_point.unexpected_order",
                "record": {
                    "signal_id": "stationary_point.unexpected_order",
                    "status": "unreplicated",
                    "values": {"observed_imaginary_modes": 1},
                },
            },
        )
        return SimpleNamespace(status="partial", analysis_status="partial")

    seen_prior = {}

    def clean(run_directory):
        seen_prior["file"] = (run_directory / PRIOR_ANOMALIES_FILE).exists()
        _engine_stream(tmp_path, run_directory, failed=False)
        return SimpleNamespace(status="completed", analysis_status="")

    result = _loop(
        tmp_path,
        sessions=[
            capture(_planning_session("live-1", review=_review_payload())),
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[failed_with_anomaly, clean],
    )

    assert result.settlement == "achieved_with_observations"
    assert "stationary_point.unexpected_order" in result.reasons[0]
    assert seen_prior["file"] is True
    assert [a["signal_id"] for a in contexts[1]["anomalies"]] == [
        "stationary_point.unexpected_order"
    ]
    ledger = GoalLedger(
        tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1"
    )
    kinds = [entry["kind"] for entry in ledger.entries()]
    assert "anomalies_observed" in kinds


def test_the_completion_receipt_certifies_not_the_terminal_word(tmp_path):
    """A session delivered its claims, passed the gate, then left an
    analysis-only draft and ended 'planned'; the goal returned saying
    the gate had not passed. The receipt decides, and the posture at
    exit is a note on the settlement."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1", terminal="planned", wake_rows=_delivery_rows()
            ),
        ],
        executes=[],
    )
    assert result.settlement == "achieved"
    assert any("ended 'planned'" in reason for reason in result.reasons)


def test_a_dispatched_run_parks_the_goal_and_resumes_at_its_outcome(
    tmp_path,
):
    """The one human decision continues in its own run directory.

    With a scheduler in the loop the driver does not wait: it records
    the job it submitted, parks, and a later process -- the job's own
    tail, or a human's `agent wake` -- rebuilds the driver from the
    ledger at the outcome phase and settles from the run's durable
    record. Nothing here creates a second decision.
    """

    from chemsmart.agent.driver import (
        EXECUTION_RESULT_FILE,
        GoalDriver,
        GoalLoopResultV1,
    )

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    submitted: dict = {}

    def dispatch_run(**kwargs):
        submitted.update(kwargs)
        return SimpleNamespace(
            scheduler="SLURM",
            job_id="191",
            submitted_at="2026-09-01T00:00:00+00:00",
            submit_script=str(kwargs["run_directory"] / "sub.sh"),
        )

    def never_execute(**_kwargs):
        raise AssertionError("a dispatched run is not executed in-process")

    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-parked",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: _planning_session(
            "live-1", review=_review_payload()
        )(workspace, kw),
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=never_execute,
        dispatch_run=dispatch_run,
    )
    parked = driver.run()
    assert isinstance(parked, GoalLoopResultV1)
    assert parked.settlement == "parked"
    assert "job 191" in parked.reasons[0]
    assert submitted["cycle"] == 1
    kinds = [entry["kind"] for entry in driver.ledger.entries()]
    assert kinds[-1] == "run_dispatched"
    assert "goal_settled" not in kinds
    assert (driver.goal_dir / "task.md").read_text() == "the goal task"

    # A fresh process may not start the goal over ...
    with pytest.raises(ContractError, match="already exists"):
        GoalDriver(
            task="the goal task",
            workspace=workspace,
            execution_envelope_file=_envelope_file(tmp_path),
            goal_id="goal-parked",
            granted_by="claude-owner-delegated-reviewer",
        )

    # ... but it resumes it once the job has written the run's record.
    run_directory = driver.run_directory
    _engine_stream(tmp_path, run_directory, failed=False)
    (run_directory / EXECUTION_RESULT_FILE).write_text(
        json.dumps({"status": "completed", "analysis_status": ""}),
        encoding="utf-8",
    )

    def no_session(**_kw):
        raise AssertionError("no new session is needed to settle")

    resumed = GoalDriver.resume(
        workspace=workspace,
        goal_id="goal-parked",
        execution_envelope_file=_envelope_file(tmp_path),
        granted_by="claude-owner-delegated-reviewer",
        plan_session=no_session,
        execute_bundle=never_execute,
    )
    assert resumed.phase == "outcome"
    assert resumed.task == "the goal task"
    assert resumed.cycles == 1
    result = resumed.run()
    assert result.settlement == "achieved"
    recorded = [
        entry
        for entry in resumed.ledger.entries()
        if entry["kind"] == "run_recorded"
    ]
    assert recorded[-1]["payload"]["cycle"] == 1
    assert recorded[-1]["payload"]["queue_wait_seconds"] > 0.0

    # Settled goals do not resume, and a goal with no parked run does not.
    with pytest.raises(ContractError, match="is settled"):
        GoalDriver.resume(
            workspace=workspace,
            goal_id="goal-parked",
            execution_envelope_file=_envelope_file(tmp_path),
            granted_by="claude-owner-delegated-reviewer",
        )


def test_a_failed_run_opens_a_typed_recovery_with_a_repair_menu(tmp_path):
    """A timeout is ordinary work: the ledger says a recovery opened and
    names how each node ended, and the next wake carries the route that
    answers it. The host names the route; the physics grades it."""

    from chemsmart.agent.driver import REPAIR_MENU, GoalDriver

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    wakes: list = []

    def plan_session(**kwargs):
        wakes.append(kwargs["goal_context"])
        name = f"live-{len(wakes)}"
        rows = list(_READ_OUTCOME_ROWS) if len(wakes) > 1 else ()
        return _planning_session(
            name, review=_review_payload(), wake_rows=rows
        )(workspace, kwargs)

    executes = iter(
        (
            _execute(tmp_path, failed=True, status="failed"),
            _execute(tmp_path, failed=False, status="completed"),
        )
    )
    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-repair",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=plan_session,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **kw: next(executes)(kw["run_directory"]),
    )
    result = driver.run()
    assert result.settlement == "achieved"
    assert result.cycles == 2
    opened = [
        entry
        for entry in driver.ledger.entries()
        if entry["kind"] == "recovery_opened"
    ]
    assert len(opened) == 1
    assert opened[0]["payload"]["cycle"] == 1
    assert set(opened[0]["payload"]["terminal_states"].values()) == {
        "timeout_terminated"
    }
    second_wake = wakes[1]
    assert second_wake["repair_menu"] == {
        "timeout_terminated": REPAIR_MENU["timeout_terminated"]
    }
    assert "repair_menu names" in second_wake["authority"]
    assert wakes[0]["deliverables"]["delivered_quantity_ids"] == ()


def test_a_run_no_revision_can_answer_returns_to_the_human(tmp_path):
    """A launch that never happened is not evidence a revision can stand
    on; re-planning over it would spend budget on the host's own gap."""

    from chemsmart.agent.driver import GoalDriver

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)

    def execute(run_directory):
        _engine_stream(
            tmp_path,
            run_directory,
            failed=True,
            findings=("execution.process.launch_failed",),
        )
        return SimpleNamespace(status="failed", analysis_status="")

    sessions = 0

    def plan_session(**kwargs):
        nonlocal sessions
        sessions += 1
        return _planning_session("live-1", review=_review_payload())(
            workspace, kwargs
        )

    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-unanswerable",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=plan_session,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **kw: execute(kw["run_directory"]),
    )
    result = driver.run()
    assert result.settlement == "returned_to_human"
    assert "no revision can answer" in result.reasons[0]
    assert "launch_failed" in result.reasons[0]
    assert sessions == 1, "no second planning session was started"


def _driver_after_run(tmp_path, *, calls, rows, failed, settle=True):
    """A driver standing at its settle step after one recorded run whose
    stream carries the delivery rows, with the engine line spent.

    ``settle=False`` leaves the goal unsettled, for a caller that must
    put something on the ledger first. A goal settles once --
    `GoalLedger.settle` is keyed and `resume` refuses a settled goal --
    so a caller that needs a different word builds the ledger it wants
    and settles once, rather than settling twice and reading the
    second."""

    workspace = tmp_path / "ws"
    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path, calls),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: None,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **_kw: None,
    )
    driver.cycles = 1
    driver.goal = driver._goal_record(
        identity="b" * 64,
        conditions={"solvents": (), "thermochemistry": ((298.15, 1.0, None),)},
        review_sha256="c" * 64,
    )
    driver.ledger.create(driver.goal)
    driver.ledger.append(
        "run_recorded",
        {
            "cycle": 1,
            "engine_calls_consumed": calls,
            "engine_wall_seconds": 10.0,
            "workflow_state": "failed" if failed else "validated",
        },
    )
    driver.run_directory = driver.goal_dir / "runs" / "cycle-1"
    driver.run_directory.mkdir(parents=True)
    (driver.run_directory / "events.jsonl").write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    driver.execute_result = SimpleNamespace(
        status="failed" if failed else "completed",
        analysis_status="partial" if failed else "completed",
    )
    if settle:
        driver._settle()
    return driver


def test_a_certified_delivery_never_says_achieved_over_an_undelivered_declared_observable(
    tmp_path,
):
    """Two live goals settled achieved_with_observations on a run stream
    whose completion passed while listing declared observables no claim
    carried by id (E4 window, 2026-09-03): the limitation branch wanted
    decisions in the same stream and the executor's stream has none. A
    declared headline nobody claimed is not a delivery, whatever word
    the chain's kernels earned; and it is not a typed refusal either."""

    result = _loop(
        tmp_path,
        sessions=[
            _planning_session(
                "live-1",
                terminal="complete",
                wake_rows=_delivery_rows(
                    limitations=("declared_observable:nh3-gibbs-298",)
                ),
            ),
        ],
        executes=[],
    )
    assert result.settlement == "returned_to_human"
    assert any("nh3-gibbs-298" in reason for reason in result.reasons)


def test_a_claim_by_id_recovery_costs_no_engine_call(tmp_path):
    """A goal with every receipt on disk and its engine line spent
    settled exhausted because the recovery test demanded an engine
    call; claiming a declared observable by id from receipts in hand
    needs none. The cycle wakes with zero engine calls, and the
    plan-time budget gate refuses any engine node it might plan."""

    driver = _driver_after_run(
        tmp_path,
        calls=1,
        failed=False,
        rows=_delivery_rows(
            limitations=("declared_observable:gibbs-rrho-298",)
        ),
    )
    assert driver.phase == "plan"
    (opened,) = [
        entry
        for entry in driver.ledger.entries()
        if entry["kind"] == "recovery_opened"
    ]
    assert opened["payload"]["undelivered_declared_observable_ids"] == [
        "gibbs-rrho-298"
    ]
    assert opened["payload"]["engine_calls_remaining"] == 0


def test_a_failed_run_with_receipts_in_hand_is_not_exhausted(tmp_path):
    """butadiene-a: the pool was spent on a saddle, a scan and a
    validated minimum, four declared observables sat unclaimed with
    their receipts on disk, and the word was exhausted. An analysis-only
    cycle can still claim them by id."""

    driver = _driver_after_run(
        tmp_path,
        calls=1,
        failed=True,
        rows=_delivery_rows(
            completion="partial",
            limitations=("declared_observable:rel-energy-scis",),
        ),
    )
    assert driver.phase == "plan"
    (opened,) = [
        entry
        for entry in driver.ledger.entries()
        if entry["kind"] == "recovery_opened"
    ]
    assert opened["payload"]["analysis_only"] is True
    assert opened["payload"]["undelivered_declared_observable_ids"] == [
        "rel-energy-scis"
    ]


def test_the_goal_keeps_every_budget_line_the_envelope_granted(tmp_path):
    """Two hand-listed copies of the goal's envelope record dropped the
    excursion line: every woken cycle of every excursion arm in two
    sealed windows read zero excursions remaining and the plan gate
    refused the line the human had granted (E4, E4', 2026-09-03)."""

    from chemsmart.agent.driver import _wake_context

    driver = GoalDriver(
        task="the goal task",
        workspace=tmp_path / "ws",
        execution_envelope_file=_envelope_file(tmp_path, 3, excursions=2),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: None,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **_kw: None,
    )
    assert driver.envelope_record["max_excursion_calls"] == 2
    goal = driver._goal_record(
        identity="b" * 64,
        conditions={"solvents": (), "thermochemistry": ((298.15, 1.0, None),)},
        review_sha256="c" * 64,
    )
    driver.ledger.create(goal)
    budgets = driver.ledger.budgets(goal)
    assert budgets.engine_calls_remaining == 3
    assert budgets.excursion_calls_remaining == 2
    wake = _wake_context(goal, driver.ledger, None, workspace=tmp_path / "ws")
    assert wake["budgets"]["excursion_calls_remaining"] == 2


def test_the_word_names_the_delivered_number_that_came_from_a_flagged_node(
    tmp_path,
):
    """Naming the anomaly is not naming the number. Five deliveries in
    two cyclohexane windows reported a structure the basin sensor had
    flagged as "the minimum", and the settlement said only that
    something somewhere had been observed. The host knows which claim
    descends from the flagged result: the anomaly names that node's
    artifacts, and the delivery walk already follows claims back to
    artifacts."""

    artifact = "9" * 64
    rows = [
        {
            "kind": "result_quantities_extracted",
            "payload": {
                "receipt_sha256": "e" * 64,
                "artifact_sha256": artifact,
                "status": "extracted",
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "a1" + "a" * 62,
                "status": "recorded",
                "record": {
                    "claims": (
                        {
                            "source_receipt_sha256": "e" * 64,
                            "quantity_id": "gibbs-298",
                        },
                    )
                },
            },
        },
        {
            "kind": "analysis_completion_evaluated",
            "payload": {
                "receipt_sha256": "c1" + "c" * 62,
                "status": "passed",
                "limitation_output_ids": [],
            },
        },
    ]
    driver = _driver_after_run(
        tmp_path, calls=1, failed=False, rows=rows, settle=False
    )
    driver.ledger.append(
        "anomalies_observed",
        {
            "cycle": 1,
            "anomalies": [
                {
                    "receipt_sha256": "d" * 64,
                    "signal_id": "geometry.heavy_atom_rmsd_ge_0.3",
                    "status": "unreplicated",
                    "node_id": "opt",
                    "flagged_artifact_sha256s": [artifact],
                }
            ],
        },
    )
    # The anomalies are on the ledger before the single settlement, so
    # the word is produced once rather than by re-settling a goal that
    # has already ended.
    driver.phase = "settle"
    driver._settle()

    settled = driver.ledger.entries()[-1]
    assert settled["payload"]["state"] == "achieved_with_observations"
    reasons = " ".join(settled["payload"]["reasons"])
    assert "geometry.heavy_atom_rmsd_ge_0.3" in reasons
    assert "delivered from the flagged result: gibbs-298" in reasons


def _failed_source_rows(*, characterised: bool):
    """A run whose delivered number was read off a node that missed its
    promise, optionally with the host asked what that structure is."""

    artifact = "7" * 64
    rows = [
        {
            "kind": "program_result_verified",
            "payload": {
                "receipt_sha256": "b" * 64,
                "node_id": "opt-a",
                "record": {
                    "state": "invalid",
                    "output_artifacts": [{"sha256": artifact}],
                },
            },
        },
        {
            "kind": "result_quantities_extracted",
            "payload": {
                "receipt_sha256": "e" * 64,
                "artifact_sha256": artifact,
                "status": "extracted",
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "a1" + "a" * 62,
                "status": "recorded",
                "record": {
                    "claims": (
                        {
                            "source_receipt_sha256": "e" * 64,
                            "quantity_id": "barrier-kcal",
                        },
                    )
                },
            },
        },
        {
            "kind": "analysis_completion_evaluated",
            "payload": {
                "receipt_sha256": "c1" + "c" * 62,
                "status": "passed",
                "limitation_output_ids": [],
            },
        },
    ]
    if characterised:
        rows.insert(
            1,
            {
                "kind": "stationary_point_characterised",
                "payload": {
                    "receipt_sha256": "d" * 64,
                    "node_id": "opt-a",
                    "order_claimed": 1,
                    "record": {"result_artifact_sha256": artifact},
                },
            },
        )
    return rows


def test_a_number_read_off_a_failed_node_is_named_by_the_settlement(tmp_path):
    """Twenty-four archived saddles were readable all along and no
    settlement ever said a delivered number came from one. The word says
    it now, and says nothing about whether that was the right thing to
    do."""

    driver = _driver_after_run(
        tmp_path,
        calls=1,
        failed=False,
        rows=_failed_source_rows(characterised=False),
    )
    settled = driver.ledger.entries()[-1]
    reasons = " ".join(settled["payload"]["reasons"])
    assert settled["payload"]["state"] == "achieved"
    assert "did not meet its promise, uncharacterised: barrier-kcal" in reasons


def test_a_characterised_source_says_so_instead(tmp_path):
    """With the host asked what the structure is, the same delivery reads
    as a checked statement rather than an unexamined one."""

    driver = _driver_after_run(
        tmp_path,
        calls=1,
        failed=False,
        rows=_failed_source_rows(characterised=True),
    )
    reasons = " ".join(driver.ledger.entries()[-1]["payload"]["reasons"])
    assert "delivered from a characterised result: barrier-kcal" in reasons
    assert "uncharacterised" not in reasons


def test_a_denial_holds_when_the_first_cycle_only_read_results(tmp_path):
    """A goal record's existence is not an execution grant.

    An analysis-only first cycle creates the goal record in `_plan`
    with an empty initial review, so a later cycle's first executable
    review found `self.goal` already set, took the revision path, and
    resolved its own review with decision="approve". The initial
    decision was never consulted. Observed live before this repair:
    `--initial-decision deny`, one engine partition launched, settled
    `achieved`. The gate now reads the grant the human gave.
    """

    from chemsmart.agent.driver import run_goal_loop

    analysis_only = [
        {
            "kind": "requested_observable_declared",
            "payload": {
                "observables": [
                    {
                        "observable_id": "dg_solv",
                        "unit": "kJ/mol",
                        "meaning": "solvation free energy",
                        "dimension": [1, 0, 0, 0, 0, 0],
                    }
                ]
            },
        },
        {
            "kind": "result_quantities_extracted",
            "payload": {"receipt_sha256": "e" * 64},
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "q",
                            "quantity_id": "q",
                            "display_value": 1.0,
                            "display_unit": "kJ/mol",
                            "dimension": [1, 0, 0, 0, 0, 0],
                            "source_receipt_sha256": "4" * 64,
                        }
                    ]
                },
            },
        },
        {"kind": "scientific_decision_recorded", "payload": {}},
    ]
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    sessions = iter(
        [
            _planning_session(
                "live-1", terminal="planned", wake_rows=analysis_only
            ),
            _planning_session(
                "live-2", review=_review_payload(), wake_rows=analysis_only
            ),
            _planning_session(
                "live-3", terminal="planned", wake_rows=analysis_only
            ),
        ]
    )
    executes = iter([_execute(tmp_path, failed=False, status="completed")])
    decisions: list[str] = []
    launched: list[str] = []

    def plan_session(**kwargs):
        return next(sessions)(workspace, kwargs)

    def resolve_review(**kwargs):
        decisions.append(str(kwargs.get("decision")))
        return ("d" * 64, tmp_path / "bundle.json")

    def execute_bundle(*, approval_file, workspace, run_directory):
        launched.append(str(run_directory))
        return next(executes)(run_directory)

    result = run_goal_loop(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path, 6),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        max_revisions=3,
        plan_session=plan_session,
        resolve_review=resolve_review,
        execute_bundle=execute_bundle,
        initial_decision="deny",
    )

    assert result.settlement == "returned_to_human"
    assert decisions == []
    assert launched == []


def test_a_typed_error_still_records_what_the_cycle_delivered(tmp_path):
    """Surviving an error and preserving what it interrupted differ.

    SUFFICIENCY-2's session recorded 57 claims, 11 declarations, a
    sufficiency assessment and a scientific decision, then hit a red
    completion gate on its terminal event. `_typed_error_settlement`
    caught it and the goal settled -- and the ledger held one line and
    the workspace record none, because the projection runs after the
    planning session returns and the error returned first. The evidence
    was never destroyed; it was made unreachable, which is the same
    thing to every later reader.
    """

    import json as _json

    from chemsmart.agent._contracts import ContractError
    from chemsmart.agent.driver import run_goal_loop

    rows = [
        {
            "kind": "requested_observable_declared",
            "payload": {
                "observables": [
                    {
                        "observable_id": "dg_solv",
                        "unit": "kJ/mol",
                        "meaning": "solvation free energy",
                        "dimension": [1, 0, 0, 0, 0, 0],
                    }
                ]
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "dg_solv",
                            "quantity_id": "dg_solv",
                            "display_value": 1.0,
                            "display_unit": "kJ/mol",
                            "dimension": [1, 0, 0, 0, 0, 0],
                            "source_receipt_sha256": "4" * 64,
                        }
                    ]
                },
            },
        },
    ]
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)

    def plan_session(**kwargs):
        _write_session_stream(workspace, "live-1", rows)
        raise ContractError("a required completion gate is red")

    result = run_goal_loop(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path, 6),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        max_revisions=3,
        plan_session=plan_session,
        resolve_review=lambda **kw: ("d" * 64, tmp_path / "bundle.json"),
        execute_bundle=lambda **kw: None,
    )
    assert result.settlement == "returned_to_human"

    ledger = (
        workspace / ".chemsmart-agent" / "goals" / "goal-t1" / "ledger.jsonl"
    )
    kinds = [
        _json.loads(line)["kind"]
        for line in ledger.read_text(encoding="utf-8").splitlines()
    ]
    # A ledger holding only its own settlement is a malformed story.
    assert "goal_created" in kinds
    assert "observables_declared" in kinds

    record = workspace / ".chemsmart-agent" / "workspace-record.jsonl"
    claims = [
        _json.loads(line)
        for line in record.read_text(encoding="utf-8").splitlines()
        if _json.loads(line).get("kind") == "claim"
    ]
    assert [row["claim_id"] for row in claims] == ["dg_solv"]


def test_a_re_woken_cycle_records_what_it_delivered(tmp_path):
    """The one cycle the mechanism exists to produce was the one lost.

    A cycle that re-wakes returned before the projection, so its
    delivery never reached the record. SUFFICIENCY-3 lost 26 rows that
    way, including its own `attested` assessment -- and had its next
    cycle delivered anything else, the goal-grain join would have
    fallen back to a staler, worse row from two cycles earlier.
    """

    import json as _json

    rows = [
        {
            "kind": "requested_observable_declared",
            "payload": {
                "observables": [
                    {
                        "observable_id": "dg_solv",
                        "unit": "kJ/mol",
                        "meaning": "solvation free energy",
                        "dimension": [1, 0, 0, 0, 0, 0],
                        "required_tolerance": 2.0,
                        "tolerance_basis": "the author asked for 2 kJ/mol",
                    }
                ]
            },
        },
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "dg_solv",
                            "quantity_id": "dg_solv",
                            "display_value": 1.0,
                            "display_unit": "kJ/mol",
                            "dimension": [1, 0, 0, 0, 0, 0],
                            "source_receipt_sha256": "4" * 64,
                            "uncertainty": 1.0,
                            "uncertainty_basis": "asserted",
                        }
                    ]
                },
                "sufficiency": [
                    {
                        "observable_id": "dg_solv",
                        "unit": "kJ/mol",
                        "required_tolerance": 2.0,
                        "uncertainty": 1.0,
                        "uncertainty_basis": "asserted",
                        "uncertainty_evidence_backed": False,
                        "meets_tolerance": True,
                        "state": "attested",
                    }
                ],
            },
        },
        {"kind": "scientific_decision_recorded", "payload": {}},
    ]
    result = _loop(
        tmp_path,
        sessions=[
            _planning_session("live-1", terminal="planned", wake_rows=rows),
            _planning_session("live-2", terminal="planned", wake_rows=rows),
            _planning_session("live-3", terminal="planned", wake_rows=rows),
        ],
        executes=[],
        max_revisions=3,
    )
    assert result.settlement

    record = tmp_path / "ws" / ".chemsmart-agent" / "workspace-record.jsonl"
    claims = [
        _json.loads(line)
        for line in record.read_text(encoding="utf-8").splitlines()
        if _json.loads(line).get("kind") == "claim"
    ]
    # The re-woken cycle's own delivery, and its assessment, are on the
    # record: one row per cycle that delivered, and each cycle recorded
    # once rather than once per path that projects.
    cycles = sorted(row["cycle"] for row in claims)
    assert cycles == sorted(set(cycles))
    assert len(cycles) > 1, "the re-woken cycle projected nothing"
    assert {row["claim_id"] for row in claims} == {"dg_solv"}
    assert all(row["sufficiency"]["state"] == "attested" for row in claims)


def test_a_transport_loss_does_not_block_the_requirement_wake(tmp_path):
    """The previous cycle's failure report belongs to the previous cycle.

    `self.failure_report` is set when a wake opens and was never
    cleared, so the guard in `_rewake` that reads it blocked every later
    wake unconditionally. A12 taught the *ledger* scan to ignore a
    transport continuation and left this reader, two lines above it,
    asking the same question by different means. Observed live:
    SUFFICIENCY-4 arm A lost cycle 1 to four inter-event timeouts,
    delivered an `attested` requirement at cycle 2, and settled with
    forty engine calls and every revision unspent -- so the arm that
    was supposed to receive the sufficiency consequence never did, and
    the window was void.
    """

    import json as _json

    from chemsmart.agent.driver import run_goal_loop
    from chemsmart.agent.terminal_states import (
        PROVIDER_TRANSPORT_TERMINAL_REASON,
    )

    declared = {
        "kind": "requested_observable_declared",
        "payload": {
            "observables": [
                {
                    "observable_id": "dg_solv",
                    "unit": "kJ/mol",
                    "meaning": "solvation free energy",
                    "dimension": [1, 0, 0, 0, 0, 0],
                    "required_tolerance": 2.0,
                    "tolerance_basis": "the author asked for 2 kJ/mol",
                }
            ]
        },
    }
    attested = [
        declared,
        {
            "kind": "analysis_claims_recorded",
            "payload": {
                "receipt_sha256": "3" * 64,
                "record": {
                    "claims": [
                        {
                            "claim_id": "dg_solv",
                            "quantity_id": "dg_solv",
                            "display_value": 1.0,
                            "display_unit": "kJ/mol",
                            "dimension": [1, 0, 0, 0, 0, 0],
                            "source_receipt_sha256": "4" * 64,
                            "uncertainty": 1.0,
                            "uncertainty_basis": "asserted",
                        }
                    ]
                },
                "sufficiency": [
                    {
                        "observable_id": "dg_solv",
                        "unit": "kJ/mol",
                        "required_tolerance": 2.0,
                        "uncertainty": 1.0,
                        "uncertainty_basis": "asserted",
                        "uncertainty_evidence_backed": False,
                        "meets_tolerance": True,
                        "state": "attested",
                    }
                ],
            },
        },
        {"kind": "scientific_decision_recorded", "payload": {}},
    ]
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    sessions = iter(
        [
            _planning_session(
                "live-1",
                terminal="failed",
                wake_rows=[
                    declared,
                    {
                        "kind": "runtime_terminated",
                        "payload": {
                            "reason": PROVIDER_TRANSPORT_TERMINAL_REASON,
                            "terminal_state": "failed",
                        },
                    },
                ],
            ),
        ]
        + [
            _planning_session(
                f"live-{index}", terminal="planned", wake_rows=attested
            )
            for index in range(2, 6)
        ]
    )
    run_goal_loop(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path, 6),
        goal_id="goal-t1",
        granted_by="claude-owner-delegated-reviewer",
        max_revisions=5,
        plan_session=lambda **kwargs: next(sessions)(workspace, kwargs),
        resolve_review=lambda **kwargs: ("d" * 64, tmp_path / "b.json"),
        execute_bundle=lambda **kwargs: None,
    )
    gates = [
        (
            (_json.loads(line)["payload"].get("failure_report") or {}).get(
                "gate"
            ),
            _json.loads(line)["payload"].get("transport_continuation"),
        )
        for line in (
            workspace
            / ".chemsmart-agent"
            / "goals"
            / "goal-t1"
            / "ledger.jsonl"
        )
        .read_text(encoding="utf-8")
        .splitlines()
        if _json.loads(line)["kind"] == "rewake_opened"
    ]
    assert ("goal.cycle_delivers_or_returns", True) in gates
    assert any(
        gate == "goal.requirement_is_resolved" for gate, _ in gates
    ), f"the requirement wake never fired: {gates}"


def test_a_refusal_in_the_decide_phase_still_settles_the_goal(tmp_path):
    """E1 of PySCF round 2 (2026-09-13): the session planned, the host
    built the review, and the resolver refused the driver's own approval
    id inside the decide phase before the ledger existed; the exception
    escaped the loop and the process exited with no goal record and no
    settlement. Every ending is a settlement: a ContractError in any phase
    settles the goal returned_to_human, naming the phase and the refusal,
    and the goal id the human typed is normalised through the one
    identifier rule so the directory, the ledger and the approval carry
    one spelling."""

    def refuse(**kwargs):
        raise ContractError("execution bundle approval IDs differ")

    result = _loop(
        tmp_path,
        sessions=[_planning_session("live-1", review=_review_payload())],
        executes=[],
        goal_id="goal-T1-Upper",
        resolve=refuse,
    )
    assert result.goal_id == "goal-t1-upper"
    assert result.settlement == "returned_to_human"
    assert any(
        "decide phase" in reason and "approval IDs differ" in reason
        for reason in result.reasons
    ), result.reasons
    goal_dir = tmp_path / "ws" / ".chemsmart-agent" / "goals" / "goal-t1-upper"
    ledger = [
        json.loads(line)
        for line in (goal_dir / "ledger.jsonl").read_text().splitlines()
        if line.strip()
    ]
    kinds = [entry["kind"] for entry in ledger]
    assert "goal_settled" in kinds, kinds
    assert (goal_dir / "goal.json").exists()
