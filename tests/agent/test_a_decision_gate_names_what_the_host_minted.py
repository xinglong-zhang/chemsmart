"""A decision gate tells the truth about a digest the host minted.

``record_scientific_decision`` refuses a postprocessing citation that is not
a receipt the host minted (``decision.receipt_is_one_the_host_minted``). It
read this session's registries and the goal runs' streams for five receipt
kinds, and said of everything else "no digest this host minted". R11
`truth` (census 3): of 71 such refusals the Agent met, 29 cited a receipt
the host had minted -- 22 in an earlier run's stream (10 result
verifications, 12 anomaly observations) and 3 in an earlier planning
session, one of them the failed validation receipt the wake tells a
session to cite -- and every one was told the digest was never minted.

Driven through the goal loop with real hosts: (a) cycle 1's planning
session walks an analysis-only chain over the archived PySCF O2 singlet
bytes whose stability criterion fails, and cycle 2's woken session cites
that session's failed receipt, as the wake prescribes; (b) cycle 1's run
records an anomaly (appended through the run stream's own event store, as
the goal loop's anomaly test does), and cycle 2's woken session cites it
in the postprocessing field, then where the host says it belongs.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent.driver import run_goal_loop
from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.runtime.alibaba import (
    Qwen38MaxConfigV1,
    Qwen38MaxToolSession,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

from .provider_fakes import _run_contracts
from .test_a_failed_criterion_is_a_finding_the_goal_can_deliver import (
    _TASK,
    _artifact,
    _stream_rows,
)
from .test_a_partial_delivery_ends_its_session import _call, _o2r_turns, _turn
from .test_the_goal_loop_recovers_or_returns import (
    _bundle_file,
    _engine_stream,
    _envelope_file,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability(
    "rule:wake.failed_validation_receipt_answers_verdict"
)

_GATE = "decision.receipt_is_one_the_host_minted"


def _agent(workspace):
    return workspace / ".chemsmart-agent"


def _session_stream(workspace, run_id):
    return _agent(workspace) / "runs" / run_id / "events.jsonl"


def _real_turns(stream, turns, *, workspace, scratch, prior_anomalies=()):
    """A real host driven through the tool loop into ``stream``."""

    artifact = _artifact()
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(stream, session_id="protocol-session"),
        task_spec_sha256s=(_TASK,),
        approved_workspace=scratch / "approved",
        execute_analysis_only_plans=True,
        analysis_only_run_directory=scratch / f"analysis-{stream.parent.name}",
        analysis_only_workspace=scratch / "approved",
        run_evidence_root=workspace,
        # live_session seeds a woken host with the goal-wide anomalies its
        # wake context carries; this is that one line.
        prior_anomaly_observations=tuple(prior_anomalies),
    )
    host.artifacts[artifact.artifact_id] = artifact
    pending = iter(turns(artifact.artifact_id))
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda payload: next(pending)(payload),
        messages=[{"role": "user", "content": "Is the reference stable?"}],
        config=config,
    )
    envelope, request_context, network = _run_contracts(host, config)
    return ToolLoopRunner(host=host, event_store=host.event_store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
        should_stop=lambda: False,
    )


def _session(run_id, turns, *, review, seeded=False):
    """A goal session: real turns into its own stream; with ``review`` it
    plans a run and selects one wave, without it the cycle ends there."""

    def step(workspace, kwargs):
        anomalies = (
            (kwargs.get("goal_context") or {}).get("anomalies") or ()
            if seeded
            else ()
        )
        _real_turns(
            _session_stream(workspace, run_id),
            turns,
            workspace=workspace,
            scratch=workspace.parent / f"scratch-{run_id}",
            prior_anomalies=anomalies,
        )
        from chemsmart.agent.cohort import build_execution_wave_decision

        selected = ()
        if review:
            review_file = kwargs["review_file"]
            review_file.parent.mkdir(parents=True, exist_ok=True)
            review_file.write_text(
                json.dumps(_review_payload()), encoding="utf-8"
            )
            selected = ("sp-initial",)
        return SimpleNamespace(
            terminal_state="waiting_for_approval" if review else "complete",
            run_id=run_id,
            task_spec_sha256=_TASK,
            selected_execution_wave=selected,
            execution_wave_decision=build_execution_wave_decision(
                state="selected" if selected else "undecided",
                workflow_id="water-workflow" if selected else "",
                ready_node_ids=selected,
                node_ids=selected,
            ),
        )

    return step


def _goal(tmp_path, *, sessions, executes, goal_id):
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    session_iter = iter(sessions)
    execute_iter = iter(executes)

    def plan_session(**kwargs):
        return next(session_iter)(workspace, kwargs)

    def resolve_review(**_kwargs):
        return ("d" * 64, _bundle_file(tmp_path))

    def execute_bundle(*, approval_file, workspace, run_directory):
        return next(execute_iter)(run_directory)

    return run_goal_loop(
        task="Is the restricted reference of singlet O2 stable?",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id=goal_id,
        granted_by="claude-researcher-truth-owner-delegated",
        max_revisions=2,
        plan_session=plan_session,
        resolve_review=resolve_review,
        execute_bundle=execute_bundle,
    )


def _decision(ordinal, *, postprocessing=(), evidence=()):
    return _call(
        ordinal,
        "record_scientific_decision",
        {
            "decision_id": f"o2-rks-read-{ordinal}",
            "assumptions": ["the restricted reference at the geometry"],
            "method_rationale": "the task fixed B3LYP/def2-SVP",
            "alternatives": ["a broken-symmetry solution, not asked"],
            "uncertainties": ["SCF convergence precision"],
            "diagnostics": ["what the earlier cycle recorded"],
            "stage_order": ["read", "decide"],
            "evidence_refs": list(evidence),
            "postprocessing_receipt_sha256s": list(postprocessing),
        },
    )


def _failed_rows(rows):
    return [
        row["payload"]
        for row in rows
        if row["kind"] == "tool_failed"
        and (row["payload"].get("failure_report") or {}).get("gate") == _GATE
    ]


def _cited(rows):
    return {
        reference
        for row in rows
        if row["kind"] == "scientific_decision_recorded"
        for reference in row["payload"]["record"]["evidence_refs"]
    }


def test_a_woken_decision_may_cite_what_an_earlier_session_minted(tmp_path):
    """r10/q22 gh2: the failed validation receipt was minted by an earlier
    planning session, the wake told the session to cite it, and the gate
    said the host had never minted it."""

    goal_id = "goal-o2r-sessions"
    first = "live-20260928T000000000000Z-truth-first"
    woken = "live-20260928T010000000000Z-truth-woken"
    workspace = tmp_path / "ws"

    def failed_receipt():
        return [
            row["payload"]["receipt_sha256"]
            for row in _stream_rows(_session_stream(workspace, first))
            if row["kind"] == "scientific_validation_evaluated"
            and not row["payload"]["all_rules_passed"]
        ][-1]

    def cite(_artifact_id):
        return [
            lambda payload: _turn(
                1,
                "The earlier session's criterion failed; that is read.",
                (_decision(1, postprocessing=(failed_receipt(),)),),
            ),
            lambda payload: _turn(2, "Recorded."),
        ]

    def judged(artifact_id):
        # o2r's acts without its decision: the criterion fails in this
        # session's own walk and nothing answers it here.
        turns = _o2r_turns(artifact_id, cite_verdict=False)
        return turns[:3] + turns[4:]

    _goal(
        tmp_path,
        goal_id=goal_id,
        sessions=[
            _session(first, judged, review=True),
            _session(woken, cite, review=False),
        ],
        executes=[
            lambda run_directory: (
                _engine_stream(tmp_path, run_directory, failed=True),
                SimpleNamespace(status="partial", analysis_status=""),
            )[1]
        ],
    )

    rows = _stream_rows(_session_stream(workspace, woken))
    digest = failed_receipt()
    assert not _failed_rows(rows), _failed_rows(rows)
    assert f"receipt:{digest}" in _cited(rows)


def test_a_refused_anomaly_citation_is_told_what_it_is_and_where_it_goes(
    tmp_path,
):
    """An anomaly receipt a run recorded, cited as postprocessing evidence:
    refused, because this field takes no anomaly -- and told so, with the
    stream that recorded it and the reference the host accepts, which the
    session then uses."""

    from chemsmart.agent.terminal_states import (
        derive_run_outcome,
        read_run_events,
    )

    goal_id = "goal-anomaly-cited"
    woken = "live-20260928T010000000000Z-truth-anomaly"
    anomaly = "c1" * 32
    workspace = tmp_path / "ws"

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
                "receipt_sha256": anomaly,
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

    def cite(_artifact_id):
        return [
            lambda payload: _turn(
                1,
                "Citing the anomaly the run recorded.",
                (_decision(1, postprocessing=(anomaly,)),),
            ),
            lambda payload: _turn(
                2,
                "Citing it where the host says it belongs.",
                (_decision(2, evidence=(f"anomaly:{anomaly}",)),),
            ),
            lambda payload: _turn(3, "Recorded."),
        ]

    _goal(
        tmp_path,
        goal_id=goal_id,
        sessions=[
            _planning_session(
                "live-20260928T000000000000Z-truth-plan",
                review=_review_payload(),
            ),
            _session(woken, cite, review=False, seeded=True),
        ],
        executes=[failed_with_anomaly],
    )

    rows = _stream_rows(_session_stream(workspace, woken))
    (refused,) = _failed_rows(rows)
    diagnosis = refused["failure_report"]["diagnosis"]
    # What it is, where it was recorded, and where it may be cited.
    assert "anomaly_observed" in diagnosis, diagnosis
    assert f"goals/{goal_id}/runs/cycle-1" in diagnosis, diagnosis
    assert f"anomaly:{anomaly}" in diagnosis, diagnosis
    assert "no digest this host minted" not in diagnosis
    # And that route is walkable.
    assert f"anomaly:{anomaly}" in _cited(rows)
