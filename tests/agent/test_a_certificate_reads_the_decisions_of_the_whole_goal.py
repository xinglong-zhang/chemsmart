"""A completion certificate reads the verdicts and decisions of the whole goal.

A certificate names every claim that stands on a failed acceptance
criterion no recorded decision has answered, and the settlement trusts
what it names. The settlement asks whether a verdict is answered of every
stream the goal's ledger names (``goal.failed_criteria`` over the goal's
streams); the certificate asked it of its own host alone -- the decisions
that host recorded and the validation receipts it minted. A verdict is
answered by a decision in one stream and judged again in another all the
time: the wake tells a woken session to cite the failed receipt of the run
it was woken with, and the next run, or the session's own analysis-only
walk, judges the same criterion over the same result again. The judging
host then held a receipt nobody in it had cited, minted a partial
certificate, and the goal returned to the human over a verdict its own
records answered (R11 truth-3, item 1: truth-2's dropped witness arm, and
the same class in the session's walk).

Driven through the goal loop with real hosts over the archived bytes of a
PySCF run of closed-shell singlet O2, whose restricted reference is
unstable to spin symmetry breaking. The woken session of the second arm is
the production session, ``run_live_agent_session``, with only the provider
transport replaced by a script, so the host it builds is the one a live
goal gets.
"""

from __future__ import annotations

import json
import shutil

import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.driver import run_goal_loop

from .test_a_failed_criterion_is_a_finding_the_goal_can_deliver import (
    _RESULT,
    _RULE,
    _run_with_the_failed_criterion,
    _stream_rows,
)
from .test_a_goal_word_reads_every_criterion_the_goal_holds import (
    _GOAL,
    _agent,
    _goal,
    _run_stream,
    _woken,
)
from .test_a_partial_delivery_ends_its_session import _o2r_turns
from .test_a_refused_input_says_why import _fenced_host, _ScriptedModel
from .test_the_goal_loop_recovers_or_returns import (
    _bundle_file,
    _envelope_file,
)

pytestmark = pytest.mark.capability(
    "rule:wake.failed_validation_receipt_answers_verdict"
)


def _failed_receipts(stream):
    return [
        row["payload"]["receipt_sha256"]
        for row in _stream_rows(stream)
        if row["kind"] == "scientific_validation_evaluated"
        and not row["payload"]["all_rules_passed"]
    ]


def _completions(stream):
    return [
        row["payload"]
        for row in _stream_rows(stream)
        if row["kind"] == "analysis_completion_evaluated"
    ]


def test_a_run_that_judges_an_answered_verdict_again_certifies_it(tmp_path):
    """Cycle 2's woken session cites the failed receipt cycle 1's run
    minted; cycle 2's run judges the same criterion over the same result
    again. The verdict is answered in the goal, so the run's certificate
    certifies the delivery and carries the verdict, and the word says so.
    """

    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    result = _goal(
        tmp_path,
        cite=True,
        second_run=_run_with_the_failed_criterion(tmp_path / "second"),
    )

    workspace = tmp_path / "ws"
    kinds = [
        row["kind"]
        for row in _stream_rows(
            _agent(workspace) / "goals" / _GOAL / "ledger.jsonl"
        )
    ]
    assert kinds.count("run_recorded") == 2, kinds
    # The run judged the verdict again: its own stream holds a failed
    # receipt of it, and no decision.
    rejudged = _failed_receipts(_run_stream(workspace, 2))
    assert rejudged
    assert not [
        row
        for row in _stream_rows(_run_stream(workspace, 2))
        if row["kind"] == "scientific_decision_recorded"
    ]
    certificate = _completions(_run_stream(workspace, 2))[-1]
    assert certificate["status"] == "passed", certificate
    assert any(
        item.startswith(f"failed_criterion:{_RULE}:answered:")
        for item in certificate["anomaly_output_ids"]
    )
    assert result.settlement == "achieved_with_observations", result.reasons
    assert any(
        f"failed_criterion:{_RULE}:answered" in r for r in result.reasons
    )


def _plan_calls(artifact_id):
    """The o2r session's plan, as the (tool, arguments) pairs it called:
    extraction, the stability criterion, the claims, and the plan itself,
    which the host of a woken session walks the moment it is planned."""

    turn = _o2r_turns(artifact_id, cite_verdict=False)[0]({})
    return [
        (call["function"]["name"], json.loads(call["function"]["arguments"]))
        for call in turn["choices"][0]["message"]["tool_calls"]
    ]


def _live_woken_session(monkeypatch, tmp_path, *, cite):
    """Cycle 2's session, run by run_live_agent_session over the goal's own
    wake: it records a decision citing the failed receipt cycle 1's run
    minted (``cite``), then plans the o2r analysis chain, whose walk judges
    the criterion again."""

    import chemsmart.agent.runtime.alibaba as alibaba
    from chemsmart.agent.live_session import run_live_agent_session

    agent, keys = _fenced_host(monkeypatch, tmp_path)
    artifact_id = f"pyscf-result-{file_sha256(_RESULT)[:16]}"
    run_stream = _run_stream(tmp_path / "ws", 1)

    def decided(_model):
        return [
            (
                "record_scientific_decision",
                {
                    "decision_id": "o2-rks-unstable",
                    "assumptions": [
                        "the restricted reference at the geometry"
                    ],
                    "method_rationale": "the task fixed B3LYP/def2-SVP",
                    "alternatives": ["a broken-symmetry solution, not asked"],
                    "uncertainties": ["SCF convergence precision"],
                    "diagnostics": ["the external eigenvalue is negative"],
                    "stage_order": ["validate", "decide"],
                    "evidence_refs": [],
                    "postprocessing_receipt_sha256s": _failed_receipts(
                        run_stream
                    )[-1:],
                },
            )
        ]

    steps = [
        *((decided,) if cite else ()),
        lambda _model: _plan_calls(artifact_id),
        lambda _model: "The energy and the eigenvalue are delivered.",
    ]
    script = _ScriptedModel(steps)
    monkeypatch.setattr(
        alibaba, "AlibabaTokenPlanHttpsTransport", script.transport()
    )

    def step(workspace, kwargs):
        # The result the run read, where a workspace keeps its results.
        shutil.copytree(_RESULT.parent, workspace / "results" / "o2")
        return run_live_agent_session(
            **{
                **kwargs,
                "provider": "alibaba-token-plan",
                "provider_config_file": agent,
                "secret_file": keys,
                "exposure_mode": "eager",
            }
        )

    return step, script


def _session_goal(monkeypatch, tmp_path, *, cite):
    (tmp_path / "first").mkdir()
    woken, script = _live_woken_session(monkeypatch, tmp_path, cite=cite)
    sessions = iter(
        [
            _woken(
                "live-20260928T000000000000Z-truth-plan",
                cite=False,
                reads=False,
            ),
            woken,
        ]
    )
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    executes = iter([_run_with_the_failed_criterion(tmp_path / "first")])

    def plan_session(**kwargs):
        return next(sessions)(workspace, kwargs)

    def resolve_review(**_kwargs):
        return ("d" * 64, _bundle_file(tmp_path))

    def execute_bundle(*, approval_file, workspace, run_directory):
        return next(executes)(run_directory)

    result = run_goal_loop(
        task="Is the restricted reference of singlet O2 stable?",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id=_GOAL,
        granted_by="claude-researcher-truth-owner-delegated",
        max_revisions=1,
        plan_session=plan_session,
        resolve_review=resolve_review,
        execute_bundle=execute_bundle,
    )
    return result, script


@pytest.mark.parametrize("cite", [True, False])
def test_a_session_walk_reads_the_verdict_its_goal_answered(
    monkeypatch, tmp_path, cite
):
    """The woken session cites the run's failed receipt, then its own
    analysis-only walk judges the criterion again. Answered in the goal,
    the walk's certificate certifies the delivery and the word carries the
    verdict; with no citing decision, nothing standing on it is certified
    and the word names the verdict."""

    result, _script = _session_goal(monkeypatch, tmp_path, cite=cite)

    workspace = tmp_path / "ws"
    ledger = _stream_rows(_agent(workspace) / "goals" / _GOAL / "ledger.jsonl")
    (woken,) = [
        row["payload"]["run_id"]
        for row in ledger
        if row["kind"] == "session_stream_recorded"
    ][1:]
    stream = _agent(workspace) / "runs" / woken / "events.jsonl"
    # The walk judged the verdict again in the woken session's own stream.
    assert _failed_receipts(stream)
    decisions = [
        row
        for row in _stream_rows(stream)
        if row["kind"] == "scientific_decision_recorded"
    ]
    assert bool(decisions) is cite
    certificate = _completions(stream)[-1]
    reasons = " ".join(result.reasons)
    if cite:
        assert certificate["status"] == "passed", certificate
        assert result.settlement == "achieved_with_observations", reasons
        assert f"failed_criterion:{_RULE}:answered" in reasons
    else:
        assert certificate["status"] == "partial", certificate
        assert result.settlement == "returned_to_human", reasons
        assert _RULE in reasons, reasons
