"""A goal's word reads every failed criterion the goal holds, on every path.

A failed acceptance criterion is a finding the goal can deliver once a
recorded decision cites a receipt that states it, in any cycle of the goal
(R10 Q19, ``goal.failed_criteria``). The settlement of a cycle that ends in
its session reads the verdicts and the decisions of every stream of the
goal; the settlement of a cycle that ends in a run read the run's own stream
-- which never holds a decision -- and carried an earlier cycle's rejections
forward as a set of results, whether or not a decision had since answered
them. R11 `truth` (Phase II item 1): the two paths answer one question and
must call one reading of it.

Driven through the goal loop with real hosts over the archived bytes of a
PySCF run of closed-shell singlet O2, whose restricted reference is unstable
to spin symmetry breaking: cycle 1's run fails the stability criterion and
claims the energy; cycle 2's woken session reads how that run ended, cites
the failed receipt or not, and approves a second run whose chain reads the
same result again -- re-claiming the energy, or re-judging it with the same
criterion, under the same two receipts. Cycle 2 is the goal's last
revision, so its run's settlement is final.

Not driven here: a second run that re-judges the verdict the woken session
answered. The settlement now reads that verdict as answered, and the word
still returns to the human, because the run's own completion receipt is
partial -- the executor's walk judges the claims standing on a failed
criterion with the decisions of its own host, which holds none. That is the
certificate's organ, not the settlement's, and it is reported rather than
repaired here.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from chemsmart.agent.driver import run_goal_loop

from .test_a_failed_criterion_is_a_finding_the_goal_can_deliver import (
    _RULE,
    _TASK,
    _engine_prefix,
    _run_turns,
    _run_with_the_failed_criterion,
    _stream_rows,
)
from .test_a_partial_delivery_ends_its_session import _WORKFLOW, _call, _turn
from .test_the_goal_loop_recovers_or_returns import (
    _bundle_file,
    _envelope_file,
    _review_payload,
)

pytestmark = pytest.mark.capability(
    "rule:wake.failed_validation_receipt_answers_verdict"
)

_GOAL = "goal-o2r-two-runs"


def _agent(workspace):
    return workspace / ".chemsmart-agent"


def _run_stream(workspace, cycle):
    return (
        _agent(workspace)
        / "goals"
        / _GOAL
        / "runs"
        / f"cycle-{cycle}"
        / "events.jsonl"
    )


def _planning(run_id):
    """Cycle 1's session: it plans the run and approves one wave."""

    return _woken(run_id, cite=False, reads=False)


def _woken(run_id, *, cite, reads=True):
    """A session that reads how cycle 1's run ended (``reads``), records a
    decision citing the run's failed validation receipt when ``cite`` --
    the route the wake prescribes -- and plans the next run."""

    def step(workspace, kwargs):
        stream = _agent(workspace) / "runs" / run_id / "events.jsonl"

        def turns(_artifact_id):
            acts = []
            if reads:
                acts.append(
                    lambda payload: _turn(
                        1,
                        "Reading how the run ended.",
                        (
                            _call(
                                1,
                                "inspect_run_outcome",
                                {"run": f"goals/{_GOAL}/runs/cycle-1"},
                            ),
                        ),
                    )
                )
            if cite:

                def decided(payload):
                    failed = [
                        row["payload"]["receipt_sha256"]
                        for row in _stream_rows(_run_stream(workspace, 1))
                        if row["kind"] == "scientific_validation_evaluated"
                        and not row["payload"]["all_rules_passed"]
                    ][-1:]
                    return _turn(
                        2,
                        "The reference is unstable; that is the finding.",
                        (
                            _call(
                                2,
                                "record_scientific_decision",
                                {
                                    "decision_id": "o2-rks-unstable",
                                    "assumptions": [
                                        "the restricted reference at the "
                                        "geometry"
                                    ],
                                    "method_rationale": "the task fixed "
                                    "B3LYP/def2-SVP",
                                    "alternatives": [
                                        "a broken-symmetry solution, not "
                                        "asked"
                                    ],
                                    "uncertainties": [
                                        "SCF convergence precision"
                                    ],
                                    "diagnostics": [
                                        "the external eigenvalue is negative"
                                    ],
                                    "stage_order": ["validate", "decide"],
                                    "evidence_refs": [],
                                    "postprocessing_receipt_sha256s": failed,
                                },
                            ),
                        ),
                    )

                acts.append(decided)
            acts.append(lambda payload: _turn(3, "The next run is planned."))
            return acts

        _run_turns(
            stream,
            turns,
            session_id="protocol-session",
            scratch=workspace.parent / f"scratch-{run_id}",
            workspace=workspace,
        )
        review_file = kwargs["review_file"]
        review_file.parent.mkdir(parents=True, exist_ok=True)
        review_file.write_text(json.dumps(_review_payload()), encoding="utf-8")
        from chemsmart.agent.cohort import build_execution_wave_decision

        return SimpleNamespace(
            terminal_state="waiting_for_approval",
            run_id=run_id,
            task_spec_sha256=_TASK,
            selected_execution_wave=("sp-initial",),
            execution_wave_decision=build_execution_wave_decision(
                state="selected",
                workflow_id="water-workflow",
                ready_node_ids=("sp-initial",),
                node_ids=("sp-initial",),
            ),
        )

    return step


def _energy_chain(artifact_id):
    """An approved chain that reads the energy from the result again and
    claims it, judging nothing."""

    extraction = {
        "workflow_id": _WORKFLOW,
        "stages": [
            {
                "node_id": "extract-o2-energy",
                "artifact_id": artifact_id,
                "dependencies": [],
                "inputs": [],
                "selectors": [
                    {"quantity_id": "ref-energy", "selector": "energy"}
                ],
                "outputs": [
                    {
                        "output_id": "ref-energy",
                        "quantity_kind": "energy",
                        "unit": "hartree",
                    }
                ],
                "support_state": "planned",
                "blocked_reason": "",
            }
        ],
    }
    claims = {
        "workflow_id": _WORKFLOW,
        "stages": [
            {
                "node_id": "claim-o2-energy",
                "dependencies": [],
                "inputs": [
                    {
                        "input_id": "ref-energy",
                        "source_kind": "analysis_output",
                        "producer_node_id": "extract-o2-energy",
                        "producer_output_id": "ref-energy",
                    }
                ],
                "outputs": [
                    {
                        "output_id": "ref-energy",
                        "quantity_kind": "energy",
                        "unit": "hartree",
                    }
                ],
                "support_state": "planned",
                "blocked_reason": "",
            }
        ],
    }
    finalise = {
        "plan_id": "o2-rks-energy-plan",
        "workflow_id": _WORKFLOW,
        "required_output_ids": ["ref-energy"],
    }
    return [
        lambda payload: _turn(
            1,
            "Reading the energy again.",
            (
                _call(1, "plan_result_extraction", extraction),
                _call(2, "plan_claim_rendering", claims),
                _call(3, "plan_scientific_workflow", finalise),
            ),
        ),
        lambda payload: _turn(2, "The energy is claimed."),
    ]


def _run_reclaiming_the_energy(tmp_path):
    """Cycle 2's run: an engine node, then a chain that reads the energy of
    the result cycle 1's criterion judged and claims it; a run's stream
    holds no decision."""

    def step(run_directory):
        _engine_prefix(tmp_path / "second", run_directory)
        _run_turns(
            run_directory / "events.jsonl",
            _energy_chain,
            session_id="protocol-session",
            scratch=tmp_path / "second",
        )
        return SimpleNamespace(status="completed", analysis_status="completed")

    return step


def _goal(tmp_path, *, cite, second_run):
    workspace = tmp_path / "ws"
    workspace.mkdir(parents=True, exist_ok=True)
    sessions = iter(
        [
            _planning("live-20260928T000000000000Z-truth-plan"),
            _woken("live-20260928T010000000000Z-truth-woken", cite=cite),
        ]
    )
    executes = iter(
        [
            _run_with_the_failed_criterion(tmp_path / "first"),
            second_run,
        ]
    )

    def plan_session(**kwargs):
        return next(sessions)(workspace, kwargs)

    def resolve_review(**_kwargs):
        return ("d" * 64, _bundle_file(tmp_path))

    def execute_bundle(*, approval_file, workspace, run_directory):
        return next(executes)(run_directory)

    return run_goal_loop(
        task="Is the restricted reference of singlet O2 stable?",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id=_GOAL,
        granted_by="claude-researcher-truth-owner-delegated",
        # Cycle 2 is the last revision, so its run's word is final.
        max_revisions=1,
        plan_session=plan_session,
        resolve_review=resolve_review,
        execute_bundle=execute_bundle,
    )


@pytest.mark.parametrize(
    "second, cite",
    [
        ("reclaims", True),
        ("reclaims", False),
        ("judges-again", False),
    ],
)
def test_a_run_word_reads_the_criteria_and_decisions_of_the_whole_goal(
    tmp_path, second, cite
):
    """Cycle 2's run reads the result cycle 1's criterion rejected: it
    claims its energy again (a number standing on that verdict), or judges
    it again with the same criterion (the same verdict). A decision of
    cycle 2's session that cites the failed receipt answers the verdict
    whichever stream typed it, and the word delivers the finding; without
    one, nothing standing on it is certified and the reason names the
    verdict."""

    (tmp_path / "first").mkdir()
    (tmp_path / "second").mkdir()
    second_run = (
        _run_with_the_failed_criterion(tmp_path / "second")
        if second == "judges-again"
        else _run_reclaiming_the_energy(tmp_path)
    )
    result = _goal(tmp_path, cite=cite, second_run=second_run)

    workspace = tmp_path / "ws"
    ledger = _stream_rows(_agent(workspace) / "goals" / _GOAL / "ledger.jsonl")
    kinds = [row["kind"] for row in ledger]
    # Both cycles ran: the first run's criterion held the goal open, the
    # woken session's plan was admitted, and the second run settled.
    assert kinds.count("run_recorded") == 2, kinds
    assert "revision_admitted" in kinds, kinds
    assert result.cycles == 2
    woken = _stream_rows(
        _agent(workspace)
        / "runs"
        / "live-20260928T010000000000Z-truth-woken"
        / "events.jsonl"
    )
    failed = {
        row["payload"]["receipt_sha256"]
        for row in _stream_rows(_run_stream(workspace, 1))
        if row["kind"] == "scientific_validation_evaluated"
        and not row["payload"]["all_rules_passed"]
    }
    cited = {
        reference.split(":")[-1]
        for row in woken
        if row["kind"] == "scientific_decision_recorded"
        for reference in row["payload"]["record"]["evidence_refs"]
    }
    # The woken decision stands in its stream exactly when it was asked
    # for, citing the receipt cycle 1's run minted.
    assert bool(cited & failed) is cite, (cited, failed)
    (settled,) = [row for row in ledger if row["kind"] == "goal_settled"]
    reasons = " ".join(settled["payload"]["reasons"])
    if cite:
        assert result.settlement == "achieved_with_observations", reasons
        assert f"failed_criterion:{_RULE}:answered" in reasons
        # The word cites the receipt the decision cites.
        assert (
            cited
            & failed
            & set(settled["payload"]["evidence"].get("receipt_sha256s") or ())
        )
    else:
        assert result.settlement == "returned_to_human", reasons
        assert _RULE in reasons
