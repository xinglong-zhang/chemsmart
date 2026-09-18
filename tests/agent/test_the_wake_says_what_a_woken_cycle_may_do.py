"""A woken cycle is told what it may do, what binds it, and how to read.

NOVEL-3 ino3 (2026-09-05): the woken cycle's system prompt said
"Execution is not exposed" while the wake context beneath it granted
execution under the goal, and two cycles ended "planned"; four
remaining budget numbers with no grant beside them read as plenty
while the engine wall was the binding line; and inspect_run refused
the run_id the wake's own outcome record had just named.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent.live_session import _system_prompt
from chemsmart.agent.rules import rules_by_id
from tests.agent.test_the_goal_loop_recovers_or_returns import (
    _engine_stream,
    _execute,
    _loop,
    _planning_session,
    _review_payload,
)

_AUTHORITY = rules_by_id()["wake.goal_authority"].text


@pytest.mark.capability("rule:wake.goal_authority")
def test_the_prompt_says_the_goal_authority_at_wake_and_not_before():
    cold = _system_prompt({})
    assert "Execution is not exposed" in cold
    assert _AUTHORITY not in cold
    woken = _system_prompt({}, goal_record={"goal_id": "goal-t1"})
    assert "Execution is not exposed" not in woken
    assert woken.count(_AUTHORITY) == 1


#: Every way the prompt has said, or could say, that the provider does
#: not get execution. The test above pinned one exact string, so a
#: paraphrase assembled somewhere else survived it.
_DENIALS = (
    "not exposed",
    "do not expose engine execution",
    "does not expose engine execution",
    "no engine execution is exposed",
    "inert exact workflow",
)


@pytest.mark.capability("rule:wake.goal_authority")
@pytest.mark.parametrize("bounded", (True, False))
def test_a_goal_prompt_never_denies_the_execution_it_grants(bounded):
    """One sentence may speak to execution authority, and under a goal it
    is the grant.

    OPEN-2's ino3-qwen (2026-09-07) spent zero of forty engine calls and
    gave "engine execution was not exposed to this session" as the
    reason. The removed sentence was gone; the bounded-review sentence
    still called the review inert and said the operating bounds "do not
    expose engine execution or human approval to the provider", and it
    sat in that run's prompt beside the grant. Assert the family, not
    one string.
    """

    woken = _system_prompt(
        {},
        bounded_review_requested=bounded,
        goal_record={"goal_id": "goal-t1"},
    )
    assert woken.count(_AUTHORITY) == 1
    remainder = woken.replace(_AUTHORITY, "")
    for phrase in _DENIALS:
        assert phrase not in remainder, phrase

    # Without a goal there is nothing to contradict, and the session is
    # still told plainly that it launches nothing.
    cold = _system_prompt({}, bounded_review_requested=bounded)
    assert "Execution is not exposed" in cold


@pytest.mark.capability("rule:stem.preview_only_is_review_authority")
def test_bounded_review_affordance_reaches_the_model_prompt():
    """The registered review-authority rule reaches the planning surface."""

    rule = rules_by_id()["stem.preview_only_is_review_authority"]
    prompt = _system_prompt(
        {},
        bounded_review_requested=True,
        goal_record={"goal_id": "goal-t1"},
    )
    assert rule.placement == "stem"
    assert prompt.count(rule.text) == 1


@pytest.mark.capability("rule:wake.failed_validation_receipt_answers_verdict")
def test_failed_validation_citation_affordance_reaches_each_goal_context(
    tmp_path,
):
    """A failed typed verdict must name its receipt at both goal entries.

    CUHK acetamide r8 delivered an honest saddle but cited only its
    characterisation and completion receipts.  The settlement correctly
    returned it to the human; this witness reaches the actual initial and
    recovery contexts rather than testing a detached policy string.
    """

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
    rule = rules_by_id()["wake.failed_validation_receipt_answers_verdict"]
    assert rule.placement == "wake"
    assert len(contexts) == 2
    assert all(
        context["authority"].count(rule.text) == 1 for context in contexts
    )


def test_the_budget_block_leads_with_the_binding_line(tmp_path):
    """Production path: the second cycle's wake context, after one run
    whose slowest node took five seconds."""

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
    assert first["budgets"]["binding_line"] == (
        "nearest exhausted: engine calls 6 of 6 remaining (100%)"
    )
    assert "slowest node" not in first["budgets"]["binding_line"]
    line = second["budgets"]["binding_line"]
    # The revision is charged when admitted, after this wake is built,
    # so engine calls at 83% are the binding line here.
    assert line.startswith(
        "nearest exhausted: engine calls 5 of 6 remaining (83%)"
    )
    assert "the slowest node this goal ran took 5 s" in line
    assert "fits 1439 such nodes" in line
    assert second["authority"].startswith(_AUTHORITY)
    assert "by its run_id" in second["authority"]


def test_inspect_run_resolves_the_outcomes_run_id(tmp_path):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    host = _host(tmp_path)
    evidence = tmp_path / "evidence"
    target = evidence / ".chemsmart-agent" / "goals" / "goal-t1" / "runs"
    _engine_stream(tmp_path, target / "cycle-1", failed=False)
    host.run_evidence_root = Path(evidence)

    by_reference = host._inspect_run_outcome(
        "t1", {"run": "goals/goal-t1/runs/cycle-1"}
    )
    assert "resolved_from" not in by_reference
    by_id = host._inspect_run_outcome("t1", {"run": "cycle-1"})
    assert by_id["run"] == "goals/goal-t1/runs/cycle-1"
    assert by_id["resolved_from"] == "cycle-1"
    assert by_id["workflow_state"] == by_reference["workflow_state"]

    with pytest.raises(Exception, match="records no run 'cycle-9'"):
        host._inspect_run_outcome("t1", {"run": "cycle-9"})
