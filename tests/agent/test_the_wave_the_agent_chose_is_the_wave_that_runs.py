"""The wave the Agent selected is the wave that reaches the scheduler.

A selection that lives only in a tool reply is a selection the array
never hears about: the model chooses in the planning session, and the
dispatcher runs in the driver afterwards, in another process for a woken
cycle. The selection therefore has to survive that gap, and the digest it
is bound by must be the one that matters -- the cohort manifest's, which
binds goal, cycle, membership and the approved bundle together, so an
element cannot resolve itself through a manifest written for a different
approval.

What it must *not* do is enter `LiveAgentSessionResultV1`'s digested
body. That digest covers every archived session result on disk, and
extending the body would change all of them -- the same mistake this
tree already paid for once, when contract v6's vocabulary grew after its
first artifacts were minted.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

from chemsmart.agent.live_session import LiveAgentSessionResultV1


def _result(**overrides):
    from chemsmart.agent._contracts import canonical_sha256

    body = {
        "schema_version": "chemsmart.live-agent-session-result.v1",
        "session_id": "s1",
        "task_spec_sha256": "a" * 64,
        "terminal_state": "planned",
        "execution_requested": True,
        "execution_profile_status": "ready",
        "final_text": "planned",
        "artifact_records": (),
        "conformance_records": (),
        "public_transcript": (),
        "successful_tool_calls": 1,
        "failed_tool_calls": 0,
        "execution_review": {},
        "event_stream_head_sha256": "",
    }
    return LiveAgentSessionResultV1(
        **body, result_sha256=canonical_sha256(body), **overrides
    )


def test_a_session_carries_the_wave_it_selected():
    result = _result(selected_execution_wave=("a2", "a1"))
    assert result.selected_execution_wave == ("a2", "a1")


def test_the_wave_does_not_change_an_archived_result_s_digest():
    """Every session result on disk must still verify."""

    without = _result()
    with_wave = _result(selected_execution_wave=("a1",))
    assert without.result_sha256 == with_wave.result_sha256, (
        "adding the wave to the digested body would invalidate every "
        "archived session result, which is a cost this field does not "
        "need to pay: the manifest is where membership is digest-bound"
    )


def test_the_driver_dispatches_the_wave_the_session_chose(tmp_path):
    """Drive the real driver: what the session chose is what is submitted."""

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import (
        _envelope_file,
        _planning_session,
        _review_payload,
    )

    workspace = tmp_path / "ws"
    workspace.mkdir()
    seen: dict = {}

    def dispatch_run(**kwargs):
        seen.update(kwargs)
        return SimpleNamespace(
            scheduler="SLURM",
            job_id="191",
            submitted_at="2026-09-01T00:00:00+00:00",
            submit_script=str(kwargs["run_directory"] / "sub.sh"),
        )

    def plan_session(**kw):
        session = _planning_session("live-1", review=_review_payload())(
            workspace, kw
        )
        object.__setattr__(
            session,
            "selected_execution_wave",
            ("conf-b-opt", "conf-a-opt", "conf-c-opt"),
        )
        return session

    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-w1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=plan_session,
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        dispatch_run=dispatch_run,
        dispatch="scheduler",
        server="canned-slurm",
    )
    assert driver.run().settlement == "parked"

    assert tuple(seen.get("cohort_node_ids") or ()) == (
        "conf-b-opt",
        "conf-a-opt",
        "conf-c-opt",
    ), (
        "the Agent's wave never reached the dispatcher, so the cohort "
        f"submitted was {seen.get('cohort_node_ids')!r}"
    )


def test_a_session_that_selected_no_wave_dispatches_as_one_job(tmp_path):
    """Nothing is invented: no selection is the single-job path."""

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import (
        _envelope_file,
        _planning_session,
        _review_payload,
    )

    workspace = tmp_path / "ws"
    workspace.mkdir()
    seen: dict = {}

    def dispatch_run(**kwargs):
        seen.update(kwargs)
        return SimpleNamespace(
            scheduler="SLURM",
            job_id="191",
            submitted_at="2026-09-01T00:00:00+00:00",
            submit_script=str(kwargs["run_directory"] / "sub.sh"),
        )

    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-w2",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: _planning_session(
            "live-1", review=_review_payload()
        )(workspace, kw),
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        dispatch_run=dispatch_run,
        dispatch="scheduler",
        server="canned-slurm",
    )
    assert driver.run().settlement == "parked"
    assert tuple(seen.get("cohort_node_ids") or ()) == ()


def test_the_ledger_records_which_wave_was_dispatched(tmp_path):
    """A reader of the goal must be able to say what ran together."""

    from chemsmart.agent.dispatch import DispatchReceiptV1

    receipt = DispatchReceiptV1(
        scheduler="SLURM",
        job_id="412",
        submitted_at="2026-09-16T00:00:00+00:00",
        submit_command="sbatch x.sh",
        submit_script="x.sh",
        run_directory="/tmp/run",
        approval_file="/tmp/bundle.json",
        goal_id="g1",
        cycle=1,
        wake_command="wake",
        wake_job_id="422",
        cohort_node_ids=("a1", "a2"),
    )
    record = json.loads(json.dumps(receipt.public_record()))
    assert record["cohort_node_ids"] == ["a1", "a2"]
    assert record["wake_job_id"] == "422"


def test_the_attribute_the_tool_sets_is_the_attribute_the_session_reads(
    tmp_path,
):
    """The hop a `getattr` default would have hidden.

    A renamed attribute read with a default returns an empty wave, which
    dispatches as a single job and is indistinguishable from a session
    that chose not to select one -- the declared-but-no-reader defect
    this round exists to remove. Both ends are driven here: a real host
    runs the real tool, and the expression the session builder uses reads
    it back.
    """

    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    # A real host declares it before any tool runs, so the session's read
    # cannot be an AttributeError on a session that selected nothing.
    assert host.selected_execution_wave == ()

    host._resolve_program_workflow = lambda workflow_id: SimpleNamespace(
        draft=SimpleNamespace(
            workflow_id="w1",
            nodes=(
                SimpleNamespace(node_id="a1", inputs=()),
                SimpleNamespace(node_id="a2", inputs=()),
            ),
        ),
        scientific_plan=SimpleNamespace(plan_sha256="d" * 64),
    )
    host._workflow_context = lambda draft, **_kw: SimpleNamespace(
        ready_node_ids=("a1", "a2"), waiting_node_ids=()
    )
    host._select_execution_wave(
        "t1", {"workflow_id": "w1", "node_ids": ["a1", "a2"]}
    )

    # Exactly what `run_live_agent_session` writes onto the result.
    assert tuple(host.selected_execution_wave or ()) == ("a1", "a2")
