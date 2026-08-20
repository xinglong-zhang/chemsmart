"""Resume reads only durable evidence and never resurrects spent authority."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("textual")

from chemsmart.agent.execution import (  # noqa: E402
    workflow_execution_review_json,
)
from chemsmart.agent.tui.resume import (  # noqa: E402
    find_latest_session,
    find_pending_review,
    has_history,
    read_session_story,
    run_rows,
)
from tests.agent.test_exact_execution_approval_chain import (  # noqa: E402
    _review,
)


def _line(kind: str, payload: dict, session_id: str = "s") -> str:
    return (
        json.dumps(
            {"kind": kind, "payload": payload, "session_id": session_id}
        )
        + "\n"
    )


def _planning_session(workspace: Path, name: str, task: str, ended: str):
    run_dir = workspace / ".chemsmart-agent" / "runs" / name
    run_dir.mkdir(parents=True)
    (run_dir / "events.jsonl").write_text(
        _line("session_started", {})
        + _line("runtime_terminated", {"terminal_state": ended, "reason": ""}),
        encoding="utf-8",
    )
    transcript = [
        {"role": "system", "content": "system prompt"},
        {"role": "user", "content": json.dumps({"task": task})},
        {"role": "assistant", "content": "The workflow is ready for review."},
    ]
    (run_dir / "public-transcript-abcd1234abcd1234.json").write_text(
        json.dumps({"transcript": transcript}), encoding="utf-8"
    )
    return run_dir


def test_the_story_is_read_from_durable_records(tmp_path: Path):
    assert not has_history(tmp_path)
    _planning_session(
        tmp_path, "live-20260820T010000000000Z-aa-1", "old task", "complete"
    )
    _planning_session(
        tmp_path,
        "live-20260820T020000000000Z-bb-2",
        "optimize the water dimer",
        "cancelled",
    )
    execution = tmp_path / ".chemsmart-agent" / "executions" / "tui-1"
    execution.mkdir(parents=True)
    (execution / "events.jsonl").write_text(
        json.dumps(
            {
                "kind": "workflow_execution_started",
                "payload": {},
                "session_id": "execute-water-sp",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    assert has_history(tmp_path)
    latest = find_latest_session(tmp_path)
    assert latest is not None
    assert latest.name.endswith("-bb-2"), "lexical UTC order, not mtime"
    story = read_session_story(latest)
    assert story is not None
    assert story.task == "optimize the water dimer"
    assert story.ended == "cancelled"
    assert story.final_prose == "The workflow is ready for review."

    rows = run_rows(tmp_path)
    assert rows[0].workflow_id == "water-sp"


def test_a_pending_review_is_found_and_a_decided_one_is_not(tmp_path: Path):
    review = _review(tmp_path)
    workspace = tmp_path / "workspace"
    sha16 = review.review_sha256[:16]
    reviews = workspace / ".chemsmart-agent" / "reviews"
    reviews.mkdir(parents=True)
    (reviews / f"{sha16}.json").write_text(
        workflow_execution_review_json(review), encoding="utf-8"
    )

    pending = find_pending_review(workspace)
    assert pending is not None and pending.status == "pending"
    assert pending.review_sha256 == review.review_sha256

    # An approve line with a bundle but no consumption: the honest recovery
    # is the exact provider-free launch, not a second decision.
    scope = workspace / ".chemsmart-agent" / "decisions" / sha16
    scope.mkdir(parents=True)
    (scope / "decisions.jsonl").write_text(
        json.dumps({"payload": {"decision": "approve"}}) + "\n",
        encoding="utf-8",
    )
    (scope / "bundle.json").write_text("{}", encoding="utf-8")

    approved = find_pending_review(workspace)
    assert approved is not None
    assert approved.status == "approved_unexecuted"
    assert "chemsmart agent run --approval-file" in approved.recovery

    # A durable consumption entry ends the pendingness entirely.
    ledger = (
        workspace
        / ".chemsmart-agent"
        / "approvals"
        / "tui-x"
        / "consumption-events.jsonl"
    )
    ledger.parent.mkdir(parents=True)
    ledger.write_text(
        json.dumps({"review_sha256": review.review_sha256}) + "\n",
        encoding="utf-8",
    )
    assert find_pending_review(workspace) is None
