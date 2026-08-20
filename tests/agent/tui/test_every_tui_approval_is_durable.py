"""One /approve, durable everywhere -- the RESTART test.

Before this round, the TUI's one-shot fences were RAM: a process restart
could re-execute an already-approved review with nothing refusing. Every
TUI approval now records the same evidence the file pipeline records --
decision log, one-shot bundle, workspace consumption ledger -- and a second
approve of the same review in the same workspace is refused across
restarts, first by the UI-thread pre-check and independently by the
deterministic resolution identity inside the shared per-review decision
log. A deliberate second decision remains `chemsmart agent review`.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("textual")

from chemsmart.agent._contracts import ContractError  # noqa: E402
from chemsmart.agent.live_session import (  # noqa: E402
    claim_workflow_execution_approval_bundle,
    load_workflow_execution_approval_bundle,
)
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
    AgentTuiPhase,
)
from tests.agent.test_exact_execution_approval_chain import (  # noqa: E402
    _review,
)


def _controller(workspace: Path) -> AgentTuiController:
    secret = workspace / "secret.env"
    if not secret.exists():
        secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    return AgentTuiController(
        AgentSessionConfigV1(workspace=workspace, secret_file=secret)
    )


def _claiming_stub(monkeypatch, workspace: Path, record: dict):
    """Stand in for the engine walk but keep the durable claim real."""

    from chemsmart.agent import executor

    def fake_execute(*, approval_file, workspace, run_directory, **_kw):
        bundle = load_workflow_execution_approval_bundle(Path(approval_file))
        claim_workflow_execution_approval_bundle(
            bundle, workspace=Path(workspace)
        )
        record["executed"] = record.get("executed", 0) + 1
        return SimpleNamespace(
            workflow_id=bundle.approved_scientific_plan.workflow_id,
            status="completed",
            nodes=(),
            analysis_nodes=(),
            analysis_status="",
            analysis_report_path="",
            analysis_completion_receipt_sha256s=(),
            non_executable_node_ids=(),
            provider_calls=0,
            run_directory=str(run_directory),
        )

    monkeypatch.setattr(executor, "execute_approved_workflow", fake_execute)


def test_a_tui_approval_leaves_durable_evidence_and_refuses_a_restart_replay(
    tmp_path, monkeypatch
):
    review = _review(tmp_path)
    workspace = tmp_path / "workspace"  # the builder's bound workspace
    record: dict = {}
    _claiming_stub(monkeypatch, workspace, record)

    first = _controller(workspace)
    first.restore_prepared_execution(review)
    assert first.phase is AgentTuiPhase.REQUEST_REVIEWED
    first.begin_execution()
    result = first.execute_begun()
    assert result.status == "completed"
    assert record["executed"] == 1

    sha16 = review.review_sha256[:16]
    scope = workspace / ".chemsmart-agent" / "decisions" / sha16
    assert (scope / "decisions.jsonl").exists()
    assert (scope / "bundle.json").exists()
    assert (
        workspace / ".chemsmart-agent" / "reviews" / f"{sha16}.json"
    ).exists()
    ledgers = list(
        (workspace / ".chemsmart-agent" / "approvals").glob(
            "tui-*/consumption-events.jsonl"
        )
    )
    assert ledgers, "the workspace consumption ledger must be durable"

    # A fresh process: the RAM fences are gone; the durable ones are not.
    second = _controller(workspace)
    second.restore_prepared_execution(review)
    with pytest.raises(ContractError, match="already approved"):
        second.begin_execution()

    # Even bypassing the pre-check, the shared decision log refuses the
    # second resolution before any bundle could be written.
    second._begun_execution = review
    second.execution_id = "tui-forced-second"
    second.execution_run_directory = (
        workspace / ".chemsmart-agent" / "executions" / "tui-forced-second"
    )
    with pytest.raises(ContractError):
        second.execute_begun()
    assert record["executed"] == 1, "no second engine walk happened"


def test_planning_exports_the_workspace_review_copy(tmp_path, monkeypatch):
    review = _review(tmp_path)
    workspace = tmp_path / "workspace"  # the builder's bound workspace

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(
        live_session,
        "run_live_agent_session",
        lambda **_kw: SimpleNamespace(
            prepared_execution=review,
            terminal_state="waiting_for_approval",
            public_transcript=(),
        ),
    )
    controller = _controller(workspace)
    controller.run_planning("bind the water input and preview a single point")

    copy = (
        workspace
        / ".chemsmart-agent"
        / "reviews"
        / f"{review.review_sha256[:16]}.json"
    )
    assert copy.exists()
    payload = json.loads(copy.read_text(encoding="utf-8"))
    assert (
        payload["workflow_execution_review"]["review_sha256"]
        == review.review_sha256
    )
    assert controller.review_copy_note == ""

    # Idempotent on a re-plan of the identical review.
    controller.run_planning("bind the water input and preview a single point")
    assert controller.review_copy_note == ""
