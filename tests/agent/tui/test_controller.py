"""Focused gates for the interactive approval chain."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tui.approval import (
    ReviewedExecutionApprovalV1,
    ReviewedExecutionReviewV1,
)
from chemsmart.agent.tui.controller import (
    AgentSessionConfigV1,
    AgentTuiController,
    AgentTuiPhase,
)


def _controller(tmp_path: Path) -> AgentTuiController:
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    envelope = tmp_path / "execution-envelope.json"
    envelope.write_text("{}\n", encoding="utf-8")
    controller = AgentTuiController(
        AgentSessionConfigV1(
            workspace=tmp_path,
            secret_file=secret,
            execution_envelope_file=envelope,
            review_file=tmp_path / "review-output.json",
        )
    )
    controller.task = "single point"
    controller.plan_result = SimpleNamespace(task_spec_sha256="a" * 64)
    controller.phase = AgentTuiPhase.PREVIEW_READY
    return controller


def test_approval_cannot_be_bound_before_review(tmp_path):
    controller = _controller(tmp_path)
    with pytest.raises(ContractError, match="review an execution review"):
        controller.bind_approval(tmp_path / "approval.json")


def test_execution_requires_full_approval_file_digest(tmp_path, monkeypatch):
    controller = _controller(tmp_path)
    review_path = tmp_path / "review.json"
    review_path.write_text("{}\n", encoding="utf-8")
    approval_path = tmp_path / "approval.json"
    approval_path.write_text("{}\n", encoding="utf-8")
    controller.reviewed_review = ReviewedExecutionReviewV1(
        path=review_path,
        artifact_sha256="b" * 64,
        review=SimpleNamespace(review_sha256="c" * 64),
        canonical_text="{}",
    )
    controller.reviewed_approval = ReviewedExecutionApprovalV1(
        path=approval_path,
        artifact_sha256="d" * 64,
        bundle=SimpleNamespace(),
        canonical_text="{}",
    )
    controller.phase = AgentTuiPhase.APPROVAL_BOUND

    monkeypatch.setattr(
        "chemsmart.agent.executor.execute_approved_workflow",
        lambda **_kwargs: pytest.fail("executor must not be reached"),
    )
    with pytest.raises(ContractError, match="complete reviewed"):
        controller.execute(
            confirmation_digest="d" * 12,
            run_directory=tmp_path / "run",
        )


def test_deny_forgets_every_grant_but_keeps_preview(tmp_path):
    controller = _controller(tmp_path)
    controller.reviewed_review = object()
    controller.reviewed_approval = object()

    controller.deny()

    assert controller.reviewed_review is None
    assert controller.reviewed_approval is None
    assert controller.phase is AgentTuiPhase.PREVIEW_READY


def test_explicit_human_approve_uses_review_digest_and_writes_bundle(
    tmp_path, monkeypatch
):
    controller = _controller(tmp_path)
    review_path = tmp_path / "review.json"
    review_path.write_text("{}\n", encoding="utf-8")
    controller.reviewed_review = ReviewedExecutionReviewV1(
        path=review_path,
        artifact_sha256="b" * 64,
        review=SimpleNamespace(review_sha256="c" * 64),
        canonical_text="{}",
    )
    captured = {}

    def resolve(**kwargs):
        captured.update(kwargs)
        return SimpleNamespace(decision="approve"), SimpleNamespace()

    monkeypatch.setattr(
        "chemsmart.agent.live_session.resolve_workflow_execution_review",
        resolve,
    )
    reviewed = ReviewedExecutionApprovalV1(
        path=tmp_path / "approval.json",
        artifact_sha256="d" * 64,
        bundle=SimpleNamespace(),
        canonical_text="{}",
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.controller.review_execution_approval",
        lambda *_args, **_kwargs: reviewed,
    )

    _resolution, observed = controller.resolve_review(
        decision="approve",
        actor="human-reviewer",
        output_file=tmp_path / "approval.json",
    )

    assert observed is reviewed
    assert captured["reviewed_sha256"] == "c" * 64
    assert captured["decision"] == "approve"
    assert captured["actor"] == "human-reviewer"
    assert captured["output_file"] == (tmp_path / "approval.json").resolve()
    assert controller.phase is AgentTuiPhase.APPROVAL_BOUND


def test_deny_resolution_requires_an_explicit_rereview(tmp_path, monkeypatch):
    controller = _controller(tmp_path)
    review_path = tmp_path / "review.json"
    review_path.write_text("{}\n", encoding="utf-8")
    controller.reviewed_review = ReviewedExecutionReviewV1(
        path=review_path,
        artifact_sha256="b" * 64,
        review=SimpleNamespace(review_sha256="c" * 64),
        canonical_text="{}",
    )
    monkeypatch.setattr(
        "chemsmart.agent.live_session.resolve_workflow_execution_review",
        lambda **_kwargs: (SimpleNamespace(decision="deny"), None),
    )

    controller.resolve_review(decision="deny", actor="human-reviewer")

    assert controller.reviewed_review is None
    assert controller.reviewed_approval is None
    assert controller.phase is AgentTuiPhase.PREVIEW_READY
