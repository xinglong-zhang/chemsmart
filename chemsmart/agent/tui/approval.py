"""Read-only TUI projections of canonical approval-chain artifacts."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
)
from chemsmart.agent.execution import (
    WorkflowExecutionApprovalBundleV1,
    WorkflowExecutionReviewV1,
)


@dataclass(frozen=True)
class ReviewedExecutionReviewV1:
    path: Path
    artifact_sha256: str
    review: WorkflowExecutionReviewV1
    canonical_text: str


@dataclass(frozen=True)
class ReviewedExecutionApprovalV1:
    path: Path
    artifact_sha256: str
    bundle: WorkflowExecutionApprovalBundleV1
    canonical_text: str


def _artifact_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def review_execution_review(
    path: str | Path,
    *,
    task_spec_sha256: str,
    workspace: str | Path,
) -> ReviewedExecutionReviewV1:
    """Load the host-produced inert review and bind it to this TUI task."""

    from chemsmart.agent.live_session import load_workflow_execution_review

    candidate = Path(path).expanduser()
    if candidate.is_symlink() or not candidate.is_file():
        raise ContractError("execution review must be a current regular file")
    candidate = candidate.resolve()
    review = load_workflow_execution_review(candidate)
    if review.scientific_plan.task_spec_sha256 != task_spec_sha256:
        raise ContractError("execution review targets another task")
    if Path(review.request.workspace).resolve() != Path(workspace).resolve():
        raise ContractError("execution review targets another workspace")
    return ReviewedExecutionReviewV1(
        path=candidate,
        artifact_sha256=_artifact_sha256(candidate),
        review=review,
        canonical_text=canonical_json(
            {"workflow_execution_review": canonical_data(review)}
        ),
    )


def review_execution_approval(
    path: str | Path,
    *,
    reviewed_review: ReviewedExecutionReviewV1,
) -> ReviewedExecutionApprovalV1:
    """Load a one-shot grant and prove it binds the reviewed packet."""

    from chemsmart.agent.live_session import (
        load_workflow_execution_approval_bundle,
    )

    candidate = Path(path).expanduser()
    if candidate.is_symlink() or not candidate.is_file():
        raise ContractError(
            "execution approval must be a current regular file"
        )
    candidate = candidate.resolve()
    artifact_sha256 = _artifact_sha256(candidate)
    bundle = load_workflow_execution_approval_bundle(
        candidate,
        expected_file_sha256=artifact_sha256,
    )
    if bundle.review_sha256 != reviewed_review.review.review_sha256:
        raise ContractError("approval bundle differs from the reviewed packet")
    if bundle.resolution.decision != "approve" or not bundle.one_shot:
        raise ContractError("approval bundle does not carry a one-shot grant")
    return ReviewedExecutionApprovalV1(
        path=candidate,
        artifact_sha256=artifact_sha256,
        bundle=bundle,
        canonical_text=canonical_json(
            {"workflow_execution_approval_bundle": canonical_data(bundle)}
        ),
    )


__all__ = [
    "ReviewedExecutionApprovalV1",
    "ReviewedExecutionReviewV1",
    "review_execution_approval",
    "review_execution_review",
]
