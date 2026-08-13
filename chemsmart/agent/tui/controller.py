"""Provider-neutral controller shared by the Textual presentation layer."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TYPE_CHECKING

from chemsmart.agent._contracts import ContractError, file_sha256

if TYPE_CHECKING:
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1

from .approval import (
    ReviewedExecutionApprovalV1,
    ReviewedExecutionReviewV1,
    review_execution_approval,
    review_execution_review,
)


class AgentTuiPhase(str, Enum):
    READY = "ready"
    PLANNING = "planning"
    PREVIEW_READY = "preview-ready"
    REQUEST_REVIEWED = "request-reviewed"
    APPROVAL_BOUND = "approval-bound"
    EXECUTING = "executing"
    COMPLETE = "complete"
    BLOCKED = "blocked"
    ERROR = "error"


@dataclass(frozen=True)
class AgentSessionConfigV1:
    workspace: Path
    secret_file: Path
    provider: str | None = None
    provider_config_file: Path | None = None
    execution_envelope_file: Path | None = None
    review_file: Path | None = None
    identity_manifest: Path | None = None
    analysis_completion_file: Path | None = None

    def __post_init__(self) -> None:
        workspace = self.workspace.expanduser().resolve()
        secret = self.secret_file.expanduser().resolve()
        if not workspace.is_dir() or workspace.is_symlink():
            raise ContractError("TUI workspace must be a current directory")
        if not secret.is_file() or secret.is_symlink():
            raise ContractError(
                "TUI secret file must be a current regular file"
            )
        review = self.review_file
        envelope = self.execution_envelope_file
        if (envelope is None) != (review is None):
            raise ContractError(
                "TUI execution review requires both an execution envelope "
                "and an absolute review output path"
            )
        if envelope is not None:
            envelope = envelope.expanduser().resolve()
            if not envelope.is_file() or envelope.is_symlink():
                raise ContractError(
                    "TUI execution envelope must be a current regular file"
                )
            review = review.expanduser()
            if not review.is_absolute():
                raise ContractError("TUI review output file must be absolute")
        object.__setattr__(self, "workspace", workspace)
        object.__setattr__(self, "secret_file", secret)
        object.__setattr__(self, "execution_envelope_file", envelope)
        object.__setattr__(self, "review_file", review)


class AgentTuiController:
    """Human-driven state machine; provider turns cannot grant approval."""

    def __init__(self, config: AgentSessionConfigV1) -> None:
        self.config = config
        self.phase = AgentTuiPhase.READY
        self.task = ""
        self.plan_result: LiveAgentSessionResultV1 | None = None
        self.reviewed_review: ReviewedExecutionReviewV1 | None = None
        self.reviewed_approval: ReviewedExecutionApprovalV1 | None = None
        self.execution_result: WorkflowExecutionResultV1 | None = None

    def plan(self, task: str) -> LiveAgentSessionResultV1:
        from chemsmart.agent.identity import (
            load_approved_molecular_input_manifest,
        )
        from chemsmart.agent.live_session import run_live_agent_session

        normalized = str(task).strip()
        if not normalized:
            raise ContractError("agent task must not be empty")
        self.phase = AgentTuiPhase.PLANNING
        approved_inputs = (
            load_approved_molecular_input_manifest(
                self.config.identity_manifest,
                workspace=self.config.workspace,
            )
            if self.config.identity_manifest is not None
            else ()
        )
        try:
            result = run_live_agent_session(
                task=normalized,
                provider=(
                    self.config.provider.lower()
                    if self.config.provider
                    else None
                ),
                provider_config_file=self.config.provider_config_file,
                secret_file=self.config.secret_file,
                workspace=self.config.workspace,
                execution_enabled=False,
                approval_file=None,
                execution_envelope_file=self.config.execution_envelope_file,
                analysis_completion_file=(
                    self.config.analysis_completion_file
                ),
                approved_molecular_inputs=approved_inputs,
                review_file=self.config.review_file,
            )
        except Exception:
            self.phase = AgentTuiPhase.ERROR
            raise
        self.task = normalized
        self.plan_result = result
        self.reviewed_review = None
        self.reviewed_approval = None
        self.execution_result = None
        if result.execution_review:
            self.phase = AgentTuiPhase.PREVIEW_READY
        elif result.terminal_state in {"complete", "planned"}:
            self.phase = AgentTuiPhase.COMPLETE
        else:
            self.phase = AgentTuiPhase.BLOCKED
        return result

    def resolve_review(
        self,
        *,
        decision: str,
        actor: str,
        output_file: str | Path | None = None,
    ):
        """Record one explicit human decision and optionally bind its bundle."""

        if self.reviewed_review is None:
            raise ContractError("review the execution packet first")
        if self.reviewed_approval is not None:
            raise ContractError(
                "this review already produced an approval bundle; it cannot "
                "be revoked or replaced in place"
            )
        from chemsmart.agent.live_session import (
            resolve_workflow_execution_review,
        )

        resolution, bundle = resolve_workflow_execution_review(
            review_file=self.reviewed_review.path,
            reviewed_sha256=self.reviewed_review.review.review_sha256,
            decision=decision,
            actor=actor,
            output_file=(
                Path(output_file).expanduser().resolve()
                if output_file is not None
                else None
            ),
        )
        if bundle is None:
            self.reviewed_review = None
            self.reviewed_approval = None
            self.phase = AgentTuiPhase.PREVIEW_READY
        else:
            self.reviewed_approval = review_execution_approval(
                Path(output_file).expanduser().resolve(),
                reviewed_review=self.reviewed_review,
            )
            self.phase = AgentTuiPhase.APPROVAL_BOUND
        return resolution, self.reviewed_approval

    def review_request(self, path: str | Path) -> ReviewedExecutionReviewV1:
        if self.plan_result is None or self.phase not in {
            AgentTuiPhase.PREVIEW_READY,
            AgentTuiPhase.REQUEST_REVIEWED,
            AgentTuiPhase.APPROVAL_BOUND,
        }:
            raise ContractError(
                "finish a safe preview before reviewing approval"
            )
        reviewed = review_execution_review(
            path,
            task_spec_sha256=self.plan_result.task_spec_sha256,
            workspace=self.config.workspace,
        )
        self.reviewed_review = reviewed
        self.reviewed_approval = None
        self.phase = AgentTuiPhase.REQUEST_REVIEWED
        return reviewed

    def bind_approval(self, path: str | Path) -> ReviewedExecutionApprovalV1:
        if self.reviewed_review is None:
            raise ContractError("review an execution review packet first")
        reviewed = review_execution_approval(
            path, reviewed_review=self.reviewed_review
        )
        self.reviewed_approval = reviewed
        self.phase = AgentTuiPhase.APPROVAL_BOUND
        return reviewed

    def deny(self) -> None:
        """Forget a locally loaded chain without recording a human decision."""

        self.reviewed_review = None
        self.reviewed_approval = None
        self.phase = (
            AgentTuiPhase.PREVIEW_READY
            if self.plan_result is not None
            else AgentTuiPhase.READY
        )

    def execute(
        self,
        *,
        confirmation_digest: str,
        run_directory: str | Path,
    ) -> WorkflowExecutionResultV1:
        from chemsmart.agent.executor import execute_approved_workflow

        reviewed = self.reviewed_approval
        if reviewed is None or self.plan_result is None:
            raise ContractError("bind a reviewed execution approval first")
        if reviewed.path.is_symlink() or not reviewed.path.is_file():
            raise ContractError(
                "reviewed execution approval is no longer a regular file"
            )
        current_artifact_sha256 = file_sha256(reviewed.path)
        if current_artifact_sha256 != reviewed.artifact_sha256:
            raise ContractError(
                "execution approval bytes changed after human review; bind "
                "and confirm the current file again"
            )
        normalized = str(confirmation_digest).strip().lower()
        if normalized != current_artifact_sha256:
            raise ContractError(
                "execution confirmation must equal the complete reviewed "
                "approval file SHA-256"
            )
        target = Path(run_directory).expanduser()
        if target.exists() and target.is_symlink():
            raise ContractError("execution run directory cannot be a symlink")
        target = target.resolve()
        self.phase = AgentTuiPhase.EXECUTING
        try:
            result = execute_approved_workflow(
                approval_file=reviewed.path,
                workspace=self.config.workspace,
                run_directory=target,
                task_spec_sha256=self.plan_result.task_spec_sha256,
                expected_approval_file_sha256=current_artifact_sha256,
            )
        except Exception:
            self.phase = AgentTuiPhase.ERROR
            raise
        self.execution_result = result
        self.phase = (
            AgentTuiPhase.COMPLETE
            if result.status == "completed"
            else AgentTuiPhase.BLOCKED
        )
        return result


__all__ = [
    "AgentSessionConfigV1",
    "AgentTuiController",
    "AgentTuiPhase",
]
