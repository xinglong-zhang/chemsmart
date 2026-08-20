"""Provider-neutral controller shared by the Textual presentation layer."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import getpass
from pathlib import Path
import threading
from typing import TYPE_CHECKING, Callable
import uuid

from chemsmart.agent._contracts import ContractError

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1
    from chemsmart.agent.executor import WorkflowExecutionResultV1
    from chemsmart.agent.live_session import LiveAgentSessionResultV1


class AgentTuiPhase(str, Enum):
    READY = "ready"
    PLANNING = "planning"
    PREVIEW_READY = "preview-ready"
    REQUEST_REVIEWED = "request-reviewed"
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
        if envelope is None and review is not None:
            raise ContractError(
                "a review export path requires an execution envelope"
            )
        if envelope is not None:
            envelope = envelope.expanduser().resolve()
            if not envelope.is_file() or envelope.is_symlink():
                raise ContractError(
                    "TUI execution envelope must be a current regular file"
                )
            if review is not None:
                review = review.expanduser()
                if not review.is_absolute():
                    raise ContractError(
                        "TUI review export file must be absolute"
                    )
        object.__setattr__(self, "workspace", workspace)
        object.__setattr__(self, "secret_file", secret)
        object.__setattr__(self, "execution_envelope_file", envelope)
        object.__setattr__(self, "review_file", review)


class AgentTuiController:
    """Human-driven state machine; provider turns cannot grant approval.

    The guards run on the UI thread (`begin_planning`, `begin_execution`) so
    a refused action refuses *before* any banner or worker starts; the heavy
    host calls run on a worker thread (`run_planning`, `execute_begun`).
    """

    def __init__(self, config: AgentSessionConfigV1) -> None:
        self.config = config
        self.phase = AgentTuiPhase.READY
        self.task = ""
        self.plan_result: LiveAgentSessionResultV1 | None = None
        self.prepared_execution: WorkflowExecutionReviewV1 | None = None
        self.execution_result: WorkflowExecutionResultV1 | None = None
        self.execution_id = ""
        self.execution_run_directory: Path | None = None
        #: Set by the view when the human double-escapes a planning session.
        self.cancel_planning = threading.Event()
        #: Set by the view to observe the planning run directory (event tail).
        self.on_run_directory: Callable[[Path], None] | None = None
        self._begun_execution: WorkflowExecutionReviewV1 | None = None

    # -- planning ----------------------------------------------------------

    def begin_planning(self, task: str) -> str:
        """UI-thread guard: validate and take the planning phase now."""

        normalized = str(task).strip()
        if not normalized:
            raise ContractError("agent task must not be empty")
        self.cancel_planning.clear()
        self.phase = AgentTuiPhase.PLANNING
        self.task = normalized
        return normalized

    def run_planning(self, task: str) -> LiveAgentSessionResultV1:
        """Worker-thread body: the live provider session."""

        from chemsmart.agent.identity import (
            load_approved_molecular_input_manifest,
        )
        from chemsmart.agent.live_session import run_live_agent_session

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
                task=task,
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
                on_run_directory=self.on_run_directory,
                should_stop=self.cancel_planning.is_set,
            )
        except Exception:
            self.phase = AgentTuiPhase.ERROR
            raise
        self.plan_result = result
        self.prepared_execution = result.prepared_execution
        self.execution_result = None
        if self.prepared_execution is not None:
            self.phase = AgentTuiPhase.REQUEST_REVIEWED
        elif result.terminal_state in {"complete", "planned"}:
            self.phase = AgentTuiPhase.COMPLETE
        elif result.terminal_state == "cancelled":
            self.phase = AgentTuiPhase.READY
        else:
            self.phase = AgentTuiPhase.BLOCKED
        return result

    # -- decision ----------------------------------------------------------

    def decline(self) -> None:
        """Decline the displayed workflow without creating run authority."""

        self.prepared_execution = None
        self.phase = (
            AgentTuiPhase.PREVIEW_READY
            if self.plan_result is not None
            else AgentTuiPhase.READY
        )

    def begin_execution(self) -> WorkflowExecutionReviewV1:
        """UI-thread guard: consume the pending authority before any launch.

        A failed run remains an observed run; rerunning it requires another
        explicit plan/review act.
        """

        prepared = self.prepared_execution
        if (
            prepared is None
            or self.phase is not AgentTuiPhase.REQUEST_REVIEWED
        ):
            raise ContractError(
                "finish planning and review the displayed ChemSmart "
                "workflow first"
            )
        self.execution_id = "tui-" + uuid.uuid4().hex
        self.execution_run_directory = (
            self.config.workspace
            / ".chemsmart-agent"
            / "executions"
            / self.execution_id
        )
        self.phase = AgentTuiPhase.EXECUTING
        self.prepared_execution = None
        self._begun_execution = prepared
        return prepared

    def execute_begun(self) -> WorkflowExecutionResultV1:
        """Worker-thread body: run the consumed review provider-free."""

        from chemsmart.agent.executor import execute_prepared_workflow

        prepared = self._begun_execution
        if prepared is None or self.execution_run_directory is None:
            raise ContractError("no begun execution to run")
        self._begun_execution = None
        try:
            result = execute_prepared_workflow(
                review=prepared,
                actor=getpass.getuser() or "local-user",
                execution_id=self.execution_id,
                workspace=self.config.workspace,
                run_directory=self.execution_run_directory,
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
