"""Stable public contracts between providers, tools, and runtime policy."""

from __future__ import annotations

from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class RuntimeV2Mode(str, Enum):
    OFF = "off"
    SHADOW = "shadow"
    ACTIVE = "active"

    @classmethod
    def parse(cls, value: str | bool | None) -> "RuntimeV2Mode":
        if value is True:
            return cls.ACTIVE
        if value in {False, None, "", "0", "false", "off"}:
            return cls.OFF
        normalized = str(value).strip().lower()
        if normalized in {"1", "true", "on", "active"}:
            return cls.ACTIVE
        if normalized == "shadow":
            return cls.SHADOW
        raise ValueError(f"Unsupported agent runtime v2 mode: {value!r}")


class ProviderRole(str, Enum):
    CONTROLLER = "controller"
    SYNTHESIS_SPECIALIST = "synthesis_specialist"


class TaskPhase(str, Enum):
    ROUTE = "route"
    PROJECT = "project"
    PROJECT_READ = "project_read"
    PROJECT_WRITE = "project_write"
    SYNTHESIS = "synthesis"
    VALIDATION = "validation"
    REPAIR = "repair"
    EXECUTION = "execution"
    WAITING_USER = "waiting_user"
    COMPLETE = "complete"
    BLOCKED = "blocked"


class AgentAction(str, Enum):
    ANSWER = "answer"
    ASK_USER = "ask_user"
    BUILD_PROJECT = "build_project"
    READ_PROJECT = "read_project"
    UPDATE_PROJECT = "update_project"
    SYNTHESIZE_COMMAND = "synthesize_command"
    REPAIR_COMMAND = "repair_command"
    VALIDATE = "validate"
    EXECUTE = "execute"


class ExecutionMode(str, Enum):
    NONE = "none"
    TEST_FAKE = "test_fake"
    LOCAL = "local"
    HPC = "hpc"


class RuntimeContract(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class ArtifactRef(RuntimeContract):
    artifact_id: str
    kind: str
    path: str
    sha256: str
    size_bytes: int = Field(ge=0)
    producer_tool: str = ""
    metadata: dict[str, Any] = Field(default_factory=dict)


class WorkspaceRef(RuntimeContract):
    name: str
    program: str = ""
    path: str
    sha256: str


class TaskEnvelope(RuntimeContract):
    task_id: str
    session_id: str
    turn_id: str
    request: str
    cwd: str
    provider_role: ProviderRole
    phase: TaskPhase
    execution_mode: ExecutionMode = ExecutionMode.NONE
    project: WorkspaceRef | None = None
    server: WorkspaceRef | None = None
    previous_command: str = ""
    unresolved_slots: tuple[str, ...] = ()


class AgentDecision(RuntimeContract):
    """Auditable provider decision; private chain-of-thought is excluded."""

    action: AgentAction
    phase: TaskPhase
    confidence: Literal["low", "medium", "high"] = "medium"
    summary: str
    evidence: tuple[str, ...] = ()
    required_slots: tuple[str, ...] = ()
    requested_tools: tuple[str, ...] = ()


class ToolReceipt(RuntimeContract):
    request_id: str
    tool_name: str
    status: Literal["ok", "error", "denied", "skipped"]
    normalized_args: dict[str, Any] = Field(default_factory=dict)
    artifacts: tuple[ArtifactRef, ...] = ()
    rule_ids: tuple[str, ...] = ()
    state_delta: dict[str, Any] = Field(default_factory=dict)


__all__ = [
    "AgentAction",
    "AgentDecision",
    "ArtifactRef",
    "ExecutionMode",
    "ProviderRole",
    "RuntimeV2Mode",
    "TaskEnvelope",
    "TaskPhase",
    "ToolReceipt",
    "WorkspaceRef",
]
