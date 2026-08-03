"""Dependency-free Runtime V2 nucleus for the v3.1.4 agent port."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from chemsmart.agent._contracts import ContractError, canonical_sha256


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
        raise ContractError(f"unsupported Runtime V2 mode: {value!r}")


class TaskPhase(str, Enum):
    ROUTE = "route"
    SPECIFY = "specify"
    PROJECT = "project"
    COMPILE = "compile"
    PREFLIGHT = "preflight"
    WAITING_APPROVAL = "waiting_approval"
    EXECUTION = "execution"
    VALIDATION = "validation"
    REVIEW = "review"
    REPORT = "report"
    TERMINAL = "terminal"


class TerminalState(str, Enum):
    COMPLETE = "complete"
    FAILED = "failed"
    BLOCKED = "blocked"
    WAITING_APPROVAL = "waiting_for_approval"


@dataclass(frozen=True)
class ResourceBudgetV1:
    max_input_tokens_per_request: int
    max_output_tokens_per_request: int
    max_tool_calls: int
    wall_time_seconds: float
    max_cost_usd: float | None = None
    chemistry_engine_calls: int = 0
    hpc_calls: int = 0

    def __post_init__(self) -> None:
        if min(
            self.max_input_tokens_per_request,
            self.max_output_tokens_per_request,
            self.max_tool_calls,
        ) < 1:
            raise ContractError("token, tool, and wall-time budgets must be positive")
        if self.wall_time_seconds <= 0:
            raise ContractError("wall_time_seconds must be positive")
        if self.max_cost_usd is not None and self.max_cost_usd < 0:
            raise ContractError("max_cost_usd must be non-negative")
        if self.chemistry_engine_calls < 0 or self.hpc_calls < 0:
            raise ContractError("compute call budgets must be non-negative")


@dataclass(frozen=True)
class TaskEnvelopeV1:
    schema_version: str
    task_id: str
    session_id: str
    turn_id: str
    request_sha256: str
    cwd_sha256: str
    phase: TaskPhase
    budget: ResourceBudgetV1
    tool_schema_sha256: str
    envelope_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.task-envelope.v1":
            raise ContractError("unsupported task envelope schema")
        body = {
            "schema_version": self.schema_version,
            "task_id": self.task_id,
            "session_id": self.session_id,
            "turn_id": self.turn_id,
            "request_sha256": self.request_sha256,
            "cwd_sha256": self.cwd_sha256,
            "phase": self.phase,
            "budget": self.budget,
            "tool_schema_sha256": self.tool_schema_sha256,
        }
        if self.envelope_sha256 != canonical_sha256(body):
            raise ContractError("task envelope digest mismatch")


@dataclass(frozen=True)
class ProviderStateRefV1:
    """Opaque, adapter-owned continuation identity with no evidence authority."""

    schema_version: str
    provider: str
    session_id: str
    turn_id: str
    state_id: str
    sanitized_history_sha256: str
    tool_schema_sha256: str
    evidentiary: bool = False

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.provider-state-ref.v1":
            raise ContractError("unsupported provider state ref schema")
        if self.evidentiary:
            raise ContractError("provider continuation state is non-evidentiary")


__all__ = [
    "ProviderStateRefV1",
    "ResourceBudgetV1",
    "RuntimeV2Mode",
    "TaskEnvelopeV1",
    "TaskPhase",
    "TerminalState",
]
