"""Stable data contracts shared by ChemSmart agent frontends and services."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

UTC = timezone.utc


def utc_now_iso() -> str:
    """Return an explicit UTC timestamp for persisted agent records."""
    return datetime.now(UTC).isoformat()


class Step(BaseModel):
    tool: str
    args: dict[str, Any] = Field(default_factory=dict)
    rationale: str = ""


class Plan(BaseModel):
    steps: list[Step]
    rationale: str = ""
    estimated_cost: str | None = None
    intent: Literal["workflow", "advisory", "chitchat"] | None = None

    def resolved_intent(self) -> Literal["workflow", "advisory", "chitchat"]:
        if self.intent in {"workflow", "advisory", "chitchat"}:
            return self.intent
        if not self.steps:
            return "advisory"
        return "workflow"

    def is_chitchat(self) -> bool:
        return self.resolved_intent() == "chitchat"


class CriticVerdict(BaseModel):
    verdict: Literal["ok", "warn", "reject"]
    confidence: float = Field(default=0.5, ge=0.0, le=1.0)
    issues: list[str] = Field(default_factory=list)
    rationale: str = ""


class SessionState(BaseModel):
    model_config = ConfigDict(arbitrary_types_allowed=True)

    session_id: str
    cwd: str
    started_at: str = Field(default_factory=utc_now_iso)
    request_started_at: str = Field(default_factory=utc_now_iso)
    turn_index: int = 1
    request_intent: str = "unknown"
    total_steps_planned: int = 0
    current_step_index: int = 0
    plan: Plan | None = None
    pending_ask_user: dict[str, Any] | None = None
    pending_messages: list[dict[str, Any]] | None = None
    request: str | None = None
    env_snapshot: dict[str, str | None] = Field(default_factory=dict)

    def save(self, path: Path) -> None:
        path.write_text(self.model_dump_json(indent=2))

    @classmethod
    def load(cls, path: Path) -> "SessionState":
        return cls.model_validate_json(path.read_text())


__all__ = [
    "CriticVerdict",
    "Plan",
    "SessionState",
    "Step",
    "utc_now_iso",
]
