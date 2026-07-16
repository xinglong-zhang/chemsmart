"""Phase state for the agent TUI."""

from __future__ import annotations

from enum import Enum


class Phase(str, Enum):
    IDLE = "idle"
    PLANNING = "planning"
    TOOL_RUNNING = "tool-running"
    APPROVAL_REQUIRED = "approval-required"
    VALIDATING = "validating"
    WAITING_USER = "waiting-user"
    DRY_RUN_READY = "dry-run-ready"
    RUNNING = "running"
    EXECUTING = "executing"
    SUBMITTING = "submitting"
    FINISHED = "finished"
    COMPLETED = "finished"
    ERROR = "error"
    FAILED = "error"
    INTERRUPTED = "interrupted"

    @property
    def label(self) -> str:
        return {
            Phase.IDLE: "◌ idle",
            Phase.PLANNING: "◐ planning",
            Phase.TOOL_RUNNING: "◐ tool-running",
            Phase.APPROVAL_REQUIRED: "! approval-required",
            Phase.VALIDATING: "◐ validating",
            Phase.WAITING_USER: "? waiting-user",
            Phase.DRY_RUN_READY: "◑ dry-run-ready",
            Phase.RUNNING: "▶ running",
            Phase.EXECUTING: "▶ executing",
            Phase.SUBMITTING: "▶ submitting",
            Phase.FINISHED: "✓ finished",
            Phase.ERROR: "✗ error",
            Phase.INTERRUPTED: "■ interrupted",
        }[self]
