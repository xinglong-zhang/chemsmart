"""Phase state for the agent TUI."""

from __future__ import annotations

from enum import Enum


class Phase(str, Enum):
    IDLE = "idle"
    PLANNING = "planning"
    DRY_RUN_READY = "dry-run-ready"
    FINISHED = "finished"
    ERROR = "error"

    @property
    def label(self) -> str:
        return {
            Phase.IDLE: "◌ idle",
            Phase.PLANNING: "◐ planning",
            Phase.DRY_RUN_READY: "◑ dry-run-ready",
            Phase.FINISHED: "✓ finished",
            Phase.ERROR: "✗ error",
        }[self]
