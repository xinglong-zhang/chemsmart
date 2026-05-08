"""Prompt shown when resuming a session from a different cwd."""

from __future__ import annotations

from dataclasses import dataclass

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static


@dataclass(slots=True)
class CwdMismatchChoice:
    choice: str


class CwdMismatchOverlay(ModalScreen[CwdMismatchChoice | None]):
    BINDINGS = [
        Binding("c", "resume_from_recorded", "Recorded cwd", show=False),
        Binding("i", "ignore_and_continue", "Ignore", show=False),
        Binding("n", "cancel_resume", "Cancel", show=False),
        Binding("escape", "cancel_resume", "Cancel", show=False),
    ]

    DEFAULT_CSS = """
    CwdMismatchOverlay {
        align: center middle;
    }

    #cwd-mismatch-modal {
        width: 92;
        height: auto;
        border: round $warning;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(self, *, recorded_cwd: str, current_cwd: str) -> None:
        super().__init__()
        self.recorded_cwd = recorded_cwd
        self.current_cwd = current_cwd

    def compose(self) -> ComposeResult:
        text = (
            "Resume from a different working directory?\n\n"
            f"recorded: {self.recorded_cwd}\n"
            f"current:  {self.current_cwd}\n\n"
            "c recorded cwd · i ignore and continue here · n cancel"
        )
        with Vertical(id="cwd-mismatch-modal"):
            summary = Static(text)
            summary.border_title = "cwd mismatch"
            yield summary

    def action_resume_from_recorded(self) -> None:
        self.dismiss(CwdMismatchChoice("c"))

    def action_ignore_and_continue(self) -> None:
        self.dismiss(CwdMismatchChoice("i"))

    def action_cancel_resume(self) -> None:
        self.dismiss(CwdMismatchChoice("n"))
