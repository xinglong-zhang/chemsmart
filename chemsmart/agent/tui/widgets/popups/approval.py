"""Approval overlay shown before risky execution."""

from __future__ import annotations

from dataclasses import dataclass

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Input, Static


@dataclass(slots=True)
class ApprovalResult:
    choice: str
    corrective_text: str | None = None


class ApprovalOverlay(ModalScreen[ApprovalResult | None]):
    BINDINGS = [
        Binding("y", "approve_once", "Yes once", show=False),
        Binding("n", "deny", "No", show=False),
        Binding("s", "approve_session", "Session", show=False),
        Binding("r", "revise", "Revise", show=False),
        Binding("escape", "cancel", "Cancel", show=False),
        Binding("enter", "submit_revision", "Submit", show=False),
    ]

    DEFAULT_CSS = """
    ApprovalOverlay {
        align: center middle;
    }

    #approval-modal {
        width: 88;
        height: auto;
        max-height: 22;
        border: round $warning;
        padding: 1 2;
        background: $surface;
    }

    #approval-revision {
        margin-top: 1;
        display: none;
    }

    ApprovalOverlay.-revising #approval-revision {
        display: block;
    }
    """

    def __init__(self, *, action: str, request: str) -> None:
        super().__init__()
        self.action = action
        self.request = request
        self._revising = False

    def compose(self) -> ComposeResult:
        summary = (
            f"Approve `{self.action}` for this request?\n\n"
            "y once · n deny · s this-session · r decline-and-revise\n\n"
            f"Request: {self.request}"
        )
        with Vertical(id="approval-modal"):
            yield Static(summary, id="approval-summary")
            yield Input(
                placeholder="Add a corrective instruction, then press Enter",
                id="approval-revision",
            )

    def on_mount(self) -> None:
        self.query_one("#approval-summary", Static).border_title = "Approval"

    def action_approve_once(self) -> None:
        self.dismiss(ApprovalResult("y"))

    def action_deny(self) -> None:
        self.dismiss(ApprovalResult("n"))

    def action_approve_session(self) -> None:
        self.dismiss(ApprovalResult("s"))

    def action_revise(self) -> None:
        self._revising = True
        self.add_class("-revising")
        self.query_one("#approval-revision", Input).focus()

    def action_submit_revision(self) -> None:
        if not self._revising:
            return
        text = self.query_one("#approval-revision", Input).value.strip()
        if not text:
            return
        self.dismiss(ApprovalResult("r", corrective_text=text))

    def action_cancel(self) -> None:
        self.dismiss(None)
