"""Approval overlay shown before risky execution."""

from __future__ import annotations

from dataclasses import dataclass

from textual import events
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
        Binding("y", "approve_once", "Yes once", show=False, priority=True),
        Binding("n", "deny", "No", show=False, priority=True),
        Binding("s", "approve_session", "Session", show=False, priority=True),
        Binding("r", "revise", "Revise", show=False, priority=True),
        Binding("escape", "cancel", "Cancel", show=False, priority=True),
        Binding(
            "enter", "submit_revision", "Submit", show=False, priority=True
        ),
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
        mode, remote_effect, rollback = _approval_copy(self.action)
        summary = (
            f"Approve `{self.action}` for this request?\n\n"
            f"mode: {mode}\n"
            f"target: {self.action}\n"
            f"remote-effect: {remote_effect}\n"
            f"rollback: {rollback}\n\n"
            "y once · n deny · s this-session · r decline-and-revise\n\n"
            f"request: {self.request}"
        )
        with Vertical(id="approval-modal"):
            yield Static(summary, id="approval-summary")
            yield Input(
                placeholder="Add a corrective instruction, then press Enter",
                id="approval-revision",
                disabled=True,
            )

    def on_mount(self) -> None:
        self.query_one("#approval-summary", Static).border_title = "Approval"

    def on_key(self, event: events.Key) -> None:
        if self._revising:
            if event.key == "enter":
                event.stop()
                self.action_submit_revision()
            elif event.key == "escape":
                event.stop()
                self.action_cancel()
            return
        if event.key == "y":
            event.stop()
            self.action_approve_once()
        elif event.key == "n":
            event.stop()
            self.action_deny()
        elif event.key == "s":
            event.stop()
            self.action_approve_session()
        elif event.key == "r":
            event.stop()
            self.action_revise()
        elif event.key == "escape":
            event.stop()
            self.action_cancel()

    def action_approve_once(self) -> None:
        self.dismiss(ApprovalResult("y"))

    def action_deny(self) -> None:
        self.dismiss(ApprovalResult("n"))

    def action_approve_session(self) -> None:
        self.dismiss(ApprovalResult("s"))

    def action_revise(self) -> None:
        self._revising = True
        self.add_class("-revising")
        revision = self.query_one("#approval-revision", Input)
        revision.disabled = False
        revision.focus()

    def action_submit_revision(self) -> None:
        if not self._revising:
            return
        text = self.query_one("#approval-revision", Input).value.strip()
        if not text:
            return
        self.dismiss(ApprovalResult("r", corrective_text=text))

    def action_cancel(self) -> None:
        self.dismiss(None)


def _approval_copy(action: str) -> tuple[str, str, str]:
    if action == "submit_hpc":
        return (
            "execute once",
            "This may submit to the remote queue.",
            "After submission, cancellation requires a separate queue cancel.",
        )
    if action == "run_local":
        return (
            "execute once",
            "A local run will start in the current working directory.",
            "Reject now to prevent it; cleanup after start is handled separately.",
        )
    return (
        "execute once",
        "The requested execution step will start.",
        "You can reject before it starts; cleanup after start is handled separately.",
    )
