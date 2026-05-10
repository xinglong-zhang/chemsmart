"""Approval and permission-mode popups."""

from __future__ import annotations

from dataclasses import dataclass

from textual import events
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static

from chemsmart.agent.permissions import ApprovalDecision, PermissionMode
from chemsmart.agent.tui.tool_meta import (
    pretty_tool_args,
    tool_risk_badge,
    tool_side_effect_summary,
)


@dataclass(slots=True, frozen=True)
class ApprovalResult:
    choice: str

    def to_decision(self) -> ApprovalDecision | None:
        return {
            "y": ApprovalDecision.ALLOW_ONCE,
            "s": ApprovalDecision.ALLOW_SESSION,
            "n": ApprovalDecision.DENY,
        }.get(self.choice)


@dataclass(slots=True, frozen=True)
class PermissionModeResult:
    mode: PermissionMode
    yolo: bool


class ApprovalOverlay(ModalScreen[ApprovalResult | None]):
    BINDINGS = [
        Binding("y", "approve_once", "Yes once", show=False, priority=True),
        Binding("s", "approve_session", "Session", show=False, priority=True),
        Binding("n", "deny", "No", show=False, priority=True),
        Binding("escape", "cancel", "Cancel", show=False, priority=True),
    ]

    DEFAULT_CSS = """
    ApprovalOverlay {
        align: center middle;
    }

    #approval-modal {
        width: 92;
        height: auto;
        max-height: 24;
        border: round $warning;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(
        self,
        *,
        action: str | None = None,
        request: str | None = None,
        tool_name: str = "",
        description: str = "",
        arguments: dict | None = None,
        risk_badge: str = "",
        side_effect_summary: str = "",
        session_rule_active: bool = False,
        queue_index: int | None = None,
        queue_total: int | None = None,
    ) -> None:
        super().__init__()
        tool_name = tool_name or action or ""
        description = description or tool_name
        arguments = dict(arguments or {})
        risk_badge = risk_badge or tool_risk_badge(tool_name)[0]
        side_effect_summary = side_effect_summary or tool_side_effect_summary(
            tool_name
        )
        self.tool_name = tool_name
        self.description = description
        self.arguments = arguments
        self.risk_badge = risk_badge
        self.side_effect_summary = side_effect_summary
        self.session_rule_active = session_rule_active
        self.queue_index = queue_index
        self.queue_total = queue_total
        self.request = request or ""

    def compose(self) -> ComposeResult:
        queue_text = ""
        if (
            self.queue_index is not None
            and self.queue_total is not None
            and self.queue_total > 0
        ):
            queue_text = f" ({self.queue_index} of {self.queue_total})"

        session_rule = "active" if self.session_rule_active else "none"
        summary = (
            f"Approve tool use{queue_text}?\n\n"
            f"tool: {self.tool_name}\n"
            f"description: {self.description}\n"
            f"risk: {self.risk_badge}\n"
            f"mode: execute once\n"
            f"remote-effect: {self.side_effect_summary}\n"
            "rollback: deny before execution; cleanup after start is separate.\n"
            f"effect: {self.side_effect_summary}\n"
            f"session rule: {session_rule}\n\n"
            "args:\n"
            f"{pretty_tool_args(self.arguments)}\n\n"
            "y once · s this-session · n deny · esc cancel"
        )
        if self.request:
            summary += f"\n\nrequest: {self.request}"
        with Vertical(id="approval-modal"):
            yield Static(summary, id="approval-summary")

    def on_mount(self) -> None:
        self.query_one("#approval-summary", Static).border_title = "Approval"

    def on_key(self, event: events.Key) -> None:
        if event.key == "y":
            event.stop()
            self.action_approve_once()
        elif event.key == "s":
            event.stop()
            self.action_approve_session()
        elif event.key == "n":
            event.stop()
            self.action_deny()
        elif event.key == "escape":
            event.stop()
            self.action_cancel()

    def action_approve_once(self) -> None:
        self.dismiss(ApprovalResult("y"))

    def action_approve_session(self) -> None:
        self.dismiss(ApprovalResult("s"))

    def action_deny(self) -> None:
        self.dismiss(ApprovalResult("n"))

    def action_cancel(self) -> None:
        self.dismiss(None)


class PermissionModeOverlay(ModalScreen[PermissionModeResult | None]):
    BINDINGS = [
        Binding("p", "set_permission", "Permission", show=False),
        Binding("d", "set_driving", "Driving", show=False),
        Binding("y", "toggle_yolo", "YOLO", show=False),
        Binding("enter", "apply", "Apply", show=False),
        Binding("escape", "cancel", "Cancel", show=False),
    ]

    DEFAULT_CSS = """
    PermissionModeOverlay {
        align: center middle;
    }

    #permissions-modal {
        width: 74;
        height: auto;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(
        self,
        *,
        mode: PermissionMode,
        yolo: bool,
    ) -> None:
        super().__init__()
        self.mode = mode
        self.yolo = yolo

    def compose(self) -> ComposeResult:
        with Vertical(id="permissions-modal"):
            yield Static(self._summary_text(), id="permissions-summary")

    def on_mount(self) -> None:
        self.query_one("#permissions-summary", Static).border_title = (
            "Permissions"
        )

    def action_set_permission(self) -> None:
        self.mode = PermissionMode.PERMISSION
        self._refresh()

    def action_set_driving(self) -> None:
        self.mode = PermissionMode.DRIVING
        self._refresh()

    def action_toggle_yolo(self) -> None:
        self.yolo = not self.yolo
        self._refresh()

    def action_apply(self) -> None:
        self.dismiss(PermissionModeResult(mode=self.mode, yolo=self.yolo))

    def action_cancel(self) -> None:
        self.dismiss(None)

    def _refresh(self) -> None:
        self.query_one("#permissions-summary", Static).update(
            self._summary_text()
        )

    def _summary_text(self) -> str:
        return (
            "Toggle permission policy.\n\n"
            f"mode: {self.mode.value}\n"
            f"yolo: {'on' if self.yolo else 'off'}\n\n"
            "p permission · d driving · y toggle yolo · enter apply"
        )


def build_approval_overlay(
    *,
    tool_name: str,
    description: str,
    arguments: dict,
    session_rule_active: bool = False,
    queue_index: int | None = None,
    queue_total: int | None = None,
) -> ApprovalOverlay:
    risk_badge, _ = tool_risk_badge(tool_name)
    return ApprovalOverlay(
        tool_name=tool_name,
        description=description,
        arguments=arguments,
        risk_badge=risk_badge,
        side_effect_summary=tool_side_effect_summary(tool_name),
        session_rule_active=session_rule_active,
        queue_index=queue_index,
        queue_total=queue_total,
    )
