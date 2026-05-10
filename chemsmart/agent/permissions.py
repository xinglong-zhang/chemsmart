from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from chemsmart.agent.provider_adapter import ToolRequest


class PermissionMode(str, Enum):
    PERMISSION = "permission"
    DRIVING = "driving"


class ApprovalDecision(str, Enum):
    ALLOW_ONCE = "allow_once"
    ALLOW_SESSION = "allow_session"
    DENY = "deny"


class ResolvedDecision(str, Enum):
    AUTO_ALLOW = "auto_allow"
    NEEDS_USER = "needs_user"
    AUTO_DENY = "auto_deny"


DRIVING_DEFAULT_DENY = {"run_local", "submit_hpc"}


@dataclass(frozen=True)
class ResolvedPermission:
    decision: ResolvedDecision
    reason: str


@dataclass
class PermissionPolicy:
    mode: PermissionMode
    yolo: bool = False
    session_allow: set[str] = field(default_factory=set)
    driving_denylist: set[str] = field(
        default_factory=lambda: set(DRIVING_DEFAULT_DENY)
    )

    def resolve(self, req: ToolRequest) -> ResolvedPermission:
        tool = req.name
        if self.mode == PermissionMode.DRIVING:
            if tool in self.driving_denylist:
                if self.yolo:
                    return ResolvedPermission(
                        decision=ResolvedDecision.AUTO_ALLOW,
                        reason="yolo",
                    )
                return ResolvedPermission(
                    decision=ResolvedDecision.AUTO_DENY,
                    reason="missing_yolo",
                )
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_ALLOW,
                reason="driving_mode",
            )

        if tool in self.session_allow:
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_ALLOW,
                reason="session_rule",
            )
        return ResolvedPermission(
            decision=ResolvedDecision.NEEDS_USER,
            reason="needs_user",
        )

    def record(self, tool: str, decision: ApprovalDecision) -> None:
        if decision == ApprovalDecision.ALLOW_SESSION:
            self.session_allow.add(tool)
