from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from chemsmart.agent.provider_adapter import ToolRequest


class PermissionMode(str, Enum):
    PERMISSION = "permission"
    DRIVING = "driving"


class RuntimePermissionMode(str, Enum):
    READ_ONLY = "read_only"
    ACCEPT_EDITS = "accept_edits"
    BYPASS = "bypass"
    PLAN = "plan"


class ApprovalDecision(str, Enum):
    ALLOW_ONCE = "allow_once"
    ALLOW_SESSION = "allow_session"
    DENY = "deny"


class ResolvedDecision(str, Enum):
    AUTO_ALLOW = "auto_allow"
    NEEDS_USER = "needs_user"
    AUTO_DENY = "auto_deny"


DRIVING_DEFAULT_DENY = {"run_local", "submit_hpc", "remote_probe"}
ALWAYS_REQUIRE_APPROVAL = {"wizard_write"}
READ_ONLY_TOOLS = {"read", "ssh_probe", "scheduler_query", "log_tail"}
EDIT_SAFE_TOOLS = {"edit", "write"}
PLAN_MODE_REASON = "plan mode active"
NEVER_AUTO_ALLOW_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("pip_install", re.compile(r"\bpip3?\s+(install|uninstall)\b")),
    ("apt_install", re.compile(r"\bapt(-get)?\s+install\b")),
    ("brew_install", re.compile(r"\bbrew\s+install\b")),
    (
        "npm_install_global",
        re.compile(r"\bnpm\s+install\s+(-g|--global)\b"),
    ),
    ("sudo", re.compile(r"\bsudo\b")),
    ("rm_root", re.compile(r"\brm\s+-[rf]+\s+/(?:\s|$)")),
    (
        "curl_pipe_shell",
        re.compile(r"\b(curl|wget)\b[^|]*\|\s*(bash|sh|zsh)\b"),
    ),
    ("chmod_777", re.compile(r"\bchmod\s+(-R\s+)?777\b")),
]

PermissionPolicyMode = PermissionMode | RuntimePermissionMode


@dataclass(frozen=True)
class ResolvedPermission:
    decision: ResolvedDecision
    reason: str


@dataclass
class PermissionPolicy:
    mode: PermissionPolicyMode
    yolo: bool = False
    session_allow: set[str] = field(default_factory=set)
    driving_denylist: set[str] = field(
        default_factory=lambda: set(DRIVING_DEFAULT_DENY)
    )

    def resolve(self, req: ToolRequest) -> ResolvedPermission:
        return resolve(
            req,
            mode=self.mode,
            yolo=self.yolo,
            session_allow=self.session_allow,
            driving_denylist=self.driving_denylist,
        )

    def record(self, tool: str, decision: ApprovalDecision) -> None:
        if tool in ALWAYS_REQUIRE_APPROVAL:
            return
        if decision == ApprovalDecision.ALLOW_SESSION:
            self.session_allow.add(tool)


def _decision_keys(req: ToolRequest) -> set[str]:
    keys = {req.name}
    if req.name == "wizard_probe":
        ssh_host_hint = req.arguments.get("ssh_host_hint")
        keys.add("remote_probe" if ssh_host_hint else "local_probe")
    return keys


def legacy_to_runtime(mode: PermissionMode) -> RuntimePermissionMode:
    return {
        PermissionMode.PERMISSION: RuntimePermissionMode.READ_ONLY,
        PermissionMode.DRIVING: RuntimePermissionMode.ACCEPT_EDITS,
    }[mode]


def _matches_never_auto_allow(
    req: ToolRequest,
) -> tuple[str, str] | None:
    def iter_string_values(value: Any) -> list[str]:
        if isinstance(value, str):
            return [value]
        if isinstance(value, dict):
            return [
                nested
                for item in value.values()
                for nested in iter_string_values(item)
            ]
        if isinstance(value, list):
            return [
                nested for item in value for nested in iter_string_values(item)
            ]
        return []

    string_values = iter_string_values(req.arguments)
    for pattern_id, pattern in NEVER_AUTO_ALLOW_PATTERNS:
        for value in string_values:
            match = pattern.search(value)
            if match is not None:
                return pattern_id, match.group(0)
    return None


def resolve(
    req: ToolRequest,
    mode: PermissionPolicyMode,
    *,
    yolo: bool = False,
    session_allow: set[str] | None = None,
    driving_denylist: set[str] | None = None,
) -> ResolvedPermission:
    never_auto_allow_match = _matches_never_auto_allow(req)
    if never_auto_allow_match is not None:
        pattern_id, _ = never_auto_allow_match
        return ResolvedPermission(
            decision=ResolvedDecision.NEEDS_USER,
            reason=f"never_auto_allow:{pattern_id}",
        )

    if isinstance(mode, RuntimePermissionMode):
        tool = req.name
        if mode == RuntimePermissionMode.BYPASS:
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_ALLOW,
                reason="bypass_mode",
            )
        if mode == RuntimePermissionMode.PLAN:
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_DENY,
                reason=PLAN_MODE_REASON,
            )
        if tool in READ_ONLY_TOOLS:
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_ALLOW,
                reason="read_only_tool",
            )
        if (
            mode == RuntimePermissionMode.ACCEPT_EDITS
            and tool in READ_ONLY_TOOLS | EDIT_SAFE_TOOLS
        ):
            return ResolvedPermission(
                decision=ResolvedDecision.AUTO_ALLOW,
                reason="edit_safe_tool",
            )
        return ResolvedPermission(
            decision=ResolvedDecision.NEEDS_USER,
            reason="needs_user",
        )

    tool = req.name
    if tool in ALWAYS_REQUIRE_APPROVAL:
        return ResolvedPermission(
            decision=ResolvedDecision.NEEDS_USER,
            reason="always_require_approval",
        )

    decision_keys = _decision_keys(req)
    denylist = (
        DRIVING_DEFAULT_DENY if driving_denylist is None else driving_denylist
    )
    if mode == PermissionMode.DRIVING:
        if denylist.intersection(decision_keys):
            if yolo:
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

    if tool in (session_allow or set()):
        return ResolvedPermission(
            decision=ResolvedDecision.AUTO_ALLOW,
            reason="session_rule",
        )
    return ResolvedPermission(
        decision=ResolvedDecision.NEEDS_USER,
        reason="needs_user",
    )
