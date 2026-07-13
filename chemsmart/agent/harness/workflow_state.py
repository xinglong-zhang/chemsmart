"""Observable workspace state shared by agent command and YAML tools.

The tool registry exposes stateless Python callables, while a research session
is stateful: reading or writing a project YAML selects it for later command
synthesis.  This module keeps that small piece of state deterministic and
scoped to the current working directory.
"""

from __future__ import annotations

import hashlib
import re
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Iterator

from chemsmart.settings.workspace_project import workspace_project_path


@dataclass(frozen=True)
class WorkspaceSelection:
    name: str
    program: str
    path: str
    sha256: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class WorkflowState:
    cwd: str
    project: WorkspaceSelection | None = None
    server: WorkspaceSelection | None = None
    previous_command: str = ""
    unresolved_slots: tuple[str, ...] = ()
    asked_slots: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "cwd": self.cwd,
            "project": self.project.to_dict() if self.project else None,
            "server": self.server.to_dict() if self.server else None,
            "previous_command": self.previous_command,
            "unresolved_slots": list(self.unresolved_slots),
            "asked_slots": list(self.asked_slots),
        }


_DEFAULT_SCOPE = "default"
_CURRENT_SCOPE: ContextVar[str] = ContextVar(
    "chemsmart_agent_workflow_scope",
    default=_DEFAULT_SCOPE,
)
_STATES: dict[tuple[str, str], WorkflowState] = {}


def current_workflow_scope() -> str:
    """Return the active agent-session scope for workspace state."""

    return _CURRENT_SCOPE.get()


@contextmanager
def workflow_state_scope(
    scope: str,
    *,
    cwd: str | Path | None = None,
    inherit_default: bool = True,
) -> Iterator[None]:
    """Isolate project/clarification state for one agent session.

    A newly-created session inherits the TUI's default workspace selection once,
    then evolves independently from other sessions in the same working directory.
    """

    normalized_scope = str(scope or _DEFAULT_SCOPE).strip() or _DEFAULT_SCOPE
    cwd_key = _cwd_key(cwd)
    target_key = (normalized_scope, cwd_key)
    if inherit_default and target_key not in _STATES:
        default_state = _STATES.get((_DEFAULT_SCOPE, cwd_key))
        if default_state is not None:
            _STATES[target_key] = replace(default_state)
    token = _CURRENT_SCOPE.set(normalized_scope)
    try:
        yield
    finally:
        _CURRENT_SCOPE.reset(token)


def current_workflow_state(
    cwd: str | Path | None = None,
) -> WorkflowState:
    cwd_key = _cwd_key(cwd)
    key = (current_workflow_scope(), cwd_key)
    state = _STATES.get(key)
    if state is None:
        state = WorkflowState(cwd=cwd_key)
        _STATES[key] = state
    return state


def reset_workflow_state(cwd: str | Path | None = None) -> None:
    if cwd is None:
        _STATES.clear()
        return
    _STATES.pop((current_workflow_scope(), _cwd_key(cwd)), None)


def select_workspace_project(
    project: str,
    program: str,
    *,
    cwd: str | Path | None = None,
) -> dict[str, Any]:
    base = Path(cwd or Path.cwd()).resolve()
    normalized_program = str(program or "").strip().lower()
    path = workspace_project_path(project, normalized_program, cwd=base)
    if not path.is_file():
        return {
            "selected": False,
            "project": str(project or ""),
            "program": normalized_program,
            "path": str(path),
            "rule_id": "workflow.project.not_found",
        }
    selection = WorkspaceSelection(
        name=path.stem,
        program=normalized_program,
        path=str(path),
        sha256=_file_sha256(path),
    )
    key = (current_workflow_scope(), str(base))
    state = current_workflow_state(base)
    _STATES[key] = replace(state, project=selection)
    return {
        "selected": True,
        "project": selection.name,
        "program": selection.program,
        "path": selection.path,
        "sha256": selection.sha256,
    }


def select_project_from_request(
    request: str,
    *,
    cwd: str | Path | None = None,
) -> dict[str, Any] | None:
    """Select an explicitly named project when it exists in the workspace."""

    name = project_name_from_request(request)
    if not name:
        return None
    base = Path(cwd or Path.cwd()).resolve()
    program_hint = _program_from_request(request)
    programs = (program_hint,) if program_hint else ("gaussian", "orca")
    for program in programs:
        result = select_workspace_project(name, program, cwd=base)
        if result.get("selected"):
            return result
    return {
        "selected": False,
        "project": name,
        "program": program_hint or "",
        "rule_id": "workflow.project.not_found",
    }


def project_name_from_request(request: str) -> str:
    text = str(request or "")
    patterns = (
        r"(?:^|\s)(?:-p|--project)(?:=|\s+)([A-Za-z0-9_.-]+)",
        r"\busing\s+(?:the\s+)?([A-Za-z0-9_.-]+)\s+project\b",
        r"\bwith\s+(?:the\s+)?([A-Za-z0-9_.-]+)\s+project\b",
        r"\b([A-Za-z0-9_.-]+)\s+project\s+settings\b",
        r"\bproject(?:\s+named|\s+name)?\s+([A-Za-z0-9_.-]+)\b",
    )
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if not match:
            continue
        candidate = Path(match.group(1)).stem.strip("._-")
        if candidate.lower() not in {"yaml", "settings", "configuration"}:
            return candidate
    return ""


def record_command(
    command: str,
    *,
    cwd: str | Path | None = None,
) -> None:
    cwd_key = _cwd_key(cwd)
    key = (current_workflow_scope(), cwd_key)
    state = current_workflow_state(cwd_key)
    _STATES[key] = replace(state, previous_command=str(command or "").strip())


def record_clarification_slots(
    slots: list[str] | tuple[str, ...],
    *,
    cwd: str | Path | None = None,
) -> tuple[str, ...]:
    """Record only newly requested slots and return those new slots."""

    cwd_key = _cwd_key(cwd)
    key = (current_workflow_scope(), cwd_key)
    state = current_workflow_state(cwd_key)
    normalized = tuple(dict.fromkeys(str(item).strip() for item in slots if str(item).strip()))
    new = tuple(item for item in normalized if item not in state.asked_slots)
    _STATES[key] = replace(
        state,
        unresolved_slots=normalized,
        asked_slots=tuple(dict.fromkeys((*state.asked_slots, *normalized))),
    )
    return new


def clear_resolved_slots(*, cwd: str | Path | None = None) -> None:
    cwd_key = _cwd_key(cwd)
    key = (current_workflow_scope(), cwd_key)
    state = current_workflow_state(cwd_key)
    _STATES[key] = replace(state, unresolved_slots=(), asked_slots=())


def _program_from_request(request: str) -> str:
    lowered = str(request or "").lower()
    if "gaussian" in lowered:
        return "gaussian"
    if "orca" in lowered:
        return "orca"
    return ""


def _cwd_key(cwd: str | Path | None) -> str:
    return str(Path(cwd or Path.cwd()).resolve())


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


__all__ = [
    "WorkflowState",
    "WorkspaceSelection",
    "clear_resolved_slots",
    "current_workflow_scope",
    "current_workflow_state",
    "project_name_from_request",
    "record_clarification_slots",
    "record_command",
    "reset_workflow_state",
    "select_project_from_request",
    "select_workspace_project",
    "workflow_state_scope",
]
