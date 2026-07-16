"""Validated user-facing TUI preferences from ``agent.yaml``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from chemsmart.agent.provider_config import _default_yaml_path

DEFAULT_KEYBINDINGS = {
    "show_shortcuts": "f1",
    "toggle_transcript": "ctrl+o",
    "show_activity": "ctrl+t",
    "show_calculations": "ctrl+b",
    "search_history": "ctrl+r",
    "show_project_yaml": "shift+tab",
}

_ALLOWED_ACTIONS = frozenset(DEFAULT_KEYBINDINGS)
_RESERVED_KEYS = frozenset(
    {"ctrl+c", "ctrl+d", "ctrl+l", "y", "s", "n", "r"}
)


@dataclass(slots=True, frozen=True)
class TuiConfig:
    keybindings: dict[str, str] = field(
        default_factory=lambda: dict(DEFAULT_KEYBINDINGS)
    )
    tool_detail: str = "compact"
    issues: tuple[str, ...] = ()


def load_tui_config(
    yaml_path: str | Path | None = None,
) -> TuiConfig:
    """Load only the non-secret ``tui`` block; invalid entries use defaults."""

    path = Path(yaml_path) if yaml_path is not None else _default_yaml_path()
    try:
        document: Any = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return TuiConfig()
    block = document.get("tui") if isinstance(document, dict) else None
    if not isinstance(block, dict):
        return TuiConfig()

    issues: list[str] = []
    configured = block.get("keybindings", {})
    bindings = dict(DEFAULT_KEYBINDINGS)
    if configured is not None and not isinstance(configured, dict):
        issues.append("tui.keybindings must be a mapping")
    elif isinstance(configured, dict):
        seen_keys: dict[str, str] = {}
        modified: set[str] = set()
        for raw_action, raw_key in configured.items():
            action = str(raw_action).strip()
            key = str(raw_key).strip().lower()
            if action not in _ALLOWED_ACTIONS:
                issues.append(f"unknown TUI action: {action}")
                continue
            if not key or key in _RESERVED_KEYS:
                issues.append(f"reserved or empty TUI key: {key or '<empty>'}")
                continue
            previous = seen_keys.get(key)
            if previous is not None and previous != action:
                issues.append(f"duplicate TUI key {key}: {previous}, {action}")
                continue
            bindings[action] = key
            modified.add(action)
            seen_keys[key] = action

        by_key: dict[str, list[str]] = {}
        for action, key in bindings.items():
            by_key.setdefault(key, []).append(action)
        for key, actions in by_key.items():
            if len(actions) < 2:
                continue
            issues.append(f"duplicate TUI key {key}: {', '.join(actions)}")
            for action in actions:
                if action in modified:
                    bindings[action] = DEFAULT_KEYBINDINGS[action]

    detail = str(block.get("tool_detail", "compact")).strip().lower()
    if detail not in {"compact", "full"}:
        issues.append("tui.tool_detail must be compact or full")
        detail = "compact"
    return TuiConfig(bindings, detail, tuple(issues))
