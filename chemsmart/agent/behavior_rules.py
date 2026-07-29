"""CHEMSMART.md agent rules: loading, merging, and read/write tools.

CHEMSMART.md is the agent's user-editable initial-rules file, in the spirit
of CLAUDE.md/AGENTS.md: plain markdown that shapes answer style, verbosity,
creativity, answer language, and the xtb_real_runs execution policy. Two
scopes are merged on every turn:

- user global: ``~/.chemsmart/CHEMSMART.md`` (personal preferences)
- workspace:   ``./CHEMSMART.md`` (project rules; wins on conflicts)

Free-form prose is passed verbatim into the system prompt (bounded). Only
lines inside a ``## Policy`` section are machine-parsed, as ``key: value``
pairs, so execution policy stays deterministic and never depends on a model
interpreting prose.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

RULES_FILENAME = "CHEMSMART.md"
MAX_RULES_CHARS_PER_FILE = 4000
TRUNCATION_MARKER = "\n… [truncated]"

_POLICY_HEADING_RE = re.compile(r"^##\s+policy\s*$", re.IGNORECASE)
_HEADING_RE = re.compile(r"^#{1,6}\s")
_POLICY_LINE_RE = re.compile(r"^[-*]?\s*([A-Za-z_][A-Za-z0-9_-]*)\s*:\s*(.+)$")


@dataclass(frozen=True)
class BehaviorRules:
    """Merged CHEMSMART.md rules for one session turn."""

    text: str = ""
    policy: dict[str, str] = field(default_factory=dict)
    workspace_path: Path | None = None
    user_path: Path | None = None
    truncated: bool = False

    @property
    def loaded(self) -> bool:
        return bool(self.text)


def workspace_rules_path(cwd: str | Path | None = None) -> Path:
    return Path(cwd or Path.cwd()).resolve() / RULES_FILENAME


def user_rules_path() -> Path:
    return Path.home() / ".chemsmart" / RULES_FILENAME


def load_behavior_rules(cwd: str | Path | None = None) -> BehaviorRules:
    """Load and merge user-global and workspace CHEMSMART.md files.

    Both files are size-bounded before prompt injection. The workspace file
    is appended after the user file (so it dominates prompt attention) and
    its ``## Policy`` keys override user-level keys.
    """

    user_path = user_rules_path()
    workspace_path = workspace_rules_path(cwd)
    sections: list[str] = []
    policy: dict[str, str] = {}
    truncated = False
    found_user: Path | None = None
    found_workspace: Path | None = None

    for scope, path in (("user", user_path), ("workspace", workspace_path)):
        text = _read_rules_file(path)
        if text is None:
            continue
        if scope == "user":
            found_user = path
        else:
            found_workspace = path
        clamped, was_truncated = _clamp(text)
        truncated = truncated or was_truncated
        sections.append(f"[{scope} rules — {path}]\n{clamped}")
        policy.update(parse_policy_block(text))

    return BehaviorRules(
        text="\n\n".join(sections),
        policy=policy,
        workspace_path=found_workspace,
        user_path=found_user,
        truncated=truncated,
    )


def parse_policy_block(text: str) -> dict[str, str]:
    """Extract ``key: value`` pairs from the ``## Policy`` section."""

    policy: dict[str, str] = {}
    in_policy = False
    for line in text.splitlines():
        stripped = line.strip()
        if _POLICY_HEADING_RE.match(stripped):
            in_policy = True
            continue
        if in_policy and _HEADING_RE.match(stripped):
            in_policy = False
            continue
        if not in_policy or not stripped:
            continue
        match = _POLICY_LINE_RE.match(stripped)
        if match is not None:
            policy[match.group(1).strip().lower()] = match.group(2).strip()
    return policy


def read_behavior_rules() -> dict[str, Any]:
    """Read the current CHEMSMART.md rules files and parsed policy."""

    rules = load_behavior_rules()
    return {
        "workspace_path": (
            str(rules.workspace_path) if rules.workspace_path else None
        ),
        "user_path": str(rules.user_path) if rules.user_path else None,
        "policy": dict(rules.policy),
        "content": rules.text,
        "truncated_for_prompt": rules.truncated,
    }


def write_behavior_rules(
    content: str,
    scope: str = "workspace",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write a CHEMSMART.md rules file after explicit user approval."""

    normalized_scope = str(scope or "workspace").strip().lower()
    if normalized_scope not in {"workspace", "user"}:
        return {
            "error": "unknown_scope",
            "message": "scope must be 'workspace' or 'user'",
        }
    path = (
        user_rules_path()
        if normalized_scope == "user"
        else workspace_rules_path()
    )
    if path.exists() and not overwrite:
        return {
            "error": "file_exists",
            "path": str(path),
            "message": "pass overwrite=true to replace the existing file",
        }
    text = str(content or "")
    if not text.strip():
        return {
            "error": "empty_content",
            "message": "refusing to write an empty CHEMSMART.md",
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return {
        "path": str(path),
        "scope": normalized_scope,
        "chars": len(text),
        "policy": parse_policy_block(text),
        "prompt_bounded": len(text) <= MAX_RULES_CHARS_PER_FILE,
    }


def _read_rules_file(path: Path) -> str | None:
    try:
        if not path.is_file():
            return None
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None


def _clamp(text: str) -> tuple[str, bool]:
    if len(text) <= MAX_RULES_CHARS_PER_FILE:
        return text, False
    kept = MAX_RULES_CHARS_PER_FILE - len(TRUNCATION_MARKER)
    return text[:kept] + TRUNCATION_MARKER, True
