"""The single registry behind slash commands, help, and key hints.

Every human-facing action is declared once: its slash name, aliases, what it
does, and whether it takes arguments. The dispatcher, the generated /help
text, the footer hints, and the command palette all read this table, so a
command cannot exist in one surface and be missing from another.

This module imports no UI framework; it is pure data and string logic.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CommandSpecV1:
    name: str
    slash: str
    title: str
    category: str
    usage: str = ""
    takes_argument: bool = False
    aliases: tuple[str, ...] = ()

    def matches(self, token: str) -> bool:
        return token == self.slash or token in self.aliases


@dataclass(frozen=True)
class KeybindSpecV1:
    key: str
    action: str
    description: str
    show: bool = True


#: Declared in display order. Handlers are bound by the app; a spec whose
#: handler is not yet mounted must not be listed here.
COMMANDS: tuple[CommandSpecV1, ...] = (
    CommandSpecV1(
        name="approve",
        slash="/approve",
        title="Approve the displayed workflow and run it once",
        category="decision",
    ),
    CommandSpecV1(
        name="deny",
        slash="/deny",
        title="Decline the displayed workflow",
        category="decision",
    ),
    CommandSpecV1(
        name="revise",
        slash="/revise",
        title="Decline and edit the request in the composer",
        category="decision",
    ),
    CommandSpecV1(
        name="status",
        slash="/status",
        title="Show the session phase, task, and pending authority",
        category="session",
    ),
    CommandSpecV1(
        name="capabilities",
        slash="/capabilities",
        title="Show the declared program/engine/job surface",
        category="session",
    ),
    CommandSpecV1(
        name="resume",
        slash="/resume",
        title="Restore this workspace's previous session and any pending review",
        category="session",
    ),
    CommandSpecV1(
        name="skills",
        slash="/skills",
        title="List the domain skills the session can consult",
        category="session",
    ),
    CommandSpecV1(
        name="skill",
        slash="/skill",
        title="Tag the next request with a domain skill to consult",
        category="session",
        usage="/skill <id>",
        takes_argument=True,
    ),
    CommandSpecV1(
        name="export",
        slash="/export",
        title="Save the transcript to a file in the workspace",
        category="view",
    ),
    CommandSpecV1(
        name="dag",
        slash="/dag",
        title="Toggle the workflow panel with every node's status",
        category="view",
    ),
    CommandSpecV1(
        name="report",
        slash="/report",
        title="Open the completed-analysis report",
        category="view",
        usage="/report [n]",
        takes_argument=True,
    ),
    CommandSpecV1(
        name="runs",
        slash="/runs",
        title="List this workspace's executions and their reports",
        category="view",
    ),
    CommandSpecV1(
        name="help",
        slash="/help",
        title="Show every command and keybinding",
        category="app",
    ),
    CommandSpecV1(
        name="quit",
        slash="/quit",
        title="Exit when no host operation is running",
        category="app",
        aliases=("/exit",),
    ),
)

KEYBINDS: tuple[KeybindSpecV1, ...] = (
    KeybindSpecV1("ctrl+p", "command_palette", "Command palette"),
    KeybindSpecV1("pageup", "scroll_transcript_up", "Scroll up", show=False),
    KeybindSpecV1(
        "pagedown", "scroll_transcript_down", "Scroll down", show=False
    ),
    KeybindSpecV1("ctrl+home", "scroll_transcript_top", "Top", show=False),
    KeybindSpecV1(
        "ctrl+end", "scroll_transcript_bottom", "Bottom", show=False
    ),
    KeybindSpecV1("ctrl+c", "safe_quit", "Quit"),
)


def command_for(token: str) -> CommandSpecV1 | None:
    lowered = token.lower()
    for spec in COMMANDS:
        if spec.matches(lowered):
            return spec
    return None


def suggest(token: str) -> str | None:
    """Nearest command for a mistyped slash token, or None."""

    lowered = token.lower().lstrip("/")
    if not lowered:
        return None
    candidates = []
    for spec in COMMANDS:
        for name in (spec.slash, *spec.aliases):
            bare = name.lstrip("/")
            if bare.startswith(lowered) or lowered.startswith(bare):
                candidates.append((abs(len(bare) - len(lowered)), name))
    if not candidates:
        shared = [
            (len(set(lowered) ^ set(spec.slash.lstrip("/"))), spec.slash)
            for spec in COMMANDS
        ]
        shared.sort()
        return shared[0][1] if shared and shared[0][0] <= 3 else None
    candidates.sort()
    return candidates[0][1]


def render_help() -> str:
    """The /help text, generated so it can never omit a command."""

    lines = [
        "Enter a scientific request to create a project-YAML/CLI plan and "
        "safe preview.",
        "",
        "Commands:",
    ]
    width = max(
        len(spec.usage or spec.slash)
        + 2
        + sum(len(a) + 2 for a in spec.aliases)
        for spec in COMMANDS
    )
    for spec in COMMANDS:
        shown = spec.usage or spec.slash
        if spec.aliases:
            shown += " (" + ", ".join(spec.aliases) + ")"
        lines.append(f"  {shown:<{width}}  {spec.title}")
    lines += [
        "",
        "Keys:",
    ]
    for bind in KEYBINDS:
        if bind.show:
            lines.append(f"  {bind.key:<12}  {bind.description}")
    lines += [
        "",
        "The provider/runtime cannot grant itself authority. The single "
        "/approve action ends planning authority and starts the displayed "
        "YAML/CLI DAG through the provider-free ChemSmart executor. Internal "
        "receipts remain provenance; no hash or approval-file token is "
        "required from the human.",
    ]
    return "\n".join(lines)


def footer_hint(phase: str) -> str:
    """One line of always-true guidance for the footer."""

    by_phase = {
        "ready": "Enter a scientific request · /help commands",
        "planning": "Planning with the selected model",
        "request-reviewed": "/approve runs once · /revise or /deny declines",
        "executing": "Executing the approved ChemSmart DAG",
    }
    return by_phase.get(phase, "/help commands") + " · ctrl+c quit"


__all__ = [
    "CommandSpecV1",
    "KeybindSpecV1",
    "COMMANDS",
    "KEYBINDS",
    "command_for",
    "suggest",
    "render_help",
    "footer_hint",
]
