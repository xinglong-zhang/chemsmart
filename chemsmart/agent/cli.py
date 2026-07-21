"""Public CLI facade with lazy optional TUI integration."""

from __future__ import annotations

from chemsmart.agent.cli_commands import agent, configure_tui_launcher
from chemsmart.agent.services.cli_presenters import (
    sanitize_inline_cli_output as sanitize_inline_cli_output,
)


def _launch_tui(*, plain: bool) -> None:
    try:
        from chemsmart.agent.tui import launch_tui
    except ImportError as exc:
        raise ImportError("agent TUI dependencies are unavailable") from exc
    launch_tui(plain=plain)


configure_tui_launcher(_launch_tui)

__all__ = ["agent", "sanitize_inline_cli_output"]
