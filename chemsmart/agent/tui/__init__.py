"""Textual TUI for ``chemsmart agent``."""

from __future__ import annotations


def launch_tui(*, plain: bool = False, session_root=None) -> None:
    from .app import launch_tui as _launch_tui

    _launch_tui(plain=plain, session_root=session_root)


__all__ = ["launch_tui"]
