"""Semantic color language of the ChemSmart terminal.

One restrained palette with a single cyan accent. State encoding follows one
rule everywhere: finished work goes muted, active work stays plain text,
pending human decision is warning, failure is error, validated success is
success. Depth comes from background steps, never from extra hues.
"""

from __future__ import annotations

from textual.theme import Theme

_VARIABLES = {
    # Transcript row states (the one rule, named once).
    "row-running": "$text",
    "row-finished": "$text-muted",
    "row-failed": "$text-error",
    "row-pending-decision": "$text-warning",
    # Panel chrome.
    "panel-border": "$primary 40%",
    "chip-accent": "$primary",
}

CHEMSMART_THEME = Theme(
    name="chemsmart",
    primary="#22b8cf",
    secondary="#5c7cfa",
    accent="#22b8cf",
    warning="#e8a33d",
    error="#e05252",
    success="#51cf66",
    foreground="#e6e6e6",
    background="#101216",
    surface="#181b20",
    panel="#1f232a",
    dark=True,
    variables=_VARIABLES,
)

__all__ = ["CHEMSMART_THEME"]
