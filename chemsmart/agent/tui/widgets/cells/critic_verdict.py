"""Critic verdict transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text

from chemsmart.agent.core import CriticVerdict

from .base import BaseCell


class CriticVerdictCell(BaseCell):
    def __init__(self, verdict: CriticVerdict) -> None:
        lines = [
            Text(f"Verdict: {verdict.verdict}", style="bold"),
            Text(f"Confidence: {verdict.confidence:.2f}"),
        ]
        if verdict.issues:
            lines.append(Text("Issues:", style="bold"))
            lines.extend(Text(f"- {issue}") for issue in verdict.issues)
        if verdict.rationale:
            lines.extend([Text(""), Text(verdict.rationale)])
        super().__init__(
            Group(*lines),
            title="Critic verdict",
            classes="critic-cell",
        )
