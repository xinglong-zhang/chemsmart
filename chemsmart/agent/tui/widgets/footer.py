"""Footer widget for phase and hint display."""

from __future__ import annotations

from textual.reactive import reactive
from textual.widgets import Static

from chemsmart.agent.tui.phase import Phase


class FooterWidget(Static):
    phase: reactive[Phase] = reactive(Phase.IDLE)
    hint: reactive[str] = reactive("Enter to submit • /help for commands")

    DEFAULT_CSS = """
    FooterWidget {
        height: 1;
        padding: 0 1;
        background: $boost;
        color: $text;
    }
    """

    def on_mount(self) -> None:
        self._refresh_text()

    def watch_phase(self, old: Phase, new: Phase) -> None:
        self._refresh_text()

    def watch_hint(self, old: str, new: str) -> None:
        self._refresh_text()

    def set_phase(self, phase: Phase) -> None:
        self.phase = phase

    def set_hint(self, hint: str) -> None:
        self.hint = hint

    def _refresh_text(self) -> None:
        self.update(f"{self.phase.label}   {self.hint}")
