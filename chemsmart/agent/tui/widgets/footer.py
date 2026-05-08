"""Footer widget for phase and status display."""

from __future__ import annotations

from collections import OrderedDict

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

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._job_counts = OrderedDict(queued=0, running=0, failed=0)

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

    def set_job_counts(self, **counts: int) -> None:
        for key in self._job_counts:
            if key in counts:
                self._job_counts[key] = int(counts[key])
        self._refresh_text()

    def reset_job_counts(self) -> None:
        self.set_job_counts(queued=0, running=0, failed=0)

    def _refresh_text(self) -> None:
        jobs = (
            f"jobs: {self._job_counts['queued']}q "
            f"{self._job_counts['running']}r "
            f"{self._job_counts['failed']}f"
        )
        self.update(f"{self.phase.label}   {jobs}   {self.hint}")
