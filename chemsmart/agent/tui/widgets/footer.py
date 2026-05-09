"""Footer widget for phase and status display."""

from __future__ import annotations

import os
from collections import OrderedDict

from rich.text import Text
from textual.reactive import reactive
from textual.widgets import Static

from chemsmart.agent.tui.phase import Phase

_DEFAULT_MODELS = {
    "anthropic": "claude-sonnet-4-6",
    "openai": "gpt-5.4",
}


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
        provider = (os.environ.get("AI_PROVIDER") or "").strip().lower()
        self._provider = provider or "offline"
        self._model = _DEFAULT_MODELS.get(provider, "auto")
        self._draft_tokens = 0

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

    def set_provider_model(
        self,
        provider: str | None = None,
        model: str | None = None,
    ) -> None:
        if provider:
            self._provider = provider
        if model:
            self._model = model
        self._refresh_text()

    def update_draft(self, text: str) -> None:
        self._draft_tokens = len([part for part in text.split() if part])
        self._refresh_text()

    def set_job_counts(self, **counts: int) -> None:
        for key in self._job_counts:
            if key in counts:
                self._job_counts[key] = int(counts[key])
        self._refresh_text()

    def reset_job_counts(self) -> None:
        self.set_job_counts(queued=0, running=0, failed=0)

    def _refresh_text(self) -> None:
        text = Text()
        text.append(self.phase.label, style="bold")
        text.append(" • ", style="dim")
        text.append(self._provider, style="dim")
        text.append("/", style="dim")
        text.append(self._model, style="dim")
        text.append(" • ", style="dim")
        text.append(f"tok {self._draft_tokens}", style="dim")
        text.append(" • ", style="dim")
        text.append("jobs ", style="dim")
        text.append(f"q{self._job_counts['queued']} ", style="dim")
        text.append(
            f"r{self._job_counts['running']} ",
            style="accent" if self._job_counts["running"] else "dim",
        )
        text.append(
            f"f{self._job_counts['failed']}",
            style="error" if self._job_counts["failed"] else "dim",
        )
        if self.hint:
            text.append(" • ", style="dim")
            text.append(self.hint, style="dim")
        self.update(text)
