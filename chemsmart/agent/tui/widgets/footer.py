"""Footer widget for phase and status display."""

from __future__ import annotations

import os
from collections import OrderedDict

from rich.spinner import Spinner
from rich.text import Text
from textual.reactive import reactive
from textual.timer import Timer
from textual.widgets import Static

from chemsmart.agent.tui.phase import Phase

_DEFAULT_MODELS = {
    "anthropic": "claude-sonnet-4-6",
    "openai": "gpt-5.4",
}
_ACTIVE_SPINNER_PHASES = frozenset({Phase.PLANNING, Phase.RUNNING})
_SPINNER_LABELS = (
    "thinking…",
    "calculating…",
    "resting…",
    "mapping orbitals…",
    "checking constraints…",
    "tuning the plan…",
    "queuing tools…",
)


class FooterWidget(Static):
    phase: reactive[Phase] = reactive(Phase.IDLE)
    hint: reactive[str] = reactive("Enter to submit • /help for commands")
    entity_status: reactive[str | None] = reactive(None)

    DEFAULT_CSS = """
    FooterWidget {
        height: 1;
        padding: 0 1;
        background: $boost;
        color: $text;
    }
    """

    SPINNER_LABELS = _SPINNER_LABELS
    _SPINNER_INTERVAL_SECONDS = 0.1
    _SPINNER_LABEL_SECONDS = 2.5
    _STATIC_WORKING_LABEL = "Working…"

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._job_counts = OrderedDict(queued=0, running=0, failed=0)
        provider = (os.environ.get("AI_PROVIDER") or "").strip().lower()
        self._provider = provider or "offline"
        self._model = _DEFAULT_MODELS.get(provider, "auto")
        self._draft_tokens = 0
        spinner = Spinner("arc")
        self._spinner_frames = tuple(spinner.frames)
        self._spinner_frame_index = 0
        self._spinner_label_index = 0
        self._spinner_tick_count = 0
        self._spinner_timer: Timer | None = None

    @property
    def spinner_visible(self) -> bool:
        return self.phase in _ACTIVE_SPINNER_PHASES

    @property
    def spinner_text(self) -> str:
        if not self.spinner_visible:
            return ""
        if self._has_reduced_motion():
            return self._STATIC_WORKING_LABEL
        return (
            f"{self._spinner_frames[self._spinner_frame_index]} "
            f"{self.SPINNER_LABELS[self._spinner_label_index]}"
        )

    def on_mount(self) -> None:
        self._spinner_timer = self.set_interval(
            self._SPINNER_INTERVAL_SECONDS,
            self._advance_spinner,
            pause=True,
        )
        self._sync_spinner_timer(reset=True)
        self._refresh_text()

    def watch_phase(self, old: Phase, new: Phase) -> None:
        self._sync_spinner_timer(
            reset=(
                old not in _ACTIVE_SPINNER_PHASES
                and new in _ACTIVE_SPINNER_PHASES
            )
        )
        self._refresh_text()

    def watch_hint(self, old: str, new: str) -> None:
        self._refresh_text()

    def watch_entity_status(
        self,
        old: str | None,
        new: str | None,
    ) -> None:
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

    def _has_reduced_motion(self) -> bool:
        if not self.is_mounted:
            return False
        return bool(getattr(self.app, "plain", False)) or (
            self.app.animation_level == "none"
        )

    def _reset_spinner_state(self) -> None:
        self._spinner_frame_index = 0
        self._spinner_label_index = 0
        self._spinner_tick_count = 0

    def _should_animate_spinner(self) -> bool:
        return self.spinner_visible and not self._has_reduced_motion()

    def _sync_spinner_timer(self, *, reset: bool = False) -> None:
        if reset:
            self._reset_spinner_state()
        if self._spinner_timer is None:
            return
        if self._should_animate_spinner():
            if reset:
                self._spinner_timer.reset()
            self._spinner_timer.resume()
            return
        self._spinner_timer.pause()

    def _advance_spinner(self) -> None:
        if not self._should_animate_spinner():
            return
        self._spinner_tick_count += 1
        self._spinner_frame_index = (self._spinner_frame_index + 1) % len(
            self._spinner_frames
        )
        label_interval = max(
            1,
            round(
                self._SPINNER_LABEL_SECONDS / self._SPINNER_INTERVAL_SECONDS
            ),
        )
        if self._spinner_tick_count % label_interval == 0:
            self._spinner_label_index = (self._spinner_label_index + 1) % len(
                self.SPINNER_LABELS
            )
        self._refresh_text()

    def _append_spinner(self, text: Text) -> None:
        if not self.spinner_visible:
            return
        text.append(" • ", style="dim")
        if self._has_reduced_motion():
            text.append(self._STATIC_WORKING_LABEL, style="dim")
            return
        text.append(
            self._spinner_frames[self._spinner_frame_index], style="accent"
        )
        text.append(" ", style="dim")
        text.append(
            self.SPINNER_LABELS[self._spinner_label_index], style="dim"
        )

    def _refresh_text(self) -> None:
        text = Text()
        text.append(self.phase.label, style="bold")
        self._append_spinner(text)
        if self.hint:
            text.append(" • ", style="dim")
            text.append(self.hint, style="dim")
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
        text.append(" • ", style="dim")
        text.append(self._provider, style="dim")
        text.append("/", style="dim")
        text.append(self._model, style="dim")
        text.append(" • ", style="dim")
        text.append(f"tok {self._draft_tokens}", style="dim")
        if self.entity_status:
            text.append(" | ", style="dim")
            text.append(self.entity_status, style="dim")
        self.update(text)
