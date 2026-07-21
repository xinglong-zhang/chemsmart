"""Responsive footer backed by observable TUI state."""

from __future__ import annotations

from rich.spinner import Spinner
from rich.text import Text
from textual import events
from textual.reactive import reactive
from textual.timer import Timer
from textual.widgets import Static

from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.state import TuiState, TuiStateReducer

_ACTIVE_SPINNER_PHASES = frozenset(
    {
        Phase.PLANNING,
        Phase.TOOL_RUNNING,
        Phase.VALIDATING,
        Phase.RUNNING,
        Phase.EXECUTING,
        Phase.SUBMITTING,
    }
)
_SPINNER_LABELS = ("working",)


class FooterWidget(Static):
    phase: reactive[Phase] = reactive(Phase.IDLE)
    hint: reactive[str] = reactive("Enter to submit · F1 shortcuts")
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
    _STATIC_WORKING_LABEL = "Working…"

    def __init__(
        self,
        *args,
        state: TuiState | None = None,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.reducer = TuiStateReducer(state)
        self.state = self.reducer.state
        spinner = Spinner("arc")
        self._spinner_frames = tuple(spinner.frames)
        self._spinner_frame_index = 0
        self._spinner_timer: Timer | None = None
        self._draft_estimate: int | None = None

    @property
    def spinner_visible(self) -> bool:
        return self.phase in _ACTIVE_SPINNER_PHASES

    @property
    def spinner_text(self) -> str:
        if not self.spinner_visible:
            return ""
        if self._has_reduced_motion():
            return self._STATIC_WORKING_LABEL
        return f"{self._spinner_frames[self._spinner_frame_index]} working"

    def on_mount(self) -> None:
        self._spinner_timer = self.set_interval(
            self._SPINNER_INTERVAL_SECONDS,
            self._advance_spinner,
            pause=True,
        )
        self._sync_spinner_timer(reset=True)
        self._refresh_text()

    def on_resize(self, _event: events.Resize) -> None:
        self._refresh_text()

    def watch_phase(self, old: Phase, new: Phase) -> None:
        self.reducer.set_phase(new)
        self._sync_spinner_timer(
            reset=old not in _ACTIVE_SPINNER_PHASES
            and new in _ACTIVE_SPINNER_PHASES
        )
        self._refresh_text()

    def watch_hint(self, _old: str, new: str) -> None:
        self.reducer.update(operation=new)
        self._refresh_text()

    def watch_entity_status(self, _old: str | None, _new: str | None) -> None:
        self._refresh_text()

    def set_phase(self, phase: Phase) -> None:
        self.phase = phase

    def set_hint(self, hint: str) -> None:
        self.hint = hint

    def set_provider_model(
        self,
        provider: str | None = None,
        model: str | None = None,
        project: str | None = None,
    ) -> None:
        if provider:
            self.reducer.update(provider=provider)
        if model:
            self.reducer.update(model=model)
        if project is not None:
            self.reducer.update(project=project)
        self._refresh_text()

    def set_yaml_status(
        self, *, loaded: bool, label: str | None = None
    ) -> None:
        self.reducer.update(
            yaml_loaded=bool(loaded),
            yaml_label=label or ("YAML OK" if loaded else "YAML MISSING"),
        )
        self._refresh_text()

    def set_server(self, server: str | None) -> None:
        self.reducer.update(server=server or "local/default")
        self._refresh_text()

    def set_permission(self, mode: str, *, yolo: bool) -> None:
        self.reducer.update(permission_mode=mode, yolo=yolo)
        self._refresh_text()

    def set_tool_progress(
        self,
        tool: str,
        *,
        step: int | None = None,
        total: int | None = None,
    ) -> None:
        self.reducer.set_tool(tool, step=step, total=total)
        self._refresh_text()

    def set_queued_prompt(self, queued: bool) -> None:
        self.reducer.update(queued_prompt=queued)
        self._refresh_text()

    def set_usage(
        self,
        *,
        input_tokens: int | None,
        output_tokens: int | None,
    ) -> None:
        self.reducer.update(
            usage_input_tokens=input_tokens,
            usage_output_tokens=output_tokens,
        )
        self._refresh_text()

    def update_draft(self, text: str) -> None:
        stripped = text.strip()
        self._draft_estimate = (
            max(1, (len(stripped) + 3) // 4) if stripped else None
        )
        self._refresh_text()

    def set_job_counts(self, **counts: int) -> None:
        self.reducer.set_jobs(**counts)
        self._refresh_text()

    def reset_job_counts(self) -> None:
        self.set_job_counts(queued=0, running=0, failed=0)

    def _has_reduced_motion(self) -> bool:
        if not self.is_mounted:
            return False
        return bool(getattr(self.app, "plain", False)) or (
            self.app.animation_level == "none"
        )

    def _sync_spinner_timer(self, *, reset: bool = False) -> None:
        if reset:
            self._spinner_frame_index = 0
        if self._spinner_timer is None:
            return
        if self.spinner_visible and not self._has_reduced_motion():
            self._spinner_timer.resume()
        else:
            self._spinner_timer.pause()

    def _advance_spinner(self) -> None:
        if not self.spinner_visible or self._has_reduced_motion():
            return
        self._spinner_frame_index = (self._spinner_frame_index + 1) % len(
            self._spinner_frames
        )
        self._refresh_text()

    @staticmethod
    def _segment(label: str, value: str, style: str = "dim") -> Text:
        return Text.assemble((label, "dim"), (value, style))

    def _refresh_text(self) -> None:
        width = self.size.width or 120
        text = Text(self.phase.label, style="bold")
        if self.spinner_visible:
            text.append(" ", style="dim")
            if self._has_reduced_motion():
                text.append(self._STATIC_WORKING_LABEL, style="dim")
            else:
                text.append(
                    self._spinner_frames[self._spinner_frame_index],
                    style="accent",
                )
                text.append(" working", style="dim")

        segments: list[tuple[int, Text]] = []
        operation = self.hint or self.state.tool_progress
        if operation:
            segments.append((0, Text(operation, style="dim")))
        yaml_style = "success" if self.state.yaml_loaded else "error"
        project = self.state.project or "none"
        segments.append(
            (
                1,
                Text.assemble(
                    ("project ", "dim"),
                    (project, "accent"),
                    (" · ", "dim"),
                    (self.state.yaml_label, yaml_style),
                    (" (S-Tab)", "dim"),
                ),
            )
        )
        segments.append(
            (
                2,
                self._segment(
                    "ask:", f"{self.state.provider}/{self.state.model}"
                ),
            )
        )
        counts = self.state.job_counts
        segments.append(
            (
                3,
                Text(
                    f"jobs q{counts['queued']} r{counts['running']} "
                    f"f{counts['failed']}",
                    style="dim",
                ),
            )
        )
        if self.state.queued_prompt:
            segments.append((1, Text("next queued", style="warning")))
        segments.append(
            (
                4,
                Text(
                    f"server {self.state.server} · permission "
                    f"{self.state.permission_mode}"
                    f"{'/yolo' if self.state.yolo else ''}",
                    style="dim",
                ),
            )
        )
        if (
            self.state.usage_input_tokens is not None
            or self.state.usage_output_tokens is not None
        ):
            segments.append(
                (
                    4,
                    Text(
                        "usage "
                        f"{self.state.usage_input_tokens or 0}/"
                        f"{self.state.usage_output_tokens or 0}",
                        style="dim",
                    ),
                )
            )
        elif self._draft_estimate is not None:
            segments.append(
                (4, Text(f"draft ~{self._draft_estimate}", style="dim"))
            )
        if self.entity_status:
            segments.append((5, Text(self.entity_status, style="dim")))

        # Preserve high-priority state on narrow terminals and progressively
        # reveal operational detail as space becomes available.
        for priority, segment in segments:
            minimum_width = (0, 52, 78, 102, 126, 150)[priority]
            if width < minimum_width:
                continue
            text.append(" · ", style="dim")
            text.append_text(segment)
        self.update(text)
