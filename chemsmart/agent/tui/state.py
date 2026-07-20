"""Observable state for the agent TUI presentation layer."""

from __future__ import annotations

from dataclasses import dataclass, field

from chemsmart.agent.tui.phase import Phase


@dataclass(slots=True)
class TuiState:
    """Small, provider-neutral state projected into the footer and overlays."""

    phase: Phase = Phase.IDLE
    operation: str = "Ready"
    provider: str = "offline"
    model: str = "auto"
    project: str = ""
    server: str = "local/default"
    yaml_loaded: bool = False
    yaml_label: str = "YAML MISSING"
    permission_mode: str = "driving"
    yolo: bool = False
    queued_prompt: bool = False
    current_tool: str = ""
    current_step: int | None = None
    total_steps: int | None = None
    usage_input_tokens: int | None = None
    usage_output_tokens: int | None = None
    job_counts: dict[str, int] = field(
        default_factory=lambda: {"queued": 0, "running": 0, "failed": 0}
    )

    def set_phase(self, phase: Phase, operation: str | None = None) -> None:
        self.phase = phase
        if operation is not None:
            self.operation = operation

    def set_tool(
        self,
        tool: str,
        *,
        step: int | None = None,
        total: int | None = None,
    ) -> None:
        self.current_tool = tool
        self.current_step = step
        self.total_steps = total

    @property
    def tool_progress(self) -> str:
        if not self.current_tool:
            return ""
        if self.current_step is None or self.total_steps is None:
            return self.current_tool
        return (
            f"tool {self.current_step}/{self.total_steps} {self.current_tool}"
        )


class TuiStateReducer:
    """Apply presentation events to one shared, observable state object."""

    def __init__(self, state: TuiState | None = None) -> None:
        self.state = state or TuiState()

    def update(self, **changes: object) -> TuiState:
        for name, value in changes.items():
            if not hasattr(self.state, name):
                raise KeyError(f"Unknown TUI state field: {name}")
            setattr(self.state, name, value)
        return self.state

    def set_phase(
        self, phase: Phase, operation: str | None = None
    ) -> TuiState:
        self.state.set_phase(phase, operation)
        return self.state

    def set_tool(
        self,
        tool: str,
        *,
        step: int | None = None,
        total: int | None = None,
    ) -> TuiState:
        self.state.set_tool(tool, step=step, total=total)
        return self.state

    def set_jobs(self, **counts: int) -> TuiState:
        for key in self.state.job_counts:
            if key in counts:
                self.state.job_counts[key] = int(counts[key])
        return self.state
