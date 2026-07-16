"""Scrollable transcript widget."""

from __future__ import annotations

from textual.containers import Vertical, VerticalScroll
from textual.widget import Widget

# Cell types that represent in-flight planning chrome for a single turn —
# removed at the start of a new turn so prior conversational
# UserMessageCell / AgentMessageCell entries stay visible.
_EPHEMERAL_CELL_TYPE_NAMES = frozenset(
    {
        "PlanCell",
        "CriticVerdictCell",
        "DryRunInputCell",
        "RuntimeValidationCell",
        "WorkflowCell",
        "SubmissionPreviewCell",
        "MethodCell",
        "MoleculeCell",
        "GeometryHandoffCell",
        "RunResultCell",
        "ToolCallCell",
        "DecisionTraceCell",
        "CommandInterpretationCell",
        "SynthesisTraceCell",
    }
)


class Transcript(VerticalScroll):
    DEFAULT_CSS = """
    Transcript {
        border: round $panel;
        padding: 0 1;
        overflow-x: hidden;
        overflow-y: auto;
    }

    #cells {
        width: 100%;
        height: auto;
    }
    """

    def compose(self):
        yield Vertical(id="cells")

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.detail_expanded = False
        self.active_turn_id = "turn-0"

    def start_turn(self, turn_id: str) -> None:
        """Select the turn assigned to subsequently mounted cells."""

        self.active_turn_id = str(turn_id or "turn-0")

    def add_cell(self, cell: Widget, *, turn_id: str | None = None) -> None:
        setattr(cell, "_chemsmart_turn_id", turn_id or self.active_turn_id)
        self.query_one("#cells", Vertical).mount(cell)
        setter = getattr(cell, "set_expanded", None)
        if callable(setter) and self.detail_expanded:
            setter(True)
        self.call_after_refresh(self.scroll_end, animate=False)

    def toggle_detail_mode(self) -> bool:
        self.detail_expanded = not self.detail_expanded
        cells = self.query_one("#cells", Vertical)
        for child in cells.children:
            setter = getattr(child, "set_expanded", None)
            if callable(setter):
                setter(self.detail_expanded)
        if self.detail_expanded:
            self.focus()
        return self.detail_expanded

    def clear_cells(self) -> None:
        cells = self.query_one("#cells", Vertical)
        for child in list(cells.children):
            child.remove()

    def clear_turn_chrome(self) -> None:
        cells = self.query_one("#cells", Vertical)
        for child in list(cells.children):
            if type(child).__name__ in _EPHEMERAL_CELL_TYPE_NAMES:
                child.remove()
                continue
            if (
                type(child).__name__ == "AgentMessageCell"
                and getattr(child, "border_title", None) == "Summary"
            ):
                child.remove()
