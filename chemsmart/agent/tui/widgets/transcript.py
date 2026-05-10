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

    def add_cell(self, cell: Widget) -> None:
        self.query_one("#cells", Vertical).mount(cell)
        self.call_after_refresh(self.scroll_end, animate=False)

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
