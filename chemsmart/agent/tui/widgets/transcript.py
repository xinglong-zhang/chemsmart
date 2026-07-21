"""Scrollable transcript widget."""

from __future__ import annotations

from textual.containers import Vertical, VerticalScroll
from textual.widget import Widget

from chemsmart.agent.tui.widgets.cells.tool_chain import ToolChainToggleCell

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

_TOOL_CHAIN_CELL_TYPE_NAMES = _EPHEMERAL_CELL_TYPE_NAMES - {
    "RunResultCell",
}


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
        self._turn_tool_chains: dict[
            str, tuple[ToolChainToggleCell, tuple[Widget, ...]]
        ] = {}

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
        for turn_id in tuple(self._turn_tool_chains):
            self.set_tool_chain_expanded(turn_id, self.detail_expanded)
        if self.detail_expanded:
            self.focus()
        return self.detail_expanded

    def clear_cells(self) -> None:
        self._turn_tool_chains.clear()
        cells = self.query_one("#cells", Vertical)
        for child in list(cells.children):
            child.remove()

    def clear_turn_chrome(self) -> None:
        cells = self.query_one("#cells", Vertical)
        preserved = {
            cell
            for _toggle, chain_cells in self._turn_tool_chains.values()
            for cell in chain_cells
        }
        for child in list(cells.children):
            if child in preserved:
                continue
            if type(child).__name__ in _EPHEMERAL_CELL_TYPE_NAMES:
                child.remove()
                continue
            if (
                type(child).__name__ == "AgentMessageCell"
                and getattr(child, "border_title", None) == "Summary"
            ):
                child.remove()

    def collapse_tool_chain(self, turn_id: str) -> bool:
        """Hide completed activity while keeping one final deliverable visible."""

        normalized = str(turn_id or self.active_turn_id)
        container = self.query_one("#cells", Vertical)
        turn_cells = [
            child
            for child in container.children
            if getattr(child, "_chemsmart_turn_id", None) == normalized
        ]
        final_cell = _final_turn_deliverable(turn_cells)
        if final_cell is None:
            return False
        chain_cells = tuple(
            child
            for child in turn_cells
            if child is not final_cell
            and (
                type(child).__name__ in _TOOL_CHAIN_CELL_TYPE_NAMES
                or type(child).__name__
                in {"FinalAnswerCell", "AgentMessageCell", "ErrorCell"}
            )
        )
        if not chain_cells:
            return False
        existing = self._turn_tool_chains.get(normalized)
        if existing is not None:
            toggle, _previous_cells = existing
            toggle.set_tool_count(len(chain_cells))
            self._turn_tool_chains[normalized] = (toggle, chain_cells)
            final_cell.display = True
            self.set_tool_chain_expanded(normalized, False)
            return True
        toggle = ToolChainToggleCell(
            turn_id=normalized,
            tool_count=len(chain_cells),
        )
        setattr(toggle, "_chemsmart_turn_id", normalized)
        container.mount(toggle, before=final_cell)
        self._turn_tool_chains[normalized] = (toggle, chain_cells)
        self.set_tool_chain_expanded(normalized, False)
        return True

    def set_tool_chain_expanded(self, turn_id: str, expanded: bool) -> bool:
        entry = self._turn_tool_chains.get(str(turn_id))
        if entry is None:
            return False
        toggle, chain_cells = entry
        for cell in chain_cells:
            cell.display = bool(expanded)
        toggle.set_expanded(expanded)
        return True

    def on_tool_chain_toggle_cell_toggled(
        self, event: ToolChainToggleCell.Toggled
    ) -> None:
        self.set_tool_chain_expanded(event.turn_id, event.expanded)


def _final_turn_deliverable(turn_cells: list[Widget]) -> Widget | None:
    explicit = next(
        (
            cell
            for cell in reversed(turn_cells)
            if getattr(cell, "_chemsmart_final_deliverable", False)
        ),
        None,
    )
    if explicit is not None:
        return explicit

    last_answer = next(
        (
            cell
            for cell in reversed(turn_cells)
            if type(cell).__name__ == "FinalAnswerCell"
        ),
        None,
    )
    last_error = next(
        (
            cell
            for cell in reversed(turn_cells)
            if type(cell).__name__ == "ErrorCell"
        ),
        None,
    )
    if last_error is not None and (
        last_answer is None
        or turn_cells.index(last_error) > turn_cells.index(last_answer)
    ):
        return last_error
    if last_answer is not None:
        return last_answer

    conversational = next(
        (
            cell
            for cell in reversed(turn_cells)
            if type(cell).__name__ == "AgentMessageCell"
            and getattr(cell, "border_title", None) != "Summary"
        ),
        None,
    )
    if conversational is not None:
        return conversational
    return next(
        (
            cell
            for cell in reversed(turn_cells)
            if type(cell).__name__ in {"AgentMessageCell", "ErrorCell"}
        ),
        None,
    )
