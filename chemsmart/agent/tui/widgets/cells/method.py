"""Method recommendation transcript cell."""

from __future__ import annotations

from rich.console import Group
from rich.text import Text
from textual.binding import Binding

from .base import BaseCell


class MethodCell(BaseCell):
    BINDINGS = [Binding("e", "edit_method", "Edit", show=False)]

    def __init__(self, recommendation: dict) -> None:
        self.recommendation = recommendation
        super().__init__(
            self._build_renderable(),
            title="Method recommendation",
            classes="method-cell",
        )

    def action_edit_method(self) -> None:
        handler = getattr(self.app.screen, "edit_method_from_cell", None)
        if handler is not None:
            handler(self.recommendation)

    def _build_renderable(self):
        functional = self.recommendation.get("functional") or "(manual pick)"
        basis = self.recommendation.get("basis") or "(manual pick)"
        match = self.recommendation.get("match") or "no automatic match"
        solvent_model = self.recommendation.get("solvent_model") or "gas phase"
        solvent_id = self.recommendation.get("solvent_id")
        solvent = solvent_model
        if solvent_id:
            solvent = f"{solvent_model} / {solvent_id}"

        lines = [
            Text(f"Method: {functional}/{basis}", style="bold"),
            Text(f"Project match: {match}"),
            Text(f"Solvent: {solvent}"),
        ]
        heavy = self.recommendation.get("heavy_elements") or []
        if heavy:
            heavy_basis = (
                self.recommendation.get("heavy_elements_basis") or "default"
            )
            lines.append(
                Text(
                    "Heavy elements: "
                    f"{', '.join(str(item) for item in heavy)} [{heavy_basis}]"
                )
            )
        rationale = self.recommendation.get("rationale")
        if rationale:
            lines.extend([Text(""), Text(str(rationale))])
        lines.append(Text(""))
        lines.append(Text("Press e to revise method/basis.", style="dim"))
        return Group(*lines)
