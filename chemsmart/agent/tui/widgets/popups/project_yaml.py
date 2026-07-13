"""Workspace project YAML preview overlay."""

from __future__ import annotations

from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static


class ProjectYamlOverlay(ModalScreen[None]):
    BINDINGS = [Binding("escape", "close", "Close", show=False)]

    DEFAULT_CSS = """
    ProjectYamlOverlay {
        align: center middle;
    }

    #project-yaml-modal {
        width: 100;
        height: auto;
        max-height: 88%;
        border: round $primary;
        padding: 1 2;
        background: $surface;
    }

    #project-yaml-content {
        margin-top: 1;
    }
    """

    def __init__(
        self,
        *,
        title: str,
        path: Path | None,
        yaml_text: str,
    ) -> None:
        super().__init__()
        self.title = title
        self.path = path
        self.yaml_text = yaml_text

    def compose(self) -> ComposeResult:
        with Vertical(id="project-yaml-modal"):
            header = Static(self._header_text(), id="project-yaml-header")
            header.border_title = self.title
            yield header
            yield Static(
                f"```yaml\n{self.yaml_text.rstrip()}\n```",
                id="project-yaml-content",
            )

    def action_close(self) -> None:
        self.dismiss(None)

    def _header_text(self) -> str:
        if self.path is None:
            return "No workspace project YAML is loaded."
        return f"Path: {self.path}\nEsc returns to the chat."
