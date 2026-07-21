"""Workspace project YAML selector and preview overlay."""

from __future__ import annotations

from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static


class ProjectYamlOverlay(ModalScreen[Path | None]):
    BINDINGS = [
        Binding("escape", "close", "Close", show=False),
        Binding("up", "previous_project", "Previous project", show=False),
        Binding("down", "next_project", "Next project", show=False),
        Binding(
            "shift+tab", "previous_project", "Previous project", show=False
        ),
        Binding("tab", "next_project", "Next project", show=False),
        Binding("enter", "select_project", "Use project", show=False),
    ]

    DEFAULT_CSS = """
    ProjectYamlOverlay {
        align: center middle;
    }

    #project-yaml-modal {
        width: 96%;
        max-width: 100;
        height: auto;
        max-height: 88%;
        border: round $primary;
        padding: 1 2;
        background: $surface;
    }

    #project-yaml-selector {
        margin-top: 1;
        color: $text-muted;
    }

    #project-yaml-content {
        margin-top: 1;
    }
    """

    def __init__(
        self,
        *,
        title: str,
        candidates: tuple[Path, ...] = (),
        selected_path: Path | None = None,
        path: Path | None = None,
        yaml_text: str = "",
    ) -> None:
        super().__init__()
        self.title = title
        normalized = tuple(
            candidate.expanduser().resolve() for candidate in candidates
        )
        fallback_path = selected_path or path
        if not normalized and fallback_path is not None:
            normalized = (fallback_path.expanduser().resolve(),)
        self.candidates = normalized
        resolved_selection = (
            fallback_path.expanduser().resolve()
            if fallback_path is not None
            else None
        )
        self._index = (
            normalized.index(resolved_selection)
            if resolved_selection in normalized
            else 0
        )
        self.path = normalized[self._index] if normalized else None
        self.yaml_text = yaml_text

    def compose(self) -> ComposeResult:
        with Vertical(id="project-yaml-modal"):
            header = Static(self._header_text(), id="project-yaml-header")
            header.border_title = self.title
            yield header
            yield Static(
                self._selector_text(),
                id="project-yaml-selector",
                markup=False,
            )
            yield Static(
                self._content_text(),
                id="project-yaml-content",
                markup=False,
            )

    def on_mount(self) -> None:
        self._render_current()

    def action_close(self) -> None:
        self.dismiss(None)

    def action_previous_project(self) -> None:
        self._move(-1)

    def action_next_project(self) -> None:
        self._move(1)

    def action_select_project(self) -> None:
        if self.path is not None:
            self.dismiss(self.path)

    def _move(self, offset: int) -> None:
        if len(self.candidates) < 2:
            return
        self._index = (self._index + offset) % len(self.candidates)
        self.path = self.candidates[self._index]
        self.yaml_text = ""
        self._render_current()

    def _render_current(self) -> None:
        if self.path is not None and not self.yaml_text:
            try:
                self.yaml_text = self.path.read_text(encoding="utf-8")
            except OSError as exc:
                self.yaml_text = f"# Failed to read {self.path}: {exc}\n"
        self.query_one("#project-yaml-header", Static).update(
            self._header_text()
        )
        self.query_one("#project-yaml-selector", Static).update(
            self._selector_text()
        )
        self.query_one("#project-yaml-content", Static).update(
            self._content_text()
        )

    def _header_text(self) -> str:
        if self.path is None:
            return (
                "No workspace project YAML is available. Esc returns to chat."
            )
        return (
            f"Path: {self.path}\n"
            "Tab/Shift+Tab or arrows switch; Enter activates; Esc returns."
        )

    def _selector_text(self) -> str:
        if not self.candidates:
            return "Projects: none"
        lines = [f"Projects ({self._index + 1}/{len(self.candidates)}):"]
        for index, candidate in enumerate(self.candidates):
            marker = ">" if index == self._index else " "
            lines.append(f"{marker} {candidate.parent.name}:{candidate.stem}")
        return "\n".join(lines)

    def _content_text(self) -> str:
        if self.path is None:
            return (
                "# YAML MISSING\n"
                "# Build one with /init, then save it with /write-project.\n"
            )
        return self.yaml_text.rstrip()
