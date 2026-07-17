"""Keyboard shortcut reference overlay."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import ModalScreen
from textual.widgets import Static


class ShortcutOverlay(ModalScreen[None]):
    BINDINGS = [Binding("escape", "close", "Close", show=False, priority=True)]

    DEFAULT_CSS = """
    ShortcutOverlay { align: center middle; }
    #shortcut-modal {
        width: 76;
        height: auto;
        max-height: 26;
        border: round $accent;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(self, bindings: dict[str, str]) -> None:
        super().__init__()
        self.bindings = bindings

    def compose(self) -> ComposeResult:
        rows = [
            (self.bindings.get("show_shortcuts", "F1"), "shortcut reference"),
            (
                self.bindings.get("toggle_transcript", "Ctrl+O"),
                "compact/full transcript",
            ),
            (self.bindings.get("show_activity", "Ctrl+T"), "tool activity"),
            (
                self.bindings.get("show_calculations", "Ctrl+B"),
                "calculation monitor and logs",
            ),
            (
                self.bindings.get("search_history", "Ctrl+R"),
                "previous request",
            ),
            (
                self.bindings.get("show_project_yaml", "Shift+Tab"),
                "workspace project YAML",
            ),
            ("/ · ↑/↓ · Tab · Enter", "command palette"),
            ("@", "workspace file picker"),
            ("Ctrl+G", "external editor"),
            ("Ctrl+J / Shift+Enter", "newline"),
            ("Click a response", "open mouse-selectable copy view"),
            ("Enter / Space on tool chain", "show or hide completed tools"),
            ("y / s / n / r", "approve once/session, deny, revise"),
            ("Esc", "close overlay and return to request"),
        ]
        width = max(len(key) for key, _ in rows)
        text = "\n".join(f"{key:<{width}}  {label}" for key, label in rows)
        with Vertical(id="shortcut-modal"):
            yield Static(text, id="shortcut-summary")

    def on_mount(self) -> None:
        self.query_one("#shortcut-summary", Static).border_title = (
            "Keyboard shortcuts"
        )

    def action_close(self) -> None:
        self.dismiss(None)
