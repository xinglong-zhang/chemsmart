"""Key bindings for the Phase 1 agent TUI."""

from textual.binding import Binding

BINDINGS = [
    Binding("ctrl+c", "soft_cancel", "Soft cancel", priority=True),
    Binding("ctrl+d", "quit_if_empty", "Quit", priority=True),
    Binding("ctrl+l", "refresh_screen", "Redraw"),
    Binding("shift+tab", "show_project_yaml", "Project YAML", priority=True),
    Binding("escape", "dismiss_overlay", "Dismiss"),
]
