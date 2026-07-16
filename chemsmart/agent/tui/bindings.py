"""Safety bindings that remain fixed regardless of TUI preferences."""

from textual.binding import Binding

BINDINGS = [
    Binding("ctrl+c", "soft_cancel", "Soft cancel", priority=True),
    Binding("ctrl+d", "quit_if_empty", "Quit", priority=True),
    Binding("ctrl+l", "refresh_screen", "Redraw"),
    Binding("escape", "dismiss_overlay", "Dismiss"),
]
