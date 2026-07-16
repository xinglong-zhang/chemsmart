"""Bottom composer widget."""

from __future__ import annotations

import os
import shlex
import subprocess
from time import monotonic
from tempfile import NamedTemporaryFile

from textual import events
from textual.binding import Binding
from textual.message import Message
from textual.widgets import TextArea


class Composer(TextArea):
    DEFAULT_CSS = """
    Composer {
        height: 5;
        min-height: 5;
        border: round $primary;
    }
    """

    BINDINGS = [
        Binding("enter", "submit", show=False, priority=True),
        Binding("tab", "context_tab", show=False, priority=True),
        Binding("ctrl+r", "history", show=False, priority=True),
        Binding("ctrl+j", "insert_newline", show=False, priority=True),
        Binding("shift+enter", "insert_newline", show=False, priority=True),
        Binding("ctrl+g", "external_editor", show=False, priority=True),
        Binding("@", "file_popup", show=False, priority=True),
    ]

    class Submitted(Message):
        def __init__(self, composer: "Composer", text: str) -> None:
            super().__init__()
            self.composer = composer
            self.text = text

    def __init__(self) -> None:
        super().__init__(
            text="",
            soft_wrap=True,
            show_line_numbers=False,
            id="composer",
        )
        self.border_title = "Request"
        self._large_paste_chunks: list[tuple[str, str]] = []
        self._submitting = False
        self._paste_guard_text = ""
        self._paste_guard_deadline = 0.0

    async def _on_paste(self, event: events.Paste) -> None:
        event.stop()
        event.prevent_default()
        if len(event.text) <= 10000:
            if result := self._replace_via_keyboard(
                event.text,
                *self.selection,
            ):
                self.move_cursor(result.end_location)
            self._arm_paste_guard(event.text)
            return
        placeholder = f"[Pasted {len(event.text)} chars]"
        self._large_paste_chunks.append((placeholder, event.text))
        if result := self._replace_via_keyboard(placeholder, *self.selection):
            self.move_cursor(result.end_location)
        self._arm_paste_guard(event.text)

    async def _on_key(self, event: events.Key) -> None:
        if self._consume_duplicate_paste_key(event):
            return
        if self._matches_configured_key(event, "show_project_yaml"):
            event.stop()
            event.prevent_default()
            handler = getattr(self.app, "action_show_project_yaml", None)
            if callable(handler):
                handler()
            return
        palette = self._slash_palette()
        if palette is not None and palette.is_open:
            if event.key == "up":
                event.stop()
                event.prevent_default()
                palette.move_highlight(-1)
                return
            if event.key == "down":
                event.stop()
                event.prevent_default()
                palette.move_highlight(1)
                return
            if event.key == "escape":
                event.stop()
                event.prevent_default()
                palette.hide()
                return
        await super()._on_key(event)

    def _matches_configured_key(
        self,
        event: events.Key,
        action: str,
    ) -> bool:
        config = getattr(self.app, "tui_config", None)
        bindings = getattr(config, "keybindings", {})
        configured = str(bindings.get(action) or "").strip().lower()
        if not configured:
            return False
        event_key = str(event.key or "").strip().lower()
        aliases = {configured}
        if configured == "shift+tab":
            # Textual/terminal backends may normalize Shift+Tab either way.
            aliases.add("backtab")
        return event_key in aliases

    def action_submit(self) -> None:
        if self._submitting:
            return
        screen = getattr(self.app, "screen", None)
        accept = getattr(screen, "accept_slash_palette", None)
        if callable(accept):
            selected = accept()
            if selected is False:
                return
            if selected:
                self.load_text(selected)
        text = self.resolve_text().strip()
        if not text:
            return
        self._submitting = True
        self.post_message(self.Submitted(self, text))

    def action_insert_newline(self) -> None:
        self.insert("\n")

    def action_context_tab(self) -> None:
        screen = getattr(self.app, "screen", None)
        handler = getattr(screen, "action_context_tab", None)
        if callable(handler):
            handler()

    def action_history(self) -> None:
        screen = getattr(self.app, "screen", None)
        handler = getattr(screen, "action_search_history", None)
        if callable(handler):
            handler()

    def action_external_editor(self) -> None:
        editor = os.environ.get("EDITOR", "vi").strip() or "vi"
        with NamedTemporaryFile(
            "w+",
            suffix=".md",
            prefix="chemsmart-agent-",
            delete=False,
            encoding="utf-8",
        ) as handle:
            handle.write(self.resolve_text())
            handle.flush()
            path = handle.name
        subprocess.run([*shlex.split(editor), path], check=False)
        with open(path, encoding="utf-8") as handle:
            self.load_text(handle.read())

    def action_file_popup(self) -> None:
        handler = getattr(self.app.screen, "open_file_picker", None)
        if handler is not None:
            handler()

    def insert_file_reference(self, value: str) -> None:
        self.insert(value)

    def resolve_text(self) -> str:
        text = self.text
        for placeholder, payload in self._large_paste_chunks:
            if placeholder not in text:
                continue
            text = text.replace(placeholder, payload, 1)
        return text

    def _arm_paste_guard(self, text: str) -> None:
        # Some macOS terminal/Textual combinations deliver Cmd+V as both a
        # bracketed paste event and an immediate printable key stream. Consume
        # only the exact duplicate stream, and only for a short window.
        self._paste_guard_text = text
        self._paste_guard_deadline = monotonic() + 0.5

    def _consume_duplicate_paste_key(self, event: events.Key) -> bool:
        if not self._paste_guard_text:
            return False
        if monotonic() > self._paste_guard_deadline:
            self._paste_guard_text = ""
            return False
        character = event.character
        if not character or not event.is_printable:
            return False
        if not self._paste_guard_text.startswith(character):
            self._paste_guard_text = ""
            return False
        self._paste_guard_text = self._paste_guard_text[len(character) :]
        event.stop()
        event.prevent_default()
        return True

    def clear_text(self) -> None:
        self._large_paste_chunks.clear()
        self._submitting = False
        self._paste_guard_text = ""
        self.load_text("")

    def load_text(self, text: str) -> None:
        self._large_paste_chunks.clear()
        self._submitting = False
        self._paste_guard_text = ""
        super().load_text(text)

    def _slash_palette(self):
        screen = getattr(self.app, "screen", None)
        if screen is None:
            return None
        try:
            from chemsmart.agent.tui.widgets.slash_palette import (
                SlashCommandPalette,
            )

            return screen.query_one(SlashCommandPalette)
        except Exception:
            return None
