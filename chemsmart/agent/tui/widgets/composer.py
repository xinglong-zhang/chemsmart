"""Bottom composer widget."""

from __future__ import annotations

import os
import shlex
import subprocess
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

    async def _on_paste(self, event: events.Paste) -> None:
        if len(event.text) <= 10000:
            await super()._on_paste(event)
            return
        placeholder = f"[Pasted {len(event.text)} chars]"
        self._large_paste_chunks.append((placeholder, event.text))
        if result := self._replace_via_keyboard(placeholder, *self.selection):
            self.move_cursor(result.end_location)

    def action_submit(self) -> None:
        if self._submitting:
            return
        text = self.resolve_text().strip()
        if not text:
            return
        self._submitting = True
        self.post_message(self.Submitted(self, text))

    def action_insert_newline(self) -> None:
        self.insert("\n")

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

    def clear_text(self) -> None:
        self._large_paste_chunks.clear()
        self._submitting = False
        self.load_text("")

    def load_text(self, text: str) -> None:
        self._large_paste_chunks.clear()
        self._submitting = False
        super().load_text(text)
