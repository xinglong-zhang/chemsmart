"""Watchdog-backed tailer for decision_log.jsonl."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

from watchdog.events import FileSystemEvent, FileSystemEventHandler
from watchdog.observers import Observer


class _TailEventHandler(FileSystemEventHandler):
    def __init__(self, tailer: "LogTailer") -> None:
        self.tailer = tailer

    def on_created(self, event: FileSystemEvent) -> None:
        self._maybe_read(event)

    def on_modified(self, event: FileSystemEvent) -> None:
        self._maybe_read(event)

    def _maybe_read(self, event: FileSystemEvent) -> None:
        if Path(event.src_path) == self.tailer.path:
            self.tailer.read_available()


class LogTailer:
    def __init__(
        self,
        path: str | Path,
        on_entry: Callable[[dict], None],
    ) -> None:
        self.path = Path(path)
        self.on_entry = on_entry
        self._offset = 0
        self._observer: Observer | None = None
        self._handler = _TailEventHandler(self)

    def start(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.read_available()
        observer = Observer()
        observer.schedule(
            self._handler, str(self.path.parent), recursive=False
        )
        observer.start()
        self._observer = observer

    def stop(self) -> None:
        if self._observer is None:
            return
        self._observer.stop()
        self._observer.join(timeout=1)
        self._observer = None

    def read_available(self) -> None:
        if not self.path.exists():
            return
        size = self.path.stat().st_size
        if size < self._offset:
            self._offset = 0
        with self.path.open(encoding="utf-8") as handle:
            handle.seek(self._offset)
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                self.on_entry(json.loads(line))
            self._offset = handle.tell()
