"""Crash-tolerant JSONL event storage with hash-chain validation."""

from __future__ import annotations

import fcntl
import json
import os
from pathlib import Path
from threading import RLock
from typing import Any, Iterable

from chemsmart.agent.runtime.events import EventKind, RuntimeEvent


class EventStoreCorruptionError(RuntimeError):
    pass


class RuntimeEventStore:
    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = RLock()

    def append(
        self,
        *,
        session_id: str,
        turn_id: str,
        kind: EventKind,
        payload: dict[str, Any] | None = None,
        idempotency_key: str = "",
    ) -> RuntimeEvent:
        with self._lock, self.path.open("a+", encoding="utf-8") as handle:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                handle.seek(0)
                events = self._parse_lines(handle)
                self._validate(events)
                if idempotency_key:
                    for event in events:
                        if event.idempotency_key == idempotency_key:
                            return event
                event = RuntimeEvent.create(
                    sequence=len(events) + 1,
                    session_id=session_id,
                    turn_id=turn_id,
                    kind=kind,
                    payload=dict(payload or {}),
                    previous_hash=events[-1].event_hash if events else "",
                    idempotency_key=idempotency_key,
                )
                handle.seek(0, os.SEEK_END)
                handle.write(event.model_dump_json() + "\n")
                handle.flush()
                os.fsync(handle.fileno())
                return event
            finally:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    def load(self) -> list[RuntimeEvent]:
        if not self.path.exists():
            return []
        with self._lock, self.path.open(encoding="utf-8") as handle:
            events = self._parse_lines(handle)
        self._validate(events)
        return events

    def write_snapshot(
        self, payload: dict[str, Any], path: str | Path
    ) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)

    @staticmethod
    def _parse_lines(lines: Iterable[str]) -> list[RuntimeEvent]:
        events: list[RuntimeEvent] = []
        for line_number, line in enumerate(lines, start=1):
            if not line.strip():
                continue
            try:
                events.append(RuntimeEvent.model_validate_json(line))
            except Exception as exc:
                raise EventStoreCorruptionError(
                    f"invalid runtime event at line {line_number}: {exc}"
                ) from exc
        return events

    @staticmethod
    def _validate(events: list[RuntimeEvent]) -> None:
        previous_hash = ""
        for expected_sequence, event in enumerate(events, start=1):
            if event.sequence != expected_sequence:
                raise EventStoreCorruptionError(
                    "runtime event sequence is not contiguous"
                )
            if event.previous_hash != previous_hash:
                raise EventStoreCorruptionError(
                    "runtime event hash chain is broken"
                )
            if not event.verify_hash():
                raise EventStoreCorruptionError(
                    f"runtime event {event.event_id} has an invalid hash"
                )
            previous_hash = event.event_hash


__all__ = ["EventStoreCorruptionError", "RuntimeEventStore"]
