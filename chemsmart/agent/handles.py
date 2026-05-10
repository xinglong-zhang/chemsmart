from __future__ import annotations

import json
import re
import secrets
from copy import deepcopy
from pathlib import Path
from typing import Any

_HANDLE_PREFIXES = {
    "mol": "mol",
    "gset": "gset",
    "oset": "oset",
    "job": "job",
    "dryrun": "dryrun",
    "runtime": "runtime",
    "runresult": "runresult",
    "geom": "geom",
    "submit": "submit",
    "recmethod": "recmethod",
}
_HANDLE_PREFIX_PATTERN = "|".join(
    sorted(_HANDLE_PREFIXES.values(), key=len, reverse=True)
)
HANDLE_ID_RE = re.compile(rf"^(?:{_HANDLE_PREFIX_PATTERN})_[0-9a-f]{{4,}}$")


class HandleStore:
    """Session-scoped artifact handle store.

    Invariants:
    - Handle ids are unique within one store instance.
    - Every persisted record contains only JSON-safe data.
    - `get()` returns only objects inserted in the current process.
    - `get_summary()` works for current-process and reloaded handles.
    """

    def __init__(self, session_dir: str | Path) -> None:
        self.session_dir = Path(session_dir)
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.path = self.session_dir / "handles.jsonl"
        self._objects: dict[str, Any] = {}
        self._summaries: dict[str, dict[str, Any]] = {}
        self._kinds: dict[str, str] = {}
        self._load()

    def put(self, kind: str, obj: Any, summary: dict[str, Any]) -> str:
        if kind not in _HANDLE_PREFIXES:
            known = ", ".join(sorted(_HANDLE_PREFIXES))
            raise ValueError(
                f"Unknown handle kind {kind!r}. Known kinds: {known}"
            )
        if not isinstance(summary, dict):
            raise TypeError("Handle summaries must be dict objects")

        handle_id = self._mint(kind)
        safe_summary = _json_safe(summary)
        self._objects[handle_id] = obj
        self._summaries[handle_id] = safe_summary
        self._kinds[handle_id] = kind
        with self.path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(
                    {
                        "handle_id": handle_id,
                        "kind": kind,
                        "summary": safe_summary,
                    },
                    sort_keys=True,
                )
                + "\n"
            )
        return handle_id

    def get(self, handle_id: str) -> Any:
        if handle_id not in self._objects:
            raise KeyError(handle_id)
        return self._objects[handle_id]

    def get_summary(self, handle_id: str) -> dict[str, Any]:
        if handle_id not in self._summaries:
            raise KeyError(handle_id)
        return deepcopy(self._summaries[handle_id])

    def _mint(self, kind: str) -> str:
        prefix = _HANDLE_PREFIXES[kind]
        while True:
            handle_id = f"{prefix}_{secrets.token_hex(2)}"
            if handle_id not in self._summaries:
                return handle_id

    def _load(self) -> None:
        if not self.path.exists():
            return
        with self.path.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                payload = json.loads(line)
                handle_id = payload["handle_id"]
                kind = payload["kind"]
                summary = payload["summary"]
                if not isinstance(summary, dict):
                    continue
                self._summaries[handle_id] = summary
                self._kinds[handle_id] = kind


def is_handle_id(value: str) -> bool:
    return HANDLE_ID_RE.match(value) is not None


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value)}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"type": value.__class__.__name__, "repr": repr(value)}
