"""Small deterministic measurements used by persisted agent receipts."""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any


def elapsed_ms(
    start_time: float | None,
    *,
    started_at: str | None = None,
    ended_at: str | None = None,
) -> int:
    if started_at is not None and ended_at is not None:
        started = datetime.fromisoformat(started_at)
        ended = datetime.fromisoformat(ended_at)
        return max(0, int(round((ended - started).total_seconds() * 1000)))
    if start_time is not None:
        return max(0, int(round((time.perf_counter() - start_time) * 1000)))
    return 0


def schema_hash(tool_defs: list[dict[str, Any]]) -> str:
    payload = json.dumps(tool_defs, sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def git_sha(module_path: str | Path) -> str | None:
    repo_root = Path(module_path).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


__all__ = ["elapsed_ms", "git_sha", "schema_hash"]
