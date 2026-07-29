"""Create a small, redacted support bundle for one explicit agent session."""

from __future__ import annotations

import io
import json
import os
import platform
import re
import sys
import tarfile
from pathlib import Path
from typing import Any

from chemsmart import __version__
from chemsmart.agent.models import utc_now_iso
from chemsmart.agent.services.runtime_metrics import git_sha

_SESSION_FILES = (
    "decision_log.jsonl",
    "runtime_events.jsonl",
    "runtime_state.json",
    "session_metadata.json",
    "session.json",
    "state.json",
)
_SECRET_KEYS = re.compile(
    r"(api[-_]?key|token|password|secret|authorization|cookie)",
    re.IGNORECASE,
)
_INLINE_SECRETS = (
    re.compile(r"\bBearer\s+[A-Za-z0-9._~+/=-]+", re.IGNORECASE),
    re.compile(r"\bsk-[A-Za-z0-9_-]{6,}\b"),
    re.compile(
        r"(?i)\b(api[-_]?key|token|password|secret|authorization|cookie)"
        r"(\s*[:=]\s*)([^\s,;]+)"
    ),
)


class DebugBundleError(RuntimeError):
    """Raised when an explicit safe debug bundle cannot be created."""


def create_debug_bundle(
    session_root: str | os.PathLike[str],
    session_id: str,
    output_path: str | os.PathLike[str],
    *,
    workspace_path: str | os.PathLike[str] | None = None,
) -> dict[str, Any]:
    """Write redacted diagnostic evidence for exactly ``session_id``."""

    if session_id.strip().lower() == "latest":
        raise DebugBundleError(
            "'latest' is not allowed; pass an explicit session ID"
        )
    if not session_id.strip() or Path(session_id).name != session_id:
        raise DebugBundleError(
            "session ID must be one explicit directory name"
        )

    root = Path(session_root).expanduser().resolve()
    session_dir = (root / session_id).resolve()
    if session_dir.parent != root or not session_dir.is_dir():
        raise DebugBundleError(f"session not found: {session_id}")

    output = Path(output_path).expanduser().resolve()
    if output.exists():
        raise DebugBundleError(f"output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    home = Path.home().resolve()
    workspace = (
        Path(workspace_path).expanduser().resolve()
        if workspace_path is not None
        else None
    )

    members: list[tuple[str, bytes]] = []
    included: list[str] = []
    for filename in _SESSION_FILES:
        source = session_dir / filename
        if not source.is_file() or source.is_symlink():
            continue
        sanitized = _sanitize_file(source, home=home, workspace=workspace)
        members.append((f"session/{filename}", sanitized))
        included.append(filename)

    if not included:
        raise DebugBundleError(
            f"session has no supported diagnostic files: {session_id}"
        )

    manifest = {
        "schema_version": 1,
        "created_at": utc_now_iso(),
        "session_id": session_id,
        "included_files": included,
        "chemsmart_version": __version__,
        "chemsmart_package": str(Path(__file__).resolve().parents[2]),
        "git_sha": git_sha(Path(__file__)) or "unknown",
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "conda_environment": os.environ.get("CONDA_DEFAULT_ENV"),
    }
    manifest = _redact_value(manifest, home=home, workspace=workspace)
    members.append(
        (
            "environment_manifest.json",
            _json_bytes(manifest),
        )
    )

    try:
        with tarfile.open(output, mode="w:gz") as archive:
            for archive_name, payload in members:
                info = tarfile.TarInfo(archive_name)
                info.size = len(payload)
                info.mtime = 0
                info.mode = 0o600
                archive.addfile(info, io.BytesIO(payload))
    except OSError as exc:
        raise DebugBundleError(f"could not write debug bundle: {exc}") from exc

    return {
        "session_id": session_id,
        "output": str(output),
        "included_files": included,
    }


def _sanitize_file(
    path: Path,
    *,
    home: Path,
    workspace: Path | None,
) -> bytes:
    text = path.read_text(encoding="utf-8", errors="replace")
    if path.suffix == ".jsonl":
        lines = []
        for line in text.splitlines():
            if not line.strip():
                continue
            try:
                value = json.loads(line)
            except json.JSONDecodeError:
                lines.append(
                    _redact_text(line, home=home, workspace=workspace)
                )
            else:
                lines.append(
                    json.dumps(
                        _redact_value(value, home=home, workspace=workspace),
                        sort_keys=True,
                    )
                )
        return (("\n".join(lines) + "\n") if lines else "").encode("utf-8")

    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        return _redact_text(text, home=home, workspace=workspace).encode(
            "utf-8"
        )
    return _json_bytes(_redact_value(value, home=home, workspace=workspace))


def _redact_value(
    value: Any,
    *,
    home: Path,
    workspace: Path | None,
    key: str | None = None,
) -> Any:
    if key is not None and _SECRET_KEYS.search(key):
        return "<REDACTED>"
    if isinstance(value, dict):
        return {
            str(item_key): _redact_value(
                item_value,
                home=home,
                workspace=workspace,
                key=str(item_key),
            )
            for item_key, item_value in value.items()
        }
    if isinstance(value, list):
        return [
            _redact_value(item, home=home, workspace=workspace)
            for item in value
        ]
    if isinstance(value, str):
        return _redact_text(value, home=home, workspace=workspace)
    return value


def _redact_text(
    text: str,
    *,
    home: Path,
    workspace: Path | None,
) -> str:
    redacted = text
    replacements: list[tuple[str, str]] = []
    if workspace is not None:
        replacements.append((str(workspace), "<WORKSPACE>"))
    replacements.append((str(home), "<HOME>"))
    for original, replacement in sorted(
        replacements,
        key=lambda item: len(item[0]),
        reverse=True,
    ):
        redacted = redacted.replace(original, replacement)
    for pattern in _INLINE_SECRETS:
        if pattern.groups == 3:
            redacted = pattern.sub(r"\1\2<REDACTED>", redacted)
        else:
            redacted = pattern.sub("<REDACTED>", redacted)
    return redacted


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode("utf-8")


__all__ = ["DebugBundleError", "create_debug_bundle"]
