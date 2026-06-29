from __future__ import annotations

from pathlib import Path
from typing import Any

_BINARY_SNIFF_BYTES = 8192


def read(
    path: str,
    start_line: int = 1,
    limit: int = 200,
) -> dict[str, Any]:
    """Read text from a local file using 1-based numbered lines."""

    cwd = Path.cwd().resolve(strict=False)
    resolved_path = Path(path).resolve(strict=False)

    if not resolved_path.is_relative_to(cwd):
        return {
            "error": "path_outside_cwd",
            "path": str(resolved_path),
            "cwd": str(cwd),
        }

    if not resolved_path.exists() or not resolved_path.is_file():
        return {
            "error": "file_not_found",
            "path": str(resolved_path),
        }

    with resolved_path.open("rb") as handle:
        if b"\x00" in handle.read(_BINARY_SNIFF_BYTES):
            return {
                "error": "binary_file_refused",
                "path": str(resolved_path),
            }

    start_index = max(start_line - 1, 0)
    line_limit = max(limit, 0)
    lines = resolved_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()
    total_lines = len(lines)
    selected_lines = lines[start_index : start_index + line_limit]
    if selected_lines:
        end_line = start_line + len(selected_lines) - 1
    else:
        end_line = start_line - 1
    truncated = start_index + len(selected_lines) < total_lines
    content = "\n".join(
        f"{line_number:>6}\t{line}"
        for line_number, line in enumerate(
            selected_lines,
            start=start_index + 1,
        )
    )

    return {
        "path": str(resolved_path),
        "start_line": start_line,
        "end_line": end_line,
        "total_lines": total_lines,
        "truncated": truncated,
        "content": content,
    }
