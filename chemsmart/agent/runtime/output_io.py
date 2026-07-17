"""Bounded file readers for live and completed calculation outputs."""

from __future__ import annotations

import mmap
import os
from pathlib import Path


def read_output_window(path: Path) -> str:
    """Read bounded metadata from the start and results from the end."""

    try:
        size = path.stat().st_size
    except OSError:
        return ""
    tail = read_tail(path, max_chars=2_000_000)
    if size <= 8_000_000:
        return tail
    try:
        with path.open("rb") as handle:
            head = handle.read(800_000)
    except OSError:
        return tail
    sections = read_last_marker_windows(
        path,
        (b"VIBRATIONAL FREQUENCIES", b"Harmonic frequencies"),
    )
    return (
        head.decode("utf-8", errors="replace")
        + "\n"
        + "\n".join(sections)
        + "\n"
        + tail
    )


def read_tail(path: Path, *, max_chars: int) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - max_chars * 4))
            payload = handle.read()
    except OSError:
        return ""
    return payload.decode("utf-8", errors="replace")[-max_chars:]


def read_last_marker_windows(
    path: Path,
    markers: tuple[bytes, ...],
    *,
    after_bytes: int = 2_000_000,
) -> list[str]:
    windows: list[str] = []
    try:
        with (
            path.open("rb") as handle,
            mmap.mmap(
                handle.fileno(),
                0,
                access=mmap.ACCESS_READ,
            ) as mapped,
        ):
            for marker in markers:
                start = mapped.rfind(marker)
                if start < 0:
                    continue
                payload = mapped[start : min(len(mapped), start + after_bytes)]
                windows.append(payload.decode("utf-8", errors="replace"))
    except (OSError, ValueError):
        return []
    return windows


__all__ = ["read_last_marker_windows", "read_output_window", "read_tail"]
