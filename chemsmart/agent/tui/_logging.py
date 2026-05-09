from __future__ import annotations

import logging
from pathlib import Path

_FORMATTER = logging.Formatter(
    "{asctime} - {levelname:6s} - [{name}] {message}",
    style="{",
)


def _silence_console_logging(
    session_root: Path,
    *,
    log_name: str = "_tui.log",
) -> Path:
    """Redirect console logging to a per-install agent log file."""
    session_root = Path(session_root)
    session_root.mkdir(parents=True, exist_ok=True)
    log_path = session_root / log_name
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    for handler in list(root.handlers):
        if isinstance(handler, logging.StreamHandler) and not isinstance(
            handler, logging.FileHandler
        ):
            root.removeHandler(handler)

    file_handlers = [
        handler
        for handler in root.handlers
        if isinstance(handler, logging.FileHandler)
        and Path(handler.baseFilename) == log_path
    ]
    for handler in file_handlers[1:]:
        root.removeHandler(handler)
        handler.close()

    existing_handler = next(iter(file_handlers), None)
    if existing_handler is not None:
        existing_handler.setLevel(logging.DEBUG)
        existing_handler.setFormatter(_FORMATTER)
        return log_path

    handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(_FORMATTER)
    root.addHandler(handler)
    return log_path


def _enable_console_logging(level: int = logging.DEBUG) -> None:
    """Ensure verbose agent runs keep console logging enabled."""
    root = logging.getLogger()
    root.setLevel(level)
    stream_handlers = [
        handler
        for handler in root.handlers
        if isinstance(handler, logging.StreamHandler)
        and not isinstance(handler, logging.FileHandler)
    ]
    if not stream_handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)
        handler.setFormatter(_FORMATTER)
        root.addHandler(handler)
        return

    for handler in stream_handlers:
        handler.setLevel(level)
        handler.setFormatter(_FORMATTER)


def _flush_logging_handlers() -> None:
    """Flush any active logging handlers before reporting log hints."""
    for handler in logging.getLogger().handlers:
        handler.flush()
