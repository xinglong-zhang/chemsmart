from __future__ import annotations

import logging
from pathlib import Path

_FORMATTER = logging.Formatter(
    "{asctime} - {levelname:6s} - [{name}] {message}",
    style="{",
)


def _silence_console_logging(session_root: Path) -> Path:
    """Redirect root console logging to a per-install log file."""
    session_root = Path(session_root)
    session_root.mkdir(parents=True, exist_ok=True)
    log_path = session_root / "_tui.log"
    root = logging.getLogger()

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

    if file_handlers:
        file_handlers[0].setFormatter(_FORMATTER)
        return log_path

    handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    handler.setFormatter(_FORMATTER)
    root.addHandler(handler)
    return log_path
