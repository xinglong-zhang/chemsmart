from __future__ import annotations

import logging
from pathlib import Path

_FORMATTER = logging.Formatter(
    "{asctime} - {levelname:6s} - [{name}] {message}",
    style="{",
)


class _AgentConsoleHandler(logging.StreamHandler):
    """Marker subclass for the agent's owned console handler.

    Used so that handler-level filters and replacements installed by the
    agent command pipeline (notably :func:`_apply_third_party_silence`)
    only touch handlers the agent itself attached, leaving pytest caplog
    and other host-attached ``StreamHandler`` instances untouched.
    """


_THIRD_PARTY_LOGGERS = (
    "numexpr.utils",
    "anthropic",
    "openai",
    "httpx",
    "httpcore",
    "urllib3",
    "chemsmart.settings",
    "chemsmart.agent.registry",
)


class _ThirdPartyConsoleFilter(logging.Filter):
    """Drop records from noisy 3rd-party + framework loggers.

    Attach to a console ``StreamHandler`` to suppress chatter without
    clamping the underlying logger's level — records still propagate to
    file handlers attached to root, so the agent log file stays complete.
    """

    def filter(self, record: logging.LogRecord) -> bool:
        for name in _THIRD_PARTY_LOGGERS:
            if record.name == name or record.name.startswith(name + "."):
                return False
        return True


def _apply_third_party_silence(level: int = logging.WARNING) -> None:
    """Suppress 3rd-party logger console output without losing records.

    Two parts:

    1. Strip any direct ``StreamHandler`` (non-file) handlers from the
       listed loggers — these bypass root and would print regardless of
       what console filtering root has installed (numexpr.utils historically
       attaches its own).
    2. Add :class:`_ThirdPartyConsoleFilter` to every root console handler.
       When ``level`` is ``DEBUG`` (i.e. ``--debug``), the filter is removed
       so records pass through to console too.

    Records continue propagating to root and reach the agent log file
    regardless of this filtering — debug forensics stay intact.
    """
    for name in _THIRD_PARTY_LOGGERS:
        log = logging.getLogger(name)
        # Reset level to NOTSET so records (DEBUG/INFO) are created and can
        # propagate to root's file handler for forensics. Modules that
        # imported chemsmart.cli.main may have set this to WARNING at import
        # time; undo that for the duration of the agent command.
        log.setLevel(logging.NOTSET)
        for handler in list(log.handlers):
            if isinstance(handler, logging.StreamHandler) and not isinstance(
                handler, logging.FileHandler
            ):
                log.removeHandler(handler)

    root = logging.getLogger()
    for handler in root.handlers:
        if not isinstance(handler, _AgentConsoleHandler):
            continue
        handler.filters = [
            f
            for f in handler.filters
            if not isinstance(f, _ThirdPartyConsoleFilter)
        ]
        if level > logging.DEBUG:
            handler.addFilter(_ThirdPartyConsoleFilter())


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
    agent_handlers = [
        handler
        for handler in root.handlers
        if isinstance(handler, _AgentConsoleHandler)
    ]
    if not agent_handlers:
        handler = _AgentConsoleHandler()
        handler.setLevel(level)
        handler.setFormatter(_FORMATTER)
        root.addHandler(handler)
        return

    for handler in agent_handlers:
        handler.setLevel(level)
        handler.setFormatter(_FORMATTER)


def _flush_logging_handlers() -> None:
    """Flush any active logging handlers before reporting log hints."""
    for handler in logging.getLogger().handlers:
        handler.flush()
