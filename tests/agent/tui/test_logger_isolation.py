from __future__ import annotations

import logging
import sys
from pathlib import Path

from chemsmart.agent.tui.app import _silence_console_logging, launch_tui


def _target_handlers(root: logging.Logger, log_path: Path):
    return [
        handler
        for handler in root.handlers
        if isinstance(handler, logging.FileHandler)
        and Path(handler.baseFilename) == log_path
    ]


def test_silence_console_logging_removes_root_stream_handlers(tmp_path: Path):
    root = logging.getLogger()
    original = list(root.handlers)
    root.handlers = [
        logging.StreamHandler(sys.stdout),
        logging.StreamHandler(sys.stderr),
    ]
    try:
        log_path = _silence_console_logging(tmp_path)
        assert all(
            not isinstance(handler, logging.StreamHandler)
            or isinstance(handler, logging.FileHandler)
            for handler in root.handlers
        )
        assert len(_target_handlers(root, log_path)) == 1
    finally:
        for handler in list(root.handlers):
            if handler not in original:
                handler.close()
        root.handlers = original


def test_silence_console_logging_adds_expected_file_handler(tmp_path: Path):
    root = logging.getLogger()
    original = list(root.handlers)
    root.handlers = [logging.StreamHandler(sys.stdout)]
    try:
        log_path = _silence_console_logging(tmp_path)
        assert log_path == tmp_path / "_tui.log"
        assert len(_target_handlers(root, log_path)) == 1
    finally:
        for handler in list(root.handlers):
            if handler not in original:
                handler.close()
        root.handlers = original


def test_silence_console_logging_is_idempotent(tmp_path: Path):
    root = logging.getLogger()
    original = list(root.handlers)
    root.handlers = [logging.StreamHandler(sys.stdout)]
    try:
        log_path = _silence_console_logging(tmp_path)
        _silence_console_logging(tmp_path)
        assert len(_target_handlers(root, log_path)) == 1
    finally:
        for handler in list(root.handlers):
            if handler not in original:
                handler.close()
        root.handlers = original


def test_launch_tui_smoke_with_patched_run(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        "chemsmart.agent.tui.app.ChemsmartTuiApp.run",
        lambda self, **kwargs: None,
    )

    launch_tui(session_root=tmp_path / "sessions")
