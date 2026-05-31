"""Tests for agent third-party console-silencing helpers."""

from __future__ import annotations

import logging

import pytest

from chemsmart.agent.tui._logging import (
    _THIRD_PARTY_LOGGERS,
    _AgentConsoleHandler,
    _ThirdPartyConsoleFilter,
    _apply_third_party_silence,
)


@pytest.fixture
def isolated_root_handlers():
    """Snapshot and restore the root logger's handler list around tests."""
    root = logging.getLogger()
    original_handlers = list(root.handlers)
    original_level = root.level
    yield
    root.handlers[:] = original_handlers
    root.setLevel(original_level)


def _make_agent_handler() -> _AgentConsoleHandler:
    handler = _AgentConsoleHandler()
    handler.setLevel(logging.DEBUG)
    return handler


def test_third_party_filter_drops_listed_logger_records() -> None:
    f = _ThirdPartyConsoleFilter()
    record = logging.LogRecord(
        name="numexpr.utils",
        level=logging.INFO,
        pathname=__file__,
        lineno=0,
        msg="hi",
        args=(),
        exc_info=None,
    )
    assert f.filter(record) is False


def test_third_party_filter_drops_sublogger_records() -> None:
    f = _ThirdPartyConsoleFilter()
    record = logging.LogRecord(
        name="openai.http",
        level=logging.DEBUG,
        pathname=__file__,
        lineno=0,
        msg="hi",
        args=(),
        exc_info=None,
    )
    assert f.filter(record) is False


def test_third_party_filter_lets_other_loggers_through() -> None:
    f = _ThirdPartyConsoleFilter()
    record = logging.LogRecord(
        name="chemsmart.agent.synthesis",
        level=logging.INFO,
        pathname=__file__,
        lineno=0,
        msg="hi",
        args=(),
        exc_info=None,
    )
    assert f.filter(record) is True


def test_apply_third_party_silence_attaches_filter_to_agent_handler(
    isolated_root_handlers,
) -> None:
    root = logging.getLogger()
    root.handlers[:] = []
    agent_handler = _make_agent_handler()
    root.addHandler(agent_handler)

    _apply_third_party_silence(logging.WARNING)

    assert any(
        isinstance(f, _ThirdPartyConsoleFilter)
        for f in agent_handler.filters
    )


def test_apply_third_party_silence_with_debug_removes_filter(
    isolated_root_handlers,
) -> None:
    root = logging.getLogger()
    root.handlers[:] = []
    agent_handler = _make_agent_handler()
    agent_handler.addFilter(_ThirdPartyConsoleFilter())
    root.addHandler(agent_handler)

    _apply_third_party_silence(logging.DEBUG)

    assert all(
        not isinstance(f, _ThirdPartyConsoleFilter)
        for f in agent_handler.filters
    )


def test_apply_third_party_silence_leaves_foreign_stream_handlers_alone(
    isolated_root_handlers,
) -> None:
    """Non-agent StreamHandlers (e.g. pytest caplog) must not be touched."""
    root = logging.getLogger()
    root.handlers[:] = []
    foreign = logging.StreamHandler()
    foreign.setLevel(logging.DEBUG)
    root.addHandler(foreign)

    _apply_third_party_silence(logging.WARNING)

    assert all(
        not isinstance(f, _ThirdPartyConsoleFilter) for f in foreign.filters
    )


def test_apply_third_party_silence_strips_direct_handlers_from_listed_loggers(
    isolated_root_handlers,
) -> None:
    targets = [logging.getLogger(name) for name in _THIRD_PARTY_LOGGERS]
    stamped = []
    for log in targets:
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        log.addHandler(handler)
        stamped.append((log, handler))

    try:
        _apply_third_party_silence(logging.WARNING)
        for log, handler in stamped:
            assert handler not in log.handlers
    finally:
        for log, handler in stamped:
            if handler in log.handlers:
                log.removeHandler(handler)
