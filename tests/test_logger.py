"""Tests for chemsmart.utils.logger."""

import logging

import pytest

from chemsmart.utils.logger import LogOnceFilter, create_logger


@pytest.fixture(autouse=True)
def restore_root_logger():
    """create_logger mutates the root logger; restore state afterwards."""
    root = logging.getLogger()
    original_handlers = list(root.handlers)
    original_filters = list(root.filters)
    original_level = root.level
    yield
    root.handlers = original_handlers
    root.filters = original_filters
    root.setLevel(original_level)


class TestLogOnceFilterFormatRecord:
    def test_formats_message_with_args(self):
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="value is %s",
            args=("42",),
            exc_info=None,
        )
        assert LogOnceFilter.format_record(record) == "value is 42"

    def test_returns_msg_unchanged_without_args(self):
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="plain message",
            args=None,
            exc_info=None,
        )
        assert LogOnceFilter.format_record(record) == "plain message"


class TestLogOnceFilterRemoveTimestamp:
    def test_strips_leading_timestamp(self):
        message = "2024-01-15 10:30:00,123 - Some log message"
        assert LogOnceFilter.remove_timestamp(message) == "Some log message"

    def test_leaves_message_without_timestamp_unchanged(self):
        message = "Some log message"
        assert LogOnceFilter.remove_timestamp(message) == message


class TestLogOnceFilterFilter:
    def _make_record(self, msg):
        return logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg=msg,
            args=None,
            exc_info=None,
        )

    def test_allows_first_occurrence(self):
        log_filter = LogOnceFilter()
        record = self._make_record("hello world")
        assert log_filter.filter(record) is True

    def test_blocks_duplicate_message(self):
        log_filter = LogOnceFilter()
        record1 = self._make_record("hello world")
        record2 = self._make_record("hello world")
        assert log_filter.filter(record1) is True
        assert log_filter.filter(record2) is False

    def test_allows_different_messages(self):
        log_filter = LogOnceFilter()
        record1 = self._make_record("hello world")
        record2 = self._make_record("goodbye world")
        assert log_filter.filter(record1) is True
        assert log_filter.filter(record2) is True


class TestCreateLogger:
    def test_returns_root_logger_configured_debug(self):
        logger = create_logger(debug=True, stream=False)
        assert logger is logging.getLogger()
        assert logger.level == logging.DEBUG

    def test_info_level_when_debug_false(self):
        logger = create_logger(debug=False, stream=False)
        assert logger.level == logging.INFO

    def test_stream_true_adds_stdout_and_stderr_handlers(self):
        logger = create_logger(debug=True, stream=True)
        stream_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]
        assert len(stream_handlers) == 2

    def test_stream_false_adds_only_stderr_handler(self):
        logger = create_logger(debug=True, stream=False)
        stream_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]
        assert len(stream_handlers) == 1
        assert stream_handlers[0].level == logging.ERROR

    def test_logfile_creates_file_handler(self, tmp_path):
        logger = create_logger(
            stream=False, folder=str(tmp_path), logfile="info.log"
        )
        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert (tmp_path / "info.log").exists()

    def test_errfile_creates_file_handler(self, tmp_path):
        logger = create_logger(
            stream=False, folder=str(tmp_path), errfile="err.log"
        )
        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert file_handlers[0].level == logging.WARNING
        assert (tmp_path / "err.log").exists()

    def test_disable_sets_module_logger_disabled(self):
        create_logger(stream=False, disable=["some.module.to.disable"])
        assert logging.getLogger("some.module.to.disable").disabled is True

    def test_adds_log_once_filter(self):
        logger = create_logger(stream=False)
        assert any(isinstance(f, LogOnceFilter) for f in logger.filters)

    def test_resets_existing_handlers(self):
        logger = logging.getLogger()
        logger.addHandler(logging.NullHandler())
        create_logger(stream=False)
        assert not any(
            isinstance(h, logging.NullHandler) for h in logger.handlers
        )
