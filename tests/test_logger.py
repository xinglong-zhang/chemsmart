"""
Direct unit tests for chemsmart.utils.logger.

create_logger() is exercised indirectly by nearly every CLI test (it is
called on every ``chemsmart`` CLI invocation), but only with its default
arguments, and LogOnceFilter is never exercised directly at all. This
file targets the branches that integration-style CLI tests don't reach:
debug/stream toggles, logfile/errfile handler creation, the ``disable``
list, and LogOnceFilter's deduplication/formatting helpers.
"""

import logging
import sys

import pytest

from chemsmart.utils.logger import LogOnceFilter, create_logger


@pytest.fixture()
def clean_root_logger():
    """Snapshot and restore the root logger around each test.

    create_logger() mutates the (global, shared) root logger's level,
    handlers, and filters, so tests must not leak that state into the
    rest of the suite.
    """
    root = logging.getLogger()
    original_handlers = root.handlers[:]
    original_filters = root.filters[:]
    original_level = root.level
    yield root
    for handler in root.handlers:
        if handler not in original_handlers:
            handler.close()
    root.handlers = original_handlers
    root.filters = original_filters
    root.level = original_level


def _make_record(msg, args=None):
    return logging.LogRecord(
        name="test",
        level=logging.INFO,
        pathname="",
        lineno=0,
        msg=msg,
        args=args,
        exc_info=None,
    )


class TestLogOnceFilterFormatRecord:
    def test_formats_message_with_args(self):
        record = _make_record("value is %s", args=("foo",))
        assert LogOnceFilter.format_record(record) == "value is foo"

    def test_formats_message_without_args(self):
        record = _make_record("plain message")
        assert LogOnceFilter.format_record(record) == "plain message"


class TestLogOnceFilterRemoveTimestamp:
    def test_strips_leading_timestamp(self):
        message = "2026-01-01 12:00:00,123 - hello"
        assert LogOnceFilter.remove_timestamp(message) == "hello"

    def test_leaves_message_without_timestamp_unchanged(self):
        message = "hello world"
        assert LogOnceFilter.remove_timestamp(message) == "hello world"


class TestLogOnceFilterFilter:
    def test_first_occurrence_passes(self):
        log_filter = LogOnceFilter()
        assert log_filter.filter(_make_record("hello")) is True

    def test_duplicate_message_is_filtered_out(self):
        log_filter = LogOnceFilter()
        log_filter.filter(_make_record("hello"))
        assert log_filter.filter(_make_record("hello")) is False

    def test_distinct_messages_both_pass(self):
        log_filter = LogOnceFilter()
        assert log_filter.filter(_make_record("hello")) is True
        assert log_filter.filter(_make_record("world")) is True

    def test_formatted_args_are_deduplicated(self):
        # Two records with the same msg template and args collapse to
        # the same formatted string and are deduplicated accordingly.
        log_filter = LogOnceFilter()
        log_filter.filter(_make_record("value is %s", args=("foo",)))
        assert (
            log_filter.filter(_make_record("value is %s", args=("foo",)))
            is False
        )


class TestCreateLogger:
    def test_defaults_set_debug_level_and_both_streams(
        self, clean_root_logger
    ):
        logger = create_logger()
        assert logger.level == logging.DEBUG
        assert any(
            isinstance(h, logging.StreamHandler) and h.stream is sys.stderr
            for h in logger.handlers
        )
        assert any(
            isinstance(h, logging.StreamHandler) and h.stream is sys.stdout
            for h in logger.handlers
        )
        assert any(isinstance(f, LogOnceFilter) for f in logger.filters)

    def test_debug_false_sets_info_level(self, clean_root_logger):
        logger = create_logger(debug=False)
        assert logger.level == logging.INFO

    def test_stream_false_omits_stdout_handler(self, clean_root_logger):
        logger = create_logger(stream=False)
        stream_handlers = [
            h for h in logger.handlers if isinstance(h, logging.StreamHandler)
        ]
        # Only the always-on stderr error handler should be present.
        assert len(stream_handlers) == 1
        assert stream_handlers[0].stream is sys.stderr

    def test_logfile_creates_file_handler_at_configured_level(
        self, clean_root_logger, tmp_path
    ):
        logfile = "test.log"
        logger = create_logger(
            debug=False, folder=str(tmp_path), logfile=logfile
        )

        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert file_handlers[0].baseFilename == str(tmp_path / logfile)
        assert file_handlers[0].level == logging.INFO

    def test_errfile_creates_file_handler_at_warning_level(
        self, clean_root_logger, tmp_path
    ):
        errfile = "test.err"
        logger = create_logger(folder=str(tmp_path), errfile=errfile)

        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert file_handlers[0].baseFilename == str(tmp_path / errfile)
        assert file_handlers[0].level == logging.WARNING

    def test_logfile_and_errfile_both_create_separate_handlers(
        self, clean_root_logger, tmp_path
    ):
        logger = create_logger(
            folder=str(tmp_path), logfile="out.log", errfile="out.err"
        )
        file_handlers = [
            h for h in logger.handlers if isinstance(h, logging.FileHandler)
        ]
        assert len(file_handlers) == 2

    def test_disable_marks_listed_loggers_disabled(self, clean_root_logger):
        module_name = "chemsmart.tests.fake_module_for_logger_test"
        try:
            create_logger(disable=[module_name])
            assert logging.getLogger(module_name).disabled is True
        finally:
            logging.getLogger(module_name).disabled = False

    def test_disable_none_defaults_to_no_modules_disabled(
        self, clean_root_logger
    ):
        # Should not raise when disable is omitted (defaults to None).
        logger = create_logger()
        assert logger is logging.getLogger()

    def test_repeated_calls_reset_handlers_and_filters(
        self, clean_root_logger
    ):
        create_logger()
        logger = logging.getLogger()
        first_handler_count = len(logger.handlers)
        first_filter_count = len(logger.filters)

        create_logger()

        assert len(logger.handlers) == first_handler_count
        assert len(logger.filters) == first_filter_count

    def test_returns_root_logger(self, clean_root_logger):
        assert create_logger() is logging.getLogger()
