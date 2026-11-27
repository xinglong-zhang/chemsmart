"""Logging fixtures (loggers, log capture, root-logger isolation)."""

import logging
import os

import pytest


@pytest.fixture(scope="session")
def tests_logger():
    """Fixture to configure the root logger for tests."""
    logger = logging.getLogger()  # Root logger
    logger.setLevel(logging.INFO)
    logger.handlers = []  # Clear handlers to avoid conflicts
    logger.propagate = True
    # Set environment variable to signal test mode
    os.environ["TEST_MODE"] = "1"
    yield logger
    # Clean up
    logger.handlers = []
    os.environ.pop("TEST_MODE", None)


# Use built-in caplog fixture for capturing log messages
@pytest.fixture()
def capture_log(caplog):
    """
    Fixture to capture log messages.

    Captures messages from the root logger at DEBUG level by default.
    """
    caplog.set_level(logging.DEBUG, logger="")  # "" for root logger
    return caplog


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
