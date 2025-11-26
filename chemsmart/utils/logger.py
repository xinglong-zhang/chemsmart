"""
Logging utilities for ChemSmart applications.

Provides customized logging functionality with deduplication capabilities
and configurable output streams. Supports both file and console logging
with different levels for information and error messages.
"""

import logging
import os
import re
import sys


class LogOnceFilter(logging.Filter):
    """
    Logging filter that prevents duplicate messages from being logged.

    Maintains a set of previously logged messages and filters out
    any messages that have already been logged (excluding timestamps).
    Useful for preventing repetitive log entries in long-running processes.
    """

    def __init__(self):
        """
        Initialize the LogOnceFilter with empty message tracking.

        Creates an empty set to store previously logged messages
        for duplicate detection.
        """
        super().__init__()
        self.logged_messages = set()

    def filter(self, record):
        """
        Filter duplicate log records.

        Checks if a log message (excluding timestamp) has been logged
        before and filters it out if it's a duplicate.

        Args:
            record (logging.LogRecord): The log record to evaluate.

        Returns:
            bool: True if message should be logged, False if duplicate.
        """
        formatted_message = self.format_record(record)
        stripped_message = self.remove_timestamp(formatted_message)

        if stripped_message in self.logged_messages:
            return False
        self.logged_messages.add(stripped_message)
        return True

    @staticmethod
    def format_record(record):
        """
        Format a log record into a string message.

        Applies string formatting to the log record message using
        any provided arguments.

        Args:
            record (logging.LogRecord): The log record to format.

        Returns:
            str: The formatted message string.
        """
        if record.args:
            return record.msg % record.args
        return record.msg

    @staticmethod
    def remove_timestamp(message):
        """
        Remove timestamp prefix from log message.

        Strips the standard timestamp format from the beginning
        of log messages for duplicate comparison.

        Args:
            message (str): The log message potentially with timestamp.

        Returns:
            str: The message with timestamp removed.
        """
        return re.sub(
            r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - ", "", message
        )


def create_logger(
    debug=True,
    folder=".",
    logfile=None,
    errfile=None,
    stream=True,
    disable=None,
):
    """
    Create and configure a logger with customizable output options.

    Sets up a root logger with file and stream handlers, duplicate message
    filtering, and configurable logging levels.

    Stream behavior:
    - Errors are always sent to stderr.
    - If `stream=True`, all messages (including errors) are also sent to stdout.

    Args:
        debug (bool, optional): Enable debug level logging. Defaults to True.
        folder (str, optional): Directory for log files. Defaults to ".".
        logfile (str, optional): Name of the info/debug log file. Defaults to None.
        errfile (str, optional): Name of the error log file. Defaults to None.
        stream (bool, optional): Enable console output to stdout (stderr is
            always used for errors). Defaults to True.
        disable (list[str], optional): Module names to disable logging for.
            Defaults to None.

    Returns:
        logging.Logger: Configured root logger instance.
    """
    if disable is None:
        disable = []

    for module in disable:
        logging.getLogger(module).disabled = True

    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logger = logging.getLogger()

    level = logging.DEBUG if debug else logging.INFO
    logger.setLevel(level)
    logger.handlers = []
    formatter = logging.Formatter(
        "{asctime} - {levelname:6s} - [{name}] {message}",
        style="{",
    )

    log_once_filter = LogOnceFilter()
    logger.addFilter(log_once_filter)

    # Stream errors always
    err_stream_handler = logging.StreamHandler(stream=sys.stderr)
    err_stream_handler.setLevel(logging.ERROR)
    err_stream_handler.setFormatter(formatter)
    logger.addHandler(err_stream_handler)

    # Stream other info only if required
    if stream:
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    # logfile
    if logfile:
        infofile_handler = logging.FileHandler(
            filename=os.path.join(folder, logfile)
        )
        infofile_handler.setLevel(level)
        infofile_handler.setFormatter(formatter)
        logger.addHandler(infofile_handler)

    # errfile
    if errfile:
        errfile_handler = logging.FileHandler(
            filename=os.path.join(folder, errfile)
        )
        errfile_handler.setLevel(logging.WARNING)
        errfile_handler.setFormatter(formatter)
        logger.addHandler(errfile_handler)

    return logger
