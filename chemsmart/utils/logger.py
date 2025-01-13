import os
import sys
import logging


class LogOnceFilter(logging.Filter):
    def __init__(self):
        super().__init__()
        self.logged_messages = set()  # To track logged messages

    def filter(self, record):
        # Only log the message if it hasn't been logged before
        if record.msg in self.logged_messages:
            return False  # Skip the log message
        self.logged_messages.add(record.msg)
        return True  # Allow the log message


def create_logger(
    debug=True,
    folder=".",
    logfile=None,
    errfile=None,
    stream=True,
    disable=None,
):
    if disable is None:
        disable = []

    for module in disable:
        logging.getLogger(module).disabled = True

    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logger = logging.getLogger()

    # Stream
    level = logging.INFO
    if debug:
        level = logging.DEBUG

    logger.setLevel(level)
    logger.handlers = []
    formatter = logging.Formatter(
        '{asctime} - {levelname:6s} - [{name}] {message}',
        style="{",

    )

    # Add the filter to the logger
    log_once_filter = LogOnceFilter()
    logger.addFilter(log_once_filter)

    # Stream errors always
    err_stream_handler = logging.StreamHandler(stream=sys.stderr)
    err_stream_handler.setLevel(logging.ERROR)
    err_stream_handler.setFormatter(formatter)
    logger.addHandler(err_stream_handler)
    logger.addFilter(log_once_filter)

    # Stream other info only if required
    if stream:
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)
        logger.addFilter(log_once_filter)

    # logfile
    if logfile is not None:
        infofile_handler = logging.FileHandler(
            filename=os.path.join(folder, logfile)
        )
        infofile_handler.setLevel(level)
        infofile_handler.setFormatter(formatter)
        logger.addHandler(infofile_handler)
        logger.addFilter(log_once_filter)


    # errfile
    if errfile is not None:
        errfile_handler = logging.FileHandler(
            filename=os.path.join(folder, errfile)
        )
        errfile_handler.setLevel(logging.WARNING)
        errfile_handler.setFormatter(formatter)
        logger.addHandler(errfile_handler)
        logger.addFilter(log_once_filter)