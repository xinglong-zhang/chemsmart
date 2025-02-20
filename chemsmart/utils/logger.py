import logging
import os
import re
import sys


class LogOnceFilter(logging.Filter):
    def __init__(self):
        super().__init__()
        self.logged_messages = set()

    def filter(self, record):
        formatted_message = self.format_record(record)
        stripped_message = self.remove_timestamp(formatted_message)

        if stripped_message in self.logged_messages:
            return False
        self.logged_messages.add(stripped_message)
        return True

    @staticmethod
    def format_record(record):
        if record.args:
            return record.msg % record.args
        return record.msg

    @staticmethod
    def remove_timestamp(message):
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
