import os
import sys
import time
import hashlib
import copy
import logging
from functools import lru_cache, wraps
from itertools import groupby


def file_cache(copy_result=True, maxsize=64):
    """
    Caches results of functions that take files as input.

    Uses the modified time of the files to see if they changed or not.
    If unchanged, will return the cached results of the function call.

    Use before @classmethod or @staticmethod. For example:

    @classmethod
    @file_cache()
    def my_class_method(cls, ...):
        ...

    Args:
        copy_result (bool): If true, will copy the result of the function call before caching.
            This is to prevent changes to the result from polluting the cache.
        maxsize (int): Maximum number of results to cache. Defaults to 64.
    """

    def decorator(func):
        if isinstance(func, staticmethod):
            func = func.__func__
        elif isinstance(func, classmethod):
            raise ValueError(
                "Unable to use this with classmethod. Use a staticmethod instead."
            )

        @lru_cache(maxsize=maxsize)
        def func_with_modified_time_arg(modified_time, *args, **kwargs):
            return func(*args, **kwargs)

        @wraps(func)
        def wrapped(*args, **kwargs):
            filenames = []
            for arg in args + tuple(kwargs.values()):
                if isinstance(arg, (int, bool)):
                    continue
                if isinstance(arg, str) and os.path.isfile(arg):
                    filenames.append(arg)

            if not filenames:
                return func(*args, **kwargs)

            modified_times = tuple(os.path.getmtime(filename) for filename in filenames)

            current_time = time.time()
            time_diff = [current_time - t for t in modified_times]
            max_time_diff = 1800
            if any(t.is_integer() for t in modified_times) and any(
                t < max_time_diff for t in time_diff
            ):
                modified_times = []
                for filename in filenames:
                    with open(filename, "rb") as f:
                        data = f.read()
                    modified_times.append(
                        hashlib.blake2b(data, digest_size=10).hexdigest()
                    )

            result = func_with_modified_time_arg(tuple(modified_times), *args, **kwargs)
            if copy_result:
                return copy.copy(result)
            return result

        return wrapped

    return decorator


class FileReadError(Exception):
    pass


def is_float(string):
    """Function to test if a given string is float or not. This should exclude string that is a digit."""
    if string.replace("-", "").isdigit():
        return False  # if test string is a digit, then it is not a float
    try:
        float(string)
        return True
    except ValueError:
        return False


def content_blocks_by_paragraph(string_list):
    return [
        list(group) for k, group in groupby(string_list, lambda x: x == "") if not k
    ]


def write_list_of_lists_as_a_string_with_empty_line_between_lists(list_of_lists):
    string = ""
    num_lists = len(list_of_lists)
    for i in range(num_lists):
        for line in list_of_lists[i]:
            string += line + "\n"
        if i < num_lists - 1:
            # ensures that only empty line is written in between lists, but not after the last list
            # as the empty line after each block is written in gaussian job settings.write_gaussian_input()
            string += "\n"
    return string


def get_list_from_string_range(string_of_range):
    """Converted to 0-indexed from 1-indexed input.

    See pyatoms/tests/GaussianUtilsTest.py::GetListFromStringRangeTest::test_get_list_from_string_range
    :param string_of_range: accepts string of range. e.g., s='[1-3,28-31,34-41]' or s='1-3,28-31,34-41'
    :return: list of range.
    """
    if "[" in string_of_range:
        string_of_range = string_of_range.replace("[", "")
        string_of_range = string_of_range.replace("]", "")

    string_ranges = string_of_range.split(",")
    indices = []
    for s in string_ranges:
        if "-" in s:
            each_range = s.split("-")
            for i in range(int(each_range[0]), int(each_range[1]) + 1):
                indices.append(i - 1)  # noqa: PERF401
        else:
            indices.append(int(s))
    return indices


def create_logger(
    debug=True, folder=".", logfile=None, errfile=None, stream=True, disable=None
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

    # Stream errors always
    err_stream_handler = logging.StreamHandler(stream=sys.stderr)
    err_stream_handler.setLevel(logging.ERROR)
    logger.addHandler(err_stream_handler)

    # Stream other info only if required
    if stream:
        stream_handler = logging.StreamHandler(stream=sys.stdout)
        logger.addHandler(stream_handler)

    # logfile
    if logfile is not None:
        infofile_handler = logging.FileHandler(filename=os.path.join(folder, logfile))
        infofile_handler.setLevel(level)
        logger.addHandler(infofile_handler)

    # errfile
    if errfile is not None:
        errfile_handler = logging.FileHandler(filename=os.path.join(folder, errfile))
        errfile_handler.setLevel(logging.WARNING)
        logger.addHandler(errfile_handler)



def get_value_by_number(num, data):
    # Iterate through all keys in the dictionary
    for key in data.keys():
        # Extract the numeric part of the key
        key_number = ''.join(filter(str.isdigit, key))
        if key_number == str(num):
            return data[key]

def get_key_by_value(value, data):
    # Iterate through all items in the dictionary
    for key, val in data.items():
        if val == value:
            return key

# Sample dictionary for demonstration purposes
"""
data = {
    'O1': 'value1',
    'O2': 'value2',
    'O3': 'value3',
    'C4': 'value4',
    'C5': 'value5',
    'O6': 'value6',
    'O7': 'value7',
    'O8': 'value8',
    'C9': 'value9',
    'C10': 'value10',
    'C43': 'value43',
    'F68': 'value68',
    'H15': 'value15',
    # ... (add all other key-value pairs)
}
# Example usage
print(get_value_by_number(43, data))  # Should print 'value43'
print(get_value_by_number(1, data))   # Should print 'value1'
print(get_value_by_number(68, data))  # Should print 'value68'
print(get_value_by_number(15, data))  # Should print 'value15'
"""