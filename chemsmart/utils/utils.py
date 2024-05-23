import os
import time
import hashlib
import copy
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
