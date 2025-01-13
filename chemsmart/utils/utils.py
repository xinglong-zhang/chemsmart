import os
import re
import time
import hashlib
import copy
from functools import lru_cache, wraps
from itertools import groupby
from ase.io import string2index


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

            modified_times = tuple(
                os.path.getmtime(filename) for filename in filenames
            )

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

            result = func_with_modified_time_arg(
                tuple(modified_times), *args, **kwargs
            )
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
        list(group)
        for k, group in groupby(string_list, lambda x: x == "")
        if not k
    ]


def write_list_of_lists_as_a_string_with_empty_line_between_lists(
    list_of_lists,
):
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
    """Use 1-indexed.
    See chemsmart/tests/GaussianUtilsTest.py::GetListFromStringRangeTest::test_get_list_from_string_range
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
                indices.append(i)
        else:
            indices.append(int(s))
    return indices

def convert_list_to_gaussian_frozen_list(list_of_indices, molecule):
    """Convert a list of indices to a list of indices for Gaussian frozen list.
    Use 1-indexed list.
    """
    num_atoms_in_molecule = len(molecule.chemical_symbols)
    masks = [0] * num_atoms_in_molecule
    for i in list_of_indices:
        masks[i - 1] = -1
    return masks


def str_indices_to_list(str_indices):
    """Convert a supplied string of indices to a list of indices. All 0-indexed.

    Supported formats: e.g.,
        '1:9' -> gives [1,2,3,4,5,6,7,8]
        '1,2,4' -> gives [1,2,4]
        '1-9' -> gives [1,2,3,4,5,6,7,8]
        '[1-9]' -> gives [1,2,3,4,5,6,7,8].
    """
    try:
        indices_slice = string2index(str_indices)
        list_indices = list(range(*indices_slice.indices(100000)))
        # randomly large number (100000) to provide an upper bound onn the total number of components in the structure
    except (TypeError, ValueError, AttributeError):
        list_indices = []
        if "[" in str_indices:
            str_indices = str_indices.replace("[", "")
        if "]" in str_indices:
            str_indices = str_indices.replace("]", "")
        if "," in str_indices:
            str_indices_split = str_indices.split(",")
            for i in str_indices_split:
                list_indices.append(int(i))
        if ":" in str_indices:
            str_indices_split = str_indices.split(":")
            start_index = int(str_indices_split[0])
            end_index = int(str_indices_split[-1])
            for i in range(start_index, end_index):
                list_indices.append(i)
        if "-" in str_indices:
            str_indices_split = str_indices.split("-")
            start_index = int(str_indices_split[0])
            end_index = int(str_indices_split[-1])
            for i in range(start_index, end_index):
                list_indices.append(i)
    return list_indices

def get_value_by_number(num, data):
    # Iterate through all keys in the dictionary
    for key in data.keys():
        # Extract the numeric part of the key
        key_number = "".join(filter(str.isdigit, key))
        if key_number == str(num):
            return data[key]


def get_key_by_value_and_number(value, number, data):
    # Iterate through all items in the dictionary
    for key, val in data.items():
        # Extract the numerical part of the key using regex
        match = re.search(r"\d+$", key)
        if match:
            key_number = int(match.group())
            # Check if both the value and the numerical part of the key match
            if val == value and key_number == number:
                return key


def two_files_have_similar_contents(
    reference_file, output_file, ignored_string=None
):
    """Checks if the reference_file has same contents as output_file, except where ignored strings.

    reference_file: reference file
    output_file: output file
    ignored_string: If a string is ignored, it will not be compared.
    """
    with open(reference_file, "r") as f:
        reference_lines = f.readlines()
    with open(output_file, "r") as f:
        output_lines = f.readlines()
    for reference_line, output_line in zip(
        reference_lines, output_lines, strict=False
    ):
        if reference_line.strip() != output_line.strip():
            if ignored_string and ignored_string in reference_line:
                continue
            return False
    return True


def two_lists_have_similar_contents(list1, list2, ignore_string=None):
    """Checks if two lists have the same contents.
    If ignored string appears in the list member, then this member is ignored.
    """
    if len(list1) != len(list2):
        return False
    for item1, item2 in zip(list1, list2):
        if item1 != item2:
            if ignore_string and ignore_string in item1:
                continue
            return False
    return True


def update_dict_with_existing_keys(dict1, dict2):
    for k, v in dict2.items():
        if k in dict1:
            dict1[k] = v
        else:
            raise ValueError(
                f"Keyword `{k}` is not in list of keywords `{dict1.keys()}`\nPlease double check and rectify!"
            )
    return dict1


# Common utils for orca and gaussian
def get_prepend_string_for_modred(list_of_string):
    list_length = len(list_of_string)
    prepend_strings = {2: "B", 3: "A", 4: "D"}
    if list_length not in prepend_strings:
        raise ValueError(
            "Supplied list of coordinates should be 2 for Bond; 3 for Angle and 4 for Dihedral!"
        )
    return prepend_strings[list_length]


def convert_modred_list_to_string(modred_list):
    list_of_string = [str(a) for a in modred_list]
    return " ".join(list_of_string)


def get_prepend_string_list_from_modred_free_format(
    input_modred, program_type="gaussian"
):
    """Get prepend string list from either a list of lists or a list.
    e.g., [[1,2],[5,6]] or [1,2] -> ['B 1 2', 'B 5 6']
    program_type: gaussian or orca, variable to decide the index to start from
    Gaussian starts from 1, while orca starts from 0
    For chemsmart applications in homogeneous catalysis, we will use 1-indexing throughout,
    since this is what we usually see when we use GaussView or Avogadro to visualise the structures.
    """
    prepend_string_list = []
    if not isinstance(input_modred, list):
        raise ValueError(
            f"Required input for modredundant coordinates should be a list of lists "
            f"e.g., [[1,2],[5,6]] or a list e.g., [1,2].\n "
            f"But the given input is {input_modred} instead!"
        )
    if isinstance(input_modred[0], list):
        # for list of lists; e.g.: [[2,3],[4,5]]
        num_list = len(input_modred)
        for i in range(num_list):
            prepend_string = get_prepend_string_for_modred(input_modred[i])
            if program_type == "gaussian":
                modred_string = convert_modred_list_to_string(input_modred[i])
            elif program_type == "orca":
                modred_string = convert_modred_list_to_string(
                    [a - 1 for a in input_modred[i]]
                )
            else:
                raise ValueError(
                    f"Program type should be either gaussian or orca, but the given type is {program_type}!"
                )
            each_frozen_string = f"{prepend_string} {modred_string}"
            prepend_string_list.append(each_frozen_string)
    elif isinstance(input_modred[0], int):
        # for a single list; e.g.: [2,3]
        prepend_string = get_prepend_string_for_modred(input_modred)
        if program_type == "gaussian":
            modred_string = convert_modred_list_to_string(input_modred)
        elif program_type == "orca":
            modred_string = convert_modred_list_to_string(
                [a - 1 for a in input_modred]
            )
        each_frozen_string = f"{prepend_string} {modred_string}"
        prepend_string_list.append(each_frozen_string)
    return prepend_string_list


def prune_list_of_elements(list_of_elements, molecule):
    """Prune the list of elements so that only the specified elements that appear in a molecule is returned.
    E.g., one can specify heavy elements to be heavy_elements = ["Pd", "Ag", "Au"], but if only "Pd" appears
    in the molecule, then only "Pd" will be returned."""
    return list(
        set(molecule.chemical_symbols).intersection(set(list_of_elements))
    )
