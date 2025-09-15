"""
General utility functions and classes for Chemsmart.

Provides a collection of utility functions for file handling, caching,
data structures, text processing, and system operations. Includes
specialized decorators for function caching and ordered set implementation.
"""

import copy
import hashlib
import logging
import os
import re
import shlex
import subprocess
import sys
import time
from functools import lru_cache, wraps
from itertools import groupby
from typing import Tuple, Union

import numpy as np

logger = logging.getLogger(__name__)


class OrderedSet:
    """
    Set-like container that maintains insertion order.

    Provides set operations while preserving the order in which items
    were added. Useful when uniqueness is required but order matters.
    """

    def __init__(self, iterable=None):
        """
        Initialize an OrderedSet with optional initial items.

        Args:
            iterable (optional): Initial items to add to the set.
        """
        self.items = []
        if iterable:
            for item in iterable:
                self.add(item)

    def add(self, item):
        """
        Add an item to the set if not already present.

        Args:
            item: Item to add to the set.
        """
        if item not in self.items:
            self.items.append(item)

    def remove(self, item):
        """
        Remove an item from the set if present.

        Args:
            item: Item to remove from the set.
        """
        if item in self.items:
            self.items.remove(item)

    def __contains__(self, item):
        """
        Check if an item is in the set.

        Args:
            item: Item to check for membership.

        Returns:
            bool: True if item is in the set, False otherwise.
        """
        return item in self.items

    def __iter__(self):
        """
        Return an iterator over the items in the set.

        Returns:
            iterator: Iterator over set items in insertion order.
        """
        return iter(self.items)

    def __len__(self):
        """
        Return the number of items in the set.

        Returns:
            int: Number of items in the set.
        """
        return len(self.items)


def file_cache(copy_result=True, maxsize=64):
    """
    Cache function results based on file modification times.

    Decorator that caches results of functions that take files as input.
    Uses the modified time of the files to see if they changed or not.
    If unchanged, will return the cached results of the function call.
    For recently modified files, uses file content hash for cache validation.

    Use before @classmethod or @staticmethod. For example:

    @classmethod
    @file_cache()
    def my_class_method(cls, ...):
        ...

    Args:
        copy_result (bool, optional): If true, will copy the result of the
                                    function call before caching. This
                                    prevents changes to the result from
                                    polluting the cache. Defaults to True.
        maxsize (int, optional): Maximum number of results to cache.
                               Defaults to 64.

    Returns:
        function: Decorated function with file-based caching.
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
    """
    Exception raised when file reading operations fail.

    Used to indicate errors during file parsing or reading operations
    in computational chemistry file processing.
    """

    pass


def is_float(string):
    """
    Test if a given string represents a float value.

    Checks if a string can be converted to a float, excluding
    strings that represent integers.

    Args:
        string (str): String to test for float representation.

    Returns:
        bool: True if string represents a float, False otherwise.
    """
    if string.replace("-", "").isdigit():
        return False  # if test string is a digit, then it is not a float
    try:
        float(string)
        return True
    except ValueError:
        return False


def strip_out_comments(string):
    """
    Remove comments from a string.

    Strips out comments by removing everything that follows after '#'
    for each line in the input string.

    Args:
        string (str): Input string potentially containing comments.

    Returns:
        str: String with comments removed.
    """
    return "\n".join(line.split("#")[0].strip() for line in string.split("\n"))


def content_blocks_by_paragraph(string_list):
    """
    Group lines into content blocks separated by empty lines.

    Takes a list of strings and groups them into blocks based on
    empty line separators.

    Args:
        string_list (list): List of strings to group into blocks.

    Returns:
        list: List of content blocks, each block is a list of strings.
    """
    return [
        list(group)
        for k, group in groupby(string_list, lambda x: x == "")
        if not k
    ]


def write_list_of_lists_as_a_string_with_empty_line_between_lists(
    list_of_lists,
):
    """
    Convert nested lists to a formatted string with separating empty lines.

    Converts a list of lists into a single string where each inner list
    becomes a block of lines, separated by empty lines between blocks.

    Args:
        list_of_lists (list): Nested list structure to convert.

    Returns:
        str: Formatted string with empty lines between list blocks.
    """
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
    """
    Use 1-indexed.
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
    """
    Convert atom indices to Gaussian frozen coordinate format.

    Creates a mask list for Gaussian input where specified atoms
    are marked for coordinate freezing. Uses 1-indexed numbering.

    Args:
        list_of_indices (list): List of atom indices to freeze (1-indexed).
        molecule: Molecule object containing chemical symbols.

    Returns:
        list: Binary mask list for Gaussian frozen coordinates.
    """
    num_atoms_in_molecule = len(molecule.chemical_symbols)
    masks = [0] * num_atoms_in_molecule
    for i in list_of_indices:
        masks[i - 1] = -1
    return masks


def convert_list_to_orca_frozen_list(list_of_indices, molecule):
    """
    Convert atom indices to ORCA frozen coordinate format.

    Creates a mask list for ORCA input where specified atoms
    are marked for coordinate freezing. Uses 0-indexed numbering.

    Args:
        list_of_indices (list): List of atom indices to freeze (0-indexed).
        molecule: Molecule object containing chemical symbols.

    Returns:
        list: Binary mask list for ORCA frozen coordinates.
    """
    num_atoms_in_molecule = len(molecule.chemical_symbols)
    masks = [0] * num_atoms_in_molecule
    for i in list_of_indices:
        masks[i] = -1
    return masks


def str_indices_range_to_list(str_indices):
    """
    Convert string index notation to a list of indices.

    Parses various string formats for specifying index ranges
    and converts them to a list of individual indices. All 1-indexed.

    Supported formats:
        '1:9' -> gives [1,2,3,4,5,6,7,8]
        '1,2,4' -> gives [1,2,4]
        '1-9' -> gives [1,2,3,4,5,6,7,8]
        '[1-9]' -> gives [1,2,3,4,5,6,7,8]

    Args:
        str_indices (str): String representation of indices in supported formats.

    Returns:
        list: List of individual integer indices.
    """
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


def string2index_1based(stridx: str) -> Union[int, slice, str]:
    """
    Convert string index to appropriate type using 1-based indexing.

    Wrapper function that converts string representations of indices
    to the appropriate Python types (int, slice) while handling
    1-based to 0-based index conversion.

    Args:
        stridx (str): String representation of index (e.g., '1', '1:5', '1:5:2').

    Returns:
        Union[int, slice, str]: Converted index in appropriate Python type.

    Raises:
        ValueError: If index is 0 or has invalid format.
    """

    def adjust_to_0based(index):
        """Adjust a 1-based index to 0-based. Handles None gracefully."""
        # raise error if index is 0, since requires 1-based indexing
        if index == 0:
            raise ValueError("Index cannot be 0 in 1-based indexing.")
        return index - 1 if index is not None else None

    # If it's not a slice, handle as a single integer or string
    if ":" not in stridx:
        try:
            if int(stridx) < 0:
                # return integer itself if negative
                return int(stridx)
            else:

                return int(stridx) - 1

        except ValueError:
            # If it's not an integer, check if it's alphanumeric
            if stridx.isnumeric():
                return stridx  # Return as-is if it is a number
            else:
                # Raise error for invalid input
                raise ValueError(f"Invalid input: {stridx}")
    else:
        # Handle slice input (e.g., "1:5" or "1:5:2")
        try:
            # Split slice into components and convert to integers or None
            i = [None if s == "" else int(s) for s in stridx.split(":")]
            # Adjust start and stop to 0-based only if start value is > 0
            if i[0] is not None and i[0] > 0:
                i[0] = adjust_to_0based(i[0])
            if i[1] is not None and i[1] > 0:
                i[1] = adjust_to_0based(i[1])
            # Return slice with adjusted start and stop values
            return slice(
                i[0],
                i[1],
                i[2] if len(i) > 2 else None,  # Step remains unchanged
            )
        except ValueError:
            # Raise error for invalid slice format
            raise ValueError(f"Invalid slice input: {stridx}")


def convert_string_index_from_1_based_to_0_based(
    index: Union[str, int, slice],
) -> Union[int, slice, list]:
    """
    Convert string index from 1-based to 0-based indexing.

    Handles single indices, slices, and negative indices. Supports
    various input formats including user-defined ranges.

    Args:
        index (Union[str, int, slice]): String representing an index or slice,
                                       e.g., "1", "2:5", "3:10:2", "-1",
                                       or user-defined ranges like "1-2,5".

    Returns:
        Union[int, slice, list]: Integer, slice, or list adjusted for
                                0-based indexing.

    Raises:
        ValueError: If index is 0 (invalid for 1-based indexing).
    """
    try:
        # Try numeric index
        index_int = int(index)
    except (TypeError, ValueError):
        try:
            index_list = string2index_1based(index)
        except ValueError:
            # Last resort: user-defined ranges
            index_list = get_list_from_string_range(index)
            # convert back to 0-based indexing
            index_list = [i - 1 for i in index_list]
    else:
        # Only runs if int() succeeded
        if index_int == 0:
            raise ValueError(
                f"Index {index_int} is out of range, as 1-indexing is used!\n "
                f"Please provide a positive integer.\n"
            )
        elif index_int < 0:
            # If negative index, return as is
            return index_int
        else:
            # Convert to 0-based indexing
            return index_int - 1  # Convert to 0-based indexing
    return index_list


def return_objects_from_string_index(list_of_objects, index):
    """
    Return objects from a list of objects based on the given index.
    The index can be a single integer, a slice, or a list of integers.
    If the index is negative, it will return the last n objects.
    If the index is a string, it will convert it to an integer or slice first.
    Index is 1-based, so inputting 0 will raise an error.
    """
    # convert index from 1-based (user input) to 0-based (python code-needed)
    index = convert_string_index_from_1_based_to_0_based(index)
    if isinstance(index, list):
        # if index is a list, use it to select objects
        objects = [list_of_objects[i] for i in index]
    elif isinstance(index, int):
        # if index is a single integer, use it to select a single object
        objects = list_of_objects[index]
    else:
        # index is a Slice
        objects = list_of_objects[index]

    # if not isinstance(objects, list):
    #     objects = [objects]

    return objects


def get_value_by_number(num, data):
    """
    Retrieve dictionary value by matching numeric part of keys.

    Searches through dictionary keys to find one containing the
    specified number and returns the corresponding value.

    Args:
        num (int): Number to search for in dictionary keys.
        data (dict): Dictionary to search through.

    Returns:
        Any: Value corresponding to the key containing the number.
    """
    # Iterate through all keys in the dictionary
    for key in data.keys():
        # Extract the numeric part of the key
        key_number = "".join(filter(str.isdigit, key))
        if key_number == str(num):
            return data[key]


def get_key_by_value_and_number(value, number, data):
    """
    Find dictionary key by matching both value and numeric suffix.

    Searches for a dictionary entry where both the value matches
    the specified value and the key ends with the specified number.

    Args:
        value (Any): Value to match in the dictionary.
        number (int): Number that should appear at the end of the key.
        data (dict): Dictionary to search through.

    Returns:
        str: Key that matches both criteria, or None if not found.
    """
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
    """
    Compare two files for similar contents with optional string exclusion.

    Checks if the reference file has the same contents as the output file,
    optionally ignoring lines containing specified strings.

    Args:
        reference_file (str): Path to the reference file.
        output_file (str): Path to the output file to compare.
        ignored_string (str, optional): If specified, lines containing this
                                       string will not be compared.

    Returns:
        bool: True if files have similar contents, False otherwise.
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
    """
    Compare two lists for identical contents with optional string exclusion.

    Checks if two lists have the same contents, optionally ignoring
    list members that contain a specified string.

    Args:
        list1 (list): First list to compare.
        list2 (list): Second list to compare.
        ignore_string (str, optional): If specified, list items containing
                                      this string will be ignored in comparison.

    Returns:
        bool: True if lists have similar contents, False otherwise.
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
    """
    Update dictionary with values from another dictionary.

    Updates dict1 with values from dict2, but only for keys that
    already exist in dict1. Raises an error if dict2 contains
    keys not present in dict1.

    Args:
        dict1 (dict): Dictionary to be updated (modified in place).
        dict2 (dict): Dictionary containing new values.

    Returns:
        dict: The updated dict1.

    Raises:
        ValueError: If dict2 contains keys not present in dict1.
    """
    for k, v in dict2.items():
        if k in dict1:
            dict1[k] = v
        else:
            raise ValueError(
                f"Keyword `{k}` is not in list of keywords `{dict1.keys()}`\n"
                f"Please double check and rectify!"
            )
    return dict1


# Common utils for orca and gaussian
def get_prepend_string_for_modred(list_of_string):
    """
    Get coordinate type prefix for molecular geometry modifications.

    Returns the appropriate prefix string based on the number of
    atoms involved in the coordinate specification.

    Args:
        list_of_string (list): List of atom indices for coordinate.

    Returns:
        str: Prefix string ('B' for bond, 'A' for angle, 'D' for dihedral).

    Raises:
        ValueError: If list length is not 2, 3, or 4.
    """
    list_length = len(list_of_string)
    prepend_strings = {2: "B", 3: "A", 4: "D"}
    if list_length not in prepend_strings:
        raise ValueError(
            "Supplied list of coordinates should be 2 for Bond; 3 for Angle and 4 for Dihedral!"
        )
    return prepend_strings[list_length]


def convert_modred_list_to_string(modred_list):
    """
    Convert a list of modification coordinates to a space-separated string.

    Args:
        modred_list (list): List of coordinate values to convert.

    Returns:
        str: Space-separated string of coordinate values.
    """
    list_of_string = [str(a) for a in modred_list]
    return " ".join(list_of_string)


def get_prepend_string_list_from_modred_free_format(
    input_modred, program="gaussian"
):
    """
    Get prepend string list from either a list of lists or a list.
    e.g., [[1,2],[5,6]] or [1,2] -> ['B 1 2', 'B 5 6']
    program: gaussian or orca, variable to decide the index to start from
    Gaussian starts from 1, while orca starts from 0
    For chemsmart applications in homogeneous catalysis, we will use 1-indexing throughout,
    since this is what we usually see when we use GaussView or Avogadro to visualise the structures.
    """
    prepend_string_list = []
    if not isinstance(input_modred, list):
        raise ValueError(
            f"Required input for modred coordinates should be a list of lists "
            f"e.g., [[1,2],[5,6]] or a list e.g., [1,2].\n "
            f"But the given input is {input_modred} instead!"
        )
    if isinstance(input_modred[0], list):
        # for list of lists; e.g.: [[2,3],[4,5]]
        num_list = len(input_modred)
        for i in range(num_list):
            prepend_string = get_prepend_string_for_modred(input_modred[i])
            if program == "gaussian" or program == "pymol":
                modred_string = convert_modred_list_to_string(input_modred[i])
            elif program == "orca":
                modred_string = convert_modred_list_to_string(
                    [a - 1 for a in input_modred[i]]
                )
            else:
                raise ValueError(
                    f"Program type should be either gaussian or orca, but the given type is {program}!"
                )
            each_frozen_string = f"{prepend_string} {modred_string}"
            prepend_string_list.append(each_frozen_string)
    elif isinstance(input_modred[0], int):
        # for a single list; e.g.: [2,3]
        prepend_string = get_prepend_string_for_modred(input_modred)
        if program == "gaussian" or program == "pymol":
            modred_string = convert_modred_list_to_string(input_modred)
        elif program == "orca":
            modred_string = convert_modred_list_to_string(
                [a - 1 for a in input_modred]
            )
        each_frozen_string = f"{prepend_string} {modred_string}"
        prepend_string_list.append(each_frozen_string)
    return prepend_string_list


def prune_list_of_elements(list_of_elements, molecule):
    """
    Filter element list to only include elements present in molecule.

    Prunes the list of elements so that only the specified elements
    that appear in a molecule are returned. Useful for filtering
    predefined element lists based on actual molecular composition.

    Example:
        Heavy elements specified as ["Pd", "Ag", "Au"], but if only "Pd"
        appears in the molecule, then only "Pd" will be returned.

    Args:
        list_of_elements (list): List of element symbols to filter.
        molecule: Molecule object with chemical_symbols attribute.

    Returns:
        list: Filtered list containing only elements present in molecule.
    """
    return list(
        set(molecule.chemical_symbols).intersection(set(list_of_elements))
    )


def sdf2molecule(sdf_lines):
    """
    Convert SDF format data to a Molecule object.

    Parses Structure Data Format (SDF) lines containing atomic coordinates
    and element information to create a Molecule object.

    Args:
        sdf_lines (Union[list, str]): SDF format data as list of lines
                                     or as a single string.

    Returns:
        Molecule: Molecule object created from SDF coordinate data.
    """
    from chemsmart.io.molecules.structure import Molecule

    list_of_symbols = []
    cart_coords = []

    # sdf line pattern containing coordinates and element type
    pattern = (
        r"\s*([\d\.-]+)\s+([\d\.-]+)\s+([\d\.-]+)\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)"
        r"\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*"
    )
    if isinstance(sdf_lines, list):
        line_elements = sdf_lines
    elif isinstance(sdf_lines, str):
        line_elements = sdf_lines.split("\n")

    for line in line_elements:
        match = re.match(pattern, line)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))
            atom_type = str(match.group(4))
            list_of_symbols.append(atom_type)
            cart_coords.append((x, y, z))

    cart_coords = np.array(cart_coords)
    return Molecule.from_symbols_and_positions_and_pbc_conditions(
        list_of_symbols=list_of_symbols, positions=cart_coords
    )


def check_charge_and_multiplicity(settings):
    """
    Validate that charge and multiplicity are set in settings.

    Checks if both charge and multiplicity are properly configured
    in the settings object, raising an error if either is missing.

    Args:
        settings: Settings object containing charge and multiplicity attributes.

    Raises:
        ValueError: If charge or multiplicity is None.
    """
    if settings.charge is None or settings.multiplicity is None:
        raise ValueError(
            "Charge and multiplicity must be set for Gaussian jobs."
        )


def cmp_with_ignore(f1, f2, ignore_string=None):
    """
    Compare two files with option to ignore lines containing specific strings.

    Compares the contents of two files while optionally ignoring lines
    that contain specified strings. Useful for comparing files that may
    have minor differences in timestamps, paths, or other variable content.

    Args:
        f1 (str): First file path to compare.
        f2 (str): Second file path to compare.
        ignore_string (Union[str, list], optional): String or list of strings
                                                   to ignore during comparison.
                                                   Lines containing any of these
                                                   strings will be skipped.

    Returns:
        bool: True if files are identical except for ignored lines,
              False otherwise.

    Raises:
        ValueError: If ignore_string is not a string or list of strings.
    """

    if ignore_string is not None:

        if isinstance(ignore_string, str):
            ignore_string = [ignore_string]
        elif isinstance(ignore_string, list):
            ignore_string = ignore_string
        else:
            raise ValueError(
                "ignore_string should be either a string or a list of strings."
            )
    with open(f1, "r") as fp1, open(f2, "r") as fp2:
        for line1, line2 in zip(fp1, fp2):
            # Skip lines containing the ignore_string
            if ignore_string and any(
                ignore in line1 or ignore in line2 for ignore in ignore_string
            ):
                continue
            # Compare the current lines
            if line1 != line2:
                return False
        # Ensure both files have reached EOF (no extra lines in one file)
        return fp1.read() == fp2.read()


def run_command(command):
    """
    Execute a shell command and capture its output.

    Runs a shell command using subprocess.Popen and captures its output.
    Handles both string commands and command lists with proper error handling.

    Args:
        command (Union[str, list]): Command to execute. Can be a string
                                  (e.g., 'ls -l') or a list
                                  (e.g., ['ls', '-l']).

    Returns:
        str: The command's stdout as a string, or None if an error occurs.
    """
    try:
        # If command is a string, split it into a list using shlex.split
        if isinstance(command, str):
            command = shlex.split(command)

        # Ensure command is a list at this point
        if not isinstance(command, list):
            logger.error(
                f"Invalid command type: {type(command)}. Expected str or list."
            )
            return None

        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.info(f"Error running {command}: {stderr.strip()}")
            return None
        return stdout.strip()
    except Exception as e:
        logger.error(f"Exception while running {command}: {e}")
        return None


def quote_path(path):
    """
    Quote file paths to handle spaces and special characters.

    Adds quotes around paths on Windows to handle spaces and backslashes
    properly. On other platforms, returns the path unchanged.

    Args:
        path (str): File path to quote.

    Returns:
        str: Quoted path on Windows, original path on other platforms.
    """
    if sys.platform == "win32":
        # Double-quote paths on Windows to preserve spaces
        return f'"{path}"'
    return path


def kabsch_align(
    P: np.ndarray, Q: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Align molecules using the Kabsch algorithm.

    Performs optimal alignment of two molecular structures using the
    Kabsch algorithm for finding the optimal rotation and translation
    that minimizes RMSD between corresponding atoms.

    Args:
        P (np.ndarray): First molecule coordinates (N x 3 array).
        Q (np.ndarray): Second molecule coordinates (N x 3 array).

    Returns:
        Tuple[np.ndarray, np.ndarray]: Aligned coordinates for both molecules.

    Raises:
        AssertionError: If matrix dimensions don't match.

    Reference:
        From https://hunterheidenreich.com/posts/kabsch_algorithm/
    """
    # Center molecules
    assert P.shape == Q.shape, "Matrix dimensions must match"

    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Optimal translation
    t = centroid_Q - centroid_P

    # Center the points
    p = P - centroid_P
    q = Q - centroid_Q

    # Compute the covariance matrix
    H = np.dot(p.T, q)

    # SVD
    U, S, Vt = np.linalg.svd(H)

    # Validate right-handed coordinate system
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0

    # Optimal rotation
    R = np.dot(Vt.T, U.T)

    # rotate p
    p = np.dot(p, R.T)

    # RMSD
    rmsd = np.sqrt(np.sum(np.square(p - q)) / P.shape[0])

    return p, q, R, t, rmsd


def kabsch_align2(
    P: np.ndarray, Q: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """
    Alternative Kabsch algorithm implementation for molecular alignment.

    Note: This implementation may not work correctly and is kept for
    reference. Use kabsch_align() for reliable molecular alignment.

    Args:
        P (np.ndarray): First molecule coordinates (N x 3 array).
        Q (np.ndarray): Second molecule coordinates (N x 3 array).

    Returns:
        Tuple containing:
        - aligned_P: P after rotation and translation to match Q's position
        - R: Optimal rotation matrix
        - t: Optimal translation vector
        - rmsd: Root Mean Square Deviation after alignment

    Raises:
        AssertionError: If input matrices have different dimensions.
    """
    assert P.shape == Q.shape, "Input matrices must have the same dimensions"

    # Compute centroids and center coordinates
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    centered_P = P - centroid_P
    centered_Q = Q - centroid_Q

    # Compute covariance matrix
    H = centered_P.T @ centered_Q

    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)

    # Ensure proper rotation (right-hand system)
    det_sign = np.linalg.det(Vt.T @ U.T)
    if det_sign < 0:
        Vt[-1, :] *= -1

    # Compute optimal rotation matrix
    R = Vt.T @ U.T

    # Apply rotation to centered_P (corrected from original)
    rotated_P = centered_P @ R  # Changed from R.T to R

    # Apply final translation to Q's coordinate system
    aligned_P = rotated_P + centroid_Q  # Translate to Q's position

    # Calculate RMSD using centered coordinates
    rmsd = np.sqrt(np.mean(np.sum((rotated_P - centered_Q) ** 2, axis=1)))

    return aligned_P, Q, R, centroid_Q - centroid_P, rmsd


def extract_number(filename):
    """
    Extract numeric part from filename for sorting purposes.

    Extracts the numeric component from filenames containing patterns
    like 'c123' or 'c456*'. Used for natural sorting of files.

    Args:
        filename (str): Filename to extract number from.

    Returns:
        int or float: Extracted number, or infinity if no number found
                     (to place files without numbers at the end during sorting).
    """
    match = re.search(r"c(\d+)\*?", filename)
    if match:
        return int(match.group(1))
    else:
        return float("inf")  # If no number is found, place it at the end


## file handling


def search_file(filename):
    """
    Search for a file in current directory and subdirectories.

    Searches for a file in the current directory and its subdirectories
    using the system's find command. Returns both absolute file path
    and directory path.

    Args:
        filename (str): Name of the file to search for.

    Returns:
        tuple: (absolute_file_path, absolute_file_dir) or (None, None)
               if file is not found or an error occurs.
    """
    try:
        # Search for the absolute file path
        result = subprocess.run(
            ["find", ".", "-name", filename],
            capture_output=True,
            text=True,
            check=True,
        )
        absolute_file_path = (
            result.stdout.strip().split("\n")[0]
            if result.stdout.strip()
            else None
        )

        # Search for the absolute directory path
        result = subprocess.run(
            ["find", ".", "-name", filename, "-exec", "dirname", "{}", ";"],
            capture_output=True,
            text=True,
            check=True,
        )
        absolute_file_dir = (
            result.stdout.strip().split("\n")[0]
            if result.stdout.strip()
            else None
        )

        if absolute_file_path and absolute_file_dir:
            return absolute_file_path, absolute_file_dir
        else:
            logger.error(f"{filename} not found! Check your Excel file.")
            return None, None
    except subprocess.CalledProcessError:
        logger.error(f"Error occurred while searching for {filename}.")
        return None, None


def iterative_compare(input_list):
    """
    Extract unique elements from a list while preserving order.

    Compares an input list and returns a list of unique elements,
    maintaining the order of first occurrence. Can handle lists
    of various data types including nested lists, strings, or dictionaries.

    Args:
        input_list (list): List to process for unique elements.

    Returns:
        list: List containing unique elements in order of first occurrence.
    """
    if not input_list:
        return []
    # Use OrderedSet to maintain order
    if isinstance(input_list, list) and len(input_list) > 0:
        return list(OrderedSet(input_list))


def naturally_sorted(lst):
    """
    Sort a list of strings in natural order with numeric awareness.

    Unlike standard alphabetical sorting, this function treats numerical
    parts as numbers rather than strings, ensuring proper ordering
    (e.g., "file10.txt" comes after "file2.txt").

    Args:
        lst (list): List of strings to sort naturally.

    Returns:
        list: Sorted list using natural ordering.
    """

    def key_func(key):
        return [
            int(c) if c.isdigit() else c.lower()
            for c in re.split("([0-9]+)", key)
        ]

    return sorted(lst, key=key_func)


def spline_data(x, y, new_length=1000, k=3):
    """
    Interpolate data points using univariate spline with even spacing.

    Interpolates data points using a univariate spline and returns
    evenly spaced points. Input points are automatically sorted by
    x-values to ensure proper interpolation.

    Args:
        x (list or array-like): X-coordinates of the input data points.
        y (list or array-like): Y-coordinates of the input data points.
        new_length (int, optional): Number of points in interpolated output.
                                   Defaults to 1000.
        k (int, optional): Degree of the spline (1 <= k <= 5).
                          Defaults to 3 (cubic spline).

    Returns:
        tuple: Two arrays (new_x, new_y) containing interpolated coordinates.

    Notes:
        - Uses scipy.interpolate.UnivariateSpline for interpolation.
        - Output new_x is evenly spaced between min and max input x-values.

    Raises:
        ValueError: If input lists x and y have different lengths or are empty.
        scipy.interpolate.InterpolationError: If spline interpolation fails
                                            (e.g., invalid k or insufficient points).
    """
    from scipy.interpolate import UnivariateSpline

    # Combine lists into list of tuples
    points = zip(x, y, strict=False)

    # Sort list of tuples by x-value
    points = sorted(points, key=lambda point: point[0])

    # Split list of tuples into two list of x values any y values
    x1, y1 = zip(*points, strict=False)
    new_x = np.linspace(min(x1), max(x1), new_length)
    new_y = UnivariateSpline(x1, y1, k=k)(new_x)
    return new_x, new_y


def safe_min_lengths(*lists):
    """Return the minimum length of any number of provided lists (ignoring None).

    Args:
        *lists: Variable number of list-like objects. None values are ignored.

    Returns:
        int: The minimum length among the provided lists, or 0 if no valid lists are given.
    """
    lengths = [len(lst) for lst in lists if lst is not None]
    return min(lengths) if lengths else 0
