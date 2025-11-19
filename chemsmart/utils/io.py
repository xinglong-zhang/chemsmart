"""
Input/output utility functions for molecular structure processing.

This module provides helper functions for creating molecule objects,
cleaning duplicate structures, and text processing operations commonly
used in computational chemistry file I/O operations.

Key functionality includes:
- Molecule object creation from coordinate data
- Duplicate structure detection and removal
- Text processing for chemical file formats
"""

import logging
import os
import re
from pathlib import Path

import numpy as np

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.repattern import float_pattern_with_exponential

logger = logging.getLogger(__name__)


def create_molecule_list(
    orientations,
    orientations_pbc,
    energies,
    forces,
    symbols,
    charge,
    multiplicity,
    frozen_atoms,
    pbc_conditions,
    num_structures=None,
):
    """
    Helper to build a list of Molecule objects from arrays.

    Creates Molecule instances from per-structure coordinates and optional
    per-structure properties.

    Args:
        orientations (list[np.ndarray]): Cartesian coordinates per structure,
            each with shape (N, 3).
        orientations_pbc (list[np.ndarray]): Per-structure translation vectors
            (periodic cell info) aligned with `orientations`.
        energies (list[float] | None): Per-structure energies (optional).
        forces (list[np.ndarray] | None): Per-structure forces (optional).
        symbols (list[str]): Chemical symbols for the atoms.
        charge (int): Molecular charge.
        multiplicity (int): Spin multiplicity.
        frozen_atoms (list[int] | None): Indices of frozen atoms.
        pbc_conditions (list | None): Periodic boundary conditions.
        num_structures (int, optional): Number of structures to create; if None
            uses `len(orientations)`.

    Returns:
        list[Molecule]: Molecule objects with specified properties.

    Notes:
        - `orientations_pbc`, `energies`, and `forces` are indexed by structure;
          when provided, they should be at least `num_structures` long.
        - This function does not validate shapes beyond basic indexing.
    """
    num_structures = num_structures or len(orientations)
    logger.debug(f"Number of structures to create: {num_structures}")
    logger.debug(f"Orientations shape: {np.array(orientations).shape}")
    logger.debug(
        f"Number of PBC orientations: {len(orientations_pbc) if orientations_pbc else 'N/A'}"
    )
    logger.debug(f"Number of energies: {len(energies) if energies else 'N/A'}")
    logger.debug(
        f"Energies shape: {np.array(energies).shape if energies else 'N/A'}"
    )
    logger.debug(f"Number of forces: {len(forces) if forces else 'N/A'}")
    logger.debug(
        f"Forces shape: {np.array(forces).shape if forces else 'N/A'}"
    )

    return [
        Molecule(
            symbols=symbols,
            positions=orientations[i],
            translation_vectors=orientations_pbc[i],
            charge=charge,
            multiplicity=multiplicity,
            frozen_atoms=frozen_atoms,
            pbc_conditions=pbc_conditions,
            energy=energies[i] if energies else None,
            forces=forces[i] if forces else None,
        )
        for i in range(num_structures)
    ]


def clean_duplicate_structure(orientations):
    """
    Remove trailing duplicate structure in place.

    If the last structure is nearly identical to the second-to-last (within
    `np.allclose` tolerance), remove the last entry. Useful for pruning
    redundant frames in optimization trajectories.

    Args:
        orientations (list[np.ndarray]): Coordinate arrays to check.

    Notes:
        - Uses `np.allclose` with `rtol=1e-5` (default `atol`).
        - Mutates `orientations`; returns None.
    """
    if orientations and len(orientations) > 1:
        if np.allclose(orientations[-1], orientations[-2], rtol=1e-5):
            orientations.pop(-1)


def increment_numbers(s, increment=1):
    """
    Increment all integers in a string by a fixed amount.

    Finds decimal integer substrings (no sign or decimal point) and adds
    `increment` to each. Useful for adjusting atom indices in text files.

    Args:
        s (str): Input string containing integers to increment.
        increment (int): Amount to add to each integer. Defaults to 1.

    Returns:
        str: String with all matched integers incremented.
    """
    # Use re.sub with a lambda to increment each number
    return re.sub(r"\d+", lambda m: str(int(m.group()) + increment), s)


def remove_keyword(text, keyword):
    """
    Remove a keyword from text (case-insensitive, word boundaries).

    Removes all occurrences of `keyword` using `\b` word-boundary matching
    to avoid partial matches.

    Args:
        text (str): Input text to process.
        keyword (str): Keyword to remove.

    Returns:
        str: Text with the keyword removed.

    Note:
        This does not trim extra whitespace that may remain after removal.
    """
    return re.sub(
        r"\b" + re.escape(keyword) + r"\b", "", text, flags=re.IGNORECASE
    )


def line_of_all_integers(line: str, allow_sign: bool = True) -> bool:
    """
    Return True iff the line has 1+ whitespace-separated tokens
    and every token is an integer.
    Examples return True: "0 1 23", "+3 -5 0" (when allow_sign=True)
    """
    tokens = line.strip().split()
    if not tokens:
        return False
    try:
        for t in tokens:
            if not allow_sign and (t.startswith(("+", "-"))):
                return False
            int(t)  # raises ValueError if not an integer literal
        return True
    except ValueError:
        return False


def line_of_integer_followed_by_floats(line) -> bool:
    """
    Return True iff the line has tokens and:
      - first token is an integer (± allowed),
      - remaining tokens are floats.
    Options:
      strict_float=True  -> require decimal point or exponent in floats
      min_floats=1       -> require at least this many float tokens after the integer
    """
    float_pattern = re.compile(float_pattern_with_exponential)
    tokens = line.split()

    if len(tokens) < 2:
        return False

    # First token: integer (allows + / -)
    try:
        int(tokens[0])
    except ValueError:
        return False

    # Remaining tokens: floats
    return all(float_pattern.fullmatch(t) for t in tokens[1:])


def load_molecules_from_paths(
    file_paths,
    index,
    add_index_suffix_for_single=False,
    check_exists=False,
):
    """
    Load molecules from a list of file paths, assigning unique names to each molecule.

    For each file in `file_paths`, this function loads one or more molecular structures
    using `Molecule.from_filepath`, assigns a unique name to each molecule based on the
    file name and structure index, and returns a list of all loaded molecules.

    Args:
        file_paths (list of str or Path): List of file paths to load molecules from.
        index (int or str): Index or slice to select specific structures from each file.
        add_index_suffix_for_single (bool, optional): If True, appends an index suffix to
            the molecule name even if only a single structure is loaded from a file.
        check_exists (bool, optional): If True, checks that each file exists before loading.

    Returns:
        list of Molecule: List of loaded Molecule objects, each with a unique name.

    Raises:
        FileNotFoundError: If `check_exists` is True and a file does not exist.
        Exception: If an error occurs during molecule loading from a file.
    """
    loaded = []

    for i, file_path in enumerate(file_paths):
        logger.debug(f"Processing file {i+1}/{len(file_paths)}: {file_path}")

        if not file_path:
            logger.warning(f"Skipping invalid file path: {file_path}")
            continue

        p = Path(file_path)

        if check_exists and not p.is_file():
            logger.error(f"File not found or not a regular file: {file_path}")
            raise FileNotFoundError(f"File not found: {file_path}")

        file_path = str(p)

        try:
            mols = Molecule.from_filepath(
                filepath=file_path,
                index=index,
                return_list=True,
            )

            # assign unique names per-structure when file contains multiple structures
            base = os.path.splitext(os.path.basename(file_path))[0]
            if isinstance(mols, list) and len(mols) > 1:
                for j, mol in enumerate(mols, start=1):
                    mol.name = f"{base}_{j}"
            else:
                # Optional suffix for single-structure files (filenames branch).
                if add_index_suffix_for_single and index not in (":", "-1"):
                    for mol in mols:
                        mol.name = f"{base}_idx{index}"
                else:
                    for mol in mols:
                        mol.name = base

            loaded += mols
            logger.debug(
                f"Successfully loaded {len(mols)} molecules from {file_path}"
            )
        except Exception as e:
            logger.error(f"Error loading molecules from {file_path}: {e}")
            raise

    return loaded
