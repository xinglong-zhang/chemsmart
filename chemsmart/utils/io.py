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


def match_outfile_pattern(line) -> str | None:
    """
    Match a line of text to known quantum chemistry program signatures.

    Args:
        line (str): Line from an output file.

    Returns:
        str | None: Program name ("gaussian", "orca", "xtb", "crest") if matched, else None.
    """
    patterns = {
        "crest": [
            "C R E S T",
            "Conformer-Rotamer Ensemble Sampling Tool",
            "https://crest-lab.github.io/crest-docs/",
            "$ crest",
        ],
        "gaussian": [
            "Entering Gaussian System",
            "Gaussian, Inc.",
            "Gaussian(R)",
        ],
        "orca": [
            "* O   R   C   A *",
            "Your ORCA version",
            "ORCA versions",
        ],
        "xtb": ["x T B", "xtb version", "xtb is free software:"],
    }
    for program, keywords in patterns.items():
        if any(keyword in line for keyword in keywords):
            return program
    return None


def get_outfile_format(filepath) -> str:
    """
    Detect the type of quantum chemistry output file.

    Reads only the first 200 lines and scans for characteristic keywords
    of major QC packages to improve efficiency on large output files.

    Args:
        filepath (str): Path to the quantum chemistry output file.

    Returns:
        str: Program name, one of: "gaussian", "orca", "xtb", "crest", or "unknown" if the format cannot be detected.
    """
    max_lines = 200
    try:
        with open(filepath, "r") as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                stripped = line.strip()
                if not stripped:
                    continue
                if program := match_outfile_pattern(stripped):
                    logger.debug(
                        f"Detected output format for '{os.path.basename(filepath)}': {program}."
                    )
                    return program
    except Exception as e:
        logger.error(f"Error reading file '{filepath}': {e}.")
        return "unknown"

    logger.debug(
        f"Could not detect output format for '{os.path.basename(filepath)}'."
    )
    return "unknown"
