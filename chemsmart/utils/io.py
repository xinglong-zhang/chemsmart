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
import re

import numpy as np

from chemsmart.io.molecules.structure import Molecule

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


def create_vibrationally_displaced_molecule(molecule, vib_mode, amp=0.5):
    """
    Create a displaced Molecule object from a given molecule and
    vibrational normal mode. Keep all other attributes from the molecule
    unchanged.
    Args:
        molecule: Molecule object
        vib_mode: vibrational mode
        amp (float): amplitude of vibrational displacement
    """
    symbols = list(molecule.symbols)
    R0 = np.asarray(molecule.positions, float)
    M = np.asarray(vib_mode)
    Rt = R0 + amp * M
    return Molecule(symbols=symbols, positions=Rt, **molecule)


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
