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

from __future__ import annotations

import logging
import os
import re
import shutil
import string
import subprocess
from io import BytesIO
from pathlib import Path
from typing import List

import numpy as np
from rdkit import Chem

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.repattern import float_pattern_with_exponential

logger = logging.getLogger(__name__)

SAFE_CHARS = set(string.ascii_letters + string.digits + "_-")

PROGRAM_INFO = {
    "crest": {
        "keywords": [
            "C R E S T",
            "Conformer-Rotamer Ensemble Sampling Tool",
            "https://crest-lab.github.io/crest-docs/",
            "$ crest",
        ],
        "suffixes": [".out"],
    },
    "gaussian": {
        "keywords": [
            "Entering Gaussian System",
            "Gaussian, Inc.",
            "Gaussian(R)",
        ],
        "suffixes": [".log", ".out"],
    },
    "orca": {
        "keywords": [
            "* O   R   C   A *",
            "Your ORCA version",
            "ORCA versions",
        ],
        "suffixes": [".out"],
    },
    "xtb": {
        "keywords": ["x T B", "xtb version", "xtb is free software:"],
        "suffixes": [".out"],
    },
}
SUPPORTED_PROGRAMS = set(PROGRAM_INFO.keys())
ALL_SUFFIXES = tuple(
    suffix for info in PROGRAM_INFO.values() for suffix in info["suffixes"]
)
# Folder-level detection is currently supported only for these programs
PROGRAMS_WITH_FOLDER_DETECTION = {"xtb", "crest"}


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
        - `orientations_pbc`, `energies`, and
        `forces` are indexed by structure;
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


def replace_word(text, old_word, new_word, case_sensitive=True):
    """
    Replace a word in text using word boundary matching.

    Replaces all occurrences of `old_word` with `new_word` using `\b`
    word-boundary matching to ensure only complete words are replaced,
    not partial matches within other words.

    Args:
        text (str): Input text to process.
        old_word (str): Word to replace.
        new_word (str): Replacement word.
        case_sensitive (bool, optional): If `False`, replace all occurrences of
            `old_word` with `new_word` in a case-insensitive manner.
            Default is `True`.

    Returns:
        str: Text with the word replaced.

    Example:
        >>> replace_word("gen test noeigentest","gen","def2svp",case_sensitive=True)
        'def2svp test noeigentest'
        >>> replace_word("opt=(gen) freq","gen","6-31g",case_sensitive=False)
        'opt=(6-31g) freq'
        >>> replace_word("Gen gen GEN","gen","X")
        # -> "X X X"
        >>> replace_word("Gen gen GEN","gen","X",case_sensitive=True)
        # -> "Gen X GEN"
    """
    pattern = r"\b" + re.escape(old_word) + r"\b"
    flags = 0 if case_sensitive else re.IGNORECASE
    return re.sub(pattern, new_word, text, flags=flags)


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
      min_floats=1 -> require at least this
      many float tokens after the integer
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


def match_outfile_pattern(line):
    """
    Match a line of text to known quantum chemistry program signatures.

    Args:
        line (str): Line from an output file.

    Returns:
        str | None: Program name ("gaussian", "orca", "xtb", "crest")
        if matched, else None.
    """
    for program, info in PROGRAM_INFO.items():
        if any(keyword in line for keyword in info["keywords"]):
            return program
    return None


def get_program_type_from_file(filepath):
    """
    Detect the type of quantum chemistry output file.

    Reads only the first 200 lines and scans for characteristic keywords
    of major QC packages to improve efficiency on large output files.

    Args:
        filepath (str): Path to the quantum chemistry output file.

    Returns:
        str: Program name ("gaussian", "orca", "xtb", "crest") or "unknown"
             if the format cannot be detected.
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


def load_molecules_from_paths(
    file_paths,
    index,
    add_index_suffix_for_single=False,
    check_exists=False,
):
    """
    Load molecules from a list of file paths,
    assigning unique names to each molecule.

    For each file in `file_paths`, this function
    loads one or more molecular structures
    using `Molecule.from_filepath`, assigns a
    unique name to each molecule based on the
    file name and structure index, and returns a list of all loaded molecules.

    Args:
        file_paths (list of str or Path): List
        of file paths to load molecules from.
        index (int or str or None): Index or slice
        to select specific structures from each file.
            If None, defaults to "-1" (last structure).
        add_index_suffix_for_single (bool, optional):
        If True, appends an index suffix to
            the molecule name even if only a
            single structure is loaded from a file.
        check_exists (bool, optional): If True,
        checks that each file exists before loading.

    Returns:
        list of Molecule: List of loaded Molecule
        objects, each with a unique name.

    Raises:
        FileNotFoundError: If `check_exists` is True and a file does not exist.
        Exception: If an error occurs during molecule loading from a file.
    """
    # Default index to "-1" (last structure) if not specified
    if index is None:
        index = "-1"

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

            # assign unique names per-structure
            # when file contains multiple structures
            base = os.path.splitext(os.path.basename(file_path))[0]
            if isinstance(mols, list) and len(mols) > 1:
                for j, mol in enumerate(mols, start=1):
                    mol.name = f"{base}_{j}"
            else:
                # Optional suffix for single-structure
                # files (filenames branch).
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


def select_items_by_index(
    items_list,
    index_spec,
    allow_duplicates=False,
    allow_out_of_range=False,
):
    """
    Select items from a list based on an index specification.

    This is a general-purpose utility for applying parse_index_specification
    results to any list.

    Args:
        items_list (list): List of items to select from.
        index_spec (int or str or slice or
        None): Index specification for selection.
            If None or ":", returns all items.
            If int: Direct integer index (0-based Python indexing).
            If slice: Direct slice object (0-based Python indexing).
            If str: String specification (1-based
            indexing, parsed by parse_index_specification).
        allow_duplicates (bool, optional): If True, allows duplicate indices.
            Only applies to string specifications.
        allow_out_of_range (bool, optional):
        If True, allows out-of-range indices.
            Only applies to string specifications.

    Returns:
        list: List of selected items.

    Raises:
        IndexError: If int index is out of range.
        ValueError: If string index specification is invalid or out of range
            (when allow_out_of_range=False).

    Note:
        - int indices raise IndexError for out-of-range access.
        - slice indices never raise errors;
          they return empty or partial results.
        - str indices raise ValueError based on allow_out_of_range parameter.
    """
    # If no filtering needed, return all
    if index_spec is None or index_spec == ":":
        return list(items_list)

    # Handle int and slice types directly
    if isinstance(index_spec, int):
        # Direct integer index - return as single-item list
        return [items_list[index_spec]]
    elif isinstance(index_spec, slice):
        # Direct slice - return sliced items as list
        return items_list[index_spec]

    # Handle string specifications via parse_index_specification
    from chemsmart.utils.utils import parse_index_specification

    selected_indices = parse_index_specification(
        index_spec,
        total_count=len(items_list),
        allow_duplicates=allow_duplicates,
        allow_out_of_range=allow_out_of_range,
    )

    # Handle different return types from parse_index_specification
    if isinstance(selected_indices, list):
        return [items_list[i] for i in selected_indices]
    elif isinstance(selected_indices, int):
        return [items_list[selected_indices]]
    elif isinstance(selected_indices, slice):
        return items_list[selected_indices]
    else:
        raise ValueError(f"Unexpected index type: {type(selected_indices)}")


def clean_label(label: str) -> str:
    """
    Make a label that is safe for filenames, DB keys, RST labels, etc.
    - Keeps letters, digits, `_` and `-`.
    - Replaces spaces, commas, dots, parentheses, etc. with `_`.
    - Encodes "'" as `_prime_` and "*" as `_star_`.
    - Collapses multiple underscores and strips leading/trailing `_`.
    """

    # Preserve your special semantics
    label = label.replace("'", "_prime_")
    label = label.replace("*", "_star_")

    out = []
    for ch in label:
        if ch in SAFE_CHARS:
            out.append(ch)
        else:
            # includes ch.isspace() or ch in {",",
            # ".", "(", ")", "[", "]", "/", "\\"}
            # drop any other weird character, or map to "_"
            out.append("_")

    cleaned = "".join(out)
    # collapse runs of underscores
    cleaned = re.sub(r"_+", "_", cleaned)
    # strip leading/trailing underscores
    cleaned = cleaned.strip("_")
    return cleaned


def convert_string_indices_to_pymol_id_indices(string_indices: str) -> str:
    """
    Convert a comma-separated list of atom
    index ranges into a PyMOL `id` selection.

    The input is expected to be a string of
    indices and/or index ranges separated
    by commas, e.g.:

        "1-10,11,14,19-30"

    This will be converted into a PyMOL selection string where each element is
    prefixed with `id` and combined with `or`, e.g.:

        "id 1-10 or id 11 or id 14 or id 19-30"

    Note: PyMOL selection:
    `select mysel, id 1 or id 2 or id 8-10`
    selects all atoms where (id == 1) OR (id == 2) OR (id is in 8-10)
    So any atom that satisfies any one of those
    conditions is included in the selection.
    That gives you atoms 1, 2, 8, 9, 10.
    This is proper boolean logic:
    or -> set union (combine atoms from all conditions)

    Parameters
    ----------
    string_indices : str
        Comma-separated atom indices and/or ranges, as understood by PyMOL.

    Returns
    -------
    str
        A PyMOL selection string using `id` and `or`.

    Raises
    ------
    ValueError
        If the input string is empty or contains
        no valid indices after stripping.
    """
    # Split on commas and normalise whitespace.
    parts = [
        part.strip() for part in string_indices.split(",") if part.strip()
    ]

    if not parts:
        raise ValueError(
            "string_indices must contain at least one index or range."
        )

    if len(parts) == 1:
        return f"id {parts[0]}"

    return " or ".join(f"id {part}" for part in parts)


def obtain_mols_from_cdx_via_obabel(filename: str) -> List[Chem.Mol]:
    """
    Use the Open Babel CLI ('obabel') to convert a CDX file to SDF and
    return a list of RDKit Mol objects.

    This implementation writes no temporary files; it streams SDF from
    stdout into RDKit, which avoids Windows path issues.

    Args:
        filename: Path to the .cdx file.

    Returns:
        List of RDKit Mol objects.

    Raises:
        ValueError: If 'obabel' is not available or no molecules can be read.
        RuntimeError: If the obabel subprocess fails.
    """
    obabel = shutil.which("obabel")
    if obabel is None:
        raise ValueError(
            "Open Babel CLI ('obabel') is not available on PATH. "
            "Install Open Babel or save the ChemDraw file as CDXML instead."
        )

    # Run: obabel -icdx input.cdx -osdf
    # This writes SDF directly to stdout.
    result = subprocess.run(
        [obabel, "-icdx", filename, "-osdf"],
        check=False,
        capture_output=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"obabel failed to convert {filename!r} to SDF "
            f"(exit code {result.returncode}). stderr:\n{result.stderr.decode(errors='replace')}"
        )

    # Feed stdout bytes directly into RDKit's ForwardSDMolSupplier
    # Use sanitize=False to avoid kekulization errors for organometallic complexes
    sdf_stream = BytesIO(result.stdout)
    suppl = Chem.ForwardSDMolSupplier(
        sdf_stream, sanitize=False, removeHs=False
    )

    mols = [mol for mol in suppl if mol is not None]

    if not mols:
        raise ValueError(
            f"Open Babel produced no valid molecules from CDX file: {filename}"
        )

    return mols


def safe_sanitize(mol, skip_kekulize=False):
    """
    Safely sanitize an RDKit molecule, handling organometallic complexes.

    For organometallic/aromatic-metal complexes, RDKit's kekulization often
    fails because aromatic rings coordinated to metals cannot be properly
    kekulized. This function can skip kekulization for such molecules.

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule to sanitize.
        skip_kekulize (bool): If True, skip kekulization step. Default False.

    Returns:
        rdkit.Chem.Mol: Sanitized molecule.

    Raises:
        Exception: If sanitization fails.
    """
    if skip_kekulize:
        # Skip kekulization for organometallic complexes
        ops = (
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
        Chem.SanitizeMol(mol, sanitizeOps=ops)
        return mol

    try:
        Chem.SanitizeMol(mol)
        return mol
    except Exception as e:
        # Common for organometallic/aromatic-metal complexes:
        # keep most sanitation but skip kekulization
        logger.debug(
            f"Standard sanitization failed ({e}), retrying without kekulization "
            "(common for organometallic complexes)."
        )
        ops = (
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        )
        Chem.SanitizeMol(mol, sanitizeOps=ops)
        return mol


def normalize_metal_bonds(mol):
    """
    Remove aromatic flags from bonds involving metal atoms.

    RDKit does not support aromatic bonds to metal atoms. This function
    identifies bonds to metals and converts any aromatic bonds to single
    bonds, which prevents kekulization and sanitization errors.

    Args:
        mol (rdkit.Chem.Mol): RDKit molecule to normalize.

    Returns:
        rdkit.Chem.Mol: Molecule with normalized metal bonds.

    Notes:
        - Identifies metals using a comprehensive classification based on
          the periodic table (excludes non-metals and metalloids).
        - Converts aromatic bonds to metals to single bonds.
        - Uses element classifications from chemsmart.utils.periodictable.
    """
    from chemsmart.utils.periodictable import NON_METALS_AND_METALLOIDS

    # Identify metal atom indices
    # Metals are elements that are neither non-metals nor metalloids
    metal_idxs = {
        a.GetIdx()
        for a in mol.GetAtoms()
        if a.GetAtomicNum() not in NON_METALS_AND_METALLOIDS
    }

    # Clear aromatic flags from bonds to metals
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()

        if begin_idx in metal_idxs or end_idx in metal_idxs:
            if bond.GetIsAromatic():
                bond.SetIsAromatic(False)
            if bond.GetBondType() == Chem.BondType.AROMATIC:
                bond.SetBondType(Chem.BondType.SINGLE)

    return mol


def fix_cyclopentadienyl_aromaticity(mol: Chem.Mol) -> Chem.Mol:
    """
    RDKit cannot sanitize a neutral aromatic 5-member carbon ring (c1cccc1).
    ChemDraw often uses that for Cp. Convert each such ring into a Cp- by
    making one atom [cH-] (formal charge -1 + explicit H).

    Also de-aromatizes the ring to prevent further processing by attach_one_bond_per_cp_ring.
    """
    rw = Chem.RWMol(mol)
    ri = mol.GetRingInfo()

    for ring in ri.AtomRings():
        if len(ring) != 5:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetSymbol() == "C" and a.GetIsAromatic() for a in atoms):
            continue
        # "Cp as aromatic": 5 aromatic C, each degree 2, each has 1 H, all charges 0
        if not all(
            a.GetDegree() == 2
            and a.GetTotalNumHs() == 1
            and a.GetFormalCharge() == 0
            for a in atoms
        ):
            continue

        i0 = ring[0]
        a0 = rw.GetAtomWithIdx(i0)
        a0.SetFormalCharge(-1)
        a0.SetNumExplicitHs(1)  # important: keep it as [cH-], not [c-]
        a0.SetNoImplicit(True)  # prevent RDKit from dropping that H

        # De-aromatize the ring atoms and bonds to prevent attach_one_bond_per_cp_ring
        # from processing this ring again
        for idx in ring:
            rw.GetAtomWithIdx(idx).SetIsAromatic(False)

        # De-aromatize the ring bonds
        for i in range(len(ring)):
            a = ring[i]
            b = ring[(i + 1) % len(ring)]
            bond = rw.GetBondBetweenAtoms(a, b)
            if bond:
                bond.SetIsAromatic(False)
                bond.SetBondType(Chem.BondType.SINGLE)

    return rw.GetMol()


def _order_ring_atoms_by_walk(
    m: Chem.Mol, ring: tuple[int, ...]
) -> list[int] | None:
    """Return ring atoms in cyclic order by walking ring neighbors."""
    ring_set = set(ring)
    start = ring[0]
    # ring neighbors within ring
    nbrs0 = [
        n.GetIdx()
        for n in m.GetAtomWithIdx(start).GetNeighbors()
        if n.GetIdx() in ring_set
    ]
    if len(nbrs0) != 2:
        return None

    order = [start, nbrs0[0]]
    prev, cur = start, nbrs0[0]
    while True:
        nbrs = [
            n.GetIdx()
            for n in m.GetAtomWithIdx(cur).GetNeighbors()
            if n.GetIdx() in ring_set
        ]
        if len(nbrs) != 2:
            return None
        nxt = nbrs[0] if nbrs[1] == prev else nbrs[1]
        if nxt == start:
            break
        if nxt in order:
            return None
        order.append(nxt)
        prev, cur = cur, nxt

    if len(order) != len(ring):
        return None
    return order


def attach_one_bond_per_cp_ring(
    mol: Chem.Mol, metal_idxs: set[int]
) -> Chem.Mol:
    """
    Approximate η5 Cp coordination by:
      - finding aromatic 5-member all-carbon rings (Cp drawn aromatic)
      - de-aromatizing ring to alternating single/double bonds (kekulizable)
      - adding ONE single bond from a metal to an anchor carbon in each ring (if none exists)
    """
    if not metal_idxs:
        return mol

    rw = Chem.RWMol(mol)
    Chem.GetSymmSSSR(rw)  # ensure rings are perceived

    metal_idx = min(metal_idxs)  # heuristic

    ring_info = rw.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue

        atoms = [rw.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetSymbol() == "C" and a.GetIsAromatic() for a in atoms):
            continue

        # reconstruct cyclic order robustly
        ordered = _order_ring_atoms_by_walk(rw, ring)
        if ordered is None:
            continue

        # already connected to any metal?
        already_connected = any(
            nb.GetIdx() in metal_idxs
            for i in ordered
            for nb in rw.GetAtomWithIdx(i).GetNeighbors()
        )

        # de-aromatize ring atoms
        for i in ordered:
            rw.GetAtomWithIdx(i).SetIsAromatic(False)

        # de-aromatize ring bonds and set alternating pattern
        cycle_bonds = []
        for i in range(5):
            a = ordered[i]
            b = ordered[(i + 1) % 5]
            bond = rw.GetBondBetweenAtoms(a, b)
            if bond is None:
                cycle_bonds = []
                break
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)
            cycle_bonds.append(bond)

        if not cycle_bonds:
            continue

        # set 2 double bonds (cyclopentadiene-like)
        cycle_bonds[0].SetBondType(Chem.BondType.DOUBLE)
        cycle_bonds[2].SetBondType(Chem.BondType.DOUBLE)

        # add one metal-anchor bond if needed
        if not already_connected:
            anchor = None
            for i in ordered:
                if rw.GetAtomWithIdx(i).GetTotalNumHs() > 0:
                    anchor = i
                    break
            if anchor is None:
                anchor = ordered[0]
            rw.AddBond(metal_idx, anchor, Chem.BondType.SINGLE)

    new_mol = rw.GetMol()
    new_mol.UpdatePropertyCache(strict=False)
    return new_mol


def update_shell_config(shell_file: Path, env_vars: list) -> None:
    """
    Append chemsmart ``export`` lines to *shell_file* (idempotent).

    This is the POSIX equivalent of writing ``$env:PATH`` lines to a
    PowerShell profile or updating the Windows registry.  The update is
    idempotent: if the marker comment ``# Added by chemsmart installer`` is
    already present the file is not modified again.

    Creates the file if it does not yet exist.

    Args:
        shell_file: Path to the shell startup file (e.g. ``~/.bashrc``).
        env_vars: List of ``export VAR=...`` lines to append.
    """
    if not shell_file.exists():
        shell_file.touch()

    with shell_file.open("r+", encoding="utf-8") as f:
        lines = f.readlines()
        if not any("Added by chemsmart installer" in line for line in lines):
            f.write("\n# Added by chemsmart installer\n")
            for var in env_vars:
                f.write(f"{var}\n")
            f.write("\n")
            logger.info(f"Updated shell config: {shell_file}")
        else:
            logger.info(f"Shell config already updated: {shell_file}")

    logger.info(f"Please restart your terminal or run 'source {shell_file}'.")


_PS_BLOCK_START = "# >>> chemsmart initialize >>>"
_PS_BLOCK_END = "# <<< chemsmart initialize <<<"
# Legacy marker written by earlier versions of the installer.
_PS_BLOCK_LEGACY = "# Added by chemsmart installer"


def update_powershell_profiles(profiles: list, ps_env_vars: list) -> None:
    """
    Write chemsmart initialisation lines to each PowerShell profile,
    replacing any previously written block.  Creates profile directories
    and files as needed.

    The block is delimited by ``# >>> chemsmart initialize >>>`` /
    ``# <<< chemsmart initialize <<<`` markers.  If an older block using
    the legacy marker ``# Added by chemsmart installer`` is found it is
    removed and replaced with the current block so that users who ran
    ``make configure`` with an earlier version are migrated automatically.

    Args:
        profiles: List of :class:`~pathlib.Path` objects pointing to the PS
            profile files to update.
        ps_env_vars: List of PowerShell assignment / alias lines to write
            inside the block.
    """
    new_block = (
        f"\n{_PS_BLOCK_START}\n"
        + "\n".join(ps_env_vars)
        + f"\n{_PS_BLOCK_END}\n"
    )

    for ps_profile in profiles:
        ps_profile.parent.mkdir(parents=True, exist_ok=True)
        if not ps_profile.exists():
            ps_profile.write_text(new_block, encoding="utf-8")
            logger.info(f"Created PowerShell profile: {ps_profile}")
            continue

        content = ps_profile.read_text(encoding="utf-8")

        # ── Remove existing chemsmart block (new-style markers) ──────────
        if _PS_BLOCK_START in content:
            start = content.find(_PS_BLOCK_START)
            end = content.find(_PS_BLOCK_END, start)
            if end != -1:
                content = content[:start] + content[end + len(_PS_BLOCK_END) :]
            else:
                # Malformed block — remove from start marker to end of file
                content = content[:start]

        # ── Remove legacy block ("# Added by chemsmart installer") ───────
        if _PS_BLOCK_LEGACY in content:
            lines = content.splitlines(keepends=True)
            filtered = []
            in_legacy = False
            for line in lines:
                if _PS_BLOCK_LEGACY in line:
                    in_legacy = True
                    continue
                if in_legacy and line.strip() == "":
                    in_legacy = False
                    continue
                if not in_legacy:
                    filtered.append(line)
            content = "".join(filtered)

        # Ensure exactly one blank line separates existing content from the
        # new block (strip trailing newlines, then add exactly one).
        ps_profile.write_text(
            content.rstrip("\n") + "\n" + new_block, encoding="utf-8"
        )
        logger.info(f"Updated PowerShell profile: {ps_profile}")

    logger.info(
        "PowerShell profiles updated.\n"
        "To apply changes in the current PowerShell session, run:\n"
        "  . $PROFILE"
    )


def update_windows_env(paths_to_add: list, pythonpath_entry: str) -> None:
    """
    Add directories to the Windows user PATH and PYTHONPATH via the registry.

    This is the Windows equivalent of appending ``export PATH=...`` lines to
    ``~/.bashrc``.  Changes take effect in **new** terminal sessions; users
    must restart their terminal or execute the PowerShell one-liner below to
    refresh the current session::

        $env:PATH = [System.Environment]::GetEnvironmentVariable('PATH', 'Machine') + ';' +
                    [System.Environment]::GetEnvironmentVariable('PATH', 'User')

    Args:
        paths_to_add: List of directory paths to append to the user PATH (only
            entries that are not already present will be added).
        pythonpath_entry: Single directory to append to the user PYTHONPATH (no-op
            if already present).
    """
    try:
        import winreg  # noqa: PLC0415 — Windows-only stdlib module
    except ImportError:
        logger.warning("winreg not available; skipping Windows PATH update.")
        return

    try:
        with winreg.OpenKey(
            winreg.HKEY_CURRENT_USER,
            "Environment",
            0,
            winreg.KEY_READ | winreg.KEY_WRITE,
        ) as key:
            # ---- PATH ----
            try:
                current_path, _ = winreg.QueryValueEx(key, "PATH")
            except FileNotFoundError:
                current_path = ""

            path_parts = [p for p in current_path.split(";") if p.strip()]
            new_paths = [p for p in paths_to_add if p not in path_parts]
            if new_paths:
                winreg.SetValueEx(
                    key,
                    "PATH",
                    0,
                    winreg.REG_EXPAND_SZ,
                    ";".join(path_parts + new_paths),
                )
                logger.info(f"Added to Windows user PATH: {new_paths}")
            else:
                logger.info(
                    "Windows user PATH already contains all required paths."
                )

            # ---- PYTHONPATH ----
            try:
                current_pypath, _ = winreg.QueryValueEx(key, "PYTHONPATH")
            except FileNotFoundError:
                current_pypath = ""

            pypath_parts = [p for p in current_pypath.split(";") if p.strip()]
            if pythonpath_entry not in pypath_parts:
                winreg.SetValueEx(
                    key,
                    "PYTHONPATH",
                    0,
                    winreg.REG_SZ,
                    ";".join(pypath_parts + [pythonpath_entry]),
                )
                logger.info(
                    f"Updated Windows PYTHONPATH to include: {pythonpath_entry}"
                )
            else:
                logger.info(
                    "Windows PYTHONPATH already contains the required path."
                )

        # Notify running applications about the environment change so
        # that new terminal windows inherit the updated PATH immediately.
        import ctypes

        ctypes.windll.user32.SendMessageTimeoutW(
            0xFFFF,  # HWND_BROADCAST
            0x001A,  # WM_SETTINGCHANGE
            0,
            "Environment",
            0x0002,  # SMTO_ABORTIFHUNG
            5000,
            None,
        )
        logger.info("Broadcasted environment change to Windows.")

    except PermissionError:
        logger.warning(
            "Permission denied accessing Windows registry. "
            "Run as administrator or add these paths to PATH manually:\n"
            + "\n".join(f"  {p}" for p in paths_to_add)
        )
    except Exception as e:
        logger.warning(f"Could not update Windows environment: {e}")
