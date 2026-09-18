"""Distance analysis of recorded Gaussian IRC geometries.
"""

import os

import numpy as np

from chemsmart.io.molecules.structure import Molecule


def single_point_distances(filepath, index):
    """Return labeled distances and a square matrix for one IRC geometry. 
    Return a dictionary with source_file, geometry_index, molecule, matrix,atom_pairs and distances. 
    """
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)):
        raise TypeError("Geometry index must be an integer.")
    if index < 1:
        raise ValueError("Geometry index must be positive and 1-based.")
    filepath = os.path.abspath(os.fspath(filepath))
    suffix = os.path.splitext(filepath)[1].lower()
    if suffix not in {".log", ".out"}:
        raise ValueError(
            "Single-point analysis supports Gaussian IRC output only."
        )
    if suffix == ".out":
        from chemsmart.utils.io import get_program_type_from_file

        if get_program_type_from_file(filepath) != "gaussian":
            raise ValueError("Single-point IRC analysis is Gaussian-only.")
    molecule = Molecule.from_irc_filepath(filepath, index=str(index))
    if molecule is None or "irc_point" not in molecule.info:
        raise ValueError("No recorded Gaussian IRC geometry was returned.")
    matrix = molecule.atom_distances()
    first, second = np.triu_indices(molecule.num_atoms, k=1)
    pairs = [(int(i + 1), int(j + 1)) for i, j in zip(first, second)]
    return {
        "source_file": filepath,
        "geometry_index": int(index),
        "molecule": molecule,
        "matrix": matrix,
        "atom_pairs": pairs,
        "distances": matrix[first, second],
    }


def compare_distances(first, second):
    """Compare distances of two Molecules in the requested A, B order.

    Return molecules, atom_pairs, distance_matrix, signed_change,
    absolute_change, comparison_matrix, comparison_columns and distance_units.
    Positive change means distance increases; negative change means it decreases.
    """
    distances = Molecule.all_distances((first, second))
    i, j = np.triu_indices(first.num_atoms, k=1)
    pairs = [(int(a + 1), int(b + 1)) for a, b in zip(i, j)]
    signed = distances[1] - distances[0]
    absolute = np.abs(signed)
    return {
        "molecules": (first, second),
        "atom_pairs": pairs,
        "distance_matrix": distances,
        "signed_change": signed,
        "absolute_change": absolute,
        "comparison_matrix": np.column_stack(
            (distances[0], distances[1], signed, absolute)
        ),
        "comparison_columns": (
            "distance_A",
            "distance_B",
            "signed_change",
            "absolute_change",
        ),
        "distance_units": "angstrom",
    }


def two_point_distances(filepath, first_index, second_index):
    """Compare two manually selected recorded points in one Gaussian IRC.

    Preserve requested order: first is A, second is B.
    """
    filepath = os.path.abspath(os.fspath(filepath))
    first, second = Molecule.select_geometries(
        filepath, first_index, second_index
    )
    result = compare_distances(first, second)
    result["source_file"] = filepath
    result["geometry_indices"] = (int(first_index), int(second_index))
    return result


def distance_profile(molecules):
    """Summarize every atom-pair distance over an ordered geometry sequence.
    At least two geometries are required. Atom correspondence is assumed.
    """
    molecules = tuple(molecules)
    if len(molecules) < 2:
        raise ValueError("A distance profile needs at least two geometries.")
    distances = Molecule.all_distances(molecules)
    first, second = np.triu_indices(molecules[0].num_atoms, k=1)
    pairs = [(int(i + 1), int(j + 1)) for i, j in zip(first, second)]
    signed = distances[-1] - distances[0]
    minimum = distances.min(axis=0)
    maximum = distances.max(axis=0)
    return {
        "molecules": molecules,
        "atom_pairs": pairs,
        "distance_matrix": distances,
        "summary_matrix": np.column_stack(
            (
                distances[0],
                distances[-1],
                signed,
                np.abs(signed),
                minimum,
                maximum,
                maximum - minimum,
            )
        ),
        "summary_columns": (
            "distance_start",
            "distance_end",
            "signed_start_to_end_change",
            "absolute_start_to_end_change",
            "minimum_distance",
            "maximum_distance",
            "distance_range",
        ),
        "minimum_indices": distances.argmin(axis=0) + 1,
        "maximum_indices": distances.argmax(axis=0) + 1,
        "distance_units": "angstrom",
    }


def irc_distance_profile(filepath):
    """Return a distance profile for all recorded points of one IRC branch.
    Start/end mean first/last recorded geometry.
    """
    filepath = os.path.abspath(os.fspath(filepath))
    suffix = os.path.splitext(filepath)[1].lower()
    if suffix not in {".log", ".out"}:
        raise ValueError("IRC distance profiles support Gaussian output only.")
    if suffix == ".out":
        from chemsmart.utils.io import get_program_type_from_file

        if get_program_type_from_file(filepath) != "gaussian":
            raise ValueError("IRC distance profiles are Gaussian-only.")
    molecules = Molecule.from_irc_filepath(filepath, index=":")
    if not molecules or any(
        "irc_point" not in molecule.info or "irc_path" not in molecule.info
        for molecule in molecules
    ):
        raise ValueError("No recorded Gaussian IRC geometries were returned.")
    if len({molecule.info["irc_path"] for molecule in molecules}) != 1:
        raise ValueError("Use a single-branch IRC log, not multiple paths.")
    result = distance_profile(molecules)
    result["source_file"] = filepath
    return result


def _filter_connected_distance_profile(profile, bond_cutoff_buffer, adjust_H):
    """Select pairs connected at least once.
    Connectivity is estimated by the existing Molecule.to_graph.
    """
    if isinstance(bond_cutoff_buffer, (bool, np.bool_)) or not isinstance(
        bond_cutoff_buffer, (int, float, np.integer, np.floating)
    ):
        raise TypeError("Bond cutoff buffer must be a number in angstroms.")
    if not np.isfinite(bond_cutoff_buffer) or bond_cutoff_buffer < 0:
        raise ValueError("Bond cutoff buffer must be finite and nonnegative.")
    if not isinstance(adjust_H, (bool, np.bool_)):
        raise TypeError("adjust_H must be a boolean.")
    molecules = profile["molecules"]
    pairs = profile["atom_pairs"]
    pair_columns = {pair: column for column, pair in enumerate(pairs)}
    connected = np.zeros((len(molecules), len(pairs)), dtype=bool)
    for row, molecule in enumerate(molecules):
        graph = molecule.to_graph(
            bond_cutoff_buffer=float(bond_cutoff_buffer),
            adjust_H=bool(adjust_H),
        )
        for first, second in graph.edges:
            pair = tuple(sorted((int(first) + 1, int(second) + 1)))
            connected[row, pair_columns[pair]] = True
    retained = connected.any(axis=0)
    result = dict(profile)
    result["all_distance_profile"] = profile
    result["all_pair_connected_matrix"] = connected
    result["retained_pair_mask"] = retained
    result["atom_pairs"] = [
        pair for pair, keep in zip(pairs, retained) if keep
    ]
    result["distance_matrix"] = profile["distance_matrix"][:, retained]
    result["summary_matrix"] = profile["summary_matrix"][retained]
    result["minimum_indices"] = profile["minimum_indices"][retained]
    result["maximum_indices"] = profile["maximum_indices"][retained]
    result["connected_matrix"] = connected[:, retained]
    result["connected_counts"] = result["connected_matrix"].sum(axis=0)
    result["connectivity_settings"] = {
        "method": "Molecule.to_graph",
        "bond_cutoff_buffer_angstrom": float(bond_cutoff_buffer),
        "adjust_H": bool(adjust_H),
        "selection_rule": "connected in at least one supplied geometry",
    }
    return result


def connected_distance_profile(molecules, bond_cutoff_buffer=0.05, adjust_H=True):
    """
    Summarize pairs estimated to be connected anywhere in the sequence.
    This is an estimated-connectivity screen, not a CV quality score or proof of bonds.
    """
    profile = distance_profile(molecules)
    return _filter_connected_distance_profile(
        profile, bond_cutoff_buffer, adjust_H
    )


def irc_connected_distance_profile(filepath, bond_cutoff_buffer=0.05, adjust_H=True):
    """Connectivity-filtered distances for one recorded Gaussian IRC branch.

    Reuse irc_distance_profile, retaining its source and termination metadata.
    The full unfiltered profile remains accessible as all_distance_profile.
    """
    profile = irc_distance_profile(filepath)
    return _filter_connected_distance_profile(
        profile, bond_cutoff_buffer, adjust_H
    )
