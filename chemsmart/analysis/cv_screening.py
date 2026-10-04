"""Collective-variable analysis of ordered Molecules using one analysis object.

Notes for downstream exports:
  *_summary files: summaries for later collective-variable screening and
  selection.
  *_values files: all individual values for each Molecule, for the atom
  combinations included in the profile. Rows represent Molecules/geometries;
  columns identify the individual distances, angles or dihedrals.
"""

import os
import re

import numpy as np

from chemsmart.io.molecules.structure import Molecule


class CVScreening:
    """Analyse candidate collective variables across ordered Molecules.
    Results belong to one analysis object.
    call refresh() after changing a supplied Molecule's coordinates.
    """

    def __init__(self, molecules):
        self.molecules = tuple(molecules)
        self._metadata = {}
        self.refresh()

    @property
    def profile(self):
        """Return the complete calculated data dictionary."""
        return self._profile

    @property
    def distances(self):
        """Distance values and summaries in angstroms."""
        return self._profile["distances"]

    @property
    def angles(self):
        """Angle values and summaries in degrees."""
        return self._profile["angles"]

    @property
    def dihedrals(self):
        """Dihedral values and circular summaries in degrees."""
        return self._profile["dihedrals"]

    @staticmethod
    def _geometry_summary(values, circular=False):
        """Use shortest signed torsion endpoint changes in [-180,180)."""
        summary = np.full((values.shape[1], 7), np.nan)
        counts = np.isfinite(values).sum(axis=0)
        for column in range(values.shape[1]):
            series = values[:, column]
            first, last = series[0], series[-1]
            change = last - first
            if circular:
                change = (change + 180.0) % 360.0 - 180.0
            summary[column, :4] = first, last, change, abs(change)
            finite = series[np.isfinite(series)]
            if not len(finite):
                continue
            if circular:
                ordered = np.sort(finite % 360.0)
                gaps = np.diff(np.r_[ordered, ordered[0] + 360.0])
                gap = int(np.argmax(gaps))
                lower = ordered[(gap + 1) % len(ordered)]
                lower = (lower + 180.0) % 360.0 - 180.0
                span = max(0.0, 360.0 - gaps[gap])
                summary[column, 4:] = lower, lower + span, span
            else:
                summary[column, 4:] = (
                    finite.min(),
                    finite.max(),
                    np.ptp(finite),
                )
        columns = (
            "start",
            "end",
            "signed_change",
            "absolute_change",
            "arc_start" if circular else "minimum",
            "arc_end" if circular else "maximum",
            "circular_range" if circular else "range",
        )
        return summary, columns, counts

    def _angular_profile(self, property_name, order):
        molecules = self.molecules
        tables = [getattr(molecule, property_name) for molecule in molecules]
        detected = [
            {tuple(int(i) for i in row[:order]): row[-1] for row in table}
            for table in tables
        ]
        combinations = sorted(set().union(*(set(rows) for rows in detected)))
        indices = np.asarray(combinations, dtype=int).reshape(-1, order)
        values = np.empty((len(molecules), len(indices)), dtype=float)
        connected = np.zeros(values.shape, dtype=bool)
        for frame, (molecule, rows) in enumerate(zip(molecules, detected)):
            missing = []
            for column, atoms in enumerate(combinations):
                if atoms in rows:
                    values[frame, column] = rows[atoms]
                    connected[frame, column] = True
                else:
                    missing.append(column)
            if missing:
                values[frame, missing] = molecule.geometry_table(
                    indices[missing]
                )[:, -1]
        summary, columns, counts = self._geometry_summary(
            values, circular=order == 4
        )
        return {
            "atom_indices": indices,
            "values_matrix": values,
            "connected_matrix": connected,
            "summary_matrix": summary,
            "summary_columns": columns,
            "valid_counts": counts,
            "units": "degree",
        }

    def refresh(self):
        """Recalculate profiles after coordinate edits and return self."""
        molecules = self.molecules
        if not molecules:
            raise ValueError("Provide at least one geometry.")
        if not all(isinstance(molecule, Molecule) for molecule in molecules):
            raise TypeError("All geometries must be Molecule objects.")
        symbols = list(molecules[0].symbols)
        for molecule in molecules:
            if list(molecule.symbols) != symbols:
                raise ValueError(
                    "Geometries have different element sequences."
                )
            molecule._validate_matrix_positions()
        first, second = np.triu_indices(len(symbols), k=1)
        indices = np.column_stack((first + 1, second + 1))
        values = np.stack(
            [
                molecule.distances_matrix[first, second]
                for molecule in molecules
            ]
        )
        lookup = {tuple(pair): column for column, pair in enumerate(indices)}
        connected = np.zeros(values.shape, dtype=bool)
        for frame, molecule in enumerate(molecules):
            for i, j in molecule.to_graph().edges:
                connected[frame, lookup[tuple(sorted((i + 1, j + 1)))]] = True
        summary, columns, counts = self._geometry_summary(values)
        distances = {
            "atom_indices": indices,
            "values_matrix": values,
            "connected_matrix": connected,
            "selected_mask": connected.any(axis=0),
            "summary_matrix": summary,
            "summary_columns": columns,
            "valid_counts": counts,
            "units": "angstrom",
        }
        self._profile = {
            "molecules": molecules,
            "distances": distances,
            "angles": self._angular_profile("angles_matrix", 3),
            "dihedrals": self._angular_profile("dihedrals_matrix", 4),
            "connectivity_method": "Molecule.to_graph(default settings)",
        }
        self._profile.update(self._metadata)
        return self

    @classmethod
    def from_filepath(cls, filepath):
        """Use the Molecule.from_filepath reader for all geometries."""
        filepath = os.path.abspath(filepath)
        molecules = Molecule.from_filepath(
            filepath, index=":", return_list=True
        )
        if not molecules:
            raise ValueError("No geometries were returned by from_filepath.")
        analysis = cls(molecules)
        result = analysis._metadata
        result["source_file"] = filepath
        result["reader"] = "Molecule.from_filepath"
        result["geometry_index_definition"] = (
            "1-based position in reader output"
        )
        analysis._profile.update(result)
        return analysis

    def select(self, selections, central_bond_only=False):
        """Return separate profile snapshots for nested atom-label groups."""
        if not isinstance(central_bond_only, bool):
            raise TypeError("central_bond_only must be True or False.")
        if not isinstance(selections, (list, tuple)) or not selections:
            raise ValueError("Provide nested groups")
        symbols = list(self.molecules[0].symbols)
        groups = []
        for group in selections:
            if not isinstance(group, (list, tuple)) or len(group) not in (
                1,
                2,
                3,
            ):
                raise ValueError(
                    "Each selection must contain one, two or three labels."
                )
            if central_bond_only and len(group) != 2:
                raise ValueError("central_bond_only requires two-atom groups.")
            atoms = []
            labels = []
            for label in group:
                match = (
                    re.fullmatch(r"([A-Za-z]{1,2})([1-9][0-9]*)", label)
                    if isinstance(label, str)
                    else None
                )
                if match is None:
                    raise ValueError("Use element plus 1-based atom number.")
                atom = int(match[2])
                if atom > len(symbols):
                    raise ValueError(
                        f"Atom {atom} is outside 1..{len(symbols)}."
                    )
                actual = f"{symbols[atom - 1]}{atom}"
                if match[1].lower() != symbols[atom - 1].lower():
                    raise ValueError(
                        f"Selection {label!r} does not match {actual}."
                    )
                atoms.append(atom)
                labels.append(actual)
            if len(set(atoms)) != len(atoms):
                raise ValueError("Each selection must contain distinct atoms.")
            result = {
                key: value
                for key, value in self.profile.items()
                if key not in ("distances", "angles", "dihedrals")
            }
            result["selection_labels"] = tuple(labels)
            result["selection_indices"] = tuple(atoms)
            result["central_bond_only"] = central_bond_only
            result["selection_note"] = (
                "Existing candidates only; all pair distances, but angles and "
                "dihedrals require a connected path in at least one geometry. "
            )
            for category in ("distances", "angles", "dihedrals"):
                data = self.profile[category]
                indices = data["atom_indices"]
                if len(atoms) == 1:
                    mask = (indices == atoms[0]).any(axis=1)
                elif len(atoms) == 3:
                    if category == "distances":
                        mask = np.isin(indices, atoms).all(axis=1)
                    else:
                        mask = np.zeros(len(indices), dtype=bool)
                        for start in range(indices.shape[1] - 2):
                            window = indices[:, start : start + 3]
                            mask |= (window == atoms).all(axis=1)
                            mask |= (window == atoms[::-1]).all(axis=1)
                else:
                    left, right = indices[:, :-1], indices[:, 1:]
                    if category == "dihedrals" and central_bond_only:
                        left, right = indices[:, 1:2], indices[:, 2:3]
                    mask = (
                        ((left == atoms[0]) & (right == atoms[1]))
                        | ((left == atoms[1]) & (right == atoms[0]))
                    ).any(axis=1)
                selected = dict(data)
                for key in (
                    "atom_indices",
                    "summary_matrix",
                    "valid_counts",
                    "selected_mask",
                ):
                    if key in data:
                        selected[key] = data[key][mask].copy()
                for key in ("values_matrix", "connected_matrix"):
                    selected[key] = data[key][:, mask].copy()
                result[category] = selected
            groups.append(result)
        return groups

    def atom_maximum_changes(self, connected_distances=True):
        """Summarise largest absolute endpoint changes involving each atom."""
        profile = self._profile
        categories = ("distances", "angles", "dihedrals")
        symbols = list(profile["molecules"][0].symbols)
        labels = tuple(f"{symbol}{i}" for i, symbol in enumerate(symbols, 1))
        values = np.full((len(symbols), 3), np.nan)
        notes = np.empty((len(symbols), 3), dtype=object)
        candidate_counts = np.zeros(values.shape, dtype=int)
        undefined_counts = np.zeros(values.shape, dtype=int)
        winner_columns = {}
        for category_column, category in enumerate(categories):
            data = profile[category]
            indices = data["atom_indices"]
            changes = data["summary_matrix"][
                :, data["summary_columns"].index("absolute_change")
            ]
            eligible = np.ones(len(indices), dtype=bool)
            if category == "distances" and connected_distances:
                eligible &= data["selected_mask"]
            winners_by_atom = []
            for atom in range(1, len(symbols) + 1):
                candidates = eligible & (indices == atom).any(axis=1)
                finite = candidates & np.isfinite(changes)
                row = atom - 1
                candidate_counts[row, category_column] = candidates.sum()
                undefined_counts[row, category_column] = (
                    candidates & ~finite
                ).sum()
                if not finite.any():
                    notes[row, category_column] = (
                        "Undefined endpoint changes"
                        if candidates.any()
                        else "No eligible coordinate"
                    )
                    winners_by_atom.append(())
                    continue
                maximum = changes[finite].max()
                winners = sorted(
                    np.flatnonzero(finite & (changes == maximum)),
                    key=lambda column: tuple(indices[column]),
                )
                values[row, category_column] = maximum
                notes[row, category_column] = "; ".join(
                    "-".join(labels[i - 1] for i in indices[column])
                    for column in winners
                )
                winners_by_atom.append(tuple(int(i) for i in winners))
            winner_columns[category] = tuple(winners_by_atom)
        return {
            "atom_labels": labels,
            "categories": categories,
            "units": ("angstrom", "degree", "degree"),
            "values_matrix": values,
            "notes_matrix": notes,
            "winner_columns": winner_columns,
            "candidate_counts": candidate_counts,
            "undefined_counts": undefined_counts,
            "distance_selection": (
                "connected_any" if connected_distances else "all"
            ),
        }
