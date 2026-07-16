"""
Database export module for exporting records from a chemsmart database.

Output format is inferred from the file extension of the output path:
* .json   – Full structured records.
* .csv    – Tabular scalar properties.
* .xyz    – Cartesian coordinates of selected structure(s).
* .extxyz – Extended XYZ with per-frame energy and per-atom forces.
"""

import csv
import json
import logging
import os

import numpy as np

from chemsmart.database.database import Database
from chemsmart.database.utils import (
    collect_energies_for_structure,
    convert_numpy,
    sort_frames_by_energy,
)

logger = logging.getLogger(__name__)


# Default columns always present in CSV output.
CSV_DEFAULT_COLUMNS = ["record_index", "record_id", "chemical_formula"]

# Optional scalar keys supported via -k/--keys.
CSV_OPTIONAL_COLUMNS = {
    "program",
    "method",
    "basis",
    "charge",
    "multiplicity",
    "smiles",
    "total_energy",
    "homo_energy",
    "lumo_energy",
    "fmo_gap",
    "zero_point_energy",
    "enthalpy",
    "entropy",
    "gibbs_free_energy",
}

SUPPORTED_FORMATS = {".json", ".csv", ".xyz", ".extxyz"}


class DatabaseExporter:
    """Export a chemsmart database to JSON, CSV, XYZ, or extended XYZ."""

    def __init__(
        self,
        db_file,
        output,
        record_index=None,
        record_id=None,
        structure_index=None,
        structure_id=None,
        molecule_id=None,
        keys=None,
        method=None,
        basis=None,
    ):
        self.db_file = db_file
        self.db = Database(db_file)
        self.output = output
        self.record_index = record_index
        self.record_id = record_id
        self.structure_index = structure_index
        self.structure_id = structure_id
        self.molecule_id = molecule_id
        self.keys = keys
        self.method = method
        self.basis = basis
        self.format = self._infer_format()
        self.parsed_keys = self._parse_csv_keys()

    @property
    def user_primary_mb(self):
        """User-specified primary (method, basis), or None if not set."""
        if self.method is None and self.basis is None:
            return None
        return (self.method, self.basis)

    def export(self):
        """Run export, dispatching to the appropriate format handler."""
        handler = {
            ".json": self.to_json,
            ".csv": self.to_csv,
            ".xyz": self.to_xyz,
            ".extxyz": self.to_extxyz,
        }
        handler[self.format]()

    def to_json(self):
        """Export the full database as structured JSON."""
        records = self.db.get_all_records()
        with open(self.output, "w") as f:
            json.dump(convert_numpy(records), f, indent=4)

    def to_csv(self):
        """Export scalar properties of the full database as CSV."""
        records = self.db.get_all_records()
        columns = list(CSV_DEFAULT_COLUMNS)
        if self.parsed_keys:
            columns.extend(self.parsed_keys)
        rows = [self.record_to_csv_row(rec, columns) for rec in records]
        with open(self.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns, restval="NaN")
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def record_to_csv_row(record, columns):
        """Flatten a record dict into a CSV row using the last structure
        as molecule/structure context."""
        meta = record.get("meta", {}) or {}
        results = record.get("results", {}) or {}
        provenance = record.get("provenance", {}) or {}
        molecules = record.get("molecules", []) or []
        last_mol = molecules[-1] if molecules else {}
        lookup = {
            "record_index": record.get("record_index"),
            "record_id": record.get("record_id"),
            "program": provenance.get("program"),
            "method": meta.get("method"),
            "basis": meta.get("basis"),
            **{k: v for k, v in results.items() if v is not None},
            "chemical_formula": last_mol.get("chemical_formula"),
            "charge": last_mol.get("charge"),
            "multiplicity": last_mol.get("multiplicity"),
            "smiles": last_mol.get("smiles"),
        }
        return {col: lookup.get(col, "NaN") for col in columns}

    def to_xyz(self):
        """Export selected structure(s) as XYZ.
        * molecule_id -> all conformers of that molecule.
        * structure_id -> that single structure.
        * record_index / record_id -> structures of the record,
          filtered by structure_index (default last only).
        """
        frames = self._collect_frames(include_forces=False)
        logger.info(
            f"Writing {len(frames)} XYZ frame(s) to "
            f"{os.path.basename(self.output)}"
        )
        with open(self.output, "w") as f:
            for frame in frames:
                self._write_xyz_frame(f, frame)

    def to_extxyz(self):
        """Export selected structure(s) as extended XYZ.
        Same selection semantics as to_xyz, but each frame also
        carries the calculation energy in the comment line and per-atom
        forces as a forces:R:3 property.
        """
        frames = self._collect_frames(include_forces=True)
        logger.info(
            f"Writing {len(frames)} Extended XYZ frame(s) to "
            f"{os.path.basename(self.output)}"
        )
        with open(self.output, "w") as f:
            for frame in frames:
                self._write_extxyz_frame(f, frame)

    def _collect_frames(self, include_forces=False):
        """Validate selectors and dispatch to the appropriate frame builder."""
        selectors = (
            self.record_index,
            self.record_id,
            self.structure_id,
            self.molecule_id,
        )
        if sum(s is not None for s in selectors) != 1:
            raise ValueError(
                "XYZ/extended XYZ export requires exactly one of "
                "--ri/--record-index, --rid/--record-id, "
                "--sid/--structure-id, or --mid/--molecule-id."
            )

        # effective_primary tracks the (method, basis) actually used for
        # energy/forces — either user-specified or auto-picked.
        effective_primary = self.user_primary_mb
        sort_primary = None

        if self.molecule_id is not None:
            frames, sort_primary = self.frames_from_molecule_id(
                include_forces=include_forces
            )
            if effective_primary is None:
                effective_primary = sort_primary
        elif self.structure_id is not None:
            frames = self.frames_from_structure_id(
                include_forces=include_forces
            )
            # For single-structure path, infer effective_primary from the frame.
            if effective_primary is None and frames:
                f0 = frames[0]
                if (
                    f0.get("primary_method") is not None
                    or f0.get("primary_basis") is not None
                ):
                    effective_primary = (
                        f0.get("primary_method"),
                        f0.get("primary_basis"),
                    )
        else:
            frames = self.frames_from_record(include_forces=include_forces)
            # Infer effective_primary from the frame for --ri/--rid path.
            if effective_primary is None and frames:
                f0 = frames[0]
                if (
                    f0.get("primary_method") is not None
                    or f0.get("primary_basis") is not None
                ):
                    effective_primary = (
                        f0.get("primary_method"),
                        f0.get("primary_basis"),
                    )

        if not frames:
            raise ValueError("No structures found for the given selection.")

        if include_forces and effective_primary is None:
            raise ValueError(
                "No structures found with both energy and forces at "
                "the given selection."
            )

        if effective_primary is not None:
            mb_str = "/".join(x for x in effective_primary if x) or "unknown"
            if include_forces:
                # extxyz: keep only structures that have BOTH energy and forces
                kept, skipped = [], []
                for f in frames:
                    if f.get(
                        "primary_energy"
                    ) is not None and self._validate_forces(f["molecule"]):
                        kept.append(f)
                    else:
                        skipped.append(f.get("structure_id"))
                if skipped:
                    max_show = 5
                    shown = ", ".join(
                        (sid[:12] if sid else "?")
                        for sid in skipped[:max_show]
                    )
                    more = (
                        f", ... and {len(skipped) - max_show} more"
                        if len(skipped) > max_show
                        else ""
                    )
                    logger.warning(
                        f"Skipped {len(skipped)} structure(s) without both "
                        f"energy and forces at {mb_str}: {shown}{more}."
                    )
                if not kept:
                    raise ValueError(
                        "No structures found with both energy and forces at "
                        "the given selection."
                    )
                frames = kept
                sort_primary = effective_primary
            elif self.user_primary_mb is not None:
                # xyz: keep only structures that have energy at the requested mb.
                kept = [
                    f
                    for f in frames
                    if any(
                        (m, b) == effective_primary
                        for m, b, _ in f.get("energies", [])
                    )
                ]
                if not kept:
                    raise ValueError(
                        "No structures found at the given method/basis."
                    )
                frames = kept
                sort_primary = effective_primary

        if self.molecule_id is not None:
            frames = sort_frames_by_energy(frames, primary=sort_primary)

        return frames

    def frames_from_record(self, include_forces=False):
        """Build frames from --ri/--rid (+ optional --si)."""
        from chemsmart.io.database import DatabaseFile
        from chemsmart.utils.utils import string2index_1based

        # Resolve record
        if self.record_index is not None:
            record = self.db.get_record(record_index=self.record_index)
            if record is None:
                raise ValueError(
                    f"No record found at index {self.record_index}."
                )
        else:
            full_rid = self.db.get_record_by_partial_id(self.record_id)
            record = self.db.get_record(record_id=full_rid)
            if record is None:
                raise ValueError(
                    f"No record found with ID '{self.record_id}'."
                )

        meta = record.get("meta", {}) or {}
        method = meta.get("method")
        basis = meta.get("basis")
        mol_dicts = record.get("molecules", []) or []

        # Apply structure_index (default: -1 = last)
        si = self.structure_index if self.structure_index is not None else "-1"
        idx = string2index_1based(str(si))
        if isinstance(idx, int):
            try:
                selected = [mol_dicts[idx]]
            except IndexError as exc:
                raise ValueError(
                    f"Structure index {si} out of range "
                    f"(record has {len(mol_dicts)} structures)."
                ) from exc
        else:
            selected = mol_dicts[idx]

        db_file = DatabaseFile(filename=self.db_file)
        primary_method, primary_basis = method, basis
        frames = []
        for mol_dict in selected:
            if include_forces:
                # For --rid path, query forces within the specific record
                forces, energy_from_forces = (
                    self.db.get_forces_for_record_structure_at(
                        record["record_id"],
                        mol_dict.get("structure_id"),
                        primary_method,
                        primary_basis,
                    )
                )
                mol_dict_with_forces = dict(mol_dict)
                mol_dict_with_forces["forces"] = forces
                molecule = db_file.build_molecule_from_database(
                    mol_dict_with_forces
                )
                # Only show energy at record's method/basis
                energy = energy_from_forces
                energies = (
                    [(method, basis, energy)] if energy is not None else []
                )
                frames.append(
                    {
                        "molecule": molecule,
                        "structure_id": mol_dict.get("structure_id"),
                        "energies": energies,
                        "primary_energy": energy,
                        "primary_method": (
                            method if energy is not None else None
                        ),
                        "primary_basis": basis if energy is not None else None,
                    }
                )
            else:
                molecule = db_file.build_molecule_from_database(mol_dict)
                energy = mol_dict.get("energy")
                energies = (
                    [(method, basis, energy)] if energy is not None else []
                )
                frames.append(
                    {
                        "molecule": molecule,
                        "structure_id": mol_dict.get("structure_id"),
                        "energies": energies,
                        "primary_energy": energy,
                        "primary_method": (
                            method if energy is not None else None
                        ),
                        "primary_basis": basis if energy is not None else None,
                    }
                )
        return frames

    def frames_from_structure_id(self, include_forces=False):
        """Build a single frame from --sid."""
        from chemsmart.io.database import DatabaseFile

        full_sid = self.db.get_structure_by_partial_id(self.structure_id)
        struct = self.db.get_structure(full_sid)
        if struct is None:
            raise ValueError(
                f"No structure found with ID '{self.structure_id}'."
            )

        db_file = DatabaseFile(filename=self.db_file)
        primary_method = primary_basis = None
        if include_forces:
            if self.user_primary_mb is not None:
                primary_method, primary_basis = self.user_primary_mb
            else:
                primary_method, primary_basis = (
                    self.db.pick_primary_forces_method_basis([full_sid])
                )
        return [
            self._build_frame(
                db_file, struct, primary_method, primary_basis, include_forces
            )
        ]

    def frames_from_molecule_id(self, include_forces=False):
        """Build frames from --mid (all conformers of a molecule).
        Frames are sorted by energy at the most frequently available
        (method, basis) key (or at the user-specified key when given).
        """
        from chemsmart.io.database import DatabaseFile

        full_mid = self.db.get_molecule_by_partial_id(self.molecule_id)
        struct_dicts = self.db.get_structures_for_molecule(full_mid)
        if not struct_dicts:
            raise ValueError(
                f"No structures found for molecule '{self.molecule_id}'."
            )

        primary_method = primary_basis = None
        if self.user_primary_mb is not None:
            primary_method, primary_basis = self.user_primary_mb
        elif include_forces:
            sids = [s.get("structure_id") for s in struct_dicts]
            primary_method, primary_basis = (
                self.db.pick_primary_forces_method_basis(sids)
            )
            if primary_method is None and primary_basis is None:
                logger.info(
                    f"No forces found for any structure of molecule "
                    f"{full_mid[:12]}; extxyz frames will be skipped."
                )
            else:
                mb = (
                    "/".join(x for x in (primary_method, primary_basis) if x)
                    or "unknown"
                )
                logger.info(
                    f"Exporting extxyz frames for molecule "
                    f"{full_mid[:12]} using forces/energy at {mb}; "
                    f"structures lacking forces at this level will be skipped."
                )

        db_file = DatabaseFile(filename=self.db_file)
        frames = [
            self._build_frame(
                db_file, struct, primary_method, primary_basis, include_forces
            )
            for struct in struct_dicts
        ]
        if self.user_primary_mb is not None:
            sort_primary = self.user_primary_mb
        elif include_forces and (
            primary_method is not None or primary_basis is not None
        ):
            sort_primary = (primary_method, primary_basis)
        else:
            sort_primary = None
        return frames, sort_primary

    def _build_frame(
        self,
        db_file,
        struct,
        primary_method,
        primary_basis,
        include_forces,
    ):
        """Build a single export frame dict."""
        sid = struct.get("structure_id")
        primary_energy = None
        method_used = basis_used = None
        if (
            include_forces
            and sid is not None
            and (primary_method is not None or primary_basis is not None)
        ):
            forces, primary_energy = self.db.get_forces_for_structure_at(
                sid, primary_method, primary_basis
            )
            struct["forces"] = forces
            if forces is not None or primary_energy is not None:
                method_used = primary_method
                basis_used = primary_basis
        molecule = db_file.build_molecule_from_database(struct)
        energies = collect_energies_for_structure(self.db_file, sid)
        return {
            "molecule": molecule,
            "structure_id": sid,
            "energies": energies,
            "primary_energy": primary_energy,
            "primary_method": method_used,
            "primary_basis": basis_used,
        }

    def _write_xyz_frame(self, f, frame):
        """Write a single XYZ frame with a custom comment line."""
        molecule = frame["molecule"]
        sid = frame.get("structure_id") or ""
        energies = frame.get("energies", [])

        parts = [os.path.basename(self.output)]
        if sid:
            parts.append(f"SID: {str(sid)[:12]}")
        if molecule.chemical_formula:
            parts.append(f"Empirical formula: {molecule.chemical_formula}")
        for method, basis, energy in energies:
            mb = "/".join(x for x in (method, basis) if x) or "unknown"
            parts.append(f"Energy({mb}): {energy:.10f} Eh")
        comment = "    ".join(parts)

        symbols = molecule.chemical_symbols
        positions = molecule.positions
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for s, (x, y, z) in zip(symbols, positions):
            f.write(f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f}\n")

    def _write_extxyz_frame(self, f, frame):
        """Write a single extended-XYZ frame directly to an open file handle."""
        molecule = frame["molecule"]
        sid = frame.get("structure_id") or ""
        energies = frame.get("energies", [])
        primary_energy = frame.get("primary_energy")
        primary_method = frame.get("primary_method")
        primary_basis = frame.get("primary_basis")

        forces = molecule.forces
        forces_list = None
        if forces is not None:
            try:
                forces_arr = np.asarray(forces, dtype=float)
                if forces_arr.ndim == 2 and forces_arr.shape == (
                    molecule.num_atoms,
                    3,
                ):
                    forces_list = forces_arr.tolist()
            except (TypeError, ValueError):
                pass

        # Build Properties / key=value header line.
        properties = "species:S:1:pos:R:3"
        if forces_list is not None:
            properties += ":forces:R:3"
        parts = [f"Properties={properties}"]
        if primary_energy is not None:
            parts.append(f"energy={float(primary_energy):.10f}")
            parts.append('energy_units="Hartree"')
        if forces_list is not None:
            parts.append('forces_units="Hartree/Bohr"')
        if primary_energy is not None:
            source = (
                "/".join(x for x in (primary_method, primary_basis) if x)
                or "unknown"
            )
            parts.append(f'source="{source}"')
        # Human-readable comment field.
        comment_parts = [os.path.basename(self.output)]
        if sid:
            comment_parts.append(f"SID: {str(sid)[:12]}")
        if molecule.chemical_formula:
            comment_parts.append(f"formula={molecule.chemical_formula}")
        for i, (method, basis, energy) in enumerate(energies):
            mb = "/".join(x for x in (method, basis) if x) or "unknown"
            tag = "E" if i == 0 else f"E{i}"
            comment_parts.append(f"{tag}({mb})={energy:.10f}Eh")
        comment_str = " | ".join(comment_parts).replace('"', "'")
        parts.append(f'comment="{comment_str}"')
        header = " ".join(parts)
        symbols = molecule.chemical_symbols
        positions = molecule.positions
        f.write(f"{len(symbols)}\n")
        f.write(f"{header}\n")
        if forces_list is None:
            for s, (x, y, z) in zip(symbols, positions):
                f.write(f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f}\n")
        else:
            for i, (s, (x, y, z)) in enumerate(zip(symbols, positions)):
                fx, fy, fz = forces_list[i]
                f.write(
                    f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f} "
                    f"{fx:15.10f} {fy:15.10f} {fz:15.10f}\n"
                )

    def _infer_format(self):
        """Infer output format from the file extension."""
        ext = os.path.splitext(self.output)[1].lower()
        if ext not in SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported output format '{ext}'. "
                f"Supported extensions: {', '.join(sorted(SUPPORTED_FORMATS))}"
            )
        return ext

    def _parse_csv_keys(self):
        """Parse and validate the comma-separated keys."""
        if self.keys is None:
            return None
        parsed = [k.strip() for k in self.keys.split(",") if k.strip()]
        invalid = set(parsed) - CSV_OPTIONAL_COLUMNS
        if invalid:
            raise ValueError(
                f"Unsupported CSV key(s): {', '.join(sorted(invalid))}. "
                f"Supported: {', '.join(sorted(CSV_OPTIONAL_COLUMNS))}"
            )
        return parsed

    def _validate_forces(self, molecule):
        """Check if molecule has valid forces with correct shape.
        Returns True only if forces exist and have shape (num_atoms, 3).
        """
        forces = molecule.forces
        if forces is None:
            return False
        try:
            forces_arr = np.asarray(forces, dtype=float)
            return forces_arr.ndim == 2 and forces_arr.shape == (
                molecule.num_atoms,
                3,
            )
        except (TypeError, ValueError):
            return False
