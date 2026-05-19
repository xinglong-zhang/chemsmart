"""
Database export module for exporting records from a chemsmart database.

Output format is inferred from the file extension of the output path:
* .json – Full structured records.
* .csv  – Tabular scalar properties.
* .xyz  – Cartesian coordinates of selected structure(s).
"""

import csv
import json
import logging
import os
import sqlite3

from chemsmart.database.database import Database
from chemsmart.database.utils import convert_numpy

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

SUPPORTED_FORMATS = {".json", ".csv", ".xyz"}


class DatabaseExporter:
    """Export a chemsmart database to JSON, CSV, or XYZ."""

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
        self.format = self._infer_format()
        self.parsed_keys = self._parse_csv_keys()

    def export(self):
        """Run export, dispatching to the appropriate format handler."""
        handler = {
            ".json": self.to_json,
            ".csv": self.to_csv,
            ".xyz": self.to_xyz,
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
        selectors = (
            self.record_index,
            self.record_id,
            self.structure_id,
            self.molecule_id,
        )
        if sum(s is not None for s in selectors) != 1:
            raise ValueError(
                "XYZ export requires exactly one of --ri/--record-index, "
                "--rid/--record-id, --sid/--structure-id, or "
                "--mid/--molecule-id."
            )

        if self.molecule_id is not None:
            frames = self.frames_from_molecule_id()
        elif self.structure_id is not None:
            frames = self._frames_from_structure_id()
        else:
            frames = self._frames_from_record()

        if not frames:
            raise ValueError("No structures found for the given selection.")

        logger.info(
            f"Writing {len(frames)} frame(s) to "
            f"{os.path.basename(self.output)}"
        )
        with open(self.output, "w") as f:
            for frame in frames:
                self._write_xyz_frame(f, frame)

    def _frames_from_record(self):
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
        frames = []
        for mol_dict in selected:
            molecule = db_file.build_molecule_from_database(mol_dict)
            energy = mol_dict.get("energy")
            energies = [(method, basis, energy)] if energy is not None else []
            frames.append(
                {
                    "molecule": molecule,
                    "structure_id": mol_dict.get("structure_id"),
                    "energies": energies,
                }
            )
        return frames

    def _frames_from_structure_id(self):
        """Build a single frame from --sid."""
        from chemsmart.io.database import DatabaseFile

        full_sid = self.db.get_structure_by_partial_id(self.structure_id)
        struct = self.db.get_structure(full_sid)
        if struct is None:
            raise ValueError(
                f"No structure found with ID '{self.structure_id}'."
            )

        db_file = DatabaseFile(filename=self.db_file)
        molecule = db_file.build_molecule_from_database(struct)
        energies = self._collect_energies_for_structure(full_sid)
        return [
            {
                "molecule": molecule,
                "structure_id": full_sid,
                "energies": energies,
            }
        ]

    def frames_from_molecule_id(self):
        """Build frames from --mid (all conformers of a molecule).
        Frames are sorted by energy at the most frequently available
        (method, basis) key.
        """
        from chemsmart.io.database import DatabaseFile

        full_mid = self.db.get_molecule_by_partial_id(self.molecule_id)
        struct_dicts = self.db.get_structures_for_molecule(full_mid)
        if not struct_dicts:
            raise ValueError(
                f"No structures found for molecule '{self.molecule_id}'."
            )

        db_file = DatabaseFile(filename=self.db_file)
        frames = []
        for struct in struct_dicts:
            sid = struct.get("structure_id")
            molecule = db_file.build_molecule_from_database(struct)
            energies = self._collect_energies_for_structure(sid)
            frames.append(
                {
                    "molecule": molecule,
                    "structure_id": sid,
                    "energies": energies,
                }
            )
        return self.sort_frames_by_energy(frames)

    @staticmethod
    def sort_frames_by_energy(frames):
        """Sort frames ascending by energy at the most-covered (method,
        basis); frames missing that key go to the end (sorted by their
        own lowest available energy)."""
        # 1) Count (method, basis) coverage across all frames
        counts = {}
        for frame in frames:
            for method, basis, _ in frame.get("energies", []):
                key = (method, basis)
                counts[key] = counts.get(key, 0) + 1
        if not counts:
            return frames  # nothing to sort by

        # 2) Pick the most frequent (method, basis) as primary key
        primary = max(counts.items(), key=lambda kv: kv[1])[0]
        logger.info(
            f"Sorting {len(frames)} frame(s) by energy at "
            f"{primary[0]}/{primary[1]}; frames missing it go to the end."
        )

        def sort_key(frame):
            energies = frame.get("energies", [])
            primary_e = next(
                (e for m, b, e in energies if (m, b) == primary), None
            )
            sid = str(frame.get("structure_id") or "")
            if primary_e is not None:
                # Bucket 0: has primary key -> sort by its energy
                return (0, primary_e, sid)
            # Bucket 1: missing primary key -> sort by lowest available
            fallback = min((e for _, _, e in energies), default=float("inf"))
            return (1, fallback, sid)

        sorted_frames = sorted(frames, key=sort_key)

        # 3) Within each frame, move the primary (method, basis) entry
        #    to the front so the comment line shows it first.
        for frame in sorted_frames:
            energies = frame.get("energies", [])
            head = [(m, b, e) for (m, b, e) in energies if (m, b) == primary]
            tail = [(m, b, e) for (m, b, e) in energies if (m, b) != primary]
            frame["energies"] = head + tail
        return sorted_frames

    def _collect_energies_for_structure(self, structure_id):
        """Return list of (method, basis, energy) for every record
        that references the given structure_id."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT r.method AS method, r.basis AS basis,
                       rs.energy AS energy
                FROM record_structures rs
                JOIN records r ON rs.record_id = r.record_id
                WHERE rs.structure_id = ?
                ORDER BY r.record_index
                """,
                (structure_id,),
            )
            rows = cursor.fetchall()
        finally:
            conn.close()
        return [
            (r["method"], r["basis"], r["energy"])
            for r in rows
            if r["energy"] is not None
        ]

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
            parts.append(f"Energy({mb}): {energy:.8f} Eh")
        comment = "    ".join(parts)

        symbols = molecule.chemical_symbols
        positions = molecule.positions
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for s, (x, y, z) in zip(symbols, positions):
            f.write(f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f}\n")

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
