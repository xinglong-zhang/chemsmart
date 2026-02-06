"""
Database module for storing assembled calculation records in SQLite.

This module provides a SQLite-based database for storing assembled records
from quantum chemistry calculations. The database uses two tables:
- records: One row per calculation record (meta, results, provenance)
- molecules: One row per molecule in each record (structure, properties)
"""

import logging
import sqlite3
from pathlib import Path

from chemsmart.assembler.records import AssembledRecord
from chemsmart.assembler.utils import convert_numpy, from_json, to_json

logger = logging.getLogger(__name__)


def record_to_dict(record):
    """Convert AssembledRecord to dictionary.

    Args:
        record: An AssembledRecord instance or dictionary.

    Returns:
        Dictionary representation of the record.
    """
    if isinstance(record, AssembledRecord):
        return {
            "record_id": record.record_id,
            "meta": record.meta,
            "results": record.results,
            "molecules": record.molecules,
            "provenance": record.provenance,
        }
    return record


class Database:
    """SQLite database for storing assembled calculation records.

    This class provides methods for creating, populating, and managing
    a database of assembled calculation records. Uses two tables:
    - records: calculation-level data (one row per record)
    - molecules: molecule-level data (one row per molecule)

    Attributes:
        db_file: Path to the SQLite database file.
    """

    def __init__(self, db_file="database.db"):
        """Initialize database connection."""
        self.db_file = db_file
        if not self.db_file.endswith(".db"):
            self.db_file = self.db_file + ".db"

    def create(self):
        """Create the database with the required schema (records and molecules tables)."""
        conn = sqlite3.connect(self.db_file)
        try:
            # Create records table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS records (
                    record_index INTEGER PRIMARY KEY AUTOINCREMENT,
                    record_id TEXT UNIQUE NOT NULL,
                    program TEXT,
                    -- Meta fields
                    functional TEXT,
                    basis TEXT,
                    num_basis_functions INTEGER,
                    spin TEXT,
                    jobtype TEXT,
                    solvent_on INTEGER,
                    route_string TEXT,
                    solvent_model TEXT,
                    solvent_id TEXT,
                    custom_solvent TEXT,
                    temperature_in_K REAL,
                    pressure_in_atm REAL,
                    -- Results fields
                    total_energy REAL,
                    num_unpaired_electrons INTEGER,
                    homo_energy REAL,
                    lumo_energy REAL,
                    alpha_homo_energy REAL,
                    beta_homo_energy REAL,
                    alpha_lumo_energy REAL,
                    beta_lumo_energy REAL,
                    somo_energies_json TEXT,
                    lowest_somo_energy REAL,
                    highest_somo_energy REAL,
                    fmo_gap REAL,
                    alpha_fmo_gap REAL,
                    beta_fmo_gap REAL,
                    total_core_hours REAL,
                    total_elapsed_walltime REAL,
                    -- Thermochemistry results
                    rotational_temperatures_json TEXT,
                    rotational_constants_json TEXT,
                    zero_point_energy REAL,
                    thermal_vibration_correction REAL,
                    thermal_rotation_correction REAL,
                    thermal_translation_correction REAL,
                    thermal_energy_correction REAL,
                    thermal_enthalpy_correction REAL,
                    thermal_free_energy_correction REAL,
                    internal_energy REAL,
                    enthalpy REAL,
                    electronic_entropy REAL,
                    vibrational_entropy REAL,
                    rotational_entropy REAL,
                    translational_entropy REAL,
                    entropy REAL,
                    entropy_times_temperature REAL,
                    gibbs_free_energy REAL,
                    -- Provenance fields
                    source_file TEXT,
                    source_file_hash TEXT,
                    source_file_size INTEGER,
                    source_file_date TEXT,
                    program_version TEXT,
                    parser TEXT,
                    chemsmart_version TEXT,
                    assembled_at TEXT
                )
            """)
            # Create molecules table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS molecules (
                    molecule_index INTEGER PRIMARY KEY AUTOINCREMENT,
                    record_id TEXT NOT NULL,
                    index_in_record INTEGER,
                    structure_index_in_file INTEGER,
                    is_optimized_structure INTEGER,
                    charge INTEGER,
                    multiplicity INTEGER,
                    chemical_symbols_json TEXT,
                    positions_json TEXT,
                    chemical_formula TEXT,
                    number_of_atoms INTEGER,
                    mass REAL,
                    elements_json TEXT,
                    element_counts_json TEXT,
                    center_of_mass_json TEXT,
                    is_chiral INTEGER,
                    is_ring INTEGER,
                    is_monoatomic INTEGER,
                    is_diatomic INTEGER,
                    is_linear INTEGER,
                    smiles TEXT,
                    moments_of_inertia_json TEXT,
                    frozen_atoms_json TEXT,
                    energy REAL,
                    forces_json TEXT,
                    mulliken_atomic_charges_json TEXT,
                    rotational_symmetry_number INTEGER,
                    num_vibrational_modes INTEGER,
                    vibrational_frequencies_json TEXT,
                    vibrational_modes_json TEXT,
                    FOREIGN KEY (record_id) REFERENCES records(record_id)
                )
            """)
            # Create indices
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_record_id ON records(record_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_program ON records(program)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_mol_record_id ON molecules(record_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_chemical_formula ON molecules(chemical_formula)"
            )
            conn.commit()
            logger.debug(f"Created database at {self.db_file}.")
        finally:
            conn.close()

    def insert_record(
        self,
        record,
        program,
        conn=None,
    ):
        """Insert a single record into the database.

        Args:
            record: The record to insert (AssembledRecord or dict).
            program: The program used for the calculation (e.g., 'gaussian', 'orca').
            conn: Optional existing database connection.
        """
        record_dict = record_to_dict(record)
        record_dict = convert_numpy(record_dict)

        close_conn = False
        if conn is None:
            conn = sqlite3.connect(self.db_file)
            close_conn = True

        try:
            self._insert_record_row(conn, record_dict, program)
            self._insert_molecule_rows(conn, record_dict)
            if close_conn:
                conn.commit()
        finally:
            if close_conn:
                conn.close()

    def _insert_record_row(self, conn, record_dict, program):
        """Insert a row into the records table."""
        meta = record_dict.get("meta", {})
        results = record_dict.get("results", {})
        provenance = record_dict.get("provenance", {})

        conn.execute(
            """
            INSERT OR REPLACE INTO records (
                record_id, program,
                functional, basis, num_basis_functions, spin, jobtype,
                solvent_on, route_string, solvent_model, solvent_id, custom_solvent,
                temperature_in_K, pressure_in_atm,
                total_energy, num_unpaired_electrons,
                homo_energy, lumo_energy, alpha_homo_energy, beta_homo_energy,
                alpha_lumo_energy, beta_lumo_energy, somo_energies_json,
                lowest_somo_energy, highest_somo_energy,
                fmo_gap, alpha_fmo_gap, beta_fmo_gap,
                total_core_hours, total_elapsed_walltime,
                rotational_temperatures_json, rotational_constants_json,
                zero_point_energy, thermal_vibration_correction,
                thermal_rotation_correction, thermal_translation_correction,
                thermal_energy_correction, thermal_enthalpy_correction,
                thermal_free_energy_correction, internal_energy, enthalpy,
                electronic_entropy, vibrational_entropy, rotational_entropy,
                translational_entropy, entropy, entropy_times_temperature,
                gibbs_free_energy,
                source_file, source_file_hash, source_file_size, source_file_date,
                program_version, parser, chemsmart_version, assembled_at
            ) VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            )
            """,
            (
                record_dict.get("record_id"),
                program,
                # Meta
                meta.get("functional"),
                meta.get("basis"),
                meta.get("num_basis_functions"),
                meta.get("spin"),
                meta.get("jobtype"),
                1 if meta.get("solvent_on") else 0,
                meta.get("route_string"),
                meta.get("solvent_model"),
                meta.get("solvent_id"),
                meta.get("custom_solvent"),
                meta.get("temperature_in_K"),
                meta.get("pressure_in_atm"),
                # Results
                results.get("total_energy"),
                results.get("num_unpaired_electrons"),
                results.get("homo_energy"),
                results.get("lumo_energy"),
                results.get("alpha_homo_energy"),
                results.get("beta_homo_energy"),
                results.get("alpha_lumo_energy"),
                results.get("beta_lumo_energy"),
                to_json(results.get("somo_energies")),
                results.get("lowest_somo_energy"),
                results.get("highest_somo_energy"),
                results.get("fmo_gap"),
                results.get("alpha_fmo_gap"),
                results.get("beta_fmo_gap"),
                results.get("total_core_hours"),
                results.get("total_elapsed_walltime"),
                to_json(results.get("rotational_temperatures_in_K")),
                to_json(results.get("rotational_constants_in_Hz")),
                results.get("zero_point_energy"),
                results.get("thermal_vibration_correction"),
                results.get("thermal_rotation_correction"),
                results.get("thermal_translation_correction"),
                results.get("thermal_energy_correction"),
                results.get("thermal_enthalpy_correction"),
                results.get("thermal_free_energy_correction"),
                results.get("internal_energy"),
                results.get("enthalpy"),
                results.get("electronic_entropy"),
                results.get("vibrational_entropy"),
                results.get("rotational_entropy"),
                results.get("translational_entropy"),
                results.get("entropy"),
                results.get("entropy_times_temperature"),
                results.get("gibbs_free_energy"),
                # Provenance
                provenance.get("source_file"),
                provenance.get("source_file_hash"),
                provenance.get("source_file_size"),
                provenance.get("source_file_date"),
                provenance.get("program_version"),
                provenance.get("parser"),
                provenance.get("chemsmart_version"),
                provenance.get("assembled_at"),
            ),
        )

    def _insert_molecule_rows(self, conn, record_dict):
        """Insert rows into the molecules table for each molecule in the record."""
        record_id = record_dict.get("record_id")
        molecules = record_dict.get("molecules", [])

        # Delete existing molecules for this record (for REPLACE behavior)
        conn.execute("DELETE FROM molecules WHERE record_id = ?", (record_id,))

        for mol in molecules:
            conn.execute(
                """
                INSERT INTO molecules (
                    record_id, index_in_record, structure_index_in_file,
                    is_optimized_structure, charge, multiplicity,
                    chemical_symbols_json, positions_json, chemical_formula,
                    number_of_atoms, mass, elements_json, element_counts_json,
                    center_of_mass_json, is_chiral, is_ring, is_monoatomic,
                    is_diatomic, is_linear, smiles, moments_of_inertia_json,
                    frozen_atoms_json, energy, forces_json,
                    mulliken_atomic_charges_json, rotational_symmetry_number,
                    num_vibrational_modes, vibrational_frequencies_json,
                    vibrational_modes_json
                ) VALUES (
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
                )
                """,
                (
                    record_id,
                    mol.get("index"),
                    mol.get("structure_index_in_file"),
                    1 if mol.get("is_optimized_structure") else 0,
                    mol.get("charge"),
                    mol.get("multiplicity"),
                    to_json(mol.get("chemical_symbols")),
                    to_json(mol.get("positions")),
                    mol.get("chemical_formula"),
                    mol.get("number_of_atoms"),
                    mol.get("mass"),
                    to_json(mol.get("elements")),
                    to_json(mol.get("element_counts")),
                    to_json(mol.get("center_of_mass")),
                    1 if mol.get("is_chiral") else 0,
                    1 if mol.get("is_ring") else 0,
                    1 if mol.get("is_monoatomic") else 0,
                    1 if mol.get("is_diatomic") else 0,
                    1 if mol.get("is_linear") else 0,
                    mol.get("smiles"),
                    to_json(mol.get("moments_of_inertia")),
                    to_json(mol.get("frozen_atoms")),
                    mol.get("energy"),
                    to_json(mol.get("forces")),
                    to_json(mol.get("mulliken_atomic_charges")),
                    mol.get("rotational_symmetry_number"),
                    mol.get("num_vibrational_modes"),
                    to_json(mol.get("vibrational_frequencies")),
                    to_json(mol.get("vibrational_modes")),
                ),
            )

    def insert_records(self, records, program):
        """Insert multiple records into the database.

        Args:
            records: List of records to insert.
            program: The program used for the calculations.

        Returns:
            Number of records successfully inserted.
        """
        conn = sqlite3.connect(self.db_file)
        count = 0
        try:
            for record in records:
                try:
                    self.insert_record(record, program, conn)
                    count += 1
                except Exception as e:
                    logger.error(f"Failed to insert record: {e}")
            conn.commit()
        finally:
            conn.close()
        return count

    def get_connection(self):
        """Get a database connection.

        Returns:
            SQLite connection object.
        """
        return sqlite3.connect(self.db_file)

    def count_records(self):
        """Count the number of records in the database.

        Returns:
            Number of records.
        """
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute("SELECT COUNT(*) FROM records")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def exists(self):
        """Check if the database file exists.

        Returns:
            True if database file exists.
        """
        return Path(self.db_file).exists()

    def get_record(self, record_index=None, record_id=None):
        """Get a full record with molecules from the database.

        Args:
            record_index: Record index (1-based).
            record_id: Record ID string.

        Returns:
            Full record dictionary with meta, results, molecules, provenance.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if record_index is not None:
                cursor = conn.execute(
                    "SELECT * FROM records WHERE record_index = ?",
                    (record_index,),
                )
            elif record_id is not None:
                cursor = conn.execute(
                    "SELECT * FROM records WHERE record_id = ?",
                    (record_id,),
                )
            else:
                return None

            row = cursor.fetchone()
            if row is None:
                return None

            return self._row_to_full_record(conn, dict(row))
        finally:
            conn.close()

    def get_all_records(self):
        """Get all records from the database.

        Returns:
            List of full record dictionaries.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                "SELECT * FROM records ORDER BY record_index"
            )
            rows = cursor.fetchall()
            return [self._row_to_full_record(conn, dict(row)) for row in rows]
        finally:
            conn.close()

    def _row_to_full_record(self, conn, record_row):
        """Convert a record row to full record dictionary with molecules."""
        record_id = record_row["record_id"]

        # Get molecules for this record
        cursor = conn.execute(
            "SELECT * FROM molecules WHERE record_id = ? ORDER BY index_in_record",
            (record_id,),
        )
        mol_rows = cursor.fetchall()

        return {
            "record_index": record_row["record_index"],
            "record_id": record_id,
            "meta": self._extract_meta(record_row),
            "results": self._extract_results(record_row),
            "molecules": [self._mol_row_to_dict(dict(m)) for m in mol_rows],
            "provenance": self._extract_provenance(record_row),
        }

    @staticmethod
    def _extract_meta(row):
        """Extract meta fields from record row."""
        meta = {
            "functional": row.get("functional"),
            "basis": row.get("basis"),
            "num_basis_functions": row.get("num_basis_functions"),
            "spin": row.get("spin"),
            "jobtype": row.get("jobtype"),
            "solvent_on": bool(row.get("solvent_on")),
            "route_string": row.get("route_string"),
        }
        if meta["solvent_on"]:
            meta["solvent_model"] = row.get("solvent_model")
            meta["solvent_id"] = row.get("solvent_id")
            meta["custom_solvent"] = row.get("custom_solvent")
        if row.get("temperature_in_K") is not None:
            meta["temperature_in_K"] = row.get("temperature_in_K")
            meta["pressure_in_atm"] = row.get("pressure_in_atm")
        return meta

    @staticmethod
    def _extract_results(row):
        """Extract results fields from record row."""
        results = {
            "total_energy": row.get("total_energy"),
            "num_unpaired_electrons": row.get("num_unpaired_electrons"),
            "homo_energy": row.get("homo_energy"),
            "lumo_energy": row.get("lumo_energy"),
            "alpha_homo_energy": row.get("alpha_homo_energy"),
            "beta_homo_energy": row.get("beta_homo_energy"),
            "alpha_lumo_energy": row.get("alpha_lumo_energy"),
            "beta_lumo_energy": row.get("beta_lumo_energy"),
            "somo_energies": from_json(row.get("somo_energies_json")),
            "lowest_somo_energy": row.get("lowest_somo_energy"),
            "highest_somo_energy": row.get("highest_somo_energy"),
            "fmo_gap": row.get("fmo_gap"),
            "alpha_fmo_gap": row.get("alpha_fmo_gap"),
            "beta_fmo_gap": row.get("beta_fmo_gap"),
            "total_core_hours": row.get("total_core_hours"),
            "total_elapsed_walltime": row.get("total_elapsed_walltime"),
        }
        # Add thermochemistry if present
        if row.get("zero_point_energy") is not None:
            results.update(
                {
                    "rotational_temperatures_in_K": from_json(
                        row.get("rotational_temperatures_json")
                    ),
                    "rotational_constants_in_Hz": from_json(
                        row.get("rotational_constants_json")
                    ),
                    "zero_point_energy": row.get("zero_point_energy"),
                    "thermal_vibration_correction": row.get(
                        "thermal_vibration_correction"
                    ),
                    "thermal_rotation_correction": row.get(
                        "thermal_rotation_correction"
                    ),
                    "thermal_translation_correction": row.get(
                        "thermal_translation_correction"
                    ),
                    "thermal_energy_correction": row.get(
                        "thermal_energy_correction"
                    ),
                    "thermal_enthalpy_correction": row.get(
                        "thermal_enthalpy_correction"
                    ),
                    "thermal_free_energy_correction": row.get(
                        "thermal_free_energy_correction"
                    ),
                    "internal_energy": row.get("internal_energy"),
                    "enthalpy": row.get("enthalpy"),
                    "electronic_entropy": row.get("electronic_entropy"),
                    "vibrational_entropy": row.get("vibrational_entropy"),
                    "rotational_entropy": row.get("rotational_entropy"),
                    "translational_entropy": row.get("translational_entropy"),
                    "entropy": row.get("entropy"),
                    "entropy_times_temperature": row.get(
                        "entropy_times_temperature"
                    ),
                    "gibbs_free_energy": row.get("gibbs_free_energy"),
                }
            )
        return results

    @staticmethod
    def _extract_provenance(row):
        """Extract provenance fields from record row."""
        return {
            "source_file": row.get("source_file"),
            "source_file_hash": row.get("source_file_hash"),
            "source_file_size": row.get("source_file_size"),
            "source_file_date": row.get("source_file_date"),
            "program": row.get("program"),
            "program_version": row.get("program_version"),
            "parser": row.get("parser"),
            "chemsmart_version": row.get("chemsmart_version"),
            "assembled_at": row.get("assembled_at"),
        }

    @staticmethod
    def _mol_row_to_dict(row):
        """Convert molecule row to dictionary."""
        mol = {
            "index": row.get("index_in_record"),
            "structure_index_in_file": row.get("structure_index_in_file"),
            "is_optimized_structure": bool(row.get("is_optimized_structure")),
            "charge": row.get("charge"),
            "multiplicity": row.get("multiplicity"),
            "chemical_symbols": from_json(row.get("chemical_symbols_json")),
            "positions": from_json(row.get("positions_json")),
            "chemical_formula": row.get("chemical_formula"),
            "number_of_atoms": row.get("number_of_atoms"),
            "mass": row.get("mass"),
            "elements": from_json(row.get("elements_json")),
            "element_counts": from_json(row.get("element_counts_json")),
            "center_of_mass": from_json(row.get("center_of_mass_json")),
            "is_chiral": bool(row.get("is_chiral")),
            "is_ring": bool(row.get("is_ring")),
            "is_monoatomic": bool(row.get("is_monoatomic")),
            "is_diatomic": bool(row.get("is_diatomic")),
            "is_linear": bool(row.get("is_linear")),
            "smiles": row.get("smiles"),
            "moments_of_inertia": from_json(
                row.get("moments_of_inertia_json")
            ),
            "frozen_atoms": from_json(row.get("frozen_atoms_json")),
            "energy": row.get("energy"),
            "forces": from_json(row.get("forces_json")),
        }
        # Optional fields
        mulliken = from_json(row.get("mulliken_atomic_charges_json"))
        if mulliken is not None:
            mol["mulliken_atomic_charges"] = mulliken
        rot_sym = row.get("rotational_symmetry_number")
        if rot_sym is not None:
            mol["rotational_symmetry_number"] = rot_sym
        if row.get("num_vibrational_modes") is not None:
            mol["num_vibrational_modes"] = row.get("num_vibrational_modes")
            mol["vibrational_frequencies"] = from_json(
                row.get("vibrational_frequencies_json")
            )
            mol["vibrational_modes"] = from_json(
                row.get("vibrational_modes_json")
            )
        return mol

    def get_record_summary(
        self,
        record_index=None,
        record_id=None,
    ):
        """Get summary of a record (without molecules).

        Args:
            record_index: Record index (1-based).
            record_id: Record ID string.

        Returns:
            Record summary dictionary.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if record_index is not None:
                cursor = conn.execute(
                    """SELECT record_index, record_id, program, functional, basis, 
                       jobtype, total_energy, homo_energy, lumo_energy, fmo_gap
                       FROM records WHERE record_index = ?""",
                    (record_index,),
                )
            elif record_id is not None:
                cursor = conn.execute(
                    """SELECT record_index, record_id, program, functional, basis,
                       jobtype, total_energy, homo_energy, lumo_energy, fmo_gap
                       FROM records WHERE record_id = ?""",
                    (record_id,),
                )
            else:
                return None

            row = cursor.fetchone()
            if row is None:
                return None

            # Get last molecule info
            mol_cursor = conn.execute(
                """SELECT chemical_formula, charge, multiplicity, smiles 
                   FROM molecules WHERE record_id = ? 
                   ORDER BY index_in_record DESC LIMIT 1""",
                (row["record_id"],),
            )
            mol_row = mol_cursor.fetchone()

            summary = dict(row)
            if mol_row:
                summary.update(dict(mol_row))
            return summary
        finally:
            conn.close()

    def get_all_record_summaries(self):
        """Get summaries for all records.

        Returns:
            List of record summary dictionaries.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """SELECT record_index, record_id, program, functional, basis,
                   jobtype, total_energy, homo_energy, lumo_energy, fmo_gap
                   FROM records ORDER BY record_index"""
            )
            rows = cursor.fetchall()

            summaries = []
            for row in rows:
                summary = dict(row)
                # Get last molecule info
                mol_cursor = conn.execute(
                    """SELECT chemical_formula, charge, multiplicity, smiles
                       FROM molecules WHERE record_id = ?
                       ORDER BY index_in_record DESC LIMIT 1""",
                    (row["record_id"],),
                )
                mol_row = mol_cursor.fetchone()
                if mol_row:
                    summary.update(dict(mol_row))
                summaries.append(summary)
            return summaries
        finally:
            conn.close()

    def get_molecules_for_record(self, record_id):
        """Get all molecules for a specific record.

        Args:
            record_id: Record ID string.

        Returns:
            List of molecule dictionaries.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                "SELECT * FROM molecules WHERE record_id = ? ORDER BY index_in_record",
                (record_id,),
            )
            rows = cursor.fetchall()
            return [self._mol_row_to_dict(dict(row)) for row in rows]
        finally:
            conn.close()

    def query(self, where_clause):
        """Query records with a SQL WHERE clause.

        Args:
            where_clause: SQL WHERE clause (without the WHERE keyword).

        Returns:
            List of full record dictionaries.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            # Join records with last molecule to allow querying by chemical_formula
            sql = f"""
                SELECT DISTINCT r.* FROM records r
                LEFT JOIN molecules m ON r.record_id = m.record_id
                WHERE {where_clause}
                ORDER BY r.record_index
            """
            cursor = conn.execute(sql)
            rows = cursor.fetchall()
            return [self._row_to_full_record(conn, dict(row)) for row in rows]
        finally:
            conn.close()

    def query_summaries(self, where_clause):
        """Query record summaries with a SQL WHERE clause.

        Args:
            where_clause: SQL WHERE clause (without the WHERE keyword).

        Returns:
            List of record summary dictionaries.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            sql = f"""
                SELECT DISTINCT r.record_index, r.record_id, r.program,
                       r.functional, r.basis, r.jobtype, r.total_energy,
                       m.chemical_formula, m.charge, m.multiplicity
                FROM records r
                LEFT JOIN molecules m ON r.record_id = m.record_id
                WHERE {where_clause}
                ORDER BY r.record_index
            """
            cursor = conn.execute(sql)
            rows = cursor.fetchall()
            return [dict(row) for row in rows]
        finally:
            conn.close()
