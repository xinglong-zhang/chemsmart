"""
Database module for storing assembled calculation records in SQLite.

This module provides a SQLite-based database for storing assembled records
from quantum chemistry calculations. The database uses four tables:
- records: One row per calculation record (meta, results, provenance)
- molecules: One row per chemical species (molecule-level identity)
- structures: One row per unique 3D geometry instance (conformer)
- record_structures: Junction table linking records to structures
"""

import logging
import sqlite3
from pathlib import Path

from chemsmart.database.records import AssembledRecord
from chemsmart.database.utils import (
    convert_numpy,
    from_json,
    to_json,
)

logger = logging.getLogger(__name__)


def record_to_dict(record):
    """Convert AssembledRecord to dictionary."""
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
    a database of assembled calculation records. Uses four tables:
    - records: calculation-level data (one row per record)
    - molecules: chemical species data (one row per unique species)
    - structures: 3D geometry data (one row per unique conformer)
    - record_structures: junction linking records to structures

    Attributes:
        db_file: Path to the SQLite database file.
    """

    def __init__(self, db_file="database.db"):
        """Initialize database connection."""
        self.db_file = db_file
        if not self.db_file.endswith(".db"):
            self.db_file = self.db_file + ".db"

    def create(self):
        """Create the database with the required schema.
        Tables:
        - records: One row per calculation record (meta, results, provenance)
        - molecules: One row per chemical species (molecule-level identity)
        - structures: One row per unique 3D geometry instance (conformer)
        - record_structures: Junction table linking records to structures
        """
        conn = sqlite3.connect(self.db_file)
        try:
            # Create records table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS records (
                    record_index INTEGER PRIMARY KEY AUTOINCREMENT,
                    record_id TEXT UNIQUE NOT NULL,
                    program TEXT,
                    -- Meta fields
                    method TEXT,
                    basis TEXT,
                    num_basis_functions INTEGER,
                    spin TEXT,
                    jobtype TEXT,
                    solvent_on INTEGER,
                    route_string TEXT,
                    solvent_model TEXT,
                    solvent_id TEXT,
                    custom_solvent TEXT,
                    custom_basis_json TEXT,
                    trajectory_id TEXT,
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
                    assembled_at TEXT,
                    normal_termination INTEGER
                )
            """)

            # Create molecules table (chemical species)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS molecules (
                    molecule_id TEXT PRIMARY KEY NOT NULL,
                    molecule_label TEXT,
                    empirical_formula TEXT,
                    chemical_formula TEXT,
                    smiles TEXT,
                    inchi TEXT,
                    chiral_centers_json TEXT,
                    number_of_atoms INTEGER,
                    mass REAL,
                    elements_json TEXT,
                    element_counts_json TEXT,
                    is_chiral INTEGER,
                    is_ring INTEGER,
                    is_aromatic INTEGER,
                    is_monoatomic INTEGER,
                    is_diatomic INTEGER,
                    is_linear INTEGER,
                    is_multicomponent INTEGER,
                    num_components INTEGER
                )
            """)

            # Create structures table (unique 3D geometry instances)
            conn.execute("""
                CREATE TABLE IF NOT EXISTS structures (
                    structure_id TEXT PRIMARY KEY,
                    molecule_id TEXT NOT NULL,
                    charge INTEGER,
                    multiplicity INTEGER,
                    structure_label TEXT,
                    chemical_symbols_json TEXT,
                    positions_json TEXT,
                    center_of_mass_json TEXT,
                    moments_of_inertia_json TEXT,
                    FOREIGN KEY (molecule_id) REFERENCES molecules(molecule_id)
                )
            """)

            # Create record_structures junction table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS record_structures (
                    record_id TEXT NOT NULL,
                    structure_id TEXT NOT NULL,
                    index_in_record INTEGER NOT NULL,
                    structure_index_in_file INTEGER,
                    is_optimized_structure INTEGER,
                    energy REAL,
                    forces_json TEXT,
                    frozen_atoms_json TEXT,
                    mulliken_atomic_charges_json TEXT,
                    mulliken_spin_densities_json TEXT,
                    rotational_symmetry_number INTEGER,
                    rotational_constants_json TEXT,
                    point_group TEXT,
                    dipole_moment_json TEXT,
                    dipole_moment_magnitude REAL,
                    num_vibrational_modes INTEGER,
                    vibrational_frequencies_json TEXT,
                    vibrational_modes_json TEXT,
                    PRIMARY KEY (record_id, index_in_record),
                    FOREIGN KEY (record_id) REFERENCES records(record_id),
                    FOREIGN KEY (structure_id) REFERENCES structures(structure_id)
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
                "CREATE INDEX IF NOT EXISTS idx_mol_empirical_formula ON molecules(empirical_formula)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_mol_chemical_formula ON molecules(chemical_formula)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_struct_molecule_id ON structures(molecule_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_rs_record_id ON record_structures(record_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_rs_structure_id ON record_structures(structure_id)"
            )
            conn.execute(
                "CREATE INDEX IF NOT EXISTS idx_records_trajectory_id ON records(trajectory_id)"
            )
            conn.commit()
            logger.debug(f"Created database at {self.db_file}.")
        finally:
            conn.close()

    def insert_record(
        self,
        record,
        conn=None,
    ):
        """Insert a single record into the database."""
        record_dict = record_to_dict(record)
        record_dict = convert_numpy(record_dict)

        program = record_dict.get("provenance", {}).get("program", "unknown")
        close_conn = False
        if conn is None:
            conn = sqlite3.connect(self.db_file)
            close_conn = True

        try:
            self._insert_record_row(conn, record_dict, program)
            self._insert_related_rows(conn, record_dict)
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
                method, basis, num_basis_functions, spin, jobtype,
                solvent_on, route_string, solvent_model, solvent_id, custom_solvent,
                custom_basis_json,
                trajectory_id,
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
                program_version, parser, chemsmart_version, assembled_at,
                normal_termination
            ) VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            )
            """,
            (
                record_dict.get("record_id"),
                program,
                # Meta
                meta.get("method"),
                meta.get("basis"),
                meta.get("num_basis_functions"),
                meta.get("spin"),
                meta.get("jobtype"),
                1 if meta.get("solvent_on") else 0,
                meta.get("route_string"),
                meta.get("solvent_model"),
                meta.get("solvent_id"),
                to_json(meta.get("custom_solvent")),
                to_json(meta.get("custom_basis")),
                meta.get("trajectory_id"),
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
                1 if provenance.get("normal_termination") else 0,
            ),
        )

    def _insert_related_rows(self, conn, record_dict):
        """Insert rows into molecules, structures, and record_structures tables.
        For each Molecule entry in the record, splits its data into:
        1. INSERT OR IGNORE into molecules (species dedup by molecule_id)
        2. INSERT OR IGNORE into structures (geometry dedup by structure_id)
        3. INSERT into record_structures (per-calculation data)
        """
        record_id = record_dict.get("record_id")
        molecules = record_dict.get("molecules", [])

        # Delete existing record_structures for this record (for REPLACE behavior)
        conn.execute(
            "DELETE FROM record_structures WHERE record_id = ?", (record_id,)
        )

        for mol in molecules:
            molecule_id = mol.get("molecule_id")
            structure_id = mol.get("structure_id")

            # Step 1: Insert species (dedup by molecule_id)
            conn.execute(
                """
                INSERT OR IGNORE INTO molecules (
                    molecule_id, molecule_label, empirical_formula,
                    chemical_formula, smiles, inchi, chiral_centers_json,
                    number_of_atoms, mass,
                    elements_json, element_counts_json, is_chiral, is_ring,
                    is_aromatic, is_monoatomic, is_diatomic, is_linear,
                    is_multicomponent, num_components
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    molecule_id,
                    mol.get("molecule_label"),
                    mol.get("empirical_formula"),
                    mol.get("chemical_formula"),
                    mol.get("smiles"),
                    mol.get("inchi"),
                    to_json(mol.get("chiral_centers")),
                    mol.get("number_of_atoms"),
                    mol.get("mass"),
                    to_json(mol.get("elements")),
                    to_json(mol.get("element_counts")),
                    1 if mol.get("is_chiral") else 0,
                    1 if mol.get("is_ring") else 0,
                    1 if mol.get("is_aromatic") else 0,
                    1 if mol.get("is_monoatomic") else 0,
                    1 if mol.get("is_diatomic") else 0,
                    1 if mol.get("is_linear") else 0,
                    1 if mol.get("is_multicomponent") else 0,
                    mol.get("num_components"),
                ),
            )

            # Step 2: Insert structure (dedup by structure_id)
            conn.execute(
                """
                INSERT OR IGNORE INTO structures (
                    structure_id, molecule_id, charge, multiplicity,
                    structure_label, chemical_symbols_json, positions_json,
                    center_of_mass_json, moments_of_inertia_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    structure_id,
                    molecule_id,
                    mol.get("charge"),
                    mol.get("multiplicity"),
                    mol.get("structure_label"),
                    to_json(mol.get("chemical_symbols")),
                    to_json(mol.get("positions")),
                    to_json(mol.get("center_of_mass")),
                    to_json(mol.get("moments_of_inertia")),
                ),
            )

            # Step 3: Insert record-structure link (per-calculation data)
            conn.execute(
                """
                INSERT INTO record_structures (
                    record_id, structure_id, index_in_record,
                    structure_index_in_file, is_optimized_structure,
                    energy, forces_json, frozen_atoms_json,
                    mulliken_atomic_charges_json, mulliken_spin_densities_json,
                    rotational_symmetry_number, rotational_constants_json,
                    point_group, dipole_moment_json, dipole_moment_magnitude,
                    num_vibrational_modes, vibrational_frequencies_json,
                    vibrational_modes_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    record_id,
                    structure_id,
                    mol.get("index"),
                    mol.get("structure_index_in_file"),
                    1 if mol.get("is_optimized_structure") else 0,
                    mol.get("energy"),
                    to_json(mol.get("forces")),
                    to_json(mol.get("frozen_atoms")),
                    to_json(mol.get("mulliken_atomic_charges")),
                    to_json(mol.get("mulliken_spin_densities")),
                    mol.get("rotational_symmetry_number"),
                    to_json(mol.get("rotational_constants")),
                    mol.get("point_group"),
                    to_json(mol.get("dipole_moment")),
                    mol.get("dipole_moment_magnitude"),
                    mol.get("num_vibrational_modes"),
                    to_json(mol.get("vibrational_frequencies")),
                    to_json(mol.get("vibrational_modes")),
                ),
            )

    def insert_records(self, records):
        """Insert multiple records into the database."""
        conn = sqlite3.connect(self.db_file)
        count = 0
        try:
            for record in records:
                try:
                    self.insert_record(record, conn)
                    count += 1
                except Exception as e:
                    logger.error(f"Failed to insert record: {e}")
            conn.commit()
        finally:
            conn.close()
        return count

    def get_connection(self):
        """Get a database connection."""
        return sqlite3.connect(self.db_file)

    def count_records(self):
        """Count the number of records in the database."""
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute("SELECT COUNT(*) FROM records")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def exists(self):
        """Check if the database file exists."""
        return Path(self.db_file).exists()

    def get_record(self, record_index=None, record_id=None):
        """Get a full record with structures from the database."""
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
        """Get all records from the database."""
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
        """Convert a record row to full record dictionary with structures.
        Joins record_structures, structures, and molecules tables to
        reconstruct the full structure data for each entry in the record.
        """
        record_id = record_row["record_id"]

        # Get structures for this record via three-table join
        cursor = conn.execute(
            """
            SELECT rs.*, s.charge, s.multiplicity,
                   s.structure_label, s.chemical_symbols_json, s.positions_json,
                   s.center_of_mass_json, s.moments_of_inertia_json,
                   m.molecule_id, m.molecule_label, m.empirical_formula,
                   m.chemical_formula, m.smiles, m.inchi, m.chiral_centers_json,
                   m.number_of_atoms, m.mass,
                   m.elements_json, m.element_counts_json, m.is_chiral,
                   m.is_ring, m.is_aromatic, m.is_monoatomic, m.is_diatomic,
                   m.is_linear, m.is_multicomponent, m.num_components
            FROM record_structures rs
            JOIN structures s ON rs.structure_id = s.structure_id
            JOIN molecules m ON s.molecule_id = m.molecule_id
            WHERE rs.record_id = ?
            ORDER BY rs.index_in_record
            """,
            (record_id,),
        )
        struct_rows = cursor.fetchall()

        return {
            "record_index": record_row["record_index"],
            "record_id": record_id,
            "meta": self._extract_meta(record_row),
            "results": self._extract_results(record_row),
            "molecules": [
                self._record_structure_row_to_dict(dict(r))
                for r in struct_rows
            ],
            "provenance": self._extract_provenance(record_row),
        }

    @staticmethod
    def _extract_meta(row):
        """Extract meta fields from record row."""
        meta = {
            "method": row.get("method"),
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
            meta["custom_solvent"] = from_json(row.get("custom_solvent"))
        if row.get("custom_basis_json") is not None:
            meta["custom_basis"] = from_json(row.get("custom_basis_json"))
        if row.get("trajectory_id") is not None:
            meta["trajectory_id"] = row.get("trajectory_id")
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
        provenance = {
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
        if row.get("normal_termination") is not None:
            provenance["normal_termination"] = bool(
                row.get("normal_termination")
            )
        return provenance

    @staticmethod
    def _record_structure_row_to_dict(row):
        """Convert a joined record_structures+structures+molecules row to dict.
        Produces a flat dict representing one structure entry within a record,
        combining per-calculation data (from record_structures), geometry data
        (from structures), and species identity (from molecules).
        """
        mol = {
            "index": row.get("index_in_record"),
            "molecule_id": row.get("molecule_id"),
            "molecule_label": row.get("molecule_label"),
            "structure_id": row.get("structure_id"),
            "structure_index_in_file": row.get("structure_index_in_file"),
            "is_optimized_structure": bool(row.get("is_optimized_structure")),
            "charge": row.get("charge"),
            "multiplicity": row.get("multiplicity"),
            "structure_label": row.get("structure_label"),
            "chemical_symbols": from_json(row.get("chemical_symbols_json")),
            "positions": from_json(row.get("positions_json")),
            "empirical_formula": row.get("empirical_formula"),
            "chemical_formula": row.get("chemical_formula"),
            "number_of_atoms": row.get("number_of_atoms"),
            "mass": row.get("mass"),
            "elements": from_json(row.get("elements_json")),
            "element_counts": from_json(row.get("element_counts_json")),
            "center_of_mass": from_json(row.get("center_of_mass_json")),
            "is_chiral": bool(row.get("is_chiral")),
            "is_ring": bool(row.get("is_ring")),
            "is_aromatic": bool(row.get("is_aromatic")),
            "is_monoatomic": bool(row.get("is_monoatomic")),
            "is_diatomic": bool(row.get("is_diatomic")),
            "is_linear": bool(row.get("is_linear")),
            "is_multicomponent": bool(row.get("is_multicomponent")),
            "num_components": row.get("num_components"),
            "smiles": row.get("smiles"),
            "inchi": row.get("inchi"),
            "chiral_centers": from_json(row.get("chiral_centers_json")),
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
        spin_densities = from_json(row.get("mulliken_spin_densities_json"))
        if spin_densities is not None:
            mol["mulliken_spin_densities"] = spin_densities
        rot_sym = row.get("rotational_symmetry_number")
        if rot_sym is not None:
            mol["rotational_symmetry_number"] = rot_sym
        rot_consts = from_json(row.get("rotational_constants_json"))
        if rot_consts is not None:
            mol["rotational_constants"] = rot_consts
        pg = row.get("point_group")
        if pg is not None:
            mol["point_group"] = pg
        dipole = from_json(row.get("dipole_moment_json"))
        if dipole is not None:
            mol["dipole_moment"] = dipole
        dipole_mag = row.get("dipole_moment_magnitude")
        if dipole_mag is not None:
            mol["dipole_moment_magnitude"] = dipole_mag
        if row.get("num_vibrational_modes") is not None:
            mol["num_vibrational_modes"] = row.get("num_vibrational_modes")
            mol["vibrational_frequencies"] = from_json(
                row.get("vibrational_frequencies_json")
            )
            mol["vibrational_modes"] = from_json(
                row.get("vibrational_modes_json")
            )
        return mol

    def get_record_by_partial_id(self, partial_id):
        """Resolve a partial record ID to the full record ID.
        Matches record IDs that start with the given prefix. Requires
        exactly one match; raises ValueError otherwise.
        """
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute(
                "SELECT record_id FROM records WHERE record_id LIKE ?",
                (partial_id + "%",),
            )
            matches = [row[0] for row in cursor.fetchall()]
        finally:
            conn.close()

        if len(matches) == 0:
            raise ValueError(
                f"No record found matching ID prefix '{partial_id}'."
            )
        if len(matches) > 1:
            preview = ", ".join(m[:12] for m in matches[:10])
            raise ValueError(
                f"Ambiguous ID prefix '{partial_id}' matches {len(matches)} records: {preview}"
            )
        return matches[0]

    def get_record_summary(
        self,
        record_index=None,
        record_id=None,
    ):
        """Get summary of a record (without structures)."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            if record_index is not None:
                cursor = conn.execute(
                    """SELECT record_index, record_id, program, method, basis, 
                       jobtype, total_energy, homo_energy, lumo_energy, fmo_gap,
                       source_file
                       FROM records WHERE record_index = ?""",
                    (record_index,),
                )
            elif record_id is not None:
                cursor = conn.execute(
                    """SELECT record_index, record_id, program, method, basis,
                       jobtype, total_energy, homo_energy, lumo_energy, fmo_gap,
                       source_file
                       FROM records WHERE record_id = ?""",
                    (record_id,),
                )
            else:
                return None

            row = cursor.fetchone()
            if row is None:
                return None

            # Get last structure's molecule info via three-table join
            struct_cursor = conn.execute(
                """SELECT m.chemical_formula, s.charge, s.multiplicity, m.smiles
                   FROM record_structures rs
                   JOIN structures s ON rs.structure_id = s.structure_id
                   JOIN molecules m ON s.molecule_id = m.molecule_id
                   WHERE rs.record_id = ?
                   ORDER BY rs.index_in_record DESC LIMIT 1""",
                (row["record_id"],),
            )
            struct_row = struct_cursor.fetchone()

            summary = dict(row)
            if struct_row:
                summary.update(dict(struct_row))
            return summary
        finally:
            conn.close()

    def get_all_record_summaries(self):
        """Get summaries for all records."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """SELECT record_index, record_id, program, method, basis,
                   jobtype, total_energy, homo_energy, lumo_energy, fmo_gap,
                   source_file
                   FROM records ORDER BY record_index"""
            )
            rows = cursor.fetchall()

            summaries = []
            for row in rows:
                summary = dict(row)
                # Get last structure's molecule info via three-table join
                struct_cursor = conn.execute(
                    """SELECT m.chemical_formula, s.charge, s.multiplicity, m.smiles
                       FROM record_structures rs
                       JOIN structures s ON rs.structure_id = s.structure_id
                       JOIN molecules m ON s.molecule_id = m.molecule_id
                       WHERE rs.record_id = ?
                       ORDER BY rs.index_in_record DESC LIMIT 1""",
                    (row["record_id"],),
                )
                struct_row = struct_cursor.fetchone()
                if struct_row:
                    summary.update(dict(struct_row))
                summaries.append(summary)
            return summaries
        finally:
            conn.close()

    def get_structures_for_record(self, record_id):
        """Get all structure entries for a specific record.
        Returns the full joined record_structures + structures + molecules
        data for each structure entry in the record.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT rs.*, s.charge, s.multiplicity,
                       s.structure_label, s.chemical_symbols_json, s.positions_json,
                       s.center_of_mass_json, s.moments_of_inertia_json,
                       m.molecule_id, m.molecule_label, m.empirical_formula,
                       m.chemical_formula, m.smiles, m.inchi, m.chiral_centers_json,
                       m.number_of_atoms, m.mass,
                       m.elements_json, m.element_counts_json, m.is_chiral,
                       m.is_ring, m.is_aromatic, m.is_monoatomic, m.is_diatomic,
                       m.is_linear, m.is_multicomponent, m.num_components
                FROM record_structures rs
                JOIN structures s ON rs.structure_id = s.structure_id
                JOIN molecules m ON s.molecule_id = m.molecule_id
                WHERE rs.record_id = ?
                ORDER BY rs.index_in_record
                """,
                (record_id,),
            )
            rows = cursor.fetchall()
            return [
                self._record_structure_row_to_dict(dict(row)) for row in rows
            ]
        finally:
            conn.close()

    def get_molecule(self, molecule_id):
        """Get a molecule (chemical species) by molecule_id."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                "SELECT * FROM molecules WHERE molecule_id = ?",
                (molecule_id,),
            )
            row = cursor.fetchone()
            if row is None:
                return None
            return self._molecule_row_to_dict(dict(row))
        finally:
            conn.close()

    def get_all_molecules(self):
        """Get all unique molecules (chemical species) from the database."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                "SELECT * FROM molecules ORDER BY chemical_formula"
            )
            rows = cursor.fetchall()
            return [self._molecule_row_to_dict(dict(row)) for row in rows]
        finally:
            conn.close()

    def get_molecule_by_partial_id(self, partial_id):
        """Resolve a partial molecule ID to the full molecule_id."""
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute(
                "SELECT molecule_id FROM molecules WHERE molecule_id LIKE ?",
                (partial_id + "%",),
            )
            matches = [row[0] for row in cursor.fetchall()]
        finally:
            conn.close()

        if len(matches) == 0:
            raise ValueError(
                f"No molecule found matching ID prefix '{partial_id}'."
            )
        if len(matches) > 1:
            preview = ", ".join(m[:12] for m in matches[:10])
            raise ValueError(
                f"Ambiguous ID prefix '{partial_id}' matches "
                f"{len(matches)} molecules: {preview}"
            )
        return matches[0]

    def count_molecules(self):
        """Count unique molecules in the database."""
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute("SELECT COUNT(*) FROM molecules")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def get_structure(self, structure_id):
        """Get a structure (conformer) by structure_id.
        Returns structure fields joined with parent molecule info for
        convenient display.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT s.*, m.chemical_formula, m.empirical_formula,
                       m.number_of_atoms, m.mass, m.smiles
                FROM structures s
                JOIN molecules m ON s.molecule_id = m.molecule_id
                WHERE s.structure_id = ?
                """,
                (structure_id,),
            )
            row = cursor.fetchone()
            if row is None:
                return None
            return self._structure_row_to_dict(dict(row))
        finally:
            conn.close()

    def get_all_structures(self):
        """Get all structures (conformers) from the database."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute("""
                SELECT s.*, m.chemical_formula, m.empirical_formula,
                       m.number_of_atoms, m.mass, m.smiles
                FROM structures s
                JOIN molecules m ON s.molecule_id = m.molecule_id
                ORDER BY s.molecule_id, s.structure_id
                """)
            rows = cursor.fetchall()
            return [self._structure_row_to_dict(dict(row)) for row in rows]
        finally:
            conn.close()

    def get_structure_by_partial_id(self, partial_id):
        """Resolve a partial structure ID to the full structure_id.
        Matches structure IDs that start with the given prefix. Requires
        exactly one match; raises ValueError otherwise.
        """
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute(
                "SELECT structure_id FROM structures WHERE structure_id LIKE ?",
                (partial_id + "%",),
            )
            matches = [row[0] for row in cursor.fetchall()]
        finally:
            conn.close()

        if len(matches) == 0:
            raise ValueError(
                f"No structure found matching ID prefix '{partial_id}'."
            )
        if len(matches) > 1:
            preview = ", ".join(m[:12] for m in matches[:10])
            raise ValueError(
                f"Ambiguous ID prefix '{partial_id}' matches "
                f"{len(matches)} structures: {preview}"
            )
        return matches[0]

    def count_structures(self):
        """Count unique structures in the database."""
        conn = sqlite3.connect(self.db_file)
        try:
            cursor = conn.execute("SELECT COUNT(*) FROM structures")
            return cursor.fetchone()[0]
        finally:
            conn.close()

    def get_structures_for_molecule(self, molecule_id):
        """Get all structures (conformers) belonging to a molecule."""
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT s.*, m.chemical_formula, m.empirical_formula,
                       m.number_of_atoms, m.mass, m.smiles
                FROM structures s
                JOIN molecules m ON s.molecule_id = m.molecule_id
                WHERE s.molecule_id = ?
                ORDER BY s.structure_id
                """,
                (molecule_id,),
            )
            rows = cursor.fetchall()
            return [self._structure_row_to_dict(dict(row)) for row in rows]
        finally:
            conn.close()

    def get_records_for_molecule(self, molecule_id):
        """Get all record summaries that reference a given molecule.
        Traverses: molecules -> structures -> record_structures -> records.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT DISTINCT r.record_index, r.record_id, r.program,
                       r.method, r.basis, r.jobtype, r.total_energy,
                       r.source_file
                FROM records r
                JOIN record_structures rs ON r.record_id = rs.record_id
                JOIN structures s ON rs.structure_id = s.structure_id
                WHERE s.molecule_id = ?
                ORDER BY r.record_index
                """,
                (molecule_id,),
            )
            rows = cursor.fetchall()
            return [dict(row) for row in rows]
        finally:
            conn.close()

    def get_records_for_structure(self, structure_id):
        """Get all record summaries that reference a given structure.
        Traverses: record_structures -> records.
        """
        conn = sqlite3.connect(self.db_file)
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute(
                """
                SELECT DISTINCT r.record_index, r.record_id, r.program,
                       r.method, r.basis, r.jobtype, r.total_energy,
                       r.source_file
                FROM records r
                JOIN record_structures rs ON r.record_id = rs.record_id
                WHERE rs.structure_id = ?
                ORDER BY r.record_index
                """,
                (structure_id,),
            )
            rows = cursor.fetchall()
            return [dict(row) for row in rows]
        finally:
            conn.close()

    @staticmethod
    def _molecule_row_to_dict(row):
        """Convert a molecules table row to a clean dictionary.
        Deserializes JSON fields and converts boolean integers to Python
        bools.
        """
        return {
            "molecule_id": row.get("molecule_id"),
            "molecule_label": row.get("molecule_label"),
            "empirical_formula": row.get("empirical_formula"),
            "chemical_formula": row.get("chemical_formula"),
            "smiles": row.get("smiles"),
            "inchi": row.get("inchi"),
            "chiral_centers": from_json(row.get("chiral_centers_json")),
            "number_of_atoms": row.get("number_of_atoms"),
            "mass": row.get("mass"),
            "elements": from_json(row.get("elements_json")),
            "element_counts": from_json(row.get("element_counts_json")),
            "is_chiral": bool(row.get("is_chiral")),
            "is_ring": bool(row.get("is_ring")),
            "is_aromatic": bool(row.get("is_aromatic")),
            "is_monoatomic": bool(row.get("is_monoatomic")),
            "is_diatomic": bool(row.get("is_diatomic")),
            "is_linear": bool(row.get("is_linear")),
            "is_multicomponent": bool(row.get("is_multicomponent")),
            "num_components": row.get("num_components"),
        }

    @staticmethod
    def _structure_row_to_dict(row):
        """Convert a structures table row (joined with molecules) to a
        clean dictionary.
        Deserializes JSON fields (chemical_symbols, positions, etc.) and
        includes molecule-level convenience fields when present in the row.
        """
        d = {
            "structure_id": row.get("structure_id"),
            "molecule_id": row.get("molecule_id"),
            "charge": row.get("charge"),
            "multiplicity": row.get("multiplicity"),
            "structure_label": row.get("structure_label"),
            "chemical_symbols": from_json(row.get("chemical_symbols_json")),
            "positions": from_json(row.get("positions_json")),
            "center_of_mass": from_json(row.get("center_of_mass_json")),
            "moments_of_inertia": from_json(
                row.get("moments_of_inertia_json")
            ),
        }
        # Include molecule-level fields when available (from JOIN)
        for key in (
            "chemical_formula",
            "empirical_formula",
            "number_of_atoms",
            "mass",
            "smiles",
        ):
            if key in row:
                d[key] = row[key]
        return d
