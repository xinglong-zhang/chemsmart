import csv
import json
import sqlite3

import numpy as np
import pytest

from chemsmart.database.assemble import SingleFileAssembler
from chemsmart.database.database import Database
from chemsmart.database.export import DatabaseExporter
from chemsmart.database.inspect import DatabaseInspector
from chemsmart.database.query import DatabaseQuery
from chemsmart.database.records import AssembledRecord
from chemsmart.database.utils import (
    SCHEMA_VERSION,
    bool_to_str,
    check_schema_version,
    compute_trajectory_id,
    convert_numpy,
    format_energy,
    format_float,
    format_kv,
    from_json,
    get_record_id,
    human_size,
    is_chemsmart_database,
    is_custom_basis,
    is_custom_solvent,
    resolve_record,
    sort_frames_by_energy,
    sort_structure_dicts_by_energy,
    standardize_basis_set,
    to_json,
    truncate_iso,
)

# Canonical InChIKeys (molecule_id) for molecules used across the test suite
INCHIKEY_H2O = "XLYOFNOQVPJJNP-UHFFFAOYSA-N"
INCHIKEY_HE = "SWQJXJOGLNCZEY-UHFFFAOYSA-N"
INCHIKEY_CO2 = "CURLTUGMZLYLDI-UHFFFAOYSA-N"
# Known stable record ID: Gaussian MP2/aug-cc-pVTZ geometry optimisation of water
RECORD_ID_GAUSSIAN_MP2_WATER = (
    "f9b9f1b79a7bff9a29ba999e0f625875f8ec3bc4162260ccf2a8117a9403e726"
)
# Known stable record ID: ORCA M062X/def2-SVP geometry optimisation of water
RECORD_ID_ORCA_M062X_WATER = (
    "57efe22eb7c24111c0e06f4b8ae3c6b8c38c7b916ba8375225d2cdd33a043a7e"
)
# Known stable structure ID: He atom at origin (0,0,0)
STRUCTURE_ID_ORIGIN_HE = (
    "4054373538ed8c659cd165de0c57822be13073206e4305ddc5dc6937fb2cd65b"
)


class TestDatabaseUtilities:
    def test_is_chemsmart_database(
        self, database_chemsmart_file, database_ase_file, database_empty_file
    ):
        assert is_chemsmart_database(database_chemsmart_file)
        assert not is_chemsmart_database(database_ase_file)
        assert not is_chemsmart_database(database_empty_file)

    def test_is_custom_basis(self):
        assert is_custom_basis(" GenECP ")
        assert is_custom_basis("gen")
        assert not is_custom_basis("6-31g")
        assert not is_custom_basis("cc-pvtz")
        assert not is_custom_basis(None)

    def test_is_custom_solvent(self):
        assert is_custom_solvent("generic,read")
        assert not is_custom_solvent("water")
        assert not is_custom_solvent("methanol")
        assert not is_custom_solvent(None)

    def test_get_record_id(self):
        id0 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
        )
        # different program -> different record
        id1 = get_record_id(
            program="orca",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
        )
        assert id0 != id1
        # different job type -> different record
        id2 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="opt",
            trajectory_id="traj_abc",
        )
        assert id0 != id2
        # different trajectory_id -> different record
        id3 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_adc",
        )
        assert id0 != id3
        # different basis -> different record
        id4 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="def2svp",
            jobtype="sp",
            trajectory_id="traj_abc",
        )
        assert id0 != id4
        # different basis (custom) -> different record
        id4_cb = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="customized_basis",
            jobtype="sp",
            trajectory_id="traj_abc",
            custom_basis_hash="a" * 64,
        )
        assert id4 != id4_cb
        # different custom_basis_hash -> different record
        id4_cb2 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="customized_basis",
            jobtype="sp",
            trajectory_id="traj_abc",
            custom_basis_hash="b" * 64,
        )
        assert id4_cb != id4_cb2
        # no solvent vs with solvent -> different record
        id5 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            solvent_model="pcm",
            solvent_id="water",
        )
        assert id0 != id5
        # different solvent model -> different record
        id5_sm = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            solvent_model="smd",
            solvent_id="water",
        )
        assert id5 != id5_sm
        # different solvent id -> different record
        id5_si = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            solvent_model="pcm",
            solvent_id="toluene",
        )
        assert id5 != id5_si
        # different solvent id (custom) -> different record
        id5_cs = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            solvent_model="pcm",
            solvent_id="customized_solvent",
            custom_solvent_hash="c" * 64,
        )
        assert id5 != id5_cs
        # different custom_solvent_hash -> different record
        id5_cs2 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            solvent_model="pcm",
            solvent_id="customized_solvent",
            custom_solvent_hash="d" * 64,
        )
        assert id5_cs != id5_cs2
        # idempotency: same inputs always produce the same id
        id0_repeat = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
        )
        assert id0 == id0_repeat

        assert all(
            len(i) == 64
            for i in [
                id0,
                id1,
                id2,
                id3,
                id4,
                id4_cb,
                id4_cb2,
                id5,
                id5_sm,
                id5_si,
                id5_cs,
                id5_cs2,
                id0_repeat,
            ]
        )

    def test_compute_trajectory_id(self):
        sid_a = "a" * 64
        sid_b = "b" * 64
        sid_c = "c" * 64
        t1 = compute_trajectory_id([sid_a, sid_b])
        # different structure_ids -> different trajectory
        t2 = compute_trajectory_id([sid_a, sid_c])
        assert t1 != t2
        # order matters: [A, B] != [B, A]
        t_reversed = compute_trajectory_id([sid_b, sid_a])
        assert t1 != t_reversed
        # idempotency: same inputs always yield the same id
        assert compute_trajectory_id([sid_a, sid_b]) == t1

    def test_truncate_iso(self):
        assert (
            truncate_iso("2026-05-13T08:47:52.847497+00:00")
            == "2026-05-13 08:47"
        )
        assert truncate_iso("2001-11-01T04:52:56+00:00") == "2001-11-01 04:52"
        assert truncate_iso(None) is None

    def test_human_size(self):
        assert human_size(0) == "0.0 B"
        assert human_size(1024) == "1.0 KB"
        assert human_size(1048576) == "1.0 MB"
        assert human_size(None) == "-"

    def test_convert_numpy(self):
        assert convert_numpy(np.int64(2)) == 2  # int
        assert convert_numpy(np.float64(1.23456789)) == 1.23456789  # float
        assert convert_numpy(np.array([1, 2])) == [1, 2]  # array
        assert convert_numpy((np.int64(3),)) == [3]  # tuple
        obj1 = {"a": np.int64(1), "b": np.float64(2.5)}
        assert convert_numpy(obj1) == {"a": 1, "b": 2.5}
        obj2 = {
            "a": np.array([1, 2]),
            "b": [np.int64(3), np.float64(4.5)],
            "c": {"d": np.int64(5)},
        }
        assert convert_numpy(obj2) == {
            "a": [1, 2],
            "b": [3, 4.5],
            "c": {"d": 5},
        }
        obj3 = [np.array([1, 2]), np.array([3, 4])]
        assert convert_numpy(obj3) == [[1, 2], [3, 4]]

    def test_convert_json(self):
        obj1 = {"a": 1, "b": 2}
        json_str = to_json(obj1)
        assert json.loads(json_str) == obj1
        assert to_json(None) is None
        json_str = '{"a": 1, "b": 2}'
        obj2 = from_json(json_str)
        assert obj2 == {"a": 1, "b": 2}
        assert from_json(None) is None

    def test_format(self):
        assert format_kv("Energy", None, 10) == "  Energy    : NULL"
        assert format_kv("Mass", 12.3456789, 10) == "  Mass      : 12.3456789"
        assert format_energy(-100.1234567800) == "-100.12345678"
        assert format_energy(None) == "NULL"
        assert format_float(1.23456789) == "1.234568"
        assert format_float(1.23456789, decimals=2) == "1.23"
        assert format_float(None) == "NULL"

    def test_bool_to_str(self):
        assert bool_to_str(True) == "Yes"
        assert bool_to_str(False) == "No"
        assert bool_to_str(1) == "Yes"
        assert bool_to_str(0) == "No"
        assert bool_to_str(None) == "NULL"

    def test_standardize_basis(self):
        assert standardize_basis_set("def2-svp") == "def2svp"
        assert standardize_basis_set("def2-tzvp") == "def2tzvp"
        assert standardize_basis_set("6-31g") == "6-31g"

    def test_sort_frames_by_energy(self):
        frames = [
            # Majority bucket: B3LYP/def2svp => sorted first, ascending energy.
            {"structure_id": "c", "energies": [("B3LYP", "def2svp", -1.0)]},
            {"structure_id": "a", "energies": [("B3LYP", "def2svp", -3.0)]},
            {"structure_id": "b", "energies": [("B3LYP", "def2svp", -2.0)]},
            # Minority bucket: no B3LYP/def2svp energy => pushed to the end.
            {"structure_id": "e", "energies": [("PBE0", "def2svp", -10.0)]},
            {"structure_id": "d", "energies": [("M062X", "def2svp", -20.0)]},
        ]
        sorted_frames_auto = sort_frames_by_energy(frames)
        # bucket 0 (ascending B3LYP): a(-3.0) -> b(-2.0) -> c(-1.0)
        # bucket 1 (ascending fallback): d(-20.0) -> e(-10.0)
        assert [f["structure_id"] for f in sorted_frames_auto] == [
            "a",
            "b",
            "c",
            "d",
            "e",
        ]
        # Primary (method, basis) entry must be moved to the front within each frame.
        assert sorted_frames_auto[0]["energies"][0] == (
            "B3LYP",
            "def2svp",
            -3.0,
        )
        sorted_frames_specified = sort_frames_by_energy(
            frames, primary=("M062X", "def2svp")
        )
        assert [f["structure_id"] for f in sorted_frames_specified] == [
            "d",
            "e",
            "a",
            "b",
            "c",
        ]

        # No-energy case: original order preserved, no IndexError.
        empty_frames = [{"structure_id": "x", "energies": []}]
        assert sort_frames_by_energy(empty_frames) == empty_frames

    def test_sort_structure_dicts_no_energy(self, tmp_path):
        db = Database(str(tmp_path / "empty.db"))
        db.create()
        no_energy_dicts = [
            {"structure_id": "aaaa"},
            {"structure_id": "bbbb"},
        ]
        result = sort_structure_dicts_by_energy(db.db_file, no_energy_dicts)
        assert [s["structure_id"] for s in result] == ["aaaa", "bbbb"]
        assert result[0]["primary_energy"] is None
        assert result[0]["primary_method_basis"] is None

    def test_resolve_record_helpers(
        self,
        tmp_path,
        gaussian_co2_opt_outfile,
        orca_co2_output,
    ):
        gaussian = SingleFileAssembler(gaussian_co2_opt_outfile).assemble_data
        orca = SingleFileAssembler(orca_co2_output).assemble_data
        assert gaussian is not None
        assert orca is not None

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([gaussian, orca]) == 2

        # resolve_record by index
        record_by_index = resolve_record(
            db,
            record_index=1,
            return_list=False,
        )
        assert record_by_index["record_id"] == gaussian.record_id
        # resolve_record by partial ID
        record_by_id = resolve_record(
            db,
            record_id=gaussian.record_id[:12],
            return_list=False,
        )
        assert record_by_id["record_id"] == gaussian.record_id

        assert record_by_index["record_id"] == record_by_id["record_id"]


class TestDatabaseSchemaAndInsertion:
    def test_create_schema(self, tmp_path):
        db = Database(str(tmp_path / "created_without_suffix"))
        db.create()
        assert db.db_file.endswith(".db")
        assert db.exists()
        assert db.count_records() == 0
        assert db.count_molecules() == 0
        assert db.count_structures() == 0
        conn = db.get_connection()
        tables = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        indexes = {
            row[0]
            for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='index'"
            )
        }
        conn.close()
        required_tables = {
            "records",
            "molecules",
            "structures",
            "record_structures",
        }
        assert required_tables <= tables
        assert "idx_record_id" in indexes
        assert "idx_struct_molecule_id" in indexes

    def test_insert_deduplicates_molecules(
        self,
        tmp_path,
        gaussian_co2_opt_outfile,
        orca_co2_output,
        orca_he_output_freq,
    ):
        # Two CO2 calculations (Gaussian + ORCA) and one He calculation.
        # Although CO2 appears in two separate records, it must be stored as
        # a single canonical molecular species (deduplicated by InChIKey).
        co2_gaussian = SingleFileAssembler(
            gaussian_co2_opt_outfile
        ).assemble_data
        co2_orca = SingleFileAssembler(orca_co2_output).assemble_data
        he_orca = SingleFileAssembler(orca_he_output_freq).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([co2_gaussian, co2_orca, he_orca]) == 3

        # Three records, but only TWO distinct molecular species: CO2 and He.
        assert db.count_records() == 3
        assert db.count_molecules() == 2

        # Retrieve both canonical molecules directly by InChIKey.
        co2_molecule = db.get_molecule(INCHIKEY_CO2)
        he_molecule = db.get_molecule(INCHIKEY_HE)
        assert co2_molecule is not None
        assert he_molecule is not None
        assert co2_molecule["chemical_formula"] == "CO2"
        assert he_molecule["chemical_formula"] == "He"

        # The two CO2 records must share exactly one molecule entry.
        # get_records_for_molecule proves that both computed records are
        # linked back to the same deduplicated species.
        co2_records = db.get_records_for_molecule(INCHIKEY_CO2)
        assert len(co2_records) == 2
        assert {r["program"] for r in co2_records} == {"gaussian", "orca"}

        # He has only one record and its own molecule entry.
        he_records = db.get_records_for_molecule(INCHIKEY_HE)
        assert len(he_records) == 1

        # All conformer structures from both CO2 calculations are preserved
        # and linked to the single shared CO2 molecule_id.
        co2_structures = db.get_structures_for_molecule(INCHIKEY_CO2)
        assert len(co2_structures) > 0
        assert all(s["molecule_id"] == INCHIKEY_CO2 for s in co2_structures)

    def test_insert_replaces_record(
        self, tmp_path, orca_he_output_freq, orca_he_output_freq_new
    ):
        # Both fixtures are the same He M062X/def2-SVP freq calculation run on
        # ORCA 6.0.1 and ORCA 6.1.1 respectively. Because the molecule,
        # method, basis, and job type are identical the assembler produces the
        # same record_id for both — no artificial ID override is needed.
        # The energies differ slightly, making the replace observable.
        he_v1 = SingleFileAssembler(orca_he_output_freq).assemble_data
        he_v2 = SingleFileAssembler(orca_he_output_freq_new).assemble_data
        assert he_v1.record_id == he_v2.record_id  # same record_id
        assert he_v1.results["total_energy"] != he_v2.results["total_energy"]

        db = Database(str(tmp_path / "replace.db"))
        db.create()

        # initial insert ORCA 6.0.1 result
        db.insert_record(he_v1)
        assert db.count_records() == 1
        energy_before = db.get_record(record_id=he_v1.record_id)["results"][
            "total_energy"
        ]
        assert energy_before == he_v1.results["total_energy"]

        # re-insert same record_id (INSERT OR REPLACE)
        db.insert_record(he_v2)

        # Record count must remain 1 — no duplicate was created.
        assert db.count_records() == 1

        # Stored energy must now reflect the newer run; old value is gone.
        energy_after = db.get_record(record_id=he_v1.record_id)["results"][
            "total_energy"
        ]
        assert energy_after == he_v2.results["total_energy"]
        assert energy_after != energy_before

        # The record-structure link table must be fully refreshed — no stale
        # entries from the first insert remain.
        structures = db.get_structures_for_record(he_v1.record_id)
        assert [s["structure_id"] for s in structures] == [
            mol["structure_id"] for mol in he_v2.molecules
        ]

    def test_insert_skips_failed_record(
        self,
        tmp_path,
        gaussian_benzene_opt_outfile,
        gaussian_link_failed_outfile,
    ):
        # A real-world import batch often mixes successful parses with files
        # that terminated abnormally. SingleFileAssembler returns None for
        # the latter; insert_records must skip None entries and still commit
        # the valid records.
        succeeded = SingleFileAssembler(
            gaussian_benzene_opt_outfile
        ).assemble_data
        failed = SingleFileAssembler(
            gaussian_link_failed_outfile
        ).assemble_data
        assert failed is None  # error-terminated log -> assembler skips it

        db = Database(str(tmp_path / "partial_insert.db"))
        db.create()
        assert db.insert_records([succeeded, failed]) == 1
        assert db.count_records() == 1

    def test_schema_version_set_on_create(self, tmp_path):
        """create() writes SCHEMA_VERSION into PRAGMA user_version."""
        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        conn = db.get_connection()
        version = conn.execute("PRAGMA user_version").fetchone()[0]
        conn.close()
        assert version == SCHEMA_VERSION

    def test_schema_version_mismatch_raises(self, tmp_path):
        """check_schema_version() raises RuntimeError for a stale database."""
        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        # Simulate a database written by a future (incompatible) chemsmart
        conn = sqlite3.connect(db.db_file)
        conn.execute(f"PRAGMA user_version = {SCHEMA_VERSION + 1}")
        conn.commit()
        conn.close()
        with pytest.raises(RuntimeError, match="schema version mismatch"):
            check_schema_version(db.db_file)

    def test_foreign_key_enforced(self, tmp_path):
        """open_connection() enables FK enforcement: orphan inserts are rejected."""
        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        conn = db.get_connection()
        with pytest.raises(sqlite3.IntegrityError):
            conn.execute(
                "INSERT INTO record_structures "
                "(record_id, structure_id, index_in_record) "
                "VALUES ('nonexistent_record', 'nonexistent_structure', 0)"
            )
        conn.close()


class TestDatabaseRecordMoleculeStructureQueries:
    def test_get_record_summaries(
        self,
        tmp_path,
        gaussian_benzene_opt_outfile,
        orca_co2_output,
        gaussian_he_opt_outfile,
    ):
        benzene_gaussian = SingleFileAssembler(
            gaussian_benzene_opt_outfile
        ).assemble_data
        co2_orca = SingleFileAssembler(orca_co2_output).assemble_data
        he_gaussian = SingleFileAssembler(
            gaussian_he_opt_outfile
        ).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert (
            db.insert_records([benzene_gaussian, co2_orca, he_gaussian]) == 3
        )

        # Summaries should handle mixed molecules and both supported programs.
        summaries = db.get_all_record_summaries()
        assert len(summaries) == 3
        assert {s["program"] for s in summaries} == {"gaussian", "orca"}
        assert {s["chemical_formula"] for s in summaries} == {
            "C6H6",
            "CO2",
            "He",
        }

        benzene_summary = next(
            s for s in summaries if s["chemical_formula"] == "C6H6"
        )
        summary = db.get_record_summary(record_id=benzene_summary["record_id"])
        assert summary["charge"] == 0
        assert summary["multiplicity"] == 1

    def test_molecule_structure_accessors(
        self,
        tmp_path,
        gaussian_mp2_outputfile,
        water_output_gas_path,
        orca_he_output_freq,
    ):
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data
        he_orca = SingleFileAssembler(orca_he_output_freq).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([water_gaussian, water_orca, he_orca]) == 3

        molecules = db.get_all_molecules()
        assert len(molecules) == 2
        assert db.get_molecule("missing") is None
        h2o_molecule = db.get_molecule(INCHIKEY_H2O)
        assert h2o_molecule["element_counts"] == {"H": 2, "O": 1}
        assert h2o_molecule["is_linear"] is False

        structures = db.get_all_structures()
        h2o_structures = db.get_structures_for_molecule(INCHIKEY_H2O)
        he_structures = db.get_structures_for_molecule(INCHIKEY_HE)
        assert len(structures) == len(h2o_structures) + len(he_structures)
        assert db.get_structure("missing") is None
        struct = db.get_structure(STRUCTURE_ID_ORIGIN_HE)
        assert struct["molecule_id"] == INCHIKEY_HE
        assert struct["positions"] == [[0.0, 0.0, 0.0]]

        assert len(db.get_records_for_molecule(INCHIKEY_H2O)) == 2
        assert (
            len(
                db.get_records_for_structure(h2o_structures[0]["structure_id"])
            )
            >= 1
        )

    def test_partial_id_resolution(
        self,
        tmp_path,
        gaussian_mp2_outputfile,
        orca_co2_output,
        orca_he_output_freq,
    ):
        """Verify that partial ID matching works for records, molecules, and structures."""
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        co2_orca = SingleFileAssembler(orca_co2_output).assemble_data
        he_orca = SingleFileAssembler(orca_he_output_freq).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([water_gaussian, co2_orca, he_orca]) == 3

        # A prefix should be unique in this small database
        assert (
            db.get_record_by_partial_id(RECORD_ID_GAUSSIAN_MP2_WATER[:12])
            == RECORD_ID_GAUSSIAN_MP2_WATER
        )
        assert db.get_molecule_by_partial_id(INCHIKEY_CO2[:16]) == INCHIKEY_CO2
        assert (
            db.get_structure_by_partial_id(STRUCTURE_ID_ORIGIN_HE[:12])
            == STRUCTURE_ID_ORIGIN_HE
        )


class TestDatabaseQuery:
    def test_parse_query_validation(self, tmp_path):
        db = Database(str(tmp_path / "empty.db"))
        db.create()

        # Test parsing compound WHERE conditions
        dq = DatabaseQuery(
            db.db_file,
            "total_energy < -3 AND program = 'orca'",
        )
        where, params = dq.parse_query()
        assert "r.total_energy < ?" in where
        assert "AND" in where
        assert params == (-3, "orca")

        # Test parsing LIKE pattern
        like = DatabaseQuery(db.db_file, 'source_file ~ "CO2"')
        where, params = like.parse_query()
        assert params == ("%CO2%",)

        with pytest.raises(ValueError, match="Invalid target"):
            DatabaseQuery(db.db_file, None, target="atoms")
        with pytest.raises(ValueError, match="No query string"):
            DatabaseQuery(db.db_file, None).parse_query()
        with pytest.raises(ValueError, match="Unknown field"):
            DatabaseQuery(db.db_file, "badfield = 1").parse_query()
        with pytest.raises(ValueError, match="Invalid condition"):
            DatabaseQuery(db.db_file, "not a condition").parse_query()

    def test_query_summaries_and_counts(
        self,
        tmp_path,
        gaussian_ozone_opt_outfile,
        orca_co2_output,
        gaussian_he_opt_outfile,
    ):
        ozone_gaussian = SingleFileAssembler(
            gaussian_ozone_opt_outfile
        ).assemble_data
        co2_orca = SingleFileAssembler(orca_co2_output).assemble_data
        he_gaussian = SingleFileAssembler(
            gaussian_he_opt_outfile
        ).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([ozone_gaussian, co2_orca, he_gaussian]) == 3

        # Test query limit
        records = DatabaseQuery(db.db_file, None, limit=2)
        assert len(records.query_summaries()) == 2
        assert records.count_total() == 3
        assert records.count_matched() == 3

        # Test filtering by program with the single ORCA record in this set.
        matched = DatabaseQuery(db.db_file, "program = 'orca'")
        assert matched.count_matched() == 1
        assert matched.query_summaries()[0]["program"] == "orca"

        # Test molecule-level queries: 3-atom molecules
        molecules = DatabaseQuery(
            db.db_file,
            "number_of_atoms = 3",
            target="molecules",
        )
        assert molecules.count_total() == 3  # O3, CO2, He
        assert molecules.count_matched() == 2  # O3 and CO2
        formulas = {m["chemical_formula"] for m in molecules.query_summaries()}
        assert formulas == {"O3", "CO2"}

        # Test structure-level queries: neutral singlets
        structures = DatabaseQuery(
            db.db_file,
            "charge = 0 AND multiplicity = 1",
            target="structures",
        )
        assert (
            structures.count_matched() == structures.count_total()
        )  # all structures are neutral singlets
        summaries = structures.query_summaries()
        assert summaries
        assert len(summaries) == structures.count_matched()
        formatted = structures.format_summary(summaries)
        assert "Query Summary (structures)" in formatted
        assert summaries[0]["structure_id"][:12] in formatted

        # Test empty result formatting
        empty = DatabaseQuery(db.db_file, "method = 'never'")
        empty_summaries = empty.query_summaries()
        assert len(empty_summaries) == 0
        assert "No records" in empty.format_summary(empty.query_summaries())


class TestDatabaseExport:
    def test_json_and_csv_export(
        self,
        tmp_path,
        gaussian_mp2_outputfile,
        water_output_gas_path,
    ):
        """Test JSON and CSV export with two water records from different programs."""
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([water_gaussian, water_orca]) == 2

        # Test JSON export
        json_path = tmp_path / "records.json"
        DatabaseExporter(db.db_file, str(json_path)).to_json()
        records = json.loads(json_path.read_text(encoding="utf-8"))
        assert len(records) == 2
        # Both records are for H2O molecule
        formulas = {
            record["molecules"][0]["chemical_formula"] for record in records
        }
        assert formulas == {"H2O"}
        # One from Gaussian, one from ORCA
        programs = {record["provenance"]["program"] for record in records}
        assert programs == {"gaussian", "orca"}

        # Test CSV export
        csv_path = tmp_path / "records.csv"
        exporter = DatabaseExporter(
            db.db_file,
            str(csv_path),
            keys="program,method,basis,total_energy,charge",
        )
        exporter.to_csv()
        with open(csv_path, newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 2
        # Verify all rows have expected structure
        csv_formulas = {row["chemical_formula"] for row in rows}
        csv_programs = {row["program"] for row in rows}
        csv_charges = {row["charge"] for row in rows}
        assert csv_formulas == {"H2O"}
        assert csv_programs == {"gaussian", "orca"}
        assert csv_charges == {"0"}  # both are neutral

    def test_export_helpers(self, tmp_path):
        # Create empty database for validation tests
        db = Database(str(tmp_path / "empty.db"))
        db.create()
        # Test unsupported output format validation
        with pytest.raises(ValueError, match="Unsupported output format"):
            DatabaseExporter(db.db_file, str(tmp_path / "unsupported.txt"))
        # Test unsupported CSV key validation
        with pytest.raises(ValueError, match="Unsupported CSV key"):
            DatabaseExporter(
                db.db_file,
                str(tmp_path / "output.csv"),
                keys="unsupported_key",
            )
        # Test record_to_csv_row helper
        row = DatabaseExporter.record_to_csv_row(
            {"record_index": 1, "record_id": "r", "molecules": []},
            ["record_index", "record_id", "chemical_formula", "method"],
        )
        assert row == {
            "record_index": 1,
            "record_id": "r",
            "chemical_formula": None,
            "method": None,
        }

    def test_xyz_export_selectors(
        self,
        tmp_path,
        gaussian_benzene_opt_outfile,
    ):
        benzene_gaussian = SingleFileAssembler(
            gaussian_benzene_opt_outfile
        ).assemble_data

        db = Database(str(tmp_path / "benzene_gaussian.db"))
        db.create()
        assert db.insert_records([benzene_gaussian]) == 1

        record = db.get_record(record_id=benzene_gaussian.record_id)
        last_structure = record["molecules"][-1]
        atom_count = len(last_structure["chemical_symbols"])

        # Export by record_id: should export the last structure of the record, starting with atom count
        record_xyz = tmp_path / "record.xyz"
        DatabaseExporter(
            db.db_file,
            str(record_xyz),
            record_id=benzene_gaussian.record_id[:12],
        ).to_xyz()
        record_lines = record_xyz.read_text(encoding="utf-8").splitlines()
        assert record_lines[0] == str(atom_count)

        # Export by structure_id: should export that specific structure, including SID in comment line
        structure_xyz = tmp_path / "structure.xyz"
        DatabaseExporter(
            db.db_file,
            str(structure_xyz),
            structure_id=last_structure["structure_id"][:12],
        ).to_xyz()
        structure_text = structure_xyz.read_text(encoding="utf-8")
        assert f"SID: {last_structure['structure_id'][:12]}" in structure_text

        # Export by molecule_id: should export all structures for the canonical molecule
        molecule_xyz = tmp_path / "molecule.xyz"
        DatabaseExporter(
            db.db_file,
            str(molecule_xyz),
            molecule_id=last_structure["molecule_id"],
        ).to_xyz()
        molecule_lines = molecule_xyz.read_text(encoding="utf-8").splitlines()
        assert (
            molecule_lines.count(str(atom_count)) >= 1
        )  # at least one structure block with correct atom count

        # No selector
        with pytest.raises(ValueError, match="requires exactly one"):
            DatabaseExporter(db.db_file, str(tmp_path / "bad.xyz")).to_xyz()
        # Multiple selectors
        with pytest.raises(ValueError, match="requires exactly one"):
            DatabaseExporter(
                db.db_file,
                str(tmp_path / "bad2.xyz"),
                record_id=benzene_gaussian.record_id[:12],
                structure_id=last_structure["structure_id"][:12],
            ).to_xyz()
        # structure_index out of range
        with pytest.raises(ValueError, match="out of range"):
            DatabaseExporter(
                db.db_file,
                str(tmp_path / "bad3.xyz"),
                record_id=benzene_gaussian.record_id[:12],
                structure_index="999",
            ).to_xyz()


class TestDatabaseInspect:
    def test_inspect_overview_and_record(
        self,
        tmp_path,
        gaussian_mp2_outputfile,
        water_output_gas_path,
    ):
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([water_gaussian, water_orca]) == 2

        # Overview should report two records but only one molecule: H2O.
        inspector = DatabaseInspector(db.db_file)
        stats = inspector.overview()
        assert stats["num_records"] == 2
        assert stats["num_molecules"] == 1
        assert ("gaussian", 1) in stats["programs"]
        assert ("orca", 1) in stats["programs"]
        assert "Database Overview" in inspector.format_overview()

        record_inspector = DatabaseInspector(
            db.db_file, record_id=water_gaussian.record_id[:12]
        )
        assert record_inspector.resolve_id() == water_gaussian.record_id
        detail = record_inspector.record_detail()
        assert (
            detail["record_id"]
            == water_gaussian.record_id
            == RECORD_ID_GAUSSIAN_MP2_WATER
        )
        formatted = record_inspector.format_record_detail()
        assert "Record Detail" in formatted
        assert "Electronic Results" in formatted

    def test_inspect_related_entities(
        self,
        tmp_path,
        gaussian_mp2_outputfile,
        water_output_gas_path,
        orca_he_output_freq,
    ):
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data
        he_orca = SingleFileAssembler(orca_he_output_freq).assemble_data

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        assert db.insert_records([water_gaussian, water_orca, he_orca]) == 3

        water_orca_record = db.get_record(record_id=water_orca.record_id)
        structure_inspector = DatabaseInspector(
            db.db_file,
            record_id=water_orca.record_id[:12],
            structure_index=1,
        )
        record, structure = structure_inspector.structure_detail()
        assert (
            record["record_id"]
            == water_orca.record_id
            == RECORD_ID_ORCA_M062X_WATER
        )
        assert (
            structure["structure_id"]
            == water_orca_record["molecules"][0]["structure_id"]
        )
        formatted_standalone_structure = (
            structure_inspector.format_structure_detail()
        )
        assert "Coordinates" in formatted_standalone_structure
        assert structure["structure_id"][:12] in formatted_standalone_structure
        # Invalid structure_index should fail with a clear error.
        with pytest.raises(ValueError, match="Structure index 999"):
            DatabaseInspector(
                db.db_file,
                record_id=water_orca.record_id[:12],
                structure_index=999,
            ).structure_detail()

        molecule_inspector = DatabaseInspector(
            db.db_file, molecule_id=INCHIKEY_H2O
        )
        molecule, structures, records = molecule_inspector.molecule_detail()
        assert molecule["molecule_id"] == INCHIKEY_H2O
        assert molecule["chemical_formula"] == "H2O"
        assert len(structures) > 0
        assert len(records) == 2
        formatted_molecule = molecule_inspector.format_molecule_detail()
        assert "Molecule Detail" in formatted_molecule
        assert "H2O" in formatted_molecule

        standalone_structure_inspector = DatabaseInspector(
            db.db_file,
            structure_id=he_orca.molecules[0]["structure_id"][:12],
        )
        standalone_structure, related_records = (
            standalone_structure_inspector.standalone_structure_detail()
        )
        assert (
            standalone_structure["structure_id"]
            == STRUCTURE_ID_ORIGIN_HE
            == he_orca.molecules[0]["structure_id"]
        )
        assert len(related_records) == 1
        assert related_records[0]["record_id"] == he_orca.record_id
        formatted_standalone_structure = (
            standalone_structure_inspector.format_standalone_structure_detail()
        )
        assert "Structure Detail" in formatted_standalone_structure
        assert (
            he_orca.molecules[0]["structure_id"][:12]
            in formatted_standalone_structure
        )


class TestDatabaseAssembler:
    def test_assembler_builds_records(
        self,
        gaussian_mp2_outputfile,
        water_output_gas_path,
    ):
        water_gaussian = SingleFileAssembler(
            gaussian_mp2_outputfile
        ).assemble_data
        assert isinstance(water_gaussian, AssembledRecord)
        assert water_gaussian.provenance["program"] == "gaussian"
        assert water_gaussian.provenance["parser"] == "Gaussian16Output"
        assert water_gaussian.meta["method"] == "mp2"
        assert water_gaussian.meta["basis"] == "aug-cc-pvtz"
        assert water_gaussian.meta["jobtype"] == "opt"
        assert water_gaussian.results["total_energy"] == -76.328992324258
        assert water_gaussian.molecules[-1]["chemical_formula"] == "H2O"

        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data
        assert isinstance(water_orca, AssembledRecord)
        assert water_orca.provenance["program"] == "orca"
        assert water_orca.provenance["parser"] == "ORCAOutput"
        assert water_orca.meta["method"] == "m062x"
        assert water_orca.meta["basis"] == "def2svp"
        assert water_orca.meta["jobtype"] == "opt"
        assert water_orca.results["total_energy"] == -76.323311011349
        assert water_orca.molecules[-1]["chemical_formula"] == "H2O"

    def test_assembler_skips_invalid_files(
        self, tmp_path, gaussian_link_failed_outfile
    ):
        assert (
            SingleFileAssembler(gaussian_link_failed_outfile).assemble_data
            is None
        )
        unknown_output = tmp_path / "unknown.out"
        unknown_output.write_text(
            "not a quantum chemistry output", encoding="utf-8"
        )
        with pytest.raises(ValueError, match="Unsupported format"):
            SingleFileAssembler(unknown_output)._get_assembler(unknown_output)


class TestStaticDatabaseContent:
    def test_database_counts(self, database_chemsmart_file):
        """The fixture holds exactly 40 records, 30 molecules, 303 structures."""
        db = Database(database_chemsmart_file)
        assert db.count_records() == 40
        assert db.count_molecules() == 30
        assert db.count_structures() == 303

    def test_program_distribution(self, database_chemsmart_file):
        """30 Gaussian records and 10 ORCA records are present."""
        db = Database(database_chemsmart_file)
        programs = {}
        for s in db.get_all_record_summaries():
            programs[s["program"]] = programs.get(s["program"], 0) + 1
        assert programs["gaussian"] == 30
        assert programs["orca"] == 10

    def test_gaussian_mp2_water_record(
        self, database_chemsmart_file, gaussian_mp2_outputfile
    ):
        """Gaussian MP2/aug-cc-pVTZ opt of water (water_mp2.log)."""
        db = Database(database_chemsmart_file)
        record = db.get_record(record_id=RECORD_ID_GAUSSIAN_MP2_WATER)
        assembled = SingleFileAssembler(gaussian_mp2_outputfile).assemble_data

        assert record["record_id"] == assembled.record_id
        assert (
            record["provenance"]["program"]
            == assembled.provenance["program"]
            == "gaussian"
        )
        assert record["meta"]["method"] == assembled.meta["method"] == "mp2"
        assert (
            record["meta"]["basis"] == assembled.meta["basis"] == "aug-cc-pvtz"
        )
        assert record["meta"]["jobtype"] == assembled.meta["jobtype"] == "opt"
        assert (
            record["meta"]["solvent_on"]
            == assembled.meta["solvent_on"]
            is False
        )
        assert np.isclose(
            record["results"]["total_energy"],
            assembled.results["total_energy"],
        )
        assert np.isclose(
            record["results"]["total_energy"], -76.328992324258, rtol=1e-6
        )
        assert np.isclose(
            record["results"]["homo_energy"], assembled.results["homo_energy"]
        )
        assert np.isclose(
            record["results"]["homo_energy"], -13.881616466470707, rtol=1e-6
        )
        assert (
            record["molecules"][-1]["chemical_formula"]
            == assembled.molecules[-1]["chemical_formula"]
            == "H2O"
        )
        assert (
            record["molecules"][-1]["chemical_symbols"]
            == assembled.molecules[-1]["chemical_symbols"]
            == ["O", "H", "H"]
        )
        assert (
            record["molecules"][-1]["charge"]
            == assembled.molecules[-1]["charge"]
            == 0
        )
        assert (
            record["molecules"][-1]["multiplicity"]
            == assembled.molecules[-1]["multiplicity"]
            == 1
        )
        assert len(record["molecules"]) == len(assembled.molecules) == 4

    def test_orca_m062x_water_record(
        self, database_chemsmart_file, water_output_gas_path
    ):
        """ORCA M062X/def2-SVP opt of water (water_opt.out)."""
        db = Database(database_chemsmart_file)
        record = db.get_record(record_id=RECORD_ID_ORCA_M062X_WATER)
        assembled = SingleFileAssembler(water_output_gas_path).assemble_data

        assert record["record_id"] == assembled.record_id
        assert (
            record["provenance"]["program"]
            == assembled.provenance["program"]
            == "orca"
        )
        assert record["meta"]["method"] == assembled.meta["method"] == "m062x"
        assert record["meta"]["basis"] == assembled.meta["basis"] == "def2svp"
        assert np.isclose(
            record["results"]["total_energy"],
            assembled.results["total_energy"],
        )
        assert np.isclose(
            record["results"]["total_energy"], -76.323311011349, rtol=1e-6
        )
        assert (
            record["molecules"][-1]["chemical_formula"]
            == assembled.molecules[-1]["chemical_formula"]
            == "H2O"
        )
        assert len(record["molecules"]) == len(assembled.molecules) == 5

    def test_h2o_molecule(self, database_chemsmart_file):
        """H2O molecule: composition, topology, and cross-record linkage.
        The database holds exactly 2 water records (Gaussian water_mp2.log
        + ORCA water_opt.out) that both map to the same canonical H2O
        molecular species.
        """
        db = Database(database_chemsmart_file)
        mol = db.get_molecule(INCHIKEY_H2O)
        assert mol is not None
        assert mol["chemical_formula"] == "H2O"
        assert mol["element_counts"] == {"O": 1, "H": 2}
        assert mol["is_linear"] is False
        assert mol["is_monoatomic"] is False
        records = db.get_records_for_molecule(INCHIKEY_H2O)
        assert len(records) == 2  # water_mp2.log and water_opt.out
        programs = {r["program"] for r in records}
        assert programs == {"gaussian", "orca"}

    def test_he_molecule(self, database_chemsmart_file):
        """He molecule: monoatomic, single geometry, two records.
        The database holds exactly 2 He calculations (Gaussian he.log + ORCA
        He.out) both use the same single-atom He geometry.
        """
        db = Database(database_chemsmart_file)
        mol = db.get_molecule(INCHIKEY_HE)
        assert mol is not None
        assert mol["chemical_formula"] == "He"
        assert mol["element_counts"] == {"He": 1}
        assert mol["number_of_atoms"] == 1
        assert mol["is_monoatomic"] is True
        records = db.get_records_for_molecule(INCHIKEY_HE)
        assert len(records) == 2  # he.log and He.out
        programs = {r["program"] for r in records}
        assert programs == {"gaussian", "orca"}

    def test_query_by_program(self, database_chemsmart_file):
        """Filtering by program='gaussian' returns exactly 30 of the 40 records."""
        dq = DatabaseQuery(database_chemsmart_file, "program = 'gaussian'")
        assert dq.count_total() == 40
        assert dq.count_matched() == 30

    def test_query_by_method_mp2(self, database_chemsmart_file):
        """Only one record uses plain MP2 — Gaussian water_mp2.log."""
        dq = DatabaseQuery(database_chemsmart_file, "method = 'mp2'")
        assert dq.count_matched() == 1
        result = dq.query_summaries()[0]
        assert result["program"] == "gaussian"
        assert result["basis"] == "aug-cc-pvtz"
        assert result["chemical_formula"] == "H2O"

    def test_query_limit(self, database_chemsmart_file):
        """limit= truncates results but count_total() still reflects the full DB."""
        dq = DatabaseQuery(database_chemsmart_file, None, limit=5)
        assert len(dq.query_summaries()) == 5
        assert dq.count_total() == 40

    def test_inspect_overview(self, database_chemsmart_file):
        """Inspector overview matches the manually verified DB-wide statistics."""
        inspector = DatabaseInspector(database_chemsmart_file)
        stats = inspector.overview()
        assert stats["num_records"] == 40
        assert stats["num_molecules"] == 30
        assert stats["num_structures"] == 303
        assert ("gaussian", 30) in stats["programs"]
        assert ("orca", 10) in stats["programs"]

    def test_inspect_mp2_water_record(self, database_chemsmart_file):
        """Inspector record detail matches the same ground truth as
        test_gaussian_mp2_water_record, but accessed via DatabaseInspector."""
        inspector = DatabaseInspector(
            database_chemsmart_file,
            record_id=RECORD_ID_GAUSSIAN_MP2_WATER[:12],
        )
        assert inspector.resolve_id() == RECORD_ID_GAUSSIAN_MP2_WATER
        detail = inspector.record_detail()
        assert detail["meta"]["method"] == "mp2"
        assert np.isclose(
            detail["results"]["total_energy"], -76.328992324258, rtol=1e-6
        )

    def test_inspect_molecule_h2o(self, database_chemsmart_file):
        """Molecule inspector for H2O reports 2 records that are directly
        verified from the source opt output files."""
        inspector = DatabaseInspector(
            database_chemsmart_file,
            molecule_id=INCHIKEY_H2O,
        )
        mol, structures, records = inspector.molecule_detail()
        assert mol["chemical_formula"] == "H2O"
        assert len(records) == 2

    def test_partial_id_for_mp2_water(self, database_chemsmart_file):
        """A 12-character prefix is sufficient to uniquely identify MP2 water
        in this 40-record database."""
        db = Database(database_chemsmart_file)
        resolved = db.get_record_by_partial_id(
            RECORD_ID_GAUSSIAN_MP2_WATER[:12]
        )
        assert resolved == RECORD_ID_GAUSSIAN_MP2_WATER
