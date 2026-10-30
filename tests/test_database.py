import csv
import json
import sqlite3
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemsmart.database.assemble import (
    BaseAssembler,
    GaussianAssembler,
    SingleFileAssembler,
    SingleFolderAssembler,
)
from chemsmart.database.database import Database
from chemsmart.database.export import DatabaseExporter
from chemsmart.database.inspect import DatabaseInspector
from chemsmart.database.query import DatabaseQuery
from chemsmart.database.records import AssembledRecord
from chemsmart.database.utils import (
    SCHEMA_VERSION,
    bool_to_str,
    canonicalize_route_string,
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
INCHIKEY_C6H4 = "AIESRBVWAFETPR-UHFFFAOYSA-N"
# Known stable record ID: Gaussian MP2/aug-cc-pVTZ geometry optimisation of water
RECORD_ID_GAUSSIAN_MP2_WATER = (
    "d69cb7535dbeff46026b5f1a786659cc1d9035ea2015c0c73404cc4f7781eb63"
)
# Known stable record ID: ORCA M062X/def2-SVP geometry optimisation of water
RECORD_ID_ORCA_M062X_WATER = (
    "645702e39834649153f4f8c83b8f979165cd61eb86c28d4507db8e01f76e42f7"
)
# Known stable record ID: xTB GFN2 ohess of water
RECORD_ID_XTB_GFN2_WATER = (
    "cfec83b1119920c5810c9989996017b1e82117b6bf8df7e5ecb5716119c1c126"
)
# Known stable structure ID: He atom at origin (0,0,0)
STRUCTURE_ID_ORIGIN_HE = (
    "4054373538ed8c659cd165de0c57822be13073206e4305ddc5dc6937fb2cd65b"
)
# Known stable record ID: xTB GFN2 ohess of CO2 (co2_ohess)
RECORD_ID_XTB_CO2 = (
    "ab59d330d3fa785dad8a23ac51f2bf4edde8569ace56347f262c716151f29598"
)
# Known stable record ID: xTB GFN2 hess of He (he_hess)
RECORD_ID_XTB_HE = (
    "5668063a53c661d8c07f6b56ee4985cf5740fd38674a7b0f563d98fa85986400"
)
# Known stable record ID: xTB GFN2 ALPB(toluene) opt of p-benzyne
RECORD_ID_XTB_P_BENZYNE_OPT = (
    "e213b34b678f7682ed077c0088dd97e9d229e721bf546321bc35003030f8558a"
)


def _make_inspector(**kwargs):
    with patch("chemsmart.database.inspect.Database") as mock_db_cls:
        mock_db_cls.return_value = MagicMock()
        inspector = DatabaseInspector("fake.db", **kwargs)
    return inspector


class FakeMolecule:
    def __init__(
        self,
        forces=None,
        num_atoms=2,
        chemical_formula="H2",
        chemical_symbols=("H", "H"),
        positions=None,
    ):
        self.forces = forces
        self.num_atoms = num_atoms
        self.chemical_formula = chemical_formula
        self.chemical_symbols = list(chemical_symbols)
        self.positions = (
            positions
            if positions is not None
            else np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])
        )


def _make_exporter(output="out.xyz", **kwargs):
    with patch("chemsmart.database.export.Database") as mock_db_cls:
        mock_db_cls.return_value = MagicMock()
        exporter = DatabaseExporter("fake.db", output, **kwargs)
    return exporter


def _patch_database_file(build_molecule=None):
    """Patch DatabaseFile at its definition site (frames_from_* methods
    import it locally, so patching chemsmart.database.export wouldn't
    affect the fresh import)."""
    mock_db_file_cls = MagicMock()
    mock_db_file_instance = MagicMock()
    if build_molecule is not None:
        mock_db_file_instance.build_molecule_from_database.side_effect = (
            build_molecule
        )
    mock_db_file_cls.return_value = mock_db_file_instance
    return patch("chemsmart.io.database.DatabaseFile", mock_db_file_cls)


def _make_base_assembler(output, include_failed=False):
    assembler = BaseAssembler.__new__(BaseAssembler)
    assembler.filename = "some_file.log"
    assembler.folder = None
    assembler.index = ":"
    assembler.include_failed = include_failed
    assembler.target = assembler.filename
    assembler.output = output
    return assembler


class TestDatabaseUtilities:
    def test_is_chemsmart_database(
        self, database_chemsmart_file, database_ase_file, database_empty_file
    ):
        assert is_chemsmart_database(database_chemsmart_file)
        assert not is_chemsmart_database(database_ase_file)
        assert not is_chemsmart_database(database_empty_file)

    def test_is_chemsmart_database_returns_false_for_non_db_extension(self):
        assert not is_chemsmart_database("not_a_database.txt")

    def test_is_chemsmart_database_returns_false_for_corrupt_file(
        self, tmp_path
    ):
        corrupt_db = tmp_path / "corrupt.db"
        corrupt_db.write_bytes(b"not a real sqlite file at all \x00\x01")
        assert not is_chemsmart_database(str(corrupt_db))

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
        # no route_hash vs with route_hash -> different record
        id6 = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            route_hash="c" * 64,
        )
        assert id0 != id6
        # different route_hash -> different record
        id6_rs = get_record_id(
            program="gaussian",
            method="b3lyp",
            basis="6-31g",
            jobtype="sp",
            trajectory_id="traj_abc",
            route_hash="d" * 64,
        )
        assert id6 != id6_rs
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
                id6,
                id6_rs,
                id0_repeat,
            ]
        )

    def test_canonicalize_route_string(self):
        route1 = "# opt=(ts,calcfc,noeigentest) b3lyp def2-svp"
        drop1 = ["b3lyp", "def2svp", "ts"]
        tokens1 = canonicalize_route_string(route1, drop_terms=drop1)
        assert "calcfc" in tokens1
        assert "noeigentest" in tokens1
        assert "opt" in tokens1
        for known in ("b3lyp", "def2svp", "ts"):
            assert known not in tokens1

        route2 = (
            "#p opt m06-2x genecp int=ultrafine scf=tight "
            "empiricaldispersion=gd3bj"
        )
        drop2 = ["m062x", "genecp", "opt"]
        tokens2 = canonicalize_route_string(route2, drop_terms=drop2)
        for modifier in (
            "int",
            "ultrafine",
            "scf",
            "tight",
            "empiricaldispersion",
            "gd3bj",
        ):
            assert modifier in tokens2

        route3 = "! PBE0 def2-TZVP TightSCF Grid5 D3BJ Opt"
        drop3 = ["pbe0", "def2tzvp", "opt"]
        tokens3 = canonicalize_route_string(route3, drop_terms=drop3)
        for modifier in ("tightscf", "grid5", "d3bj"):
            assert modifier in tokens3
        for known in ("pbe0", "def2tzvp", "opt"):
            assert known not in tokens3

    def test_canonicalize_route_string_returns_none_for_falsy_input(self):
        assert canonicalize_route_string("") is None
        assert canonicalize_route_string(None) is None

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
        assert human_size(1024**4) == "1.0 TB"

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
        assert standardize_basis_set(None) is None
        assert standardize_basis_set("") == ""

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

    def test_sort_structure_dicts_by_energy_empty_input(self, tmp_path):
        db = Database(str(tmp_path / "empty.db"))
        db.create()
        assert sort_structure_dicts_by_energy(db.db_file, []) == []

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

        with pytest.raises(ValueError, match="No record found at index"):
            resolve_record(db, record_index=999)

    def test_resolve_record_by_id_not_found_raises(self):
        # get_record_by_partial_id resolves to a full ID, but get_record
        # then finds nothing for it (e.g. deleted between the two calls).
        mock_db = MagicMock()
        mock_db.get_record_by_partial_id.return_value = "full-record-id"
        mock_db.get_record.return_value = None

        with pytest.raises(ValueError, match="No record found with ID"):
            resolve_record(mock_db, record_id="partial")

    def test_resolve_record_requires_index_or_id(self):
        mock_db = MagicMock()
        with pytest.raises(
            ValueError, match="Either record_index or record_id"
        ):
            resolve_record(mock_db)


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

    def test_get_record_without_index_or_id_returns_none(self, tmp_path):
        db = Database(str(tmp_path / "empty2.db"))
        db.create()
        assert db.get_record() is None

    def test_get_record_missing_id_returns_none(self, tmp_path):
        db = Database(str(tmp_path / "empty3.db"))
        db.create()
        assert db.get_record(record_id="no-such-record") is None

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
        assert {r["program"] for r in co2_records} == {"Gaussian", "ORCA"}

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

    def test_schema_version_matches_does_not_raise(self, tmp_path):
        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        check_schema_version(db.db_file)  # should not raise

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

    def test_upsert_fills_missing_perception_columns(
        self, tmp_path, gaussian_co2_opt_outfile
    ):
        """molecules UPSERT: a later record supplying SMILES / InChI fills in
        the NULLs left by an earlier record with the same molecule_id."""
        record = SingleFileAssembler(gaussian_co2_opt_outfile).assemble_data
        molecule_id = record.molecules[0]["molecule_id"]
        original_smiles = record.molecules[0]["smiles"]
        original_inchi = record.molecules[0]["inchi"]
        # sanity: parser produces both perception fields
        assert original_smiles is not None
        assert original_inchi is not None

        # First insertion: blank out perception fields on every frame so the
        # row lands with NULL across the board
        for m in record.molecules:
            m["smiles"] = None
            m["inchi"] = None
        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        db.insert_records([record])
        row = db.get_molecule(molecule_id)
        assert row["smiles"] is None
        assert row["inchi"] is None

        # Second insertion: restore values; UPSERT must fill both NULLs
        for m in record.molecules:
            m["smiles"] = original_smiles
            m["inchi"] = original_inchi
        db.insert_records([record])
        row = db.get_molecule(molecule_id)
        assert row["smiles"] == original_smiles
        assert row["inchi"] == original_inchi

    def test_upsert_does_not_overwrite_perception_columns(
        self, tmp_path, gaussian_co2_opt_outfile
    ):
        """molecules UPSERT: existing non-NULL SMILES / InChI are preserved
        when a later record arrives with these fields NULL."""
        record = SingleFileAssembler(gaussian_co2_opt_outfile).assemble_data
        molecule_id = record.molecules[0]["molecule_id"]
        original_smiles = record.molecules[0]["smiles"]
        original_inchi = record.molecules[0]["inchi"]
        assert original_smiles is not None
        assert original_inchi is not None

        db = Database(str(tmp_path / "chemsmart.db"))
        db.create()
        db.insert_records([record])

        # Second insertion: same molecule_id but no perception fields
        for m in record.molecules:
            m["smiles"] = None
            m["inchi"] = None
        db.insert_records([record])
        row = db.get_molecule(molecule_id)
        assert row["smiles"] == original_smiles
        assert row["inchi"] == original_inchi


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
        assert {s["program"] for s in summaries} == {"Gaussian", "ORCA"}
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
            db.get_record_by_partial_id(water_gaussian.record_id[:12])
            == water_gaussian.record_id
        )
        assert db.get_molecule_by_partial_id(INCHIKEY_CO2[:16]) == INCHIKEY_CO2
        assert (
            db.get_structure_by_partial_id(STRUCTURE_ID_ORIGIN_HE[:12])
            == STRUCTURE_ID_ORIGIN_HE
        )

    def test_partial_id_resolution_error_branches(self, tmp_path):
        """No-match and ambiguous-match errors for
        get_structure_by_partial_id, exercised via raw rows rather
        than real assembled data (real content-hash structure_ids
        don't collide on a short prefix in a small test database)."""
        db = Database(str(tmp_path / "partial_errors.db"))
        db.create()
        conn = db.get_connection()
        conn.execute("INSERT INTO molecules (molecule_id) VALUES (?)", ("m1",))
        conn.execute(
            "INSERT INTO structures (structure_id, molecule_id) VALUES (?, ?)",
            ("abc111", "m1"),
        )
        conn.execute(
            "INSERT INTO structures (structure_id, molecule_id) VALUES (?, ?)",
            ("abc222", "m1"),
        )
        conn.commit()
        conn.close()

        with pytest.raises(ValueError, match="No structure found"):
            db.get_structure_by_partial_id("xyz")

        with pytest.raises(ValueError, match="Ambiguous ID prefix"):
            db.get_structure_by_partial_id("abc")

    @pytest.mark.parametrize(
        ("lookup", "prefix", "message"),
        (
            ("get_record_by_partial_id", "zzz", "No record found"),
            ("get_record_by_partial_id", "rec", "matches 2 records"),
            ("get_molecule_by_partial_id", "zzz", "No molecule found"),
            ("get_molecule_by_partial_id", "mol", "matches 2 molecules"),
        ),
    )
    def test_record_and_molecule_partial_id_errors(
        self, raw_rows_db, lookup, prefix, message
    ):
        with pytest.raises(ValueError, match=message):
            getattr(raw_rows_db, lookup)(prefix)

    def test_record_summary_lookups_without_structures(self, raw_rows_db):
        by_index = raw_rows_db.get_record_summary(record_index=1)
        assert by_index["record_id"] == "rec111"
        assert by_index["program"] == "Gaussian"
        # No linked structure, so molecule fields are not merged in.
        assert "chemical_formula" not in by_index

        assert raw_rows_db.get_record_summary() is None
        assert raw_rows_db.get_record_summary(record_index=99) is None
        assert raw_rows_db.get_record_summary(record_id="missing") is None

        summaries = raw_rows_db.get_all_record_summaries()
        assert [s["record_id"] for s in summaries] == ["rec111", "rec222"]
        assert all("chemical_formula" not in s for s in summaries)

    def test_row_helpers_handle_optional_fields(self):
        meta = Database._extract_meta(
            {"custom_basis_json": '{"Pd": "def2-TZVP"}'}
        )
        assert meta["custom_basis"] == {"Pd": "def2-TZVP"}
        assert meta["solvent_on"] is False

        provenance = Database._extract_provenance({})
        assert "normal_termination" not in provenance

        structure = Database._record_structure_row_to_dict(
            {"mulliken_spin_densities_json": "[0.5, -0.5]"}
        )
        assert structure["mulliken_spin_densities"] == [0.5, -0.5]
        assert "mulliken_atomic_charges" not in structure

    def test_pick_primary_forces_method_basis_without_forces(
        self, raw_rows_db
    ):
        assert raw_rows_db.pick_primary_forces_method_basis([]) == (
            None,
            None,
        )
        assert raw_rows_db.pick_primary_forces_method_basis(["missing"]) == (
            None,
            None,
        )


class TestDatabaseForcesQueries:
    """pick_primary_forces_method_basis, resolve_method_basis,
    get_forces_for_structure_at, and get_forces_for_record_structure_at
    had no test coverage at all. Populated via raw SQL rows rather than
    a full SingleFileAssembler pipeline, since these queries only care
    about the records/structures/record_structures columns directly."""

    def _seed_forces_db(self, tmp_path):
        db = Database(str(tmp_path / "forces.db"))
        db.create()
        conn = db.get_connection()
        conn.execute("INSERT INTO molecules (molecule_id) VALUES (?)", ("m1",))
        conn.execute(
            "INSERT INTO records (record_id, method, basis, record_index) "
            "VALUES (?, ?, ?, ?)",
            ("rec1", "b3lyp", "def2svp", 1),
        )
        conn.execute(
            "INSERT INTO structures (structure_id, molecule_id) VALUES (?, ?)",
            ("struct1", "m1"),
        )
        conn.execute(
            "INSERT INTO record_structures "
            "(record_id, structure_id, index_in_record, energy, forces_json) "
            "VALUES (?, ?, ?, ?, ?)",
            ("rec1", "struct1", 0, -100.0, json.dumps([[0.1, 0.2, 0.3]])),
        )
        conn.commit()
        conn.close()
        return db

    def test_pick_primary_forces_method_basis_returns_best_match(
        self, tmp_path
    ):
        db = self._seed_forces_db(tmp_path)
        assert db.pick_primary_forces_method_basis(["struct1"]) == (
            "b3lyp",
            "def2svp",
        )

    def test_pick_primary_forces_method_basis_empty_input(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.pick_primary_forces_method_basis([]) == (None, None)

    def test_resolve_method_basis_case_insensitive_match(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.resolve_method_basis("B3LYP", "DEF2SVP") == (
            "b3lyp",
            "def2svp",
        )

    def test_resolve_method_basis_none_input_returns_none(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.resolve_method_basis(None, "def2svp") is None

    def test_resolve_method_basis_no_match_returns_none(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.resolve_method_basis("nomatch", "nomatch") is None

    def test_get_forces_for_structure_at_match(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        forces, energy = db.get_forces_for_structure_at(
            "struct1", "b3lyp", "def2svp"
        )
        assert forces == [[0.1, 0.2, 0.3]]
        assert energy == -100.0

    def test_get_forces_for_structure_at_no_match(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.get_forces_for_structure_at(
            "nomatch", "b3lyp", "def2svp"
        ) == (None, None)

    def test_get_forces_for_record_structure_at_match(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        forces, energy = db.get_forces_for_record_structure_at(
            "rec1", "struct1", "b3lyp", "def2svp"
        )
        assert forces == [[0.1, 0.2, 0.3]]
        assert energy == -100.0

    def test_get_forces_for_record_structure_at_no_match(self, tmp_path):
        db = self._seed_forces_db(tmp_path)
        assert db.get_forces_for_record_structure_at(
            "rec1", "nomatch", "b3lyp", "def2svp"
        ) == (None, None)


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

        # '==' should parse as an equality operator as well
        eq_alias = DatabaseQuery(db.db_file, "program == 'orca'")
        where_eq, params_eq = eq_alias.parse_query()
        assert "r.program == ?" in where_eq
        assert params_eq == ("orca",)

        # Test parsing LIKE pattern
        like = DatabaseQuery(db.db_file, 'source ~ "CO2"')
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
        matched = DatabaseQuery(db.db_file, "program = 'ORCA'")
        assert matched.count_matched() == 1
        assert matched.query_summaries()[0]["program"] == "ORCA"

        matched_eq_alias = DatabaseQuery(db.db_file, "program == 'ORCA'")
        assert matched_eq_alias.count_matched() == matched.count_matched()

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

        # Unquoted numeric-looking values parse as float when not a
        # plain int, and fall back to the raw string when neither.
        float_query = DatabaseQuery(db.db_file, "total_energy = -3.5")
        _, float_params = float_query.parse_query()
        assert float_params == (-3.5,)

        bareword_query = DatabaseQuery(db.db_file, "program = ORCA")
        _, bareword_params = bareword_query.parse_query()
        assert bareword_params == ("ORCA",)

        # format_summary with a limit but no query string covers the
        # "no query" skip alongside the "limit set" branch together.
        limited_no_query = DatabaseQuery(db.db_file, None, limit=1)
        formatted_limited = limited_no_query.format_summary(
            limited_no_query.query_summaries()
        )
        assert "Query   :" not in formatted_limited
        assert "Limit   : 1" in formatted_limited

    def test_format_table_empty_and_missing_values(self, tmp_path):
        db = Database(str(tmp_path / "empty.db"))
        db.create()
        dq = DatabaseQuery(db.db_file, None)

        assert dq._format_table([]) == []

        rows = dq._format_table([{"record_index": 1}])
        assert any("" in row for row in rows[2:])

    def test_parse_query_empty_clause_parts_is_unreachable_defensively(
        self, tmp_path, monkeypatch
    ):
        """See BUGS_FOUND.md #60: `if not clause_parts:` can never be
        true through any real query string, since every token either
        appends to clause_parts or raises. Force an empty token list
        directly to exercise the guard for coverage."""
        import types

        import chemsmart.database.query as query_module

        db = Database(str(tmp_path / "empty.db"))
        db.create()
        dq = DatabaseQuery(db.db_file, "program = 'ORCA'")

        fake_splitter = types.SimpleNamespace(split=lambda s: [])
        monkeypatch.setattr(query_module, "_LOGIC_SPLIT_RE", fake_splitter)
        with pytest.raises(ValueError, match="Empty query string"):
            dq.parse_query()

    def test_parse_query_unsupported_operator_is_unreachable_defensively(
        self, tmp_path, monkeypatch
    ):
        """See BUGS_FOUND.md #60: query_condition_pattern's operator
        group only ever matches operators already in
        SUPPORTED_OPERATORS, so `if operator not in SUPPORTED_OPERATORS`
        can never be true through any real query string. Force a fake
        match object to exercise the guard for coverage."""
        import types

        import chemsmart.database.query as query_module

        db = Database(str(tmp_path / "empty.db"))
        db.create()
        dq = DatabaseQuery(db.db_file, "total_energy !! -3")

        fake_match = types.SimpleNamespace(
            groups=lambda: ("total_energy", "!!", "-3")
        )
        fake_condition_re = types.SimpleNamespace(
            fullmatch=lambda token: fake_match
        )
        monkeypatch.setattr(query_module, "_CONDITION_RE", fake_condition_re)
        with pytest.raises(ValueError, match="Unsupported operator"):
            dq.parse_query()


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
        assert programs == {"Gaussian", "ORCA"}

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
        assert csv_programs == {"Gaussian", "ORCA"}
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
        assert ("Gaussian", 1) in stats["programs"]
        assert ("ORCA", 1) in stats["programs"]
        assert "Database Overview" in inspector.format_overview()

        record_inspector = DatabaseInspector(
            db.db_file, record_id=water_gaussian.record_id[:12]
        )
        assert record_inspector.resolve_id() == water_gaussian.record_id
        detail = record_inspector.record_detail()
        assert detail["record_id"] == water_gaussian.record_id
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
        assert record["record_id"] == water_orca.record_id
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
        assert water_gaussian.provenance["program"] == "Gaussian"
        assert water_gaussian.provenance["parser"] == "Gaussian16Output"
        assert water_gaussian.meta["method"] == "mp2"
        assert water_gaussian.meta["basis"] == "aug-cc-pvtz"
        assert water_gaussian.meta["jobtype"] == "opt"
        assert water_gaussian.results["total_energy"] == -76.328992324258
        assert water_gaussian.molecules[-1]["chemical_formula"] == "H2O"

        water_orca = SingleFileAssembler(water_output_gas_path).assemble_data
        assert isinstance(water_orca, AssembledRecord)
        assert water_orca.provenance["program"] == "ORCA"
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
        with pytest.raises(ValueError, match="Unsupported file"):
            SingleFileAssembler(unknown_output)._get_assembler(unknown_output)


class TestStaticDatabaseContent:
    def test_database_counts(self, database_chemsmart_file):
        """The fixture holds exactly 47 records, 33 molecules, 314 structures."""
        db = Database(database_chemsmart_file)
        assert db.count_records() == 47
        assert db.count_molecules() == 33
        assert db.count_structures() == 314

    def test_program_distribution(self, database_chemsmart_file):
        """36 Gaussian records and 11 ORCA records are present."""
        db = Database(database_chemsmart_file)
        programs = {}
        for s in db.get_all_record_summaries():
            programs[s["program"]] = programs.get(s["program"], 0) + 1
        assert programs["Gaussian"] == 36
        assert programs["ORCA"] == 11

    def test_gaussian_mp2_water_record(
        self, database_chemsmart_file, gaussian_mp2_outputfile
    ):
        """Gaussian MP2/aug-cc-pVTZ opt of water (water_mp2.log)."""
        db = Database(database_chemsmart_file)
        record = db.get_record(record_id=RECORD_ID_GAUSSIAN_MP2_WATER)
        assembled = SingleFileAssembler(gaussian_mp2_outputfile).assemble_data

        assert (
            record["provenance"]["program"]
            == assembled.provenance["program"]
            == "Gaussian"
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

        assert (
            record["provenance"]["program"]
            == assembled.provenance["program"]
            == "ORCA"
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
        assert programs == {"Gaussian", "ORCA"}

    def test_he_molecule(self, database_chemsmart_file):
        """He molecule: monoatomic, single geometry, 3 records across Gaussian and ORCA."""
        db = Database(database_chemsmart_file)
        mol = db.get_molecule(INCHIKEY_HE)
        assert mol is not None
        assert mol["chemical_formula"] == "He"
        assert mol["element_counts"] == {"He": 1}
        assert mol["number_of_atoms"] == 1
        assert mol["is_monoatomic"] is True
        records = db.get_records_for_molecule(INCHIKEY_HE)
        assert len(records) == 3
        programs = {r["program"] for r in records}
        assert programs == {"Gaussian", "ORCA"}

    def test_query_by_program(self, database_chemsmart_file):
        """Filtering by program='gaussian' returns exactly 36 of the 47 records."""
        dq = DatabaseQuery(database_chemsmart_file, "program = 'Gaussian'")
        assert dq.count_total() == 47
        assert dq.count_matched() == 36

    def test_query_by_method_mp2(self, database_chemsmart_file):
        """Only one record uses plain MP2 — Gaussian water_mp2.log."""
        dq = DatabaseQuery(database_chemsmart_file, "method = 'mp2'")
        assert dq.count_matched() == 1
        result = dq.query_summaries()[0]
        assert result["program"] == "Gaussian"
        assert result["basis"] == "aug-cc-pvtz"
        assert result["chemical_formula"] == "H2O"

    def test_query_limit(self, database_chemsmart_file):
        """limit= truncates results but count_total() still reflects the full DB."""
        dq = DatabaseQuery(database_chemsmart_file, None, limit=5)
        assert len(dq.query_summaries()) == 5
        assert dq.count_total() == 47

    def test_inspect_overview(self, database_chemsmart_file):
        """Inspector overview matches the manually verified DB-wide statistics."""
        inspector = DatabaseInspector(database_chemsmart_file)
        stats = inspector.overview()
        assert stats["num_records"] == 47
        assert stats["num_molecules"] == 33
        assert stats["num_structures"] == 314
        assert ("Gaussian", 36) in stats["programs"]
        assert ("ORCA", 11) in stats["programs"]

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
        in this 47-record database."""
        db = Database(database_chemsmart_file)
        resolved = db.get_record_by_partial_id(
            RECORD_ID_GAUSSIAN_MP2_WATER[:12]
        )
        assert resolved == RECORD_ID_GAUSSIAN_MP2_WATER

    def test_xtb_database_counts(self, database_chemsmart_xtb_file):
        """The xTB fixture holds exactly 7 records, 6 molecules, 41 structures."""
        db = Database(database_chemsmart_xtb_file)
        assert db.count_records() == 7
        assert db.count_molecules() == 6
        assert db.count_structures() == 41

    def test_xtb_program_distribution(self, database_chemsmart_xtb_file):
        """All 7 records are xTB."""
        db = Database(database_chemsmart_xtb_file)
        programs = {}
        for s in db.get_all_record_summaries():
            programs[s["program"]] = programs.get(s["program"], 0) + 1
        assert programs["xTB"] == 7

    def test_xtb_gfn2_water_record(
        self, database_chemsmart_xtb_file, xtb_water_outfolder
    ):
        """xTB GFN2 ohess of water (water_ohess/)."""
        db = Database(database_chemsmart_xtb_file)
        record = db.get_record(record_id=RECORD_ID_XTB_GFN2_WATER)
        assembled = SingleFolderAssembler(xtb_water_outfolder).assemble_data

        assert (
            record["provenance"]["program"]
            == assembled.provenance["program"]
            == "xTB"
        )
        assert record["meta"]["method"] == assembled.meta["method"] == "gfn2"
        assert record["meta"]["basis"] == assembled.meta["basis"] == "default"
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
            record["results"]["total_energy"], -5.0705444, rtol=1e-6
        )
        assert np.isclose(
            record["results"]["homo_energy"], assembled.results["homo_energy"]
        )
        assert np.isclose(
            record["results"]["homo_energy"], -12.1467, rtol=1e-6
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

    def test_xtb_h2o_molecule(self, database_chemsmart_xtb_file):
        """H2O molecule in the xTB DB links to the single water_ohess record."""
        db = Database(database_chemsmart_xtb_file)
        mol = db.get_molecule(INCHIKEY_H2O)
        assert mol is not None
        assert mol["chemical_formula"] == "H2O"
        assert mol["element_counts"] == {"O": 1, "H": 2}
        assert mol["is_linear"] is False
        assert mol["is_monoatomic"] is False
        records = db.get_records_for_molecule(INCHIKEY_H2O)
        assert len(records) == 1
        assert records[0]["program"] == "xTB"
        assert records[0]["record_id"] == RECORD_ID_XTB_GFN2_WATER

    def test_xtb_he_molecule(self, database_chemsmart_xtb_file):
        """He molecule: monoatomic, single hess record in the xTB DB."""
        db = Database(database_chemsmart_xtb_file)
        mol = db.get_molecule(INCHIKEY_HE)
        assert mol is not None
        assert mol["chemical_formula"] == "He"
        assert mol["element_counts"] == {"He": 1}
        assert mol["number_of_atoms"] == 1
        assert mol["is_monoatomic"] is True
        records = db.get_records_for_molecule(INCHIKEY_HE)
        assert len(records) == 1
        assert records[0]["record_id"] == RECORD_ID_XTB_HE
        assert records[0]["jobtype"] == "hess"

    def test_xtb_c6h4_molecule(self, database_chemsmart_xtb_file):
        """C6H4 (p-benzyne) links to both ALPB opt and sp records."""
        db = Database(database_chemsmart_xtb_file)
        mol = db.get_molecule(INCHIKEY_C6H4)
        assert mol is not None
        assert mol["chemical_formula"] == "C6H4"
        assert mol["element_counts"] == {"C": 6, "H": 4}
        records = db.get_records_for_molecule(INCHIKEY_C6H4)
        assert len(records) == 2
        assert {r["jobtype"] for r in records} == {"opt", "sp"}
        assert {r["program"] for r in records} == {"xTB"}

    def test_xtb_query_by_program(self, database_chemsmart_xtb_file):
        """Filtering by program='xTB' returns all 7 records."""
        dq = DatabaseQuery(database_chemsmart_xtb_file, "program = 'xTB'")
        assert dq.count_total() == 7
        assert dq.count_matched() == 7

    def test_xtb_query_by_jobtype_opt(self, database_chemsmart_xtb_file):
        """One of the seven xTB records are single point calculations."""
        dq = DatabaseQuery(database_chemsmart_xtb_file, "jobtype = 'sp'")
        assert dq.count_matched() == 1
        result = dq.query_summaries()[0]
        assert result["chemical_formula"] == "C6H4"
        assert np.isclose(result["total_energy"], -14.66185696)

    def test_xtb_query_limit(self, database_chemsmart_xtb_file):
        """limit= truncates results but count_total() still reflects the full DB."""
        dq = DatabaseQuery(database_chemsmart_xtb_file, None, limit=3)
        assert len(dq.query_summaries()) == 3
        assert dq.count_total() == 7

    def test_xtb_inspect_overview(self, database_chemsmart_xtb_file):
        """Inspector overview matches the xTB fixture statistics."""
        inspector = DatabaseInspector(database_chemsmart_xtb_file)
        stats = inspector.overview()
        assert stats["num_records"] == 7
        assert stats["num_molecules"] == 6
        assert stats["num_structures"] == 41
        assert ("xTB", 7) in stats["programs"]
        assert ("gfn2", 7) in stats["methods"]
        assert ("default", 7) in stats["basis_sets"]
        assert stats["solvated_records"] == 2
        assert stats["gas_phase_records"] == 5

    def test_xtb_inspect_water_record(self, database_chemsmart_xtb_file):
        """Inspector record detail for xTB water matches the stored ground truth."""
        inspector = DatabaseInspector(
            database_chemsmart_xtb_file,
            record_id=RECORD_ID_XTB_GFN2_WATER[:12],
        )
        assert inspector.resolve_id() == RECORD_ID_XTB_GFN2_WATER
        detail = inspector.record_detail()
        assert detail["meta"]["method"] == "gfn2"
        assert np.isclose(
            detail["results"]["total_energy"], -5.07054444, rtol=1e-6
        )

    def test_xtb_inspect_molecule_h2o(self, database_chemsmart_xtb_file):
        """Molecule inspector for H2O reports the single water_ohess record."""
        inspector = DatabaseInspector(
            database_chemsmart_xtb_file,
            molecule_id=INCHIKEY_H2O,
        )
        mol, structures, records = inspector.molecule_detail()
        assert mol["chemical_formula"] == "H2O"
        assert len(records) == 1

    def test_xtb_partial_id_for_water(self, database_chemsmart_xtb_file):
        """A 12-character prefix uniquely identifies xTB water in this DB."""
        db = Database(database_chemsmart_xtb_file)
        resolved = db.get_record_by_partial_id(RECORD_ID_XTB_GFN2_WATER[:12])
        assert resolved == RECORD_ID_XTB_GFN2_WATER


class TestAssembleFailedTermination:
    def test_failed_termination_without_include_failed_returns_none(self):
        output = MagicMock()
        output.normal_termination = False
        assembler = _make_base_assembler(output, include_failed=False)
        assert assembler.assemble() is None

    def test_failed_termination_with_include_failed_still_assembles(self):
        """include_failed=True logs a different warning but still
        proceeds to assemble partial data instead of bailing out."""
        output = MagicMock()
        output.normal_termination = False
        output.route_string = "opt"
        output.method = "b3lyp"
        output.basis = "def2svp"
        output.jobtype = "opt"
        output.solvent_model = None
        output.solvent_id = None
        output.solvent_on = False
        output.freq = False
        output.num_basis_functions = 10
        output.spin = "singlet"
        molecule = SimpleNamespace(structure_id="mol-1")
        assembler = _make_base_assembler(output, include_failed=True)
        with (
            patch.object(
                BaseAssembler,
                "molecules_list",
                new_callable=property,
                fget=lambda self: [molecule],
            ),
            patch.object(
                BaseAssembler,
                "get_calculation_results",
                return_value={},
            ),
            patch.object(
                BaseAssembler,
                "get_molecule_info",
                return_value={},
            ),
            patch.object(
                BaseAssembler,
                "build_provenance",
                return_value={"program": "Gaussian"},
            ),
        ):
            record = assembler.assemble()
        assert record is not None


class TestAssembleEmptyMolecules:
    def test_empty_molecules_list_returns_none(self):
        output = MagicMock()
        output.normal_termination = True
        assembler = _make_base_assembler(output)
        with patch.object(
            BaseAssembler,
            "molecules_list",
            new_callable=property,
            fget=lambda self: [],
        ):
            assert assembler.assemble() is None


class TestGetMetaDataCustomBasisAndSolvent:
    def _make_output_for_meta(self, **overrides):
        output = MagicMock()
        output.basis = "6-31g(d)"
        output.method = "b3lyp"
        output.num_basis_functions = 50
        output.spin = "singlet"
        output.jobtype = "opt"
        output.solvent_on = False
        output.solvent_model = None
        output.solvent_id = None
        output.custom_solvent = None
        output.route_string = "opt b3lyp/6-31g(d)"
        output.freq = False
        for key, value in overrides.items():
            setattr(output, key, value)
        return output

    def test_custom_basis_uses_customized_basis_label(self):
        output = self._make_output_for_meta(basis="gen")
        assembler = _make_base_assembler(output)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=True
        ):
            meta = assembler.get_meta_data()
        assert meta["basis"] == "customized_basis"

    def test_standard_basis_is_standardized_not_customized(self):
        output = self._make_output_for_meta(basis="6-31g(d)")
        assembler = _make_base_assembler(output)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert meta["basis"] != "customized_basis"

    def test_custom_solvent_uses_customized_solvent_label(self):
        output = self._make_output_for_meta(
            basis="6-31g(d)",
            solvent_on=True,
            solvent_model="smd",
            solvent_id="mysolvent",
            custom_solvent="eps=5.0",
        )
        assembler = _make_base_assembler(output)
        with (
            patch(
                "chemsmart.database.assemble.is_custom_basis",
                return_value=False,
            ),
            patch(
                "chemsmart.database.assemble.is_custom_solvent",
                return_value=True,
            ),
        ):
            meta = assembler.get_meta_data()
        assert meta["solvent_id"] == "customized_solvent"

    def test_standard_solvent_keeps_solvent_id(self):
        output = self._make_output_for_meta(
            basis="6-31g(d)",
            solvent_on=True,
            solvent_model="smd",
            solvent_id="water",
            custom_solvent=None,
        )
        assembler = _make_base_assembler(output)
        with (
            patch(
                "chemsmart.database.assemble.is_custom_basis",
                return_value=False,
            ),
            patch(
                "chemsmart.database.assemble.is_custom_solvent",
                return_value=False,
            ),
        ):
            meta = assembler.get_meta_data()
        assert meta["solvent_id"] == "water"


class TestGetMoleculeInfoMullikenSpinDensities:
    def test_non_none_spin_densities_are_included(self):
        assembler = _make_base_assembler(MagicMock())
        mol = MagicMock()
        mol.mulliken_atomic_charges = None
        mol.mulliken_spin_densities = [0.1, -0.1]
        mol.rotational_symmetry_number = None
        mol.rotational_constants = None
        mol.point_group = None
        mol.vibrational_frequencies = None
        mol.dipole_moment = None
        mol.dipole_moment_magnitude = None

        info = assembler.get_molecule_info(mol)
        assert info["mulliken_spin_densities"] == [0.1, -0.1]


class TestGaussianAssemblerMetaDataBranches:
    def _make_gaussian_assembler(self, **overrides):
        assembler = GaussianAssembler.__new__(GaussianAssembler)
        output = MagicMock()
        output.basis = "6-31g(d)"
        output.method = "b3lyp"
        output.num_basis_functions = 50
        output.spin = "singlet"
        output.jobtype = "opt"
        output.solvent_on = False
        output.solvent_model = None
        output.solvent_id = None
        output.custom_solvent = None
        output.route_string = "opt b3lyp/6-31g(d)"
        output.freq = False
        output.num_primitive_gaussians = 100
        output.num_cartesian_basis_functions = 60
        output.modredundant_group = None
        for key, value in overrides.items():
            setattr(output, key, value)
        assembler.output = output
        assembler.filename = "some_file.log"
        assembler.target = assembler.filename
        assembler.include_failed = False
        return assembler

    def test_modredundant_group_included_when_present(self):
        assembler = self._make_gaussian_assembler(modredundant_group="1 2 3 F")
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert meta["modredundant_group"] == "1 2 3 F"

    def test_modredundant_group_absent_when_none(self):
        assembler = self._make_gaussian_assembler(modredundant_group=None)
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=False
        ):
            meta = assembler.get_meta_data()
        assert "modredundant_group" not in meta

    def test_custom_basis_includes_heavy_and_light_element_details(self):
        assembler = self._make_gaussian_assembler(
            basis="gen",
            heavy_elements="Fe",
            heavy_elements_basis="lanl2dz",
            heavy_elements_ecp="lanl2dz",
            light_elements="C,H,O",
            light_elements_basis="6-31g(d)",
        )
        with patch(
            "chemsmart.database.assemble.is_custom_basis", return_value=True
        ):
            meta = assembler.get_meta_data()
        assert meta["custom_basis"]["heavy_elements"] == "Fe"
        assert meta["custom_basis"]["light_elements"] == "C,H,O"


class TestSingleFileAssemblerErrorHandling:
    def test_assemble_data_returns_none_and_logs_on_exception(self):
        single = SingleFileAssembler(filename="broken.log")
        broken_assembler = MagicMock()
        broken_assembler.assemble.side_effect = RuntimeError("parse failed")
        with patch.object(
            SingleFileAssembler,
            "_get_assembler",
            return_value=broken_assembler,
        ):
            assert single.assemble_data is None

    def test_unsupported_file_type_raises(self):
        single = SingleFileAssembler(filename="mystery.txt")
        with patch(
            "chemsmart.database.assemble.get_program_type_from_file",
            return_value="unknown",
        ):
            with pytest.raises(ValueError, match="Unsupported file"):
                single._get_assembler("mystery.txt")


class TestSingleFolderAssemblerErrorHandling:
    def test_assemble_data_returns_none_and_logs_on_exception(self):
        single = SingleFolderAssembler(folder="broken_folder")
        broken_assembler = MagicMock()
        broken_assembler.assemble.side_effect = RuntimeError("parse failed")
        with patch.object(
            SingleFolderAssembler,
            "_get_assembler",
            return_value=broken_assembler,
        ):
            assert single.assemble_data is None

    def test_mixed_program_folder_raises(self):
        single = SingleFolderAssembler(folder="mixed_folder")
        with patch(
            "chemsmart.database.assemble.BaseFolder"
        ) as mock_base_folder_cls:
            mock_base_folder_cls.return_value.get_program_type_from_folder.return_value = (
                "mixed"
            )
            with pytest.raises(ValueError, match="multiple programs"):
                single._get_assembler("mixed_folder")

    def test_unsupported_folder_type_raises(self):
        single = SingleFolderAssembler(folder="gaussian_folder")
        with patch(
            "chemsmart.database.assemble.BaseFolder"
        ) as mock_base_folder_cls:
            mock_base_folder_cls.return_value.get_program_type_from_folder.return_value = (
                "gaussian"
            )
            with pytest.raises(ValueError, match="Unsupported folder"):
                single._get_assembler("gaussian_folder")


class TestValidateForces:
    def test_none_forces_invalid(self):
        exporter = _make_exporter()
        assert exporter._validate_forces(FakeMolecule(forces=None)) is False

    def test_correct_shape_valid(self):
        exporter = _make_exporter()
        forces = [[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]]
        assert (
            exporter._validate_forces(FakeMolecule(forces=forces, num_atoms=2))
            is True
        )

    def test_wrong_shape_invalid(self):
        exporter = _make_exporter()
        forces = [[0.0, 0.0, 0.0]]  # only 1 row, but num_atoms=2
        assert (
            exporter._validate_forces(FakeMolecule(forces=forces, num_atoms=2))
            is False
        )

    def test_unconvertible_forces_invalid(self):
        exporter = _make_exporter()
        forces = "not-an-array"
        assert exporter._validate_forces(FakeMolecule(forces=forces)) is False


class TestWriteExtxyzFrame:
    def test_writes_forces_and_energy_when_present(self):
        exporter = _make_exporter(output="out.extxyz")
        mol = FakeMolecule(
            forces=[[0.0, 0.0, 0.0], [0.1, 0.2, 0.3]], num_atoms=2
        )
        frame = {
            "molecule": mol,
            "structure_id": "abc123def456",
            "energies": [("b3lyp", "def2-svp", -1.234567890123)],
            "primary_energy": -1.234567890123,
            "primary_method": "b3lyp",
            "primary_basis": "def2-svp",
        }
        import io

        buf = io.StringIO()
        exporter._write_extxyz_frame(buf, frame)
        content = buf.getvalue()
        assert "forces:R:3" in content
        assert 'energy_units="Hartree"' in content
        assert "SID: abc123def456" in content
        assert "formula=H2" in content
        assert "E(b3lyp/def2-svp)=" in content

    def test_writes_without_forces_or_energy(self):
        exporter = _make_exporter(output="out.extxyz")
        mol = FakeMolecule(forces=None, chemical_formula=None)
        frame = {
            "molecule": mol,
            "structure_id": None,
            "energies": [],
            "primary_energy": None,
            "primary_method": None,
            "primary_basis": None,
        }
        import io

        buf = io.StringIO()
        exporter._write_extxyz_frame(buf, frame)
        content = buf.getvalue()
        assert "forces:R:3" not in content
        assert "energy=" not in content
        assert "SID:" not in content

    def test_multiple_energies_use_indexed_tags(self):
        exporter = _make_exporter(output="out.extxyz")
        mol = FakeMolecule()
        frame = {
            "molecule": mol,
            "structure_id": "sid1",
            "energies": [
                ("b3lyp", "def2-svp", -1.0),
                ("wb97xd", "def2-tzvp", -1.1),
            ],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "def2-svp",
        }
        import io

        buf = io.StringIO()
        exporter._write_extxyz_frame(buf, frame)
        content = buf.getvalue()
        assert "E(b3lyp/def2-svp)=" in content
        assert "E1(wb97xd/def2-tzvp)=" in content

    def test_malformed_forces_shape_falls_back_to_no_forces(self):
        exporter = _make_exporter(output="out.extxyz")
        # Forces present but wrong shape for the molecule's atom count.
        mol = FakeMolecule(forces=[[0.0, 0.0, 0.0]], num_atoms=2)
        frame = {
            "molecule": mol,
            "structure_id": "sid",
            "energies": [],
            "primary_energy": None,
            "primary_method": None,
            "primary_basis": None,
        }
        import io

        buf = io.StringIO()
        exporter._write_extxyz_frame(buf, frame)
        assert "forces:R:3" not in buf.getvalue()

    def test_unconvertible_forces_falls_back_to_no_forces(self):
        exporter = _make_exporter(output="out.extxyz")
        # Ragged/non-numeric forces raise inside np.asarray(..., dtype=float).
        mol = FakeMolecule(forces=[[0.0, "bad", 0.0], [0.1, 0.2]])
        frame = {
            "molecule": mol,
            "structure_id": "sid",
            "energies": [],
            "primary_energy": None,
            "primary_method": None,
            "primary_basis": None,
        }
        import io

        buf = io.StringIO()
        exporter._write_extxyz_frame(buf, frame)
        assert "forces:R:3" not in buf.getvalue()


class TestToExtxyzEndToEnd:
    def test_to_extxyz_writes_frames_via_collect_frames(self, tmp_path):
        exporter = _make_exporter(
            output=str(tmp_path / "out.extxyz"), record_index=1
        )
        mol = FakeMolecule(
            forces=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]], num_atoms=2
        )
        frame = {
            "molecule": mol,
            "structure_id": "sid1",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with patch.object(exporter, "_collect_frames", return_value=[frame]):
            exporter.to_extxyz()
        content = (tmp_path / "out.extxyz").read_text()
        assert "forces:R:3" in content


class TestExportDispatch:
    @pytest.mark.parametrize(
        "output,method_name",
        [
            ("out.json", "to_json"),
            ("out.csv", "to_csv"),
            ("out.xyz", "to_xyz"),
            ("out.extxyz", "to_extxyz"),
        ],
    )
    def test_export_dispatches_to_correct_handler(self, output, method_name):
        exporter = _make_exporter(output=output)
        with patch.object(exporter, method_name) as mock_handler:
            exporter.export()
        mock_handler.assert_called_once()


class TestCollectFramesSelectorValidation:
    def test_no_selector_raises(self):
        exporter = _make_exporter()
        with pytest.raises(ValueError, match="requires exactly one"):
            exporter._collect_frames()

    def test_empty_frames_raises(self):
        exporter = _make_exporter(record_index=1)
        with patch.object(exporter, "frames_from_record", return_value=[]):
            with pytest.raises(ValueError, match="No structures found"):
                exporter._collect_frames()

    def test_include_forces_without_effective_primary_raises(self):
        exporter = _make_exporter(record_index=1)
        frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid",
            "energies": [],
            "primary_energy": None,
            "primary_method": None,
            "primary_basis": None,
        }
        with patch.object(
            exporter, "frames_from_record", return_value=[frame]
        ):
            with pytest.raises(ValueError, match="energy and forces at"):
                exporter._collect_frames(include_forces=True)

    def test_include_forces_filters_frames_lacking_valid_forces(self):
        exporter = _make_exporter(record_index=1)
        good_mol = FakeMolecule(
            forces=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]], num_atoms=2
        )
        bad_mol = FakeMolecule(forces=None)
        frames = [
            {
                "molecule": good_mol,
                "structure_id": "sid-good",
                "energies": [("b3lyp", "sto-3g", -1.0)],
                "primary_energy": -1.0,
                "primary_method": "b3lyp",
                "primary_basis": "sto-3g",
            },
            {
                "molecule": bad_mol,
                "structure_id": "sid-bad",
                "energies": [],
                "primary_energy": None,
                "primary_method": None,
                "primary_basis": None,
            },
        ]
        with patch.object(exporter, "frames_from_record", return_value=frames):
            result = exporter._collect_frames(include_forces=True)
        assert len(result) == 1
        assert result[0]["structure_id"] == "sid-good"

    def test_include_forces_raises_when_all_frames_filtered_out(self):
        exporter = _make_exporter(record_index=1)
        bad_mol = FakeMolecule(forces=None)
        frame = {
            "molecule": bad_mol,
            "structure_id": "sid-bad",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with patch.object(
            exporter, "frames_from_record", return_value=[frame]
        ):
            with pytest.raises(ValueError, match="energy and forces at"):
                exporter._collect_frames(include_forces=True)

    def test_include_forces_reports_more_than_five_skipped(self):
        exporter = _make_exporter(record_index=1)
        good_mol = FakeMolecule(
            forces=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]], num_atoms=2
        )
        good_frame = {
            "molecule": good_mol,
            "structure_id": "sid-good",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        bad_frames = [
            {
                "molecule": FakeMolecule(forces=None),
                "structure_id": f"sid-bad-{i}",
                "energies": [],
                "primary_energy": None,
                "primary_method": None,
                "primary_basis": None,
            }
            for i in range(7)
        ]
        with patch.object(
            exporter,
            "frames_from_record",
            return_value=[good_frame] + bad_frames,
        ):
            result = exporter._collect_frames(include_forces=True)
        assert len(result) == 1

    def test_xyz_user_primary_mb_filters_frames_without_matching_energy(
        self,
    ):
        exporter = _make_exporter(
            record_index=1, method="b3lyp", basis="sto-3g"
        )
        matching_frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid1",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        non_matching_frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid2",
            "energies": [("wb97xd", "def2-tzvp", -2.0)],
            "primary_energy": -2.0,
            "primary_method": "wb97xd",
            "primary_basis": "def2-tzvp",
        }
        with patch.object(
            exporter,
            "frames_from_record",
            return_value=[matching_frame, non_matching_frame],
        ):
            result = exporter._collect_frames(include_forces=False)
        assert len(result) == 1
        assert result[0]["structure_id"] == "sid1"

    def test_xyz_user_primary_mb_raises_when_none_match(self):
        exporter = _make_exporter(
            record_index=1, method="b3lyp", basis="sto-3g"
        )
        frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid1",
            "energies": [("wb97xd", "def2-tzvp", -2.0)],
            "primary_energy": -2.0,
            "primary_method": "wb97xd",
            "primary_basis": "def2-tzvp",
        }
        with patch.object(
            exporter, "frames_from_record", return_value=[frame]
        ):
            with pytest.raises(
                ValueError, match="No structures found at the given"
            ):
                exporter._collect_frames(include_forces=False)

    def test_molecule_id_path_sorts_frames_by_energy(self):
        exporter = _make_exporter(molecule_id="mid123")
        frame1 = {
            "molecule": FakeMolecule(),
            "structure_id": "sid1",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        frame2 = {
            "molecule": FakeMolecule(),
            "structure_id": "sid2",
            "energies": [("b3lyp", "sto-3g", -2.0)],
            "primary_energy": -2.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with (
            patch.object(
                exporter,
                "frames_from_molecule_id",
                return_value=([frame1, frame2], ("b3lyp", "sto-3g")),
            ),
            patch(
                "chemsmart.database.export.sort_frames_by_energy",
                return_value=[frame2, frame1],
            ) as mock_sort,
        ):
            result = exporter._collect_frames(include_forces=False)
        mock_sort.assert_called_once()
        assert result == [frame2, frame1]

    def test_structure_id_path_infers_effective_primary_from_frame(self):
        exporter = _make_exporter(structure_id="sid123")
        frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid123",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with patch.object(
            exporter, "frames_from_structure_id", return_value=[frame]
        ):
            result = exporter._collect_frames(include_forces=False)
        assert result == [frame]

    def test_structure_id_path_skips_inference_when_user_primary_set(self):
        # With method/basis already user-specified, effective_primary is
        # non-None before the structure_id branch runs, so the frame-based
        # inference block is skipped entirely.
        exporter = _make_exporter(
            structure_id="sid123", method="b3lyp", basis="sto-3g"
        )
        frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid123",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with patch.object(
            exporter, "frames_from_structure_id", return_value=[frame]
        ):
            result = exporter._collect_frames(include_forces=False)
        assert result == [frame]

    def test_molecule_id_path_skips_inference_when_user_primary_set(self):
        exporter = _make_exporter(
            molecule_id="mid123", method="b3lyp", basis="sto-3g"
        )
        frame = {
            "molecule": FakeMolecule(),
            "structure_id": "sid1",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with (
            patch.object(
                exporter,
                "frames_from_molecule_id",
                return_value=([frame], None),
            ),
            patch(
                "chemsmart.database.export.sort_frames_by_energy",
                return_value=[frame],
            ),
        ):
            result = exporter._collect_frames(include_forces=False)
        assert result == [frame]

    def test_include_forces_no_skipped_frames_omits_warning(self, caplog):
        exporter = _make_exporter(record_index=1)
        good_mol = FakeMolecule(
            forces=[[0.0, 0.0, 0.0], [0.1, 0.1, 0.1]], num_atoms=2
        )
        frame = {
            "molecule": good_mol,
            "structure_id": "sid-good",
            "energies": [("b3lyp", "sto-3g", -1.0)],
            "primary_energy": -1.0,
            "primary_method": "b3lyp",
            "primary_basis": "sto-3g",
        }
        with patch.object(
            exporter, "frames_from_record", return_value=[frame]
        ):
            result = exporter._collect_frames(include_forces=True)
        assert len(result) == 1
        assert "Skipped" not in caplog.text


class TestToCsvWithoutKeys:
    def test_to_csv_without_parsed_keys_uses_default_columns_only(
        self, tmp_path
    ):
        exporter = _make_exporter(output=str(tmp_path / "out.csv"))
        exporter.db.get_all_records.return_value = [
            {
                "record_index": 1,
                "record_id": "r1",
                "molecules": [{"chemical_formula": "H2O"}],
            }
        ]
        exporter.to_csv()
        import csv

        with open(tmp_path / "out.csv", newline="") as f:
            rows = list(csv.DictReader(f))
        assert list(rows[0].keys()) == [
            "record_index",
            "record_id",
            "chemical_formula",
        ]


class TestParseCsvKeys:
    def test_keys_that_are_only_commas_and_whitespace_yield_empty_list(self):
        exporter = _make_exporter(keys=" , , ")
        assert exporter.parsed_keys == []


class TestFramesFromRecord:
    def test_record_index_not_found_raises(self):
        exporter = _make_exporter(record_index=5)
        exporter.db.get_record.return_value = None
        with pytest.raises(ValueError, match="No record found at index"):
            exporter.frames_from_record()

    def test_record_id_not_found_raises(self):
        exporter = _make_exporter(record_id="abc123")
        exporter.db.get_record_by_partial_id.return_value = "abc123full"
        exporter.db.get_record.return_value = None
        with pytest.raises(ValueError, match="No record found with ID"):
            exporter.frames_from_record()

    def test_record_id_path_resolves_full_id(self):
        exporter = _make_exporter(record_id="abc123")
        exporter.db.get_record_by_partial_id.return_value = "abc123full"
        exporter.db.get_record.return_value = {
            "record_id": "abc123full",
            "meta": {"method": "b3lyp", "basis": "sto-3g"},
            "molecules": [{"structure_id": "sid1", "energy": -1.0}],
        }
        with _patch_database_file(build_molecule=lambda d: FakeMolecule()):
            frames = exporter.frames_from_record()
        assert len(frames) == 1
        assert frames[0]["primary_method"] == "b3lyp"
        exporter.db.get_record.assert_called_with(record_id="abc123full")

    def test_slice_structure_index_selects_multiple_structures(self):
        exporter = _make_exporter(record_index=1, structure_index=":")
        exporter.db.get_record.return_value = {
            "record_id": "rid1",
            "meta": {"method": "b3lyp", "basis": "sto-3g"},
            "molecules": [
                {"structure_id": "sid1", "energy": -1.0},
                {"structure_id": "sid2", "energy": -2.0},
            ],
        }
        with _patch_database_file(build_molecule=lambda d: FakeMolecule()):
            frames = exporter.frames_from_record()
        assert len(frames) == 2

    def test_include_forces_queries_forces_for_record_structure(self):
        exporter = _make_exporter(record_index=1)
        exporter.db.get_record.return_value = {
            "record_id": "rid1",
            "meta": {"method": "b3lyp", "basis": "sto-3g"},
            "molecules": [{"structure_id": "sid1", "energy": -1.0}],
        }
        exporter.db.get_forces_for_record_structure_at.return_value = (
            [[0.0, 0.0, 0.0]],
            -1.5,
        )
        with _patch_database_file(build_molecule=lambda d: FakeMolecule()):
            frames = exporter.frames_from_record(include_forces=True)
        assert frames[0]["primary_energy"] == -1.5
        assert frames[0]["primary_method"] == "b3lyp"

    def test_include_forces_with_no_energy_from_forces_omits_primary(self):
        exporter = _make_exporter(record_index=1)
        exporter.db.get_record.return_value = {
            "record_id": "rid1",
            "meta": {"method": "b3lyp", "basis": "sto-3g"},
            "molecules": [{"structure_id": "sid1", "energy": -1.0}],
        }
        exporter.db.get_forces_for_record_structure_at.return_value = (
            None,
            None,
        )
        with _patch_database_file(build_molecule=lambda d: FakeMolecule()):
            frames = exporter.frames_from_record(include_forces=True)
        assert frames[0]["primary_energy"] is None
        assert frames[0]["primary_method"] is None
        assert frames[0]["energies"] == []


class TestFramesFromStructureId:
    def test_struct_not_found_raises(self):
        exporter = _make_exporter(structure_id="sid123")
        exporter.db.get_structure_by_partial_id.return_value = "sid123full"
        exporter.db.get_structure.return_value = None
        with pytest.raises(ValueError, match="No structure found with ID"):
            exporter.frames_from_structure_id()

    def test_include_forces_uses_user_primary_mb(self):
        exporter = _make_exporter(
            structure_id="sid123", method="b3lyp", basis="sto-3g"
        )
        exporter.db.get_structure_by_partial_id.return_value = "sid123full"
        exporter.db.get_structure.return_value = {"structure_id": "sid123full"}
        exporter.db.get_forces_for_structure_at.return_value = (
            [[0.0, 0.0, 0.0]],
            -1.0,
        )

        def fake_energies(db_file, sid):
            return [("b3lyp", "sto-3g", -1.0)]

        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                side_effect=fake_energies,
            ),
        ):
            frames = exporter.frames_from_structure_id(include_forces=True)
        assert len(frames) == 1
        exporter.db.pick_primary_forces_method_basis.assert_not_called()

    def test_include_forces_auto_picks_primary_mb(self):
        exporter = _make_exporter(structure_id="sid123")
        exporter.db.get_structure_by_partial_id.return_value = "sid123full"
        exporter.db.get_structure.return_value = {"structure_id": "sid123full"}
        exporter.db.pick_primary_forces_method_basis.return_value = (
            "wb97xd",
            "def2-tzvp",
        )
        exporter.db.get_forces_for_structure_at.return_value = (None, None)
        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                return_value=[],
            ),
        ):
            exporter.frames_from_structure_id(include_forces=True)
        exporter.db.pick_primary_forces_method_basis.assert_called_once_with(
            ["sid123full"]
        )


class TestFramesFromMoleculeId:
    def test_no_structures_raises(self):
        exporter = _make_exporter(molecule_id="mid123")
        exporter.db.get_molecule_by_partial_id.return_value = "mid123full"
        exporter.db.get_structures_for_molecule.return_value = []
        with pytest.raises(ValueError, match="No structures found for"):
            exporter.frames_from_molecule_id()

    def test_user_primary_mb_used_directly(self):
        exporter = _make_exporter(
            molecule_id="mid123", method="b3lyp", basis="sto-3g"
        )
        exporter.db.get_molecule_by_partial_id.return_value = "mid123full"
        exporter.db.get_structures_for_molecule.return_value = [
            {"structure_id": "sid1"}
        ]
        exporter.db.get_forces_for_structure_at.return_value = (None, None)
        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                return_value=[],
            ),
        ):
            frames, sort_primary = exporter.frames_from_molecule_id(
                include_forces=True
            )
        assert sort_primary == ("b3lyp", "sto-3g")
        exporter.db.pick_primary_forces_method_basis.assert_not_called()

    def test_auto_pick_finds_no_forces_logs_info(self):
        exporter = _make_exporter(molecule_id="mid123")
        exporter.db.get_molecule_by_partial_id.return_value = "mid123full"
        exporter.db.get_structures_for_molecule.return_value = [
            {"structure_id": "sid1"}
        ]
        exporter.db.pick_primary_forces_method_basis.return_value = (
            None,
            None,
        )
        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                return_value=[],
            ),
        ):
            frames, sort_primary = exporter.frames_from_molecule_id(
                include_forces=True
            )
        assert sort_primary is None

    def test_auto_pick_finds_forces_logs_info_and_sets_sort_primary(self):
        exporter = _make_exporter(molecule_id="mid123")
        exporter.db.get_molecule_by_partial_id.return_value = "mid123full"
        exporter.db.get_structures_for_molecule.return_value = [
            {"structure_id": "sid1"}
        ]
        exporter.db.pick_primary_forces_method_basis.return_value = (
            "wb97xd",
            "def2-tzvp",
        )
        exporter.db.get_forces_for_structure_at.return_value = (None, None)
        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                return_value=[],
            ),
        ):
            frames, sort_primary = exporter.frames_from_molecule_id(
                include_forces=True
            )
        assert sort_primary == ("wb97xd", "def2-tzvp")

    def test_without_forces_sort_primary_is_none(self):
        exporter = _make_exporter(molecule_id="mid123")
        exporter.db.get_molecule_by_partial_id.return_value = "mid123full"
        exporter.db.get_structures_for_molecule.return_value = [
            {"structure_id": "sid1"}
        ]
        with (
            _patch_database_file(build_molecule=lambda d: FakeMolecule()),
            patch(
                "chemsmart.database.export.collect_energies_for_structure",
                return_value=[],
            ),
        ):
            frames, sort_primary = exporter.frames_from_molecule_id(
                include_forces=False
            )
        assert sort_primary is None


class TestBuildFrame:
    def test_build_frame_with_forces(self):
        exporter = _make_exporter()
        exporter.db.get_forces_for_structure_at.return_value = (
            [[0.0, 0.0, 0.0]],
            -1.0,
        )
        mock_db_file = MagicMock()
        mock_db_file.build_molecule_from_database.return_value = FakeMolecule()
        with patch(
            "chemsmart.database.export.collect_energies_for_structure",
            return_value=[],
        ):
            frame = exporter._build_frame(
                mock_db_file,
                {"structure_id": "sid1"},
                "b3lyp",
                "sto-3g",
                True,
            )
        assert frame["primary_energy"] == -1.0
        assert frame["primary_method"] == "b3lyp"
        assert frame["primary_basis"] == "sto-3g"

    def test_build_frame_without_forces_flag(self):
        exporter = _make_exporter()
        mock_db_file = MagicMock()
        mock_db_file.build_molecule_from_database.return_value = FakeMolecule()
        with patch(
            "chemsmart.database.export.collect_energies_for_structure",
            return_value=[],
        ):
            frame = exporter._build_frame(
                mock_db_file, {"structure_id": "sid1"}, None, None, False
            )
        assert frame["primary_energy"] is None
        exporter.db.get_forces_for_structure_at.assert_not_called()


class TestWriteXyzFrame:
    def test_writes_sid_and_formula_and_energies(self):
        exporter = _make_exporter(output="out.xyz")
        mol = FakeMolecule(chemical_formula="H2O")
        frame = {
            "molecule": mol,
            "structure_id": "sid123456789",
            "energies": [("b3lyp", "sto-3g", -1.5)],
        }
        import io

        buf = io.StringIO()
        exporter._write_xyz_frame(buf, frame)
        content = buf.getvalue()
        assert "SID: sid123456789" in content
        assert "Empirical formula: H2O" in content
        assert "Energy(b3lyp/sto-3g):" in content

    def test_writes_without_sid_or_formula_or_energies(self):
        exporter = _make_exporter(output="out.xyz")
        mol = FakeMolecule(chemical_formula=None)
        frame = {"molecule": mol, "structure_id": None, "energies": []}
        import io

        buf = io.StringIO()
        exporter._write_xyz_frame(buf, frame)
        content = buf.getvalue()
        assert "SID:" not in content
        assert "Empirical formula:" not in content
        assert "Energy(" not in content


class TestResolveId:
    def test_index_not_found_raises(self):
        inspector = _make_inspector(index=3)
        inspector.db.get_record.return_value = None
        with pytest.raises(ValueError, match="No record found with index"):
            inspector.resolve_id()

    def test_index_found_returns_record_id(self):
        inspector = _make_inspector(index=3)
        inspector.db.get_record.return_value = {"record_id": "rid123"}
        assert inspector.resolve_id() == "rid123"

    def test_record_id_path_delegates_to_partial_id_lookup(self):
        inspector = _make_inspector(record_id="abc")
        inspector.db.get_record_by_partial_id.return_value = "abcfull"
        assert inspector.resolve_id() == "abcfull"

    def test_neither_index_nor_record_id_raises(self):
        inspector = _make_inspector()
        with pytest.raises(ValueError, match="Either index or record_id"):
            inspector.resolve_id()


class TestRecordDetail:
    def test_record_not_found_raises(self):
        inspector = _make_inspector(index=1)
        inspector.db.get_record.side_effect = [
            {"record_id": "rid1"},  # resolve_id's lookup
            None,  # record_detail's own lookup
        ]
        with pytest.raises(ValueError, match="Record not found"):
            inspector.record_detail()


class TestFormatOverview:
    def _base_stats(self, **overrides):
        stats = {
            "db_file": "my.db",
            "db_size": 2048,
            "num_records": 2,
            "num_molecules": 1,
            "num_structures": 3,
            "programs": [("gaussian", 2)],
            "methods": [],
            "basis_sets": [],
            "assembled_first": "2024-01-01T00:00:00",
            "assembled_last": "2024-01-02T00:00:00",
            "gas_phase_records": 2,
            "solvated_records": 0,
            "solvents": [],
            "jobtypes": [("opt", 2)],
            "energy_count": 2,
            "forces_count": 0,
            "thermo_count": 1,
        }
        stats.update(overrides)
        return stats

    def test_empty_methods_and_basis_sets_show_null(self):
        inspector = _make_inspector()
        with patch.object(
            inspector, "overview", return_value=self._base_stats()
        ):
            text = inspector.format_overview()
        assert "Methods" in text
        assert "NULL" in text
        assert "Basis Sets" in text

    def test_methods_basis_and_solvents_present(self):
        inspector = _make_inspector()
        stats = self._base_stats(
            methods=[("b3lyp", 2)],
            basis_sets=[("sto-3g", 2)],
            solvents=[("water", 1)],
            solvated_records=1,
        )
        with patch.object(inspector, "overview", return_value=stats):
            text = inspector.format_overview()
        assert "b3lyp (2)" in text
        assert "sto-3g (2)" in text
        assert "water (1)" in text


class TestFormatRecordDetail:
    def _base_record(self, **overrides):
        record = {
            "record_index": 1,
            "record_id": "rid123456789",
            "meta": {
                "method": "b3lyp",
                "basis": "sto-3g",
                "spin": "restricted",
                "jobtype": "opt",
                "solvent_on": False,
                "route_string": "#p b3lyp/sto-3g opt",
            },
            "results": {
                "total_energy": -1.0,
                "num_unpaired_electrons": 0,
                "homo_energy": -5.0,
                "lumo_energy": 1.0,
                "fmo_gap": 6.0,
            },
            "provenance": {
                "source": "/tmp/job.log",
                "program": "gaussian",
                "program_version": "16",
                "parser": "cclib",
                "chemsmart_version": "1.0",
                "source_file_hash": "abc",
                "source_size": 1024,
                "source_date": "2024-01-01",
                "assembled_at": "2024-01-01T00:00:00",
            },
            "molecules": [],
        }
        record["meta"].update(overrides.pop("meta", {}))
        record["results"].update(overrides.pop("results", {}))
        record["provenance"].update(overrides.pop("provenance", {}))
        record.update(overrides)
        return record

    def _run(self, inspector, record):
        with patch.object(inspector, "record_detail", return_value=record):
            return inspector.format_record_detail()

    def test_closed_shell_basic_fields(self):
        inspector = _make_inspector()
        text = self._run(inspector, self._base_record())
        assert "HOMO Energy (eV)" in text
        assert "LUMO Energy (eV)" in text
        assert "α-HOMO" not in text

    def test_open_shell_shows_alpha_beta_and_somo(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "num_unpaired_electrons": 1,
                "alpha_homo_energy": -5.0,
                "beta_homo_energy": -4.5,
                "alpha_lumo_energy": 1.0,
                "beta_lumo_energy": 1.5,
                "alpha_fmo_gap": 6.0,
                "beta_fmo_gap": 6.0,
                "somo_energies": [-4.8],
            }
        )
        text = self._run(inspector, record)
        assert "α-HOMO Energy (eV)" in text
        assert "β-HOMO Energy (eV)" in text
        assert "SOMO Energies (eV)" in text

    def test_open_shell_without_somo_energies(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "num_unpaired_electrons": 1,
                "alpha_homo_energy": -5.0,
                "beta_homo_energy": -4.5,
                "alpha_lumo_energy": 1.0,
                "beta_lumo_energy": 1.5,
                "alpha_fmo_gap": 6.0,
                "beta_fmo_gap": 6.0,
            }
        )
        text = self._run(inspector, record)
        assert "SOMO Energies (eV)" not in text

    def test_customized_basis_without_detail_shows_unavailable_note(self):
        inspector = _make_inspector()
        record = self._base_record(meta={"basis": "customized_basis"})
        text = self._run(inspector, record)
        assert "unavailable" in text

    def test_customized_basis_with_full_detail(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": ["H", "C"],
                    "light_elements_basis": "6-31G",
                    "heavy_elements": ["Fe"],
                    "heavy_elements_basis": {
                        "Fe": [{"shell": "S"}, {"shell": "P"}]
                    },
                    "heavy_elements_ecp": {"Fe": "def2-ecp"},
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" in text
        assert "[6-31G]" in text
        assert "Heavy Elements" in text
        assert "Fe Shells" in text
        assert "S, P" in text
        assert "Heavy Elements ECP" in text

    def test_customized_basis_light_elements_only_no_light_basis(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": ["H", "C"],
                    "light_elements_basis": None,
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" in text
        assert "[" not in text.split("Light Elements")[1].split("\n")[0]
        assert "Heavy Elements" not in text

    def test_customized_basis_heavy_elements_only_no_light_elements(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "basis": "customized_basis",
                "custom_basis": {
                    "light_elements": None,
                    "heavy_elements": ["Fe"],
                },
            }
        )
        text = self._run(inspector, record)
        assert "Light Elements" not in text
        assert "Heavy Elements" in text

    def test_solvent_on_with_custom_solvent_detail(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": {
                    "SolventName": "MySolvent",
                    "Eps": 78.0,
                    "EpsInf": 1.8,
                    "HbondAcidity": 0.1,
                    "HbondBasicity": 0.2,
                    "SurfaceTensionAtInterface": 0.3,
                    "CarbonAromaticity": 0.4,
                    "ElectronegativeHalogenicity": 0.5,
                },
            }
        )
        text = self._run(inspector, record)
        assert "Solvent Model" in text
        assert "Custom Solvent Name" in text
        assert "Eps" in text

    def test_solvent_on_custom_solvent_generic_name_is_skipped(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": {"SolventName": "Generic"},
            }
        )
        text = self._run(inspector, record)
        assert "Custom Solvent Name" not in text

    def test_solvent_on_without_custom_solvent_dict(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={
                "solvent_on": True,
                "solvent_model": "PCM",
                "solvent_id": "water",
                "custom_solvent": None,
            }
        )
        text = self._run(inspector, record)
        assert "Solvent Model" in text
        assert "Custom Solvent Name" not in text

    def test_performance_section_shown_when_core_hours_present(self):
        inspector = _make_inspector()
        record = self._base_record(
            results={
                "total_core_hours": 12.5,
                "total_elapsed_walltime": 3.2,
            }
        )
        text = self._run(inspector, record)
        assert "Core Hours" in text
        assert "Elapsed Walltime" in text

    def test_thermochemistry_section_shown_when_gibbs_present(self):
        inspector = _make_inspector()
        record = self._base_record(
            meta={"temperature_in_K": 298.15, "pressure_in_atm": 1.0},
            results={
                "gibbs_free_energy": -1.5,
                "zero_point_energy": 0.01,
                "internal_energy": -1.4,
                "enthalpy": -1.3,
                "entropy": 0.002,
            },
        )
        text = self._run(inspector, record)
        assert "Thermochemistry" in text
        assert "Gibbs Free Energy" in text

    def test_normal_termination_shown_when_present(self):
        inspector = _make_inspector()
        record = self._base_record(provenance={"normal_termination": True})
        text = self._run(inspector, record)
        assert "Normal Termination" in text

    def test_uniform_structures_table(self):
        inspector = _make_inspector()
        record = self._base_record(
            molecules=[
                {
                    "index": 1,
                    "structure_id": "sid1",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": -1.0,
                    "is_optimized_structure": True,
                },
                {
                    "index": 2,
                    "structure_id": "sid2",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": None,
                    "is_optimized_structure": None,
                },
            ]
        )
        text = self._run(inspector, record)
        assert "Formula" in text
        assert "H2O" in text

    def test_non_uniform_structures_table(self):
        inspector = _make_inspector()
        record = self._base_record(
            molecules=[
                {
                    "index": 1,
                    "structure_id": "sid1",
                    "chemical_formula": "H2O",
                    "charge": 0,
                    "multiplicity": 1,
                    "number_of_atoms": 3,
                    "energy": -1.0,
                },
                {
                    "index": 2,
                    "structure_id": "sid2",
                    "chemical_formula": "H2O+",
                    "charge": 1,
                    "multiplicity": 2,
                    "number_of_atoms": 3,
                    "energy": None,
                },
            ]
        )
        text = self._run(inspector, record)
        assert "Idx" in text
        assert "H2O+" in text


class TestFormatStructureDetail:
    def _base_record_and_struct(self, **struct_overrides):
        struct = {
            "structure_index_in_file": 5,
            "energy": -1.0,
            "is_optimized_structure": True,
            "structure_id": "sid123",
            "charge": 0,
            "multiplicity": 1,
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "smiles": "O",
            "chemical_symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        }
        struct.update(struct_overrides)
        record = {
            "record_index": 1,
            "record_id": "rid123",
            "provenance": {"source": "/tmp/job.log"},
        }
        return record, struct

    def _run(self, inspector, record, struct):
        with patch.object(
            inspector, "structure_detail", return_value=(record, struct)
        ):
            return inspector.format_structure_detail()

    def test_full_detail_with_all_optional_sections(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            frozen_atoms=[1, 2],
            center_of_mass=[0.0, 0.0, 0.1],
            moments_of_inertia=[1.0, 2.0, 3.0],
            rotational_symmetry_number=2,
            rotational_constants=[1.0e9, 2.0e9, 3.0e9],
            point_group="C2v",
            dipole_moment=[0.0, 0.0, 1.5],
            dipole_moment_magnitude=1.5,
            mulliken_atomic_charges={"O1": -0.5, "H2": 0.25, "H3": 0.25},
            mulliken_spin_densities={"O1": 0.0, "H2": 0.0, "H3": 0.0},
            num_vibrational_modes=3,
            vibrational_frequencies=[100.0, 200.0, 300.0],
        )
        text = self._run(inspector, record, struct)
        assert "Optimized" in text
        assert "Coordinates" in text
        assert "Frozen Atoms" in text
        assert "Point Group" in text
        assert "Rotational Symmetry Number" in text
        assert "Center of Mass" in text
        assert "Moments of Inertia" in text
        assert "Rotational Constants" in text
        assert "Dipole Moment" in text
        assert "Mulliken Population Analysis" in text
        assert "Vibrational Frequencies" in text

    def test_physical_properties_point_group_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(point_group="C2v")
        text = self._run(inspector, record, struct)
        assert "Point Group" in text
        assert "Center of Mass" not in text
        assert "Moments of Inertia" not in text
        assert "Rotational Constants" not in text
        assert "Rotational Symmetry Number" not in text

    def test_mulliken_charges_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            mulliken_atomic_charges={"O1": -0.5}
        )
        text = self._run(inspector, record, struct)
        assert "Charge" in text

    def test_mulliken_spin_only(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            mulliken_spin_densities={"O1": 0.1}
        )
        text = self._run(inspector, record, struct)
        assert "Spin" in text

    def test_vibrational_modes_without_frequencies_list(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            num_vibrational_modes=0, vibrational_frequencies=[]
        )
        text = self._run(inspector, record, struct)
        assert "Vibrational Frequencies" in text

    def test_minimal_fields_skip_optional_sections(self):
        inspector = _make_inspector()
        record, struct = self._base_record_and_struct(
            is_optimized_structure=None,
            chemical_symbols=[],
            positions=[],
        )
        text = self._run(inspector, record, struct)
        assert "Coordinates" not in text
        assert "Frozen Atoms" not in text
        assert "Physical Properties" not in text
        assert "Dipole Moment" not in text
        assert "Mulliken Population Analysis" not in text


class TestFormatMoleculeDetail:
    def _base_molecule(self, **overrides):
        molecule = {
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "is_aromatic": False,
            "is_chiral": False,
            "is_linear": False,
            "is_multicomponent": False,
            "smiles": "O",
            "inchi": "InChI=1S/H2O/h1H2",
        }
        molecule.update(overrides)
        return molecule

    def test_molecule_not_found_raises(self):
        inspector = _make_inspector(molecule_id="mid123")
        inspector.db.get_molecule_by_partial_id.return_value = "mid123full"
        inspector.db.get_molecule.return_value = None
        with pytest.raises(ValueError, match="Molecule not found"):
            inspector.molecule_detail()

    def test_no_structures_shows_none(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "(none)" in text

    def test_element_counts_shown(self):
        inspector = _make_inspector()
        molecule = self._base_molecule(
            element_counts={"O": 1, "H": 2}, number_of_atoms=3
        )
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "Composition" in text
        assert "H 2" in text

    def test_all_tags_active(self):
        inspector = _make_inspector()
        molecule = self._base_molecule(
            is_aromatic=True,
            is_chiral=True,
            is_linear=True,
            is_multicomponent=True,
        )
        with patch.object(
            inspector, "molecule_detail", return_value=(molecule, [], [])
        ):
            text = inspector.format_molecule_detail()
        assert "aromatic" in text
        assert "chiral" in text

    def test_uniform_structures_with_primary_mb(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": ("b3lyp", "sto-3g"),
                "primary_energy": -2.0,
            },
            {
                "structure_id": "sid2",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": ("b3lyp", "sto-3g"),
                "primary_energy": -1.0,
            },
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy[b3lyp/sto-3g]" in text
        assert "ΔE (kcal/mol)" in text

    def test_structures_without_primary_mb_uses_default_header(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": None,
            }
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy (Eh)" in text
        assert "N/A" in text

    def test_structures_present_but_sorted_structures_empty(self):
        # sort_structure_dicts_by_energy could in principle return an
        # empty list even though `structures` itself is non-empty; the
        # `if sorted_structures:` guard exists specifically for that.
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": None,
            }
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=[],
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "Energy (Eh)" in text

    def test_mixed_electronic_states_note(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        structures = [
            {
                "structure_id": "sid1",
                "charge": 0,
                "multiplicity": 1,
                "primary_method_basis": None,
                "primary_energy": -1.0,
            },
            {
                "structure_id": "sid2",
                "charge": 1,
                "multiplicity": 2,
                "primary_method_basis": None,
                "primary_energy": -0.5,
            },
        ]
        with (
            patch.object(
                inspector,
                "molecule_detail",
                return_value=(molecule, structures, []),
            ),
            patch(
                "chemsmart.database.inspect.sort_structure_dicts_by_energy",
                return_value=structures,
            ),
        ):
            text = inspector.format_molecule_detail()
        assert "ΔE (—)" in text
        assert "Mixed electronic states detected" in text

    def test_related_records_table_populated(self):
        inspector = _make_inspector()
        molecule = self._base_molecule()
        records = [
            {
                "record_id": "rid123456789",
                "jobtype": "opt",
                "program": "gaussian",
                "method": "b3lyp",
                "basis": "sto-3g",
                "total_energy": -1.0,
            }
        ]
        with patch.object(
            inspector,
            "molecule_detail",
            return_value=(molecule, [], records),
        ):
            text = inspector.format_molecule_detail()
        assert "rid123456789"[:12] in text
        assert "gaussian" in text


class TestStandaloneStructureDetail:
    def test_structure_not_found_raises(self):
        inspector = _make_inspector(structure_id="sid123")
        inspector.db.get_structure_by_partial_id.return_value = "sid123full"
        inspector.db.get_structure.return_value = None
        with pytest.raises(ValueError, match="Structure not found"):
            inspector.standalone_structure_detail()


class TestFormatStandaloneStructureDetail:
    def _base_struct(self, **overrides):
        struct = {
            "structure_id": "sid123",
            "charge": 0,
            "multiplicity": 1,
            "molecule_id": "mid123",
            "chemical_formula": "H2O",
            "mass": 18.0153,
            "smiles": "O",
            "chemical_symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        }
        struct.update(overrides)
        return struct

    def _run(self, inspector, struct, records=None):
        with patch.object(
            inspector,
            "standalone_structure_detail",
            return_value=(struct, records or []),
        ):
            return inspector.format_standalone_structure_detail()

    def test_coordinates_and_no_physical_properties(self):
        inspector = _make_inspector()
        text = self._run(inspector, self._base_struct())
        assert "Coordinates" in text
        assert "Physical Properties" not in text

    def test_com_only(self):
        inspector = _make_inspector()
        struct = self._base_struct(center_of_mass=[0.0, 0.0, 0.1])
        text = self._run(inspector, struct)
        assert "Physical Properties" in text
        assert "Center of Mass" in text
        assert "Moments of Inertia" not in text

    def test_moi_only(self):
        inspector = _make_inspector()
        struct = self._base_struct(moments_of_inertia=[1.0, 2.0, 3.0])
        text = self._run(inspector, struct)
        assert "Moments of Inertia" in text
        assert "Center of Mass" not in text

    def test_no_coordinates(self):
        inspector = _make_inspector()
        struct = self._base_struct(chemical_symbols=[], positions=[])
        text = self._run(inspector, struct)
        assert "Coordinates" not in text

    def test_related_records_included(self):
        inspector = _make_inspector()
        records = [
            {
                "record_id": "rid987654321",
                "jobtype": "sp",
                "program": "orca",
                "method": "wb97xd",
                "basis": "def2tzvp",
                "total_energy": -2.0,
            }
        ]
        text = self._run(inspector, self._base_struct(), records)
        assert "orca" in text


class TestFormatRelatedRecordsTable:
    def test_empty_records_shows_none(self):
        inspector = _make_inspector()
        lines = inspector._format_related_records_table([])
        assert any("(none)" in line for line in lines)

    def test_records_with_missing_optional_fields(self):
        inspector = _make_inspector()
        lines = inspector._format_related_records_table(
            [{"record_id": "rid1"}]
        )
        text = "\n".join(lines)
        assert "rid1" in text
