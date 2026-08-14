"""Tests for chemsmart.io.database.DatabaseFile."""

import numpy as np
import pytest

from chemsmart.io.database import DatabaseFile


@pytest.fixture
def minimal_struct_dict():
    return {
        "chemical_symbols": ["H", "H"],
        "positions": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]],
        "charge": 0,
        "multiplicity": 1,
        "energy": -1.0,
    }


class TestBuildMoleculeFromDatabase:
    """Tests for DatabaseFile.build_molecule_from_database."""

    def test_minimal_struct_dict(self, minimal_struct_dict):
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert list(molecule.chemical_symbols) == ["H", "H"]
        assert molecule.charge == 0
        assert molecule.multiplicity == 1
        assert molecule.energy == -1.0

    def test_with_vibrational_modes(self, minimal_struct_dict):
        minimal_struct_dict["vibrational_modes"] = [
            [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]
        ]
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert isinstance(molecule.vibrational_modes[0], np.ndarray)

    def test_without_vibrational_modes(self, minimal_struct_dict):
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert not molecule.vibrational_modes

    def test_with_forces(self, minimal_struct_dict):
        minimal_struct_dict["forces"] = [[0.0, 0.0, 0.1], [0.0, 0.0, -0.1]]
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert isinstance(molecule.forces, np.ndarray)

    def test_without_forces(self, minimal_struct_dict):
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert molecule.forces is None

    def test_with_dipole_moment_and_rotational_constants(
        self, minimal_struct_dict
    ):
        minimal_struct_dict["dipole_moment"] = [0.1, 0.2, 0.3]
        minimal_struct_dict["rotational_constants"] = [1.0, 2.0, 3.0]
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.build_molecule_from_database(minimal_struct_dict)
        assert isinstance(molecule.dipole_moment, np.ndarray)
        assert isinstance(molecule.rotational_constants, np.ndarray)

    def test_without_positions_raises(self, minimal_struct_dict):
        minimal_struct_dict.pop("positions")
        db_file = DatabaseFile("dummy.db")
        with pytest.raises(ValueError, match="positions"):
            db_file.build_molecule_from_database(minimal_struct_dict)


class TestGetMoleculeByStructureId:
    """Tests for DatabaseFile.get_molecule_by_structure_id."""

    def test_returns_single_molecule(self, mocker, minimal_struct_dict):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_structure_by_partial_id.return_value = "full_sid_123"
        mock_db.get_structure.return_value = minimal_struct_dict

        db_file = DatabaseFile("dummy.db")
        molecule = db_file.get_molecule_by_structure_id("sid_123")

        mock_db.get_structure_by_partial_id.assert_called_once_with("sid_123")
        mock_db.get_structure.assert_called_once_with("full_sid_123")
        assert list(molecule.chemical_symbols) == ["H", "H"]

    def test_returns_list_when_requested(self, mocker, minimal_struct_dict):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_structure_by_partial_id.return_value = "full_sid_123"
        mock_db.get_structure.return_value = minimal_struct_dict

        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecule_by_structure_id(
            "sid_123", return_list=True
        )
        assert isinstance(molecules, list)
        assert len(molecules) == 1

    def test_raises_when_structure_not_found(self, mocker):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_structure_by_partial_id.return_value = "full_sid_123"
        mock_db.get_structure.return_value = None

        db_file = DatabaseFile("dummy.db")
        with pytest.raises(ValueError, match="No structure found"):
            db_file.get_molecule_by_structure_id("sid_123")


class TestGetMoleculesByMoleculeId:
    """Tests for DatabaseFile.get_molecules_by_molecule_id."""

    def test_returns_single_molecule_when_one_structure(
        self, mocker, minimal_struct_dict
    ):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_molecule_by_partial_id.return_value = "full_mid_123"
        mock_db.get_structures_for_molecule.return_value = [
            minimal_struct_dict
        ]
        mocker.patch(
            "chemsmart.io.database.sort_structure_dicts_by_energy",
            return_value=[minimal_struct_dict],
        )

        db_file = DatabaseFile("dummy.db")
        molecule = db_file.get_molecules_by_molecule_id("mid_123")
        assert list(molecule.chemical_symbols) == ["H", "H"]

    def test_returns_list_when_multiple_structures(
        self, mocker, minimal_struct_dict
    ):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_molecule_by_partial_id.return_value = "full_mid_123"
        struct_dicts = [minimal_struct_dict, dict(minimal_struct_dict)]
        mock_db.get_structures_for_molecule.return_value = struct_dicts
        mocker.patch(
            "chemsmart.io.database.sort_structure_dicts_by_energy",
            return_value=struct_dicts,
        )

        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecules_by_molecule_id("mid_123")
        assert isinstance(molecules, list)
        assert len(molecules) == 2

    def test_returns_list_explicitly_when_single_structure(
        self, mocker, minimal_struct_dict
    ):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_molecule_by_partial_id.return_value = "full_mid_123"
        mock_db.get_structures_for_molecule.return_value = [
            minimal_struct_dict
        ]
        mocker.patch(
            "chemsmart.io.database.sort_structure_dicts_by_energy",
            return_value=[minimal_struct_dict],
        )

        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecules_by_molecule_id(
            "mid_123", return_list=True
        )
        assert isinstance(molecules, list)
        assert len(molecules) == 1

    def test_raises_when_no_structures_found(self, mocker):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_molecule_by_partial_id.return_value = "full_mid_123"
        mock_db.get_structures_for_molecule.return_value = []

        db_file = DatabaseFile("dummy.db")
        with pytest.raises(ValueError, match="No structures found"):
            db_file.get_molecules_by_molecule_id("mid_123")


class TestGetMoleculesByRecord:
    """Tests for DatabaseFile.get_molecules_by_record."""

    def test_integer_index_selects_single_molecule(
        self, mocker, minimal_struct_dict
    ):
        mocker.patch("chemsmart.io.database.Database")
        mocker.patch(
            "chemsmart.io.database.resolve_record",
            return_value=[
                {"molecules": [minimal_struct_dict, dict(minimal_struct_dict)]}
            ],
        )
        db_file = DatabaseFile("dummy.db")
        molecule = db_file.get_molecules_by_record(
            record_index=1, structure_index="1"
        )
        assert list(molecule.chemical_symbols) == ["H", "H"]

    def test_slice_index_selects_multiple_molecules(
        self, mocker, minimal_struct_dict
    ):
        mocker.patch("chemsmart.io.database.Database")
        mocker.patch(
            "chemsmart.io.database.resolve_record",
            return_value=[
                {"molecules": [minimal_struct_dict, dict(minimal_struct_dict)]}
            ],
        )
        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecules_by_record(
            record_index=1, structure_index=":", return_list=True
        )
        assert len(molecules) == 2

    def test_skips_records_without_molecules(self, mocker):
        mocker.patch("chemsmart.io.database.Database")
        mocker.patch(
            "chemsmart.io.database.resolve_record",
            return_value=[{"molecules": []}],
        )
        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecules_by_record(
            record_index=1, structure_index=":", return_list=True
        )
        assert molecules == []

    def test_multiple_records_with_slice_index(
        self, mocker, minimal_struct_dict
    ):
        mocker.patch("chemsmart.io.database.Database")
        mocker.patch(
            "chemsmart.io.database.resolve_record",
            return_value=[
                {"molecules": [minimal_struct_dict]},
                {"molecules": [dict(minimal_struct_dict)]},
            ],
        )
        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_molecules_by_record(
            record_index=1, structure_index=":", return_list=True
        )
        assert len(molecules) == 2

    def test_out_of_range_integer_index_raises(
        self, mocker, minimal_struct_dict
    ):
        mocker.patch("chemsmart.io.database.Database")
        mocker.patch(
            "chemsmart.io.database.resolve_record",
            return_value=[{"molecules": [minimal_struct_dict]}],
        )
        db_file = DatabaseFile("dummy.db")
        with pytest.raises(ValueError, match="out of"):
            db_file.get_molecules_by_record(
                record_index=1, structure_index="5"
            )


class TestGetAllMolecules:
    """Tests for DatabaseFile.get_all_molecules."""

    def test_returns_single_molecule(self, mocker, minimal_struct_dict):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_all_structures.return_value = [minimal_struct_dict]

        db_file = DatabaseFile("dummy.db")
        molecule = db_file.get_all_molecules()
        assert list(molecule.chemical_symbols) == ["H", "H"]

    def test_returns_list_for_multiple_structures(
        self, mocker, minimal_struct_dict
    ):
        mock_db_cls = mocker.patch("chemsmart.io.database.Database")
        mock_db = mock_db_cls.return_value
        mock_db.get_all_structures.return_value = [
            minimal_struct_dict,
            dict(minimal_struct_dict),
        ]

        db_file = DatabaseFile("dummy.db")
        molecules = db_file.get_all_molecules(return_list=True)
        assert isinstance(molecules, list)
        assert len(molecules) == 2
