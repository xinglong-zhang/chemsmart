"""
Direct unit tests for chemsmart.database.export.DatabaseExporter.

Complements tests/test_database.py::TestDatabaseExport (which exercises
full JSON/CSV/XYZ export through a real assembled+inserted database) by
mocking out the Database/DatabaseFile layer to reach selector-resolution
branches, forces handling, and extxyz formatting that are hard to set up
with real quantum-chemistry fixtures.
"""

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from chemsmart.database.export import DatabaseExporter


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
