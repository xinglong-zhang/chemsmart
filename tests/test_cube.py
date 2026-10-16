"""Tests for chemsmart.io.gaussian.cube."""

import pytest

from chemsmart.io.gaussian.cube import CubeFileOperator, GaussianCubeFile


class TestGaussianCubeFileProperties:
    """Tests for GaussianCubeFile parsing properties."""

    def test_num_atoms(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        assert cube.num_atoms == 2

    def test_chemical_symbols(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        assert list(cube.chemical_symbols) == ["N", "N"]

    def test_positions_length_matches_num_atoms(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        assert len(cube.positions) == cube.num_atoms

    def test_structure_returns_molecule(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        molecule = cube.structure
        assert list(molecule.chemical_symbols) == ["N", "N"]

    def test_coordinate_block_as_list_length(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        assert len(cube.coordinate_block_as_list) == cube.num_atoms

    def test_values_by_lines_is_nonempty(self, spin_cube_file):
        cube = GaussianCubeFile(filename=spin_cube_file)
        values = cube.values_by_lines
        assert len(values) > 0
        assert all(isinstance(v, float) for line in values for v in line)


class TestCubeFileOperatorInit:
    """Tests for CubeFileOperator initialization."""

    def test_default_output_filename_subtract(
        self, spin_cube_file, esp_cube_file
    ):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator.operation == "subtract"
        assert operator.output_cubefile.endswith("_subtract.cube")

    def test_default_output_filename_add(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(
            spin_cube_file, esp_cube_file, operation="add"
        )
        assert operator.output_cubefile.endswith("_add.cube")

    def test_explicit_output_filename(
        self, spin_cube_file, esp_cube_file, tmp_path
    ):
        out_file = str(tmp_path / "custom_output.cube")
        operator = CubeFileOperator(
            spin_cube_file, esp_cube_file, output_cubefile=out_file
        )
        assert operator.output_cubefile == out_file


class TestCubeFileOperatorChecks:
    """Tests for CubeFileOperator compatibility checks."""

    def test_natoms_matched_for_same_geometry(
        self, spin_cube_file, esp_cube_file
    ):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator._check_natoms_matched() is True

    def test_coordinate_origin_matched(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator._check_coordinate_origin_matched() is True

    def test_grid_points_matched(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator._check_grid_points_matched() is True

    def test_geometries_matched(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator._check_geometries_matched()

    def test_all_checked_true_for_compatible_files(
        self, spin_cube_file, esp_cube_file
    ):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        assert operator._all_checked()


class TestCubeFileOperatorArithmetic:
    """Tests for CubeFileOperator add/subtract value computation."""

    def test_add_values(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(
            spin_cube_file, esp_cube_file, operation="add"
        )
        added = operator.add_values()
        cube1_values = operator.cube1.values_by_lines
        cube2_values = operator.cube2.values_by_lines
        assert added[0][0] == pytest.approx(
            cube1_values[0][0] + cube2_values[0][0]
        )

    def test_subtract_values(self, spin_cube_file, esp_cube_file):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        subtracted = operator.subtract_values()
        cube1_values = operator.cube1.values_by_lines
        cube2_values = operator.cube2.values_by_lines
        assert subtracted[0][0] == pytest.approx(
            cube1_values[0][0] - cube2_values[0][0]
        )

    def test_check_value_lines_raises_on_mismatched_line_count(
        self, spin_cube_file, esp_cube_file, mocker
    ):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        operator.cube2 = mocker.Mock(values_by_lines=[[1.0]])
        with pytest.raises(ValueError, match="same number of lines"):
            operator._check_value_lines()

    def test_check_value_lines_raises_on_mismatched_line_length(
        self, spin_cube_file, esp_cube_file, mocker
    ):
        operator = CubeFileOperator(spin_cube_file, esp_cube_file)
        original_values = operator.cube1.values_by_lines
        operator.cube2 = mocker.Mock(
            values_by_lines=[[1.0]] * len(original_values)
        )
        with pytest.raises(ValueError, match="different lengths"):
            operator._check_value_lines()


class TestCubeFileOperatorWriteResults:
    """Tests for CubeFileOperator.write_results end-to-end."""

    def test_write_results_subtract(
        self, spin_cube_file, esp_cube_file, tmp_path
    ):
        out_file = str(tmp_path / "result_subtract.cube")
        operator = CubeFileOperator(
            spin_cube_file,
            esp_cube_file,
            operation="subtract",
            output_cubefile=out_file,
        )
        operator.write_results()

        result_cube = GaussianCubeFile(filename=out_file)
        assert result_cube.num_atoms == 2
        assert result_cube.cube_job_title == "Gaussian job density"

    def test_write_results_add(self, spin_cube_file, esp_cube_file, tmp_path):
        out_file = str(tmp_path / "result_add.cube")
        operator = CubeFileOperator(
            spin_cube_file,
            esp_cube_file,
            operation="add",
            output_cubefile=out_file,
        )
        operator.write_results()

        result_cube = GaussianCubeFile(filename=out_file)
        assert result_cube.num_atoms == 2

    def test_write_results_unknown_operation_raises(
        self, spin_cube_file, esp_cube_file, tmp_path
    ):
        out_file = str(tmp_path / "result_unknown.cube")
        operator = CubeFileOperator(
            spin_cube_file,
            esp_cube_file,
            operation="multiply",
            output_cubefile=out_file,
        )
        with pytest.raises(ValueError, match="Unknown operation"):
            operator.write_results()

    def test_write_results_raises_when_geometries_mismatched(
        self, spin_cube_file, esp_cube_file, tmp_path, mocker
    ):
        out_file = str(tmp_path / "result_mismatch.cube")
        operator = CubeFileOperator(
            spin_cube_file,
            esp_cube_file,
            output_cubefile=out_file,
        )
        mocker.patch.object(operator, "_all_checked", return_value=False)
        with pytest.raises(ValueError, match="geometries"):
            operator.write_results()
