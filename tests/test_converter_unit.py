"""
Direct unit tests for chemsmart.io.converter.FileConverter.

Complements tests/test_converter.py (which exercises the common
log/com/xyz conversion paths end-to-end) by covering error branches,
the less-common directory/file type dispatches (gjf, out with
--program, inp, sdf, pdb), include_intermediate_structures handling,
and xyz_to_pdb (Open Babel-based conversion).
"""

from shutil import copyfile

import pytest

from chemsmart.io.converter import FileConverter
from chemsmart.io.molecules.structure import Molecule

# Minimal V2000 molfile (water) for SDF-based tests. There is no .sdf
# fixture file checked into tests/data, so build one directly.
WATER_SDF_CONTENT = """water
  Generated

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7600    0.5900 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7600    0.5900 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
"""


class TestConvertFilesValidation:
    def test_directory_out_type_without_program_raises(self, tmp_path):
        converter = FileConverter(directory=str(tmp_path), type="out")
        with pytest.raises(ValueError, match="Both --filetype out"):
            converter.convert_files()

    def test_neither_directory_nor_filename_raises(self):
        converter = FileConverter()
        with pytest.raises(ValueError, match="Either directory or filename"):
            converter.convert_files()


class TestConvertAllFilesDirectoryDispatch:
    def test_gjf_directory_conversion(self, tmp_path, hf_com_filepath):
        gjf_path = tmp_path / "hf.gjf"
        copyfile(hf_com_filepath, gjf_path)
        converter = FileConverter(
            directory=str(tmp_path), type="gjf", output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "hf.xyz").is_file()

    def test_out_directory_conversion_gaussian_program(
        self, tmp_path, gaussian_singlet_opt_outfile
    ):
        out_path = tmp_path / "nhc.out"
        copyfile(gaussian_singlet_opt_outfile, out_path)
        converter = FileConverter(
            directory=str(tmp_path),
            type="out",
            program="gaussian",
            output_filetype="xyz",
        )
        converter.convert_files()
        assert (tmp_path / "nhc.xyz").is_file()

    def test_out_directory_conversion_orca_program(
        self, tmp_path, water_output_gas_path
    ):
        out_path = tmp_path / "water.out"
        copyfile(water_output_gas_path, out_path)
        converter = FileConverter(
            directory=str(tmp_path),
            type="out",
            program="orca",
            output_filetype="xyz",
        )
        converter.convert_files()
        assert (tmp_path / "water.xyz").is_file()

    def test_inp_directory_conversion(self, tmp_path, water_sp_input_path):
        inp_path = tmp_path / "water_sp.inp"
        copyfile(water_sp_input_path, inp_path)
        converter = FileConverter(
            directory=str(tmp_path), type="inp", output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "water_sp.xyz").is_file()

    def test_sdf_directory_conversion(self, tmp_path):
        target = tmp_path / "structure.sdf"
        target.write_text(WATER_SDF_CONTENT)
        converter = FileConverter(
            directory=str(tmp_path), type="sdf", output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "structure.xyz").is_file()

    def test_pdb_directory_conversion(self, tmp_path):
        pdb_content = (
            "HETATM    1  O   HOH A   7       0.000   0.000   0.000  1.00  0.00           O\n"
            "HETATM    2  H1  HOH A   7       0.960   0.000   0.000  1.00  0.00           H\n"
            "HETATM    3  H2  HOH A   7      -0.240   0.930   0.000  1.00  0.00           H\n"
            "END\n"
        )
        pdb_path = tmp_path / "water.pdb"
        pdb_path.write_text(pdb_content)
        converter = FileConverter(
            directory=str(tmp_path), type="pdb", output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "water.xyz").is_file()

    def test_unsupported_directory_type_raises(self, tmp_path):
        converter = FileConverter(directory=str(tmp_path), type="bogus")
        with pytest.raises(ValueError, match="is not supported"):
            converter.convert_files()

    def test_include_intermediate_structures_directory(
        self, tmp_path, gaussian_singlet_opt_outfile
    ):
        log_path = tmp_path / "nhc.log"
        copyfile(gaussian_singlet_opt_outfile, log_path)
        converter = FileConverter(
            directory=str(tmp_path),
            type="log",
            output_filetype="xyz",
            include_intermediate_structures=True,
        )
        converter.convert_files()
        assert (tmp_path / "nhc.xyz").is_file()


class TestConvertSingleFileDispatch:
    def test_out_could_not_detect_program_raises(self, tmp_path):
        fake_out = tmp_path / "unknown.out"
        fake_out.write_text("this is not a real quantum chemistry output\n")
        converter = FileConverter(
            filename=str(fake_out), output_filetype="xyz"
        )
        with pytest.raises(ValueError, match="Could not detect program"):
            converter.convert_files()

    def test_out_gaussian_single_file(
        self, tmp_path, gaussian_singlet_opt_outfile
    ):
        out_path = tmp_path / "nhc.out"
        copyfile(gaussian_singlet_opt_outfile, out_path)
        converter = FileConverter(
            filename=str(out_path), output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "nhc.xyz").is_file()

    def test_out_orca_single_file(self, tmp_path, water_output_gas_path):
        out_path = tmp_path / "water.out"
        copyfile(water_output_gas_path, out_path)
        converter = FileConverter(
            filename=str(out_path), output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "water.xyz").is_file()

    def test_gjf_single_file(self, tmp_path, hf_com_filepath):
        gjf_path = tmp_path / "hf.gjf"
        copyfile(hf_com_filepath, gjf_path)
        converter = FileConverter(
            filename=str(gjf_path), output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "hf.xyz").is_file()

    def test_inp_single_file(self, tmp_path, water_sp_input_path):
        inp_path = tmp_path / "water_sp.inp"
        copyfile(water_sp_input_path, inp_path)
        converter = FileConverter(
            filename=str(inp_path), output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "water_sp.xyz").is_file()

    def test_sdf_single_file(self, tmp_path):
        target = tmp_path / "structure.sdf"
        target.write_text(WATER_SDF_CONTENT)
        converter = FileConverter(filename=str(target), output_filetype="xyz")
        converter.convert_files()
        assert (tmp_path / "structure.xyz").is_file()

    def test_pdb_single_file(self, tmp_path):
        pdb_content = (
            "HETATM    1  O   HOH A   7       0.000   0.000   0.000  1.00  0.00           O\n"
            "HETATM    2  H1  HOH A   7       0.960   0.000   0.000  1.00  0.00           H\n"
            "HETATM    3  H2  HOH A   7      -0.240   0.930   0.000  1.00  0.00           H\n"
            "END\n"
        )
        pdb_path = tmp_path / "water.pdb"
        pdb_path.write_text(pdb_content)
        converter = FileConverter(
            filename=str(pdb_path), output_filetype="xyz"
        )
        converter.convert_files()
        assert (tmp_path / "water.xyz").is_file()

    def test_xyz_single_file_to_com(self, tmp_path):
        xyz_path = tmp_path / "water.xyz"
        xyz_path.write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n"
        )
        converter = FileConverter(
            filename=str(xyz_path), output_filetype="com"
        )
        converter.convert_files()
        assert (tmp_path / "water.com").is_file()

    def test_unsupported_single_file_type_raises(self, tmp_path):
        bogus = tmp_path / "file.bogus"
        bogus.write_text("data")
        converter = FileConverter(filename=str(bogus))
        with pytest.raises(ValueError, match="is not supported"):
            converter.convert_files()

    def test_include_intermediate_structures_single_file(
        self, tmp_path, gaussian_singlet_opt_outfile
    ):
        log_path = tmp_path / "nhc.log"
        copyfile(gaussian_singlet_opt_outfile, log_path)
        converter = FileConverter(
            filename=str(log_path),
            output_filetype="xyz",
            include_intermediate_structures=True,
        )
        converter.convert_files()
        assert (tmp_path / "nhc.xyz").is_file()


class TestXyzToPdb:
    def test_converts_molecule_via_temp_xyz(self, tmp_path):
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        pdb_path = tmp_path / "water.pdb"
        FileConverter.xyz_to_pdb(water, str(pdb_path))
        assert pdb_path.is_file()
        content = pdb_path.read_text()
        assert "ATOM" in content or "HETATM" in content

    def test_uses_existing_xyz_file_without_cleanup(self, tmp_path):
        xyz_path = tmp_path / "water.xyz"
        xyz_path.write_text(
            "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.76 0.59\nH 0.0 -0.76 0.59\n"
        )
        pdb_path = tmp_path / "water.pdb"
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        FileConverter.xyz_to_pdb(
            water, str(pdb_path), xyz_filename=str(xyz_path), cleanup=False
        )
        assert pdb_path.is_file()
        # Since cleanup=False and the xyz file already existed
        # (auto_xyz stays False), the source xyz file must remain.
        assert xyz_path.is_file()

    def test_missing_provided_xyz_file_is_written_first(self, tmp_path):
        xyz_path = tmp_path / "missing.xyz"
        pdb_path = tmp_path / "out.pdb"
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        FileConverter.xyz_to_pdb(
            water, str(pdb_path), xyz_filename=str(xyz_path)
        )
        assert pdb_path.is_file()
        assert xyz_path.is_file()

    def test_unreadable_xyz_raises_value_error(self, tmp_path):
        xyz_path = tmp_path / "empty.xyz"
        xyz_path.write_text("")
        pdb_path = tmp_path / "out.pdb"
        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        with pytest.raises(ValueError, match="Unable to read molecule"):
            FileConverter.xyz_to_pdb(
                water,
                str(pdb_path),
                xyz_filename=str(xyz_path),
                cleanup=False,
            )

    def test_unreadable_auto_generated_xyz_is_cleaned_up(self, tmp_path):
        # xyz_filename=None forces auto_xyz=True; mocking pybel.readfile
        # to fail simulates an unreadable auto-generated temp file, which
        # should still be removed since cleanup defaults to True.
        from unittest.mock import patch

        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        pdb_path = tmp_path / "out.pdb"
        with patch("openbabel.pybel.readfile", return_value=iter([])):
            with pytest.raises(ValueError, match="Unable to read molecule"):
                FileConverter.xyz_to_pdb(water, str(pdb_path))

    def test_cleanup_swallows_os_remove_error(self, tmp_path):
        from unittest.mock import patch

        water = Molecule(
            symbols=["O", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [0.0, 0.76, 0.59],
                [0.0, -0.76, 0.59],
            ],
        )
        pdb_path = tmp_path / "out.pdb"
        with patch("os.remove", side_effect=OSError("cannot remove")):
            # Should not raise despite os.remove failing during final
            # temp-file cleanup.
            FileConverter.xyz_to_pdb(water, str(pdb_path))
        assert pdb_path.is_file()
