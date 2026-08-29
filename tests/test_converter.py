import inspect
import os.path
from shutil import copy, copytree, rmtree
from unittest.mock import patch

import click
import numpy as np
import pytest
from click.testing import CliRunner

from chemsmart.cli.convert import convert
from chemsmart.cli.main import entry_point
from chemsmart.io.converter import FileConverter
from chemsmart.io.gaussian.folder import (
    GaussianInputFolder,
    GaussianOutputFolder,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pdb.pdbfile import PDBFile
from chemsmart.io.xyz.folder import XYZFolder


class TestConverter:

    def test_convert_log_folder_to_xyz(
        self, tmpdir, gaussian_outputs_test_directory
    ):
        # copy whole directory gaussian_outputs_test_directory to tmpdir
        tmp_log_folder = os.path.join(
            tmpdir, "gaussian_outputs_test_directory"
        )
        copytree(gaussian_outputs_test_directory, tmp_log_folder)

        # remove link folder in tmp_log_folder if exists
        link_folder = os.path.join(tmp_log_folder, "link")
        if os.path.exists(link_folder):
            rmtree(link_folder)
        ###### TODO: when the test for link jobs are fixed, this should be
        ###### removed and the test should pass for link jobs too

        file_converter = FileConverter(
            directory=tmp_log_folder, type="log", output_filetype="xyz"
        )
        file_converter.convert_files()

        # check if the files are converted
        g16_folder = GaussianOutputFolder(folder=tmp_log_folder)
        all_logfiles = g16_folder.all_log_files

        # check all .log files have been converted to .xyz files
        for file in all_logfiles:
            assert os.path.exists(file.replace(".log", ".xyz"))

        ozone_xyz = os.path.join(tmp_log_folder, "ozone.xyz")

        assert os.path.exists(ozone_xyz)
        with open(ozone_xyz, "r") as f:
            lines = f.readlines()
            assert len(lines) == 5  # 5 lines in the log file
            assert lines[0] == "3\n"  # first line is number of atoms

    def test_convert_log_folder_to_com(
        self, tmpdir, gaussian_outputs_test_directory
    ):
        # copy whole directory gaussian_outputs_test_directory to tmpdir
        tmp_log_folder = os.path.join(
            tmpdir, "gaussian_outputs_test_directory"
        )
        copytree(gaussian_outputs_test_directory, tmp_log_folder)

        # remove link folder in tmp_log_folder if exists
        link_folder = os.path.join(tmp_log_folder, "link")
        if os.path.exists(link_folder):
            rmtree(link_folder)
        ###### TODO: when the test for link jobs are fixed, this should be
        ###### removed and the test should pass for link jobs too

        file_converter = FileConverter(
            directory=tmp_log_folder, type="log", output_filetype="com"
        )

        file_converter.convert_files()

        # check all .log files have been converted to .com files
        g16_folder = GaussianOutputFolder(folder=tmp_log_folder)
        all_logfiles = g16_folder.all_log_files
        for file in all_logfiles:
            assert os.path.exists(file.replace(".log", ".com"))

        ozone_com = os.path.join(tmp_log_folder, "ozone.com")
        assert os.path.exists(ozone_com)
        with open(ozone_com, "r") as f:
            lines = f.readlines()
            assert len(lines) == 12
            assert lines[5].startswith("Generated from")

    def test_convert_com_folder_to_xyz(
        self, tmpdir, gaussian_inputs_test_directory
    ):
        # copy whole directory gaussian_pbc_inputs_test_directory to tmpdir
        tmp_com_folder = os.path.join(tmpdir, "gaussian_inputs_test_directory")
        copytree(gaussian_inputs_test_directory, tmp_com_folder)

        file_converter = FileConverter(
            directory=tmp_com_folder, type="com", output_filetype="xyz"
        )
        file_converter.convert_files()

        # check all .com files have been converted to .xyz files
        g16_folder = GaussianInputFolder(folder=tmp_com_folder)
        all_comfiles = g16_folder.all_com_files
        for file in all_comfiles:
            assert os.path.exists(file.replace(".com", ".xyz"))

        hf_xyz = os.path.join(tmp_com_folder, "hf.xyz")
        assert os.path.exists(hf_xyz)
        with open(hf_xyz, "r") as f:
            lines = f.readlines()
            assert len(lines) == 16
            assert lines[0] == "14\n"

        # files in subfolders
        genecp_xyz = os.path.join(tmp_com_folder, "genecp", "opt_genecp.xyz")
        assert os.path.exists(genecp_xyz)
        with open(genecp_xyz, "r") as f:
            lines = f.readlines()
            assert len(lines) == 17
            assert lines[0] == "15\n"

        additional_xyz = os.path.join(
            tmp_com_folder, "additional", "model_sp_input.xyz"
        )
        assert os.path.exists(additional_xyz)
        with open(additional_xyz, "r") as f:
            lines = f.readlines()
            assert len(lines) == 16
            assert lines[0] == "14\n"

    def test_convert_xyz_folder_to_com(self, tmpdir, xyz_directory):
        # copy whole directory xyz_directory to tmpdir
        tmp_xyz_folder = os.path.join(tmpdir, "xyz_directory")
        copytree(xyz_directory, tmp_xyz_folder)

        file_converter = FileConverter(
            directory=tmp_xyz_folder, type="xyz", output_filetype="com"
        )
        file_converter.convert_files()

        # check all .xyz files have been converted to .com files
        xyz_folder = XYZFolder(folder=tmp_xyz_folder)
        all_xyzfiles = xyz_folder.all_xyzfiles
        for file in all_xyzfiles:
            assert os.path.exists(file.replace(".xyz", ".com"))

    def test_convert_single_logfile_to_com(
        self, tmpdir, gaussian_singlet_opt_outfile
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "gaussian_singlet_opt.log")
        copy(gaussian_singlet_opt_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="com"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".com"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".com"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 40
        assert mol.chemical_formula == "C19H12F3I2N3O"
        assert np.isclose(mol.mass, 609.128, rtol=1e-4)  # in thermo branch

    def test_convert_single_link_opt_logfile_to_com(
        self, tmpdir, gaussian_link_opt_outputfile
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "gaussian_singlet_opt.log")
        copy(gaussian_link_opt_outputfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="com"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".com"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".com"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 2
        assert mol.chemical_formula == "O2"

    def test_convert_single_link_sp_logfile_to_xyz(
        self, tmpdir, gaussian_dna_link_sp_outputfile
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "dna_link_sp.log")
        copy(gaussian_dna_link_sp_outputfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 603
        assert mol.chemical_formula == "C191H241Cu2N59O96P14"
        assert mol.energy == -25900.214629

    def test_convert_single_link_opt_logfile_to_xyz(
        self,
        tmpdir,
        gaussian_dppeFeCl2_link_opt_outputfile,
        gaussian_dppeFeCl2_link_opt_failed_outputfile,
    ):
        # copy file to tmpdir
        tmp_path_normal_termination = os.path.join(
            tmpdir, "dppeFeCl2_opt_quintet_link_opt_link.log"
        )
        copy(
            gaussian_dppeFeCl2_link_opt_outputfile, tmp_path_normal_termination
        )
        assert os.path.exists(tmp_path_normal_termination)
        file_converter = FileConverter(
            filename=tmp_path_normal_termination, output_filetype="xyz"
        )

        file_converter.convert_files()
        assert os.path.exists(
            tmp_path_normal_termination.replace(".log", ".xyz")
        )
        mol = Molecule.from_filepath(
            tmp_path_normal_termination.replace(".log", ".xyz")
        )
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 55
        assert mol.chemical_formula == "C26H24Cl2FeP2"
        assert mol.energy == -3869.013518

        tmp_path_error_termination = os.path.join(
            tmpdir,
            "dppeFeCl2_phenyldioxazolone_opt_triplet_opt_error_termination_link.log",
        )
        copy(
            gaussian_dppeFeCl2_link_opt_failed_outputfile,
            tmp_path_error_termination,
        )
        assert os.path.exists(tmp_path_error_termination)
        file_converter = FileConverter(
            filename=tmp_path_error_termination, output_filetype="xyz"
        )

        file_converter.convert_files()
        assert os.path.exists(
            tmp_path_error_termination.replace(".log", ".xyz")
        )
        mol = Molecule.from_filepath(
            tmp_path_error_termination.replace(".log", ".xyz")
        )
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 72
        assert mol.chemical_formula == "C34H29Cl2FeNO3P2"
        assert mol.energy == -4456.134472

    def test_convert_single_link_ts_logfile_to_xyz(
        self, tmpdir, gaussian_link_ts_outputfile
    ):  # copy file to tmpdir
        tmp_path_ts_error_termination = os.path.join(
            tmpdir, "dppeFeCl2_opt_quintet_link_opt_link.log"
        )
        copy(gaussian_link_ts_outputfile, tmp_path_ts_error_termination)
        assert os.path.exists(tmp_path_ts_error_termination)
        file_converter = FileConverter(
            filename=tmp_path_ts_error_termination, output_filetype="xyz"
        )

        file_converter.convert_files()
        assert os.path.exists(
            tmp_path_ts_error_termination.replace(".log", ".xyz")
        )
        mol = Molecule.from_filepath(
            tmp_path_ts_error_termination.replace(".log", ".xyz")
        )
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 2
        assert mol.chemical_formula == "O2"
        assert mol.energy == -150.116584

    def test_convert_single_link_logfile_to_xyz(
        self, tmpdir, gaussian_link_sp_outfile
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "intervening_n_Ap_A.log")
        copy(gaussian_link_sp_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 603
        assert mol.chemical_formula == "C191H241Cu2N59O96P14"
        assert mol.energy == -25900.214629

    def test_convert_single_comfile_to_xyz(
        self, tmpdir, gaussian_opt_inputfile
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "gaussian_opt.com")
        copy(gaussian_opt_inputfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".com", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".com", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 14
        assert mol.chemical_formula == "C7H5ClO"

        with open(tmp_path.replace(".com", ".xyz"), "r") as f:
            lines = f.readlines()
            assert len(lines) == 16
            assert lines[0] == "14\n"
        # assert np.isclose(mol.mass, 609.128, rtol=1e-4)  # in thermo branch

    def test_convert_single_sp_log_file_to_xyz(
        self, gaussian_benzene_opt_outfile, tmpdir
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "benzene_sp.log")
        copy(gaussian_benzene_opt_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 12
        assert mol.chemical_formula == "C6H6"
        assert mol.energy == -231.977725

    def test_convert_single_opt_log_file_to_xyz(
        self, gaussian_acetone_opt_outfile, tmpdir
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "acetone_opt.log")
        copy(gaussian_acetone_opt_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 10
        assert mol.chemical_formula == "C3H6O"
        assert mol.energy == -192.919416

    def test_convert_single_wbi_log_file_to_xyz(self, wbi_outputfile, tmpdir):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "TS_5coord_XIII_wbi.log")
        copy(wbi_outputfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 128
        assert mol.chemical_formula == "C51H63NNiO9P2Si"
        assert mol.energy == -5189.249707

    def test_convert_single_failed_modred_log_file_to_xyz(
        self, gaussian_failed_modred_outfile, tmpdir
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "cage_free_failed_modred.log")
        copy(gaussian_failed_modred_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 10
        assert mol.chemical_formula == "C3H4O3"
        assert mol.energy == -341.883317

    def test_convert_single_failed_oniom_log_file_to_xyz(
        self, gaussian_oniom_outputfile, tmpdir
    ):
        # copy file to tmpdir
        tmp_path = os.path.join(tmpdir, "cation_failed_scan.log")
        copy(gaussian_oniom_outputfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 483
        assert mol.chemical_formula == "C155H180CuN53O82P12"
        assert mol.energy == -5300.535128


class TestPDBFile:

    # -------------------------------------------------------------------
    # Initialisation and representation
    # -------------------------------------------------------------------

    def test_init_stores_filename(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.filename == single_model_pdb_file

    def test_repr(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert repr(pdb) == f"PDBFile({single_model_pdb_file})"

    def test_str(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert "PDBFile object" in str(pdb)
        assert single_model_pdb_file in str(pdb)

    def test_filepath_resolves_absolute(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert os.path.isabs(pdb.filepath)

    # -------------------------------------------------------------------
    # Raw line access
    # -------------------------------------------------------------------

    def test_raw_lines_preserves_column_whitespace(
        self, single_model_pdb_file
    ):
        """raw_lines must not strip leading spaces (fixed-width PDB format)."""
        pdb = PDBFile(filename=single_model_pdb_file)
        for line in pdb.raw_lines:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                assert len(line) >= 54  # at least through z-coordinate

    def test_raw_lines_strips_trailing_newlines(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        for line in pdb.raw_lines:
            assert not line.endswith("\n")

    # -------------------------------------------------------------------
    # Single-model parsing
    # -------------------------------------------------------------------

    def test_molecule_returns_molecule_object(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule
        assert isinstance(mol, Molecule)

    def test_num_atoms(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.num_atoms == 3

    def test_symbols(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.symbols == ["O", "H", "H"]

    def test_positions(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        expected = np.array(
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]
        )
        assert np.allclose(pdb.molecule.positions, expected)

    def test_atom_names(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.atom_names == ["O", "H1", "H2"]

    def test_residue_names(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.residue_names == ["HOH", "HOH", "HOH"]

    def test_residue_numbers(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.residue_numbers == [7, 7, 7]

    def test_chain_ids(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.chain_ids == ["A", "A", "A"]

    def test_record_types(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        assert pdb.molecule.record_type == ["HETATM", "HETATM", "HETATM"]

    def test_info_dict_populated(self, single_model_pdb_file):
        pdb = PDBFile(filename=single_model_pdb_file)
        info = pdb.molecule.info
        assert "atom_name" in info
        assert "residue_name" in info
        assert "residue_number" in info
        assert "chain_id" in info
        assert "record_type" in info

    # -------------------------------------------------------------------
    # Multi-model parsing
    # -------------------------------------------------------------------

    def test_get_molecules_all(self, multi_model_pdb_file):
        pdb = PDBFile(filename=multi_model_pdb_file)
        models = pdb.get_molecules(index=":", return_list=True)
        assert isinstance(models, list)
        assert len(models) == 2

    def test_get_molecules_first(self, multi_model_pdb_file):
        pdb = PDBFile(filename=multi_model_pdb_file)
        mol = pdb.get_molecules(index="1")
        assert mol.chain_ids == ["A", "A"]

    def test_get_molecules_last(self, multi_model_pdb_file):
        pdb = PDBFile(filename=multi_model_pdb_file)
        mol = pdb.get_molecules(index="-1")
        assert mol.chain_ids == ["B", "B"]
        assert np.allclose(mol.positions[0], np.array([1.5, 2.5, 3.5]))

    def test_molecule_property_returns_last_model(self, multi_model_pdb_file):
        pdb = PDBFile(filename=multi_model_pdb_file)
        assert pdb.molecule.chain_ids == ["B", "B"]

    def test_return_list_wraps_single(self, multi_model_pdb_file):
        pdb = PDBFile(filename=multi_model_pdb_file)
        result = pdb.get_molecules(index="-1", return_list=True)
        assert isinstance(result, list)
        assert len(result) == 1

    # -------------------------------------------------------------------
    # Element inference
    # -------------------------------------------------------------------

    def test_blank_element_columns_infer_two_letter_elements(
        self, blank_element_pdb_file
    ):
        pdb = PDBFile(filename=blank_element_pdb_file)
        mol = pdb.molecule
        assert list(mol.symbols) == ["Fe", "Zn", "Cl", "C"]

    def test_infer_element_fe(self):
        assert PDBFile._infer_element_from_atom_name("FE") == "Fe"

    def test_infer_element_zn(self):
        assert PDBFile._infer_element_from_atom_name("ZN") == "Zn"

    def test_infer_element_cl(self):
        assert PDBFile._infer_element_from_atom_name("CL") == "Cl"

    def test_infer_element_ca_is_carbon(self):
        """CA is a biomolecular atom label (C-alpha), should resolve to C."""
        assert PDBFile._infer_element_from_atom_name("CA") == "C"

    def test_infer_element_leading_digit(self):
        """Leading digits should be stripped: 1H -> H."""
        assert PDBFile._infer_element_from_atom_name("1H") == "H"

    def test_infer_element_empty_raises(self):
        with pytest.raises(ValueError, match="Unable to infer"):
            PDBFile._infer_element_from_atom_name("")

    def test_infer_element_digits_only_raises(self):
        with pytest.raises(ValueError, match="Unable to infer"):
            PDBFile._infer_element_from_atom_name("123")

    # -------------------------------------------------------------------
    # Error handling
    # -------------------------------------------------------------------

    def test_empty_file_raises_value_error(self, empty_pdb_file):
        pdb = PDBFile(filename=empty_pdb_file)
        with pytest.raises(ValueError, match="No ATOM/HETATM records"):
            pdb.get_molecules()

    def test_molecule_property_raises_on_empty(self, empty_pdb_file):
        pdb = PDBFile(filename=empty_pdb_file)
        with pytest.raises(ValueError, match="No ATOM/HETATM records"):
            _ = pdb.molecule

    # -------------------------------------------------------------------
    # Writing
    # -------------------------------------------------------------------

    def test_write_creates_file(self, single_model_pdb_file, tmpdir):
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule

        output_path = os.path.join(str(tmpdir), "output.pdb")
        PDBFile.write(mol, output_path)

        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0

    def test_write_round_trip_preserves_atom_count(
        self, single_model_pdb_file, tmpdir
    ):
        """Write then re-read should give the same number of atoms."""
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule

        output_path = os.path.join(str(tmpdir), "round_trip.pdb")
        PDBFile.write(mol, output_path)

        pdb2 = PDBFile(filename=output_path)
        assert pdb2.num_atoms == pdb.num_atoms

    def test_write_round_trip_preserves_symbols(
        self, single_model_pdb_file, tmpdir
    ):
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule

        output_path = os.path.join(str(tmpdir), "round_trip.pdb")
        PDBFile.write(mol, output_path)

        pdb2 = PDBFile(filename=output_path)
        assert pdb2.molecule.symbols == mol.symbols

    def test_write_round_trip_preserves_positions(
        self, single_model_pdb_file, tmpdir
    ):
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule

        output_path = os.path.join(str(tmpdir), "round_trip.pdb")
        PDBFile.write(mol, output_path)

        pdb2 = PDBFile(filename=output_path)
        assert np.allclose(pdb2.molecule.positions, mol.positions, atol=1e-3)

    def test_write_output_contains_atom_records(
        self, single_model_pdb_file, tmpdir
    ):
        pdb = PDBFile(filename=single_model_pdb_file)
        mol = pdb.molecule

        output_path = os.path.join(str(tmpdir), "records.pdb")
        PDBFile.write(mol, output_path)

        with open(output_path, "r") as f:
            content = f.read()
        assert "HETATM" in content or "ATOM" in content
        assert "END" in content

    # -------------------------------------------------------------------
    # Backward compatibility (Molecule delegates to PDBFile)
    # -------------------------------------------------------------------

    def test_molecule_from_filepath_uses_pdbfile(self, single_model_pdb_file):
        """Molecule.from_filepath for .pdb should produce identical results."""
        mol_via_molecule = Molecule.from_filepath(single_model_pdb_file)
        pdb = PDBFile(filename=single_model_pdb_file)
        mol_via_pdbfile = pdb.molecule

        assert mol_via_molecule.symbols == mol_via_pdbfile.symbols
        assert np.allclose(
            mol_via_molecule.positions, mol_via_pdbfile.positions
        )
        assert mol_via_molecule.atom_names == mol_via_pdbfile.atom_names
        assert mol_via_molecule.residue_names == mol_via_pdbfile.residue_names
        assert (
            mol_via_molecule.residue_numbers == mol_via_pdbfile.residue_numbers
        )
        assert mol_via_molecule.chain_ids == mol_via_pdbfile.chain_ids

    def test_pdb_infer_pdb_element(self):
        assert PDBFile._infer_element_from_atom_name("FE") == "Fe"
        assert PDBFile._infer_element_from_atom_name("CA") == "C"

    def test_pdb_parse_pdb_models(self, multi_model_pdb_file):
        models = PDBFile(multi_model_pdb_file)._parse_models()
        assert len(models) == 2

    def test_pdb_molecule_from_pdb_atom_lines(self):
        atom_line = (
            "HETATM    1  O   HOH A   7"
            "       0.000   0.000   0.000"
            "  1.00  0.00           O"
        )
        mol = PDBFile._get_molecule_from_atom_lines([atom_line])
        assert mol.symbols == ["O"]

    # ------------------------------------------------------------------
    # CDXML / CDX conversion tests
    # ------------------------------------------------------------------

    def test_convert_single_cdxml_to_xyz(
        self, tmpdir, single_molecule_cdxml_file_methane
    ):
        tmp_path = os.path.join(tmpdir, "methane.cdxml")
        copy(single_molecule_cdxml_file_methane, tmp_path)

        FileConverter(filename=tmp_path, output_filetype="xyz").convert_files()

        output = tmp_path.replace(".cdxml", ".xyz")
        assert os.path.exists(output)
        mol = Molecule.from_filepath(output)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 5
        assert mol.chemical_formula == "CH4"

    def test_convert_single_cdxml_to_com(
        self, tmpdir, single_molecule_cdxml_file_benzene
    ):
        tmp_path = os.path.join(tmpdir, "benzene.cdxml")
        copy(single_molecule_cdxml_file_benzene, tmp_path)

        FileConverter(filename=tmp_path, output_filetype="com").convert_files()

        output = tmp_path.replace(".cdxml", ".com")
        assert os.path.exists(output)
        mol = Molecule.from_filepath(output)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 12
        assert mol.chemical_formula == "C6H6"

    def test_convert_single_cdx_to_xyz(
        self, tmpdir, single_molecule_cdx_file_imidazole
    ):
        tmp_path = os.path.join(tmpdir, "imidazole.cdx")
        copy(single_molecule_cdx_file_imidazole, tmp_path)

        FileConverter(filename=tmp_path, output_filetype="xyz").convert_files()

        output = tmp_path.replace(".cdx", ".xyz")
        assert os.path.exists(output)
        mol = Molecule.from_filepath(output)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 21
        assert mol.chemical_formula == "C8H10N2O"

    def test_convert_multi_molecule_cdxml_to_xyz_splits_files(
        self, tmpdir, multi_molecule_cdxml_file
    ):
        # Multi-molecule cdxml should produce basename_1.xyz, basename_2.xyz
        # when include_intermediate_structures=True.
        tmp_path = os.path.join(tmpdir, "two_molecules.cdxml")
        copy(multi_molecule_cdxml_file, tmp_path)

        FileConverter(
            filename=tmp_path,
            output_filetype="xyz",
            include_intermediate_structures=True,
        ).convert_files()

        output_1 = os.path.join(tmpdir, "two_molecules_1.xyz")
        output_2 = os.path.join(tmpdir, "two_molecules_2.xyz")
        assert os.path.exists(output_1)
        assert os.path.exists(output_2)

        mol1 = Molecule.from_filepath(output_1)
        assert isinstance(mol1, Molecule)
        assert mol1.chemical_formula == "CH2O"

        mol2 = Molecule.from_filepath(output_2)
        assert isinstance(mol2, Molecule)
        assert mol2.chemical_formula == "N2"
        assert mol2.num_atoms == 2

    def test_convert_multi_molecule_cdxml_default_writes_last_only(
        self, tmpdir, multi_molecule_cdxml_file
    ):
        """Without -z, multi-fragment ChemDraw writes only the last molecule."""
        tmp_path = os.path.join(tmpdir, "two_molecules.cdxml")
        copy(multi_molecule_cdxml_file, tmp_path)

        FileConverter(
            filename=tmp_path,
            output_filetype="xyz",
            include_intermediate_structures=False,
        ).convert_files()

        output = os.path.join(tmpdir, "two_molecules.xyz")
        output_1 = os.path.join(tmpdir, "two_molecules_1.xyz")
        output_2 = os.path.join(tmpdir, "two_molecules_2.xyz")
        assert os.path.exists(output)
        assert not os.path.exists(output_1)
        assert not os.path.exists(output_2)

        mol = Molecule.from_filepath(output)
        assert mol.chemical_formula == "N2"
        assert mol.num_atoms == 2

    def test_convert_cdxml_folder_to_xyz(self, tmpdir, chemdraw_directory):
        from shutil import copytree

        tmp_cdxml_folder = os.path.join(tmpdir, "chemdraw")
        copytree(chemdraw_directory, tmp_cdxml_folder)

        FileConverter(
            directory=tmp_cdxml_folder,
            type="cdxml",
            output_filetype="xyz",
            include_intermediate_structures=True,
        ).convert_files()

        # Single-molecule cdxml files produce basename.xyz
        for fname in (
            "benzene.cdxml",
            "methane.cdxml",
            "complex_molecule.cdxml",
        ):
            assert os.path.exists(
                os.path.join(tmp_cdxml_folder, fname.replace(".cdxml", ".xyz"))
            )

        # two_molecules.cdxml contains 2 molecules → split into _1.xyz and _2.xyz
        assert os.path.exists(
            os.path.join(tmp_cdxml_folder, "two_molecules_1.xyz")
        )
        assert os.path.exists(
            os.path.join(tmp_cdxml_folder, "two_molecules_2.xyz")
        )

    def test_convert_cdxml_folder_to_com(self, tmpdir, chemdraw_directory):
        from shutil import copytree

        tmp_cdxml_folder = os.path.join(tmpdir, "chemdraw")
        copytree(chemdraw_directory, tmp_cdxml_folder)

        FileConverter(
            directory=tmp_cdxml_folder,
            type="cdxml",
            output_filetype="com",
            include_intermediate_structures=True,
        ).convert_files()

        # Single-molecule cdxml files produce basename.com
        for fname in (
            "benzene.cdxml",
            "methane.cdxml",
            "complex_molecule.cdxml",
        ):
            assert os.path.exists(
                os.path.join(tmp_cdxml_folder, fname.replace(".cdxml", ".com"))
            )

        # two_molecules.cdxml contains 2 molecules → split into _1.com and _2.com
        assert os.path.exists(
            os.path.join(tmp_cdxml_folder, "two_molecules_1.com")
        )
        assert os.path.exists(
            os.path.join(tmp_cdxml_folder, "two_molecules_2.com")
        )


class TestConverterTwoWay:
    """Tests for the generic two-way ``FileConverter.convert_file`` method."""

    def test_convert_file_pdb_to_xyz(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".xyz")
        FileConverter.convert_file(single_model_pdb_file, output_path)

        assert os.path.exists(output_path)
        mol = Molecule.from_filepath(output_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_convert_file_xyz_to_pdb(self, single_model_pdb_file):
        xyz_path = single_model_pdb_file.replace(".pdb", ".xyz")
        FileConverter.convert_file(single_model_pdb_file, xyz_path)

        pdb_path = single_model_pdb_file.replace(".pdb", "_roundtrip.pdb")
        FileConverter.convert_file(xyz_path, pdb_path)

        assert os.path.exists(pdb_path)
        mol = Molecule.from_filepath(pdb_path)
        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_convert_file_pdb_to_com(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".com")
        FileConverter.convert_file(single_model_pdb_file, output_path)

        assert os.path.exists(output_path)
        mol = Molecule.from_filepath(output_path)
        assert mol.num_atoms == 3
        assert mol.chemical_formula == "H2O"

    def test_convert_file_com_to_xyz(self, tmpdir, gaussian_opt_inputfile):
        tmp_path = os.path.join(tmpdir, "gaussian_opt.com")
        copy(gaussian_opt_inputfile, tmp_path)

        output_path = os.path.join(tmp_path.replace(".com", ".xyz"))
        FileConverter.convert_file(tmp_path, output_path)

        assert os.path.exists(output_path)
        mol = Molecule.from_filepath(output_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 14
        assert mol.chemical_formula == "C7H5ClO"

    def test_convert_file_xyz_to_com(self, tmpdir, gaussian_opt_inputfile):
        tmp_path = os.path.join(tmpdir, "gaussian_opt.com")
        copy(gaussian_opt_inputfile, tmp_path)
        xyz_path = tmp_path.replace(".com", ".xyz")
        FileConverter.convert_file(tmp_path, xyz_path)

        com_path = tmp_path.replace(".com", "_from_xyz.com")
        FileConverter.convert_file(xyz_path, com_path)

        assert os.path.exists(com_path)
        mol = Molecule.from_filepath(com_path)
        assert mol.num_atoms == 14
        assert mol.chemical_formula == "C7H5ClO"

    def test_convert_file_multi_model_pdb_splits(self, multi_model_pdb_file):
        output_path = multi_model_pdb_file.replace(".pdb", ".xyz")
        FileConverter.convert_file(
            multi_model_pdb_file,
            output_path,
            include_intermediate_structures=True,
        )

        assert not os.path.exists(output_path)
        output_dir = os.path.dirname(output_path)
        output_basename = os.path.splitext(os.path.basename(output_path))[0]
        output_1 = os.path.join(output_dir, f"{output_basename}_1.xyz")
        output_2 = os.path.join(output_dir, f"{output_basename}_2.xyz")
        assert os.path.exists(output_1)
        assert os.path.exists(output_2)

        mol1 = Molecule.from_filepath(output_1)
        mol2 = Molecule.from_filepath(output_2)
        assert mol1.num_atoms == 2
        assert mol2.num_atoms == 2
        assert not np.allclose(mol1.positions, mol2.positions)

    def test_convert_file_single_structure_without_splitting(
        self, multi_model_pdb_file
    ):
        """Default mode writes a single output for multi-model inputs."""
        output_path = multi_model_pdb_file.replace(".pdb", ".xyz")
        FileConverter.convert_file(
            multi_model_pdb_file,
            output_path,
            include_intermediate_structures=False,
        )

        assert os.path.exists(output_path)
        output_dir = os.path.dirname(output_path)
        output_basename = os.path.splitext(os.path.basename(output_path))[0]
        assert not os.path.exists(
            os.path.join(output_dir, f"{output_basename}_1.xyz")
        )
        mol = Molecule.from_filepath(output_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 2

    def test_convert_files_routes_to_convert_file_with_output_filepath(
        self, tmpdir, single_model_pdb_file
    ):
        output_path = os.path.join(tmpdir, "water.xyz")
        fc = FileConverter(
            filename=single_model_pdb_file, output_filepath=output_path
        )
        fc.convert_files()

        assert os.path.exists(output_path)
        mol = Molecule.from_filepath(output_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 3

    def test_convert_files_routes_with_include_intermediate(
        self, multi_model_pdb_file, tmpdir
    ):
        output_path = os.path.join(tmpdir, "multi.xyz")
        fc = FileConverter(
            filename=multi_model_pdb_file,
            output_filepath=output_path,
            include_intermediate_structures=True,
        )
        fc.convert_files()

        assert not os.path.exists(output_path)
        assert os.path.exists(os.path.join(tmpdir, "multi_1.xyz"))
        assert os.path.exists(os.path.join(tmpdir, "multi_2.xyz"))

    def test_convert_file_missing_input_raises(self, tmpdir):
        missing_input = os.path.join(tmpdir, "missing.xyz")
        output_path = os.path.join(tmpdir, "out.xyz")
        with pytest.raises(FileNotFoundError):
            FileConverter.convert_file(missing_input, output_path)

    def test_convert_file_missing_extension_raises(
        self, single_model_pdb_file
    ):
        output_path = single_model_pdb_file.replace(".pdb", "")
        with pytest.raises(ValueError, match="file extension"):
            FileConverter.convert_file(single_model_pdb_file, output_path)

    def test_convert_file_creates_missing_output_dir(
        self, single_model_pdb_file, tmpdir
    ):
        output_dir = os.path.join(tmpdir, "nested", "out")
        output_path = os.path.join(output_dir, "water.xyz")
        assert not os.path.isdir(output_dir)

        FileConverter.convert_file(single_model_pdb_file, output_path)

        assert os.path.exists(output_path)
        mol = Molecule.from_filepath(output_path)
        assert mol.num_atoms == 3

    def test_convert_file_none_molecules_raises(
        self, single_model_pdb_file, tmpdir
    ):
        output_path = os.path.join(tmpdir, "out.xyz")
        with patch(
            "chemsmart.io.converter.Molecule.from_filepath",
            return_value=None,
        ):
            with pytest.raises(ValueError, match="No molecule could be read"):
                FileConverter.convert_file(single_model_pdb_file, output_path)

    def test_convert_file_none_molecules_with_intermediates_raises(
        self, single_model_pdb_file, tmpdir
    ):
        output_path = os.path.join(tmpdir, "out.xyz")
        with patch(
            "chemsmart.io.converter.Molecule.from_filepath",
            return_value=None,
        ):
            with pytest.raises(ValueError, match="No molecule could be read"):
                FileConverter.convert_file(
                    single_model_pdb_file,
                    output_path,
                    include_intermediate_structures=True,
                )

    def test_convert_file_empty_list_raises(
        self, single_model_pdb_file, tmpdir
    ):
        output_path = os.path.join(tmpdir, "out.xyz")
        with patch(
            "chemsmart.io.converter.Molecule.from_filepath",
            return_value=[],
        ):
            with pytest.raises(ValueError, match="No molecules found"):
                FileConverter.convert_file(
                    single_model_pdb_file,
                    output_path,
                    include_intermediate_structures=True,
                )

    def test_convert_file_non_list_molecule_wrapped(
        self, single_model_pdb_file, tmpdir
    ):
        """Hit the defensive ``not isinstance(molecules, list)`` branch."""
        output_path = os.path.join(tmpdir, "out.xyz")
        mol = Molecule.from_filepath(single_model_pdb_file)
        with patch(
            "chemsmart.io.converter.Molecule.from_filepath",
            return_value=mol,
        ):
            FileConverter.convert_file(
                single_model_pdb_file,
                output_path,
                include_intermediate_structures=True,
            )
        assert os.path.exists(output_path)
        written = Molecule.from_filepath(output_path)
        assert written.num_atoms == 3

    def test_convert_file_single_molecule_with_intermediates(
        self, single_model_pdb_file, tmpdir
    ):
        """include_intermediate_structures=True with one molecule writes
        directly to the output path (no _1 suffix)."""
        output_path = os.path.join(tmpdir, "water.xyz")
        FileConverter.convert_file(
            single_model_pdb_file,
            output_path,
            include_intermediate_structures=True,
        )
        assert os.path.exists(output_path)
        assert not os.path.exists(os.path.join(tmpdir, "water_1.xyz"))
        mol = Molecule.from_filepath(output_path)
        assert mol.num_atoms == 3


class TestResolveOutputPath:
    """Tests for ``FileConverter.resolve_output_path``."""

    def test_resolve_output_path_default(self):
        assert FileConverter.resolve_output_path(
            "/path/to/molecule.pdb"
        ) == os.path.join("/path/to", "molecule.xyz")

    def test_resolve_output_path_output_filetype(self):
        assert FileConverter.resolve_output_path(
            "/path/to/molecule.pdb", output_filetype="com"
        ) == os.path.join("/path/to", "molecule.com")

    def test_resolve_output_path_extension_with_dot(self):
        assert FileConverter.resolve_output_path(
            "/path/to/molecule.pdb", output=".mol2"
        ) == os.path.join("/path/to", "molecule.mol2")

    def test_resolve_output_path_extension_without_dot(self):
        assert FileConverter.resolve_output_path(
            "/path/to/molecule.pdb", output="pdb"
        ) == os.path.join("/path/to", "molecule.pdb")

    def test_resolve_output_path_full_path(self):
        assert (
            FileConverter.resolve_output_path(
                "/path/to/molecule.pdb", output="/tmp/foo.xyz"
            )
            == "/tmp/foo.xyz"
        )

    def test_resolve_output_path_relative_path_with_extension(self):
        assert (
            FileConverter.resolve_output_path(
                "/path/to/molecule.pdb", output="othername.xyz"
            )
            == "othername.xyz"
        )

    def test_resolve_output_path_subdirectory_path(self):
        assert (
            FileConverter.resolve_output_path(
                "/path/to/molecule.pdb", output="subdir/out.xyz"
            )
            == "subdir/out.xyz"
        )


class TestConvertCLI:
    """Tests for ``chemsmart run convert``."""

    def test_cli_convert_pdb_to_xyz(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                output_path,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

        mol = Molecule.from_filepath(output_path)
        assert mol.num_atoms == 3

    def test_cli_convert_short_flags(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".com")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "-i",
                single_model_pdb_file,
                "-o",
                output_path,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_missing_input_fails(self, tmpdir):
        missing_input = os.path.join(tmpdir, "missing.pdb")
        output_path = os.path.join(tmpdir, "out.xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                missing_input,
                "--output",
                output_path,
            ],
        )
        assert result.exit_code != 0

    def test_cli_convert_with_debug_and_stream(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", "_dbg.xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                output_path,
                "--debug",
                "--stream",
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_batch_directory(
        self, tmpdir, gaussian_inputs_test_directory
    ):
        tmp_com_folder = os.path.join(tmpdir, "gaussian_inputs_test_directory")
        copytree(gaussian_inputs_test_directory, tmp_com_folder)

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--directory",
                tmp_com_folder,
                "--filetype",
                "com",
                "--output-filetype",
                "xyz",
            ],
        )
        assert result.exit_code == 0, result.output

        g16_folder = GaussianInputFolder(folder=tmp_com_folder)
        for file in g16_folder.all_com_files:
            assert os.path.exists(file.replace(".com", ".xyz"))

    def test_cli_convert_include_intermediate_structures(
        self, multi_model_pdb_file, tmpdir
    ):
        output_path = os.path.join(tmpdir, "multi.xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                multi_model_pdb_file,
                "--output",
                output_path,
                "--include-intermediate-structures",
            ],
        )
        assert result.exit_code == 0, result.output
        assert not os.path.exists(output_path)
        assert os.path.exists(os.path.join(tmpdir, "multi_1.xyz"))
        assert os.path.exists(os.path.join(tmpdir, "multi_2.xyz"))

    def test_cli_convert_mutually_exclusive_input_and_directory(
        self, single_model_pdb_file, tmpdir
    ):
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                os.path.join(tmpdir, "out.xyz"),
                "--directory",
                str(tmpdir),
                "--filetype",
                "pdb",
            ],
        )
        assert result.exit_code != 0
        assert "either --input" in result.output.lower() or (
            "not both" in result.output.lower()
        )

    def test_cli_convert_input_defaults_to_xyz(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["run", "convert", "--input", single_model_pdb_file],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

        mol = Molecule.from_filepath(output_path)
        assert mol.num_atoms == 3

    def test_cli_convert_input_defaults_to_output_filetype(
        self, single_model_pdb_file
    ):
        output_path = single_model_pdb_file.replace(".pdb", ".com")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output-filetype",
                "com",
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_output_extension_only(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", ".com")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                ".com",
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_output_extension_without_dot(
        self, single_model_pdb_file
    ):
        output_path = single_model_pdb_file.replace(".pdb", ".pdb")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                "pdb",
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_output_full_path(self, single_model_pdb_file, tmpdir):
        output_path = os.path.join(str(tmpdir), "othername.xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                output_path,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)

    def test_cli_convert_directory_requires_filetype(self, tmpdir):
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            ["run", "convert", "--directory", str(tmpdir)],
        )
        assert result.exit_code != 0
        assert "--filetype" in result.output

    def test_cli_convert_requires_mode(self):
        runner = CliRunner()
        result = runner.invoke(entry_point, ["run", "convert"])
        assert result.exit_code != 0
        assert "either --input" in result.output.lower() or (
            "--directory/--filetype" in result.output.lower()
        )

    def _call_convert_callback(self, **kwargs):
        """Call the convert command function, bypassing Click wrapping."""
        convert_fn = inspect.unwrap(convert.callback)
        defaults = {
            "ctx": click.Context(convert),
            "input_file": None,
            "output_file": None,
            "directory": None,
            "filetype": None,
            "program": None,
            "output_filetype": "xyz",
            "include_intermediate_structures": False,
            "debug": False,
            "stream": False,
        }
        defaults.update(kwargs)
        return convert_fn(**defaults)

    def test_convert_callback_resolves_default_output(
        self, single_model_pdb_file
    ):
        """Direct callback call covers resolve_output_path in convert.py."""
        output_path = single_model_pdb_file.replace(".pdb", ".xyz")
        result = self._call_convert_callback(
            input_file=single_model_pdb_file,
        )
        assert result is None
        assert os.path.exists(output_path)

    def test_convert_callback_uses_output_file(self, single_model_pdb_file):
        output_path = single_model_pdb_file.replace(".pdb", "_named.com")
        result = self._call_convert_callback(
            input_file=single_model_pdb_file,
            output_file=output_path,
        )
        assert result is None
        assert os.path.exists(output_path)

    def test_convert_callback_rejects_input_and_directory(
        self, single_model_pdb_file, tmpdir
    ):
        with pytest.raises(click.UsageError, match="not both"):
            self._call_convert_callback(
                input_file=single_model_pdb_file,
                directory=str(tmpdir),
            )

    def test_convert_callback_directory_requires_filetype(self, tmpdir):
        with pytest.raises(click.UsageError, match="--filetype is required"):
            self._call_convert_callback(directory=str(tmpdir))

    def test_convert_callback_requires_mode(self):
        with pytest.raises(click.UsageError, match="either --input"):
            self._call_convert_callback()

    def test_convert_callback_batch_directory(self, tmpdir):
        with patch("chemsmart.cli.convert.FileConverter") as mock_fc:
            mock_fc.return_value.convert_files.return_value = None
            result = self._call_convert_callback(
                directory=str(tmpdir),
                filetype="pdb",
                program="gaussian",
                output_filetype="xyz",
                include_intermediate_structures=True,
            )
        assert result is None
        mock_fc.assert_called_once_with(
            directory=str(tmpdir),
            type="pdb",
            program="gaussian",
            output_filetype="xyz",
            include_intermediate_structures=True,
        )
        mock_fc.return_value.convert_files.assert_called_once_with()

    def test_cli_convert_command_object_defaults_to_xyz(
        self, single_model_pdb_file
    ):
        """Invoke convert directly, without going through ``run``."""
        output_path = single_model_pdb_file.replace(".pdb", ".xyz")
        result = CliRunner().invoke(
            convert, ["--input", single_model_pdb_file]
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)


class TestOpenBabelWriteFallback:
    """Tests for Open Babel fallback in ``Molecule.write``."""

    def test_write_mol2_via_openbabel(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")
        from openbabel import pybel

        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")
        mol.write(output_path, format="mol2")

        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0
        ob_mol = next(pybel.readfile("mol2", output_path), None)
        assert ob_mol is not None
        assert len(ob_mol.atoms) == mol.num_atoms

    def test_write_cml_via_openbabel(self, single_model_pdb_file, tmpdir):
        """Second non-native format: proves no OB format hardcoding."""
        pytest.importorskip("openbabel")
        from openbabel import pybel

        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.cml")
        mol.write(output_path, format="cml")

        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0
        ob_mol = next(pybel.readfile("cml", output_path), None)
        assert ob_mol is not None
        assert len(ob_mol.atoms) == mol.num_atoms

    def test_convert_file_to_mol2(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")
        from openbabel import pybel

        output_path = os.path.join(str(tmpdir), "water.mol2")
        FileConverter.convert_file(single_model_pdb_file, output_path)

        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0
        ob_mol = next(pybel.readfile("mol2", output_path), None)
        assert ob_mol is not None
        assert len(ob_mol.atoms) == 3

    def test_cli_convert_to_mol2(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")

        output_path = os.path.join(str(tmpdir), "water.mol2")
        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "--input",
                single_model_pdb_file,
                "--output",
                output_path,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0

    def test_native_formats_still_work(self, single_model_pdb_file, tmpdir):
        mol = Molecule.from_filepath(single_model_pdb_file)

        xyz_path = os.path.join(str(tmpdir), "water.xyz")
        mol.write(xyz_path, format="xyz")
        assert os.path.exists(xyz_path)
        assert Molecule.from_filepath(xyz_path).num_atoms == 3

        pdb_path = os.path.join(str(tmpdir), "water_out.pdb")
        mol.write(pdb_path, format="pdb")
        assert os.path.exists(pdb_path)
        assert Molecule.from_filepath(pdb_path).num_atoms == 3

        com_path = os.path.join(str(tmpdir), "water.com")
        mol.write(com_path, format="com")
        assert os.path.exists(com_path)
        assert Molecule.from_filepath(com_path).num_atoms == 3

    def test_missing_openbabel_raises_import_error(
        self, single_model_pdb_file, tmpdir
    ):
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        import builtins

        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "openbabel" or name.startswith("openbabel."):
                raise ImportError("No module named 'openbabel'")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_fake_import):
            with pytest.raises(ImportError, match="openbabel"):
                mol.write(output_path, format="mol2")

    def test_nonsense_format_raises_value_error(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.notarealformat")

        with pytest.raises(ValueError, match="notarealformat"):
            mol.write(output_path, format="notarealformat")

    def test_temp_xyz_removed_after_success(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        created_temps = []
        real_ntf = __import__("tempfile").NamedTemporaryFile

        def _tracking_ntf(*args, **kwargs):
            tmp = real_ntf(*args, **kwargs)
            created_temps.append(tmp.name)
            return tmp

        with patch("tempfile.NamedTemporaryFile", side_effect=_tracking_ntf):
            mol.write(output_path, format="mol2")

        assert created_temps
        for path in created_temps:
            assert not os.path.exists(path)

    def test_temp_xyz_removed_after_write_failure(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        created_temps = []
        real_ntf = __import__("tempfile").NamedTemporaryFile

        def _tracking_ntf(*args, **kwargs):
            tmp = real_ntf(*args, **kwargs)
            created_temps.append(tmp.name)
            return tmp

        with patch("tempfile.NamedTemporaryFile", side_effect=_tracking_ntf):
            with patch(
                "openbabel.pybel.Molecule.write",
                side_effect=RuntimeError("write boom"),
            ):
                with pytest.raises(ValueError, match="mol2"):
                    mol.write(output_path, format="mol2")

        assert created_temps
        for path in created_temps:
            assert not os.path.exists(path)

    def test_empty_openbabel_read_raises_and_cleans_up(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        created_temps = []
        real_ntf = __import__("tempfile").NamedTemporaryFile

        def _tracking_ntf(*args, **kwargs):
            tmp = real_ntf(*args, **kwargs)
            created_temps.append(tmp.name)
            return tmp

        with patch("tempfile.NamedTemporaryFile", side_effect=_tracking_ntf):
            with patch(
                "openbabel.pybel.readfile",
                return_value=iter([]),
            ):
                with pytest.raises(ValueError, match="Unable to read"):
                    mol.write(output_path, format="mol2")

        assert created_temps
        for path in created_temps:
            assert not os.path.exists(path)

    def test_write_pdb_pybabel_delegates(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water_ob.pdb")
        mol.write_pdb_pybabel(output_path)
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0

    def test_cleanup_false_keeps_temp_xyz(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        created_temps = []
        real_ntf = __import__("tempfile").NamedTemporaryFile

        def _tracking_ntf(*args, **kwargs):
            tmp = real_ntf(*args, **kwargs)
            created_temps.append(tmp.name)
            return tmp

        with patch("tempfile.NamedTemporaryFile", side_effect=_tracking_ntf):
            mol._write_via_openbabel(output_path, "mol2", cleanup=False)

        assert created_temps
        assert os.path.exists(created_temps[0])
        os.remove(created_temps[0])

    def test_native_extxyz_via_write(self, single_model_pdb_file, tmpdir):
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.extxyz")
        mol.write(output_path, format="extxyz")
        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0

    def test_missing_openbabel_cleanup_false_keeps_temp(
        self, single_model_pdb_file, tmpdir
    ):
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        created_temps = []
        real_ntf = __import__("tempfile").NamedTemporaryFile

        def _tracking_ntf(*args, **kwargs):
            tmp = real_ntf(*args, **kwargs)
            created_temps.append(tmp.name)
            return tmp

        import builtins

        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "openbabel" or name.startswith("openbabel."):
                raise ImportError("No module named 'openbabel'")
            return real_import(name, *args, **kwargs)

        with patch("tempfile.NamedTemporaryFile", side_effect=_tracking_ntf):
            with patch("builtins.__import__", side_effect=_fake_import):
                with pytest.raises(ImportError, match="openbabel"):
                    mol._write_via_openbabel(
                        output_path, "mol2", cleanup=False
                    )

        assert created_temps
        assert os.path.exists(created_temps[0])
        os.remove(created_temps[0])

    def test_cleanup_oserror_on_success_is_warned(
        self, single_model_pdb_file, tmpdir, caplog
    ):
        pytest.importorskip("openbabel")
        import logging

        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        real_remove = os.remove

        def _remove_then_fail(path):
            # Allow the first remove attempts from other code; fail on temp xyz
            if path.endswith(".xyz"):
                raise OSError("permission denied")
            return real_remove(path)

        with patch("os.remove", side_effect=_remove_then_fail):
            with caplog.at_level(logging.WARNING):
                mol.write(output_path, format="mol2")

        assert os.path.exists(output_path)
        assert any(
            "Failed to remove temporary file" in r.message
            for r in caplog.records
        )

    def test_cleanup_oserror_on_import_error(
        self, single_model_pdb_file, tmpdir
    ):
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "water.mol2")

        import builtins

        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "openbabel" or name.startswith("openbabel."):
                raise ImportError("No module named 'openbabel'")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_fake_import):
            with patch("os.remove", side_effect=OSError("busy")):
                with pytest.raises(ImportError, match="openbabel"):
                    mol.write(output_path, format="mol2")

    def test_write_cosmorsxyz_is_native(self, single_model_pdb_file, tmpdir):
        """cosmorsxyz uses the native writer (no Open Babel required)."""
        mol = Molecule.from_filepath(single_model_pdb_file)
        mol.charge = 0
        mol.multiplicity = 1
        output_path = os.path.join(str(tmpdir), "water.cosmorsxyz")

        with patch(
            "chemsmart.io.molecules.structure.Molecule._write_via_openbabel"
        ) as mock_ob:
            mol.write(output_path, format="cosmorsxyz")
            mock_ob.assert_not_called()

        assert os.path.exists(output_path)
        with open(output_path, "r") as f:
            lines = f.read().strip().splitlines()
        assert lines[0] == "3"
        assert lines[1] == "0 1"
        assert len(lines) == 5

    def test_xyz_to_pdb_none_delegates_to_write_pdb_pybabel(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        output_path = os.path.join(str(tmpdir), "from_xyz_to_pdb.pdb")

        with patch.object(
            mol, "write_pdb_pybabel", wraps=mol.write_pdb_pybabel
        ) as mock_delegate:
            FileConverter.xyz_to_pdb(mol, output_path, xyz_filename=None)
            mock_delegate.assert_called_once()

        assert os.path.exists(output_path)
        assert os.path.getsize(output_path) > 0

    def test_xyz_to_pdb_existing_file_branch(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        xyz_path = os.path.join(str(tmpdir), "water.xyz")
        mol.write_xyz(xyz_path, mode="w")
        pdb_path = os.path.join(str(tmpdir), "from_existing.pdb")

        FileConverter.xyz_to_pdb(mol, pdb_path, xyz_filename=xyz_path)

        assert os.path.exists(pdb_path)
        assert os.path.exists(xyz_path)  # caller-supplied path kept


class TestConverterCoverageGaps:
    """Hit remaining FileConverter branches for ≥90% file coverage."""

    def test_convert_files_neither_directory_nor_filename_raises(self):
        with pytest.raises(ValueError, match="Either directory or filename"):
            FileConverter().convert_files()

    def test_convert_out_without_program_raises(self, tmpdir):
        with pytest.raises(ValueError, match="--program"):
            FileConverter(
                directory=str(tmpdir), type="out", program=None
            ).convert_files()

    def test_convert_unsupported_batch_type_raises(self, tmpdir):
        with pytest.raises(ValueError, match="not supported"):
            FileConverter(directory=str(tmpdir), type="xyzzy").convert_files()

    def test_convert_unsupported_single_type_raises(self, tmpdir):
        path = os.path.join(str(tmpdir), "x.xyzzy")
        with open(path, "w") as f:
            f.write("dummy\n")
        fc = FileConverter(filename=path, output_filetype="xyz")
        # After the unified convert_file path, unsupported extensions fall
        # through Open Babel then ASE; ASE raises UnknownFileTypeError.
        with pytest.raises(Exception):
            fc.convert_files()

    def test_batch_gjf_to_xyz(self, tmpdir, gaussian_opt_inputfile):
        from shutil import copy

        gjf = os.path.join(str(tmpdir), "opt.gjf")
        copy(gaussian_opt_inputfile, gjf)
        FileConverter(
            directory=str(tmpdir), type="gjf", output_filetype="xyz"
        ).convert_files()
        assert os.path.exists(gjf.replace(".gjf", ".xyz"))

    def test_batch_pdb_to_xyz(self, single_model_pdb_file):
        FileConverter(
            directory=os.path.dirname(single_model_pdb_file),
            type="pdb",
            output_filetype="xyz",
        ).convert_files()
        assert os.path.exists(single_model_pdb_file.replace(".pdb", ".xyz"))

    def test_batch_sdf_to_xyz(self, tmpdir):
        # Minimal one-molecule SDF
        sdf_path = os.path.join(str(tmpdir), "methane.sdf")
        with open(sdf_path, "w") as f:
            f.write(
                "methane\n"
                "  ChemSmart\n"
                "\n"
                "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "M  END\n"
                "$$$$\n"
            )
        FileConverter(
            directory=str(tmpdir), type="sdf", output_filetype="xyz"
        ).convert_files()
        assert os.path.exists(os.path.join(str(tmpdir), "methane.xyz"))

    def test_batch_orca_out_to_xyz(self, tmpdir, orca_co2_output):
        from shutil import copy

        copy(orca_co2_output, os.path.join(str(tmpdir), "CO2.out"))
        FileConverter(
            directory=str(tmpdir),
            type="out",
            program="orca",
            output_filetype="xyz",
        ).convert_files()
        assert os.path.exists(os.path.join(str(tmpdir), "CO2.xyz"))

    def test_batch_gaussian_out_to_xyz(
        self, tmpdir, gaussian_ozone_opt_outfile
    ):
        from shutil import copy

        # Gaussian outputs are usually .log; exercise type=out + program=gaussian
        dest = os.path.join(str(tmpdir), "ozone.out")
        copy(gaussian_ozone_opt_outfile, dest)
        FileConverter(
            directory=str(tmpdir),
            type="out",
            program="gaussian",
            output_filetype="xyz",
        ).convert_files()
        assert os.path.exists(os.path.join(str(tmpdir), "ozone.xyz"))

    def test_batch_orca_inp_to_xyz(self, tmpdir, water_opt_input_path):
        from shutil import copy

        copy(water_opt_input_path, os.path.join(str(tmpdir), "water_opt.inp"))
        FileConverter(
            directory=str(tmpdir), type="inp", output_filetype="xyz"
        ).convert_files()
        assert os.path.exists(os.path.join(str(tmpdir), "water_opt.xyz"))

    def test_batch_include_intermediate_list_write(
        self, tmpdir, gaussian_ozone_opt_outfile
    ):
        from shutil import copy

        copy(
            gaussian_ozone_opt_outfile, os.path.join(str(tmpdir), "ozone.log")
        )
        FileConverter(
            directory=str(tmpdir),
            type="log",
            output_filetype="xyz",
            include_intermediate_structures=True,
        ).convert_files()
        # Multi-structure log → numbered files, not a single overwrite path.
        assert not os.path.exists(os.path.join(str(tmpdir), "ozone.xyz"))
        assert os.path.exists(os.path.join(str(tmpdir), "ozone_1.xyz"))
        assert os.path.exists(os.path.join(str(tmpdir), "ozone_7.xyz"))
        mol = Molecule.from_filepath(os.path.join(str(tmpdir), "ozone_7.xyz"))
        assert mol.num_atoms == 3
        assert mol.chemical_formula == "O3"

    def test_single_orca_out_to_xyz(self, tmpdir, orca_co2_output):
        from shutil import copy

        path = os.path.join(str(tmpdir), "CO2.out")
        copy(orca_co2_output, path)
        FileConverter(filename=path, output_filetype="xyz").convert_files()
        assert os.path.exists(path.replace(".out", ".xyz"))

    def test_single_inp_to_xyz(self, tmpdir, water_opt_input_path):
        from shutil import copy

        path = os.path.join(str(tmpdir), "water_opt.inp")
        copy(water_opt_input_path, path)
        FileConverter(filename=path, output_filetype="xyz").convert_files()
        assert os.path.exists(path.replace(".inp", ".xyz"))

    def test_single_sdf_to_xyz(self, tmpdir):
        path = os.path.join(str(tmpdir), "c.sdf")
        with open(path, "w") as f:
            f.write(
                "c\n\n\n"
                "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                "M  END\n$$$$\n"
            )
        FileConverter(filename=path, output_filetype="xyz").convert_files()
        assert os.path.exists(path.replace(".sdf", ".xyz"))

    def test_single_pdb_to_xyz(self, single_model_pdb_file):
        FileConverter(
            filename=single_model_pdb_file, output_filetype="xyz"
        ).convert_files()
        assert os.path.exists(single_model_pdb_file.replace(".pdb", ".xyz"))

    def test_single_gjf_to_xyz(self, tmpdir, gaussian_opt_inputfile):
        from shutil import copy

        path = os.path.join(str(tmpdir), "opt.gjf")
        copy(gaussian_opt_inputfile, path)
        FileConverter(filename=path, output_filetype="xyz").convert_files()
        assert os.path.exists(path.replace(".gjf", ".xyz"))

    def test_single_unknown_program_out_raises(self, tmpdir):
        path = os.path.join(str(tmpdir), "mystery.out")
        with open(path, "w") as f:
            f.write("not a real quantum chemistry output\n")
        fc = FileConverter(filename=path, output_filetype="xyz")
        with pytest.raises(
            ValueError, match="Unsupported .out file program type"
        ):
            fc.convert_files()

    def test_single_include_intermediate_list(
        self, tmpdir, gaussian_ozone_opt_outfile
    ):
        from shutil import copy

        path = os.path.join(str(tmpdir), "ozone.log")
        copy(gaussian_ozone_opt_outfile, path)
        FileConverter(
            filename=path,
            output_filetype="xyz",
            include_intermediate_structures=True,
        ).convert_files()
        # Multi-structure log → numbered files, not a single overwrite path.
        assert not os.path.exists(path.replace(".log", ".xyz"))
        assert os.path.exists(os.path.join(str(tmpdir), "ozone_1.xyz"))
        assert os.path.exists(os.path.join(str(tmpdir), "ozone_7.xyz"))
        mol = Molecule.from_filepath(os.path.join(str(tmpdir), "ozone_7.xyz"))
        assert mol.num_atoms == 3
        assert mol.chemical_formula == "O3"

    def test_xyz_to_pdb_missing_path_writes_then_converts(
        self, single_model_pdb_file, tmpdir
    ):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        xyz_path = os.path.join(str(tmpdir), "missing.xyz")
        pdb_path = os.path.join(str(tmpdir), "out.pdb")
        assert not os.path.exists(xyz_path)
        FileConverter.xyz_to_pdb(mol, pdb_path, xyz_filename=xyz_path)
        assert os.path.exists(xyz_path)
        assert os.path.exists(pdb_path)

    def test_xyz_to_pdb_empty_read_raises(self, single_model_pdb_file, tmpdir):
        pytest.importorskip("openbabel")
        mol = Molecule.from_filepath(single_model_pdb_file)
        xyz_path = os.path.join(str(tmpdir), "emptyish.xyz")
        mol.write_xyz(xyz_path, mode="w")
        pdb_path = os.path.join(str(tmpdir), "fail.pdb")
        with patch("openbabel.pybel.readfile", return_value=iter([])):
            with pytest.raises(ValueError, match="Unable to read"):
                FileConverter.xyz_to_pdb(mol, pdb_path, xyz_filename=xyz_path)


class TestOpenBabelReadFallback:
    """Tests for Open Babel-first input fallback in ``Molecule._read_other``."""

    def _write_water_mol2(self, path):
        pytest.importorskip("openbabel")
        from openbabel import pybel

        mol = pybel.readstring(
            "xyz",
            "3\nwater\n"
            "O  0.000  0.000  0.000\n"
            "H  0.960  0.000  0.000\n"
            "H -0.240  0.930  0.000\n",
        )
        mol.write("mol2", path, overwrite=True)

    def test_read_mol2_via_openbabel(self, tmpdir):
        pytest.importorskip("openbabel")
        mol2_path = os.path.join(str(tmpdir), "water.mol2")
        self._write_water_mol2(mol2_path)

        mol = Molecule.from_filepath(mol2_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_convert_file_mol2_to_xyz(self, tmpdir):
        pytest.importorskip("openbabel")
        mol2_path = os.path.join(str(tmpdir), "water.mol2")
        self._write_water_mol2(mol2_path)
        xyz_path = os.path.join(str(tmpdir), "water_from_mol2.xyz")

        FileConverter.convert_file(mol2_path, xyz_path)

        assert os.path.exists(xyz_path)
        mol = Molecule.from_filepath(xyz_path)
        assert mol.num_atoms == 3
        assert mol.chemical_formula == "H2O"

    def test_cli_convert_mol2_to_xyz(self, tmpdir):
        pytest.importorskip("openbabel")
        mol2_path = os.path.join(str(tmpdir), "water.mol2")
        self._write_water_mol2(mol2_path)
        xyz_path = os.path.join(str(tmpdir), "water_cli.xyz")

        runner = CliRunner()
        result = runner.invoke(
            entry_point,
            [
                "run",
                "convert",
                "-i",
                mol2_path,
                "-o",
                xyz_path,
            ],
        )
        assert result.exit_code == 0, result.output
        assert os.path.exists(xyz_path)
        assert Molecule.from_filepath(xyz_path).num_atoms == 3

    def test_read_smi_via_openbabel_generates_3d(self, tmpdir):
        pytest.importorskip("openbabel")
        smi_path = os.path.join(str(tmpdir), "benzene.smi")
        with open(smi_path, "w") as f:
            f.write("c1ccccc1 benzene\n")

        mol = Molecule.from_filepath(smi_path)
        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 12  # C6H6 after make3D + implicit H
        assert mol.chemical_formula == "C6H6"
        # Coordinates should not all be zero after make3D.
        assert np.linalg.norm(mol.positions) > 0.1

    def test_ase_only_cif_skips_openbabel(self, tmpdir):
        """Periodic .cif must go straight to ASE (skip-list)."""
        from ase.build import bulk
        from ase.io import write

        cif_path = os.path.join(str(tmpdir), "nacl.cif")
        write(cif_path, bulk("NaCl", "rocksalt", a=5.64))

        with patch.object(
            Molecule, "_read_via_openbabel", side_effect=AssertionError
        ) as mock_ob:
            mol = Molecule.from_filepath(cif_path)
            mock_ob.assert_not_called()

        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 2
        assert mol.pbc_conditions is not None

    def test_ase_only_traj_skips_openbabel(self, tmpdir):
        from ase import Atoms
        from ase.io import write

        traj_path = os.path.join(str(tmpdir), "water.traj")
        write(
            traj_path,
            Atoms(
                "OHH",
                positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
            ),
        )

        with patch.object(
            Molecule, "_read_via_openbabel", side_effect=AssertionError
        ) as mock_ob:
            mol = Molecule.from_filepath(traj_path)
            mock_ob.assert_not_called()

        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_ase_only_cfg_skips_openbabel(self, tmpdir):
        """Periodic / ASE .cfg must go straight to ASE (skip-list)."""
        from ase import Atoms
        from ase.io import write

        cfg_path = os.path.join(str(tmpdir), "water.cfg")
        write(
            cfg_path,
            Atoms(
                "OHH",
                positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
            ),
        )

        with patch.object(
            Molecule, "_read_via_openbabel", side_effect=AssertionError
        ) as mock_ob:
            mol = Molecule.from_filepath(cfg_path)
            mock_ob.assert_not_called()

        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_obabel_failure_falls_back_to_ase(self, tmpdir):
        """When Open Babel fails, ASE is tried for non-skip-list formats."""
        from ase import Atoms
        from ase.io import write

        # ASE JSON is not in ASE_ONLY_EXTENSIONS; Open Babel typically cannot
        # read it, so a forced OB failure should fall through to ASE.
        json_path = os.path.join(str(tmpdir), "water.json")
        write(
            json_path,
            Atoms(
                "OHH",
                positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
            ),
            format="json",
        )

        with patch.object(
            Molecule,
            "_read_via_openbabel",
            side_effect=ValueError("forced openbabel failure"),
        ):
            mol = Molecule.from_filepath(json_path)

        assert isinstance(mol, Molecule)
        assert mol.num_atoms == 3
        assert list(mol.symbols) == ["O", "H", "H"]

    def test_missing_openbabel_falls_back_to_ase(self, tmpdir):
        from ase import Atoms
        from ase.io import write

        json_path = os.path.join(str(tmpdir), "water.json")
        write(
            json_path,
            Atoms(
                "OHH",
                positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
            ),
            format="json",
        )

        import builtins

        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "openbabel" or name.startswith("openbabel."):
                raise ImportError("No module named 'openbabel'")
            return real_import(name, *args, **kwargs)

        with patch("builtins.__import__", side_effect=_fake_import):
            mol = Molecule.from_filepath(json_path)

        assert mol.num_atoms == 3

    def test_multi_molecule_mol2_index_and_return_list(self, tmpdir):
        pytest.importorskip("openbabel")
        from openbabel import pybel

        mol2_path = os.path.join(str(tmpdir), "two.mol2")
        m1 = pybel.readstring(
            "xyz",
            "1\nmethane-C\nC  0.0  0.0  0.0\n",
        )
        m2 = pybel.readstring(
            "xyz",
            "2\nnitrogen\nN  0.0  0.0  0.0\nN  1.1  0.0  0.0\n",
        )
        m1.write("mol2", mol2_path, overwrite=True)
        # Open Babel refuses overwrite=False on an existing path; append text.
        with open(mol2_path, "a") as f:
            f.write(m2.write("mol2"))

        all_mols = Molecule.from_filepath(
            mol2_path, index=":", return_list=True
        )
        assert isinstance(all_mols, list)
        assert len(all_mols) == 2
        assert all_mols[0].num_atoms == 1
        assert all_mols[1].num_atoms == 2

        first = Molecule.from_filepath(mol2_path, index="1")
        assert first.num_atoms == 1
        last = Molecule.from_filepath(mol2_path, index="-1")
        assert last.num_atoms == 2

    def test_ase_only_extensions_constant(self):
        from chemsmart.io.molecules.structure import ASE_ONLY_EXTENSIONS

        for ext in ("traj", "db", "cif", "cfg", "vasp", "extxyz"):
            assert ext in ASE_ONLY_EXTENSIONS
        assert "mol2" not in ASE_ONLY_EXTENSIONS
        assert "smi" not in ASE_ONLY_EXTENSIONS

    def test_read_charged_smiles_propagates_charge(self, tmpdir):
        pytest.importorskip("openbabel")
        smi_path = os.path.join(str(tmpdir), "acetate.smi")
        with open(smi_path, "w") as f:
            f.write("CC(=O)[O-] acetate\n")

        mol = Molecule.from_filepath(smi_path)
        assert mol.charge == -1
        assert mol.num_atoms >= 7  # C2H3O2 after implicit H

    def test_read_neutral_mol2_charge_is_none(self, tmpdir):
        pytest.importorskip("openbabel")
        mol2_path = os.path.join(str(tmpdir), "water.mol2")
        self._write_water_mol2(mol2_path)

        mol = Molecule.from_filepath(mol2_path)
        assert mol.charge is None

    def test_convert_cif_to_xyz_drops_pbc(self, tmpdir):
        from ase.build import bulk
        from ase.io import write

        cif_path = os.path.join(str(tmpdir), "nacl.cif")
        xyz_path = os.path.join(str(tmpdir), "nacl.xyz")
        write(cif_path, bulk("NaCl", "rocksalt", a=5.64))

        mol_cif = Molecule.from_filepath(cif_path)
        assert mol_cif.pbc_conditions is not None

        FileConverter.convert_file(cif_path, xyz_path)
        mol_xyz = Molecule.from_filepath(xyz_path)
        assert mol_xyz.num_atoms == mol_cif.num_atoms
        # Standard XYZ write does not persist cell / PBC.
        assert mol_xyz.pbc_conditions is None or not any(
            mol_xyz.pbc_conditions
        )

    def test_numbered_output_warns_on_collision(
        self, tmpdir, multi_model_pdb_file, caplog
    ):
        import logging

        tmp_path = os.path.join(str(tmpdir), "ensemble.pdb")
        copy(multi_model_pdb_file, tmp_path)
        preexisting = os.path.join(str(tmpdir), "ensemble_1.xyz")
        with open(preexisting, "w") as f:
            f.write("1\nold\nHe 0 0 0\n")

        with caplog.at_level(logging.WARNING):
            FileConverter.convert_file(
                tmp_path,
                os.path.join(str(tmpdir), "ensemble.xyz"),
                include_intermediate_structures=True,
            )

        assert any(
            "Overwriting existing numbered output file" in r.message
            for r in caplog.records
        )
        # Converted content replaced the stub.
        mol = Molecule.from_filepath(preexisting)
        assert mol.chemical_formula != "He"

    def test_batch_xyz_folder_excludes_extxyz(self, tmpdir):
        from chemsmart.io.xyz.folder import XYZFolder

        xyz_path = os.path.join(str(tmpdir), "water.xyz")
        extxyz_path = os.path.join(str(tmpdir), "crystal.extxyz")
        with open(xyz_path, "w") as f:
            f.write(
                "3\nwater\n"
                "O  0.000  0.000  0.000\n"
                "H  0.960  0.000  0.000\n"
                "H -0.240  0.930  0.000\n"
            )
        with open(extxyz_path, "w") as f:
            f.write(
                '3\nLattice="1 0 0 0 1 0 0 0 1" Properties=species:S:1:pos:R:3\n'
                "O  0.000  0.000  0.000\n"
                "H  0.960  0.000  0.000\n"
                "H -0.240  0.930  0.000\n"
            )

        files = XYZFolder(folder=str(tmpdir)).all_xyzfiles
        assert xyz_path in files
        assert extxyz_path not in files

    def test_make3d_failure_logs_warning(self, tmpdir, caplog):
        pytest.importorskip("openbabel")
        import logging

        from openbabel import pybel

        smi_path = os.path.join(str(tmpdir), "benzene.smi")
        with open(smi_path, "w") as f:
            f.write("c1ccccc1 benzene\n")

        with patch.object(
            pybel.Molecule,
            "make3D",
            side_effect=RuntimeError("forced make3D failure"),
        ):
            with caplog.at_level(logging.WARNING):
                mol = Molecule.from_filepath(smi_path)

        assert isinstance(mol, Molecule)
        assert any("make3D failed" in r.message for r in caplog.records)
