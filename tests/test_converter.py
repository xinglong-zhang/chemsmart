import os.path
from shutil import copy, copytree, rmtree

import numpy as np
import pytest

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
        tmp_path = os.path.join(tmpdir, "two_molecules.cdxml")
        copy(multi_molecule_cdxml_file, tmp_path)

        FileConverter(filename=tmp_path, output_filetype="xyz").convert_files()

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

    def test_convert_cdxml_folder_to_xyz(self, tmpdir, chemdraw_directory):
        from shutil import copytree

        tmp_cdxml_folder = os.path.join(tmpdir, "chemdraw")
        copytree(chemdraw_directory, tmp_cdxml_folder)

        FileConverter(
            directory=tmp_cdxml_folder, type="cdxml", output_filetype="xyz"
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
            directory=tmp_cdxml_folder, type="cdxml", output_filetype="com"
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
