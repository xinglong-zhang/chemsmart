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

    # ------------------------------------------------------------------
    # CDXML / CDX conversion tests
    # ------------------------------------------------------------------

    def test_convert_single_cdxml_to_xyz(
        self, tmpdir, single_molecule_cdxml_file_methane, expected_methane_xyz
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

        ref = Molecule.from_filepath(expected_methane_xyz)
        assert mol.num_atoms == ref.num_atoms
        assert mol.chemical_formula == ref.chemical_formula
        assert mol.positions == pytest.approx(ref.positions, abs=1e-5)

    def test_convert_single_cdxml_to_com(
        self, tmpdir, single_molecule_cdxml_file_benzene, expected_benzene_com
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

        ref = Molecule.from_filepath(expected_benzene_com)
        assert mol.num_atoms == ref.num_atoms
        assert mol.chemical_formula == ref.chemical_formula
        assert mol.positions == pytest.approx(ref.positions, abs=1e-5)

    def test_convert_single_cdx_to_xyz(
        self,
        tmpdir,
        single_molecule_cdx_file_imidazole,
        expected_imidazole_xyz,
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

        ref = Molecule.from_filepath(expected_imidazole_xyz)
        assert mol.num_atoms == ref.num_atoms
        assert mol.chemical_formula == ref.chemical_formula
        assert mol.positions == pytest.approx(ref.positions, abs=1e-5)

    def test_convert_multi_molecule_cdxml_to_xyz_splits_files(
        self,
        tmpdir,
        multi_molecule_cdxml_file,
        expected_two_molecules_1_xyz,
        expected_two_molecules_2_xyz,
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

        ref1 = Molecule.from_filepath(expected_two_molecules_1_xyz)
        assert mol1.num_atoms == ref1.num_atoms
        assert mol1.chemical_formula == ref1.chemical_formula
        assert mol1.positions == pytest.approx(ref1.positions, abs=1e-5)

        mol2 = Molecule.from_filepath(output_2)
        assert isinstance(mol2, Molecule)
        assert mol2.chemical_formula == "N2"
        assert mol2.num_atoms == 2

        ref2 = Molecule.from_filepath(expected_two_molecules_2_xyz)
        assert mol2.num_atoms == ref2.num_atoms
        assert mol2.chemical_formula == ref2.chemical_formula
        assert mol2.positions == pytest.approx(ref2.positions, abs=1e-5)

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

    def test_convert_cdxml_folder_to_com(
        self,
        tmpdir,
        chemdraw_directory,
        expected_benzene_com,
        expected_methane_com,
        expected_complex_molecule_com,
        expected_two_molecules_1_com,
        expected_two_molecules_2_com,
    ):
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

        # Compare contents against reference files
        for out_name, ref_path in (
            ("benzene.com", expected_benzene_com),
            ("methane.com", expected_methane_com),
            ("complex_molecule.com", expected_complex_molecule_com),
            ("two_molecules_1.com", expected_two_molecules_1_com),
            ("two_molecules_2.com", expected_two_molecules_2_com),
        ):
            mol = Molecule.from_filepath(
                os.path.join(tmp_cdxml_folder, out_name)
            )
            ref = Molecule.from_filepath(ref_path)
            assert mol.num_atoms == ref.num_atoms
            assert mol.chemical_formula == ref.chemical_formula
            assert mol.positions == pytest.approx(ref.positions, abs=1e-5)
