import os.path
from shutil import copy, copytree, rmtree

from chemsmart.io.converter import FileConverter
from chemsmart.io.gaussian.folder import GaussianComFolder, GaussianLogFolder
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xyz.folder import XYZFolder


class TestConverter:

    def test_convert_log_foler_to_xyz(
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
        g16_folder = GaussianLogFolder(folder=tmp_log_folder)
        all_logfiles = g16_folder.all_logfiles

        # check all .log files have been converted to .xyz files
        for file in all_logfiles:
            assert os.path.exists(file.replace(".log", ".xyz"))

        ozone_xyz = os.path.join(tmp_log_folder, "ozone.xyz")

        assert os.path.exists(ozone_xyz)
        with open(ozone_xyz, "r") as f:
            lines = f.readlines()
            assert len(lines) == 5  # 5 lines in the log file
            assert lines[0] == "3\n"  # first line is number of atoms

    def test_convert_log_foler_to_com(
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
        g16_folder = GaussianLogFolder(folder=tmp_log_folder)
        all_logfiles = g16_folder.all_logfiles
        for file in all_logfiles:
            assert os.path.exists(file.replace(".log", ".com"))

        ozone_com = os.path.join(tmp_log_folder, "ozone.com")
        assert os.path.exists(ozone_com)
        with open(ozone_com, "r") as f:
            lines = f.readlines()
            assert len(lines) == 12
            assert lines[5].startswith("Generated from")

    def test_convert_com_foler_to_xyz(
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
        g16_folder = GaussianComFolder(folder=tmp_com_folder)
        all_comfiles = g16_folder.all_com_files
        print(all_comfiles)
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
        # assert np.isclose(mol.mass, 609.128, rtol=1e-4)  # in thermo branch

    def test_convert_single_link__opt_logfile_to_com(
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
        tmp_path = os.path.join(tmpdir, "benzene_opt.log")
        copy(gaussian_acetone_opt_outfile, tmp_path)
        assert os.path.exists(tmp_path)
        file_converter = FileConverter(
            filename=tmp_path, output_filetype="xyz"
        )

        file_converter.convert_files()

        assert os.path.exists(tmp_path.replace(".log", ".xyz"))
        mol = Molecule.from_filepath(tmp_path.replace(".log", ".xyz"))
        assert isinstance(mol, Molecule)
        print(mol.num_atoms, mol.chemical_formula, mol.energy)
        assert mol.num_atoms == 10
        assert mol.chemical_formula == "C3H6O"
        assert mol.energy == -192.919416
