import os
from filecmp import cmp
from shutil import copy

import pytest
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.gaussian import GaussianOptJob, GaussianSinglePointJob
from chemsmart.jobs.gaussian import GaussianModredJob
from chemsmart.jobs.gaussian import GaussianScanJob
from chemsmart.jobs.gaussian import GaussianTSJob
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


class TestGaussianInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file
        )  # writes input file as expected

    def test_write_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_no_scratch,
        gaussian_written_opt_file_with_route,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.route_to_be_written = "#p pbepbe/6-31g(d) opt"
        settings.title = "Optimisation job with supplied route"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_opt",
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_opt_file_with_route)

    def test_write_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_no_scratch,
        gaussian_written_modred_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.modred_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = [[1, 2], [3, 4, 5]]
        job = GaussianModredJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_modred",
        )
        assert isinstance(job, GaussianModredJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_modred.com")
        print(g16_file)
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_modred_file)

    def test_write_scan_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_no_scratch,
        gaussian_written_scan_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4, 5]],
            "num_steps": 10,
            "step_size": 0.1,
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan",
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_scan.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_scan_file)

    def test_write_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        jobrunner_no_scratch,
        gaussian_written_ts_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.ts_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianTSJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_ts",
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_ts_file)

    def test_write_opt_input_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        jobrunner_no_scratch,
        gaussian_written_ts_from_nhc_singlet_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log output
        and write aldehyde_opt.com using the settings from the .log file."""
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        # also merge the title keywords
        keywords = ("charge", "multiplicity", "title")
        ts_settings = ts_settings.merge(job_settings, keywords=keywords)
        job = GaussianTSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="gaussian_ts_from_log",
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_ts_from_nhc_singlet_log_file)

    def test_write_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
    ):
        """Test writing simple .com input file using settings from .log file,
        including solvation."""

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_solvent",
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_solvent.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
        )

    def test_it_writes_sp_with_custom_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        jobrunner_no_scratch,
        smd_TBME_solvent_parameters_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
    ):
        """Test writing input file from log file.
        Simply taking the Gaussian nhc_neutral_singlet.log output and write
        gaussian_sp_custom_solv.com using the settings from the .log
        file and including custom solvation parameters from file smd_TBME.
        """
        smd_TBME_tmp_path = os.path.join(tmpdir, "smd_TBME")
        copy(smd_TBME_solvent_parameters_text_file, smd_TBME_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        sp_settings.solvent_model = "smd"
        sp_settings.solvent_id = "generic,read"
        sp_settings.custom_solvent = smd_TBME_tmp_path

        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_custom_solv",
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(
            job=job, jobrunner=jobrunner_no_scratch
        )
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_sp_custom_solv.com")
        assert os.path.isfile(g16_file)
        print(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
        )

    def test_it_writes_single_point_from_logfile(
        self, tmpdir, normal_min_log_filepath
    ):
        """Test writing input file from log file.

        Taking the Gaussian aldehyde.log output and write sp with different functional/basis and with solvation
        aldehyde_sp.com using the settings from the .log file and updating those settings.
        """
        atoms = AtomsWrapper.from_file(normal_min_log_filepath)
        settings = GaussianJobSettings.from_logfile(normal_min_log_filepath)
        settings.functional = "b3lyp empiricaldispersion=gd3bj"
        settings.basis = "def2tzvp"
        settings.solvent_model = "smd"
        settings.solvent_id = "toluene"
        settings.job_type = "sp"
        settings.freq = False

        written_file = settings.write_gaussian_input(
            output_dir=tmpdir, job_label="input", atoms=atoms
        )
        written_settings = GaussianJobSettings.from_comfile(written_file)
        written_atoms = AtomsWrapper.from_file(written_file)
        assert written_settings == settings
        assert written_atoms == atoms

    def test_it_writes_single_point_with_diff_functional_and_basis_and_custom_solvent(
        self,
        tmpdir,
        smd_TBME_solvent_parameters_txt_file,
        normal_min_log_filepath,
    ):
        atoms = AtomsWrapper.from_file(normal_min_log_filepath)
        settings = GaussianJobSettings.from_logfile(normal_min_log_filepath)
        settings.functional = "b3lyp empiricaldispersion=gd3bj"
        settings.basis = "def2tzvp"
        settings.solvent_model = "smd"
        settings.solvent_id = "generic,read"
        settings.job_type = "sp"
        settings.freq = False

        # copy solvent parameters file:
        settings.set_custom_solvent_via_file(
            smd_TBME_solvent_parameters_txt_file
        )

        written_file = settings.write_gaussian_input(
            output_dir=tmpdir, job_label="input", atoms=atoms
        )
        written_settings = GaussianJobSettings.from_comfile(written_file)
        written_atoms = AtomsWrapper.from_file(written_file)

        # fails as written_settings.solvent_id is `generic,read`, whereas settings.solvent_id is `TBME`
        # <-- should be `generic,read`, as `TBME` is not in a list of solvents already parametrized by Gaussian
        assert written_settings == settings
        assert written_atoms == atoms

    def test_it_writes_modredundant_file_from_logfile(
        self, tmpdir, normal_ts_filepath
    ):
        """Test writing modredundant input file using settings from log file output."""
        atoms = AtomsWrapper.from_file(normal_ts_filepath)
        settings = GaussianJobSettings.from_logfile(normal_ts_filepath)
        settings.job_type = "modred"
        settings.modred = [[2, 12], [9, 2]]

        written_file = settings.write_gaussian_input(
            output_dir=tmpdir, job_label="input", atoms=atoms
        )
        written_settings = GaussianJobSettings.from_comfile(written_file)
        written_atoms = AtomsWrapper.from_file(written_file)
        assert written_settings == settings
        assert written_atoms == atoms

    def test_it_writes_modredundant_file_with_solvation_from_logfile(
        self, tmpdir, normal_ts_filepath, model_solv_modred_inputfile
    ):
        """Test writing input file from log file.

        Taking the Gaussian normal_ts.log output and write modredundant ts search file
        normal_ts_solv_modred.com using the settings from the .log file and with solvation.
        """
        gaussian_settings = GaussianJobSettings.from_logfile(
            normal_ts_filepath
        )
        gaussian_settings.job_type = "modred"
        gaussian_settings.modred = [[2, 12], [9, 2]]
        gaussian_settings.update_solvent(
            solvent_model="smd", solvent_id="toluene"
        )
        gaussian_settings.freq = False
        atoms = AtomsWrapper.from_filepath(normal_ts_filepath)

        gaussian_settings.write_gaussian_input(
            atoms=atoms,
            output_dir=tmpdir,
            job_label="normal_ts_solv_modred",
            num_cores=32,
            mem_gigs=30,
        )
        written_file = os.path.join(tmpdir, "normal_ts_solv_modred.com")
        written_atoms = AtomsWrapper.from_filepath(filepath=written_file)
        written_settings = GaussianJobSettings.from_comfile(written_file)
        assert written_settings == gaussian_settings
        assert written_atoms == atoms

        # model_file_tmp_path = os.path.join(tmpdir, 'model_solv_modred_input.com')
        # copy(model_solv_modred_inputfile, model_file_tmp_path)
        #
        # compare_result = compare_file(model_file_tmp_path, written_file)
        # assert compare_result is True

    def test_it_writes_scan_file_from_logfile(
        self, tmpdir, normal_ts_filepath
    ):
        """Test writing scan input file from log file.

        Taking the Gaussian normal_ts.log output and write scan file
        normal_ts_scan.com using the settings from the .log file.
        """
        atoms = AtomsWrapper.from_file(normal_ts_filepath)
        settings = GaussianJobSettings.from_logfile(normal_ts_filepath)
        settings.job_type = "scan"
        settings.modred = {
            "num_steps": 10,
            "step_size": 0.05,
            "coords": [[2, 12], [9, 2]],
        }
        written_file = settings.write_gaussian_input(
            output_dir=tmpdir, job_label="scan", atoms=atoms
        )

        written_atoms = AtomsWrapper.from_filepath(written_file)
        written_settings = GaussianJobSettings.from_comfile(written_file)

        assert written_settings == settings
        assert written_atoms == atoms

    def test_it_writes_opt_with_custom_basis_from_logfile(
        self,
        tmpdir,
        normal_ts_filepath,
        Ni_def2tzvp_PCHOSi_svp_txt_file,
        model_custom_basis_opt_input,
    ):
        # supports writing custom basis from file path
        atoms = AtomsWrapper.from_filepath(
            filepath=model_custom_basis_opt_input
        )
        gaussian_settings = GaussianJobSettings.from_logfile(
            normal_ts_filepath
        )

        # update functional and basis
        gaussian_settings.functional = "b3lyp empiricaldispersion=gd3bj"
        gaussian_settings.basis = "gen"

        # copy gen basis file:
        gen_basis_tmp_path = os.path.join(tmpdir, "Ni_def2tzvp_PCHOSi_svp.txt")
        copy(Ni_def2tzvp_PCHOSi_svp_txt_file, gen_basis_tmp_path)

        gaussian_settings.gen_genecp = gen_basis_tmp_path

        gaussian_settings.job_type = "opt"

        written_file = gaussian_settings.write_gaussian_input(
            atoms=atoms,
            output_dir=tmpdir,
            job_label="normal_ts_custom_basis_opt",
        )
        written_settings = GaussianJobSettings.from_filepath(
            filepath=written_file
        )
        written_atoms = AtomsWrapper.from_filepath(filepath=written_file)

        model_input_tmpdir = os.path.join(
            tmpdir, "model_custom_basis_opt_input.com"
        )
        copy(model_custom_basis_opt_input, model_input_tmpdir)
        model_settings = GaussianJobSettings.from_filepath(
            filepath=model_input_tmpdir
        )
        model_atoms = AtomsWrapper.from_filepath(filepath=model_input_tmpdir)

        assert written_settings == model_settings
        assert written_atoms == atoms == model_atoms

    @pytest.mark.slow
    def test_it_writes_opt_with_custom_basis_using_api(
        self,
        tmpdir,
        genecp_log_filepath,
        model_custom_basis_from_api_opt_input,
    ):
        # supports writing custom basis from file path
        atoms = AtomsWrapper.from_filepath(filepath=genecp_log_filepath)
        gaussian_settings = GaussianJobSettings.from_logfile(
            genecp_log_filepath
        )

        # update functional and basis
        gaussian_settings.functional = "b3lyp empiricaldispersion=gd3bj"
        gaussian_settings.job_type = "opt"
        gaussian_settings.basis = "genecp"
        gaussian_settings.heavy_elements = ["Pd"]
        gaussian_settings.heavy_elements_basis = "def2-TZVPPD"
        gaussian_settings.light_elements_basis = "def2-SVP"
        written_file = gaussian_settings.write_gaussian_input(
            atoms=atoms, output_dir=tmpdir, job_label="custom_basis"
        )
        written_settings = GaussianJobSettings.from_filepath(
            filepath=written_file
        )
        written_atoms = AtomsWrapper.from_filepath(filepath=written_file)

        model_custom_basis_from_api_opt_input_tmpdir = os.path.join(
            tmpdir, "model_custom_basis_from_api_opt_input.com"
        )
        copy(
            model_custom_basis_from_api_opt_input,
            model_custom_basis_from_api_opt_input_tmpdir,
        )
        model_settings = GaussianJobSettings.from_filepath(
            model_custom_basis_from_api_opt_input_tmpdir
        )
        model_atoms = AtomsWrapper.from_filepath(
            filepath=model_custom_basis_from_api_opt_input_tmpdir
        )

        assert written_settings == model_settings
        assert written_atoms == atoms == model_atoms

    @pytest.mark.slow
    def test_it_writes_opt_with_custom_basis_for_all_elements_in_structure_using_api(
        self, tmpdir, com_filepath, model_custom_basis_for_all_elements_input
    ):
        # supports writing custom basis from file path
        atoms = AtomsWrapper.from_filepath(filepath=com_filepath)
        gaussian_settings = GaussianJobSettings.from_comfile(com_filepath)

        # update basis settings (or given such basis set settings)
        gaussian_settings.basis = "gen"
        gaussian_settings.heavy_elements = [
            "C",
            "H",
            "O",
            "N",
            "Cl",
            "P",
            "S",
            "F",
            "Br",
            "I",
        ]
        # more than all elements in the system but will be filtered to only those in the system for input preparation
        gaussian_settings.heavy_elements_basis = "def2-TZVPPD"
        gaussian_settings.light_elements_basis = (
            None  # light element basis not specified as
        )
        # all use custom basis from heavy_elements_basis

        written_file = gaussian_settings.write_gaussian_input(
            atoms=atoms,
            output_dir=tmpdir,
            job_label="custom_basis_for_all_elements",
        )
        written_settings = GaussianJobSettings.from_filepath(
            filepath=written_file
        )
        written_atoms = AtomsWrapper.from_filepath(filepath=written_file)

        # compare to the model file
        model_file_settings = GaussianJobSettings.from_filepath(
            filepath=model_custom_basis_for_all_elements_input
        )
        model_file_atoms = AtomsWrapper.from_filepath(
            filepath=model_custom_basis_for_all_elements_input
        )

        assert written_settings == model_file_settings
        assert written_atoms == atoms == model_file_atoms

        # also check that the correct basis name is in the input file
        with open(written_file) as f:
            lines = f.readlines()
            all_lines_string = "".join(lines)
            assert "def2-TZVPPD" in all_lines_string
            assert (
                "H     0" in all_lines_string
            )  # all these elements have been specified in custom basis
            assert "C     0" in all_lines_string
            assert "O     0" in all_lines_string
            assert "Cl     0" in all_lines_string
            assert (
                "N     0" not in all_lines_string
            )  # these are in the list of heavy elements but not in the system
            assert "P     0" not in all_lines_string
            assert "S     0" not in all_lines_string
            assert "F     0" not in all_lines_string
            assert "Br     0" not in all_lines_string
            assert "I     0" not in all_lines_string

        # can further check that the written file is the same as the model file
        assert cmp(written_file, model_custom_basis_for_all_elements_input)

    def test_writes_gaussian_input_from_pbc_logfile(
        self, tmpdir, graphite_sheet_logfile
    ):
        atoms = AtomsWrapper.from_filepath(graphite_sheet_logfile)

        settings = GaussianJobSettings.from_filepath(
            filepath=graphite_sheet_logfile
        )
        settings.functional = "pbepbe"
        settings.basis = "6-31g(d,p)/Auto"
        settings.job_type = "sp"
        written_file = settings.write_gaussian_input(
            output_dir=tmpdir, job_label="graphite_2d_opt_out", atoms=atoms
        )

        with open(written_file) as f:
            lines = f.readlines()
            assert lines[-4].split()[0] == "C", "Last element is C"
            assert lines[-3].split()[0] == "TV", "PBC should be written."
            assert lines[-2].split()[0] == "TV", "PBC should be written."
            assert len(lines[-1].strip()) == 0, "Last line is empty line."
