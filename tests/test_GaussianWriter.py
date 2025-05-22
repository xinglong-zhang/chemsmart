import os
from filecmp import cmp
from shutil import copy

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian import (
    GaussianModredJob,
    GaussianOptJob,
    GaussianQMMMJob,
    GaussianScanJob,
    GaussianSinglePointJob,
    GaussianTSJob,
)
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.utils.utils import cmp_with_ignore
from tests.conftest import (
    gaussian_yaml_settings_qmmm_project_name,
)


class TestGaussianInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file_with_route, shallow=False
        )

    def test_write_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_modred.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_modred_file, shallow=False)

    def test_write_scan_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_scan.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_scan_file, shallow=False)

    def test_write_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_ts_file, shallow=False)

    def test_write_qmmm_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_written_qmmm_file,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.real_charge = 0
        qmmm_settings.real_multiplicity = 1
        qmmm_settings.high_level_atoms = [1,2,3]
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm",
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(
            job=job
        )
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm.com")
        assert os.path.isfile(g16_file)
        assert cmp(g16_file, gaussian_written_qmmm_file, shallow=False)


    def test_write_qmmm_input_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_qmmm_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_qmmm_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log output and write qmmm .com"""
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.real_charge = 0
        qmmm_settings.real_multiplicity = 1
        qmmm_settings.high_level_atoms = [3, 12, 14, 7, 9]
        qmmm_settings.medium_level_atoms = [8, 17, 19 - 25]
        qmmm_settings.bonded_atoms = [(1, 3)]
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        keywords = ("charge", "multiplicity", "title")
        qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)
        job = GaussianQMMMJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=qmmm_settings,
            label="gaussian_qmmm_from_log",
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(
            job=job
        )
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_qmmm_log_file,
            shallow=False,
        )

    def test_write_opt_input_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_ts_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_ts_from_nhc_singlet_log_file,
            shallow=False,
        )

    def test_write_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_solvent.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_solvent_file,
            shallow=False,
        )

    def test_write_sp_with_custom_solvation_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
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
        sp_settings.solvent_model = "smd"
        sp_settings.custom_solvent = smd_TBME_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)

        job = GaussianSinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="gaussian_sp_from_log_with_custom_solvent",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_solvent.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_solvent_file,
            shallow=False,
        )

    def test_write_ts_with_custom_basis_from_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        gaussian_jobrunner_no_scratch,
        Ni_def2tzvp_PCHOSi_svp_text_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file,
    ):
        custom_basis_tmp_path = os.path.join(tmpdir, "custom_basis.txt")
        copy(Ni_def2tzvp_PCHOSi_svp_text_file, custom_basis_tmp_path)

        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        ts_settings.gen_genecp_file = custom_basis_tmp_path
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        job = GaussianTSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_file,
            shallow=False,
        )

    def test_write_ts_with_custom_basis_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_ts_genecp_outfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = GaussianJobSettings.from_logfile(
            gaussian_ts_genecp_outfile
        )
        ts_settings = ts_settings.merge(job_settings)

        # update settings for heavy elements
        ts_settings.heavy_elements = ["Pd"]
        ts_settings.heavy_elements_basis = "def2-TZVPPD"
        ts_settings.light_elements_basis = "def2-SVP"

        molecule = Molecule.from_filepath(gaussian_ts_genecp_outfile)

        job = GaussianTSJob(
            molecule=molecule,
            settings=ts_settings,
            label="gaussian_sp_from_log_with_custom_basis_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianTSJob)
        g16_writer = GaussianInputWriter(job=job)
        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_sp_from_log_with_custom_basis_from_api.com"
        )

        # compare the written input file with the expected input file, except
        # the line containing the version number (basis set exchange api may be different)
        assert cmp_with_ignore(
            g16_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
            ignore_string="Version",
        )

    def test_write_modred_with_custom_basis_for_all_elements_in_structure_using_api(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        modred_genecp_inputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_modred_with_custom_basis_for_all_atoms_from_api,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        modred_settings = project_settings.modred_settings()
        job_settings = GaussianJobSettings.from_filepath(
            modred_genecp_inputfile
        )
        modred_settings = modred_settings.merge(job_settings)

        modred_settings.modred = [[1, 2], [3, 4, 5]]
        modred_settings.basis = "genecp"
        modred_settings.heavy_elements = [
            "C",
            "H",
            "O",
            "N",
            "Pd",
            "P",
            "S",
            "F",
            "Br",
            "I",
        ]
        # more than all elements in the system but will be filtered to only those in the system for input preparation
        modred_settings.heavy_elements_basis = "def2-TZVPPD"
        modred_settings.light_elements_basis = None  # light element basis not specified as all use custom basis from heavy_elements_basis

        molecule = Molecule.from_filepath(modred_genecp_inputfile)
        job = GaussianModredJob(
            molecule=molecule,
            settings=modred_settings,
            label="gaussian_modred_with_custom_basis_for_all_atoms_from_api",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianModredJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir,
            "gaussian_modred_with_custom_basis_for_all_atoms_from_api.com",
        )
        assert os.path.isfile(g16_file)

        # compare the written input file with the expected input file, except
        # the line containing the version number (basis set exchange api may be different)
        assert cmp_with_ignore(
            g16_file,
            gaussian_modred_with_custom_basis_for_all_atoms_from_api,
            ignore_string="Version",
        )

    def test_write_gaussian_input_from_pbc_logfile(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_pbc_2d_outputfile,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_from_graphite_2d_pbc_log,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        file_settings = GaussianJobSettings.from_filepath(
            filepath=gaussian_pbc_2d_outputfile
        )
        opt_settings = opt_settings.merge(file_settings)

        molecule = Molecule.from_filepath(gaussian_pbc_2d_outputfile)
        g16_outputfile = Gaussian16Output(gaussian_pbc_2d_outputfile)
        assert g16_outputfile.list_of_pbc_conditions == [1, 1, 0]
        assert g16_outputfile.input_translation_vectors == [
            [2.47532, 0.0, 0.0],
            [-1.21995, 2.13345, 0.0],
        ]

        job = GaussianSinglePointJob(
            molecule=molecule,
            settings=opt_settings,
            label="graphite_2d_opt_from_log",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianSinglePointJob)

        g16_writer = GaussianInputWriter(job=job)
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "graphite_2d_opt_from_log.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_opt_from_graphite_2d_pbc_log,
            shallow=False,
        )
