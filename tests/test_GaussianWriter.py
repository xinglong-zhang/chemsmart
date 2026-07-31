import io
import os
from filecmp import cmp
from shutil import copy
from types import SimpleNamespace

import pytest

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

        # job run will result in the job being run and
        # the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_semiempirical_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_pm6_opt_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.semiempirical = "PM6"
        job = GaussianOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_pm6_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_pm6_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_pm6_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and
        # the output file copied back to run folder
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

    def test_write_scan_job_multiple_degrees_of_freedom(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_multiple_degrees_of_freedom_file,
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
            "num_steps": [10, 18],
            "step_size": [0.1, 5.0],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_multiple_degrees_of_freedom",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_scan_multiple_degrees_of_freedom.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_multiple_degrees_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_job_single_degree_of_freedom(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_single_degree_of_freedom_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [1, 2],
            "num_steps": [10],
            "step_size": [0.1],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_single_degree_of_freedom",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir, "gaussian_scan_single_degree_of_freedom.com"
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_single_degree_of_freedom_file,
            shallow=False,
        )

    def test_write_scan_job_multiple_degrees_of_freedom_with_constraints(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_scan_multiple_degrees_of_freedom_with_constraints_file,
    ):
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [
                [1, 2],
                [3, 4],
                [5, 6, 7],
                [8, 9, 10],
                [11, 12, 13, 14],
            ],
            "num_steps": [10, 12, 18, 15, 36],
            "step_size": [0.1, 0.2, 3.0, 4.0, 10.0],
        }
        job = GaussianScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_scan_multiple_degrees_of_freedom_with_constraints",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianScanJob)
        g16_writer = GaussianInputWriter(job=job)

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(
            tmpdir,
            "gaussian_scan_multiple_degrees_of_freedom_with_constraints.com",
        )
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file,
            gaussian_written_scan_multiple_degrees_of_freedom_with_constraints_file,
            shallow=False,
        )

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
        gaussian_jobrunner_no_scratch,
    ):
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.mult_total = 1
        qmmm_settings.high_level_atoms = [1, 2, 3]
        job = GaussianQMMMJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=qmmm_settings,
            label="gaussian_qmmm",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(job=job)
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_qmmm.com")
        assert qmmm_settings.high_level_atoms == [1, 2, 3]
        with open(g16_file, "r") as f:
            g16_content = f.read()
            # check that the qmmm section is present in the written file
        with open(gaussian_written_qmmm_file, "r") as f:
            expected_content = f.read()
        assert g16_content == expected_content
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
        """Taking the Gaussian nhc_neutral_singlet.log
        output and write qmmm .com"""
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_qmmm_project_name
        )
        qmmm_settings = project_settings.qmmm_settings()
        qmmm_settings.charge = 0
        qmmm_settings.multiplicity = 1
        qmmm_settings.charge_total = 0
        qmmm_settings.real_multiplicity = 1
        qmmm_settings.high_level_atoms = [3, 12, 14, 7, 9]
        qmmm_settings.medium_level_atoms = [8, 17, 19, 20, 21, 22, 23, 24, 25]
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
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianQMMMJob)
        g16_writer = GaussianInputWriter(job=job)
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
        # the line containing the version number
        # (basis set exchange api may be different)
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
        # more than all elements in the system but will be filtered
        # to only those in the system for input preparation
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
        # the line containing the version number
        # (basis set exchange api may be different)
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


def _make_writer(settings, molecule=None, label="testjob"):
    """Build a GaussianInputWriter around a lightweight fake job,
    bypassing full Job construction for direct private-method tests."""
    job = SimpleNamespace(
        settings=settings,
        jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
        label=label,
        molecule=molecule,
    )
    return GaussianInputWriter(job=job)


class TestGaussianInputWriterJobSpecificInfo:
    """Direct tests for _append_job_specific_info's per-jobtype
    branches (nci/wbi/resp/other), bypassing the full write()
    pipeline since only settings.jobtype and job.label are used."""

    def test_nci_jobtype_writes_wfn_file(self):
        settings = GaussianJobSettings.default()
        settings.jobtype = "nci"
        writer = _make_writer(settings, label="mymol")
        buf = io.StringIO()
        writer._append_job_specific_info(buf)
        assert buf.getvalue() == "mymol.wfn\n\n"

    def test_wbi_jobtype_writes_nbo_directive(self):
        settings = GaussianJobSettings.default()
        settings.jobtype = "wbi"
        writer = _make_writer(settings, label="mymol")
        buf = io.StringIO()
        writer._append_job_specific_info(buf)
        assert buf.getvalue() == "$nbo bndidx $end\n\n"

    def test_resp_jobtype_writes_gesp_file(self):
        settings = GaussianJobSettings.default()
        settings.jobtype = "resp"
        writer = _make_writer(settings, label="mymol")
        buf = io.StringIO()
        writer._append_job_specific_info(buf)
        assert buf.getvalue() == "mymol.gesp\n\n"

    def test_other_jobtype_writes_nothing(self):
        settings = GaussianJobSettings.default()
        settings.jobtype = "opt"
        writer = _make_writer(settings, label="mymol")
        buf = io.StringIO()
        writer._append_job_specific_info(buf)
        assert buf.getvalue() == ""


class TestGaussianInputWriterModredundantValidation:
    """Direct tests for _append_modredundant's format validation."""

    def test_invalid_modredundant_format_raises(self):
        settings = GaussianJobSettings.default()
        settings.modred = "not-a-list-or-dict"
        writer = _make_writer(settings)
        buf = io.StringIO()
        with pytest.raises(ValueError, match="modredundant must be"):
            writer._append_modredundant(buf)

    def test_none_modredundant_writes_nothing(self):
        settings = GaussianJobSettings.default()
        settings.modred = None
        writer = _make_writer(settings)
        buf = io.StringIO()
        writer._append_modredundant(buf)
        assert buf.getvalue() == ""


class TestGaussianInputWriterGenecpNoHeavyElements:
    """Direct test for _append_gen_genecp_basis's "no heavy elements
    and no user-specified gen_genecp_file" no-op branch."""

    def test_no_heavy_elements_and_no_gen_genecp_file_writes_nothing(
        self, single_molecule_xyz_file
    ):
        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        settings = GaussianJobSettings.default()
        settings.heavy_elements = ["Pd"]
        settings.heavy_elements_basis = "def2-tzvp"
        settings.light_elements_basis = "def2-svp"
        settings.gen_genecp_file = None
        # settings.genecp is a derived property that is truthy when a
        # basis-related gen/genecp keyword is used; force it directly
        # via the underlying flag if exposed, otherwise skip.
        settings.basis = "gen"
        writer = _make_writer(settings, molecule=molecule)
        buf = io.StringIO()
        writer._append_gen_genecp_basis(buf)
        assert buf.getvalue() == ""


class TestGaussianInputWriterWriteDispatch:
    """Direct tests for _write's target_directory/input_string
    branches."""

    def test_uses_job_folder_when_no_target_directory(self, tmp_path):
        settings = GaussianJobSettings.default()
        settings.charge = 0
        settings.multiplicity = 1
        settings.chk = False
        settings.functional = "b3lyp"
        settings.basis = "sto-3g"
        molecule = Molecule(
            symbols=["Ar"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        job = SimpleNamespace(
            settings=settings,
            jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
            label="nodir",
            molecule=molecule,
            folder=str(tmp_path),
        )
        writer = GaussianInputWriter(job=job)
        writer._write(target_directory=None)
        written = tmp_path / "nodir.com"
        assert written.is_file()
        content = written.read_text()
        assert "%chk=" not in content

    def test_write_self_used_when_input_string_set(self, tmp_path):
        settings = GaussianJobSettings.default()
        settings.input_string = "raw gaussian input content\n"
        job = SimpleNamespace(
            settings=settings,
            jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
            label="rawjob",
            molecule=None,
        )
        writer = GaussianInputWriter(job=job)
        writer._write(target_directory=str(tmp_path))
        written = tmp_path / "rawjob.com"
        assert written.read_text() == "raw gaussian input content\n"

    def test_creates_target_directory_if_missing(self, tmp_path):
        settings = GaussianJobSettings.default()
        settings.input_string = "content\n"
        job = SimpleNamespace(
            settings=settings,
            jobrunner=SimpleNamespace(num_cores=4, mem_gb=8),
            label="newdir",
            molecule=None,
        )
        writer = GaussianInputWriter(job=job)
        new_dir = tmp_path / "does_not_exist_yet"
        writer._write(target_directory=str(new_dir))
        assert (new_dir / "newdir.com").is_file()


class TestGaussianInputWriterModredundantScanConstrainedCoordinates:
    """Direct test for _append_modredundant's scan-job
    'constrained_coordinates' branch, in addition to the plain
    num_steps/step_size coords already covered elsewhere."""

    def test_scan_with_constrained_coordinates_writes_both_sections(self):
        settings = GaussianJobSettings.default()
        settings.modred = {
            "coords": [[1, 2]],
            "num_steps": [10],
            "step_size": [0.05],
            "constrained_coordinates": [[3, 4]],
        }
        writer = _make_writer(settings)
        buf = io.StringIO()
        writer._append_modredundant(buf)
        content = buf.getvalue()
        assert "B 1 2 S 10 0.05\n" in content
        assert "B 3 4 F\n" in content


class TestGaussianInputWriterGenecpAlwaysAppendsNewline:
    """Regression test for BUGS_FOUND.md #33: the
    `genecp_section.string_list[-1] != "\\n"` guard in
    _append_gen_genecp_basis is always true (str.split("\\n") never
    yields a literal "\\n" element), so a trailing blank line is
    always written after the genecp section regardless of whether
    the source content already ends with one."""

    def test_extra_blank_line_always_appended(self, tmp_path):
        genecp_file = tmp_path / "genecp.txt"
        genecp_file.write_text("C 0\nsto-3g\n****\n")
        settings = GaussianJobSettings.default()
        settings.basis = "genecp"
        settings.gen_genecp_file = str(genecp_file)
        settings.heavy_elements = None
        molecule = Molecule(
            symbols=["C"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        writer = _make_writer(settings, molecule=molecule)
        buf = io.StringIO()
        writer._append_gen_genecp_basis(buf)
        assert buf.getvalue().endswith("\n\n")


class TestGaussianInputWriterOtherAdditionalInfo:
    """Direct tests for _append_other_additional_info's file-path vs.
    free-string branches."""

    def test_reads_from_existing_file_path(self, tmp_path):
        info_file = tmp_path / "extra_info.txt"
        info_file.write_text("extra line one\nextra line two\n")
        settings = GaussianJobSettings.default()
        settings.append_additional_info = str(info_file)
        writer = _make_writer(settings)
        buf = io.StringIO()
        writer._append_other_additional_info(buf)
        assert buf.getvalue() == "extra line one\nextra line two\n"

    def test_writes_free_string_when_not_a_file_path(self):
        settings = GaussianJobSettings.default()
        settings.append_additional_info = "some free\nform text"
        writer = _make_writer(settings)
        buf = io.StringIO()
        writer._append_other_additional_info(buf)
        assert buf.getvalue() == "some free\nform text\n"


class TestGaussianInputWriterLinkSection:
    """Direct test for the full link-job write path: _write_all
    skipping modredundant for link jobs (line 112->115), the link
    section dispatch (line 124), and _write_link_section /
    _write_link_route (lines 585-600, 613-615)."""

    def test_write_all_writes_link_section_for_link_job(self):
        from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

        settings = GaussianLinkJobSettings.default()
        settings.functional = "b3lyp"
        settings.basis = "sto-3g"
        settings.charge = 0
        settings.multiplicity = 1
        settings.chk = False
        settings.jobtype = "opt"
        molecule = Molecule(
            symbols=["Ar"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        writer = _make_writer(settings, molecule=molecule, label="linkjob")
        buf = io.StringIO()
        writer._write_all(buf)
        content = buf.getvalue()
        assert "--Link1--" in content
        # link route line combines the link jobtype's route with
        # geom=check/guess=read continuation directives
        assert "geom=check" in content
        assert "guess=read" in content

    def test_write_link_section_no_op_when_link_is_false(self):
        from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

        settings = GaussianLinkJobSettings.default()
        settings.functional = "b3lyp"
        settings.basis = "sto-3g"
        settings.charge = 0
        settings.multiplicity = 1
        settings.link = False
        writer = _make_writer(settings)
        buf = io.StringIO()
        writer._write_link_section(buf)
        assert buf.getvalue() == ""


class TestGaussianInputWriterRouteSectionHeavyElementsBranches:
    """Direct tests for _write_route_section's heavy-elements basis
    replacement branches: light-basis substitution when the molecule
    has no heavy elements present, and gen/genecp keyword replacement
    when the determined basis differs from the configured one."""

    def test_no_heavy_elements_in_structure_uses_light_basis(self):
        settings = GaussianJobSettings.default()
        settings.functional = "b3lyp"
        settings.basis = "genecp"
        settings.charge = 0
        settings.multiplicity = 1
        settings.heavy_elements = ["Pd"]
        settings.heavy_elements_basis = "def2-tzvp"
        settings.light_elements_basis = "def2-svp"
        molecule = Molecule(
            symbols=["Ar"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        writer = _make_writer(settings, molecule=molecule)
        buf = io.StringIO()
        writer._write_route_section(buf)
        assert "def2svp" in buf.getvalue()
        assert "genecp" not in buf.getvalue()

    def test_heavy_elements_present_replaces_with_determined_basis(self):
        settings = GaussianJobSettings.default()
        settings.functional = "b3lyp"
        settings.basis = "genecp"
        settings.charge = 0
        settings.multiplicity = 2
        settings.heavy_elements = ["Ca"]
        settings.heavy_elements_basis = "def2-tzvp"
        settings.light_elements_basis = "def2-svp"
        molecule = Molecule(
            symbols=["Ca"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=2,
        )
        writer = _make_writer(settings, molecule=molecule)
        buf = io.StringIO()
        writer._write_route_section(buf)
        content = buf.getvalue()
        assert "gen " in content or content.strip().endswith("gen")
        assert "genecp" not in content
