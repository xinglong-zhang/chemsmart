import os
from filecmp import cmp

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca import (
    ORCAModredJob,
    ORCAOptJob,
    ORCAScanJob,
    ORCASinglePointJob,
    ORCATSJob,
)
from chemsmart.jobs.orca.writer import ORCAInputWriter
from chemsmart.settings.orca import ORCAJobSettings, ORCAProjectSettings


class TestORCAInputWriter:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_orca_project_name,
        orca_jobrunner_no_scratch,
        orca_written_opt_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_orca_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_opt",
            jobrunner=orca_jobrunner_no_scratch,
        )

        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_opt.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file, orca_written_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and the output file copied back to run folder
        # job.run()
        # assert job.is_complete()

    def test_write_opt_job_with_route(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_opt_file_with_route,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.route_to_be_written = "pbepbe/6-31g(d) opt"
        job = ORCAOptJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_opt_with_route",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_opt_with_route.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_opt_file_with_route, shallow=False)

    def test_write_modred_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_modred_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.modred_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = [[1, 2], [3, 4, 5]]
        job = ORCAModredJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_modred",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAModredJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_modred.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_modred_file, shallow=False)

    def test_write_scan_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_scan_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.scan_settings()
        settings.charge = 0
        settings.multiplicity = 1
        settings.modred = {
            "coords": [[1, 2], [3, 4]],
            "num_steps": 10,
            "dist_start": 1.5,
            "dist_end": 3.5,
        }
        job = ORCAScanJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_scan",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCAScanJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_scan.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_scan_file, shallow=False)

    def test_write_ts_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_ts_file,
    ):
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.ts_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = ORCATSJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="orca_ts",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCATSJob)
        orca_writer = ORCAInputWriter(job=job)

        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_ts.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_ts_file, shallow=False)

    def test_write_opt_input_from_logfile(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        orca_jobrunner_no_scratch,
        orca_written_ts_from_nhc_singlet_log_file,
    ):
        """Taking the Gaussian nhc_neutral_singlet.log output
        and write orca input file."""
        # get project settings
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        ts_settings = project_settings.ts_settings()
        job_settings = ORCAJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        # also merge the title keywords
        keywords = ("charge", "multiplicity", "title")
        ts_settings = ts_settings.merge(job_settings, keywords=keywords)
        job = ORCATSJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=ts_settings,
            label="orca_ts_from_log",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCATSJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_ts_from_log.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file,
            orca_written_ts_from_nhc_singlet_log_file,
            shallow=False,
        )

    def test_write_sp_input_with_solvation_from_logfile(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        gaussian_singlet_opt_outfile,
        orca_jobrunner_no_scratch,
        orca_written_sp_from_nhc_singlet_log_with_solvent_file,
    ):
        """Test writing simple .inp input file using settings from .out file,
        including solvation."""

        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        sp_settings = project_settings.sp_settings()
        job_settings = ORCAJobSettings.from_logfile(
            gaussian_singlet_opt_outfile
        )
        sp_settings = sp_settings.merge(job_settings)
        job = ORCASinglePointJob.from_filename(
            filename=gaussian_singlet_opt_outfile,
            settings=sp_settings,
            label="orca_sp_from_log_with_solvent",
            jobrunner=orca_jobrunner_no_scratch,
        )
        assert isinstance(job, ORCASinglePointJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_sp_from_log_with_solvent.inp")
        assert os.path.isfile(orca_file)
        assert cmp(
            orca_file, orca_written_sp_from_nhc_singlet_log_with_solvent_file
        )

    def test_write_opt_input_on_monoatomic_species(
        self,
        tmpdir,
        orca_yaml_settings_gas_solv_project_name,
        orca_jobrunner_no_scratch,
        orca_written_he_monoatomic_opt_file,
    ):

        helium = Molecule(
            symbols=["He"],
            positions=[(0.0, 0.0, 0.0)],
        )
        project_settings = ORCAProjectSettings.from_project(
            orca_yaml_settings_gas_solv_project_name
        )
        opt_settings = project_settings.opt_settings()
        opt_settings.charge = 0
        opt_settings.multiplicity = 1
        job = ORCAOptJob.from_molecule(
            molecule=helium,
            settings=opt_settings,
            label="orca_he_monoatomic_opt",
        )
        assert isinstance(job, ORCAOptJob)
        orca_writer = ORCAInputWriter(job=job)
        # write input file
        orca_writer.write(target_directory=tmpdir)
        orca_file = os.path.join(tmpdir, "orca_he_monoatomic_opt.inp")
        assert os.path.isfile(orca_file)
        assert cmp(orca_file, orca_written_he_monoatomic_opt_file)
