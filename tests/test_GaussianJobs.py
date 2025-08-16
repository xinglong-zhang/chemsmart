import os
from filecmp import cmp

from chemsmart.jobs.gaussian import (
    GaussianOptJob,
)
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.settings.gaussian import GaussianProjectSettings


class TestGaussianJobs:
    def test_write_opt_job(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        gaussian_written_opt_file,
    ):
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir
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
        g16_writer = GaussianInputWriter(
            job=job,
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_opt.com")
        assert os.path.isfile(g16_file)
        assert cmp(
            g16_file, gaussian_written_opt_file, shallow=False
        )  # writes input file as expected

        # job run will result in the job being run and the output file copied back to run folder
        # job.run(jobrunner=jobrunner_no_scratch)
        # assert job.is_complete()

    def test_gaussian_job_from_db_with_pbc_and_constraints(
        self,
        tmpdir,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
        constrained_pbc_db_file,
    ):
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir
        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1
        job = GaussianOptJob.from_filename(
            filename=constrained_pbc_db_file,
            settings=settings,
            label="gaussian_pbc_constraint_opt",
            index="0",  # assuming the first index is used
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job,
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_pbc_constraint_opt.com")
        print(g16_file)
        assert os.path.isfile(g16_file)
        # check that the written file has constraints and PBC conditions
        # with open(g16_file) as f:
        #     for line in f.readlines():
        #         if line.split()[0] in
        #             assert "fix" in line.lower()
        #         if "pbc" in line.lower():
        #             assert "pbc" in line.lower()
