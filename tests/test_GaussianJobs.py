import os
from filecmp import cmp

import pytest

from chemsmart.jobs.gaussian import GaussianOptJob
from chemsmart.jobs.gaussian.link import GaussianLinkJob
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
            index="-1",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job,
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file = os.path.join(tmpdir, "gaussian_pbc_constraint_opt.com")
        assert os.path.isfile(g16_file)
        # check that the written file has constraints and PBC conditions
        lines = open(g16_file).readlines()
        for i in range(8):
            assert (
                " -1 " in lines[i + 8]
            )  # 8 lines before the coordinates in Gaussian input file
            # this structure has frozen atoms at these positions,
            # see test_structure.py::TestMoleculeAdvanced::test_molecule_from_db_with_pbc_and_constraints.py
            # mol object (last structure/image).

        job2 = GaussianOptJob.from_filename(
            filename=constrained_pbc_db_file,
            settings=settings,
            label="gaussian_pbc_constraint_opt2",
            index="1",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job2, GaussianOptJob)
        g16_writer = GaussianInputWriter(
            job=job2,
        )

        # write input file
        g16_writer.write(target_directory=tmpdir)
        g16_file2 = os.path.join(tmpdir, "gaussian_pbc_constraint_opt2.com")
        assert os.path.isfile(g16_file2)

        # check that the written file has constraints and PBC conditions
        lines = open(g16_file2).readlines()
        for i in range(7):
            assert (
                " -1 " in lines[i + 8]
            )  # 8 lines before the coordinates in Gaussian input file
            # this structure has frozen atoms at these positions,
            # see test_structure.py::TestMoleculeAdvanced::test_molecule_from_db_with_pbc_and_constraints.py
            # mol2 object (first structure/image).

        # job3 will fail to be created because the index is not valid (1-indexed)
        with pytest.raises(ValueError):
            GaussianOptJob.from_filename(
                filename=constrained_pbc_db_file,
                settings=settings,
                label="gaussian_pbc_constraint_opt3",
                index="0",
                jobrunner=gaussian_jobrunner_no_scratch,
            )


class TestGaussianlinkIRCJobs:
    def test_gaussian_link_irc_job_creation(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test creation of Gaussian Link IRC job."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.direction = None  # Both forward and reverse IRC

        # create link IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianLinkJob)
        assert job.settings.job_type == "irc"
        assert job._is_irc_job()

    def test_gaussian_link_irc_subjob_creation(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test creation of IRC subjobs (forward and reverse)."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"

        # create main IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc_test",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test forward IRC subjob creation
        ircf_job = job._ircf_link_job()
        assert isinstance(ircf_job, GaussianLinkJob)
        assert ircf_job.settings.job_type == "ircf"
        assert "irc_test_f" in ircf_job.label

        # test reverse IRC subjob creation
        ircr_job = job._ircr_link_job()
        assert isinstance(ircr_job, GaussianLinkJob)
        assert ircr_job.settings.job_type == "ircr"
        assert "irc_test_r" in ircr_job.label

    def test_gaussian_link_irc_job_label_naming(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test correct naming of IRC subjob labels."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.flat_irc = False

        # create main IRC job with standard link naming
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="ts1_step1_bcbond_link_irc_link",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test forward IRC subjob naming
        ircf_job = job._ircf_link_job()
        expected_ircf_label = "ts1_step1_bcbond_link_ircf_link"
        assert ircf_job.label == expected_ircf_label

        # test reverse IRC subjob naming
        ircr_job = job._ircr_link_job()
        expected_ircr_label = "ts1_step1_bcbond_link_ircr_link"
        assert ircr_job.label == expected_ircr_label

    def test_gaussian_link_irc_job_flat_naming(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test correct naming of IRC subjob labels with flat_irc option."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.flat_irc = True

        # create main IRC job with flat IRC option
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="ts1_step1_bcbond_link_irc_link",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test forward IRC subjob naming with flat
        ircf_job = job._ircf_link_job()
        expected_ircf_label = "ts1_step1_bcbond_link_ircf_flat_link"
        assert ircf_job.label == expected_ircf_label

        # test reverse IRC subjob naming with flat
        ircr_job = job._ircr_link_job()
        expected_ircr_label = "ts1_step1_bcbond_link_ircr_flat_link"
        assert ircr_job.label == expected_ircr_label

    def test_gaussian_link_irc_job_forward_only(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test IRC job with direction='forward' (forward only)."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.direction = "forward"

        # create forward-only IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc_forward_only",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test that only forward IRC job is returned
        irc_jobs = job._get_irc_jobs()
        assert len(irc_jobs) == 1
        assert irc_jobs[0].settings.job_type == "ircf"

    def test_gaussian_link_irc_job_reverse_only(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test IRC job with direction='reverse' (reverse only)."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.direction = "reverse"

        # create reverse-only IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc_reverse_only",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test that only reverse IRC job is returned
        irc_jobs = job._get_irc_jobs()
        assert len(irc_jobs) == 1
        assert irc_jobs[0].settings.job_type == "ircr"

    def test_gaussian_link_irc_job_both_directions(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test IRC job with direction=None (both directions)."""
        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.job_type = "irc"
        settings.direction = None

        # create both-directions IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc_both",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # test that both forward and reverse IRC jobs are returned
        irc_jobs = job._get_irc_jobs()
        assert len(irc_jobs) == 2
        job_types = [j.settings.job_type for j in irc_jobs]
        assert "ircf" in job_types
        assert "ircr" in job_types
