import os
from filecmp import cmp
from unittest.mock import MagicMock, PropertyMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian import GaussianOptJob
from chemsmart.jobs.gaussian.crest import GaussianCrestJob
from chemsmart.jobs.gaussian.custom import GaussianCustomJob
from chemsmart.jobs.gaussian.dias import GaussianDIASJob
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.job import (
    GaussianComJob,
    GaussianGeneralJob,
    GaussianJob,
)
from chemsmart.jobs.gaussian.link import GaussianLinkJob
from chemsmart.jobs.gaussian.nci import GaussianNCIJob
from chemsmart.jobs.gaussian.pka import (
    GaussianpKaAnalyzeJob,
    GaussianpKaJob,
    GaussianpKaThermoJob,
)
from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
from chemsmart.jobs.gaussian.resp import GaussianRESPJob
from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianLinkJobSettings,
    GaussianpKaJobSettings,
    GaussianTDDFTJobSettings,
)
from chemsmart.jobs.gaussian.tddft import GaussianTDDFTJob
from chemsmart.jobs.gaussian.traj import GaussianTrajJob
from chemsmart.jobs.gaussian.uvvis import GaussianUVVISJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.gaussian import GaussianProjectSettings


def _make_job(settings=None):
    job = GaussianUVVISJob.__new__(GaussianUVVISJob)
    job.folder = "/tmp/uvvis_job"
    job.molecule = MagicMock(name="molecule")
    job.settings = settings or GaussianJobSettings.default()
    return job


class _FakePKaSettings(GaussianpKaJobSettings):
    """Settings double with deterministic sub-job wiring.

    Overrides the pair/settings-generation methods to return plain
    sentinel strings instead of doing real molecule-copying or
    basis/functional resolution, since ``GaussianOptJob`` and
    ``GaussianSinglePointJob`` are mocked out in these tests and never
    actually inspect these values.
    """

    def conjugate_pair_molecules(self, molecule):
        return "PROT_MOL", "CONJ_MOL"

    def _create_gas_phase_job_settings(self, molecule):
        return "PROT_OPT_SETTINGS", "CONJ_OPT_SETTINGS"

    def reference_pair_molecules(self):
        return "REF_ACID_MOL", "REF_CONJ_MOL"

    def reference_pair_job_settings(self):
        return "REF_ACID_OPT_SETTINGS", "REF_CONJ_OPT_SETTINGS"

    def _create_solution_phase_sp_settings(self, molecule):
        return "PROT_SP_SETTINGS", "CONJ_SP_SETTINGS"

    def reference_pair_sp_job_settings(self):
        return "REF_ACID_SP_SETTINGS", "REF_CONJ_SP_SETTINGS"


@pytest.fixture
def _clear_shared_reference_cache():
    # `_shared_reference_molecule_cache` is a class attribute shared by
    # every GaussianpKaJob instance; reset it so tests don't leak state.
    GaussianpKaJob._shared_reference_molecule_cache.clear()
    yield
    GaussianpKaJob._shared_reference_molecule_cache.clear()


@pytest.fixture()
def pka_settings():
    return _FakePKaSettings()


@pytest.fixture()
def pka_settings_with_reference():
    return _FakePKaSettings(
        reference_file="href.xyz",
        reference_proton_index=1,
        reference_charge=0,
        reference_multiplicity=1,
    )


def _make_job_factory():
    """Return a side_effect factory for GaussianOptJob/GaussianSinglePointJob.

    Records every created mock job in ``factory.created`` (in creation
    order) so tests can reach in and configure ``is_complete``/
    ``_output``/``outputfile`` per job.
    """

    created = []

    def factory(*, molecule, settings, label, jobrunner, skip_completed):
        job = MagicMock()
        job.molecule = molecule
        job.settings = settings
        job.label = label
        job.jobrunner = jobrunner
        job.skip_completed = skip_completed
        job.is_complete.return_value = False
        job._output.return_value = None
        job.outputfile = f"{label}.log"

        # By default, simulate a successful run: once `.run()` is
        # called the job reports itself complete. Tests that need to
        # simulate a job staying incomplete after running should clear
        # this with `job.run.side_effect = None`.
        def _mark_complete_on_run():
            job.is_complete.return_value = True

        job.run.side_effect = _mark_complete_on_run
        created.append(job)
        return job

    factory.created = created
    return factory


@pytest.fixture()
def opt_factory():
    return _make_job_factory()


@pytest.fixture()
def sp_factory():
    return _make_job_factory()


@pytest.fixture()
def patched_sub_jobs(opt_factory, sp_factory):
    with (
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianOptJob",
            side_effect=opt_factory,
        ),
        patch(
            "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob",
            side_effect=sp_factory,
        ),
    ):
        yield opt_factory, sp_factory


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

        # job run will result in the job being run and
        # the output file copied back to run folder
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
            # see test_structure.py::TestMoleculeAdvanced::
            # test_molecule_from_db_with_pbc_and_constraints.py
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
            # see test_structure.py::TestMoleculeAdvanced::
            # test_molecule_from_db_with_pbc_and_constraints.py
            # mol2 object (first structure/image).

        # job3 will fail to be created because
        # the index is not valid (1-indexed)
        with pytest.raises(ValueError):
            GaussianOptJob.from_filename(
                filename=constrained_pbc_db_file,
                settings=settings,
                label="gaussian_pbc_constraint_opt3",
                index="0",
                jobrunner=gaussian_jobrunner_no_scratch,
            )


class TestGaussianQRCJobs:

    def test_init_raises_if_no_vibrations(
        self, mock_molecule, gaussian_qrc_settings
    ):
        mock_molecule.has_vibrations = False
        with pytest.raises(ValueError, match="no vibrational modes"):
            GaussianQRCJob(
                molecule=mock_molecule, settings=gaussian_qrc_settings
            )

    def test_run_both_jobs_runs_forward_and_reverse_jobs(
        self, mock_molecule, gaussian_qrc_settings, mock_jobrunner
    ):
        job = GaussianQRCJob(
            molecule=mock_molecule,
            settings=gaussian_qrc_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )

        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianGeneralJob"
        ) as mock_general_job:
            forward_job = MagicMock()
            reverse_job = MagicMock()
            mock_general_job.side_effect = [forward_job, reverse_job]

            job._run_both_jobs()

            assert mock_general_job.call_count == 2
            forward_job.run.assert_called_once()
            reverse_job.run.assert_called_once()

    def test_run_delegates_to_run_both_jobs(
        self, mock_molecule, gaussian_qrc_settings, mock_jobrunner
    ):
        job = GaussianQRCJob(
            molecule=mock_molecule,
            settings=gaussian_qrc_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )
        with patch.object(job, "_run_both_jobs") as mock_run_both:
            job._run()
        mock_run_both.assert_called_once()

    def test_is_complete_checks_both_qrc_jobs(
        self, mock_molecule, gaussian_qrc_settings, mock_jobrunner
    ):
        job = GaussianQRCJob(
            molecule=mock_molecule,
            settings=gaussian_qrc_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )
        complete_job = MagicMock()
        complete_job.is_complete.return_value = True
        incomplete_job = MagicMock()
        incomplete_job.is_complete.return_value = False
        with patch.object(
            type(job),
            "both_qrc_jobs",
            new_callable=PropertyMock,
            return_value=[complete_job, incomplete_job],
        ):
            assert job.is_complete() is False

    def test_label_no_jobtype_suffix_when_jobtype_none(
        self, mock_molecule, gaussian_qrc_settings, mock_jobrunner
    ):
        gaussian_qrc_settings.jobtype = None
        job = GaussianQRCJob(
            molecule=mock_molecule,
            settings=gaussian_qrc_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )
        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianGeneralJob"
        ) as mock_general_job:
            job._prepare_both_qrc_jobs()

        labels = [
            call.kwargs["label"] for call in mock_general_job.call_args_list
        ]
        assert labels == ["test_qrcf", "test_qrcr"]


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
        settings.jobtype = "irc"
        settings.direction = None  # Both forward and reverse IRC

        # create link IRC job
        job = GaussianLinkJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="gaussian_link_irc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert isinstance(job, GaussianLinkJob)
        assert job.settings.jobtype == "irc"
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
        settings.jobtype = "irc"

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
        assert ircf_job.settings.jobtype == "ircf"
        assert "irc_test_f" in ircf_job.label

        # test reverse IRC subjob creation
        ircr_job = job._ircr_link_job()
        assert isinstance(ircr_job, GaussianLinkJob)
        assert ircr_job.settings.jobtype == "ircr"
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
        settings.jobtype = "irc"
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
        settings.jobtype = "irc"
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
        settings.jobtype = "irc"
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
        assert irc_jobs[0].settings.jobtype == "ircf"

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
        settings.jobtype = "irc"
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
        assert irc_jobs[0].settings.jobtype == "ircr"

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
        settings.jobtype = "irc"
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
        jobtypes = [j.settings.jobtype for j in irc_jobs]
        assert "ircf" in jobtypes
        assert "ircr" in jobtypes


class TestGaussianIRCJobs:
    def test_gaussian_irc_job_flat_naming(
        self,
        tmpdir,
        single_molecule_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test correct naming of IRC sub-job labels with flat_irc option.

        For a flat IRC job with no direction specified (both forward and
        reverse), the sub-job labels should be ``{label}_ircf_flat`` and
        ``{label}_ircr_flat``.

        - Auto-generated label (``file_irc``, already ends with ``_irc``):
          sub-jobs become ``file_ircf_flat`` / ``file_ircr_flat``.
        - Custom label via ``-l label`` (no ``_irc`` suffix):
          sub-jobs become ``label_ircf_flat`` / ``label_ircr_flat``.
        """
        from chemsmart.jobs.gaussian.irc import GaussianIRCJob

        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.irc_settings()
        settings.charge = -2
        settings.multiplicity = 1
        settings.jobtype = "irc"
        settings.flat_irc = True

        # --- auto-generated label case (label ends with _irc) ---
        # simulates `gaussian -f file.xyz irc --flat`
        # where the gaussian group sets label="file_irc"
        job_auto = GaussianIRCJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="job_label_irc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert job_auto._ircf_job().label == "job_label_ircf_flat"
        assert job_auto._ircr_job().label == "job_label_ircr_flat"

        # --- custom label case (label does NOT end with _irc) ---
        # simulates `gaussian -l label irc --flat`
        job_custom = GaussianIRCJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=settings,
            label="label",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert job_custom._ircf_job().label == "label_ircf_flat"
        assert job_custom._ircr_job().label == "label_ircr_flat"


class TestGaussianCrestGroupingForwarding:
    def test_crest_grouping_forwards_num_procs_and_representative(
        self,
        tmpdir,
        multiple_molecules_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.settings.gaussian import GaussianProjectSettings

        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1

        molecules = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index=":", return_list=True
        )

        mock_grouper = MagicMock()
        mock_grouper.group.return_value = None
        mock_grouper.unique.return_value = molecules[:1]

        with patch(
            "chemsmart.jobs.gaussian.crest.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            GaussianCrestJob(
                molecules=molecules,
                settings=settings,
                label="crest_forwarding_test",
                jobrunner=gaussian_jobrunner_no_scratch,
                grouping_strategy="rmsd",
                num_procs=4,
                representative_strategy="center",
            )
        mock_create.assert_called_once()
        create_kwargs = mock_create.call_args.kwargs
        assert create_kwargs["num_procs"] == 4
        assert create_kwargs["representative_strategy"] == "center"


class TestGaussianTrajJobs:
    def test_traj_job_uses_all_molecules(
        self,
        tmpdir,
        multiple_molecules_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        """Test that GaussianTrajJob uses all molecules from the file."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob
        from chemsmart.settings.gaussian import GaussianProjectSettings

        # set scratch directory for jobrunner
        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir

        # get project settings
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1

        # Read all molecules from the file
        molecules = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index=":", return_list=True
        )
        num_molecules_in_file = len(molecules)

        # Create Traj job with all molecules
        # Using proportion_structures_to_use=1.0 to use all structures
        job = GaussianTrajJob(
            molecules=molecules,
            settings=settings,
            label="traj_test_opt",
            jobrunner=gaussian_jobrunner_no_scratch,
            proportion_structures_to_use=1.0,
        )

        # Check that the job uses all molecules (no grouping by default)
        assert job.num_structures == num_molecules_in_file

        # When no grouping strategy is provided, all structures are unique
        assert job.num_unique_structures == num_molecules_in_file

    def test_traj_grouping_forwards_num_procs_and_representative(
        self,
        tmpdir,
        multiple_molecules_xyz_file,
        gaussian_yaml_settings_gas_solv_project_name,
        gaussian_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob
        from chemsmart.settings.gaussian import GaussianProjectSettings

        gaussian_jobrunner_no_scratch.scratch_dir = tmpdir
        project_settings = GaussianProjectSettings.from_project(
            gaussian_yaml_settings_gas_solv_project_name
        )
        settings = project_settings.opt_settings()
        settings.charge = 0
        settings.multiplicity = 1

        molecules = Molecule.from_filepath(
            filepath=multiple_molecules_xyz_file, index=":", return_list=True
        )

        mock_grouper = MagicMock()
        mock_grouper.group.return_value = None

        with patch(
            "chemsmart.jobs.gaussian.traj.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            GaussianTrajJob(
                molecules=molecules,
                settings=settings,
                label="traj_forwarding_test",
                jobrunner=gaussian_jobrunner_no_scratch,
                grouping_strategy="rmsd",
                num_procs=4,
                representative_strategy="top3",
                proportion_structures_to_use=1.0,
            )
        mock_create.assert_called_once()
        create_kwargs = mock_create.call_args.kwargs
        assert create_kwargs["num_procs"] == 4
        assert create_kwargs["representative_strategy"] == "top3"


class TestGaussianJobConstruction:
    def test_rejects_invalid_settings_type(self, a_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            GaussianJob(molecule=a_molecule, settings="not-settings")

    def test_rejects_invalid_molecule_type(self, gaussian_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianJob(molecule="not-a-molecule", settings=gaussian_settings)

    def test_label_defaults_to_chemical_formula(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob(molecule=a_molecule, settings=gaussian_settings)
        assert job.label == a_molecule.get_chemical_formula(empirical=True)

    def test_custom_label_used(self, a_molecule, gaussian_settings):
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job.label == "mylabel"

    def test_settings_class_returns_expected_type(self):
        assert GaussianJob.settings_class() is GaussianJobSettings


class TestGaussianJobPaths:
    def test_input_output_chk_err_file_paths(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job.inputfile == os.path.join(job.folder, "mylabel.com")
        assert job.outputfile == os.path.join(job.folder, "mylabel.log")
        assert job.chkfile == os.path.join(job.folder, "mylabel.chk")
        assert job.errfile == os.path.join(job.folder, "mylabel.err")


class TestGaussianJobBackupAndOutput:
    def test_backup_files_copies_input_and_output(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.inputfile, "w") as f:
            f.write("input")
        with open(job.outputfile, "w") as f:
            f.write("output")

        job._backup_files()

        backups = job._previous_backup_folders()
        assert len(backups) == 1
        assert os.path.exists(os.path.join(backups[0], "mylabel.com"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.log"))

    def test_backup_files_includes_chk_when_requested(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        for f in [job.inputfile, job.outputfile, job.chkfile]:
            with open(f, "w") as fh:
                fh.write("data")

        job._backup_files(backup_chk=True)

        backups = job._previous_backup_folders()
        assert os.path.exists(os.path.join(backups[0], "mylabel.chk"))

    def test_output_none_when_outputfile_missing(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job._output() is None

    def test_output_returns_gaussian_output_when_present(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_output = MagicMock()
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        ) as mock_cls:
            result = job._output()

        mock_cls.assert_called_once_with(filename=job.outputfile)
        assert result is mock_output

    def test_output_none_on_attribute_error(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            side_effect=AttributeError("boom"),
        ):
            assert job._output() is None

    def test_output_falls_back_to_pbc_output_on_value_error(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_pbc_output = MagicMock()
        with (
            patch(
                "chemsmart.io.gaussian.output.Gaussian16Output",
                side_effect=ValueError("not a standard output"),
            ),
            patch(
                "chemsmart.io.gaussian.output.Gaussian16OutputWithPBC",
                return_value=mock_pbc_output,
            ) as mock_pbc_cls,
        ):
            result = job._output()

        mock_pbc_cls.assert_called_once_with(filename=job.outputfile)
        assert result is mock_pbc_output


class TestGaussianJobRun:
    def test_run_delegates_to_jobrunner(self, a_molecule, gaussian_settings):
        mock_runner = MagicMock()
        job = GaussianJob(
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestGaussianJobFactories:
    def test_from_filename_builds_job(
        self, single_molecule_xyz_file, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GaussianJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_filename_without_explicit_jobrunner_creation_attempted(
        self, single_molecule_xyz_file, gaussian_settings
    ):
        """The base GaussianJob has no TYPE registered with any runner,
        so JobRunner.from_job() genuinely can't find one -- but this
        still exercises the "if jobrunner is None:" branch itself
        before failing."""
        with pytest.raises(ValueError, match="Could not find any runners"):
            GaussianJob.from_filename(
                filename=single_molecule_xyz_file,
                settings=gaussian_settings,
            )

    def test_from_pubchem_builds_job(self, gaussian_settings):
        pubchem_molecule = MagicMock(spec=Molecule)
        pubchem_molecule.get_chemical_formula.return_value = "H2O"
        pubchem_molecule.copy.return_value = pubchem_molecule
        mock_jobrunner = MagicMock()
        with patch(
            "chemsmart.jobs.gaussian.job.Molecule.from_pubchem",
            return_value=pubchem_molecule,
        ) as mock_from_pubchem:
            job = GaussianJob.from_pubchem(
                identifier="water",
                settings=gaussian_settings,
                jobrunner=mock_jobrunner,
            )

        mock_from_pubchem.assert_called_once_with(identifier="water")
        assert isinstance(job, GaussianJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_pubchem_without_explicit_jobrunner_creation_attempted(
        self, gaussian_settings
    ):
        pubchem_molecule = MagicMock(spec=Molecule)
        pubchem_molecule.get_chemical_formula.return_value = "H2O"
        pubchem_molecule.copy.return_value = pubchem_molecule
        with patch(
            "chemsmart.jobs.gaussian.job.Molecule.from_pubchem",
            return_value=pubchem_molecule,
        ):
            with pytest.raises(ValueError, match="Could not find any runners"):
                GaussianJob.from_pubchem(
                    identifier="water", settings=gaussian_settings
                )

    def test_from_jobtype_opt(self, a_molecule, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_jobtype(
            jobtype="opt",
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        from chemsmart.jobs.gaussian.opt import GaussianOptJob

        assert isinstance(job, GaussianOptJob)

    def test_from_jobtype_opt_without_explicit_jobrunner_creates_one(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob.from_jobtype(
            jobtype="opt", molecule=a_molecule, settings=gaussian_settings
        )
        from chemsmart.jobs.gaussian.opt import GaussianOptJob

        assert isinstance(job, GaussianOptJob)
        assert job.jobrunner is not None

    def test_from_jobtype_com(self, a_molecule, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_jobtype(
            jobtype="com",
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GaussianComJob)

    def test_from_jobtype_com_without_explicit_jobrunner_creates_one(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob.from_jobtype(
            jobtype="com", molecule=a_molecule, settings=gaussian_settings
        )
        assert isinstance(job, GaussianComJob)
        assert job.jobrunner is not None

    def test_from_jobtype_g16_without_explicit_jobrunner_works(
        self, a_molecule, gaussian_settings
    ):
        """When ``jobrunner`` is *not* supplied, the "g16" branch creates
        one and correctly returns a ``GaussianGeneralJob``."""
        job = GaussianJob.from_jobtype(
            jobtype="g16",
            molecule=a_molecule,
            settings=gaussian_settings,
        )
        assert isinstance(job, GaussianGeneralJob)

    def test_from_jobtype_g16_with_explicit_jobrunner_is_buggy(
        self, a_molecule, gaussian_settings
    ):
        """Documents a real bug in ``GaussianJob.from_jobtype``: the
        ``return GaussianGeneralJob(...)`` for the "g16" branch is nested
        inside ``if jobrunner is None:`` instead of the ``elif
        jobtype.lower() == "g16":`` block. When a caller *already* has a
        jobrunner (the common case — see e.g. the ``qmmm`` subcommands
        that build parent jobs with a jobrunner in hand), this falls
        through to the ``else: raise ValueError(f"Invalid job type:
        {jobtype}")`` clause even though "g16" is a valid job type."""
        mock_jobrunner = MagicMock()
        with pytest.raises(ValueError, match="Invalid job type: g16"):
            GaussianJob.from_jobtype(
                jobtype="g16",
                molecule=a_molecule,
                settings=gaussian_settings,
                jobrunner=mock_jobrunner,
            )

    def test_from_jobtype_invalid_silently_returns_none(
        self, a_molecule, gaussian_settings
    ):
        """Documents a second bug in ``from_jobtype``: the "invalid job
        type" ``raise ValueError`` is only reachable from inside the
        "g16" branch (see the "buggy" test above). There is no top-level
        ``else`` for jobtypes that are not "opt"/"com"/"g16" at all, so a
        genuinely unrecognized jobtype silently falls through and returns
        ``None`` instead of raising."""
        result = GaussianJob.from_jobtype(
            jobtype="bogus",
            molecule=a_molecule,
            settings=gaussian_settings,
        )
        assert result is None


class TestGaussianComJob:
    def test_from_filename_crashes_on_none_molecule(
        self, gaussian_opt_inputfile
    ):
        """Documents a third bug: ``GaussianComJob.from_filename`` always
        constructs the job with ``molecule=None`` (it only reads route
        info, not coordinates, from the ``.com`` file), but the parent
        ``GaussianJob.__init__`` unconditionally requires
        ``isinstance(molecule, Molecule)``. So this factory method
        currently cannot succeed at all."""
        mock_jobrunner = MagicMock()
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianComJob.from_filename(
                filename=gaussian_opt_inputfile, jobrunner=mock_jobrunner
            )

    def test_from_filename_with_explicit_label_still_crashes_on_molecule(
        self, gaussian_opt_inputfile
    ):
        """Same underlying bug as above, but confirms the label-is-not-None
        branch (an explicit label skips deriving one from the filename)
        is reached before the same molecule=None crash."""
        mock_jobrunner = MagicMock()
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianComJob.from_filename(
                filename=gaussian_opt_inputfile,
                label="custom_label",
                jobrunner=mock_jobrunner,
            )

    def test_from_filename_without_jobrunner_still_crashes_on_molecule(
        self, gaussian_opt_inputfile
    ):
        """Without an explicit jobrunner, the factory tries to build one
        via JobRunner.from_job(cls(molecule=None, ...)) first -- which
        hits the same molecule=None crash even earlier."""
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianComJob.from_filename(filename=gaussian_opt_inputfile)


class TestGaussianGeneralJob:
    def test_construction(self, a_molecule, gaussian_settings):
        job = GaussianGeneralJob(
            molecule=a_molecule, settings=gaussian_settings
        )
        assert job.TYPE == "g16"
        assert isinstance(job, GaussianJob)


class TestGaussianIRCJobConstruction:
    def test_freq_disabled(self, a_molecule, gaussian_irc_settings):
        gaussian_irc_settings.freq = True
        job = GaussianIRCJob(
            molecule=a_molecule, settings=gaussian_irc_settings
        )
        assert job.settings.freq is False

    def test_settings_class_returns_expected_type(self):
        assert GaussianIRCJob.settings_class() is GaussianIRCJobSettings


class TestSubJobLabelDerivation:
    def test_ircf_label_auto_generated_irc_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf"

    def test_ircr_label_auto_generated_irc_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircr_job().label == "mymol_ircr"

    def test_ircf_label_custom_label_gets_underscored_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="custom_label",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "custom_label_ircf"

    def test_ircr_label_custom_label_gets_underscored_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="custom_label",
            jobrunner=MagicMock(),
        )
        assert job._ircr_job().label == "custom_label_ircr"

    def test_ircf_label_flat_irc_appends_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.flat_irc = True
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf_flat"

    def test_label_unchanged_when_direction_already_set(
        self, a_molecule, gaussian_irc_settings
    ):
        # When direction is already embedded (e.g. by update_irc_label
        # upstream), the sub-job label is reused verbatim.
        gaussian_irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_ircf",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf"

    def test_label_none_defaults_to_chemical_formula_based_label(
        self, a_molecule, gaussian_irc_settings
    ):
        # GaussianJob.__init__ defaults a None label to the molecule's
        # chemical formula, so job.label is never actually None by the
        # time _ircf_job() runs; it just goes through the "custom
        # label" branch (no "_irc" suffix to detect) like any other
        # non-None label.
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label=None,
            jobrunner=MagicMock(),
        )
        expected_label = a_molecule.get_chemical_formula(empirical=True)
        assert job.label == expected_label
        assert job._ircf_job().label == f"{expected_label}_ircf"

    def test_ircf_label_none_skips_label_derivation(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.label = None
        # Sub-job label derivation is skipped entirely; the sub-job's own
        # GaussianJob.__init__ then defaults the None label to the
        # molecule's chemical formula (no "_ircf" suffix applied).
        expected_label = a_molecule.get_chemical_formula(empirical=True)
        assert job._ircf_job().label == expected_label

    def test_ircr_label_none_skips_label_derivation(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.label = None
        expected_label = a_molecule.get_chemical_formula(empirical=True)
        assert job._ircr_job().label == expected_label

    def test_ircr_label_flat_irc_appends_suffix(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.flat_irc = True
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircr_job().label == "mymol_ircr_flat"

    def test_ircf_jobtype_set_on_sub_settings(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().settings.jobtype == "ircf"
        assert job._ircr_job().settings.jobtype == "ircr"


class TestRunAndCompleteHelpers:
    def test_run_forward_runs_the_ircf_job(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        mock_ircf_job = MagicMock()
        job._ircf_job = MagicMock(return_value=mock_ircf_job)

        job._run_forward()

        mock_ircf_job.run.assert_called_once()

    def test_run_reverse_runs_the_ircr_job(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        mock_ircr_job = MagicMock()
        job._ircr_job = MagicMock(return_value=mock_ircr_job)

        job._run_reverse()

        mock_ircr_job.run.assert_called_once()

    def test_run_forward_is_complete_delegates_to_ircf_job(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        mock_ircf_job = MagicMock()
        mock_ircf_job.is_complete.return_value = True
        job._ircf_job = MagicMock(return_value=mock_ircf_job)

        assert job._run_forward_is_complete() is True

    def test_run_reverse_is_complete_delegates_to_ircr_job(
        self, a_molecule, gaussian_irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        mock_ircr_job = MagicMock()
        mock_ircr_job.is_complete.return_value = False
        job._ircr_job = MagicMock(return_value=mock_ircr_job)

        assert job._run_reverse_is_complete() is False


class TestRunDirectionDispatch:
    def test_forward_direction_runs_only_forward(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward = MagicMock()
        job._run_reverse = MagicMock()

        job._run()

        job._run_forward.assert_called_once()
        job._run_reverse.assert_not_called()

    def test_reverse_direction_runs_only_reverse(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward = MagicMock()
        job._run_reverse = MagicMock()

        job._run()

        job._run_forward.assert_not_called()
        job._run_reverse.assert_called_once()

    def test_no_direction_runs_both(self, a_molecule, gaussian_irc_settings):
        gaussian_irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward = MagicMock()
        job._run_reverse = MagicMock()

        job._run()

        job._run_forward.assert_called_once()
        job._run_reverse.assert_called_once()


class TestJobIsCompleteDirectionDispatch:
    def test_forward_direction_checks_only_forward(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=True)
        job._run_reverse_is_complete = MagicMock(return_value=False)

        assert job._job_is_complete() is True
        job._run_reverse_is_complete.assert_not_called()

    def test_reverse_direction_checks_only_reverse(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=False)
        job._run_reverse_is_complete = MagicMock(return_value=True)

        assert job._job_is_complete() is True
        job._run_forward_is_complete.assert_not_called()

    def test_no_direction_requires_both_complete(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=True)
        job._run_reverse_is_complete = MagicMock(return_value=False)

        assert job._job_is_complete() is False


class TestBackupFiles:
    def test_forward_direction_backs_up_only_forward_files(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=True)
        # inputfile, outputfile, chkfile -> 3 calls for forward only
        assert job.backup_file.call_count == 3

    def test_forward_direction_without_chk_skips_chkfile(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=False)
        # inputfile, outputfile only (no chk) -> 2 calls
        assert job.backup_file.call_count == 2

    def test_reverse_direction_backs_up_only_reverse_files(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=False)
        # inputfile, outputfile only (no chk) -> 2 calls
        assert job.backup_file.call_count == 2

    def test_reverse_direction_with_chk_backs_up_chkfile(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=True)
        # inputfile, outputfile, chkfile -> 3 calls for reverse only
        assert job.backup_file.call_count == 3

    def test_no_direction_backs_up_both_directions(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=True)
        # input+output for both directions (4) + chk for both (2) = 6
        assert job.backup_file.call_count == 6

    def test_no_direction_without_chk_skips_chkfiles(
        self, a_molecule, gaussian_irc_settings
    ):
        gaussian_irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=gaussian_irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=False)
        # input+output for both directions only -> 4 calls
        assert job.backup_file.call_count == 4


class TestSettingsClass:
    def test_returns_gaussian_link_job_settings(self):
        assert GaussianLinkJob.settings_class() is GaussianLinkJobSettings


class TestIsIrcJob:
    def test_true_for_irc_jobtype(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        assert job._is_irc_job() is True

    def test_false_for_non_irc_jobtype(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        assert job._is_irc_job() is False


class TestCreateIrcSubjob:
    def test_raises_for_non_irc_job(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        with pytest.raises(ValueError, match="only for IRC link jobs"):
            job._create_irc_subjob("f")

    def test_forward_subjob_has_ircf_jobtype(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        subjob = job._create_irc_subjob("f")
        assert subjob.settings.jobtype == "ircf"
        assert subjob.label == "link_test_f"

    def test_reverse_subjob_has_ircr_jobtype(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        subjob = job._create_irc_subjob("r")
        assert subjob.settings.jobtype == "ircr"
        assert subjob.label == "link_test_r"

    def test_irc_link_label_pattern_is_rewritten_in_place(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="mymol_irc_link",
        )
        subjob = job._create_irc_subjob("f")
        assert subjob.label == "mymol_ircf_link"

    def test_flat_irc_appends_flat_suffix(self, ethanol_molecule):
        settings = GaussianLinkJobSettings(
            jobtype="irc", direction=None, flat_irc=True
        )
        job = GaussianLinkJob(
            molecule=ethanol_molecule, settings=settings, label="link_test"
        )
        subjob = job._create_irc_subjob("f")
        assert subjob.label == "link_test_f_flat"

    def test_flat_irc_with_irc_link_pattern(self, ethanol_molecule):
        settings = GaussianLinkJobSettings(
            jobtype="irc", direction=None, flat_irc=True
        )
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=settings,
            label="mymol_irc_link",
        )
        subjob = job._create_irc_subjob("r")
        assert subjob.label == "mymol_ircr_flat_link"

    def test_none_label_skips_label_derivation(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        job.label = None
        subjob = job._create_irc_subjob("f")
        assert subjob.label == ethanol_molecule.get_chemical_formula(
            empirical=True
        )


class TestGetIrcJobs:
    def test_non_irc_job_returns_empty_list(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        assert job._get_irc_jobs() == []

    def test_forward_only_direction(self, ethanol_molecule):
        settings = GaussianLinkJobSettings(jobtype="irc", direction="forward")
        job = GaussianLinkJob(
            molecule=ethanol_molecule, settings=settings, label="link_test"
        )
        jobs = job._get_irc_jobs()
        assert len(jobs) == 1
        assert jobs[0].settings.jobtype == "ircf"

    def test_reverse_only_direction(self, ethanol_molecule):
        settings = GaussianLinkJobSettings(jobtype="irc", direction="reverse")
        job = GaussianLinkJob(
            molecule=ethanol_molecule, settings=settings, label="link_test"
        )
        jobs = job._get_irc_jobs()
        assert len(jobs) == 1
        assert jobs[0].settings.jobtype == "ircr"

    def test_no_direction_returns_both(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        jobs = job._get_irc_jobs()
        assert [j.settings.jobtype for j in jobs] == ["ircf", "ircr"]


class TestRunAndCompletion:
    def test_run_delegates_to_super_for_non_irc(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        with patch(
            "chemsmart.jobs.gaussian.job.GaussianJob._run"
        ) as mock_super_run:
            job._run()
        mock_super_run.assert_called_once()

    def test_run_runs_each_irc_subjob(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        fwd, rev = MagicMock(), MagicMock()
        fwd.settings.jobtype = "ircf"
        rev.settings.jobtype = "ircr"
        with patch.object(job, "_get_irc_jobs", return_value=[fwd, rev]):
            job._run()
        fwd.run.assert_called_once()
        rev.run.assert_called_once()

    def test_job_is_complete_non_irc_delegates_to_super(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        with patch(
            "chemsmart.jobs.gaussian.job.GaussianJob._job_is_complete",
            return_value=True,
        ) as mock_super_complete:
            assert job._job_is_complete() is True
        mock_super_complete.assert_called_once()

    def test_job_is_complete_irc_requires_all_subjobs_complete(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        fwd, rev = MagicMock(), MagicMock()
        fwd._job_is_complete.return_value = True
        rev._job_is_complete.return_value = False
        with patch.object(job, "_get_irc_jobs", return_value=[fwd, rev]):
            assert job._job_is_complete() is False

        fwd._job_is_complete.return_value = True
        rev._job_is_complete.return_value = True
        with patch.object(job, "_get_irc_jobs", return_value=[fwd, rev]):
            assert job._job_is_complete() is True


class TestBackupFilesIsDeadCode:
    def test_backup_files_is_not_the_method_invoked_by_backup(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        """Documents a bug: ``GaussianLinkJob`` defines a public
        ``backup_files(self, backup_chk=False)`` method with IRC-aware
        subjob backup logic, but ``Job.backup()`` (the only caller in
        the codebase) invokes ``self._backup_files(**kwargs)`` — the
        *private*, underscore-prefixed name. ``GaussianLinkJob`` never
        overrides ``_backup_files``, so it silently inherits
        ``GaussianJob._backup_files`` unchanged, and the IRC-subjob-aware
        ``backup_files`` override is unreachable dead code."""
        assert GaussianLinkJob._backup_files is GaussianJob._backup_files
        assert "backup_files" in vars(GaussianLinkJob)

    def test_backup_files_method_itself_works_when_called_directly(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        fwd, rev = MagicMock(), MagicMock()
        fwd.inputfile = "fwd.com"
        fwd.outputfile = "fwd.log"
        fwd.chkfile = "fwd.chk"
        rev.inputfile = "rev.com"
        rev.outputfile = "rev.log"
        rev.chkfile = "rev.chk"

        with (
            patch.object(job, "_get_irc_jobs", return_value=[fwd, rev]),
            patch.object(job, "backup_file") as mock_backup_file,
        ):
            job.backup_files(backup_chk=True)

        mock_backup_file.assert_any_call("fwd.com")
        mock_backup_file.assert_any_call("fwd.log")
        mock_backup_file.assert_any_call("fwd.chk")
        mock_backup_file.assert_any_call("rev.com")
        mock_backup_file.assert_any_call("rev.log")
        mock_backup_file.assert_any_call("rev.chk")

    def test_backup_files_with_backup_chk_false_skips_chk(
        self, ethanol_molecule, gaussian_link_irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_irc_settings,
            label="link_test",
        )
        fwd = MagicMock()
        fwd.inputfile = "fwd.com"
        fwd.outputfile = "fwd.log"
        fwd.chkfile = "fwd.chk"

        with (
            patch.object(job, "_get_irc_jobs", return_value=[fwd]),
            patch.object(job, "backup_file") as mock_backup_file,
        ):
            job.backup_files(backup_chk=False)

        called_paths = [
            call.args[0] for call in mock_backup_file.call_args_list
        ]
        assert called_paths == ["fwd.com", "fwd.log"]

    def test_backup_files_non_irc_delegates_to_super_backup_files(
        self, ethanol_molecule, gaussian_link_opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=gaussian_link_opt_settings,
            label="link_test",
        )
        with patch.object(
            GaussianJob, "_backup_files"
        ) as mock_super_backup_files:
            job.backup_files(backup_chk=True)
        mock_super_backup_files.assert_called_once_with(True)


class TestGaussianCrestJobConstruction:
    def test_defaults_num_confs_to_opt_to_all_conformers(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        assert job.num_conformers == 4
        assert job.num_confs_to_opt == 4
        assert job.all_conformers == ethanol_conformers
        assert job.grouper is None

    def test_num_confs_to_run_is_respected(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            num_confs_to_run=2,
        )
        assert job.num_confs_to_opt == 2
        # all_conformers still holds every molecule; only the "opt" count
        # is restricted
        assert job.num_conformers == 4

    def test_uses_first_molecule_for_base_job_construction(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        assert (
            job.molecule.get_chemical_formula()
            == ethanol_conformers[0].get_chemical_formula()
        )

    def test_empty_molecules_list_bypasses_validation_and_crashes(self):
        """Documents a bug (see BUGS_FOUND.md): the guard
        ``if not isinstance(molecules, list) and len(molecules) == 0``
        uses ``and`` where the intent (per the raised message) was
        clearly ``or``. For an empty list, ``not isinstance([], list)``
        is ``False``, so the whole condition short-circuits to ``False``
        and the intended ``ValueError`` is never raised. Execution falls
        through to ``molecules[0]`` a few lines later and crashes with an
        unrelated ``IndexError`` instead of the documented, actionable
        error."""
        with pytest.raises(IndexError):
            GaussianCrestJob(molecules=[])

    def test_non_list_non_empty_molecules_also_bypasses_validation(
        self, gaussian_settings, ethanol_molecule
    ):
        """Same root cause as above: a tuple is ``not isinstance(x,
        list)`` (True) but has nonzero length, so ``True and False`` is
        ``False`` and the "must be a list" check never fires even though
        a tuple (not a list) was passed. Construction proceeds normally
        instead of raising the documented ``ValueError``."""
        molecules = (ethanol_molecule.copy(), ethanol_molecule.copy())
        job = GaussianCrestJob(
            molecules=molecules,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        assert job.num_conformers == 2

    def test_empty_tuple_molecules_raises_the_intended_valueerror(self):
        """The only input that actually satisfies the buggy ``and``
        guard: an empty, non-list iterable (e.g. an empty tuple) is
        both ``not isinstance(x, list)`` (True) and ``len(x) == 0``
        (True), so this is the one case where the documented
        ``ValueError`` fires as intended."""
        with pytest.raises(ValueError, match="must be a list"):
            GaussianCrestJob(molecules=())


class TestGaussianCrestJobProperties:
    def _job(self, ethanol_conformers, gaussian_settings, **kwargs):
        return GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            **kwargs,
        )

    def test_all_conformers_jobs_builds_one_job_per_conformer_with_indexed_labels(
        self, ethanol_conformers, gaussian_settings
    ):
        job = self._job(ethanol_conformers, gaussian_settings)
        sub_jobs = job.all_conformers_jobs
        assert len(sub_jobs) == 4
        assert [j.label for j in sub_jobs] == [
            "crest_test_c1",
            "crest_test_c2",
            "crest_test_c3",
            "crest_test_c4",
        ]

    def test_incomplete_conformers_jobs_filters_by_is_complete(
        self, ethanol_conformers, gaussian_settings
    ):
        job = self._job(ethanol_conformers, gaussian_settings)
        statuses = [True, False, True, False]
        fake_jobs = []
        for status in statuses:
            fake_job = MagicMock()
            fake_job.is_complete.return_value = status
            fake_jobs.append(fake_job)

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            incomplete = job.incomplete_conformers_jobs
        assert incomplete == [fake_jobs[1], fake_jobs[3]]

    def test_last_run_job_index_returns_first_incomplete(
        self, ethanol_conformers, gaussian_settings
    ):
        job = self._job(ethanol_conformers, gaussian_settings)
        statuses = [True, True, False, False]
        fake_jobs = []
        for status in statuses:
            fake_job = MagicMock()
            fake_job.is_complete.return_value = status
            fake_jobs.append(fake_job)

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.last_run_job_index == 2

    def test_last_run_job_index_is_total_count_when_all_complete(
        self, ethanol_conformers, gaussian_settings
    ):
        job = self._job(ethanol_conformers, gaussian_settings)
        fake_jobs = [MagicMock(is_complete=lambda: True) for _ in range(4)]

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.last_run_job_index == 4


class TestGaussianCrestJobRunAndCompletion:
    def test_run_all_jobs_only_runs_up_to_num_confs_to_opt(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            num_confs_to_run=2,
        )
        fake_jobs = [MagicMock() for _ in range(4)]

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            job._run_all_jobs()

        fake_jobs[0].run.assert_called_once()
        fake_jobs[1].run.assert_called_once()
        fake_jobs[2].run.assert_not_called()
        fake_jobs[3].run.assert_not_called()

    def test_run_delegates_to_run_all_jobs(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        with patch.object(job, "_run_all_jobs") as mock_run_all:
            job._run()
        mock_run_all.assert_called_once()

    def test_is_complete_true_only_when_all_within_limit_are_complete(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            num_confs_to_run=2,
        )
        # jobs 3 and 4 (index 2, 3) are incomplete but outside the limit
        fake_jobs = [
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: False),
            MagicMock(is_complete=lambda: False),
        ]

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.is_complete() is True

    def test_is_complete_false_when_a_required_job_is_incomplete(
        self, ethanol_conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=ethanol_conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            num_confs_to_run=3,
        )
        fake_jobs = [
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: False),
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: True),
        ]

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.is_complete() is False


class TestGaussianCrestJobGrouping:
    def test_grouping_strategy_replaces_all_conformers_with_uniques(
        self, ethanol_conformers, gaussian_settings
    ):
        unique_subset = ethanol_conformers[:2]
        mock_grouper = MagicMock()
        mock_grouper.unique.return_value = unique_subset

        with patch(
            "chemsmart.utils.grouper.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            job = GaussianCrestJob(
                molecules=ethanol_conformers,
                settings=gaussian_settings,
                label="crest_test",
                jobrunner=None,
                grouping_strategy="rmsd",
                num_groups=2,
            )

        mock_create.assert_called_once()
        _, kwargs = mock_create.call_args
        assert kwargs["strategy"] == "rmsd"
        assert kwargs["num_groups"] == 2
        mock_grouper.group.assert_called_once()
        assert job.all_conformers == unique_subset
        assert job.grouper is mock_grouper

    def test_grouping_strategy_without_num_groups_uses_threshold_mode(
        self, ethanol_conformers, gaussian_settings
    ):
        unique_subset = ethanol_conformers[:3]
        mock_grouper = MagicMock()
        mock_grouper.unique.return_value = unique_subset
        with patch(
            "chemsmart.utils.grouper.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            job = GaussianCrestJob(
                molecules=ethanol_conformers,
                settings=gaussian_settings,
                label="crest_test",
                jobrunner=None,
                grouping_strategy="rmsd",
            )
        _, kwargs = mock_create.call_args
        assert kwargs["num_groups"] is None
        assert job.all_conformers == unique_subset


class TestGaussianTrajJobValidation:
    def test_non_list_molecules_raises(self, gaussian_settings):
        with pytest.raises(ValueError, match="must be a list"):
            GaussianTrajJob(
                molecules="not-a-list",
                settings=gaussian_settings,
                label="traj",
                jobrunner=None,
            )

    def test_empty_list_raises(self, gaussian_settings):
        with pytest.raises(ValueError, match="cannot be empty"):
            GaussianTrajJob(
                molecules=[],
                settings=gaussian_settings,
                label="traj",
                jobrunner=None,
            )

    def test_non_molecule_elements_raise(self, gaussian_settings):
        with pytest.raises(ValueError, match="instances of Molecule"):
            GaussianTrajJob(
                molecules=[1, 2, 3],
                settings=gaussian_settings,
                label="traj",
                jobrunner=None,
            )


class TestGaussianTrajJobTailExtraction:
    def test_default_proportion_extracts_tail_fraction(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        # 20 * 0.1 = 2.0 -> last 2 structures
        assert job.num_structures == 2
        assert job.molecules == ethanol_traj_molecules[-2:]

    def test_small_trajectory_clamps_to_at_least_one(
        self, ethanol_molecule, gaussian_settings
    ):
        mols = [ethanol_molecule.copy() for _ in range(3)]
        for i, mol in enumerate(mols):
            mol.energy = float(i)
        job = GaussianTrajJob(
            molecules=mols,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=0.1,
        )
        # 3 * 0.1 = 0.3 -> would round/truncate to 0, clamped to 1
        assert job.num_structures == 1
        assert job.molecules == mols[-1:]

    def test_proportion_one_uses_entire_trajectory(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=1.0,
        )
        assert job.num_structures == 20
        assert job.molecules == ethanol_traj_molecules

    def test_original_conformer_ids_track_1_based_source_indices(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        # 20 structures, last 2 kept -> original indices 19, 20 (1-based)
        assert job._original_conformer_ids == ["19", "20"]


class TestGaussianTrajJobEnergySorting:
    def _job(self, ethanol_traj_molecules, gaussian_settings):
        return GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=1.0,
        )

    def test_all_energies_matches_molecule_order(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job(ethanol_traj_molecules, gaussian_settings)
        assert job.all_energies == [float(i) for i in range(20)]

    def test_sorted_energies_ascending(
        self, ethanol_molecule, gaussian_settings
    ):
        mols = [ethanol_molecule.copy() for _ in range(4)]
        energies = [3.0, 1.0, 4.0, 2.0]
        for mol, e in zip(mols, energies):
            mol.energy = e
        job = GaussianTrajJob(
            molecules=mols,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=1.0,
        )
        assert job.sorted_energies == [1.0, 2.0, 3.0, 4.0]
        assert [m.energy for m in job.sorted_molecules] == [
            1.0,
            2.0,
            3.0,
            4.0,
        ]


class TestGaussianTrajJobUniqueStructures:
    def test_no_grouper_returns_all_molecules(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        assert job.unique_structures == job.molecules
        assert job.num_unique_structures == job.num_structures

    def test_grouping_strategy_delegates_to_grouper(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        unique_subset = ethanol_traj_molecules[-2:][:1]
        mock_grouper = MagicMock()
        mock_grouper.unique.return_value = unique_subset

        with patch(
            "chemsmart.utils.grouper.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            job = GaussianTrajJob(
                molecules=ethanol_traj_molecules,
                settings=gaussian_settings,
                label="traj",
                jobrunner=None,
                grouping_strategy="rmsd",
            )

        mock_create.assert_called_once()
        _, kwargs = mock_create.call_args
        assert kwargs["strategy"] == "rmsd"
        assert kwargs["conformer_ids"] == job._original_conformer_ids
        mock_grouper.group.assert_called_once()
        assert job.grouper is mock_grouper
        assert job.unique_structures == unique_subset
        assert job.num_unique_structures == 1
        assert job.unique_structures_energies == [
            m.energy for m in unique_subset
        ]


class TestGaussianTrajJobPrepareAndRun:
    def _job_with_fake_uniques(
        self, ethanol_traj_molecules, gaussian_settings, n_unique
    ):
        job = GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        fake_unique = job.molecules[:n_unique]
        job.grouper = MagicMock()
        job.grouper.unique.return_value = fake_unique
        return job

    def test_prepare_all_jobs_builds_indexed_labels(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            ethanol_traj_molecules, gaussian_settings, n_unique=2
        )
        jobs = job._prepare_all_jobs()
        assert [j.label for j in jobs] == ["traj_c1", "traj_c2"]

    def test_run_all_jobs_runs_everything_when_no_limit(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            ethanol_traj_molecules, gaussian_settings, n_unique=2
        )
        fake_jobs = [MagicMock(), MagicMock()]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            job._run_all_jobs()
        for fj in fake_jobs:
            fj.run.assert_called_once()

    def test_run_all_jobs_limits_to_incomplete_subset(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            ethanol_traj_molecules, gaussian_settings, n_unique=2
        )
        job.num_structures_to_run = 1
        complete_job = MagicMock()
        complete_job.is_complete.return_value = True
        incomplete_job_1 = MagicMock()
        incomplete_job_1.is_complete.return_value = False
        incomplete_job_2 = MagicMock()
        incomplete_job_2.is_complete.return_value = False
        fake_jobs = [complete_job, incomplete_job_1, incomplete_job_2]

        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            job._run_all_jobs()

        complete_job.run.assert_not_called()
        incomplete_job_1.run.assert_called_once()
        incomplete_job_2.run.assert_not_called()

    def test_run_delegates_to_run_all_jobs(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            ethanol_traj_molecules, gaussian_settings, n_unique=2
        )
        with patch.object(job, "_run_all_jobs") as mock_run_all:
            job._run()
        mock_run_all.assert_called_once()


class TestGaussianTrajJobCompletion:
    def _job(self, ethanol_traj_molecules, gaussian_settings):
        return GaussianTrajJob(
            molecules=ethanol_traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )

    def test_is_complete_checks_all_jobs_when_no_limit(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job(ethanol_traj_molecules, gaussian_settings)
        fake_jobs = [
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: True),
        ]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.is_complete() is True

    def test_is_complete_false_when_a_required_job_incomplete(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job(ethanol_traj_molecules, gaussian_settings)
        job.num_structures_to_run = 1
        fake_jobs = [
            MagicMock(is_complete=lambda: False),
            MagicMock(is_complete=lambda: True),
        ]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.is_complete() is False

    def test_last_run_job_index_returns_first_incomplete(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job(ethanol_traj_molecules, gaussian_settings)
        fake_jobs = [
            MagicMock(is_complete=lambda: True),
            MagicMock(is_complete=lambda: False),
        ]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.last_run_job_index == 1

    def test_incomplete_structure_run_jobs_filters(
        self, ethanol_traj_molecules, gaussian_settings
    ):
        job = self._job(ethanol_traj_molecules, gaussian_settings)
        complete = MagicMock(is_complete=lambda: True)
        incomplete = MagicMock(is_complete=lambda: False)
        fake_jobs = [complete, incomplete]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.incomplete_structure_run_jobs == [incomplete]


class TestGaussianDIASJobConstruction:
    def test_uses_first_molecule_as_placeholder(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        assert job.num_molecules == 5
        assert job.all_molecules == ethanol_irc_molecules

    def test_freq_disabled_on_settings(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        gaussian_settings.freq = True
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        assert job.settings.freq is False

    def test_fragment_charge_and_multiplicity_overrides(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
            charge_of_fragment1=1,
            multiplicity_of_fragment1=2,
            charge_of_fragment2=-1,
            multiplicity_of_fragment2=3,
        )
        assert job.fragment1_settings.charge == 1
        assert job.fragment1_settings.multiplicity == 2
        assert job.fragment2_settings.charge == -1
        assert job.fragment2_settings.multiplicity == 3

    def test_fragment_settings_default_to_base_charge(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        gaussian_settings.charge = 0
        gaussian_settings.multiplicity = 1
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        assert job.fragment1_settings.charge == 0
        assert job.fragment2_settings.multiplicity == 1


class TestFragmentStructure:
    def test_fragment_structure_splits_by_indices(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        f1, f2 = job._fragment_structure(ethanol_irc_molecules[0])
        assert len(f1) == 3
        assert len(f2) == 6
        assert len(f1) + len(f2) == len(ethanol_irc_molecules[0])

    def test_fragment1_atoms_and_fragment2_atoms_properties(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        assert len(job.fragment1_atoms) == len(ethanol_irc_molecules)
        assert len(job.fragment2_atoms) == len(ethanol_irc_molecules)
        assert all(len(f1) == 3 for f1 in job.fragment1_atoms)
        assert all(len(f2) == 6 for f2 in job.fragment2_atoms)


class TestSampleMolecules:
    def test_samples_every_n_points_and_appends_last_again(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        """``_sample_molecules`` slices every ``every_n_points``-th
        molecule (``[0::n]``), then *unconditionally* appends the last
        molecule again whenever ``(num_molecules - 1) / every_n_points
        != 0`` — which is true for essentially any real trajectory. So
        for 5 molecules sampled every 2 points, the last molecule
        (index 4) appears twice: once from the slice ``[0, 2, 4]`` and
        once more from the unconditional append."""
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        sampled = job._sample_molecules(ethanol_irc_molecules)
        assert sampled == [
            ethanol_irc_molecules[0],
            ethanol_irc_molecules[2],
            ethanol_irc_molecules[4],
            ethanol_irc_molecules[4],
        ]


class TestFragmentJobsIRCMode:
    def test_all_molecules_jobs_irc_mode_creates_one_per_sample(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        jobs = job.all_molecules_jobs
        # 4 sampled points due to the duplicate-last-molecule behavior
        # documented in TestSampleMolecules.
        assert len(jobs) == 4
        labels = [j.label for j in jobs]
        assert labels == [
            "dias_test_p0",
            "dias_test_p1",
            "dias_test_p2",
            "dias_test_p3",
        ]

    def test_fragment1_jobs_irc_mode(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        jobs = job.fragment1_jobs
        assert len(jobs) == 4
        assert jobs[0].label == "dias_test_p0_f1"

    def test_fragment2_jobs_irc_mode(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        jobs = job.fragment2_jobs
        assert len(jobs) == 4
        assert jobs[0].label == "dias_test_p0_f2"


class TestFragmentJobsTSMode:
    def test_all_molecules_jobs_ts_mode_uses_last_only(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        jobs = job.all_molecules_jobs
        assert len(jobs) == 1
        assert jobs[0].label == "dias_test_p1"

    def test_fragment1_jobs_ts_mode(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        jobs = job.fragment1_jobs
        assert len(jobs) == 1
        assert jobs[0].label == "dias_test_p1_f1"

    def test_fragment2_jobs_ts_mode(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=mock_jobrunner,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        jobs = job.fragment2_jobs
        assert len(jobs) == 1
        assert jobs[0].label == "dias_test_p1_f2"


class TestFragmentJobsInvalidMode:
    def test_fragment1_jobs_invalid_mode_raises(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="bogus",
        )
        with pytest.raises(ValueError, match="Invalid mode"):
            _ = job.fragment1_jobs

    def test_fragment2_jobs_invalid_mode_raises(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="bogus",
        )
        with pytest.raises(ValueError, match="Invalid mode"):
            _ = job.fragment2_jobs


class TestRunAndIsComplete:
    def test_run_executes_all_three_job_sets_in_order(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        job._run_all_molecules_jobs = MagicMock()
        job._run_fragment1_jobs = MagicMock()
        job._run_fragment2_jobs = MagicMock()

        job._run()

        job._run_all_molecules_jobs.assert_called_once()
        job._run_fragment1_jobs.assert_called_once()
        job._run_fragment2_jobs.assert_called_once()

    def test_is_complete_true_when_all_job_sets_complete(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        job._run_all_molecules_jobs_are_complete = MagicMock(return_value=True)
        job._run_fragment1_jobs_are_complete = MagicMock(return_value=True)
        job._run_fragment2_jobs_are_complete = MagicMock(return_value=True)

        assert job.is_complete() is True

    def test_is_complete_false_when_one_set_incomplete(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        job._run_all_molecules_jobs_are_complete = MagicMock(return_value=True)
        job._run_fragment1_jobs_are_complete = MagicMock(return_value=False)
        job._run_fragment2_jobs_are_complete = MagicMock(return_value=True)

        assert job.is_complete() is False


class TestRunHelpersAndCompletionChecksDirect:
    """The tests above mock out _run_*_jobs/_run_*_jobs_are_complete
    entirely, so their real bodies (each a one-line loop or `all(...)`
    call) are never exercised. These call the real methods directly,
    mocking only GaussianGeneralJob so no actual calculation runs."""

    def _make_job(self, ethanol_irc_molecules, gaussian_settings, mode="ts"):
        return GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode=mode,
        )

    def test_run_all_molecules_jobs_runs_each_job(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            job._run_all_molecules_jobs()
        mock_job.run.assert_called_once()

    def test_run_fragment1_jobs_runs_each_job(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            job._run_fragment1_jobs()
        mock_job.run.assert_called_once()

    def test_run_fragment2_jobs_runs_each_job(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            job._run_fragment2_jobs()
        mock_job.run.assert_called_once()

    def test_run_all_molecules_jobs_are_complete_reflects_job_state(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        mock_job.is_complete.return_value = True
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            assert job._run_all_molecules_jobs_are_complete() is True

    def test_run_fragment1_jobs_are_complete_reflects_job_state(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        mock_job.is_complete.return_value = False
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            assert job._run_fragment1_jobs_are_complete() is False

    def test_run_fragment2_jobs_are_complete_reflects_job_state(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        job = self._make_job(ethanol_irc_molecules, gaussian_settings)
        mock_job = MagicMock()
        mock_job.is_complete.return_value = True
        with patch(
            "chemsmart.jobs.gaussian.dias.GaussianGeneralJob",
            return_value=mock_job,
        ):
            assert job._run_fragment2_jobs_are_complete() is True


class TestSampleMoleculesSingleMolecule:
    def test_single_molecule_does_not_duplicate_last(
        self, ethanol_molecule, gaussian_settings
    ):
        """With only one molecule, (num_molecules - 1) / every_n_points
        is 0, so the "always include the last point" append must be
        skipped -- otherwise the single molecule would be duplicated."""
        job = GaussianDIASJob(
            molecules=[ethanol_molecule],
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        sampled = job._sample_molecules([ethanol_molecule])
        assert sampled == [ethanol_molecule]


class TestAllMoleculesJobsInvalidMode:
    def test_invalid_mode_silently_returns_none(
        self, ethanol_irc_molecules, gaussian_settings
    ):
        """Unlike fragment1_jobs/fragment2_jobs, which raise
        ValueError for an unrecognized mode, all_molecules_jobs has no
        matching `else` clause and silently falls through, returning
        None instead of raising -- an inconsistency worth flagging
        (see BUGS_FOUND.md)."""
        job = GaussianDIASJob(
            molecules=ethanol_irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=MagicMock(),
            fragment_indices="1-3",
            every_n_points=2,
            mode="bogus",
        )
        assert job.all_molecules_jobs is None


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestInitValidation:
    def test_rejects_non_pka_settings(self, ethanol_molecule):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        with pytest.raises(ValueError, match="GaussianpKaJobSettings"):
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=GaussianJobSettings.default(),
                label="pka_test",
            )


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestPrepareJobsNoReference:
    def test_creates_only_target_opt_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        opt_factory, _ = patched_sub_jobs
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings,
            label="pka_test",
        )
        assert job.has_reference_jobs is False
        assert job.ref_opt_jobs == []
        assert job.ref_acid_job is None
        assert job.ref_conjugate_base_job is None
        assert len(opt_factory.created) == 2
        assert job.protonated_job.label == "pka_test_HA_opt"
        assert job.conjugate_base_job.label == "pka_test_A_opt"
        assert job.opt_jobs == [job.protonated_job, job.conjugate_base_job]


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestPrepareJobsWithReference:
    def test_creates_reference_opt_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        opt_factory, _ = patched_sub_jobs
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.has_reference_jobs is True
        assert len(opt_factory.created) == 4
        assert job.ref_acid_job.label == "pka_test_HRef_opt"
        assert job.ref_conjugate_base_job.label == "pka_test_Ref_opt"
        assert job.ref_opt_jobs == [
            job.ref_acid_job,
            job.ref_conjugate_base_job,
        ]

    def test_reference_pair_molecules_cached_across_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        with patch.object(
            _FakePKaSettings,
            "reference_pair_molecules",
            return_value=("REF_ACID_MOL", "REF_CONJ_MOL"),
        ) as mocked:
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference,
                label="pka_job_1",
            )
            GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings_with_reference.copy(),
                label="pka_job_2",
            )
        # Same cache key (scheme, reference_file, indices/charges) ->
        # the expensive geometry-building call only happens once.
        assert mocked.call_count == 1


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestReferenceCacheHelpers:
    def test_reference_cache_key_none_for_none_settings(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._reference_cache_key(None) is None

    def test_reference_cache_key_none_without_reference_file(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._reference_cache_key(pka_settings) is None

    def test_get_cached_reference_pair_none_without_reference_file(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._get_cached_reference_pair(pka_settings) is None

    def test_prepare_pka_jobs_noop_when_settings_none(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.settings = None
        job._prepare_pka_jobs()  # must not raise

    def test_run_ref_opt_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_opt_jobs()  # must not raise; nothing to run

    def test_create_ref_sp_jobs_falls_back_when_cache_miss(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            job._create_ref_sp_jobs()
        assert job.ref_sp_jobs is not None

    def test_prepare_ref_opt_jobs_falls_back_when_cache_miss(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        with patch.object(
            job, "_get_cached_reference_pair", return_value=None
        ):
            ref_acid_job, ref_conjugate_base_job = job._prepare_ref_opt_jobs()
        assert ref_acid_job.label == "pka_test_HRef_opt"
        assert ref_conjugate_base_job.label == "pka_test_Ref_opt"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestOptimizedMoleculeFromJob:
    def test_returns_output_molecule_on_normal_termination(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        fake_output = MagicMock(normal_termination=True, molecule="OPT_MOL")
        child_job = MagicMock()
        child_job._output.return_value = fake_output
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "OPT_MOL"

    def test_falls_back_when_not_normally_terminated(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        fake_output = MagicMock(normal_termination=False, molecule="OPT_MOL")
        child_job = MagicMock()
        child_job._output.return_value = fake_output
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "FALLBACK"

    def test_falls_back_when_no_output(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        child_job = MagicMock()
        child_job._output.return_value = None
        result = job._optimized_molecule_from_job(child_job, "FALLBACK")
        assert result == "FALLBACK"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestCompletionTracking:
    def test_opt_jobs_are_complete_false_when_empty(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.opt_jobs = []
        assert job._opt_jobs_are_complete() is False

    def test_opt_jobs_are_complete_true_when_all_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        assert job._opt_jobs_are_complete() is True

    def test_ref_opt_jobs_are_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_opt_jobs_are_complete() is True

    def test_ref_opt_jobs_are_complete_false_when_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job._ref_opt_jobs_are_complete() is False
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        assert job._ref_opt_jobs_are_complete() is True

    def test_sp_jobs_are_complete_false_when_none(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._sp_jobs_are_complete() is False

    def test_ref_sp_jobs_are_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job._ref_sp_jobs_are_complete() is True

    def test_is_complete_requires_all_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
        sp_factory,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.is_complete() is False

        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        assert job.is_complete() is False  # SP jobs not even created yet

        job._create_sp_jobs()
        for sp_job in job.sp_jobs:
            sp_job.is_complete.return_value = True
        assert job.is_complete() is False  # reference opt still incomplete

        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        assert job.is_complete() is False  # reference SP not created yet

        job._create_ref_sp_jobs()
        for ref_sp_job in job.ref_sp_jobs:
            ref_sp_job.is_complete.return_value = True
        assert job.is_complete() is True

    def test_is_complete_true_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()
        for sp_job in job.sp_jobs:
            sp_job.is_complete.return_value = True
        assert job.is_complete() is True


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestRunOptAndSpJobs:
    def test_run_sp_jobs_warns_and_skips_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_sp_jobs()
        assert job.sp_jobs is None

    def test_run_sp_jobs_creates_and_runs_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, sp_factory
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._run_sp_jobs()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()

    def test_run_ref_sp_jobs_noop_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is None

    def test_run_ref_sp_jobs_warns_when_ref_opt_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is None

    def test_run_ref_sp_jobs_creates_and_runs_when_ref_opt_complete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._run_ref_sp_jobs()
        assert job.ref_sp_jobs is not None
        for ref_sp_job in job.ref_sp_jobs:
            ref_sp_job.run.assert_called_once()

    def test_run_sp_jobs_reuses_existing_sp_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._run_sp_jobs()
        existing_sp_jobs = job.sp_jobs

        job._run_sp_jobs()  # second call must not recreate sp_jobs
        assert job.sp_jobs is existing_sp_jobs
        for sp_job in existing_sp_jobs:
            assert sp_job.run.call_count == 2

    def test_run_ref_sp_jobs_reuses_existing_ref_sp_jobs(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._run_ref_sp_jobs()
        existing_ref_sp_jobs = job.ref_sp_jobs

        job._run_ref_sp_jobs()  # second call must not recreate ref_sp_jobs
        assert job.ref_sp_jobs is existing_ref_sp_jobs
        for ref_sp_job in existing_ref_sp_jobs:
            assert ref_sp_job.run.call_count == 2


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestRunPhaseStateMachine:
    def test_stops_after_incomplete_opt_phase(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # Both opt jobs stay incomplete after running (clear the
        # factory's default auto-complete-on-run behaviour).
        for opt_job in job.opt_jobs:
            opt_job.run.side_effect = None
        job._run()
        # `_execute_phase_jobs` stops the loop as soon as one job is
        # found incomplete, so only the first job actually ran.
        job.opt_jobs[0].run.assert_called_once()
        job.opt_jobs[1].run.assert_not_called()
        assert job.sp_jobs is None

    def test_runs_sp_phase_when_opt_completes(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        # Default factory behaviour: each job completes once run().
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        assert job.sp_jobs is not None
        for sp_job in job.sp_jobs:
            sp_job.run.assert_called_once()

    def test_reference_opt_failure_blocks_target_sp_phase(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        # Target opt jobs succeed (default behaviour); reference opt
        # jobs are run but never mark themselves complete.
        for ref_job in job.ref_opt_jobs:
            ref_job.run.side_effect = None
        job._run()
        for opt_job in job.opt_jobs:
            opt_job.run.assert_called_once()
        job.ref_opt_jobs[0].run.assert_called_once()
        # Target SP jobs must not have been created/run: reference opt
        # incompleteness halts the whole workflow before that phase.
        assert job.sp_jobs is None

    def test_full_run_with_reference_completes_all_phases(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job._run()
        assert job.is_complete() is True

    def test_stops_when_sp_phase_fails(self, ethanol_molecule, pka_settings):
        opt_factory = _make_job_factory()

        def _sp_factory_never_completes(
            *, molecule, settings, label, jobrunner, skip_completed
        ):
            job = MagicMock()
            job.label = label
            job.outputfile = f"{label}.log"
            job.is_complete.return_value = False
            job._output.return_value = None
            return job

        with (
            patch(
                "chemsmart.jobs.gaussian.pka.GaussianOptJob",
                side_effect=opt_factory,
            ),
            patch(
                "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob",
                side_effect=_sp_factory_never_completes,
            ),
        ):
            job = GaussianpKaJob(
                molecule=ethanol_molecule,
                settings=pka_settings,
                label="pka_test",
            )
            job._run()

        assert job.sp_jobs is not None
        # `_execute_phase_jobs` stops as soon as one job is found
        # incomplete, so only the first SP job actually ran.
        job.sp_jobs[0].run.assert_called_once()
        job.sp_jobs[1].run.assert_not_called()
        assert job._sp_jobs_are_complete() is False


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestOutputProperties:
    def test_target_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.protonated_job._output.return_value = "HA_OUT"
        job.conjugate_base_job._output.return_value = "A_OUT"
        assert job.protonated_output == "HA_OUT"
        assert job.conjugate_base_output == "A_OUT"

    def test_reference_output_properties_none_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.ref_acid_output is None
        assert job.ref_conjugate_base_output is None
        assert job.ref_acid_sp_output is None
        assert job.ref_conjugate_base_sp_output is None

    def test_reference_output_properties_delegate_when_present(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job.ref_acid_job._output.return_value = "HREF_OUT"
        job.ref_conjugate_base_job._output.return_value = "REF_OUT"
        assert job.ref_acid_output == "HREF_OUT"
        assert job.ref_conjugate_base_output == "REF_OUT"

    def test_target_sp_output_properties_delegate_to_jobs(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()
        job.protonated_sp_job._output.return_value = "HA_SP_OUT"
        job.conjugate_base_sp_job._output.return_value = "A_SP_OUT"
        assert job.protonated_sp_output == "HA_SP_OUT"
        assert job.conjugate_base_sp_output == "A_SP_OUT"

    def test_reference_sp_output_properties_delegate_when_present(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._create_ref_sp_jobs()
        job.ref_acid_sp_job._output.return_value = "HREF_SP_OUT"
        job.ref_conjugate_base_sp_job._output.return_value = "REF_SP_OUT"
        assert job.ref_acid_sp_output == "HREF_SP_OUT"
        assert job.ref_conjugate_base_sp_output == "REF_SP_OUT"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestCompletionCheckEdgeCases:
    def test_ref_opt_jobs_are_complete_false_when_list_empty(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        job.ref_opt_jobs = []
        assert job._ref_opt_jobs_are_complete() is False

    def test_ref_sp_jobs_are_complete_false_when_none(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        assert job.ref_sp_jobs is None
        assert job._ref_sp_jobs_are_complete() is False


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestPkaOutputFiles:
    def test_excludes_reference_paths_without_reference(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        ha_file, a_file, href_file, ref_file = job._pka_output_files()
        assert ha_file == "pka_test_HA_opt.log"
        assert a_file == "pka_test_A_opt.log"
        assert href_file is None
        assert ref_file is None

    def test_excludes_reference_paths_when_ref_opt_incomplete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        _, _, href_file, ref_file = job._pka_output_files()
        assert href_file is None
        assert ref_file is None

    def test_includes_reference_paths_when_ref_opt_complete(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        _, _, href_file, ref_file = job._pka_output_files()
        assert href_file == "pka_test_HRef_opt.log"
        assert ref_file == "pka_test_Ref_opt.log"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestThermochemistryDelegation:
    def test_compute_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        with pytest.raises(ValueError, match="not complete"):
            job.compute_thermochemistry()

    def test_compute_thermochemistry_delegates_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".compute_pka_thermochemistry",
            return_value={"pKa": 4.2},
        ) as mocked:
            result = job.compute_thermochemistry()

        assert result == {"pKa": 4.2}
        mocked.assert_called_once()
        _, kwargs = mocked.call_args
        assert kwargs["ha_file"] == "pka_test_HA_opt.log"
        assert kwargs["a_file"] == "pka_test_A_opt.log"
        assert kwargs["href_file"] is None
        assert kwargs["temperature"] == pka_settings.temperature

    def test_print_thermochemistry_raises_when_opt_incomplete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        with pytest.raises(ValueError, match="not complete"):
            job.print_thermochemistry()

    def test_print_thermochemistry_delegates_when_opt_complete(
        self, ethanol_molecule, pka_settings, patched_sub_jobs, sp_factory
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        job._create_sp_jobs()

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        mocked.assert_called_once()
        _, kwargs = mocked.call_args
        assert kwargs["ha_gas_file"] == "pka_test_HA_opt.log"
        assert kwargs["ha_solv_file"] == "pka_test_HA_sp.log"
        assert kwargs["href_solv_file"] is None

    def test_print_thermochemistry_includes_reference_solv_files(
        self,
        ethanol_molecule,
        pka_settings_with_reference,
        patched_sub_jobs,
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule,
            settings=pka_settings_with_reference,
            label="pka_test",
        )
        for opt_job in job.opt_jobs:
            opt_job.is_complete.return_value = True
        for ref_job in job.ref_opt_jobs:
            ref_job.is_complete.return_value = True
        job._create_sp_jobs()
        job._create_ref_sp_jobs()

        with patch(
            "chemsmart.io.gaussian.output.Gaussian16pKaOutput"
            ".print_pka_summary"
        ) as mocked:
            job.print_thermochemistry()

        _, kwargs = mocked.call_args
        assert kwargs["href_solv_file"] == "pka_test_HRef_sp.log"
        assert kwargs["ref_solv_file"] == "pka_test_Ref_sp.log"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestMoleculeProperties:
    def test_original_mol_is_input_molecule(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        assert job.original_mol is job.molecule

    def test_conjugate_base_mol_from_job(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.conjugate_base_job.molecule = "CONJ_FROM_JOB"
        assert job.conjugate_base_mol == "CONJ_FROM_JOB"

    def test_conjugate_base_mol_falls_back_without_job(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaJob(
            molecule=ethanol_molecule, settings=pka_settings, label="pka_test"
        )
        job.conjugate_base_job = None
        assert job.conjugate_base_mol == "CONJ_MOL"


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestGaussianpKaAnalyzeJob:
    def test_run_calls_print_thermochemistry(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaAnalyzeJob(
            input_file=ethanol_molecule,
            settings=pka_settings,
            label="pka_analyze",
        )
        with patch.object(job, "print_thermochemistry") as mocked:
            job._run()
        mocked.assert_called_once()

    def test_run_swallows_exceptions(
        self, ethanol_molecule, pka_settings, patched_sub_jobs
    ):
        job = GaussianpKaAnalyzeJob(
            input_file=ethanol_molecule,
            settings=pka_settings,
            label="pka_analyze",
        )
        with patch.object(
            job, "print_thermochemistry", side_effect=RuntimeError("boom")
        ):
            job._run()  # must not raise


@pytest.mark.usefixtures("_clear_shared_reference_cache")
class TestGaussianpKaThermoJob:
    def test_type_identifier(self):
        assert GaussianpKaThermoJob.TYPE == "g16pka_thermo"
        assert issubclass(GaussianpKaThermoJob, GaussianpKaAnalyzeJob)


class TestGaussianUVVISJobInitCrash:
    def test_init_always_raises_typeerror(self):
        """Documents a bug (see BUGS_FOUND.md #11): __init__ forwards
        ``atoms=`` to ``GaussianJob.__init__``, whose signature requires
        ``molecule=`` (with no default), so construction always fails."""
        settings = GaussianJobSettings.default()
        molecule = MagicMock(name="molecule")
        with pytest.raises(
            TypeError, match="missing 1 required positional argument"
        ):
            GaussianUVVISJob(
                folder="/tmp/uvvis_job", atoms=molecule, settings=settings
            )


class TestGsGeomOpt:
    def test_returns_gaussian_opt_job_with_copied_settings(self):
        job = _make_job()
        with patch(
            "chemsmart.jobs.gaussian.uvvis.GaussianOptJob"
        ) as mock_opt_cls:
            mock_opt_cls.return_value = "opt-job"
            result = job._gs_geom_opt()

        assert result == "opt-job"
        mock_opt_cls.assert_called_once()
        _, kwargs = mock_opt_cls.call_args
        assert kwargs["folder"] == job.folder
        assert kwargs["atoms"] == job.molecule
        # settings passed in is a deepcopy, not the same object
        assert kwargs["settings"] is not job.settings

    def test_run_gs_geom_opt_delegates_to_run(self):
        job = _make_job()
        mock_opt_job = MagicMock()
        with patch.object(job, "_gs_geom_opt", return_value=mock_opt_job):
            job._run_gs_geom_opt(jobrunner="runner", queue_manager="qm")

        mock_opt_job.run.assert_called_once_with(
            jobrunner="runner", queue_manager="qm"
        )

    def test_gs_geom_opt_complete_delegates_to_is_complete(self):
        job = _make_job()
        mock_opt_job = MagicMock()
        mock_opt_job.is_complete.return_value = True
        with patch.object(job, "_gs_geom_opt", return_value=mock_opt_job):
            assert job._gs_geom_opt_complete() is True
        mock_opt_job.is_complete.assert_called_once()


class TestVerticalExcitation:
    def test_returns_none_when_gs_not_complete(self):
        job = _make_job()
        with patch.object(job, "_gs_geom_opt_complete", return_value=False):
            assert job._vertical_excitation() is None

    def test_returns_sp_job_with_jobtype_sp_and_freq_false_when_complete(
        self,
    ):
        job = _make_job()
        job.settings.jobtype = "opt"
        job.settings.freq = True

        with (
            patch.object(job, "_gs_geom_opt_complete", return_value=True),
            patch(
                "chemsmart.jobs.gaussian.uvvis.GaussianSinglePointJob"
            ) as mock_sp_cls,
        ):
            mock_sp_cls.return_value = "sp-job"
            result = job._vertical_excitation()

        assert result == "sp-job"
        _, kwargs = mock_sp_cls.call_args
        assert kwargs["folder"] == job.folder
        assert kwargs["atoms"] == job.molecule
        sp_settings = kwargs["settings"]
        assert sp_settings.jobtype == "sp"
        assert sp_settings.freq is False
        # original settings untouched (deepcopy, not aliased)
        assert job.settings.jobtype == "opt"
        assert job.settings.freq is True

    def test_run_vertical_excitation_runs_then_reports_completion(self):
        job = _make_job()
        mock_sp_job = MagicMock()
        mock_sp_job.is_complete.return_value = True
        with patch.object(
            job, "_vertical_excitation", return_value=mock_sp_job
        ):
            result = job._run_vertical_excitation(
                jobrunner="runner", queue_manager="qm"
            )

        mock_sp_job.run.assert_called_once_with(
            jobrunner="runner", queue_manager="qm"
        )
        mock_sp_job.is_complete.assert_called_once()
        assert result is True


class TestRunWorkflow:
    def test_run_executes_both_steps_in_order(self):
        job = _make_job()
        calls = []

        def record_gs(jobrunner, queue_manager=None):
            calls.append(("gs", jobrunner, queue_manager))

        def record_ve(jobrunner, queue_manager=None):
            calls.append(("ve", jobrunner, queue_manager))

        with (
            patch.object(job, "_run_gs_geom_opt", side_effect=record_gs),
            patch.object(
                job, "_run_vertical_excitation", side_effect=record_ve
            ),
        ):
            job._run(jobrunner="runner", queue_manager="qm")

        assert calls == [
            ("gs", "runner", "qm"),
            ("ve", "runner", "qm"),
        ]


class TestGaussianWBIJob:
    def test_type_identifier(self):
        assert GaussianWBIJob.TYPE == "g16wbi"

    def test_freq_disabled_on_settings(
        self, ethanol_molecule, gaussian_settings
    ):
        gaussian_settings.freq = True
        job = GaussianWBIJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="wbi_test",
            jobrunner=None,
        )
        assert job.settings.freq is False

    def test_jobrunner_is_respected(self, ethanol_molecule, gaussian_settings):
        sentinel_runner = object()
        job = GaussianWBIJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="wbi_test",
            jobrunner=sentinel_runner,
        )
        assert job.jobrunner is sentinel_runner


class TestGaussianNCIJob:
    def test_type_identifier(self):
        assert GaussianNCIJob.TYPE == "g16nci"

    def test_freq_disabled_on_settings(
        self, ethanol_molecule, gaussian_settings
    ):
        gaussian_settings.freq = True
        job = GaussianNCIJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="nci_test",
        )
        assert job.settings.freq is False

    def test_jobrunner_argument_is_ignored_and_class_used_instead(
        self, ethanol_molecule, gaussian_settings
    ):
        """Documents a bug (see BUGS_FOUND.md #16): __init__ hardcodes
        ``jobrunner=JobRunner`` (the class itself) when calling
        ``super().__init__``, instead of forwarding the ``jobrunner``
        parameter it was actually given. Any caller-supplied jobrunner
        instance is silently discarded."""
        sentinel_runner = object()
        job = GaussianNCIJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="nci_test",
            jobrunner=sentinel_runner,
        )
        assert job.jobrunner is not sentinel_runner
        assert job.jobrunner is JobRunner


class TestGaussianRESPJob:
    def test_type_identifier(self):
        assert GaussianRESPJob.TYPE == "g16resp"

    def test_jobrunner_is_respected(self, ethanol_molecule, gaussian_settings):
        sentinel_runner = object()
        job = GaussianRESPJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="resp_test",
            jobrunner=sentinel_runner,
        )
        assert job.jobrunner is sentinel_runner

    def test_settings_not_mutated(self, ethanol_molecule, gaussian_settings):
        gaussian_settings.freq = True
        job = GaussianRESPJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="resp_test",
        )
        # unlike WBI/NCI, RESP does not force-disable freq
        assert job.settings.freq is True


class TestGaussianCustomJob:
    """The CLI tests in test_gaussian_cli.py mock GaussianCustomJob
    entirely, so its real __init__ (a thin pass-through to GaussianJob)
    was never actually exercised."""

    def test_type_identifier(self):
        assert GaussianCustomJob.TYPE == "g16job"

    def test_jobrunner_is_respected(self, ethanol_molecule, gaussian_settings):
        sentinel_runner = object()
        job = GaussianCustomJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="custom_test",
            jobrunner=sentinel_runner,
        )
        assert job.jobrunner is sentinel_runner

    def test_settings_not_mutated(self, ethanol_molecule, gaussian_settings):
        gaussian_settings.freq = True
        job = GaussianCustomJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="custom_test",
        )
        assert job.settings.freq is True


class TestGaussianTDDFTJob:
    def test_type_identifier(self):
        assert GaussianTDDFTJob.TYPE == "g16td"

    def test_settings_class_is_tddft_specific(self):
        assert GaussianTDDFTJob.settings_class() is GaussianTDDFTJobSettings

    def test_jobrunner_is_respected(self, ethanol_molecule, gaussian_settings):
        sentinel_runner = object()
        job = GaussianTDDFTJob(
            molecule=ethanol_molecule,
            settings=gaussian_settings,
            label="td_test",
            jobrunner=sentinel_runner,
        )
        assert job.jobrunner is sentinel_runner
