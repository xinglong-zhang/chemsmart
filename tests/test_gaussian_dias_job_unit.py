"""
Direct unit tests for :class:`GaussianDIASJob` in
``chemsmart.jobs.gaussian.dias``.

CLI-level invocation (through ``chemsmart sub gaussian ... dias``) is
covered by ``TestGaussianCLIDiasCommand`` in ``test_gaussian_cli.py``,
but that mocks the job class entirely. These tests exercise the class
directly: fragment splitting, IRC-mode point sampling, TS-mode endpoint
selection, job list generation, and run/is_complete orchestration.
"""

from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.gaussian.dias import GaussianDIASJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


@pytest.fixture()
def irc_molecules(ethanol_molecule):
    # A short synthetic "IRC trajectory": 5 copies of the same
    # 9-atom ethanol geometry, standing in for points along a path.
    return [ethanol_molecule.copy() for _ in range(5)]


class TestGaussianDIASJobConstruction:
    def test_uses_first_molecule_as_placeholder(
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="irc",
        )
        assert job.num_molecules == 5
        assert job.all_molecules == irc_molecules

    def test_freq_disabled_on_settings(self, irc_molecules, gaussian_settings):
        gaussian_settings.freq = True
        job = GaussianDIASJob(
            molecules=irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        assert job.settings.freq is False

    def test_fragment_charge_and_multiplicity_overrides(
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        gaussian_settings.charge = 0
        gaussian_settings.multiplicity = 1
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        f1, f2 = job._fragment_structure(irc_molecules[0])
        assert len(f1) == 3
        assert len(f2) == 6
        assert len(f1) + len(f2) == len(irc_molecules[0])

    def test_fragment1_atoms_and_fragment2_atoms_properties(
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        assert len(job.fragment1_atoms) == len(irc_molecules)
        assert len(job.fragment2_atoms) == len(irc_molecules)
        assert all(len(f1) == 3 for f1 in job.fragment1_atoms)
        assert all(len(f2) == 6 for f2 in job.fragment2_atoms)


class TestSampleMolecules:
    def test_samples_every_n_points_and_appends_last_again(
        self, irc_molecules, gaussian_settings
    ):
        """``_sample_molecules`` slices every ``every_n_points``-th
        molecule (``[0::n]``), then *unconditionally* appends the last
        molecule again whenever ``(num_molecules - 1) / every_n_points
        != 0`` — which is true for essentially any real trajectory. So
        for 5 molecules sampled every 2 points, the last molecule
        (index 4) appears twice: once from the slice ``[0, 2, 4]`` and
        once more from the unconditional append."""
        job = GaussianDIASJob(
            molecules=irc_molecules,
            settings=gaussian_settings,
            label="dias_test",
            jobrunner=None,
            fragment_indices="1-3",
            every_n_points=2,
            mode="ts",
        )
        sampled = job._sample_molecules(irc_molecules)
        assert sampled == [
            irc_molecules[0],
            irc_molecules[2],
            irc_molecules[4],
            irc_molecules[4],
        ]


class TestFragmentJobsIRCMode:
    def test_all_molecules_jobs_irc_mode_creates_one_per_sample(
        self, irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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

    def test_fragment1_jobs_irc_mode(self, irc_molecules, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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

    def test_fragment2_jobs_irc_mode(self, irc_molecules, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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

    def test_fragment1_jobs_ts_mode(self, irc_molecules, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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

    def test_fragment2_jobs_ts_mode(self, irc_molecules, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
        self, irc_molecules, gaussian_settings
    ):
        job = GaussianDIASJob(
            molecules=irc_molecules,
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
