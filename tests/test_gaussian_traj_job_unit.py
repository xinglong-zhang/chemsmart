"""
Direct unit tests for :class:`GaussianTrajJob` in
``chemsmart.jobs.gaussian.traj``.

Covers input validation, the tail-fraction extraction/clamping logic
for ``proportion_structures_to_use``, energy sorting, grouping
delegation, job preparation, and run/is_complete orchestration.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.traj import GaussianTrajJob


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


@pytest.fixture()
def traj_molecules(ethanol_molecule):
    mols = []
    for i in range(20):
        mol = ethanol_molecule.copy()
        mol.energy = float(i)
        mols.append(mol)
    return mols


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
        self, traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        # 20 * 0.1 = 2.0 -> last 2 structures
        assert job.num_structures == 2
        assert job.molecules == traj_molecules[-2:]

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
        self, traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=1.0,
        )
        assert job.num_structures == 20
        assert job.molecules == traj_molecules

    def test_original_conformer_ids_track_1_based_source_indices(
        self, traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        # 20 structures, last 2 kept -> original indices 19, 20 (1-based)
        assert job._original_conformer_ids == ["19", "20"]


class TestGaussianTrajJobEnergySorting:
    def _job(self, traj_molecules, gaussian_settings):
        return GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
            proportion_structures_to_use=1.0,
        )

    def test_all_energies_matches_molecule_order(
        self, traj_molecules, gaussian_settings
    ):
        job = self._job(traj_molecules, gaussian_settings)
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
        self, traj_molecules, gaussian_settings
    ):
        job = GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        assert job.unique_structures == job.molecules
        assert job.num_unique_structures == job.num_structures

    def test_grouping_strategy_delegates_to_grouper(
        self, traj_molecules, gaussian_settings
    ):
        unique_subset = traj_molecules[-2:][:1]
        mock_grouper = MagicMock()
        mock_grouper.unique.return_value = unique_subset

        with patch(
            "chemsmart.utils.grouper.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            job = GaussianTrajJob(
                molecules=traj_molecules,
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
        self, traj_molecules, gaussian_settings, n_unique
    ):
        job = GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )
        fake_unique = job.molecules[:n_unique]
        job.grouper = MagicMock()
        job.grouper.unique.return_value = fake_unique
        return job

    def test_prepare_all_jobs_builds_indexed_labels(
        self, traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            traj_molecules, gaussian_settings, n_unique=2
        )
        jobs = job._prepare_all_jobs()
        assert [j.label for j in jobs] == ["traj_c1", "traj_c2"]

    def test_run_all_jobs_runs_everything_when_no_limit(
        self, traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            traj_molecules, gaussian_settings, n_unique=2
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
        self, traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            traj_molecules, gaussian_settings, n_unique=2
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
        self, traj_molecules, gaussian_settings
    ):
        job = self._job_with_fake_uniques(
            traj_molecules, gaussian_settings, n_unique=2
        )
        with patch.object(job, "_run_all_jobs") as mock_run_all:
            job._run()
        mock_run_all.assert_called_once()


class TestGaussianTrajJobCompletion:
    def _job(self, traj_molecules, gaussian_settings):
        return GaussianTrajJob(
            molecules=traj_molecules,
            settings=gaussian_settings,
            label="traj",
            jobrunner=None,
        )

    def test_is_complete_checks_all_jobs_when_no_limit(
        self, traj_molecules, gaussian_settings
    ):
        job = self._job(traj_molecules, gaussian_settings)
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
        self, traj_molecules, gaussian_settings
    ):
        job = self._job(traj_molecules, gaussian_settings)
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
        self, traj_molecules, gaussian_settings
    ):
        job = self._job(traj_molecules, gaussian_settings)
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
        self, traj_molecules, gaussian_settings
    ):
        job = self._job(traj_molecules, gaussian_settings)
        complete = MagicMock(is_complete=lambda: True)
        incomplete = MagicMock(is_complete=lambda: False)
        fake_jobs = [complete, incomplete]
        with patch.object(
            GaussianTrajJob,
            "all_structures_run_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.incomplete_structure_run_jobs == [incomplete]
