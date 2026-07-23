"""
Direct unit tests for :class:`GaussianCrestJob` in
``chemsmart.jobs.gaussian.crest``.

Exercises construction/validation, the conformer-job preparation and
completion-tracking properties, and the run/is_complete orchestration
directly against the class rather than through the CLI.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.crest import GaussianCrestJob
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


@pytest.fixture()
def conformers(ethanol_molecule):
    return [ethanol_molecule.copy() for _ in range(4)]


class TestGaussianCrestJobConstruction:
    def test_defaults_num_confs_to_opt_to_all_conformers(
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        assert job.num_conformers == 4
        assert job.num_confs_to_opt == 4
        assert job.all_conformers == conformers
        assert job.grouper is None

    def test_num_confs_to_run_is_respected(
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
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
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        assert (
            job.molecule.get_chemical_formula()
            == conformers[0].get_chemical_formula()
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


class TestGaussianCrestJobProperties:
    def _job(self, conformers, gaussian_settings, **kwargs):
        return GaussianCrestJob(
            molecules=conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
            **kwargs,
        )

    def test_all_conformers_jobs_builds_one_job_per_conformer_with_indexed_labels(
        self, conformers, gaussian_settings
    ):
        job = self._job(conformers, gaussian_settings)
        sub_jobs = job.all_conformers_jobs
        assert len(sub_jobs) == 4
        assert [j.label for j in sub_jobs] == [
            "crest_test_c1",
            "crest_test_c2",
            "crest_test_c3",
            "crest_test_c4",
        ]

    def test_incomplete_conformers_jobs_filters_by_is_complete(
        self, conformers, gaussian_settings
    ):
        job = self._job(conformers, gaussian_settings)
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
        self, conformers, gaussian_settings
    ):
        job = self._job(conformers, gaussian_settings)
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
        self, conformers, gaussian_settings
    ):
        job = self._job(conformers, gaussian_settings)
        fake_jobs = [MagicMock(is_complete=lambda: True) for _ in range(4)]

        with patch.object(
            GaussianCrestJob,
            "all_conformers_jobs",
            new_callable=lambda: property(lambda self: fake_jobs),
        ):
            assert job.last_run_job_index == 4


class TestGaussianCrestJobRunAndCompletion:
    def test_run_all_jobs_only_runs_up_to_num_confs_to_opt(
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
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
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
            settings=gaussian_settings,
            label="crest_test",
            jobrunner=None,
        )
        with patch.object(job, "_run_all_jobs") as mock_run_all:
            job._run()
        mock_run_all.assert_called_once()

    def test_is_complete_true_only_when_all_within_limit_are_complete(
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
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
        self, conformers, gaussian_settings
    ):
        job = GaussianCrestJob(
            molecules=conformers,
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
        self, conformers, gaussian_settings
    ):
        unique_subset = conformers[:2]
        mock_grouper = MagicMock()
        mock_grouper.unique.return_value = unique_subset

        with patch(
            "chemsmart.utils.grouper.StructureGrouperFactory.create",
            return_value=mock_grouper,
        ) as mock_create:
            job = GaussianCrestJob(
                molecules=conformers,
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
