"""
Direct unit tests for ``chemsmart.jobs.gaussian.uvvis.GaussianUVVISJob``.

The class is explicitly marked "incomplete and under development" in its
own docstring, and indeed its ``__init__`` always crashes (see
``BUGS_FOUND.md``) because it forwards an ``atoms=`` kwarg to
``GaussianJob.__init__``, which requires ``molecule=`` instead. Since no
real instance can be constructed through the public API, the workflow
helper methods below are tested by building a bare instance with
``object.__new__`` and only setting the attributes those methods touch,
mirroring the pattern used for other stub/incomplete classes elsewhere in
this test suite (e.g. ``test_thermochemistry_job_unit.py``).
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.uvvis import GaussianUVVISJob


def _make_job(settings=None):
    job = GaussianUVVISJob.__new__(GaussianUVVISJob)
    job.folder = "/tmp/uvvis_job"
    job.molecule = MagicMock(name="molecule")
    job.settings = settings or GaussianJobSettings.default()
    return job


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
