"""
Direct unit tests for :class:`IterateJob` and :class:`IterateJobSettings`
in ``chemsmart.jobs.iterate``.
"""

from unittest.mock import MagicMock, patch

from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.settings import IterateJobSettings


class TestIterateJobSettings:
    def test_defaults(self):
        settings = IterateJobSettings()
        assert settings.config_file is None
        assert settings.skeleton_list == []
        assert settings.substituent_list == []
        assert settings.method == "lagrange_multipliers"
        assert settings.sphere_direction_samples_num == 96
        assert settings.axial_rotations_sample_num == 6

    def test_custom_values(self):
        settings = IterateJobSettings(
            config_file="config.toml",
            method="custom_method",
            sphere_direction_samples_num=50,
            axial_rotations_sample_num=12,
        )
        assert settings.config_file == "config.toml"
        assert settings.method == "custom_method"
        assert settings.sphere_direction_samples_num == 50
        assert settings.axial_rotations_sample_num == 12

    def test_copy_produces_independent_deep_copy(self):
        settings = IterateJobSettings()
        settings.skeleton_list.append({"file": "a.xyz"})
        copied = settings.copy()

        assert copied is not settings
        assert copied.skeleton_list == settings.skeleton_list
        copied.skeleton_list.append({"file": "b.xyz"})
        assert copied.skeleton_list != settings.skeleton_list
        assert len(settings.skeleton_list) == 1


class TestIterateJobConstruction:
    def test_defaults(self):
        job = IterateJob()
        assert job.label == "iterate_job"
        assert job.molecule is None
        assert job.skip_completed is False
        assert isinstance(job.settings, IterateJobSettings)
        assert job.nprocs == 1
        assert job.timeout == 120
        assert job.separate_outputs is False
        assert job.output_directory is None

    def test_custom_settings_used_directly(self):
        settings = IterateJobSettings(method="custom")
        job = IterateJob(settings=settings)
        assert job.settings is settings

    def test_settings_class_returns_expected_type(self):
        assert IterateJob.settings_class() is IterateJobSettings


class TestIterateJobOutputfile:
    def test_adds_xyz_extension_when_missing(self):
        job = IterateJob(outputfile="my_results")
        assert job.outputfile == "my_results.xyz"

    def test_preserves_existing_xyz_extension(self):
        job = IterateJob(outputfile="my_results.xyz")
        assert job.outputfile == "my_results.xyz"


class TestIterateJobRun:
    def test_run_uses_existing_jobrunner(self):
        mock_runner = MagicMock()
        job = IterateJob(jobrunner=mock_runner)
        result = job.run()

        mock_runner.run.assert_called_once_with(job)
        assert result == job.outputfile

    def test_run_creates_default_runner_when_none_given(self):
        mock_runner_instance = MagicMock()
        with patch(
            "chemsmart.jobs.iterate.runner.IterateJobRunner",
            return_value=mock_runner_instance,
        ) as mock_runner_cls:
            job = IterateJob(jobrunner=None)
            result = job.run()

        mock_runner_cls.assert_called_once_with()
        mock_runner_instance.run.assert_called_once_with(job)
        assert result == job.outputfile
