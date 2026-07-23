"""
Direct unit tests for four small, thin Gaussian job wrapper classes that
had little to no direct coverage: ``GaussianWBIJob``, ``GaussianNCIJob``,
``GaussianRESPJob``, and ``GaussianTDDFTJob``.

Also documents a bug in ``GaussianNCIJob.__init__`` (see BUGS_FOUND.md):
it forwards the ``JobRunner`` *class* to the parent constructor instead
of the caller-supplied ``jobrunner`` argument.
"""

import pytest

from chemsmart.jobs.gaussian.nci import GaussianNCIJob
from chemsmart.jobs.gaussian.resp import GaussianRESPJob
from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianTDDFTJobSettings,
)
from chemsmart.jobs.gaussian.tddft import GaussianTDDFTJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob
from chemsmart.jobs.runner import JobRunner


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


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
