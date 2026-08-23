"""
Direct unit tests for :class:`ORCAQRCJob` (``chemsmart.jobs.orca.qrc``).

Mirrors the pattern used for ``GaussianQRCJob`` in
``TestGaussianQRCJobs`` (test_GaussianJobs.py): a mocked ``Molecule``
with ``has_vibrations``/``vibrationally_displaced`` and a real
``ORCAJobSettings`` instance to pass the parent's ``isinstance`` checks.
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca.qrc import ORCAQRCJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.runner import JobRunner


@pytest.fixture()
def mock_molecule():
    mol = MagicMock(spec=Molecule)
    mol.has_vibrations = True
    mol.copy.return_value = mol
    mol.vibrationally_displaced.return_value = mol
    mol.get_chemical_formula.return_value = "C1H4"
    return mol


@pytest.fixture()
def real_settings():
    settings = ORCAJobSettings.default()
    settings.jobtype = "qrc"
    return settings


@pytest.fixture()
def mock_jobrunner():
    return MagicMock(spec=JobRunner)


class TestORCAQRCJobConstruction:
    def test_init_raises_if_no_vibrations(self, mock_molecule, real_settings):
        mock_molecule.has_vibrations = False
        with pytest.raises(ValueError, match="no vibrational modes"):
            ORCAQRCJob(molecule=mock_molecule, settings=real_settings)

    def test_init_stores_qrc_parameters(self, mock_molecule, real_settings):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            label="test_qrc",
            mode_idx=2,
            amp=1.2,
            nframes=10,
            normalize=True,
            return_xyz=True,
        )
        assert job.mode_idx == 2
        assert job.amp == 1.2
        assert job.nframes == 10
        assert job.normalize is True
        assert job.return_xyz is True
        assert job.jobtype == "qrc"


class TestORCAQRCJobDisplacedMolecules:
    def test_qrcf_molecule_uses_positive_amplitude(
        self, mock_molecule, real_settings
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule, settings=real_settings, amp=0.5
        )
        _ = job.qrcf_molecule
        mock_molecule.vibrationally_displaced.assert_called_with(
            mode_idx=job.mode_idx,
            amp=0.5,
            nframes=job.nframes,
            phase=job.phase,
            normalize=job.normalize,
            return_xyz=job.return_xyz,
        )

    def test_qrcr_molecule_uses_negative_amplitude(
        self, mock_molecule, real_settings
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule, settings=real_settings, amp=0.5
        )
        _ = job.qrcr_molecule
        mock_molecule.vibrationally_displaced.assert_called_with(
            mode_idx=job.mode_idx,
            amp=-0.5,
            nframes=job.nframes,
            phase=job.phase,
            normalize=job.normalize,
            return_xyz=job.return_xyz,
        )


class TestORCAQRCJobBothJobs:
    def test_prepare_both_qrc_jobs_creates_forward_and_reverse(
        self, mock_molecule, real_settings, mock_jobrunner
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )

        with patch("chemsmart.jobs.orca.qrc.ORCAGeneralJob") as mock_general:
            forward_job = MagicMock()
            reverse_job = MagicMock()
            mock_general.side_effect = [forward_job, reverse_job]

            jobs = job._prepare_both_qrc_jobs()

        assert jobs == [forward_job, reverse_job]
        assert mock_general.call_count == 2
        labels = [c.kwargs["label"] for c in mock_general.call_args_list]
        assert labels == ["test_qrcf_qrc", "test_qrcr_qrc"]

    def test_run_both_jobs_runs_forward_and_reverse(
        self, mock_molecule, real_settings, mock_jobrunner
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            jobrunner=mock_jobrunner,
            label="test_qrc",
        )

        with patch("chemsmart.jobs.orca.qrc.ORCAGeneralJob") as mock_general:
            forward_job = MagicMock()
            reverse_job = MagicMock()
            mock_general.side_effect = [forward_job, reverse_job]

            job._run_both_jobs()

        forward_job.run.assert_called_once()
        reverse_job.run.assert_called_once()

    def test_run_delegates_to_run_both_jobs(
        self, mock_molecule, real_settings, mock_jobrunner
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            jobrunner=mock_jobrunner,
        )
        job._run_both_jobs = MagicMock()
        job._run()
        job._run_both_jobs.assert_called_once()

    def test_is_complete_true_when_both_jobs_complete(
        self, mock_molecule, real_settings, mock_jobrunner
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            jobrunner=mock_jobrunner,
        )
        job1 = MagicMock()
        job1.is_complete.return_value = True
        job2 = MagicMock()
        job2.is_complete.return_value = True
        job._prepare_both_qrc_jobs = MagicMock(return_value=[job1, job2])

        assert job.is_complete() is True

    def test_is_complete_false_when_one_job_incomplete(
        self, mock_molecule, real_settings, mock_jobrunner
    ):
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=real_settings,
            jobrunner=mock_jobrunner,
        )
        job1 = MagicMock()
        job1.is_complete.return_value = True
        job2 = MagicMock()
        job2.is_complete.return_value = False
        job._prepare_both_qrc_jobs = MagicMock(return_value=[job1, job2])

        assert job.is_complete() is False
