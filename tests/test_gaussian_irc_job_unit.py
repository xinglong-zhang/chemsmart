"""
Direct unit tests for :class:`GaussianIRCJob` in
``chemsmart.jobs.gaussian.irc``.

CLI-level propagation is already covered by
``TestGaussianCLIIrcCommand`` in ``test_gaussian_cli.py``, but that
mocks the job class entirely. These tests exercise the class directly:
sub-job label derivation, direction-based run/completion dispatch, and
backup file selection.
"""

from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.settings import GaussianIRCJobSettings


@pytest.fixture()
def irc_settings():
    return GaussianIRCJobSettings.default()


class TestGaussianIRCJobConstruction:
    def test_freq_disabled(self, a_molecule, irc_settings):
        irc_settings.freq = True
        job = GaussianIRCJob(molecule=a_molecule, settings=irc_settings)
        assert job.settings.freq is False

    def test_settings_class_returns_expected_type(self):
        assert GaussianIRCJob.settings_class() is GaussianIRCJobSettings


@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    from chemsmart.io.molecules.structure import Molecule

    return Molecule.from_filepath(single_molecule_xyz_file)


class TestSubJobLabelDerivation:
    def test_ircf_label_auto_generated_irc_suffix(
        self, a_molecule, irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf"

    def test_ircr_label_auto_generated_irc_suffix(
        self, a_molecule, irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircr_job().label == "mymol_ircr"

    def test_ircf_label_custom_label_gets_underscored_suffix(
        self, a_molecule, irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="custom_label",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "custom_label_ircf"

    def test_ircr_label_custom_label_gets_underscored_suffix(
        self, a_molecule, irc_settings
    ):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="custom_label",
            jobrunner=MagicMock(),
        )
        assert job._ircr_job().label == "custom_label_ircr"

    def test_ircf_label_flat_irc_appends_suffix(
        self, a_molecule, irc_settings
    ):
        irc_settings.flat_irc = True
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf_flat"

    def test_label_unchanged_when_direction_already_set(
        self, a_molecule, irc_settings
    ):
        # When direction is already embedded (e.g. by update_irc_label
        # upstream), the sub-job label is reused verbatim.
        irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_ircf",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().label == "mymol_ircf"

    def test_label_none_defaults_to_chemical_formula_based_label(
        self, a_molecule, irc_settings
    ):
        # GaussianJob.__init__ defaults a None label to the molecule's
        # chemical formula, so job.label is never actually None by the
        # time _ircf_job() runs; it just goes through the "custom
        # label" branch (no "_irc" suffix to detect) like any other
        # non-None label.
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label=None,
            jobrunner=MagicMock(),
        )
        expected_label = a_molecule.get_chemical_formula(empirical=True)
        assert job.label == expected_label
        assert job._ircf_job().label == f"{expected_label}_ircf"

    def test_ircf_jobtype_set_on_sub_settings(self, a_molecule, irc_settings):
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        assert job._ircf_job().settings.jobtype == "ircf"
        assert job._ircr_job().settings.jobtype == "ircr"


class TestRunDirectionDispatch:
    def test_forward_direction_runs_only_forward(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward = MagicMock()
        job._run_reverse = MagicMock()

        job._run()

        job._run_forward.assert_called_once()
        job._run_reverse.assert_not_called()

    def test_reverse_direction_runs_only_reverse(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward = MagicMock()
        job._run_reverse = MagicMock()

        job._run()

        job._run_forward.assert_not_called()
        job._run_reverse.assert_called_once()

    def test_no_direction_runs_both(self, a_molecule, irc_settings):
        irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
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
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=True)
        job._run_reverse_is_complete = MagicMock(return_value=False)

        assert job._job_is_complete() is True
        job._run_reverse_is_complete.assert_not_called()

    def test_reverse_direction_checks_only_reverse(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=False)
        job._run_reverse_is_complete = MagicMock(return_value=True)

        assert job._job_is_complete() is True
        job._run_forward_is_complete.assert_not_called()

    def test_no_direction_requires_both_complete(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job._run_forward_is_complete = MagicMock(return_value=True)
        job._run_reverse_is_complete = MagicMock(return_value=False)

        assert job._job_is_complete() is False


class TestBackupFiles:
    def test_forward_direction_backs_up_only_forward_files(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "forward"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=True)
        # inputfile, outputfile, chkfile -> 3 calls for forward only
        assert job.backup_file.call_count == 3

    def test_reverse_direction_backs_up_only_reverse_files(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = "reverse"
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=False)
        # inputfile, outputfile only (no chk) -> 2 calls
        assert job.backup_file.call_count == 2

    def test_no_direction_backs_up_both_directions(
        self, a_molecule, irc_settings
    ):
        irc_settings.direction = None
        job = GaussianIRCJob(
            molecule=a_molecule,
            settings=irc_settings,
            label="mymol_irc",
            jobrunner=MagicMock(),
        )
        job.backup_file = MagicMock()
        job.backup_files(backup_chk=True)
        # input+output for both directions (4) + chk for both (2) = 6
        assert job.backup_file.call_count == 6
