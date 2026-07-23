"""
Direct unit tests for :class:`GaussianLinkJob` in
``chemsmart.jobs.gaussian.link``.

Covers IRC-subjob construction (label rewriting, forward/reverse
selection), and the IRC-aware ``_run``/``_job_is_complete`` overrides.
Also documents a dead-code bug in ``backup_files`` (see
``BUGS_FOUND.md``).
"""

from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.link import GaussianLinkJob
from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings


@pytest.fixture()
def irc_settings():
    return GaussianLinkJobSettings(jobtype="irc", direction=None)


@pytest.fixture()
def opt_settings():
    return GaussianLinkJobSettings(jobtype="opt")


class TestIsIrcJob:
    def test_true_for_irc_jobtype(self, ethanol_molecule, irc_settings):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
            label="link_test",
        )
        assert job._is_irc_job() is True

    def test_false_for_non_irc_jobtype(self, ethanol_molecule, opt_settings):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=opt_settings,
            label="link_test",
        )
        assert job._is_irc_job() is False


class TestCreateIrcSubjob:
    def test_raises_for_non_irc_job(self, ethanol_molecule, opt_settings):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=opt_settings,
            label="link_test",
        )
        with pytest.raises(ValueError, match="only for IRC link jobs"):
            job._create_irc_subjob("f")

    def test_forward_subjob_has_ircf_jobtype(
        self, ethanol_molecule, irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
            label="link_test",
        )
        subjob = job._create_irc_subjob("f")
        assert subjob.settings.jobtype == "ircf"
        assert subjob.label == "link_test_f"

    def test_reverse_subjob_has_ircr_jobtype(
        self, ethanol_molecule, irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
            label="link_test",
        )
        subjob = job._create_irc_subjob("r")
        assert subjob.settings.jobtype == "ircr"
        assert subjob.label == "link_test_r"

    def test_irc_link_label_pattern_is_rewritten_in_place(
        self, ethanol_molecule, irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
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


class TestGetIrcJobs:
    def test_non_irc_job_returns_empty_list(
        self, ethanol_molecule, opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=opt_settings,
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

    def test_no_direction_returns_both(self, ethanol_molecule, irc_settings):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
            label="link_test",
        )
        jobs = job._get_irc_jobs()
        assert [j.settings.jobtype for j in jobs] == ["ircf", "ircr"]


class TestRunAndCompletion:
    def test_run_delegates_to_super_for_non_irc(
        self, ethanol_molecule, opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=opt_settings,
            label="link_test",
        )
        with patch(
            "chemsmart.jobs.gaussian.job.GaussianJob._run"
        ) as mock_super_run:
            job._run()
        mock_super_run.assert_called_once()

    def test_run_runs_each_irc_subjob(self, ethanol_molecule, irc_settings):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
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
        self, ethanol_molecule, opt_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=opt_settings,
            label="link_test",
        )
        with patch(
            "chemsmart.jobs.gaussian.job.GaussianJob._job_is_complete",
            return_value=True,
        ) as mock_super_complete:
            assert job._job_is_complete() is True
        mock_super_complete.assert_called_once()

    def test_job_is_complete_irc_requires_all_subjobs_complete(
        self, ethanol_molecule, irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
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
        self, ethanol_molecule, irc_settings
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
        self, ethanol_molecule, irc_settings
    ):
        job = GaussianLinkJob(
            molecule=ethanol_molecule,
            settings=irc_settings,
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
