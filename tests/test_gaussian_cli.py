"""
Tests for Gaussian CLI option propagation and subcommand behaviour.

This module verifies that solvent-related options
 (``-sm``/``--solvent-model``, ``-si``/``--solvent-id``,
 ``-so``/``--solvent-options``, and ``--remove-solvent``) on the
 ``gaussian`` CLI *group* are correctly propagated to every relevant
 subcommand (``opt``, ``td``, ``sp``, …) via the ``merge()`` mechanism.
 It also exercises non-solvent Gaussian CLI functionality for various
 subcommands (such as ``sp``, ``ts``, ``irc``, ``scan``, and ``crest``),
 including job type flags, directions, scan/QRC settings, and related
 options.

Each test uses :class:`click.testing.CliRunner` to invoke the ``gaussian``
group and :mod:`unittest.mock` to intercept the job constructor so that
the merged :class:`~chemsmart.jobs.gaussian.settings.GaussianJobSettings`
can be inspected without running an actual calculation.
"""

from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.main import entry_point


class TestGaussianCLIPubChemOptCommand:
    """CLI tests for ``--pubchem`` structure loading without ``-f``."""

    def test_pubchem_only_opt_passes_molecule_to_job(
        self,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--pubchem <cid> opt`` loads a structure and passes it to the job."""
        pubchem_molecule = MagicMock(name="pubchem_molecule")

        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[pubchem_molecule],
        ) as mock_from_pubchem:
            result, settings = run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                [
                    "-p",
                    "gas_solv",
                    "--pubchem",
                    "222",
                    "-l",
                    "ammonia",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )

        assert result.exit_code == 0, result.output
        mock_from_pubchem.assert_called_once_with(
            identifier="222", return_list=True
        )
        assert settings is not None, "GaussianOptJob was never instantiated"


class TestGaussianSolventCLIOptCommand:
    """CLI solvent options propagated to the ``opt`` subcommand."""

    def test_solvent_model_and_id_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-sm smd -si water`` sets solvent on the opt job settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert settings is not None, "GaussianOptJob was never instantiated"
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_solvent_options_iterative_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-sm smd -si water -so iterative`` sets iterative solvent on opt."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "iterative",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "iterative"

    def test_solvent_route_keyword_for_opt(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Route string for opt job contains ``scrf=(smd,solvent=water,iterative)``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "iterative",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert "scrf=(smd,solvent=water,iterative)" in settings.route_string

    def test_remove_solvent_clears_solvent_from_opt(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        # The ``solv`` project sets solvent_model=smd and solvent_id=toluene
        # for every job type (including opt).  ``--remove-solvent`` must strip
        # these from the merged settings.
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianSolventCLITdCommand:
    """CLI solvent options propagated to the ``td`` subcommand."""

    def test_solvent_model_and_id_override_td_project_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-sm smd -si water`` overrides project td solvent (toluene→water)."""
        # ``solv`` project has smd/toluene for td; CLI overrides solvent_id.
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "td",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_td_route_keyword_with_smd_water_iterative(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """TD route string contains ``scrf=(smd,solvent=water,iterative)``."""
        # ``solv`` project has smd/toluene for td; CLI overrides to water+iterative.
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "iterative",
                "td",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert "scrf=(smd,solvent=water,iterative)" in settings.route_string

    def test_remove_solvent_clears_solvent_from_td(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--remove-solvent`` nulls solvent settings for a td job."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "td",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_no_solvent_options_leaves_project_settings_unchanged(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Without solvent CLI options the project solvent settings are kept."""
        # ``solv`` project has smd/toluene for td; no CLI override → kept.
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "td",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"


class TestGaussianCLISinglePointCommand:
    """CLI tests for the ``sp`` (single point) subcommand."""

    def test_basic_sp_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``sp`` subcommand creates a ``GaussianSinglePointJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.basis == "def2tzvp"
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_sp_subcommand_solvent_override(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """sp-level ``-sm``/``-si`` options override project solvent settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "sp",
                "-sm",
                "pcm",
                "-si",
                "acetonitrile",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "pcm"
        assert settings.solvent_id == "acetonitrile"

    def test_sp_subcommand_remove_solvent(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """sp-level ``--remove-solvent`` strips solvent from project settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "sp",
                "--remove-solvent",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_sp_group_level_solvent_propagated(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level ``-sm``/``-si`` options are merged into sp settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "thf",
                "sp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "thf"


class TestGaussianCLITsCommand:
    """CLI tests for the ``ts`` (transition state) subcommand."""

    def test_basic_ts_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``ts`` subcommand creates a ``GaussianTSJob`` with gas settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.ts.GaussianTSJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "ts",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        # gas_solv project ts settings use ``gas`` config: def2svp, no solvent
        assert settings.basis == "def2svp"
        assert settings.solvent_model is None

    def test_ts_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``ts`` with ``solv`` project inherits solvent settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.ts.GaussianTSJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "ts",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_ts_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``ts`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.ts.GaussianTSJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "dmso",
                "ts",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "dmso"


class TestGaussianCLIIrcCommand:
    """CLI tests for the ``irc`` (Intrinsic Reaction Coordinate) subcommand."""

    def test_basic_irc_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``irc`` subcommand creates a ``GaussianIRCJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "irc",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert settings.basis == "def2svp"

    def test_irc_direction_forward_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-d forward`` sets the IRC direction to ``forward``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "irc",
                "-d",
                "forward",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.direction == "forward"

    def test_irc_direction_reverse_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-d reverse`` sets the IRC direction to ``reverse``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "irc",
                "-d",
                "reverse",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.direction == "reverse"

    def test_irc_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``irc`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "methanol",
                "irc",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "methanol"


class TestGaussianCLIScanCommand:
    """CLI tests for the ``scan`` (potential energy surface scan) subcommand."""

    def test_basic_scan_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``scan`` subcommand creates a ``GaussianScanJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "scan",
                "-c",
                "[[1,2]]",
                "-s",
                "0.1",
                "-n",
                "10",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.basis == "def2svp"


class TestGaussianBatchTriggeringGate:
    """Regression tests for Gaussian multi-molecule fan-out gating."""

    @pytest.mark.parametrize(
        "subcommand,job_class_path,extra_args",
        [
            ("opt", "chemsmart.jobs.gaussian.opt.GaussianOptJob", []),
            (
                "sp",
                "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
                [],
            ),
        ],
    )
    def test_fanout_when_multiple_indices(
        self,
        subcommand,
        job_class_path,
        extra_args,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Multiple selected molecule indices fan out into one job each."""
        runner = CliRunner()
        with patch(job_class_path) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    subcommand,
                    *extra_args,
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2

    @pytest.mark.parametrize(
        "subcommand,job_class_path,extra_args",
        [
            ("opt", "chemsmart.jobs.gaussian.opt.GaussianOptJob", []),
            (
                "sp",
                "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob",
                [],
            ),
        ],
    )
    def test_fanout_still_happens_with_no_run_in_parallel(
        self,
        subcommand,
        job_class_path,
        extra_args,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """jobrunner.no_run_in_parallel controls execution mode, not fan-out."""
        gaussian_jobrunner_no_scratch.no_run_in_parallel = True

        runner = CliRunner()
        with patch(job_class_path) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    subcommand,
                    *extra_args,
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        for call in mock_job_cls.call_args_list:
            assert call.kwargs["jobrunner"].no_run_in_parallel is True


class TestGaussianRunSubNoParallelIntegration:
    """Integration-style tests for `run/sub --no-run-in-parallel`."""

    def test_run_no_run_in_parallel_builds_serial_batch_job(
        self,
        multiple_molecules_xyz_file,
        pbs_server,
    ):
        """`run --no-run-in-parallel` builds Gaussian batch with serial execution."""
        from chemsmart.jobs.job import Job

        class DummyBatchJob(Job):
            TYPE = "g16opt"

            def run(self, **kwargs):
                return None

        runner = CliRunner()
        dummy_batch_job = DummyBatchJob(
            molecule=None,
            label="opt_batch",
            jobrunner=None,
        )
        with (
            patch(
                "chemsmart.cli.run.Server.from_servername",
                return_value=pbs_server,
            ),
            patch(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob"
            ) as mock_job_cls,
            patch(
                "chemsmart.jobs.gaussian.batch.GaussianBatchJob",
                return_value=dummy_batch_job,
            ) as mock_batch_cls,
            patch.object(dummy_batch_job, "run") as mock_batch_run,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "run",
                    "-s",
                    "PBS",
                    "-N",
                    "1",
                    "--no-run-in-parallel",
                    "--no-scratch",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        assert mock_batch_cls.call_count == 1
        assert mock_batch_cls.call_args.kwargs["no_run_in_parallel"] is True
        assert len(mock_batch_cls.call_args.kwargs["jobs"]) == 2
        assert mock_batch_run.call_count == 1

    def test_sub_no_run_in_parallel_submits_array_batch_job(
        self,
        multiple_molecules_xyz_file,
        pbs_server,
    ):
        """`sub --no-run-in-parallel` submits BatchJob as array with %1."""
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob

        runner = CliRunner()
        mock_batch_job = GaussianBatchJob(
            jobs=[MagicMock(label="j1"), MagicMock(label="j2")],
            label="mols_batch",
        )
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob"
            ) as mock_job_cls,
            patch(
                "chemsmart.jobs.gaussian.batch.GaussianBatchJob",
                return_value=mock_batch_job,
            ) as mock_batch_cls,
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "-N",
                    "1",
                    "--no-run-in-parallel",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        assert mock_batch_cls.call_count == 1
        assert mock_batch_cls.call_args.kwargs["no_run_in_parallel"] is True
        assert mock_submit_array.call_count == 1
        assert mock_submit_array.call_args.kwargs["test"] is True
        assert mock_submit_array.call_args.kwargs["num_nodes"] == 1
        assert (
            mock_submit_array.call_args.kwargs["batch_label"] == "mols_batch"
        )
        assert len(mock_submit_array.call_args.kwargs["jobs"]) == 2
        assert (
            "--no-run-in-parallel"
            in mock_submit_array.call_args.kwargs["cli_args"]
        )

    def test_sub_default_is_serial_array_throttle(
        self,
        multiple_molecules_xyz_file,
        pbs_server,
    ):
        """Default ``sub`` (no parallel flag) uses %1 array throttle."""
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob

        runner = CliRunner()
        mock_batch_job = GaussianBatchJob(
            jobs=[MagicMock(label="j1"), MagicMock(label="j2")],
            label="mols_batch",
        )
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch("chemsmart.jobs.gaussian.opt.GaussianOptJob"),
            patch(
                "chemsmart.jobs.gaussian.batch.GaussianBatchJob",
                return_value=mock_batch_job,
            ) as mock_batch_cls,
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "-N",
                    "4",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_batch_cls.call_args.kwargs["no_run_in_parallel"] is True
        assert mock_submit_array.call_count == 1
        # -N 4 ignored for throttle when serial default applies
        assert mock_submit_array.call_args.kwargs["num_nodes"] == 1

    def test_sub_run_in_parallel_uses_num_nodes_as_array_throttle(
        self,
        multiple_molecules_xyz_file,
        pbs_server,
    ):
        """`sub --run-in-parallel -N 2` throttles array concurrency to 2."""
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob

        runner = CliRunner()
        mock_batch_job = GaussianBatchJob(
            jobs=[
                MagicMock(label="j1"),
                MagicMock(label="j2"),
                MagicMock(label="j3"),
            ],
            label="mols_batch",
        )
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch("chemsmart.jobs.gaussian.opt.GaussianOptJob"),
            patch(
                "chemsmart.jobs.gaussian.batch.GaussianBatchJob",
                return_value=mock_batch_job,
            ),
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "-N",
                    "2",
                    "--run-in-parallel",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1,2,3",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock_submit_array.call_count == 1
        assert mock_submit_array.call_args.kwargs["num_nodes"] == 2
        assert len(mock_submit_array.call_args.kwargs["jobs"]) == 3

    def test_sub_run_in_parallel_expands_qrc_to_array(
        self,
        single_molecule_xyz_file,
        pbs_server,
    ):
        """`sub --run-in-parallel ... qrc` expands forward/reverse to an array."""
        child_f = MagicMock(label="ts_qrcf", PROGRAM="Gaussian")
        child_r = MagicMock(label="ts_qrcr", PROGRAM="Gaussian")
        mock_qrc = MagicMock()
        mock_qrc.PROGRAM = "Gaussian"
        mock_qrc.label = "ts_qrc"
        mock_qrc.get_array_child_jobs.return_value = [child_f, child_r]

        runner = CliRunner()
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch(
                "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
                return_value=mock_qrc,
            ),
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
            patch.object(pbs_server, "submit") as mock_submit,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "--run-in-parallel",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "qrc",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_submit.assert_not_called()
        assert mock_submit_array.call_count == 1
        assert len(mock_submit_array.call_args.kwargs["jobs"]) == 2
        assert (
            mock_submit_array.call_args.kwargs["batch_label"] == "ts_qrc_array"
        )
        cli_args = mock_submit_array.call_args.kwargs["cli_args"]
        assert "qrc" in cli_args

    def test_sub_no_run_in_parallel_submits_qrc_as_single_parent(
        self,
        single_molecule_xyz_file,
        pbs_server,
    ):
        """`sub --no-run-in-parallel ... qrc` submits one nestable parent job."""
        child_f = MagicMock(label="ts_qrcf", PROGRAM="Gaussian")
        child_r = MagicMock(label="ts_qrcr", PROGRAM="Gaussian")
        mock_qrc = MagicMock()
        mock_qrc.PROGRAM = "Gaussian"
        mock_qrc.label = "ts_qrc"
        mock_qrc.get_array_child_jobs.return_value = [child_f, child_r]

        runner = CliRunner()
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch(
                "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
                return_value=mock_qrc,
            ),
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
            patch.object(pbs_server, "submit") as mock_submit,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "--no-run-in-parallel",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "qrc",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_submit_array.assert_not_called()
        assert mock_submit.call_count == 1
        assert mock_submit.call_args.kwargs["job"] is mock_qrc
        assert "qrc" in mock_submit.call_args.kwargs["cli_args"]

    def test_sub_default_submits_qrc_as_single_parent(
        self,
        single_molecule_xyz_file,
        pbs_server,
    ):
        """Default ``sub ... qrc`` (no parallel flag) is one nestable parent."""
        child_f = MagicMock(label="ts_qrcf", PROGRAM="Gaussian")
        child_r = MagicMock(label="ts_qrcr", PROGRAM="Gaussian")
        mock_qrc = MagicMock()
        mock_qrc.PROGRAM = "Gaussian"
        mock_qrc.label = "ts_qrc"
        mock_qrc.get_array_child_jobs.return_value = [child_f, child_r]

        runner = CliRunner()
        with (
            patch(
                "chemsmart.cli.sub.Server.from_servername",
                return_value=pbs_server,
            ),
            patch(
                "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
                return_value=mock_qrc,
            ),
            patch(
                "chemsmart.settings.server.Server.submit_array_job"
            ) as mock_submit_array,
            patch.object(pbs_server, "submit") as mock_submit,
        ):
            result = runner.invoke(
                entry_point,
                [
                    "sub",
                    "-s",
                    "PBS",
                    "--test",
                    "gaussian",
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "qrc",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_submit_array.assert_not_called()
        assert mock_submit.call_count == 1
        assert mock_submit.call_args.kwargs["job"] is mock_qrc

    def test_scan_settings_from_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``scan`` with ``solv`` project inherits solvent settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "scan",
                "-c",
                "[[1,2]]",
                "-s",
                "0.1",
                "-n",
                "10",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_scan_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``scan`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "scan",
                "-c",
                "[[1,2]]",
                "-s",
                "0.1",
                "-n",
                "10",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_scan_multiple_coords_single_step_size_and_num_steps(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Multiple scan coordinates with a single step_size and num_steps broadcasts them."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "scan",
                "-c",
                "[[1,2],[2,3]]",
                "-s",
                "-0.1",
                "-n",
                "10",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.modred["step_size"] == [-0.1, -0.1]
        assert settings.modred["num_steps"] == [10, 10]

    def test_scan_multiple_coords_explicit_step_size_and_num_steps(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Multiple scan coordinates with explicit per-coordinate step_size and num_steps."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "scan",
                "-c",
                "[[1,2],[2,3]]",
                "-s",
                "[-0.1,-0.2]",
                "-n",
                "[10,15]",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.modred["step_size"] == [-0.1, -0.2]
        assert settings.modred["num_steps"] == [10, 15]


class TestGaussianCLICrestCommand:
    """CLI tests for the ``crest`` (conformer search) subcommand."""

    def test_basic_crest_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``crest -j opt`` subcommand creates a ``GaussianCrestJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "crest",
                "-j",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.basis == "def2svp"

    def test_crest_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``crest`` with ``solv`` project inherits solvent settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "crest",
                "-j",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_crest_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``crest`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "crest",
                "-j",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestGaussianCLIQrcCommand:
    """CLI tests for the ``qrc`` (Quadratic Reaction Coordinate) subcommand."""

    def test_basic_qrc_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``qrc`` subcommand creates a ``GaussianQRCJob`` (default jobtype=opt)."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "qrc",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.basis == "def2svp"

    def test_qrc_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``qrc`` with ``solv`` project inherits solvent settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "qrc",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_qrc_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``qrc`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "qrc",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_qrc_explicit_ts_jobtype(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``qrc -j ts`` uses TS settings from the project for the QRC job."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "qrc",
                "-j",
                "ts",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.basis == "def2svp"
