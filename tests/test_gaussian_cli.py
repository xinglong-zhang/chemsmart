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

from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian


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
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "iterative"

    def test_solvent_route_keyword_for_opt(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Route string for opt job contains ``scrf=(smd,solvent=water,iterative)``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert "settings" in captured
        route_string = captured["settings"].route_string
        assert "scrf=(smd,solvent=water,iterative)" in route_string

    def test_remove_solvent_clears_solvent_from_opt(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        # The ``solv`` project sets solvent_model=smd and solvent_id=toluene
        # for every job type (including opt).  ``--remove-solvent`` must strip
        # these from the merged settings.
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestGaussianSolventCLITdCommand:
    """CLI solvent options propagated to the ``td`` subcommand."""

    def test_solvent_model_and_id_override_td_project_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``-sm smd -si water`` overrides project td solvent (toluene→water)."""
        # ``solv`` project has smd/toluene for td; CLI overrides solvent_id.
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_td_route_keyword_with_smd_water_iterative(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """TD route string contains ``scrf=(smd,solvent=water,iterative)``."""
        # ``solv`` project has smd/toluene for td; CLI overrides to water+iterative.
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert "settings" in captured
        route_string = captured["settings"].route_string
        assert "scrf=(smd,solvent=water,iterative)" in route_string

    def test_remove_solvent_clears_solvent_from_td(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``--remove-solvent`` nulls solvent settings for a td job."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_no_solvent_options_leaves_project_settings_unchanged(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Without solvent CLI options the project solvent settings are kept."""
        # ``solv`` project has smd/toluene for td; no CLI override → kept.
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"


class TestGaussianCLISinglePointCommand:
    """CLI tests for the ``sp`` (single point) subcommand."""

    def test_basic_sp_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``sp`` subcommand creates a ``GaussianSinglePointJob``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert (
            "settings" in captured
        ), "GaussianSinglePointJob was never instantiated"
        # gas_solv project sp settings use ``solv`` config: def2tzvp + smd/toluene
        settings = captured["settings"]
        assert settings.basis == "def2tzvp"
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_sp_subcommand_solvent_override(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """sp-level ``-sm``/``-si`` options override project solvent settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "pcm"
        assert settings.solvent_id == "acetonitrile"

    def test_sp_subcommand_remove_solvent(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """sp-level ``--remove-solvent`` strips solvent from project settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_sp_group_level_solvent_propagated(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level ``-sm``/``-si`` options are merged into sp settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "thf"


class TestGaussianCLITsCommand:
    """CLI tests for the ``ts`` (transition state) subcommand."""

    def test_basic_ts_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``ts`` subcommand creates a ``GaussianTSJob`` with gas settings."""
        runner = CliRunner()
        captured = {}

        with patch("chemsmart.jobs.gaussian.ts.GaussianTSJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured, "GaussianTSJob was never instantiated"
        # gas_solv project ts settings use ``gas`` config: def2svp, no solvent
        settings = captured["settings"]
        assert settings.basis == "def2svp"
        assert settings.solvent_model is None

    def test_ts_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``ts`` with ``solv`` project inherits solvent settings."""
        runner = CliRunner()
        captured = {}

        with patch("chemsmart.jobs.gaussian.ts.GaussianTSJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_ts_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level solvent options are propagated to ``ts`` settings."""
        runner = CliRunner()
        captured = {}

        with patch("chemsmart.jobs.gaussian.ts.GaussianTSJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "dmso"


class TestGaussianCLIIrcCommand:
    """CLI tests for the ``irc`` (Intrinsic Reaction Coordinate) subcommand."""

    def test_basic_irc_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``irc`` subcommand creates a ``GaussianIRCJob``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured, "GaussianIRCJob was never instantiated"
        settings = captured["settings"]
        assert settings.basis == "def2svp"

    def test_irc_direction_forward_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``-d forward`` sets the IRC direction to ``forward``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        assert captured["settings"].direction == "forward"

    def test_irc_direction_reverse_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``-d reverse`` sets the IRC direction to ``reverse``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        assert captured["settings"].direction == "reverse"

    def test_irc_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level solvent options are propagated to ``irc`` settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.irc.GaussianIRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "methanol"


class TestGaussianCLIScanCommand:
    """CLI tests for the ``scan`` (potential energy surface scan) subcommand."""

    def test_basic_scan_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``scan`` subcommand creates a ``GaussianScanJob``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured, "GaussianScanJob was never instantiated"
        settings = captured["settings"]
        assert settings.basis == "def2svp"

    def test_scan_settings_from_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``scan`` with ``solv`` project inherits solvent settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_scan_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level solvent options are propagated to ``scan`` settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestGaussianCLICrestCommand:
    """CLI tests for the ``crest`` (conformer search) subcommand."""

    def test_basic_crest_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``crest -j opt`` subcommand creates a ``GaussianCrestJob``."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert (
            "settings" in captured
        ), "GaussianCrestJob was never instantiated"
        settings = captured["settings"]
        assert settings.basis == "def2svp"

    def test_crest_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``crest`` with ``solv`` project inherits solvent settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_crest_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level solvent options are propagated to ``crest`` settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.crest.GaussianCrestJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestGaussianCLIQrcCommand:
    """CLI tests for the ``qrc`` (Quadratic Reaction Coordinate) subcommand."""

    def test_basic_qrc_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``qrc`` subcommand creates a ``GaussianQRCJob`` (default jobtype=opt)."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured, "GaussianQRCJob was never instantiated"
        settings = captured["settings"]
        assert settings.basis == "def2svp"

    def test_qrc_settings_from_solv_project(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``qrc`` with ``solv`` project inherits solvent settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "toluene"

    def test_qrc_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """Group-level solvent options are propagated to ``qrc`` settings."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_qrc_explicit_ts_jobtype(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``qrc -j ts`` uses TS settings from the project for the QRC job."""
        runner = CliRunner()
        captured = {}

        with patch(
            "chemsmart.jobs.gaussian.qrc.GaussianQRCJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured
        # TS settings from gas_solv use the ``gas`` config: def2svp
        settings = captured["settings"]
        assert settings.basis == "def2svp"
