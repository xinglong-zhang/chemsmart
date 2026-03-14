"""
Tests for Gaussian CLI solvent options.

Verifies that the ``-sm``/``--solvent-model``, ``-si``/``--solvent-id``,
``-so``/``--solvent-options``, and ``--remove-solvent`` options on the
``gaussian`` CLI *group* are correctly propagated to every subcommand
(``opt``, ``td``, ``sp``, …) via the ``merge()`` mechanism.

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
    ):
        """``-sm smd -si water`` sets solvent on the opt job settings."""
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
                    "opt",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
            # Capture settings from the call
            if mock_job_cls.call_args is not None:
                captured["settings"] = mock_job_cls.call_args[1]["settings"]

        assert result.exit_code == 0, result.output
        assert "settings" in captured, "GaussianOptJob was never instantiated"
        settings = captured["settings"]
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"

    def test_solvent_options_iterative_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
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
