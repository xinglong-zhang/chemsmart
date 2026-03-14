"""
Tests for ORCA CLI option propagation and subcommand behaviour.

This module verifies that solvent-related options
 (``-sm``/``--solvent-model``, ``-si``/``--solvent-id``,
 ``-so``/``--solvent-options``, and ``--remove-solvent``) on the
 ``orca`` CLI *group* and on individual subcommands (``sp``, ``opt``,
 ``ts``) are correctly propagated to the merged
 :class:`~chemsmart.jobs.orca.settings.ORCAJobSettings`.

Each test uses :class:`click.testing.CliRunner` to invoke the ``orca``
group and :mod:`unittest.mock` to intercept the job constructor so that
the merged settings can be inspected without running an actual calculation.
"""


class TestORCASolventCLISpCommand:
    """CLI solvent options propagated to the ``sp`` subcommand."""

    def test_solvent_model_and_id_injected_into_sp_settings_group_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` at group level sets solvent on sp settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "water",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCASinglePointJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_model_and_id_injected_into_sp_settings_subcommand_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``sp -sm cpcm -si water`` at subcommand level sets solvent on sp settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "cpcm",
                "-si",
                "water",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCASinglePointJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_options_injected_into_sp_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water -so 'Epsilon 78.36'`` sets additional options on sp."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "water",
                "-so",
                "Epsilon 78.36",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_remove_solvent_clears_solvent_from_sp(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        # The ``solv`` project sets solvent_model=smd and solvent_id=cyclohexane
        # for every job type.  ``--remove-solvent`` must strip these from the
        # merged settings.
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_no_solvent_options_leaves_project_sp_settings_unchanged(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """No solvent CLI flags leave the project sp solvent settings intact."""
        # ``gas_solv`` project sp has smd/cyclohexane; no CLI flags → preserved.
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "cyclohexane"


class TestORCASolventCLIOptCommand:
    """CLI solvent options propagated to the ``opt`` subcommand."""

    def test_solvent_model_and_id_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` sets solvent on the opt job settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
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
                "water",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAOptJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_options_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water -so 'Epsilon 78.36'`` propagates to opt."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
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
                "water",
                "-so",
                "Epsilon 78.36",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_remove_solvent_clears_solvent_from_opt(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
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
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_subcommand_level_solvent_overrides_group_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Subcommand-level ``-sm``/``-si`` overrides group-level solvent."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
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
                "toluene",
                "opt",
                "-sm",
                "smd",
                "-si",
                "water",
            ],
        )

        assert result.exit_code == 0, result.output
        # Subcommand-level smd/water overrides group-level cpcm/toluene
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestORCASolventCLITsCommand:
    """CLI solvent options propagated to the ``ts`` subcommand."""

    def test_solvent_model_and_id_injected_into_ts_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` sets solvent on the ts job settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.ts.ORCATSJob",
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
                "water",
                "ts",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCATSJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_remove_solvent_clears_solvent_from_ts(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` removes solvent from ts job."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.ts.ORCATSJob",
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
                "ts",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestORCACpcmBlockOptions:
    """Tests for the ORCA-specific ``%cpcm`` block options via CLI ``-so``."""

    def test_custom_epsilon_no_solvent_id(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Custom dielectric via ``--remove-solvent`` + ``sp -sm cpcm -so 'Epsilon 78.36'``.

        The project-level solvent (cyclohexane) is cleared by ``--remove-solvent``
        at the group level; the subcommand-level flags then set the custom
        dielectric without a named solvent.
        """
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "sp",
                "-sm",
                "cpcm",
                "-so",
                "Epsilon 78.36",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id is None
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_custom_epsilon_and_refrac(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Custom Epsilon + Refrac via ``--remove-solvent`` + ``sp -sm cpcm -so '...'``.

        Both ``Epsilon`` and ``Refrac`` are passed as a newline-separated
        string to ``-so``; each should appear in ``additional_solvent_options``.
        """
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "sp",
                "-sm",
                "cpcm",
                "-so",
                "Epsilon 78.36\nRefrac 1.33",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id is None
        assert "Epsilon 78.36" in settings.additional_solvent_options
        assert "Refrac 1.33" in settings.additional_solvent_options

    def test_smd_with_surface_type(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm smd -si water -so 'SurfaceType gepol_ses'`` stores all options."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "SurfaceType gepol_ses",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "SurfaceType gepol_ses"

    def test_smd_with_rsolv(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm smd -si water -so 'Rsolv 1.30'`` stores Rsolv option."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
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
                "Rsolv 1.30",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Rsolv 1.30"
