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
        assert (
            settings is not None
        ), "ORCASinglePointJob was never instantiated"
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
        assert (
            settings is not None
        ), "ORCASinglePointJob was never instantiated"
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

    def test_solventfilename_injected_into_sp_settings(
        self,
        tmp_path,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """-sf /path/water.cosmorsxyz stores solventfilename on sp settings."""
        # Create a dummy .cosmorsxyz file so click.Path(exists=True) is satisfied
        sf = tmp_path / "water.cosmorsxyz"
        sf.write_text("")

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
                "cosmors",
                "-si",
                "water",
                "-sf",
                str(sf),
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cosmors"
        assert settings.solvent_id == "water"
        assert settings.solventfilename == str(sf)

    def test_solventfilename_group_level_injected_into_sp_settings(
        self,
        tmp_path,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """-sf at the orca group level propagates solventfilename to sp settings."""
        sf = tmp_path / "custom.cosmorsxyz"
        sf.write_text("")

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
                "cosmors",
                "-sf",
                str(sf),
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solventfilename == str(sf)


class TestORCALabelAndAuxBasisOptions:
    def test_short_a_sets_append_label(self, single_molecule_xyz_file):
        from os.path import basename, splitext
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock:
            mock.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-a",
                    "tag",
                    "sp",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock.call_args is not None
        base = splitext(basename(single_molecule_xyz_file))[0]
        assert mock.call_args.kwargs["label"].startswith(f"{base}_tag")
        assert mock.call_args.kwargs["settings"].aux_basis is None

    def test_short_B_sets_aux_basis(self, single_molecule_xyz_file):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock:
            mock.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-B",
                    "def2/J",
                    "sp",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock.call_args is not None


class TestORCACLIInpCommand:
    """CLI tests for the ``inp`` (run input file as-is) subcommand."""

    def test_inp_job_creation_from_inp_file(self, water_sp_input_path):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch("chemsmart.jobs.orca.job.ORCAInpJob") as mock_job_cls:
            mock_job_cls.from_filename.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    water_sp_input_path,
                    "inp",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        mock_job_cls.from_filename.assert_called_once()
        assert (
            mock_job_cls.from_filename.call_args.kwargs["filename"]
            == water_sp_input_path
        )


class TestORCACLIIrcCommand:
    """CLI tests for the ``irc`` subcommand."""

    def test_irc_basic_job_creation(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``irc`` subcommand creates an ``ORCAIRCJob``."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.irc.ORCAIRCJob",
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
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAIRCJob was never instantiated"

    def test_irc_direction_and_maxiter_options(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``-d forward --maxiter 50`` set direction and maxiter on IRC settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.irc.ORCAIRCJob",
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
                "--maxiter",
                "50",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.direction == "forward"
        assert settings.maxiter == 50


class TestORCACLIScanCommand:
    """CLI tests for the ``scan`` subcommand group."""

    def test_scan_requires_full_coordinate_spec(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``scan`` without ``-x``/``-y``/``-n`` raises an assertion error."""
        import pytest

        with pytest.raises(AssertionError, match="Scanning coordinates"):
            run_orca_and_capture_settings(
                "chemsmart.jobs.orca.scan.ORCAScanJob",
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
                    "[[2,3]]",
                ],
            )

    def test_scan_basic_job_creation(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``scan -c [[2,3]] -x 3.0 -y 1.2 -n 15`` creates an ``ORCAScanJob``."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.scan.ORCAScanJob",
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
                "[[2,3]]",
                "-x",
                "3.0",
                "-y",
                "1.2",
                "-n",
                "15",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAScanJob was never instantiated"

    def test_scan_constrained_coordinates_option(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``-cc`` adds additional modredundant constraints to the scan settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.scan.ORCAScanJob",
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
                "[[2,3]]",
                "-x",
                "3.0",
                "-y",
                "1.2",
                "-n",
                "15",
                "-cc",
                "[[1,2,3]]",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.modred["constrained_coordinates"] == [[1, 2, 3]]


class TestORCACLIModredCommand:
    """CLI tests for the ``modred`` subcommand group."""

    def test_modred_requires_coordinates(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``modred`` without coordinates raises an assertion error."""
        import pytest

        with pytest.raises(
            AssertionError, match="Coordinates must be provided"
        ):
            run_orca_and_capture_settings(
                "chemsmart.jobs.orca.modred.ORCAModredJob",
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "modred",
                ],
            )

    def test_modred_basic_job_creation(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``modred -c "[[1,2]]"`` creates an ``ORCAModredJob``."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.modred.ORCAModredJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "modred",
                "-c",
                "[[1,2]]",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAModredJob was never instantiated"


class TestORCACLIQmmmSubcommand:
    """CLI tests for the ``qmmm`` subcommand attached to ``opt``."""

    def test_opt_qmmm_job_creation(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``opt qmmm`` creates an ``ORCAQMMMJob`` inheriting opt jobtype."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "qmmm",
                "-hx",
                "b3lyp",
                "-hb",
                "def2-svp",
            ],
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAQMMMJob was never instantiated"


class TestORCACLINebCommand:
    """CLI tests for the ``neb`` subcommand group."""

    def test_neb_basic_job_creation(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``neb -j NEB-TS -e <file>`` creates an ``ORCANEBJob``."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.neb.ORCANEBJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "neb",
                "-j",
                "NEB-TS",
                "-e",
                single_molecule_xyz_file,
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCANEBJob was never instantiated"
        assert settings.joboption == "NEB-TS"
        assert settings.ending_xyzfile == single_molecule_xyz_file


class TestORCACLIQrcCommand:
    """CLI tests for the ``qrc`` subcommand group."""

    def test_qrc_default_jobtype_opt(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``qrc`` with no ``-j`` defaults to the ``opt`` jobtype settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qrc.ORCAQRCJob",
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
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAQRCJob was never instantiated"

    def test_qrc_explicit_ts_jobtype(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``qrc -j ts`` uses TS settings from the project for the QRC job."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qrc.ORCAQRCJob",
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
            {"jobrunner": None},
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAQRCJob was never instantiated"
