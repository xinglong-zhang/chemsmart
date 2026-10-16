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

from unittest.mock import MagicMock


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
    def test_default_label_appends_subcommand_once(
        self, single_molecule_xyz_file
    ):
        from os.path import basename, splitext
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch("chemsmart.jobs.orca.opt.ORCAOptJob") as mock:
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
                    "opt",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        base = splitext(basename(single_molecule_xyz_file))[0]
        assert mock.call_args.kwargs["label"] == f"{base}_opt"

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
        assert mock.call_args.kwargs["settings"].aux_basis == "def2/J"


class TestORCAQMMMCLIjobtypeOverride:
    def test_cli_qm_qm2_uses_three_layer_yaml_lot(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """CLI -j QM/QM2 keeps YAML per-layer theory and drops the MM layer."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "test_qmmm_3layer",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-j",
                "QM/QM2",
                "-ha",
                "1-3",
                "-ct",
                "0",
                "-mt",
                "2",
                "-ch",
                "0",
                "-mh",
                "2",
            ],
            ctx_obj={"jobrunner": MagicMock()},
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.jobtype == "QM/QM2"
        assert settings.high_level_functional == "B3LYP"
        assert settings.high_level_basis == "def2-SVP"
        assert settings.intermediate_level_functional == "HF"
        assert settings.intermediate_level_basis == "STO-3G"
        assert settings.intermediate_level_atoms is None
        assert settings.low_level_method is None
        route = settings.qmmm_route_string
        assert "QM/QM2" in route
        assert "/MM" not in route
        qmmm_block = settings.qmmm_block
        assert "QM2Atoms" not in qmmm_block
        assert "ORCAFFFilename" not in qmmm_block
        assert "Charge_Total 0" in qmmm_block
        assert "Mult_Total 2" in qmmm_block

    def test_cli_qm_qm2_keeps_cli_intermediate_atoms(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Explicit -ia is kept when overriding QM/QM2/MM YAML to QM/QM2."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "test_qmmm_3layer",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-j",
                "QM/QM2",
                "-ha",
                "1-3",
                "-ia",
                "4-6",
                "-ct",
                "0",
                "-mt",
                "2",
                "-ch",
                "0",
                "-mh",
                "2",
            ],
            ctx_obj={"jobrunner": MagicMock()},
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.jobtype == "QM/QM2"
        assert settings.intermediate_level_atoms == "4-6"
        assert settings.low_level_method == "system.ORCAFF.prms"


class TestORCAQMMMCLIHighLevelHBondLength:
    def test_cli_accepts_dict_string(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "opt",
                "qmmm",
                "-j",
                "QMMM",
                "-lm",
                "system.ORCAFF.prms",
                "-ha",
                "1-3",
                "-ch",
                "0",
                "-mh",
                "1",
                "-h",
                "{'C_H': 1.09, 'N_H': 1.01}",
            ],
            ctx_obj={"jobrunner": MagicMock()},
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.low_level_method == "system.ORCAFF.prms"
        assert settings.high_level_h_bond_length == {
            "C_H": 1.09,
            "N_H": 1.01,
        }
        h_block = settings._get_h_bond_length()
        assert "Dist_C_HLA 1.09" in h_block
        assert "Dist_N_HLA 1.01" in h_block
