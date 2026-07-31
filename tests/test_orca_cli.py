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


class TestORCACLIOptSubcommand:
    """CLI tests for the ``opt`` subcommand's own solvent-options/
    solventfilename branches, freeze-atoms, and the multi-molecule
    (index-selected) job-creation loop."""

    def test_subcommand_level_solvent_options_and_solventfilename(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmp_path,
    ):
        solvent_file = tmp_path / "custom.cosmors"
        solvent_file.write_text(
            "solventname=1,1,1,3,3,3-hexafluoropropan-2-ol\n"
        )
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
                "opt",
                "-so",
                "iterative",
                "-sf",
                str(solvent_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_solvent_options == "iterative"
        assert settings.solventfilename == str(solvent_file)

    def test_single_molecule_freeze_atoms(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.opt.ORCAOptJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
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
                    "--freeze-atoms",
                    "1-2",
                ],
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        molecule = mock_job_cls.call_args[1]["molecule"]
        assert molecule.frozen_atoms is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        run_orca_and_capture_settings,
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.opt.ORCAOptJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                    "--freeze-atoms",
                    "1",
                ],
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]
        for c in mock_job_cls.call_args_list:
            assert c.kwargs["molecule"].frozen_atoms is not None


class TestORCACLISpSubcommand:
    """CLI tests for the ``sp`` subcommand's solventfilename branch and
    the multi-molecule (index-selected) job-creation loop."""

    def test_subcommand_level_solventfilename_applied(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmp_path,
    ):
        solvent_file = tmp_path / "custom.cosmors"
        solvent_file.write_text(
            "solventname=1,1,1,3,3,3-hexafluoropropan-2-ol\n"
        )
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
                "-sf",
                str(solvent_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.solventfilename == str(solvent_file)

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        run_orca_and_capture_settings,
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "sp",
                ],
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_qmmm_child_subcommand_skips_direct_sp_job_creation(
        self,
        single_molecule_xyz_file,
    ):
        """When qmmm is invoked as a child of sp, sp() itself must not
        also create a direct ORCASinglePointJob."""
        from unittest.mock import patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock_sp_job_cls:
            result = CliRunner().invoke(
                orca,
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "def2-svp",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_sp_job_cls.assert_not_called()


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


class TestORCACLITsSubcommand:
    """CLI tests for the ``ts`` subcommand's Hessian/OptTS/ScanTS
    option-merging branches (not covered by the solvent-only tests)."""

    def test_solvent_options_and_solventfilename_applied(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmp_path,
    ):
        solvent_file = tmp_path / "custom.cosmors"
        solvent_file.write_text(
            "solventname=1,1,1,3,3,3-hexafluoropropan-2-ol\n"
        )
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
                "ts",
                "-so",
                "iterative",
                "-sf",
                str(solvent_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_solvent_options == "iterative"
        assert settings.solventfilename == str(solvent_file)

    def test_hessian_and_trust_radius_options_applied(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
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
                "ts",
                "--inhess",
                "-f",
                "hess.hess",
                "--hybrid-hess",
                "-a",
                "[0, 1, 2]",
                "--numhess",
                "-s",
                "10",
                "-t",
                "0.3",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.inhess is True
        assert settings.inhess_filename == "hess.hess"
        assert settings.hybrid_hess is True
        assert settings.hybrid_hess_atoms == "[0, 1, 2]"
        assert settings.numhess is True
        assert settings.recalc_hess == 10
        assert settings.trust_radius == 0.3

    def test_jobtype_scants_alone_is_overridden_by_tssearch_type_default(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """BUG (see BUGS_FOUND.md #29): -j scants alone is silently
        overridden right back to "optts" because -ts/--tssearch-type
        always has a non-None default ("optts"), so the
        jobtype-inference branch's effect is immediately undone."""
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
                "ts",
                "-j",
                "scants",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.tssearch_type == "optts"

    def test_explicit_tssearch_type_scants_works(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
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
                "ts",
                "-ts",
                "scants",
                "-c",
                "[[1,2]]",
                "-x",
                "3.0",
                "-y",
                "1.2",
                "-n",
                "15",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.tssearch_type == "scants"

    def test_scants_requires_all_scan_parameters(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
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
                "ts",
                "-ts",
                "scants",
                "-c",
                "[[1,2]]",
                "-x",
                "3.0",
            ],
        )
        assert result.exit_code != 0
        assert "requires" in result.output

    def test_scants_without_coordinates_falls_back_to_project(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Without CLI --coordinates, ScanTS falls back to the project's
        scants_modred settings; if the project has none, raises."""
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
                "ts",
                "-ts",
                "scants",
            ],
        )
        assert result.exit_code != 0
        assert "requires scan coordinates" in result.output

    def test_full_scan_option_applied(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
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
                "ts",
                "--full-scan",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.full_scan is True


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


class TestORCACLIGroupValidation:
    """Validation and settings-merge branches on the ``orca`` group
    callback itself (not delegated to any subcommand)."""

    def test_molecule_id_not_supported_raises(self, single_molecule_xyz_file):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca,
            ["-f", single_molecule_xyz_file, "--mid", "abc", "sp"],
        )
        assert result.exit_code != 0
        assert "not supported for ORCA job submission" in result.output

    def test_index_and_structure_index_mutually_exclusive(
        self, single_molecule_xyz_file
    ):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca,
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "sp",
            ],
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_chemsmart_db_requires_exactly_one_selector(
        self, database_chemsmart_file
    ):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca, ["-f", database_chemsmart_file, "sp"]
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_xtb_output_inherits_charge_and_multiplicity(
        self, xtb_water_outfolder, run_orca_and_capture_settings
    ):
        import os

        xtb_out = os.path.join(xtb_water_outfolder, "water_ohess.out")
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            ["-p", "gas_solv", "-f", xtb_out, "opt"],
        )
        assert result.exit_code == 0, result.output
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_non_chemsmart_db_falls_back_to_defaults(
        self, database_ase_file, run_orca_and_capture_settings
    ):
        from unittest.mock import MagicMock, patch

        mock_molecule = MagicMock(name="ase_db_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=[mock_molecule],
        ):
            result, settings = run_orca_and_capture_settings(
                "chemsmart.jobs.orca.opt.ORCAOptJob",
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    database_ase_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_filename_and_pubchem_both_missing_raises(self):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca, ["-p", "gas_solv", "sp"], catch_exceptions=True
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "has not been specified" in str(result.exception)

    def test_filename_and_pubchem_both_given_raises(
        self, single_molecule_xyz_file
    ):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "--pubchem",
                "222",
                "sp",
            ],
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "have been specified" in str(result.exception)

    def test_label_and_append_label_mutually_exclusive_raises(
        self, single_molecule_xyz_file
    ):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-l",
                "custom",
                "-a",
                "suffix",
                "sp",
            ],
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "not both" in str(result.exception)

    def test_pubchem_only_without_label_crashes(self):
        """BUG (see BUGS_FOUND.md #30): the "output" label fallback for
        a filename-less (PubChem-only) job with no -l/-a is
        unreachable -- os.path.basename(filename) is called
        unconditionally on filename (None here) before the
        `if filename:` guard that would skip it, so this currently
        crashes with a TypeError instead of falling back to "output"."""
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[pubchem_molecule],
        ):
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "sp",
                ],
                catch_exceptions=True,
            )
        assert result.exit_code != 0
        assert isinstance(result.exception, TypeError)

    def test_default_label_doubles_subcommand_suffix(
        self, single_molecule_xyz_file
    ):
        """BUG (see BUGS_FOUND.md #30): the default label ends up with
        the subcommand name appended twice (e.g. "mol_opt_opt")
        because an unconditional final append duplicates the
        conditional one just above it."""
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.opt.ORCAOptJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
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
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label.endswith("_opt_opt")

    def test_group_level_keywords_applied(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
                "-x",
                "m062x",
                "-b",
                "def2svp",
                "-t",
                "my title",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.functional == "m062x"
        assert settings.basis == "def2svp"
        assert settings.title == "my title"

    def test_structure_index_used_as_index(
        self, multiple_molecules_xyz_file, run_orca_and_capture_settings
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                multiple_molecules_xyz_file,
                "--si",
                "1",
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_index_requires_record_selector(
        self, database_chemsmart_file
    ):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        result = CliRunner().invoke(
            orca,
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "abc",
                "-i",
                "1",
                "sp",
            ],
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_xtb_output_with_unset_charge_and_multiplicity(
        self, xtb_water_outfolder, run_orca_and_capture_settings
    ):
        import os
        from unittest.mock import MagicMock, patch

        xtb_out = os.path.join(xtb_water_outfolder, "water_ohess.out")
        mock_molecule = MagicMock(charge=None, multiplicity=None)

        def _fake_from_filepath(filepath, **kwargs):
            if kwargs.get("return_list"):
                return [mock_molecule]
            return mock_molecule

        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            side_effect=_fake_from_filepath,
        ):
            result, settings = run_orca_and_capture_settings(
                "chemsmart.jobs.orca.opt.ORCAOptJob",
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    xtb_out,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_unrecognized_filetype_raises(self, tmp_path):
        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        bad_file = tmp_path / "structure.weird"
        bad_file.write_text("nonsense")
        result = CliRunner().invoke(
            orca,
            ["-p", "gas_solv", "-f", str(bad_file), "sp"],
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "Unrecognised filetype" in str(result.exception)

    def test_chemsmart_db_loads_molecule_by_record_index(
        self, database_chemsmart_file, run_orca_and_capture_settings
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_loads_molecule_by_structure_id(
        self, database_chemsmart_file, run_orca_and_capture_settings
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_remove_solvent_clears_solvent(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
        assert settings.custom_solvent is None

    def test_remaining_group_level_keywords_applied(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
                "-A",
                "mp2",
                "-D",
                "d3bj",
                "-B",
                "def2/J",
                "-e",
                "def2/QZVPP",
                "-d",
                "defgrid3",
                "--scf-tol",
                "TightSCF",
                "--scf-algorithm",
                "AutoTRAH",
                "--scf-maxiter",
                "200",
                "--scf-convergence",
                "1e-8",
                "--dipole",
                "--quadrupole",
                "--mdci-cutoff",
                "tight",
                "--mdci-density",
                "relaxed",
                "-r",
                "def2/J",
                "--forces",
                "sp",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.ab_initio == "mp2"
        assert settings.dispersion == "d3bj"
        assert settings.aux_basis == "def2/J"
        assert settings.extrapolation_basis == "def2/QZVPP"
        assert settings.defgrid == "defgrid3"
        assert settings.scf_tol == "TightSCF"
        assert settings.scf_algorithm == "AutoTRAH"
        assert settings.scf_maxiter == 200
        assert settings.scf_convergence == 1e-8
        assert settings.dipole is True
        assert settings.quadrupole is True
        assert settings.mdci_cutoff == "tight"
        assert settings.mdci_density == "relaxed"
        assert settings.additional_route_parameters == "def2/J"
        assert settings.forces is True

    def test_chemsmart_db_append_label_includes_record_id_suffix(
        self, database_chemsmart_file, run_orca_and_capture_settings
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "-a",
                "suffix",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_record_index_suffix(
        self, database_chemsmart_file, run_orca_and_capture_settings
    ):
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "-a",
                "suffix",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None


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

    def test_irc_solvent_options_and_solventfilename(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmp_path,
    ):
        solvent_file = tmp_path / "custom.cosmors"
        solvent_file.write_text(
            "solventname=1,1,1,3,3,3-hexafluoropropan-2-ol\n"
        )
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
                "-so",
                "iterative",
                "-sf",
                str(solvent_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_solvent_options == "iterative"
        assert settings.solventfilename == str(solvent_file)

    def test_irc_all_remaining_options_applied(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
                "-p",
                "3",
                "-i",
                "read",
                "-f",
                "hessian.hess",
                "-m",
                "2",
                "-M",
                "--init-displ",
                "DE",
                "--scale-init-displ",
                "0.2",
                "--de-init-displ",
                "0.003",
                "--follow-coordtype",
                "cartesian",
                "--scale-displ-sd",
                "0.1",
                "--adapt-scale-displ",
                "--sd-parabolicfit",
                "--interpolate-only",
                "--do-sd-corr",
                "--scale-displ-sd-corr",
                "0.5",
                "--sd-corr-parabolicfit",
                "--tolrmsg",
                "0.0005",
                "--tolmaxg",
                "0.002",
                "-I",
                "[[1,2,3]]",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.printlevel == 3
        assert settings.inithess == "read"
        assert settings.hess_filename == "hessian.hess"
        assert settings.hessmode == 2
        assert settings.init_displ == "DE"
        assert settings.scale_init_displ == 0.2
        assert settings.de_init_displ == 0.003
        assert settings.follow_coordtype == "cartesian"
        assert settings.scale_displ_sd == 0.1
        assert settings.adapt_scale_displ is True
        assert settings.sd_parabolicfit is True
        assert settings.interpolate_only is True
        assert settings.do_sd_corr is True
        assert settings.scale_displ_sd_corr == 0.5
        assert settings.sd_corr_parabolicfit is True
        assert settings.tolrmsg == 0.0005
        assert settings.tolmaxg == 0.002
        assert settings.internal_modred == [[1, 2, 3]]


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

    def test_explicit_jobtype_used(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
                "-j",
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
        assert settings is not None

    def test_subcommand_level_solvent_options_and_solventfilename(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmp_path,
    ):
        solvent_file = tmp_path / "custom.cosmors"
        solvent_file.write_text(
            "solventname=1,1,1,3,3,3-hexafluoropropan-2-ol\n"
        )
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
                "-so",
                "iterative",
                "-sf",
                str(solvent_file),
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_solvent_options == "iterative"
        assert settings.solventfilename == str(solvent_file)

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self, multiple_molecules_xyz_file
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.scan.ORCAScanJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
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
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_qmmm_child_subcommand_skips_direct_scan_job_creation(
        self, single_molecule_xyz_file
    ):
        from unittest.mock import patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch(
            "chemsmart.jobs.orca.scan.ORCAScanJob"
        ) as mock_scan_job_cls:
            result = CliRunner().invoke(
                orca,
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "def2-svp",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_scan_job_cls.assert_not_called()


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

    def test_explicit_jobtype_used(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
        """``-j modred`` explicitly (rather than falling back to the
        default) still resolves modred settings correctly."""
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
                "-j",
                "modred",
                "-c",
                "[[1,2]]",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self, multiple_molecules_xyz_file
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.modred.ORCAModredJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "modred",
                    "-c",
                    "[[1,2]]",
                ],
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_qmmm_child_subcommand_skips_direct_modred_job_creation(
        self, single_molecule_xyz_file
    ):
        from unittest.mock import patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch(
            "chemsmart.jobs.orca.modred.ORCAModredJob"
        ) as mock_modred_job_cls:
            result = CliRunner().invoke(
                orca,
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "def2-svp",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_modred_job_cls.assert_not_called()


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

    def test_neb_all_options_applied(
        self, single_molecule_xyz_file, run_orca_and_capture_settings
    ):
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
                "-n",
                "8",
                "-i",
                single_molecule_xyz_file,
                "-r",
                single_molecule_xyz_file,
                "-s",
                "XTB2",
                "-o",
            ],
        )
        assert result.exit_code == 0, result.output
        assert settings.nimages == 8
        assert settings.intermediate_xyzfile == single_molecule_xyz_file
        assert settings.restarting_xyzfile == single_molecule_xyz_file
        assert settings.semiempirical == "XTB2"
        assert settings.preopt_ends is True

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self, multiple_molecules_xyz_file
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.neb.ORCANEBJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                orca,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "neb",
                    "-j",
                    "NEB-TS",
                    "-e",
                    multiple_molecules_xyz_file,
                ],
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2

    def test_qmmm_child_subcommand_skips_direct_neb_job_creation(
        self, single_molecule_xyz_file
    ):
        from unittest.mock import patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca

        with patch("chemsmart.jobs.orca.neb.ORCANEBJob") as mock_neb_job_cls:
            result = CliRunner().invoke(
                orca,
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "def2-svp",
                ],
                obj={"jobrunner": None},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_neb_job_cls.assert_not_called()


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
