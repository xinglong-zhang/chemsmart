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


class TestGaussianCLIMecpCommand:
    """CLI tests for the ``mecp`` subcommand."""

    def test_mecp_defaults_state_b_to_plus_two_multiplicity(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Default ``multiplicity-b`` is inferred as ``multiplicity-a + 2``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.mecp.GaussianMECPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "mecp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.multiplicity_a == 1
        assert settings.multiplicity_b == 3
        assert settings.charge_a == 0
        assert settings.charge_b == 0

    def test_mecp_custom_spin_and_charge_options_are_forwarded(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """MECP CLI options are forwarded to the GaussianMECPJob constructor."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.mecp.GaussianMECPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "1",
                "-m",
                "2",
                "mecp",
                "--multiplicity-a",
                "2",
                "--multiplicity-b",
                "4",
                "--charge-a",
                "1",
                "--charge-b",
                "1",
                "--max-steps",
                "120",
                "--step-size",
                "0.08",
                "--trust-radius",
                "0.12",
                "--energy-diff-tol",
                "2e-4",
                "--force-max-tol",
                "8e-4",
                "--force-rms-tol",
                "6e-4",
                "--disp-max-tol",
                "2.0e-3",
                "--disp-rms-tol",
                "1.5e-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.multiplicity_a == 2
        assert settings.multiplicity_b == 4
        assert settings.charge_a == 1
        assert settings.charge_b == 1
        assert settings.max_steps == 120
        assert settings.step_size == 0.08
        assert settings.trust_radius == 0.12
        assert settings.energy_diff_tol == 2e-4
        assert settings.force_max_tol == 8e-4
        assert settings.force_rms_tol == 6e-4
        assert settings.disp_max_tol == 2e-3
        assert settings.disp_rms_tol == 1.5e-3

    def test_mecp_adaptive_step_size_options_are_forwarded(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Adaptive step size CLI options are forwarded to GaussianMECPJob settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.mecp.GaussianMECPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "mecp",
                "--no-adaptive-step-size",
                "--step-size-grow",
                "1.2",
                "--step-size-shrink",
                "0.6",
                "--step-size-min",
                "1e-3",
                "--step-size-max",
                "0.5",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.adaptive_step_size is False
        assert settings.step_size_grow == 1.2
        assert settings.step_size_shrink == 0.6
        assert settings.step_size_min == 1e-3
        assert settings.step_size_max == 0.5

    def test_mecp_adaptive_step_size_enabled_by_default(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Adaptive step size is enabled by default with Barzilai-Borwein method."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.mecp.GaussianMECPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "mecp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.adaptive_step_size is True
        assert settings.step_size_method == "bb"
        assert settings.step_size_grow == 1.2
        assert settings.step_size_shrink == 0.833
        assert settings.step_size_min == 1e-4
        assert settings.step_size_max == 1.0

    def test_mecp_step_size_method_option_is_forwarded(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--step-size-method`` CLI option is forwarded to GaussianMECPJob settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.mecp.GaussianMECPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "mecp",
                "--step-size-method",
                "grow_shrink",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings.step_size_method == "grow_shrink"
