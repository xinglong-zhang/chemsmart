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

import sys
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

import chemsmart.cli.gaussian.traj  # noqa: F401 - ensures module is cached
from chemsmart.cli.gaussian.gaussian import gaussian

# ``chemsmart.cli.gaussian.__init__`` does ``from .traj import traj``,
# which shadows the ``traj`` submodule attribute on the package with the
# click Command object. Fetch the real submodule from ``sys.modules``
# directly so we can patch ``GaussianTrajJob`` on it.
traj_module = sys.modules["chemsmart.cli.gaussian.traj"]


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


class TestGaussianCLIComCommand:
    """CLI tests for the ``com`` (run input file as-is) subcommand."""

    def test_com_job_creation_from_com_file(
        self,
        gaussian_opt_inputfile,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``com`` subcommand creates a ``GaussianComJob`` from a ``.com`` file."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.job.GaussianComJob",
            [
                "-p",
                "gas_solv",
                "-f",
                gaussian_opt_inputfile,
                "com",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianComJob was never instantiated"
        assert settings.input_string is not None


class TestGaussianCLIRespCommand:
    """CLI tests for the ``resp`` subcommand."""

    def test_resp_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``resp`` subcommand creates a ``GaussianRESPJob`` with fixed route."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.resp.GaussianRESPJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "resp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianRESPJob was never instantiated"
        assert "Pop=MK" in settings.route_to_be_written


class TestGaussianCLIUserjobCommand:
    """CLI tests for the ``userjob`` (custom route) subcommand."""

    def test_userjob_requires_route_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``userjob`` fails without the required ``-r`` route option."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.custom.GaussianCustomJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "userjob",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code != 0

    def test_userjob_job_creation_with_route(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``userjob -r <route>`` creates a ``GaussianCustomJob`` with that route."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.custom.GaussianCustomJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "userjob",
                "-r",
                "opt freq",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.route_to_be_written == "opt freq"

    def test_userjob_append_info_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-a`` appends decoded additional info to the job settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.custom.GaussianCustomJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "userjob",
                "-r",
                "opt",
                "-a",
                "extra info",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.append_additional_info == "extra info"


class TestGaussianCLINciCommand:
    """CLI tests for the ``nci`` (non-covalent interaction) subcommand."""

    def test_nci_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``nci`` subcommand creates a ``GaussianNCIJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.nci.GaussianNCIJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "nci",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianNCIJob was never instantiated"

    def test_nci_group_level_solvent_injected(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """Group-level solvent options are propagated to ``nci`` settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.nci.GaussianNCIJob",
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
                "nci",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestGaussianCLIWbiCommand:
    """CLI tests for the ``wbi`` (Wiberg bond index) subcommand."""

    def test_wbi_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``wbi`` subcommand creates a ``GaussianWBIJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.wbi.GaussianWBIJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "wbi",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianWBIJob was never instantiated"


class TestGaussianCLILinkCommand:
    """CLI tests for the ``link`` subcommand."""

    def test_link_requires_jobtype(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``link`` without ``-j`` raises since jobtype is required."""
        with pytest.raises(ValueError, match="Jobtype must be provided"):
            run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.link.GaussianLinkJob",
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "link",
                ],
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )

    def test_link_opt_jobtype_creates_job(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``link -j opt`` creates a ``GaussianLinkJob`` with unrestricted functional."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.link.GaussianLinkJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "link",
                "-j",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.functional.lower().startswith("u")
        assert settings.stable == "opt"
        assert settings.guess == "mix"

    def test_link_irc_jobtype_sets_irc_params_and_label(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``link -j irc -d forward`` sets IRC parameters on the link settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.link.GaussianLinkJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "link",
                "-j",
                "irc",
                "-d",
                "forward",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.direction == "forward"
        assert settings.maxpoints == 512

    def test_link_custom_route_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--route`` sets a custom route for the link section."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.link.GaussianLinkJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "link",
                "-j",
                "opt",
                "--route",
                "opt=(calcfc)",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.link_route == "opt=(calcfc)"


class TestGaussianCLIModredCommand:
    """CLI tests for the ``modred`` subcommand group."""

    def test_modred_requires_coordinates(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``modred`` without coordinates raises an assertion error."""
        with pytest.raises(
            AssertionError, match="Coordinates must be provided"
        ):
            run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.modred.GaussianModredJob",
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
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )

    def test_modred_basic_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``modred -c "[[1,2]]"`` creates a ``GaussianModredJob``."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.modred.GaussianModredJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianModredJob was never instantiated"

    def test_modred_qmmm_subcommand_creates_qmmm_job(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``modred -c "[[1,2]]" qmmm`` creates a ``GaussianQMMMJob`` inheriting modred info."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob",
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
                "6-31g(d)",
                "-ha",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianQMMMJob was never instantiated"
        assert settings.high_level_functional == "b3lyp"
        assert settings.high_level_basis == "6-31g(d)"
        assert settings.jobtype == "modred"


class TestGaussianCLITrajCommand:
    """CLI tests for the ``traj`` subcommand."""

    def test_traj_basic_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``traj -j opt`` creates a ``GaussianTrajJob`` for structures from a trajectory."""
        runner = CliRunner()
        with patch.object(traj_module, "GaussianTrajJob") as mock_job_cls:
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
                    "traj",
                    "-j",
                    "opt",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()

    def test_traj_mutually_exclusive_options_raise(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        """``-g`` and ``-ns`` together raise a usage error."""
        runner = CliRunner()
        with patch.object(traj_module, "GaussianTrajJob") as mock_job_cls:
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
                    "traj",
                    "-j",
                    "opt",
                    "-g",
                    "rmsd",
                    "-ns",
                    "5",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )
        assert result.exit_code != 0


class TestGaussianCLIDiasCommand:
    """CLI tests for the ``dias`` (distortion/interaction) subcommand."""

    def test_dias_basic_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``dias -i 1-3`` creates a ``GaussianDIASJob`` with solvent removed by default."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.dias.GaussianDIASJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "dias",
                "-i",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None, "GaussianDIASJob was never instantiated"

    def test_dias_solv_flag_keeps_solvent(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``--solv`` keeps solvent on the DI-AS job settings."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.dias.GaussianDIASJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "dias",
                "-i",
                "1-3",
                "--solv",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"

    def test_dias_ts_mode_option(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-m ts`` selects the TS mode for DI-AS analysis."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.dias.GaussianDIASJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "dias",
                "-i",
                "1-3",
                "-m",
                "ts",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
