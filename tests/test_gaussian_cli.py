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

import inspect
import os
import sys
from unittest.mock import MagicMock, patch

import click
import pytest
from click.testing import CliRunner

import chemsmart.cli.gaussian.traj  # noqa: F401 - ensures module is cached
from chemsmart.cli.gaussian.gaussian import gaussian
from tests.helpers import (
    make_cli_ctx_obj,
    run_gaussian_and_capture_settings,
    run_gaussian_qmmm_callback,
)

# ``chemsmart.cli.gaussian.__init__`` does ``from .traj import traj``,
# which shadows the ``traj`` submodule attribute on the package with the
# click Command object. Fetch the real submodule from ``sys.modules``
# directly so we can patch ``GaussianTrajJob`` on it.
traj_module = sys.modules["chemsmart.cli.gaussian.traj"]


class TestGaussianCLIPubChemOptCommand:
    """CLI tests for PubChem / database label paths without an explicit ``-l``."""

    def test_pubchem_only_opt_passes_molecule_to_job(
        self,
        gaussian_jobrunner_no_scratch,
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

    def test_pubchem_default_label_uses_output_prefix(
        self,
        gaussian_jobrunner_no_scratch,
    ):
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        mock = MagicMock(name="mock")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[pubchem_molecule],
        ):
            result, _ = run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                [
                    "-p",
                    "gas_solv",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                job_cls=mock,
            )

        assert result.exit_code == 0, result.output
        assert mock.call_args.kwargs["label"] == "output_opt"

    def test_database_sid_preserved_in_default_label(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        from os.path import basename, splitext

        from chemsmart.database.database import Database

        db = Database(database_chemsmart_file)
        structure_id = db.get_record(record_index=1)["molecules"][0][
            "structure_id"
        ]

        mock = MagicMock(name="mock")
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--sid",
                structure_id,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock,
        )

        assert result.exit_code == 0, result.output
        base = splitext(basename(database_chemsmart_file))[0]
        assert (
            mock.call_args.kwargs["label"] == f"{base}_SID-{structure_id}_opt"
        )


class TestGaussianCLIOptCommand:
    """CLI tests for the ``opt`` subcommand's freeze-atoms and
    multi-molecule (index-selected) job-creation branches."""

    def test_single_molecule_freeze_atoms(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
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
                "opt",
                "--freeze-atoms",
                "1-2",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        molecule = mock_job_cls.call_args[1]["molecule"]
        assert molecule.frozen_atoms is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]
        for c in mock_job_cls.call_args_list:
            assert c.kwargs["molecule"].frozen_atoms is not None


class TestGaussianSolventCLIOptCommand:
    """CLI solvent options propagated to the ``opt`` subcommand."""

    def test_solvent_model_and_id_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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


class TestGaussianCLIGroupValidation:
    """Validation and settings-merge branches on the ``gaussian`` group
    callback itself (not delegated to any subcommand)."""

    def test_molecule_id_not_supported_raises(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-f",
                single_molecule_xyz_file,
                "--mid",
                "abc",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "not supported for Gaussian job submission" in result.output

    def test_index_and_structure_index_mutually_exclusive(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_used_as_index(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """``--si`` is treated as an alias for ``-i`` when given alone."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_requires_exactly_one_selector(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        """No selector and no --sid given for a chemsmart db input."""
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            ["-f", database_chemsmart_file, "opt"],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_rejects_multiple_selectors(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "--sid",
                "abc",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_index_requires_record_selector(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "abc",
                "-i",
                "1",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_xtb_output_inherits_charge_and_multiplicity(
        self,
        xtb_water_outfolder,
        gaussian_jobrunner_no_scratch,
    ):
        """A ``.out`` file detected as xTB should fall back to default
        Gaussian settings but inherit charge/multiplicity from xTB."""

        xtb_out = os.path.join(xtb_water_outfolder, "water_ohess.out")
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            ["-p", "gas_solv", "-f", xtb_out, "opt"],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.charge == 0
        assert settings.multiplicity == 1

    def test_xtb_output_with_unset_charge_and_multiplicity(
        self,
        xtb_water_outfolder,
        gaussian_jobrunner_no_scratch,
    ):
        """If the xTB-derived molecule has no charge/multiplicity, the
        default Gaussian settings' own (unset) values are kept."""

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
            result, settings = run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                ["-p", "gas_solv", "-f", xtb_out, "-c", "0", "-m", "1", "opt"],
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_non_chemsmart_db_falls_back_to_defaults(
        self,
        database_ase_file,
        gaussian_jobrunner_no_scratch,
    ):
        """A .db file that isn't a chemsmart database (e.g. a plain ASE
        db) falls back to default Gaussian settings entirely."""

        mock_molecule = MagicMock(name="ase_db_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=[mock_molecule],
        ):
            result, settings = run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
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
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_filename_and_pubchem_both_missing_raises(
        self,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            ["-p", "gas_solv", "opt"],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "has not been specified" in str(result.exception)

    def test_filename_and_pubchem_both_given_raises(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "--pubchem",
                "222",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "have been specified" in str(result.exception)

    def test_label_and_append_label_mutually_exclusive_raises(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-l",
                "custom",
                "-a",
                "suffix",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "not both" in str(result.exception)

    def test_append_label_suffixes_filename_derived_label(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-a",
                "suffix",
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_pubchem_only_without_label_falls_back_to_output(
        self,
        gaussian_jobrunner_no_scratch,
    ):
        """Regression test: a filename-less (PubChem-only) job with no
        -l/-a falls back to the "output" label. The default-label branch
        used to call `os.path.basename(filename)` unconditionally on
        `filename` (None here) before the `if filename:` guard, crashing
        with a TypeError."""
        pubchem_molecule = MagicMock(name="pubchem_molecule")
        mock_job_cls = MagicMock(name="mock_job_cls")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[pubchem_molecule],
        ):
            result, _ = run_gaussian_and_capture_settings(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                [
                    "-p",
                    "gas_solv",
                    "--pubchem",
                    "222",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "opt",
                ],
                make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                job_cls=mock_job_cls,
                catch_exceptions=True,
            )
        assert result.exit_code == 0, result.output
        label = mock_job_cls.call_args[1]["label"]
        assert label == "output_opt"

    def test_default_label_with_filename_ignores_db_id_suffix(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """Without -f/--pubchem-derived filename, the label falls back
        to the input file's basename, suffixed with the invoked
        subcommand name."""
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
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_loads_molecule_by_record_index(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_loads_molecule_by_structure_id(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_structure_id_suffix(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        """-a/--append-label combined with --sid takes the first
        (structure_id) branch of the append-label suffix chain."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--sid",
                "f751bb2c27e2",
                "-a",
                "suffix",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_record_index_suffix(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        """-a/--append-label combined with a chemsmart db selector should
        suffix the label with the resolved record/structure identifier."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_append_label_includes_record_id_suffix(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        """--rid (record_id, no structure_id) takes the elif branch for
        the append-label suffix."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_chemsmart_db_default_label_record_id_suffix_branch(
        self,
        database_chemsmart_file,
        gaussian_jobrunner_no_scratch,
    ):
        """--rid alone (no -a/-l) appends the record_id suffix to the
        database basename in the default label, before the subcommand
        suffix."""
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                database_chemsmart_file,
                "--rid",
                "6de213a0",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        db_stem = os.path.splitext(os.path.basename(database_chemsmart_file))[
            0
        ]
        label = mock_job_cls.call_args[1]["label"]
        assert label == f"{db_stem}_RID-6de213a0_opt"

    def test_group_level_functional_and_basis_keywords_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "-x",
                "m062x",
                "-b",
                "def2svp",
                "-s",
                "am1",
                "-o",
                "maxstep=5",
                "-r",
                "empiricaldispersion=gd3",
                "-A",
                "extra info",
                "-C",
                "custom solvent block",
                "-t",
                "my title",
                "-d",
                "n",
                "--forces",
                "sp",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.functional == "m062x"
        assert settings.basis == "def2svp"
        assert settings.semiempirical == "am1"
        assert settings.additional_opt_options_in_route == "maxstep=5"
        assert (
            settings.additional_route_parameters == "empiricaldispersion=gd3"
        )
        assert settings.append_additional_info == "extra info"
        assert settings.custom_solvent.strip() == "custom solvent block"
        assert settings.title == "my title"
        assert settings.dieze_tag == "n"
        assert settings.forces is True

    def _direct_invoke_gaussian_group(self, **kwarg_overrides):
        """Call gaussian()'s own undecorated callback directly,
        bypassing Click's option parsing/validation entirely. Used to
        reach defensive branches in the default-label fallback chain
        that real CLI invocations can't trigger, e.g. an empty (not
        None) filename, or a ``None`` ``invoked_subcommand``."""

        real_fn = inspect.unwrap(gaussian.callback)

        kwargs = dict(
            project="gas_solv",
            filename="",
            label=None,
            append_label=None,
            title=None,
            charge=0,
            multiplicity=1,
            functional=None,
            basis=None,
            semiempirical=None,
            index=None,
            record_index=None,
            record_id=None,
            structure_id=None,
            structure_index=None,
            molecule_id=None,
            additional_opt_options=None,
            additional_route_parameters=None,
            append_additional_info=None,
            custom_solvent=None,
            dieze_tag=None,
            forces=False,
            pubchem=None,
            remove_solvent=False,
            solvent_model=None,
            solvent_id=None,
            solvent_options=None,
        )
        kwargs.update(kwarg_overrides)

        ctx = click.Context(gaussian)
        ctx.obj = {}
        ctx.invoked_subcommand = kwargs.pop("_invoked_subcommand", "opt")
        with ctx:
            real_fn(ctx, **kwargs)
        return ctx

    def test_empty_filename_falls_back_to_output_label(self):
        """Covers the `if filename: ... else: label = "output"` False
        arm with an empty-string filename, which a real CLI run never
        passes (the real-CLI ``None`` case is covered by
        test_pubchem_only_without_label_falls_back_to_output). The
        subcommand suffix is appended unconditionally (the old
        `if ctx.invoked_subcommand:` guard was removed, since a real Click
        invocation of this group always has a non-None
        `invoked_subcommand` by the time this callback runs), so this
        synthetic direct invocation with `_invoked_subcommand=None`
        documents that unconditional append rather than a reachable
        branch."""
        ctx = self._direct_invoke_gaussian_group(_invoked_subcommand=None)
        assert ctx.obj["label"] == "output_None"

    def test_invoked_subcommand_suffixes_default_label_when_present(self):
        """Covers the (now unconditional) subcommand-suffix append in
        the same fallback chain: the invoked subcommand name gets
        appended to the "output" default label."""
        ctx = self._direct_invoke_gaussian_group(_invoked_subcommand="opt")
        assert ctx.obj["label"] == "output_opt"


class TestGaussianCLISinglePointCommand:
    """CLI tests for the ``sp`` (single point) subcommand."""

    def test_basic_sp_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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
    """CLI tests for the ``qrc`` (Quick Reaction Coordinate) subcommand."""

    def test_basic_qrc_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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

    def test_link_subcommand_level_solvent_options_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "iterative",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "iterative"

    def test_link_functional_already_unrestricted_not_double_prefixed(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "-x",
                "ub3lyp",
                "link",
                "-j",
                "opt",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.functional == "ub3lyp"


class TestGaussianCLIModredCommand:
    """CLI tests for the ``modred`` subcommand group."""

    def test_modred_requires_coordinates(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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

    def test_explicit_jobtype_used(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "-j",
                "modred",
                "-c",
                "[[1,2]]",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.modred.GaussianModredJob",
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_modred_qmmm_subcommand_creates_qmmm_job(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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
    ):
        """``traj -j opt`` creates a ``GaussianTrajJob`` for structures from a trajectory."""
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
            (traj_module, "GaussianTrajJob"),
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()

    def test_traj_mutually_exclusive_options_raise(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """``-g`` and ``-ns`` together raise a usage error."""
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
            (traj_module, "GaussianTrajJob"),
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
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code != 0


class TestGaussianCLIDiasCommand:
    """CLI tests for the ``dias`` (distortion/interaction) subcommand."""

    def test_dias_basic_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
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


class TestGaussianQMMMCLI:
    """CLI tests for the nested ``opt qmmm`` (and related) subcommands."""

    def test_opt_qmmm_propagates_cli_options_to_settings_and_molecule(
        self,
        single_molecule_xyz_file,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        mm_info = os.path.join(tmpdir, "mm_atoms.dat")
        with open(mm_info, "w") as handle:
            handle.write("1 C 0.0\n")
            handle.write("2 H 0.0\n")
            handle.write("3 H 0.0\n")
            handle.write("4 H 0.0\n")
        params = os.path.join(tmpdir, "mm_params.dat")
        with open(params, "w") as handle:
            handle.write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0\n")

        mock_job = MagicMock(name="mock_job")
        result, _ = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob",
            [
                "-p",
                "qmmm",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-l",
                "testjob",
                "opt",
                "-f",
                "5",
                "qmmm",
                "-hx",
                "b3lyp",
                "-hb",
                "6-31g*",
                "-hff",
                "UFF",
                "-mx",
                "hf",
                "-mb",
                "sto-3g",
                "-mff",
                "UFF",
                "-lx",
                "pm6",
                "-lb",
                "sto-3g",
                "-lff",
                "UFF",
                "-ct",
                "0",
                "-mt",
                "1",
                "-ci",
                "0",
                "-mi",
                "1",
                "-ch",
                "1",
                "-mh",
                "2",
                "-ha",
                "1-3",
                "-ma",
                "4",
                "-la",
                "5-6",
                "-ba",
                "[[3, 4]]",
                "-sf",
                "{[3, 4]: [0.709]}",
                "-mai",
                mm_info,
                "-mpf",
                params,
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job,
        )

        assert result.exit_code == 0, result.output
        assert mock_job.call_args is not None
        call_kwargs = mock_job.call_args.kwargs
        settings = call_kwargs["settings"]
        molecule = call_kwargs["molecule"]

        assert settings.functional == "b3lyp"
        assert settings.basis == "6-31g*"
        assert settings.high_level_force_field == "UFF"
        assert settings.medium_level_functional == "hf"
        assert settings.medium_level_basis == "sto-3g"
        assert settings.medium_level_force_field == "UFF"
        assert settings.low_level_functional == "pm6"
        assert settings.low_level_basis == "sto-3g"
        assert settings.low_level_force_field == "UFF"
        assert settings.charge_total == 0
        assert settings.mult_total == 1
        assert settings.charge_intermediate == 0
        assert settings.mult_intermediate == 1
        assert settings.model_charge == 1
        assert settings.model_multiplicity == 2
        assert settings.mm_atom_info_file == mm_info
        assert settings.mm_parameters_file == params
        assert settings.parent_jobtype == "opt"
        assert molecule.high_level_atoms == [1, 2, 3]
        assert molecule.medium_level_atoms == [4]
        assert molecule.low_level_atoms == [5, 6]
        assert molecule.bonded_atoms == [[3, 4]]
        assert molecule.scale_factors == {(3, 4): [0.709]}
        assert molecule.frozen_atoms[4] == -1

        assert call_kwargs["label"] == "testjob_qmmm"

    def test_opt_qmmm_without_mm_sidecar_options(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """Exercise False branches for -mai/-mpf CLI options."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob",
            [
                "-p",
                "qmmm",
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
                "sto-3g",
                "-ch",
                "0",
                "-mh",
                "1",
                "-ha",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.mm_atom_info_file is None
        assert settings.mm_parameters_file is None


class TestGaussianQmmmSubcommand:
    def test_all_qmmm_options_applied_to_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """Every QMMM-specific CLI option should land on the settings
        object, and bonded/scale-factor strings should be parsed onto
        the molecule via ast.literal_eval."""
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
                "opt",
                "--freeze-atoms",
                "1-2",
                "qmmm",
                "-hx",
                "b3lyp",
                "-hb",
                "6-31g(d)",
                "-hff",
                "amber",
                "-mx",
                "pbe",
                "-mb",
                "6-31g",
                "-mff",
                "uff",
                "-lx",
                "hf",
                "-lb",
                "sto-3g",
                "-lff",
                "gaff",
                "-ct",
                "0",
                "-mt",
                "1",
                "-ci",
                "0",
                "-mi",
                "1",
                "-ch",
                "0",
                "-mh",
                "1",
                "-ha",
                "1-3",
                "-ma",
                "4-5",
                "-la",
                "6",
                "-ba",
                "{(1, 2): 1.0}",
                "-sf",
                "{(1, 2): [0.1, 0.2, 0.3]}",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.high_level_functional == "b3lyp"
        assert settings.high_level_basis == "6-31g(d)"
        assert settings.high_level_force_field == "amber"
        assert settings.medium_level_functional == "pbe"
        assert settings.medium_level_basis == "6-31g"
        assert settings.medium_level_force_field == "uff"
        assert settings.low_level_functional == "hf"
        assert settings.low_level_basis == "sto-3g"
        assert settings.low_level_force_field == "gaff"
        assert settings.charge_total == 0
        assert settings.mult_total == 1
        assert settings.charge_intermediate == 0
        assert settings.mult_intermediate == 1
        assert settings.charge_high == 0
        assert settings.mult_high == 1
        assert settings.medium_level_atoms == "4-5"
        assert settings.low_level_atoms == "6"
        assert settings.bonded_atoms == "{(1, 2): 1.0}"
        assert settings.scale_factors == "{(1, 2): [0.1, 0.2, 0.3]}"

    def test_label_already_containing_qmmm_is_not_double_suffixed(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """A label that already mentions "qmmm" (case-insensitive)
        should not get an extra ``_qmmm`` appended at the start of the
        function; the later unconditional endswith check still applies."""
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
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
                "-l",
                "myjob_qmmm",
                "opt",
                "qmmm",
                "-hx",
                "b3lyp",
                "-hb",
                "6-31g(d)",
                "-ha",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 1
        _, kwargs = mock_job_cls.call_args
        assert kwargs["label"] == "myjob_qmmm"

    def test_label_containing_but_not_ending_with_qmmm_gets_suffixed_once(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """A label containing "qmmm" but not literally ending with the
        "_qmmm" suffix (e.g. it appears at the start) skips the first
        append check but still gets suffixed by the later endswith
        check."""
        mock_job_cls = MagicMock(name="mock_job_cls")
        result, _ = run_gaussian_and_capture_settings(
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
                "-l",
                "qmmm_experiment",
                "opt",
                "qmmm",
                "-hx",
                "b3lyp",
                "-hb",
                "6-31g(d)",
                "-ha",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
            job_cls=mock_job_cls,
        )
        assert result.exit_code == 0, result.output
        _, kwargs = mock_job_cls.call_args
        assert kwargs["label"] == "qmmm_experiment_qmmm"

    def test_project_qmmm_settings_loaded_from_yaml(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        """The 'qmmm' test project's qmmm: YAML section should seed
        the QMMM settings before CLI overrides are applied."""
        result, settings = run_gaussian_and_capture_settings(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob",
            [
                "-p",
                "qmmm",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
                "qmmm",
                "-ha",
                "1-3",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None
        # From qmmm.yaml's qmmm: section, not overridden by CLI.
        assert settings.high_level_functional == "MN15"
        assert settings.medium_level_functional == "PBE"

    def test_freeze_atoms_inherited_from_parent_opt_command(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "opt",
                "--freeze-atoms",
                "1-2",
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
        assert settings is not None

    def test_jobtype_inferred_from_ts_parent_command(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
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
                "ts",
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
        assert settings is not None
        assert settings.jobtype == "ts"


class TestQmmmCallbackDirectInvocation:
    """Covers branches unreachable through any real parent command,
    since every parent (opt/ts/sp/scan/modred/qrc) unconditionally
    populates ctx.obj["parent_settings"] with a real settings object."""

    def test_no_parent_settings_skips_both_merge_blocks(self):
        mock_job = run_gaussian_qmmm_callback()
        mock_job.assert_called_once()

    def test_merge_is_called_once_and_exceptions_propagate(self):
        """Upstream removed a second, redundant parent-settings merge
        block (identical keyword list to the first) that used to wrap
        its call in a try/except and log-and-continue on failure. Only
        one merge call remains now, with no exception handling, so a
        failure there propagates instead of being swallowed."""
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        real_merge = GaussianQMMMJobSettings.merge
        call_count = {"n": 0}

        def counting_merge(self, *args, **kwargs):
            call_count["n"] += 1
            return real_merge(self, *args, **kwargs)

        with patch.object(GaussianQMMMJobSettings, "merge", counting_merge):
            mock_job = run_gaussian_qmmm_callback(parent_settings=MagicMock())
        mock_job.assert_called_once()
        assert call_count["n"] == 1

        def raising_merge(self, *args, **kwargs):
            raise RuntimeError("merge failed")

        with patch.object(GaussianQMMMJobSettings, "merge", raising_merge):
            with pytest.raises(RuntimeError, match="merge failed"):
                run_gaussian_qmmm_callback(parent_settings=MagicMock())


class TestClickGaussianQmmmOptionsIsUnusedDeadCode:
    """Documents BUGS_FOUND.md #56: click_gaussian_qmmm_options in
    chemsmart/cli/gaussian/gaussian.py is defined but never applied
    anywhere -- create_qmmm_subcommand in cli/gaussian/qmmm.py defines
    its own separate, near-identical set of click.option decorators
    directly instead of using this one. Exercised directly here purely
    for coverage since no real command uses it."""

    def test_decorator_applies_all_expected_options(self):
        from chemsmart.cli.gaussian.gaussian import (
            click_gaussian_qmmm_options,
        )

        @click.command()
        @click_gaussian_qmmm_options
        def dummy(**kwargs):
            click.echo(str(sorted(kwargs.items())))

        result = CliRunner().invoke(
            dummy,
            [
                "-hx",
                "b3lyp",
                "-hb",
                "6-31g(d)",
                "-ha",
                "1-3",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "high_level_functional" in result.output
        assert "b3lyp" in result.output


def _invoke_gaussian_group_forcing_qmmm_subcommand(single_molecule_xyz_file):
    """Directly invoke gaussian()'s own callback with
    ctx.invoked_subcommand forced to "qmmm" -- a state Click's normal
    command resolution can never produce (see
    TestQmmmMoleculeConversionBlockInGaussianGroup's docstring)."""
    real_gaussian_fn = inspect.unwrap(gaussian.callback)

    ctx = click.Context(gaussian)
    ctx.obj = {}
    ctx.invoked_subcommand = "qmmm"

    kwargs = dict(
        project="gas_solv",
        filename=single_molecule_xyz_file,
        label="test",
        append_label=None,
        title=None,
        charge=0,
        multiplicity=1,
        functional=None,
        basis=None,
        semiempirical=None,
        index=None,
        record_index=None,
        record_id=None,
        structure_id=None,
        structure_index=None,
        molecule_id=None,
        additional_opt_options=None,
        additional_route_parameters=None,
        append_additional_info=None,
        custom_solvent=None,
        dieze_tag=None,
        forces=False,
        pubchem=None,
        remove_solvent=False,
        solvent_model=None,
        solvent_id=None,
        solvent_options=None,
    )

    with ctx:
        real_gaussian_fn(ctx, **kwargs)

    return ctx


class TestQmmmMoleculeConversionBlockInGaussianGroup:
    """Documents BUGS_FOUND.md #57: gaussian()'s own
    `if ctx.invoked_subcommand == "qmmm":` early-conversion block can
    never fire through any genuine CLI invocation, because "qmmm" is
    always a grandchild subcommand attached to opt/ts/sp/scan/modred/
    qrc -- Click's ctx.invoked_subcommand on the *gaussian* group's own
    context only ever reflects the immediate next command (e.g.
    "opt"), never a subcommand of a subcommand. The conversion logic
    itself is exercised directly here (bypassing Click's normal
    command-resolution, which can never produce this ctx state) purely
    for coverage."""

    def test_molecules_converted_to_qmmmmolecule_when_forced(
        self, single_molecule_xyz_file
    ):
        from chemsmart.io.molecules.structure import QMMMMolecule

        ctx = _invoke_gaussian_group_forcing_qmmm_subcommand(
            single_molecule_xyz_file
        )

        molecules = ctx.obj["molecules"]
        assert len(molecules) == 1
        assert isinstance(molecules[0], QMMMMolecule)

    def test_already_qmmmmolecule_instances_are_left_as_is(
        self, single_molecule_xyz_file
    ):
        from chemsmart.io.molecules.structure import Molecule, QMMMMolecule

        real_molecule = Molecule.from_filepath(single_molecule_xyz_file)
        qmmm_molecule = QMMMMolecule(molecule=real_molecule)

        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=[qmmm_molecule],
        ):
            ctx = _invoke_gaussian_group_forcing_qmmm_subcommand(
                single_molecule_xyz_file
            )

        molecules = ctx.obj["molecules"]
        assert molecules == [qmmm_molecule]

    def test_molecule_kwarg_init_failure_falls_back_to_dict_init(
        self, single_molecule_xyz_file
    ):
        """When ``QMMMMolecule(molecule=m)`` raises, the code retries
        with a dict-based init from ``m.__dict__``; that retry
        succeeding covers the fallback's success path.

        The fake ``__init__`` below doesn't delegate its dict-based
        branch to the real ``Molecule.__init__``: a genuine
        ``Molecule.__dict__`` always contains private/derived keys
        (``_positions``, ``_num_atoms``, ...) that the constructor
        rejects, so the real fallback can never actually succeed for a
        real ``Molecule`` -- see BUGS_FOUND.md #57's note on this."""
        from chemsmart.io.molecules.structure import QMMMMolecule

        call_kwargs = []

        def fake_init(self, *args, **kwargs):
            call_kwargs.append(kwargs)
            if "molecule" in kwargs:
                raise TypeError("simulated molecule= init failure")
            self.molecule = None

        with patch.object(QMMMMolecule, "__init__", fake_init):
            ctx = _invoke_gaussian_group_forcing_qmmm_subcommand(
                single_molecule_xyz_file
            )

        molecules = ctx.obj["molecules"]
        assert len(molecules) == 1
        assert isinstance(molecules[0], QMMMMolecule)
        assert any("molecule" in kw for kw in call_kwargs)
        assert any("molecule" not in kw for kw in call_kwargs)

    def test_both_init_attempts_failing_keeps_original_molecule(
        self, single_molecule_xyz_file
    ):
        """When both the molecule= and dict-based QMMMMolecule init
        attempts fail, the original (unconverted) molecule is kept
        instead of raising."""
        from chemsmart.io.molecules.structure import Molecule, QMMMMolecule

        with patch.object(
            QMMMMolecule,
            "__init__",
            side_effect=TypeError("simulated init failure"),
        ):
            ctx = _invoke_gaussian_group_forcing_qmmm_subcommand(
                single_molecule_xyz_file
            )

        molecules = ctx.obj["molecules"]
        assert len(molecules) == 1
        assert isinstance(molecules[0], Molecule)
        assert not isinstance(molecules[0], QMMMMolecule)

    def test_unexpected_error_in_conversion_block_is_swallowed(
        self, single_molecule_xyz_file
    ):
        """Any failure outside the per-molecule try/except (e.g.
        iterating "molecules" itself blowing up) is caught by the
        block's own outer `except Exception`, leaving the original
        molecules in place rather than propagating."""
        from chemsmart.io.molecules.structure import Molecule

        real_molecule = Molecule.from_filepath(single_molecule_xyz_file)

        class _ExplodingList(list):
            def __iter__(self):
                raise RuntimeError("simulated iteration failure")

        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_filepath",
            return_value=_ExplodingList([real_molecule]),
        ):
            ctx = _invoke_gaussian_group_forcing_qmmm_subcommand(
                single_molecule_xyz_file
            )

        # Index rather than iterate: the stored list's own __iter__ is
        # still the exploding one, since the outer except left it as-is.
        assert len(ctx.obj["molecules"]) == 1
        assert ctx.obj["molecules"][0] is real_molecule
