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

import os
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
    """CLI tests for PubChem / database label paths without an explicit ``-l``."""

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


class TestGaussianCLIOptCommand:
    """CLI tests for the ``opt`` subcommand's freeze-atoms and
    multi-molecule (index-selected) job-creation branches."""

    def test_single_molecule_freeze_atoms(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
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
                    "opt",
                    "--freeze-atoms",
                    "1-2",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        molecule = mock_job_cls.call_args[1]["molecule"]
        assert molecule.frozen_atoms is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
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

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.tddft.GaussianTDDFTJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
                [
                    "-p",
                    "solv",
                    "-f",
                    multiple_molecules_xyz_file,
                    "-i",
                    "1-2",
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "td",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]


class TestGaussianCLIGroupValidation:
    """Validation and settings-merge branches on the ``gaussian`` group
    callback itself (not delegated to any subcommand)."""

    def test_molecule_id_not_supported_raises(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-f",
                single_molecule_xyz_file,
                "--mid",
                "abc",
                "opt",
            ],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
        )
        assert result.exit_code != 0
        assert "not supported for Gaussian job submission" in result.output

    def test_index_and_structure_index_mutually_exclusive(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-f",
                single_molecule_xyz_file,
                "-i",
                "1",
                "--si",
                "1",
                "opt",
            ],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_structure_index_used_as_index(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        self, database_chemsmart_file, gaussian_jobrunner_no_scratch
    ):
        """No selector and no --sid given for a chemsmart db input."""
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            ["-f", database_chemsmart_file, "opt"],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_rejects_multiple_selectors(
        self, database_chemsmart_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-f",
                database_chemsmart_file,
                "--ri",
                "1",
                "--sid",
                "abc",
                "opt",
            ],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
        )
        assert result.exit_code != 0
        assert "select exactly one of" in result.output

    def test_chemsmart_db_index_requires_record_selector(
        self, database_chemsmart_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-f",
                database_chemsmart_file,
                "--sid",
                "abc",
                "-i",
                "1",
                "opt",
            ],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
        )
        assert result.exit_code != 0
        assert "can only be used together with" in result.output

    def test_xtb_output_inherits_charge_and_multiplicity(
        self,
        xtb_water_outfolder,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """A ``.out`` file detected as xTB should fall back to default
        Gaussian settings but inherit charge/multiplicity from xTB."""
        import os

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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """If the xTB-derived molecule has no charge/multiplicity, the
        default Gaussian settings' own (unset) values are kept."""
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """A .db file that isn't a chemsmart database (e.g. a plain ASE
        db) falls back to default Gaussian settings entirely."""
        from unittest.mock import MagicMock, patch

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
        self, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            ["-p", "gas_solv", "opt"],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "has not been specified" in str(result.exception)

    def test_filename_and_pubchem_both_given_raises(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "--pubchem",
                "222",
                "opt",
            ],
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "have been specified" in str(result.exception)

    def test_label_and_append_label_mutually_exclusive_raises(
        self, single_molecule_xyz_file, gaussian_jobrunner_no_scratch
    ):
        runner = CliRunner()
        result = runner.invoke(
            gaussian,
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
            obj={"jobrunner": gaussian_jobrunner_no_scratch},
            catch_exceptions=True,
        )
        assert result.exit_code != 0
        assert isinstance(result.exception, ValueError)
        assert "not both" in str(result.exception)

    def test_append_label_suffixes_filename_derived_label(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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

    def test_pubchem_only_without_label_crashes(
        self,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """BUG (see BUGS_FOUND.md #25): the intended "output" label
        fallback for a filename-less (PubChem-only) job with no -l/-a
        is unreachable -- `os.path.basename(filename)` at gaussian.py's
        default-label branch is called unconditionally on `filename`
        (None here) before the `if filename:` guard that would skip it,
        so this currently crashes with a TypeError instead of falling
        back to "output"."""
        from unittest.mock import MagicMock, patch

        pubchem_molecule = MagicMock(name="pubchem_molecule")
        with patch(
            "chemsmart.io.molecules.structure.Molecule.from_pubchem",
            return_value=[pubchem_molecule],
        ):
            with pytest.raises(TypeError):
                run_gaussian_and_capture_settings(
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
                )

    def test_default_label_with_filename_ignores_db_id_suffix(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """--rid alone (no -a/-l) takes the default-label elif branch
        for the record_id suffix (later overwritten by the basename
        recompute -- see BUGS_FOUND.md #25 -- but the branch itself
        must still execute)."""
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
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_group_level_functional_and_basis_keywords_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        reach a couple of defensive branches in the default-label
        fallback chain that real CLI invocations can't trigger: an
        empty (not None) filename sidesteps the unconditional
        `os.path.basename(filename)` crash documented in
        BUGS_FOUND.md #25, letting execution reach the `if filename:
        ... else: label = "output"` branch below it."""
        import inspect

        import click

        from chemsmart.cli.gaussian.gaussian import gaussian

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
        arm directly, since a real CLI run can't reach it (filename is
        either a real path or None, and None crashes earlier per
        BUGS_FOUND.md #25)."""
        ctx = self._direct_invoke_gaussian_group(_invoked_subcommand=None)
        assert ctx.obj["label"] == "output"

    def test_invoked_subcommand_suffixes_default_label_when_present(self):
        """Covers the `if ctx.invoked_subcommand:` True arm of the
        same fallback chain: the invoked subcommand name gets appended
        to the "output" default label."""
        ctx = self._direct_invoke_gaussian_group(_invoked_subcommand="opt")
        assert ctx.obj["label"] == "output_opt"


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

    def test_sp_subcommand_level_solvent_options_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
                "sp",
                "-so",
                "iterative",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.additional_solvent_options == "iterative"

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_qmmm_child_subcommand_skips_direct_sp_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.singlepoint.GaussianSinglePointJob"
        ) as mock_sp_job_cls:
            result = CliRunner().invoke(
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "6-31g(d)",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_sp_job_cls.assert_not_called()

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

    def test_single_molecule_freeze_atoms(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch("chemsmart.jobs.gaussian.ts.GaussianTSJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
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
                    "--freeze-atoms",
                    "1-2",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_job_cls.assert_called_once()
        molecule = mock_job_cls.call_args[1]["molecule"]
        assert molecule.frozen_atoms is not None

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch("chemsmart.jobs.gaussian.ts.GaussianTSJob") as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
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
                    "ts",
                    "--freeze-atoms",
                    "1",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]
        for c in mock_job_cls.call_args_list:
            assert c.kwargs["molecule"].frozen_atoms is not None

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

    def test_irc_predictor_and_recorrect_options(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
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
                "-pt",
                "HPC",
                "-rc",
                "Always",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.predictor == "HPC"
        assert settings.recorrect == "Always"

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

    def test_explicit_jobtype_modred_used(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """``-j modred`` explicitly selects modred settings resolution."""
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
                "-j",
                "modred",
                "-c",
                "[[1,2]]",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings is not None

    def test_constrained_coordinates_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
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
                "-cc",
                "[[3,4]]",
            ],
            make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
        )
        assert result.exit_code == 0, result.output
        assert settings.modred["constrained_coordinates"] == [[3, 4]]

    def test_multiple_molecules_with_indices_creates_one_job_each(
        self,
        multiple_molecules_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
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
                    "[[1,2]]",
                    "-s",
                    "0.1",
                    "-n",
                    "10",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

    def test_qmmm_child_subcommand_skips_direct_scan_job_creation(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.scan.GaussianScanJob"
        ) as mock_scan_job_cls:
            result = CliRunner().invoke(
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
                    "qmmm",
                    "-hx",
                    "b3lyp",
                    "-hb",
                    "6-31g(d)",
                ],
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        mock_scan_job_cls.assert_not_called()

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
    """CLI tests for the ``qrc`` (Quick Reaction Coordinate) subcommand."""

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

    def test_link_subcommand_level_solvent_options_applied(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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

    def test_explicit_jobtype_used(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
    ):
        with patch(
            "chemsmart.jobs.gaussian.modred.GaussianModredJob"
        ) as mock_job_cls:
            mock_job_cls.return_value = MagicMock()
            result = CliRunner().invoke(
                gaussian,
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_job_cls.call_count == 2
        labels = [c.kwargs["label"] for c in mock_job_cls.call_args_list]
        assert labels[0] != labels[1]

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


class TestGaussianQMMMCLI:
    """CLI tests for the nested ``opt qmmm`` (and related) subcommands."""

    def test_opt_qmmm_propagates_cli_options_to_settings_and_molecule(
        self,
        single_molecule_xyz_file,
        tmpdir,
        gaussian_jobrunner_no_scratch,
    ):
        from click.testing import CliRunner

        from chemsmart.cli.gaussian.gaussian import gaussian

        mm_info = os.path.join(tmpdir, "mm_atoms.dat")
        with open(mm_info, "w") as handle:
            handle.write("1 C 0.0\n")
            handle.write("2 H 0.0\n")
            handle.write("3 H 0.0\n")
            handle.write("4 H 0.0\n")
        params = os.path.join(tmpdir, "mm_params.dat")
        with open(params, "w") as handle:
            handle.write("NonBon 3 1 0 0 0.0 0.0 0.5 0.0 0.0 0.0\n")

        runner = CliRunner()
        with patch("chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob") as mock_job:
            mock_job.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj={"jobrunner": gaussian_jobrunner_no_scratch},
                catch_exceptions=False,
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
        from click.testing import CliRunner

        from chemsmart.cli.gaussian.gaussian import gaussian

        runner = CliRunner()
        with patch("chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob") as mock_job:
            mock_job.return_value = MagicMock()
            result = runner.invoke(
                gaussian,
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
                obj={"jobrunner": gaussian_jobrunner_no_scratch},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        settings = mock_job.call_args.kwargs["settings"]
        assert settings.mm_atom_info_file is None
        assert settings.mm_parameters_file is None
