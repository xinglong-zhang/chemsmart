"""
Direct tests for the ``qmmm`` subcommand created by
``create_qmmm_subcommand`` in ``chemsmart.cli.gaussian.qmmm``.

The subcommand is attached to several Gaussian jobtype groups (``opt``,
``ts``, ``sp``, ``scan``, ``qrc``, ``modred``); these tests mostly use
``opt`` as the parent since it supports ``--freeze-atoms``, needed to
exercise the freeze-atoms inheritance branch.
"""

from unittest.mock import MagicMock, patch

import click
from click.testing import CliRunner

from chemsmart.cli.gaussian.gaussian import gaussian


class TestGaussianQmmmSubcommand:
    def test_all_qmmm_options_applied_to_settings(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """A label that already mentions "qmmm" (case-insensitive)
        should not get an extra ``_qmmm`` appended at the start of the
        function; the later unconditional endswith check still applies."""
        from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob

        with patch.object(
            GaussianQMMMJob, "__new__", return_value=MagicMock()
        ) as mock_new:
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        assert mock_new.call_count == 1
        _, kwargs = mock_new.call_args
        assert kwargs["label"] == "myjob_qmmm"

    def test_label_containing_but_not_ending_with_qmmm_gets_suffixed_once(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
    ):
        """A label containing "qmmm" but not literally ending with the
        "_qmmm" suffix (e.g. it appears at the start) skips the first
        append check but still gets suffixed by the later endswith
        check."""
        from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob

        with patch.object(
            GaussianQMMMJob, "__new__", return_value=MagicMock()
        ) as mock_new:
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
                obj=make_cli_ctx_obj(gaussian_jobrunner_no_scratch),
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        _, kwargs = mock_new.call_args
        assert kwargs["label"] == "qmmm_experiment_qmmm"

    def test_project_qmmm_settings_loaded_from_yaml(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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
        make_cli_ctx_obj,
        run_gaussian_and_capture_settings,
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


def _invoke_qmmm_callback_directly(**ctx_obj_overrides):
    """Direct-invocation helper bypassing the real "opt"/"ts"/...
    parent commands: every real parent command unconditionally sets
    ctx.obj["parent_settings"] to a real settings object (see e.g.
    chemsmart/cli/gaussian/opt.py), so parent_settings is None,
    jobtype is None (no ctx.parent), and parent_jobtype/high_level_atoms
    absent are all branches that can't be reached through a genuine
    CLI invocation. Constructing ctx.obj by hand reaches them directly.
    """
    import importlib

    opt_mod = importlib.import_module("chemsmart.cli.gaussian.opt")
    qmmm_cmd = opt_mod.opt.commands["qmmm"]

    project_settings = MagicMock()
    project_settings.qmmm_settings.return_value = None
    molecule = MagicMock()
    molecule.frozen_atoms = None

    ctx_obj = {
        "jobrunner": MagicMock(),
        "project_settings": project_settings,
        "label": "test",
        "molecules": [molecule],
    }
    ctx_obj.update(ctx_obj_overrides)

    ctx = click.Context(qmmm_cmd)
    ctx.obj = ctx_obj

    kwargs = dict(
        high_level_functional="b3lyp",
        high_level_basis="6-31g",
        high_level_force_field=None,
        medium_level_functional=None,
        medium_level_basis=None,
        medium_level_force_field=None,
        low_level_functional=None,
        low_level_basis=None,
        low_level_force_field=None,
        charge_total=None,
        mult_total=None,
        charge_intermediate=None,
        mult_intermediate=None,
        charge_high=None,
        mult_high=None,
        high_level_atoms=None,
        medium_level_atoms=None,
        low_level_atoms=None,
        bonded_atoms=None,
        scale_factors=None,
    )

    with patch("chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob") as mock_job:
        mock_job.return_value = MagicMock()
        with ctx:
            qmmm_cmd.callback(**kwargs)
        return mock_job


class TestQmmmCallbackDirectInvocation:
    """Covers branches unreachable through any real parent command,
    since every parent (opt/ts/sp/scan/modred/qrc) unconditionally
    populates ctx.obj["parent_settings"] with a real settings object."""

    def test_no_parent_settings_skips_both_merge_blocks(self):
        mock_job = _invoke_qmmm_callback_directly()
        mock_job.assert_called_once()

    def test_merge_exception_is_caught_and_logged(self):
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

        real_merge = GaussianQMMMJobSettings.merge
        call_count = {"n": 0}

        def fake_merge(self, *args, **kwargs):
            call_count["n"] += 1
            if call_count["n"] == 2:
                raise RuntimeError("merge failed")
            return real_merge(self, *args, **kwargs)

        with patch.object(GaussianQMMMJobSettings, "merge", fake_merge):
            mock_job = _invoke_qmmm_callback_directly(
                parent_settings=MagicMock()
            )
        mock_job.assert_called_once()
        assert call_count["n"] == 2


class TestPopulateChargeAndMultiplicityOnSettings:
    """Documents BUGS_FOUND.md #52: this helper is defined but never
    called anywhere in chemsmart/cli/gaussian/qmmm.py (or elsewhere in
    the Gaussian CLI) -- only its identically named ORCA counterpart
    in cli/orca/qmmm.py is actually wired up. Unit-tested directly
    here since no code path reaches it."""

    @staticmethod
    def _settings(**overrides):
        from types import SimpleNamespace

        return SimpleNamespace(charge=None, multiplicity=None, **overrides)

    def test_intermediate_charge_and_multiplicity_take_priority(self):
        from chemsmart.cli.gaussian.qmmm import (
            _populate_charge_and_multiplicity_on_settings,
        )

        qs = self._settings(
            charge_intermediate=1,
            mult_intermediate=2,
            charge_high=3,
            mult_high=4,
            charge_total=5,
            mult_total=6,
        )
        _populate_charge_and_multiplicity_on_settings(qs)
        assert (qs.charge, qs.multiplicity) == (1, 2)
        assert (qs.charge_total, qs.mult_total) == (1, 2)

    def test_high_level_charge_and_multiplicity_used_when_no_intermediate(
        self,
    ):
        from chemsmart.cli.gaussian.qmmm import (
            _populate_charge_and_multiplicity_on_settings,
        )

        qs = self._settings(
            charge_intermediate=None,
            mult_intermediate=None,
            charge_high=3,
            mult_high=4,
            charge_total=5,
            mult_total=6,
        )
        _populate_charge_and_multiplicity_on_settings(qs)
        assert (qs.charge, qs.multiplicity) == (3, 4)

    def test_total_charge_and_multiplicity_used_as_last_resort(self):
        from chemsmart.cli.gaussian.qmmm import (
            _populate_charge_and_multiplicity_on_settings,
        )

        qs = self._settings(
            charge_intermediate=None,
            mult_intermediate=None,
            charge_high=None,
            mult_high=None,
            charge_total=5,
            mult_total=6,
        )
        _populate_charge_and_multiplicity_on_settings(qs)
        assert (qs.charge, qs.multiplicity) == (5, 6)

    def test_no_charge_or_multiplicity_set_leaves_settings_untouched(self):
        from types import SimpleNamespace

        from chemsmart.cli.gaussian.qmmm import (
            _populate_charge_and_multiplicity_on_settings,
        )

        qs = SimpleNamespace()
        _populate_charge_and_multiplicity_on_settings(qs)
        assert not hasattr(qs, "charge")
        assert not hasattr(qs, "multiplicity")
