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
        with patch(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob"
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
        assert mock_job_cls.call_count == 1
        _, kwargs = mock_job_cls.call_args
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
        with patch(
            "chemsmart.jobs.gaussian.qmmm.GaussianQMMMJob"
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
        _, kwargs = mock_job_cls.call_args
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


class TestClickGaussianQmmmOptionsIsUnusedDeadCode:
    """Documents BUGS_FOUND.md #56: click_gaussian_qmmm_options in
    chemsmart/cli/gaussian/gaussian.py is defined but never applied
    anywhere -- create_qmmm_subcommand in cli/gaussian/qmmm.py defines
    its own separate, near-identical set of click.option decorators
    directly instead of using this one. Exercised directly here purely
    for coverage since no real command uses it."""

    def test_decorator_applies_all_expected_options(self):
        import click
        from click.testing import CliRunner

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
    import inspect

    import click

    from chemsmart.cli.gaussian.gaussian import gaussian

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
