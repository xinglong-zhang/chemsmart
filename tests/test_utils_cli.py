"""Tests for chemsmart.utils.cli."""

from unittest.mock import MagicMock

import click
import pytest
from click.testing import CliRunner

from chemsmart.utils.cli import (
    CtxObjArguments,
    MyCommand,
    MyGroup,
    check_scan_coordinates_gaussian,
    check_scan_coordinates_orca,
    check_scan_parameters_consistency_gaussian,
    check_scan_parameters_consistency_orca,
    create_sp_label,
    get_setting_from_jobtype_for_gaussian,
    get_setting_from_jobtype_for_orca,
    update_irc_label,
)


class TestMyGroupAndCommand:
    """Tests that MyGroup/MyCommand record subcommand metadata on ctx.obj."""

    def test_group_and_command_record_subcommand_chain(self):
        @click.group(cls=MyGroup)
        @click.option("--verbose", is_flag=True)
        @click.pass_context
        def cli(ctx, verbose):
            pass

        @cli.command(cls=MyCommand)
        @click.option("--name", default="x")
        @click.pass_context
        def sub(ctx, name):
            pass

        runner = CliRunner()
        result = runner.invoke(
            cli, ["--verbose", "sub", "--name", "abc"], obj={}
        )
        assert result.exit_code == 0

    def test_ctx_obj_none_falls_back_to_dict(self, mocker):
        """Covers the ensure_object exception fallback path."""
        from chemsmart.utils.cli import _add_subcommand_info_to_ctx

        ctx = MagicMock()
        ctx.ensure_object.side_effect = Exception("no ensure_object")
        ctx.obj = None
        ctx.info_name = "mycmd"
        ctx.command.params = []
        ctx.params = {}
        ctx.parent = None

        _add_subcommand_info_to_ctx(ctx)
        assert ctx.obj["subcommand"][0]["name"] == "mycmd"


class TestCtxObjArgumentsValue:
    """Tests for CtxObjArguments._value."""

    def test_true_returns_empty_string(self):
        assert CtxObjArguments._value(True) == ""

    def test_false_returns_empty_string(self):
        assert CtxObjArguments._value(False) == ""

    def test_other_values_stringified(self):
        assert CtxObjArguments._value(42) == "42"
        assert CtxObjArguments._value("abc") == "abc"


class TestCtxObjArgumentsKeyword:
    """Tests for CtxObjArguments.argument_keyword."""

    def test_long_option_name(self):
        keyword = CtxObjArguments.argument_keyword(
            "my_arg", "value", False, []
        )
        assert keyword == "--my-arg"

    def test_short_option_name(self):
        keyword = CtxObjArguments.argument_keyword("n", "5", False, [])
        assert keyword == "-n"

    def test_false_simple_flag_returns_empty(self):
        keyword = CtxObjArguments.argument_keyword("verbose", False, True, [])
        assert keyword == ""

    def test_false_value_with_one_secondary_opt(self):
        keyword = CtxObjArguments.argument_keyword(
            "flag", False, False, ["--no-flag"]
        )
        assert keyword == "--no-flag"

    def test_false_value_with_multiple_secondary_opts(self):
        keyword = CtxObjArguments.argument_keyword(
            "flag", False, False, ["--flag/--no-flag", "--other"]
        )
        assert keyword == "--other"


class TestCtxObjArgumentsEntryPoint:
    """Tests for CtxObjArguments._entry_point."""

    def test_finds_entry_point_with_no_parent(self):
        commands = [
            {"name": "root", "parent": None, "kwargs": {}},
            {"name": "child", "parent": "root", "kwargs": {}},
        ]
        obj = CtxObjArguments(commands)
        assert obj._entry_point()["name"] == "root"

    def test_finds_explicit_entry_point(self):
        commands = [
            {"name": "root", "parent": None, "kwargs": {}},
            {"name": "child", "parent": "root", "kwargs": {}},
        ]
        obj = CtxObjArguments(commands, entry_point="child")
        assert obj._entry_point()["name"] == "child"

    def test_raises_value_error_for_missing_explicit_entry_point(self):
        commands = [{"name": "root", "parent": None, "kwargs": {}}]
        obj = CtxObjArguments(commands, entry_point="missing")
        with pytest.raises(ValueError, match="could not be found"):
            obj._entry_point()

    def test_raises_assertion_error_when_no_entry_point_found(self):
        commands = [{"name": "child", "parent": "root", "kwargs": {}}]
        obj = CtxObjArguments(commands)
        with pytest.raises(AssertionError, match="Could not find entry"):
            obj._entry_point()


class TestCtxObjArgumentsReconstruct:
    """Tests for command-line reconstruction."""

    @staticmethod
    def _kw(
        value,
        nargs=1,
        is_multiple=False,
        type_name="string",
        is_flag=False,
        secondary_opts=None,
    ):
        return {
            "value": value,
            "nargs": nargs,
            "is_multiple": is_multiple,
            "type": MagicMock(name=type_name),
            "is_flag": is_flag,
            "secondary_opts": secondary_opts or [],
        }

    def test_reconstruct_command_skips_none_values(self):
        commands = [
            {
                "name": "root",
                "parent": None,
                "kwargs": {"skip_me": self._kw(None)},
            }
        ]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root"]

    def test_reconstruct_command_skips_empty_tuple(self):
        commands = [
            {
                "name": "root",
                "parent": None,
                "kwargs": {"skip_me": self._kw(())},
            }
        ]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root"]

    def test_reconstruct_command_simple_value(self):
        kw = self._kw("val")
        kw["type"].name = "string"
        commands = [{"name": "root", "parent": None, "kwargs": {"name": kw}}]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "--name", "val"]

    def test_reconstruct_command_literal_eval_type(self):
        kw = self._kw("[1, 2]")
        kw["type"].name = "literal_eval"
        commands = [{"name": "root", "parent": None, "kwargs": {"coords": kw}}]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "--coords", "[1, 2]"]

    def test_reconstruct_command_list_value(self):
        kw = self._kw([1, 2])
        kw["type"].name = "string"
        commands = [{"name": "root", "parent": None, "kwargs": {"nums": kw}}]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "--nums", "1", "2"]

    def test_reconstruct_command_multiple_single_nargs(self):
        kw = self._kw(("a", "b"), nargs=1, is_multiple=True)
        kw["type"].name = "string"
        commands = [{"name": "root", "parent": None, "kwargs": {"item": kw}}]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "--item", "a", "--item", "b"]

    def test_reconstruct_command_multiple_multi_nargs(self):
        kw = self._kw((("a", "b"), ("c", "d")), nargs=2, is_multiple=True)
        kw["type"].name = "string"
        commands = [{"name": "root", "parent": None, "kwargs": {"pair": kw}}]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "--pair", "a", "b", "--pair", "c", "d"]

    def test_reconstruct_family_with_children(self):
        commands = [
            {"name": "root", "parent": None, "kwargs": {}},
            {"name": "child", "parent": "root", "kwargs": {}},
        ]
        obj = CtxObjArguments(commands)
        result = obj.reconstruct_command_line()
        assert result == ["root", "child"]


class TestGetSettingFromJobtypeForGaussian:
    """Tests for get_setting_from_jobtype_for_gaussian."""

    def test_none_jobtype_raises(self):
        with pytest.raises(ValueError, match="Jobtype must be provided"):
            get_setting_from_jobtype_for_gaussian(
                MagicMock(), None, None, None, None
            )

    @pytest.mark.parametrize(
        "jobtype,settings_attr",
        [
            ("opt", "opt_settings"),
            ("ts", "ts_settings"),
            ("irc", "irc_settings"),
            ("sp", "sp_settings"),
            ("td", "td_settings"),
            ("wbi", "wbi_settings"),
            ("nci", "nci_settings"),
            ("qmmm", "qmmm_settings"),
            ("neb", "neb_settings"),
        ],
    )
    def test_dispatches_to_expected_settings_method(
        self, jobtype, settings_attr
    ):
        project_settings = MagicMock()
        get_setting_from_jobtype_for_gaussian(
            project_settings, jobtype, None, None, None
        )
        getattr(project_settings, settings_attr).assert_called_once()

    def test_modred_without_coordinates_raises_assertion(self):
        project_settings = MagicMock()
        with pytest.raises(AssertionError):
            get_setting_from_jobtype_for_gaussian(
                project_settings, "modred", None, None, None
            )

    def test_modred_with_coordinates_sets_modred(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_gaussian(
            project_settings, "modred", "[2,3]", None, None
        )
        assert settings.modred == [2, 3]

    def test_scan_single_coordinate(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_gaussian(
            project_settings, "scan", "[2,3]", "0.1", "10"
        )
        assert settings.modred == {
            "coords": [2, 3],
            "num_steps": [10],
            "step_size": [0.1],
        }

    def test_scan_broadcast_multiple_coordinates(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_gaussian(
            project_settings, "scan", "[[2,3],[6,7]]", "0.1", "10"
        )
        assert settings.modred == {
            "coords": [[2, 3], [6, 7]],
            "num_steps": [10, 10],
            "step_size": [0.1, 0.1],
        }

    def test_unsupported_jobtype_returns_none_settings(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_gaussian(
            project_settings, "bogus", None, None, None
        )
        assert settings is None


class TestCheckScanCoordinatesGaussian:
    def test_all_present_does_not_raise(self):
        check_scan_coordinates_gaussian("[1,2]", "0.1", "10")

    def test_missing_raises_assertion(self):
        with pytest.raises(AssertionError):
            check_scan_coordinates_gaussian(None, "0.1", "10")


class TestCheckScanParametersConsistencyGaussian:
    def test_single_coordinate_consistent(self):
        check_scan_parameters_consistency_gaussian([1, 2], [0.1], [10])

    def test_multi_coordinate_consistent(self):
        check_scan_parameters_consistency_gaussian(
            [[1, 2], [3, 4]], [0.1, 0.2], [10, 15]
        )

    def test_invalid_coordinate_format_raises(self):
        with pytest.raises(ValueError, match="Invalid format"):
            check_scan_parameters_consistency_gaussian(["a"], [0.1], [10])

    def test_mismatched_counts_raises(self):
        with pytest.raises(ValueError, match="Mismatch"):
            check_scan_parameters_consistency_gaussian(
                [[1, 2], [3, 4]], [0.1], [10, 15]
            )


class TestGetSettingFromJobtypeForOrca:
    """Tests for get_setting_from_jobtype_for_orca."""

    def test_none_jobtype_raises(self):
        with pytest.raises(ValueError, match="Jobtype must be provided"):
            get_setting_from_jobtype_for_orca(
                MagicMock(), None, None, None, None, None
            )

    @pytest.mark.parametrize(
        "jobtype,settings_attr",
        [
            ("opt", "opt_settings"),
            ("ts", "ts_settings"),
            ("irc", "irc_settings"),
            ("sp", "sp_settings"),
            ("td", "td_settings"),
            ("wbi", "wbi_settings"),
            ("nci", "nci_settings"),
            ("qmmm", "qmmm_settings"),
            ("neb", "neb_settings"),
        ],
    )
    def test_dispatches_to_expected_settings_method(
        self, jobtype, settings_attr
    ):
        project_settings = MagicMock()
        get_setting_from_jobtype_for_orca(
            project_settings, jobtype, None, None, None, None
        )
        getattr(project_settings, settings_attr).assert_called_once()

    def test_modred_without_coordinates_raises_assertion(self):
        project_settings = MagicMock()
        with pytest.raises(AssertionError):
            get_setting_from_jobtype_for_orca(
                project_settings, "modred", None, None, None, None
            )

    def test_modred_with_coordinates_sets_modred(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_orca(
            project_settings, "modred", "[2,3]", None, None, None
        )
        assert settings.modred == [2, 3]

    def test_scan_single_coordinate(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_orca(
            project_settings, "scan", "[2,3]", "3.0", "1.2", "10"
        )
        assert settings.modred == {
            "coords": [2, 3],
            "dist_start": [3.0],
            "dist_end": [1.2],
            "num_steps": [10],
        }

    def test_scan_broadcast_multiple_coordinates(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_orca(
            project_settings,
            "scan",
            "[[2,3],[6,7]]",
            "3.0",
            "1.2",
            "10",
        )
        assert settings.modred == {
            "coords": [[2, 3], [6, 7]],
            "dist_start": [3.0, 3.0],
            "dist_end": [1.2, 1.2],
            "num_steps": [10, 10],
        }

    def test_unsupported_jobtype_returns_none_settings(self):
        project_settings = MagicMock()
        settings = get_setting_from_jobtype_for_orca(
            project_settings, "bogus", None, None, None, None
        )
        assert settings is None


class TestCheckScanCoordinatesOrca:
    def test_all_present_does_not_raise(self):
        check_scan_coordinates_orca("[1,2]", "3.0", "1.2", "10")

    def test_missing_raises_assertion(self):
        with pytest.raises(AssertionError):
            check_scan_coordinates_orca(None, "3.0", "1.2", "10")


class TestCheckScanParametersConsistencyOrca:
    def test_single_coordinate_consistent(self):
        check_scan_parameters_consistency_orca([1, 2], [3.0], [1.2], [10])

    def test_multi_coordinate_consistent(self):
        check_scan_parameters_consistency_orca(
            [[1, 2], [3, 4]], [3.0, 2.5], [1.2, 1.0], [10, 15]
        )

    def test_invalid_coordinate_format_raises(self):
        with pytest.raises(ValueError, match="Invalid format"):
            check_scan_parameters_consistency_orca(["a"], [3.0], [1.2], [10])

    def test_mismatched_counts_raises(self):
        with pytest.raises(ValueError, match="Mismatch"):
            check_scan_parameters_consistency_orca(
                [[1, 2], [3, 4]], [3.0], [1.2, 1.0], [10, 15]
            )


class TestUpdateIrcLabel:
    def test_forward_direction(self):
        assert update_irc_label("mol", "forward", False) == "molf"

    def test_reverse_direction(self):
        assert update_irc_label("mol", "reverse", False) == "molr"

    def test_forward_with_flat_irc(self):
        assert update_irc_label("mol", "forward", True) == "molf_flat"

    def test_none_direction_unchanged(self):
        assert update_irc_label("mol", None, True) == "mol"

    def test_invalid_direction_raises(self):
        with pytest.raises(ValueError, match="Invalid direction"):
            update_irc_label("mol", "sideways", False)


class TestCreateSpLabel:
    def test_with_solvent_model_and_id(self):
        sp_settings = MagicMock(
            solvent_model="smd",
            solvent_id="Di-Methyl,Formamide",
            custom_solvent=None,
        )
        label = create_sp_label("mol", sp_settings)
        assert label == "mol_smd_Di_Methyl_Formamide"

    def test_gas_phase_when_no_solvent(self):
        sp_settings = MagicMock(
            solvent_model=None, solvent_id=None, custom_solvent=None
        )
        label = create_sp_label("mol", sp_settings)
        assert label == "mol_gas_phase"

    def test_custom_solvent_when_no_solvent_model_or_id(self):
        sp_settings = MagicMock(
            solvent_model=None, solvent_id=None, custom_solvent="myslv"
        )
        label = create_sp_label("mol", sp_settings)
        assert label == "mol_custom_solvent"

    def test_partial_solvent_info_leaves_label_unchanged(self):
        sp_settings = MagicMock(
            solvent_model="smd", solvent_id=None, custom_solvent=None
        )
        label = create_sp_label("mol", sp_settings)
        assert label == "mol"
