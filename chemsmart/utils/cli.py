"""Utility functions for command line interface."""

import copy
import logging
import typing as t

import click
from click import Context

logger = logging.getLogger(__name__)


def _add_subcommand_info_to_ctx(ctx):
    if "subcommand" not in ctx.obj:
        ctx.obj["subcommand"] = []

    subcommand = {}
    subcommand["name"] = ctx.info_name

    # Order of ctx.params may not be the same as order of ctx.command.params!
    nargs = {param.name: param.nargs for param in ctx.command.params}
    is_multiple = {param.name: param.multiple for param in ctx.command.params}
    types = {param.name: param.type for param in ctx.command.params}
    is_flag = {param.name: param.is_flag for param in ctx.command.params}
    secondary_opts = {
        param.name: param.secondary_opts for param in ctx.command.params
    }

    subcommand["kwargs"] = {}
    for param, value in ctx.params.items():
        subcommand["kwargs"][param] = {}
        subcommand["kwargs"][param]["value"] = value
        subcommand["kwargs"][param]["nargs"] = nargs[param]
        subcommand["kwargs"][param]["is_multiple"] = is_multiple[param]
        subcommand["kwargs"][param]["type"] = types[param]
        subcommand["kwargs"][param]["is_flag"] = is_flag[param]
        subcommand["kwargs"][param]["secondary_opts"] = secondary_opts[param]

    parent = ctx.parent.info_name if ctx.parent is not None else None

    subcommand["parent"] = parent

    ctx.obj["subcommand"] += [subcommand]


class MyGroup(click.Group):
    """For click commands to store subcommand information."""

    def invoke(self, ctx: Context) -> t.Any:
        """Add subcommand information to context before invoking."""
        _add_subcommand_info_to_ctx(ctx)
        return super().invoke(ctx)


class MyCommand(click.Command):
    def invoke(self, ctx):
        _add_subcommand_info_to_ctx(ctx)
        return super().invoke(ctx)


class CtxObjArguments:
    def __init__(self, commands, entry_point=None):
        self.commands = copy.copy(commands)
        self.entry_point = entry_point

    @staticmethod
    def _value(v):
        if v is True or v is False:
            return ""
        return str(v)

    @staticmethod
    def argument_keyword(arg, value, is_flag, secondary_opts):
        arg = arg.replace("_", "-")

        if value is False:
            if is_flag and len(secondary_opts) == 0:
                # click "is_flag" option. Click param class has no way to check if flag is True or not
                # only can determine by checking arg_secondary_opts
                return ""

            # click "bool" variable e.g. '--with-flag/--no-with-flag'
            assert len(secondary_opts) == 1
            arg = secondary_opts[0].strip("-")

        return "-" + arg if len(arg) == 1 else "--" + arg

    @property
    def _subcommands(self):
        return self.commands

    def _entry_point(self):
        if self.entry_point is not None:
            subcommand = [
                s for s in self._subcommands if s["name"] == self.entry_point
            ]
            if len(subcommand) == 0:
                raise ValueError(
                    f"Specified entry point: {self.entry_point} could not be found"
                )
            assert len(subcommand) == 1
            return subcommand[0]

        for subcommand in self._subcommands:
            if subcommand["parent"] is None:
                return subcommand
        raise AssertionError("Could not find entry point")

    def _reconstruct_command(self, command):
        subcommand_args = command["kwargs"]
        subcommand_name = command["name"]
        # logger.debug(f'Subcommand args {pprint.pformat(subcommand_args)}')

        command_line_string = [subcommand_name]
        for k, subdict in subcommand_args.items():
            v = subdict["value"]
            multiple_arg = subdict["is_multiple"]
            arg_type = subdict["type"]
            arg_nargs = subdict["nargs"]
            arg_is_flag = subdict["is_flag"]
            arg_secondary_opts = subdict["secondary_opts"]

            if v is None or (isinstance(v, tuple) and len(v) == 0):
                continue

            keyword = self.argument_keyword(
                k, v, arg_is_flag, arg_secondary_opts
            )
            if multiple_arg:
                logger.debug(f"Keyword {keyword}, value {v}")
                for i in v:
                    command_line_string += [keyword]
                    if arg_nargs == 1:
                        command_line_string += [self._value(i)]
                    else:
                        command_line_string += [self._value(j) for j in i]
            else:
                command_line_string += [keyword]
                if arg_type.name == "literal_eval":
                    command_line_string += [self._value(v)]
                elif isinstance(v, (list, tuple)):
                    command_line_string += [self._value(i) for i in v]
                else:
                    command_line_string += [self._value(v)]
        return command_line_string

    def _reconstruct_family(self, parent):
        children = self._subcommands_of(parent)
        logger.debug(
            f'Parent: {parent["name"]};  children: {[c["name"] for c in children]}'
        )

        parent_command_line = self._reconstruct_command(parent)
        children_command_line = []
        for child in children:
            child_command_line = self._reconstruct_family(child)
            children_command_line.extend(child_command_line)
        return parent_command_line + children_command_line

    def _subcommands_of(self, parent):
        return [
            subcommand
            for subcommand in self._subcommands
            if subcommand["parent"] == parent["name"]
        ]

    def reconstruct_command_line(self):
        parent = self._entry_point()
        command_line = self._reconstruct_family(parent)
        # logger.info(f'Cmd: {command_line}')
        return [i for i in command_line if len(i) != 0]


def get_setting_from_jobtype_for_gaussian(
    project_settings, jobtype, coordinates, step_size, num_steps
):
    if jobtype is None:
        raise ValueError("Jobtype must be provided for Crest and Link job.")

    settings = None

    if jobtype.lower() == "opt":
        settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for modred job."
        settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        assert all(
            v is not None for v in [coordinates, step_size, num_steps]
        ), (
            "Scanning coordinates, step size and number of steps of scan required!\n"
            "Use the flags `-c -s -n` for coordinates, step-size and num-steps respectively.\n"
            "Example usage: `-c [[2,3],[6,7]] -s 0.1 -n 15`"
        )
        settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        settings = project_settings.nci_settings()

    if coordinates is not None:
        modred_info = eval(coordinates)
        if jobtype == "modred":
            settings.modred = modred_info
        elif jobtype == "scan":
            scan_info = {
                "coords": modred_info,
                "num_steps": int(num_steps),
                "step_size": float(step_size),
            }
            settings.modred = scan_info

    return settings

def get_setting_from_jobtype_for_orca(
    project_settings, jobtype, coordinates, step_size, num_steps
):
    if jobtype is None:
        raise ValueError("Jobtype must be provided for Crest and Link job.")

    settings = None

    if jobtype.lower() == "opt":
        settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for modred job."
        settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        assert all(
            v is not None for v in [coordinates, step_size, num_steps]
        ), (
            "Scanning coordinates, step size and number of steps of scan required!\n"
            "Use the flags `-c -s -n` for coordinates, step-size and num-steps respectively.\n"
            "Example usage: `-c [[2,3],[6,7]] -s 0.1 -n 15`"
        )
        settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        settings = project_settings.nci_settings()

    if coordinates is not None:
        modred_info = eval(coordinates)
        if jobtype == "modred":
            settings.modred = modred_info
        elif jobtype == "scan":
            scan_info = {
                "coords": modred_info,
                "num_steps": int(num_steps),
                "step_size": float(step_size),
            }
            settings.modred = scan_info

    return settings



