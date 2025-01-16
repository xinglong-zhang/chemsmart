"""Utility functions for command line interface."""

import typing as t

import click
from click import Context


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


def get_setting_from_jobtype(
    project_settings, jobtype, coordinates, step_size, num_steps
):
    if jobtype is None:
        raise ValueError("Jobtype must be provided for Crest and Link job.")

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
