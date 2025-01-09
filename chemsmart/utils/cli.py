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
