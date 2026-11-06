"""
CLI pieces shared by ``gaussian ... external`` and ``orca ... external``.

Both groups take the same backend options, store the backend on the
parent's job settings, and rename auto-generated labels after the real
job and backend.
"""

import functools
import logging

import click
import yaml

from chemsmart.io.external import available_backends, get_backend
from chemsmart.utils.io import clean_label

logger = logging.getLogger(__name__)


def parse_backend_options(options):
    """Turn ``("key=value", ...)`` into a dict with YAML-typed values."""
    parsed = {}
    for item in options:
        key, sep, value = item.partition("=")
        if not sep or not key.strip():
            raise click.BadParameter(
                f"{item!r} is not of the form key=value.",
                param_hint="'-o/--option'",
            )
        parsed[key.strip()] = yaml.safe_load(value)
    return parsed


class ExternalBackendType(click.ParamType):
    """Backend name checked against the registry when parsed, so backends
    registered by plugins are accepted without changing the CLI."""

    name = "backend"

    def convert(self, value, param, ctx):
        try:
            return get_backend(value).name
        except ValueError as exc:
            self.fail(str(exc), param, ctx)

    def get_metavar(self, param, *args, **kwargs):
        return f"[{'|'.join(available_backends())}]"


def click_external_backend_options(f):
    """Backend options common to every external level of theory."""

    @click.option(
        "-b",
        "--backend",
        type=ExternalBackendType(),
        required=True,
        help="Level of theory evaluated outside the host program.",
    )
    @click.option(
        "--calculator",
        default=None,
        help="ASE calculator class or factory as a dotted path, e.g. "
        "ase.calculators.emt.EMT (required for --backend ase).",
    )
    @click.option(
        "-o",
        "--option",
        "options",
        multiple=True,
        help="Backend option key=value (repeatable), e.g. -o model=medium "
        "-o device=cuda, or for xtb -o method=gfn2 -o alpb=water.",
    )
    @click.option(
        "--log-backend-output/--no-log-backend-output",
        default=False,
        show_default=True,
        help="Also copy the backend program's own output (e.g. xtb's) "
        "into <label>_external.log. One line per request is always "
        "logged there.",
    )
    @click.option(
        "--python",
        "python_executable",
        default=None,
        help="Python interpreter that runs the external launcher "
        "(default: the interpreter running chemsmart).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def build_external_calculator(calculator_cls, options, **kwargs):
    """Instantiate ``calculator_cls``, turning bad input into CLI errors."""
    try:
        return calculator_cls(options=parse_backend_options(options), **kwargs)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def use_external_level_of_theory(ctx, spec):
    """Put ``spec`` on the parent job settings for the chosen subcommand.

    Auto-generated labels end in ``_external``; they are renamed after the
    real job and backend (e.g. ``mol_ts_xtb``) so external results never
    overwrite results of the same structure at the host's own level of
    theory. Labels given with ``-l``/``-a`` are kept.
    """
    logger.info(f"External level of theory: {spec}")
    ctx.obj["job_settings"].external = spec
    ctx.obj["keywords"] = tuple(ctx.obj["keywords"]) + ("external",)

    parent_params = ctx.parent.params
    if (
        parent_params.get("label") is None
        and parent_params.get("append_label") is None
        and ctx.invoked_subcommand is not None
    ):
        stem = ctx.obj["label"].removesuffix("_external")
        ctx.obj["label"] = clean_label(
            f"{stem}_{ctx.invoked_subcommand}_{spec.backend}"
        )
