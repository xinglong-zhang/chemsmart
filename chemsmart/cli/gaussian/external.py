"""
``chemsmart run/sub gaussian ... external`` — use Gaussian's optimisers
(opt, TS, modred, scan, IRC, sp) on an External level of theory.

Examples::

    chemsmart run gaussian -p proj -f ts_guess.xyz -c 0 -m 1 \\
        external -b xtb -o method=gfn2 -o alpb=water ts

    chemsmart run gaussian -p proj -f mol.xyz -c 0 -m 1 \\
        external -b mace-off -o model=medium opt

    chemsmart run gaussian -p proj -f mol.xyz -c 0 -m 1 \\
        external -b ase --calculator ase.calculators.emt.EMT opt
"""

import logging

import click
import yaml

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.gaussian.irc import irc
from chemsmart.cli.gaussian.modred import modred
from chemsmart.cli.gaussian.opt import opt
from chemsmart.cli.gaussian.scan import scan
from chemsmart.cli.gaussian.singlepoint import sp
from chemsmart.cli.gaussian.ts import ts
from chemsmart.io.gaussian.external import DEFAULT_FD_STEP, HESSIAN_MODES
from chemsmart.jobs.gaussian.external import (
    EXTERNAL_BACKENDS,
    GaussianExternalCalculator,
)
from chemsmart.utils.cli import MyGroup
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


@gaussian.group("external", cls=MyGroup)
@click.option(
    "-b",
    "--backend",
    type=click.Choice(EXTERNAL_BACKENDS, case_sensitive=False),
    required=True,
    help="Level of theory evaluated outside Gaussian.",
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
    "--hessian",
    type=click.Choice(HESSIAN_MODES),
    default="auto",
    show_default=True,
    help="Hessian source when Gaussian requests one (TS calcfc, freq, "
    "IRC): auto uses the backend's analytic Hessian when available and "
    "consistent, otherwise finite differences.",
)
@click.option(
    "--fd-step",
    type=float,
    default=DEFAULT_FD_STEP,
    show_default=True,
    help="Finite-difference displacement for Hessians, in Angstrom.",
)
@click.option(
    "--python",
    "python_executable",
    default=None,
    help="Python interpreter Gaussian uses to run the External launcher "
    "(default: the interpreter running chemsmart).",
)
@click.pass_context
def external(
    ctx, backend, calculator, options, hessian, fd_step, python_executable
):
    """Run a Gaussian job on an External level of theory (ML potential,
    xTB, or any ASE calculator)."""
    try:
        spec = GaussianExternalCalculator(
            backend=backend,
            calculator=calculator,
            options=parse_backend_options(options),
            hessian=hessian,
            fd_step=fd_step,
            python_executable=python_executable,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    logger.info(f"External level of theory: {spec}")
    ctx.obj["job_settings"].external = spec
    ctx.obj["keywords"] = tuple(ctx.obj["keywords"]) + ("external",)

    # Auto-generated labels end in "_external"; name them after the real
    # job and backend instead (e.g. mol_ts_xtb), so External results never
    # overwrite DFT results for the same structure.
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


for _command in (opt, ts, modred, scan, irc, sp):
    external.add_command(_command)
