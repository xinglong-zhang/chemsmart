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

import click

from chemsmart.cli.external import (
    build_external_calculator,
    click_external_backend_options,
    use_external_level_of_theory,
)
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.gaussian.irc import irc
from chemsmart.cli.gaussian.modred import modred
from chemsmart.cli.gaussian.opt import opt
from chemsmart.cli.gaussian.scan import scan
from chemsmart.cli.gaussian.singlepoint import sp
from chemsmart.cli.gaussian.ts import ts
from chemsmart.io.external import DEFAULT_FD_STEP, HESSIAN_MODES
from chemsmart.jobs.gaussian.external import GaussianExternalCalculator
from chemsmart.utils.cli import MyGroup


@gaussian.group("external", cls=MyGroup)
@click_external_backend_options
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
@click.pass_context
def external(
    ctx,
    backend,
    calculator,
    options,
    log_backend_output,
    python_executable,
    hessian,
    fd_step,
):
    """Run a Gaussian job on an External level of theory (ML potential,
    xTB, or any ASE calculator)."""
    spec = build_external_calculator(
        GaussianExternalCalculator,
        options,
        backend=backend,
        calculator=calculator,
        hessian=hessian,
        fd_step=fd_step,
        python_executable=python_executable,
        log_backend_output=log_backend_output,
    )
    use_external_level_of_theory(ctx, spec)


for _command in (opt, ts, modred, scan, irc, sp):
    external.add_command(_command)
