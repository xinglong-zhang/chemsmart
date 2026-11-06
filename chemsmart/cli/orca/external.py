"""
``chemsmart run/sub orca ... external`` — use ORCA's optimisers (Opt,
OptTS/ScanTS, modred, scan, IRC, sp) on an external level of theory via
ORCA's ``ExtOpt`` interface.

Examples::

    chemsmart run orca -p proj -f ts_guess.xyz -c 0 -m 1 \\
        external -b xtb -o method=gfn2 -o alpb=water ts

    chemsmart run orca -p proj -f mol.xyz -c 0 -m 1 \\
        external -b mace-off -o model=medium opt

    chemsmart run orca -p proj -f mol.xyz -c 0 -m 1 \\
        external -b ase --calculator ase.calculators.emt.EMT opt

ORCA computes the Hessians it needs (TS searches, frequencies, IRC
start) numerically from the external gradients.
"""

import click

from chemsmart.cli.external import (
    build_external_calculator,
    click_external_backend_options,
    use_external_level_of_theory,
)
from chemsmart.cli.orca.irc import irc
from chemsmart.cli.orca.modred import modred
from chemsmart.cli.orca.opt import opt
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.orca.scan import scan
from chemsmart.cli.orca.singlepoint import sp
from chemsmart.cli.orca.ts import ts
from chemsmart.jobs.orca.external import ORCAExternalCalculator
from chemsmart.utils.cli import MyGroup


@orca.group("external", cls=MyGroup)
@click_external_backend_options
@click.pass_context
def external(
    ctx, backend, calculator, options, keep_backend_output, python_executable
):
    """Run an ORCA job on an external level of theory (ML potential,
    xTB, or any ASE calculator) through ORCA's ExtOpt interface."""
    spec = build_external_calculator(
        ORCAExternalCalculator,
        options,
        backend=backend,
        calculator=calculator,
        python_executable=python_executable,
        keep_backend_output=keep_backend_output,
    )
    use_external_level_of_theory(ctx, spec)


for _command in (opt, ts, modred, scan, irc, sp):
    external.add_command(_command)
