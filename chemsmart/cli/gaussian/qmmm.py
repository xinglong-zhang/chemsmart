import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.gaussian.irc import irc
from chemsmart.cli.gaussian.modred import modred
from chemsmart.cli.gaussian.opt import opt
from chemsmart.cli.gaussian.scan import scan
from chemsmart.cli.gaussian.ts import ts
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


# @click.group(cls=MyGroup)
@gaussian.group("qmmm", cls=MyGroup)
@click.pass_context
def qmmm(
    ctx,
    # additional options for qmmm job
):
    """CLI for running Gaussian jobs using the chemsmart framework."""

    print("QMMM command group initialized")


qmmm.add_command(opt)
qmmm.add_command(ts)
qmmm.add_command(modred)
qmmm.add_command(irc)
qmmm.add_command(scan)

# @qmmm.result_callback()
# @click.pass_context
# def qmmm_process_pipeline(ctx, *args, **kwargs):
#     qmmm.add_command(ctx.invoked_subcommand)
#     kwargs.update({"subcommand": ctx.invoked_subcommand})
#     ctx.obj[ctx.info_name] = kwargs
#     return args[0]
