import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@gaussian.command("qmmm", cls=MyCommand)
@click_job_options
@click.option(
    "-fh",
    "--functional-high",
    type=str,
    help="functional of high-level layer.",
)
@click.option(
    "-bh",
    "--basis-high",
    type=str,
    help="functional of high-level layer.",
)
@click.option(
    "-fm",
    "--functional-medium",
    type=str,
    help="functional of medium-level layer.",
)
@click.option(
    "-bm",
    "--basis-medium",
    type=str,
    help="basis of medium-level layer.",
)
@click.option(
    "-fl",
    "--functional-low",
    type=str,
    help="functional of low-level layer.",
)
@click.option(
    "-bl",
    "--basis-low",
    type=str,
    help="basis of low-level layer.",
)
@click.option(
    "-ff",
    "--force field",
    type=str,
    help="name of force field.",
)
@click.option(
    "-hlc",
    "--model-high-level-charges",
    type=str,
    help="charge for model system, high-level.",
)
@click.option(
    "-hlm",
    "--model-high-level-multiplicity",
    type=str,
    help="multiplicity for model system, high-level.",
)
@click.option(
    "-llc",
    "--model-low-level-charges",
    type=str,
    help="charge for model system, low-level.",
)
@click.option(
    "-mlm",
    "--model-low-level-multiplicity",
    type=str,
    help="multiplicity for model system, low-level.",
)
@click.option(
    "-rlc",
    "--real-low-level-charges",
    type=str,
    help="charge for real system, low-level.",
)
@click.option(
    "-llm",
    "--real-low-level-multiplicity",
    type=str,
    help="multiplicity for model system, low-level.",
)
@click.option(
    "-ha",
    "--high-level-atom",
    type=str,
    help="atom index for high level.",
)
@click.option(
    "-ma",
    "--medium-level-atom",
    type=str,
    help="atom index for medium level.",
)
@click.option(
    "-la",
    "--low-level-atom",
    type=str,
    help="atom index for low level.",
)
@click.option(
    "-b",
    "--bonded-atom",
    type=str,
    help="the bond to be cut, specified by two atomic indexes .",
)
@click.option(
    "-sf1",
    "--scale-factor1",
    type=str,
    help="Scale factor for bonds between QM and MM region,default=1.0.",
)
@click.option(
    "-sf2",
    "--scale-factor2",
    type=str,
    help="Scale factor for angles involving  MM and MM region,default=1.0.",
)
@click.option(
    "-sf3",
    "--scale-factor3",
    type=str,
    help="Scale factor for torsions, default=1.0",
)
@click.pass_context
def cli(ctx, **kwargs):
    """Main CLI Command"""
    ctx.ensure_object(dict)
    ctx.obj["settings"] = kwargs


@click.pass_context
def qmmm(ctx):
    """Run the QMMM Job"""
    settings = ctx.obj["settings"]

    from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings

    job_settings = GaussianQMMMJobSettings(**settings)

    click.echo(f"Running QM/MM job with settings: {job_settings.__dict__}")
