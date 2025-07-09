import logging

import click

from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@thermochemistry.command(cls=MyCommand)
@click.pass_context
@click.option(
    "-w",
    "--energy-type-for-weighting",
    default="gibbs",
    type=click.Choice(["gibbs", "electronic"]),
    show_default=True,
    help="Type of energy to use for Boltzmann weighting.",
)
def boltzmann(
    ctx,
    energy_type_for_weighting="gibbs",
    outputfile=None,
):
    """Run the Boltzmann weighted averaging for thermochemistry jobs."""
    jobs = ctx.obj.get("jobs", [])
    files = ctx.obj.get("filenames", [])
    job_settings = ctx.obj.get("job_settings", {})
    from chemsmart.jobs.thermochemistry.boltzmann import (
        BoltzmannAverageThermochemistryJob,
    )

    boltzmann_thermochemistry = BoltzmannAverageThermochemistryJob(
        files=files,
        energy_type=energy_type_for_weighting,
        outputfile=outputfile,
        settings=job_settings.copy(),
    )

    boltzmann_thermochemistry.compute_boltzmann_averages()

    logger.info(
        f"Boltzmann-averaged thermochemistry calculation completed for {jobs}."
    )

    boltzmann_thermochemistry.show_results()
