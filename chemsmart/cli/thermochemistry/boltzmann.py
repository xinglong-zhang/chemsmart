import logging

import click

from chemsmart.cli.thermochemistry.thermochemistry import thermochemistry
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@thermochemistry.command(cls=MyCommand)
@click.option(
    "-w",
    "--energy-type-for-weighting",
    default="gibbs",
    type=click.Choice(["gibbs", "electronic"]),
    show_default=True,
    help="Type of energy to use for Boltzmann weighting.",
)
@click.option(
    "-o",
    "--outputfile",
    default=None,
    type=str,
    help="Output file to save the Boltzmann averaged thermochemistry results.",
)
@click.pass_context
def boltzmann(
    ctx,
    energy_type_for_weighting="gibbs",
    outputfile=None,
):
    """Run the Boltzmann weighted averaging for thermochemistry jobs."""
    jobs = ctx.obj.get("jobs", [])
    files = ctx.obj.get("filenames", [])
    job_settings = ctx.obj.get("job_settings", {})
    from chemsmart.analysis.thermochemistry import (
        BoltzmannAverageThermochemistry,
    )

    boltzmann_thermochemistry = BoltzmannAverageThermochemistry(
        files=files,
        energy_type=energy_type_for_weighting,
        outputfile=outputfile,
        temperature=job_settings.temperature,
        concentration=job_settings.concentration,
        pressure=job_settings.pressure,
        use_weighted_mass=job_settings.use_weighted_mass,
        alpha=job_settings.alpha,
        s_freq_cutoff=job_settings.s_freq_cutoff,
        h_freq_cutoff=job_settings.h_freq_cutoff,
        energy_units=job_settings.energy_units,
    )

    boltzmann_thermochemistry.compute_boltzmann_averages()

    logger.info(
        f"Boltzmann-averaged thermochemistry calculation completed for {jobs}."
    )
