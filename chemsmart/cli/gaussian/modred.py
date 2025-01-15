import logging
import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.gaussian import gaussian
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.cli import get_setting_from_jobtype
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("modred", cls=MyCommand)
@click_job_options
@click.option('-c', '--coordinates', default=None, help='list of coordinates to be fixed for modredundant job')
@click.pass_context
def modred(ctx, coordinates, **kwargs):
    from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    modred_settings = project_settings.modred_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.gaussian.py subcommands
    modred_settings = modred_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(modred_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    assert coordinates is not None, (
        'Fixed coordinates for modredundant job required!\n'
        'Use the flags `-c` for coordinates.\n'
        'Example usage: `-c "[[2,3],[6,7]]"`'
    )
    modred_info = eval(coordinates)
    modred_settings.modred = modred_info

    logger.info(f"IRC settings from project: {modred_settings.__dict__}")

    from chemsmart.jobs.gaussian import GaussianModredJob

    return GaussianModredJob(
        molecule=molecule, settings=modred_settings, label=label, **kwargs
    )
