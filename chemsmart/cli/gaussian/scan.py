import ast
import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("scan", cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click.option(
    "-cc",
    "--const-coord",
    default=None,
    help="Additional modredundant constraints for scan job. "
         "Format: List of constraints separated by semicolons. "
         "Example: [[1,2],[3,4,5],[1,2,3,4]]. "
         "1-indexed.",
)
@click.pass_context
def scan(ctx, jobtype, coordinates, step_size, num_steps, const_coord=None, **kwargs):
    """CLI for running Gaussian scan jobs."""

    # get jobrunner for running Gaussian scan jobs
    jobrunner = ctx.obj["jobrunner"]

    if jobtype is None:
        jobtype = "scan"

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    scan_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    scan_settings = scan_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(scan_settings)

    if const_coord is not None:
        const_coord_info = ast.literal_eval(const_coord)
        scan_settings.modred['const_coords'] = const_coord_info
    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"Scan job settings from project: {scan_settings.__dict__}")

    from chemsmart.jobs.gaussian.scan import GaussianScanJob

    return GaussianScanJob(
        molecule=molecule,
        settings=scan_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
