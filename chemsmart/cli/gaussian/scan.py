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
    "--constrained-coordinates",
    default=None,
    help="Additional modredundant constraints for scan jobs. "
    "Format: List of constraints separated by semicolons. "
    "Example: [[1,2],[3,4,5],[1,2,3,4]]. 1-indexed.",
)
@click.pass_context
def scan(
    ctx,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    constrained_coordinates=None,
    **kwargs,
):
    """CLI subcommand for running Gaussian scan jobs."""

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

    if constrained_coordinates is not None:
        constrained_coordinates_info = ast.literal_eval(
            constrained_coordinates
        )
        scan_settings.modred["constrained_coordinates"] = (
            constrained_coordinates_info
        )
    # get molecules
    molecules = ctx.obj["molecules"]

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    logger.info(f"Scan job settings from project: {scan_settings.__dict__}")

    from chemsmart.jobs.gaussian.scan import GaussianScanJob

    # Get the original molecule indices from context
    molecule_indices = ctx.obj["molecule_indices"]

    # Handle multiple molecules: create one job per molecule
    if len(molecules) > 1 and molecule_indices is not None:
        logger.info(f"Creating {len(molecules)} scan jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            logger.info(
                f"Running scan for molecule {idx}: {molecule} with label {molecule_label}"
            )

            job = GaussianScanJob(
                molecule=molecule,
                settings=scan_settings,
                label=molecule_label,
                jobrunner=jobrunner,
                **kwargs,
            )
            jobs.append(job)
        return jobs
    else:
        # Single molecule case
        molecule = molecules[-1]
        return GaussianScanJob(
            molecule=molecule,
            settings=scan_settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
