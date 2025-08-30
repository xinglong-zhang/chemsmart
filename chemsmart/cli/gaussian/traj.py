import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.jobs.gaussian import GaussianTrajJob
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command(cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click.option(
    "-N",  # avoid conflict with num_steps if scan
    "--num-structures-to-run",
    type=int,
    default=None,
    help="Number of structures from the list of unique structures to run the job on.",
)
@click.option(
    "-g",
    "--grouping-strategy",
    type=click.Choice(
        [
            "rmsd",
            "rcm",
            "fingerprint",
            "isomorphism",
            "formula",
            "connectivity",
        ],
        case_sensitive=False,
    ),
    default="rmsd",
    help="Grouping strategy to use for grouping. \n"
    "Available options are 'rmsd', 'tanimoto', 'isomorphism', 'formula', 'connectivity'",
)
@click.option(
    "-i/",
    "--ignore-hydrogens/--no-ignore-hydrogens",
    type=bool,
    default=False,
    help="Ignore H atoms in the grouping.",
)
@click.option(
    "-p",
    "--num-procs",
    type=int,
    default=4,
    help="Number of processors to use for grouper.",
)
@click.option(
    "-x",
    "--proportion-structures-to-use",
    type=float,
    default=0.1,
    help="Proportion of structures from the end of trajectory to use. \n"
    "Values ranges from 0.0 < x <=1.0. Defaults to 0.1 (last 10% of structures).",
)
@click.pass_context
def traj(
    ctx,
    skip_completed,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    num_structures_to_run,
    grouping_strategy,
    ignore_hydrogens,
    num_procs,
    proportion_structures_to_use,
    **kwargs,
):
    """CLI for running Gaussian set jobs."""

    # get jobrunner for running Gaussian set jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    structure_set_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    structure_set_settings = structure_set_settings.merge(
        job_settings, keywords=keywords
    )

    check_charge_and_multiplicity(structure_set_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]
    label = f"{label}_{jobtype}"
    logger.debug(f"Label for job: {label}")

    logger.info(
        f"Simulated annealing {type} settings from project: {structure_set_settings.__dict__}"
    )

    return GaussianTrajJob(
        molecules=molecules,
        settings=structure_set_settings,
        label=label,
        jobrunner=jobrunner,
        grouping_strategy=grouping_strategy,
        num_procs=num_procs,
        proportion_structures_to_use=proportion_structures_to_use,
        num_structures_to_run=num_structures_to_run,
        ignore_hydrogens=ignore_hydrogens,
        skip_completed=skip_completed,
        **kwargs,
    )
