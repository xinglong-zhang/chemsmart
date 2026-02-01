import logging

import click

from chemsmart.cli.gaussian.crest import click_crest_grouper_options
from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.grouper.grouper import click_grouper_common_options
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
@click_grouper_common_options
@click_crest_grouper_options
@click.option(
    "-x",
    "--proportion-structures-to-use",
    type=float,
    default=0.1,
    help="Proportion of structures from the end of trajectory to use. \n"
    "Values ranges from 0.0 < x <=1.0. Defaults to 0.1 (last 10% of "
    "structures).",
)
@click.pass_context
def traj(
    ctx,
    skip_completed,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    # grouper common options
    ignore_hydrogens,
    num_procs,
    threshold,
    num_groups,
    # crest grouper options
    grouping_strategy,
    check_stereo,
    fingerprint_type,
    use_weights,
    proportion_structures_to_use,
    **kwargs,
):
    """CLI subcommand for running Gaussian set jobs."""

    # get jobrunner for running Gaussian set jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    structure_set_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
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
        f"Simulated annealing {jobtype} settings from project: "
        f"{structure_set_settings.__dict__}"
    )

    return GaussianTrajJob(
        molecules=molecules,
        settings=structure_set_settings,
        label=label,
        jobrunner=jobrunner,
        grouping_strategy=grouping_strategy,
        num_procs=num_procs,
        proportion_structures_to_use=proportion_structures_to_use,
        num_groups=num_groups,
        ignore_hydrogens=ignore_hydrogens,
        use_weights=use_weights,
        threshold=threshold,
        fingerprint_type=fingerprint_type,
        check_stereo=check_stereo,
        skip_completed=skip_completed,
        **kwargs,
    )
