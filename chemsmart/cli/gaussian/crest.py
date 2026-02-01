import functools
import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.grouper.grouper import click_grouper_common_options
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


def click_crest_grouper_options(f):
    """Grouper options specific to crest command (strategy selection and strategy-specific options)."""

    @click.option(
        "-g",
        "--grouping-strategy",
        type=click.Choice(
            [
                "rmsd",
                "hrmsd",
                "spyrmsd",
                "irmsd",
                "pymolrmsd",
                "tanimoto",
                "torsion",
                "isomorphism",
                "formula",
                "connectivity",
            ],
            case_sensitive=False,
        ),
        default=None,
        help="Grouping strategy to use for conformer grouping.",
    )
    @click.option(
        "--check-stereo",
        type=click.Choice(["auto", "on", "off"], case_sensitive=False),
        default="auto",
        help="Control stereochemistry/inversion checking in iRMSD grouper. "
        "'auto' (default): automatically detect, 'on': force check, 'off': disable.",
    )
    @click.option(
        "-ft",
        "--fingerprint-type",
        type=click.Choice(
            [
                "rdkit",
                "rdk",
                "morgan",
                "maccs",
                "atompair",
                "torsion",
                "usr",
                "usrcat",
            ],
            case_sensitive=False,
        ),
        default="rdkit",
        help="Fingerprint type for tanimoto grouping.",
    )
    @click.option(
        "--use-weights/--no-use-weights",
        type=bool,
        default=True,
        help="Whether to use torsion weights in TFD calculation for torsion grouping.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@gaussian.command(cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click_grouper_common_options
@click_crest_grouper_options
@click.pass_context
def crest(
    ctx,
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
    skip_completed,
    **kwargs,
):
    """CLI subcommand for running Gaussian CREST jobs."""

    # get jobrunner for running Gaussian crest jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    crest_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    crest_settings = crest_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(crest_settings)

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs

    # get label for the job
    label = ctx.obj["label"]
    label = f"{label}_{jobtype}"
    logger.debug(f"Label for job: {label}")

    logger.info(
        f"Crest {jobtype} settings from project: {crest_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.crest import GaussianCrestJob

    logger.debug(f"Creating GaussianCrestJob with {len(molecules)} molecules")

    return GaussianCrestJob(
        molecules=molecules,
        settings=crest_settings,
        label=label,
        jobrunner=jobrunner,
        grouping_strategy=grouping_strategy,
        ignore_hydrogens=ignore_hydrogens,
        threshold=threshold,
        num_groups=num_groups,
        num_procs=num_procs,
        check_stereo=check_stereo,
        use_weights=use_weights,
        fingerprint_type=fingerprint_type,
        skip_completed=skip_completed,
        **kwargs,
    )
