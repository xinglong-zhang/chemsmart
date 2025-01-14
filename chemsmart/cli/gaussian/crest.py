import click
import logging
from chemsmart.utils.utils import check_charge_and_multiplicity
from chemsmart.cli.job import click_job_options
from chemsmart.cli.job import click_gaussian_link_jobs_options
from chemsmart.utils.cli import MyCommand
from chemsmart.cli.gaussian.gaussian import gaussian

logger = logging.getLogger(__name__)


@gaussian.command(cls=MyCommand)
@click_job_options
@click_gaussian_link_jobs_options
@click.option(
    "-j",
    "--jobtype",
    type=str,
    required=True,
    help="Type of job to run for crest.",
)
@click.option(
    "-n",
    "--num-confs-to-run",
    type=int,
    default=None,
    help="Number of conformers to optimize.",
)
@click.pass_context
def crest(
    ctx,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    num_confs_to_run,
    skip_completed,
    **kwargs,
):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    crest_settings = project_settings.main_settings()

    if jobtype.lower() == "opt":
        crest_settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        crest_settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for modred job"
        crest_settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        crest_settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for scan job"
        assert step_size is not None, "Step size must be provided for scan job"
        assert (
            num_steps is not None
        ), "Number of steps must be provided for scan job"
        crest_settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        crest_settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        crest_settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        crest_settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        crest_settings = project_settings.nci_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
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

    if coordinates is not None:
        modred_info = eval(coordinates)
        if jobtype == "modred":
            crest_settings.modred = modred_info
        elif jobtype == "scan":
            scan_info = {
                "coords": modred_info,
                "num_steps": int(num_steps),
                "step_size": float(step_size),
            }
            crest_settings.modred = scan_info

    logger.info(
        f"Crest {type} settings from project: {crest_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian import GaussianCrestJob

    return GaussianCrestJob(
        molecules=molecules,
        settings=crest_settings,
        label=label,
        num_confs_to_run=num_confs_to_run,
        skip_completed=skip_completed,
        **kwargs,
    )
