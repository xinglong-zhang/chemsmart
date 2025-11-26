import logging

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command("userjob", cls=MyCommand)
@click_job_options
@click.option(
    "-r", "--route", required=True, type=str, help="user-defined route"
)
@click.option(
    "-a",
    "--append-info",
    type=str,
    default=None,
    help="information to be appended at the end of the file",
)
@click.pass_context
def userjob(ctx, route, append_info, **kwargs):
    """CLI for running Gaussian custom jobs."""

    # get jobrunner for running Gaussian custom jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    opt_settings.route_to_be_written = route

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from
    # cli.gaussian.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    check_charge_and_multiplicity(opt_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[
        -1
    ]  # get last molecule from list of molecules from cli.gaussian.py subcommands
    # index = '-1' would access the right structure from the list of molecule
    # returned from cli.gaussian.py subcommands
    # user specified index was used there to return the right molecule and
    # store it as a list of single element/itself

    # get label for the job
    label = ctx.obj["label"]

    if append_info is not None:
        opt_settings.append_additional_info = append_info.encode().decode(
            "unicode-escape"
        )

    from chemsmart.jobs.gaussian.custom import GaussianCustomJob

    return GaussianCustomJob(
        molecule=molecule,
        settings=opt_settings,
        label=label,
        jobrunner=jobrunner,
        **kwargs,
    )
