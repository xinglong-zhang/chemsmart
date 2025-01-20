import click
import logging
from chemsmart.cli.job import click_job_options
from chemsmart.io.orca.input import ORCAInput
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity
from chemsmart.cli.orca.orca import orca

logger = logging.getLogger(__name__)


@orca.command("inp", cls=MyCommand)
@click_job_options
@click.pass_context
def inp(ctx, skip_completed, **kwargs):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.gaussian.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(opt_settings)
    filename = ctx.obj["filename"]
    if filename is not None:
        opt_settings.input_string = ORCAInput(
            filename=filename
        ).content_lines_string

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]

    # get label for the job
    label = ctx.obj["label"]

    from chemsmart.jobs.orca.job import ORCAInpJob

    return ORCAInpJob(
        molecule=molecule,
        settings=opt_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
