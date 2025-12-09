import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_jobtype_options, orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("neb", cls=MyCommand)
@click_job_options
@click_orca_jobtype_options
@click.option(
    "-j",
    "--job-type",
    type=click.Choice(
        [
            "NEB",
            "NEB-CI",
            "NEB-TS",
            "FAST-NEB-TS",
            "TIGHT-NEB-TS",
            "LOOSE-NEB",
            "ZOOM-NEB",
            "ZOOM-NEB-CI",
            "ZOOM-NEB-TS",
            "NEB-IDPP",
        ],
        case_sensitive=False,
    ),
    help="Option to run NEB jobs.",
)
@click.option(
    "-s",
    "--starting-xyzfile",
    type=str,
    help="Filename of starting geometry.",
)
@click.option(
    "-e",
    "--ending-xyzfile",
    type=str,
    help="Filename of ending geometry.",
)
@click.option(
    "-i",
    "--intermediate-xyzfile",
    type=str,
    help="Filename of intermediate geometry.",
)
@click.option(
    "-r",
    "--restarting-xyzfile",
    type=str,
    help="Filename of geometry for restarting.",
)
@click.option(
    "-o",
    "--pre-optimization",
    type=bool,
    default=False,
    help="Whether to optimize the input geometries.",
)
@click.option(
    "-c",
    "--charge",
    type=int,
    help="Charge of the starting geometries.",
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    help="Multiplicity of the starting geometries.",
)
@click.option(
    "-A",
    "--semiempirical",
    type=click.Choice(
        [
            "XTB0",
            "XTB1",
            "XTB2",
        ],
        case_sensitive=False,
    ),
    help="Option to run NEB jobs.",
)
@click.pass_context
def neb(
    ctx,
    job_type,
    starting_xyzfile,
    ending_xyzfile,
    intermediate_xyzfile,
    restarting_xyzfile,
    pre_optimization,
    semiempirical,
    **kwargs,
):
    """
    NEB calculation using ORCA.
    """
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    neb_project_settings = project_settings.neb_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py orca -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project irc settings with job settings from cli keywords from cli.orca.py subcommands
    neb_settings = neb_project_settings.merge(job_settings, keywords=keywords)

    # get label for the job
    label = ctx.obj["label"]

    # update neb_settings if any attribute is specified in cli options
    # suppose project has a non None value, and user does not specify a value (None),
    # then the project value should be used and unmodified, ie, should not be merged.
    # update value only if user specifies a value for the attribute:
    neb_settings.jobtype = job_type
    if starting_xyzfile:
        neb_settings.neb_start_xyz = starting_xyzfile
    if ending_xyzfile:
        neb_settings.ending_xyzfile = ending_xyzfile
    if intermediate_xyzfile:
        neb_settings.intermediate_xyzfile = intermediate_xyzfile
    if restarting_xyzfile:
        neb_settings.restarting_xyzfile = restarting_xyzfile
    if pre_optimization:
        neb_settings.preopt_ends = pre_optimization
    if semiempirical:
        neb_settings.semiempirical = semiempirical

    check_charge_and_multiplicity(neb_settings)

    logger.debug(f"Label for job: {label}")

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs
    molecule = molecules[-1]  # get last molecule from list of molecules

    logger.info(f"NEB job settings from project: {neb_settings.__dict__}")

    from chemsmart.jobs.orca.neb import ORCANEBJob

    return ORCANEBJob(
        molecule=molecule,
        settings=neb_settings,
        label=label,
        **kwargs,
    )
