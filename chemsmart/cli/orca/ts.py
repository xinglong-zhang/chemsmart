import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_jobtype_options, orca
from chemsmart.utils.cli import MyCommand, check_scan_coordinates_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("ts", cls=MyCommand)
@click_job_options
@click_orca_jobtype_options
@click.option(
    "-i/",
    "--inhess/--no-inhess",
    type=bool,
    default=False,
    help="Option to read in Hessian file.",
)
@click.option(
    "-f",
    "--inhess-filename",
    type=str,
    default=None,
    help="Filename of Hessian file.",
)
@click.option(
    "-h/",
    "--hybrid-hess/--no-hybrid-hess",
    type=bool,
    default=False,
    help="Option to use hybrid Hessian.",
)
@click.option(
    "-a",
    "--hybrid-hess-atoms",
    default=None,
    help="List of atoms to use for hybrid Hessian.\nzero-indexed, e.g. [0, 1, 2, 3]",
)
@click.option(
    "--numhess/--no-numhess",
    type=bool,
    default=False,
    help="Option to use numerical Hessian.",
)
@click.option(
    "-s",
    "--recalc-hess",
    type=int,
    default=5,
    help="Number of steps to recalculate Hessian.",
)
@click.option(
    "-t",
    "--trust-radius",
    type=float,
    default=None,
    help="Trust radius for TS optimization.",
)
@click.option(
    "-ts",
    "--tssearch-type",
    type=str,
    default="optts",
    help='Type of TS search to perform. Options are ["optts", "scants"]',
)
@click.option(
    "-fs/",
    "--full-scan/--no-full-scan",
    type=bool,
    default=False,
    help="Option to perform a full scan.",
)
@click.pass_context
def ts(
    ctx,
    jobtype,
    coordinates,
    dist_start,
    dist_end,
    num_steps,
    inhess,
    inhess_filename,
    hybrid_hess,
    hybrid_hess_atoms,
    numhess,
    recalc_hess,
    trust_radius,
    tssearch_type,
    full_scan,
    **kwargs,
):
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    ts_project_settings = project_settings.ts_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py orca -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project irc settings with job settings from cli keywords from cli.orca.py subcommands
    ts_settings = ts_project_settings.merge(job_settings, keywords=keywords)

    # get label for the job
    label = ctx.obj["label"]

    # update irc_settings if any attribute is specified in cli options
    # suppose project has a non None value, and user does not specify a value (None),
    # then the project value should be used and unmodified, ie, should not be merged.
    # update value only if user specifies a value for the attribute:
    if inhess is True:
        ts_settings.inhess = inhess
    if inhess_filename is not None:
        ts_settings.inhess_filename = inhess_filename
    if hybrid_hess is True:
        ts_settings.hybrid_hess = hybrid_hess
    if hybrid_hess_atoms is not None:
        ts_settings.hybrid_hess_atoms = hybrid_hess_atoms
    if numhess is True:
        ts_settings.numhess = numhess
    if recalc_hess is not None:
        ts_settings.recalc_hess = recalc_hess
    if trust_radius is not None:
        ts_settings.trust_radius = trust_radius
    if tssearch_type is not None:
        ts_settings.tssearch_type = tssearch_type

    if tssearch_type.lower() == "scants":
        check_scan_coordinates_orca(
            coordinates, dist_start, dist_end, num_steps
        )
        coordinates = ast.literal_eval(coordinates)
        scan_info = {
            "coordinates": coordinates,
            "dist_start": dist_start,
            "dist_end": dist_end,
            "num_steps": num_steps,
        }
        ts_settings.scants_modred = scan_info
        label = label.replace("ts", "scants")
    else:
        label = label.replace("ts", "optts")

    if full_scan is True:
        ts_settings.full_scan = full_scan

    check_charge_and_multiplicity(ts_settings)

    logger.debug(f"Label for job: {label}")

    # get molecule
    molecules = ctx.obj[
        "molecules"
    ]  # use all molecules as a list for crest jobs
    molecule = molecules[-1]  # get last molecule from list of molecules

    logger.info(f"TS job settings from project: {ts_settings.__dict__}")

    from chemsmart.jobs.orca.ts import ORCATSJob

    return ORCATSJob(
        molecule=molecule,
        settings=ts_settings,
        label=label,
        **kwargs,
    )
