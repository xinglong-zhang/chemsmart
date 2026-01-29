"""
ORCA Transition State Search CLI Module

This module provides the command-line interface for ORCA transition state
(TS) search calculations. It supports various TS optimization methods
including OptTS and ScanTS approaches with comprehensive Hessian handling
options.
"""

import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.orca.qmmm_helper import create_orca_qmmm_subcommand
from chemsmart.utils.cli import MyGroup, check_scan_coordinates_orca
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.group("ts", cls=MyGroup, invoke_without_command=True)
@click_job_options
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
    help="List of atoms to use for hybrid Hessian.\n"
    "zero-indexed, e.g. [0, 1, 2, 3]",
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
    jobtype=None,
    coordinates=None,
    dist_start=None,
    dist_end=None,
    num_steps=None,
    inhess=False,
    inhess_filename=None,
    hybrid_hess=False,
    hybrid_hess_atoms=None,
    numhess=False,
    recalc_hess=5,
    trust_radius=None,
    tssearch_type=None,
    full_scan=False,
    skip_completed=True,
    **kwargs,
):
    """
    Run ORCA transition state search calculations.

    This command performs transition state searches using ORCA with support
    for multiple search strategies and Hessian handling options. It can
    perform both direct TS optimization (OptTS) and coordinate scanning
    approaches (ScanTS).

    The calculation uses settings from the project configuration merged
    with command-line overrides. Various Hessian options are available
    including analytical, numerical, hybrid, and read-in from files.
    """
    # get transition state settings from project configuration
    project_settings = ctx.obj["project_settings"]
    ts_project_settings = project_settings.ts_settings()
    logger.debug(f"Loaded TS settings from project: {ts_project_settings}")

    # job setting from filename or default, with updates from user in cli
    # specified in keywords
    # e.g., `chemsmart orca -c <user_charge> -m <user_multiplicity> ts`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project TS settings with job settings from cli keywords
    ts_settings = ts_project_settings.merge(job_settings, keywords=keywords)

    # get label for the job output files
    label = ctx.obj["label"]

    # update ts_settings if any attribute is specified in cli options
    # note: only update value if user explicitly specifies a value for
    # the attribute to preserve project defaults
    if inhess is True:
        ts_settings.inhess = inhess
        logger.debug("Enabled reading Hessian from file")
    if inhess_filename is not None:
        ts_settings.inhess_filename = inhess_filename
        logger.debug(f"Set Hessian filename: {inhess_filename}")
    if hybrid_hess is True:
        ts_settings.hybrid_hess = hybrid_hess
        logger.debug("Enabled hybrid Hessian calculation")
    if hybrid_hess_atoms is not None:
        ts_settings.hybrid_hess_atoms = hybrid_hess_atoms
        logger.debug(f"Set hybrid Hessian atoms: {hybrid_hess_atoms}")
    if numhess is True:
        ts_settings.numhess = numhess
        logger.debug("Enabled numerical Hessian calculation")
    if recalc_hess is not None:
        ts_settings.recalc_hess = recalc_hess
        logger.debug(f"Set Hessian recalculation interval: {recalc_hess}")
    if trust_radius is not None:
        ts_settings.trust_radius = trust_radius
        logger.debug(f"Set trust radius: {trust_radius}")

    jobtype_normalized = (jobtype or "").lower()
    cli_tssearch_type = tssearch_type.lower() if tssearch_type else None
    effective_tssearch_type = ts_settings.tssearch_type or "optts"

    if jobtype_normalized == "scants":
        effective_tssearch_type = "scants"
    if cli_tssearch_type is not None:
        effective_tssearch_type = cli_tssearch_type

    ts_settings.tssearch_type = effective_tssearch_type
    logger.debug(f"Using TS search type: {ts_settings.tssearch_type}")

    is_scants = ts_settings.tssearch_type.lower() == "scants"

    if is_scants:
        label = label.replace("ts", "scants")

        if coordinates is not None:
            missing = [
                name
                for name, value in (
                    ("dist_start", dist_start),
                    ("dist_end", dist_end),
                    ("num_steps", num_steps),
                )
                if value is None
            ]
            if missing:
                raise click.BadParameter(
                    "ScanTS (--tssearch-type scants or -j scants) requires "
                    "--coordinates, --dist-start, --dist-end, and --num-steps."
                )
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
            logger.info(f"Configured ScanTS with scan info: {scan_info}")
        elif ts_settings.scants_modred is None:
            raise click.BadParameter(
                "ScanTS requires scan coordinates via CLI options or project settings."
            )
        else:
            logger.debug(
                "Using ScanTS coordinate settings inherited from the project configuration."
            )
    else:
        label = label.replace("ts", "optts")
        logger.debug("Using OptTS approach")

    if full_scan is True:
        ts_settings.full_scan = full_scan
        logger.debug("Enabled full coordinate scan")

    logger.debug(f"Final job label: {label}")

    # get molecule from context (use the last molecule if multiple)
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]  # get last molecule from list of molecules
    logger.info(f"Running TS search on molecule: {molecule}")

    logger.info(f"Final TS job settings: {ts_settings.__dict__}")

    ctx.obj["parent_skip_completed"] = skip_completed
    ctx.obj["parent_kwargs"] = kwargs
    ctx.obj["parent_settings"] = ts_settings
    ctx.obj["parent_jobtype"] = "ts"

    if ctx.invoked_subcommand is not None:
        return

    # validate charge and multiplicity consistency only for direct ts jobs
    check_charge_and_multiplicity(ts_settings)

    from chemsmart.jobs.orca.ts import ORCATSJob

    return ORCATSJob(
        molecule=molecule,
        settings=ts_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )


create_orca_qmmm_subcommand(ts)
