import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("neb", cls=MyCommand)
@click_job_options
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
    help="NEB calculation type (e.g., NEB-TS for transition state search).",
)
@click.option(
    "-e",
    "--ending-xyzfile",
    type=str,
    help="Product geometry file path.",
)
@click.option(
    "-i",
    "--intermediate-xyzfile",
    type=str,
    help="Initial TS geometry guess file path.",
)
@click.option(
    "-r",
    "--restarting-xyzfile",
    type=str,
    help="Restart file path for continuation.",
)
@click.option(
    "-o",
    "--pre-optimization",
    type=bool,
    default=False,
    show_default=True,
    help="Whether to optimize the input geometries (accepts True/False).",
)
@click.option(
    "-a",
    "--semiempirical",
    type=click.Choice(
        [
            "XTB0",
            "XTB1",
            "XTB2",
        ],
        case_sensitive=False,
    ),
    help="Semiempirical method for NEB calculation.",
)
@click.option(
    "-n",
    "--nimages",
    type=int,
    help="Number of images used for NEB calculation.",
)
@click.pass_context
def neb(
    ctx,
    job_type,
    ending_xyzfile,
    intermediate_xyzfile,
    restarting_xyzfile,
    pre_optimization,
    semiempirical,
    nimages,
    **kwargs,
):
    """
    Run ORCA Nudged Elastic Band (NEB) calculations.

    Finds minimum energy pathways and transition states by optimizing a chain
    of molecular structures connecting reactant and product geometries.

    The reactant geometry is provided via the main -f option. At minimum, you
    need both reactant (-f) and product (-e) geometries.

    Examples:
        Standard NEB calculation:
        $ chemsmart sub orca -f reactant.xyz neb -j NEB-TS -e product.xyz -A XTB2

        NEB with climbing image for precise TS location:
        $ chemsmart sub orca -f reactant.xyz neb -j NEB-CI -e product.xyz -A XTB1

        NEB with intermediate TS guess:
        $ chemsmart sub orca -f reactant.xyz neb -j NEB-CI -e product.xyz -i ts_guess.xyz

        Restart from previous calculation:
        $ chemsmart sub orca -f reactant.xyz neb -j NEB -r restart.allxyz

        Pre-optimize endpoint geometries:
        $ chemsmart sub orca -f reactant.xyz neb -j NEB-TS -e product.xyz -o -A XTB2

    Args:
        job_type: NEB calculation type (NEB, NEB-CI, NEB-TS, etc.)
        ending_xyzfile: Product geometry file
        intermediate_xyzfile: Intermediate/TS geometry guess file
        restarting_xyzfile: Restart file from previous calculation
        pre_optimization: Pre-optimize endpoint geometries before NEB
        semiempirical: Semiempirical method (XTB0, XTB1, XTB2)

    Note:
        Use NEB-TS or NEB-CI for transition state searches.
        Semiempirical methods (XTB) recommended for initial exploration.
    """
    # get settings from project
    project_settings = ctx.obj["project_settings"]
    neb_project_settings = project_settings.neb_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py orca -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project NEB settings with job settings from cli keywords from cli.orca.py subcommands
    neb_settings = neb_project_settings.merge(job_settings, keywords=keywords)

    # get label for the job
    label = ctx.obj["label"]

    # update neb_settings if any attribute is specified in cli options
    # suppose project has a non None value, and user does not specify a value (None),
    # then the project value should be used and unmodified, ie, should not be merged.
    # update value only if user specifies a value for the attribute:
    neb_settings.jobtype = job_type
    if ending_xyzfile:
        neb_settings.ending_xyzfile = ending_xyzfile
    if intermediate_xyzfile:
        neb_settings.intermediate_xyzfile = intermediate_xyzfile
    if restarting_xyzfile:
        neb_settings.restarting_xyzfile = restarting_xyzfile
    # Always set pre_optimization since it's now a proper boolean flag
    neb_settings.preopt_ends = pre_optimization
    if semiempirical:
        neb_settings.semiempirical = semiempirical
    if nimages:
        neb_settings.nimages = nimages

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
