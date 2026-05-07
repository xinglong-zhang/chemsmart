import logging
import os

import click

from chemsmart.cli.gromacs.gromacs import gromacs
from chemsmart.cli.job import click_job_options
from chemsmart.jobs.gromacs.job import GromacsEMJob
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@gromacs.group("em", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click.option(
    "--mdp",
    "mdp_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default=None,
    help="GROMACS .mdp file for energy minimization",
)
@click.option(
    "--structure",
    "-s",
    "structure_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default=None,
    help="Input structure file,such as .gro or .pdb.",
)
@click.option(
    "--top",
    "top_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default=None,
    help="GROMACS topology .top file.",
)
@click.option(
    "--index",
    "index_file",
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    default=None,
    help="Optional GROMACS index .ndx file.",
)
@click.pass_context
def em(
    ctx,
    skip_completed,
    mdp_file,
    structure_file,
    top_file,
    itp_files,
    index_file,
    **kwargs,
):
    """CLI subcommand for running GROMACS energy minimization."""

    jobrunner = ctx.obj.get("jobrunner", None)
    project = ctx.obj.get("project", None)
    filename = ctx.obj.get("filename", None)
    if structure_file is None:
        structure_file = filename
    label = "gromacs_em"
    if structure_file is not None:
        label = os.path.splitext(os.path.basename(structure_file))[0] + "_em"

    logger.info(f"GROMACS EM project: {project}")
    logger.info(f"GROMACS EM structure file: {structure_file}")
    logger.info(f"GROMACS EM mdp file: {mdp_file}")
    logger.info(f"GROMACS EM topology file: {top_file}")
    logger.info(f"GROMACS EM itp files: {itp_files}")
    logger.info(f"GROMACS EM index file: {index_file}")

    if ctx.invoked_subcommand is None:
        return GromacsEMJob(
            molecule=None,
            label=label,
            jobrunner=jobrunner,
            mdp_file=mdp_file,
            structure_file=structure_file,
            top_file=top_file,
            itp_files=itp_files,
            index_file=index_file,
            skip_completed=skip_completed,
            **kwargs,
        )
