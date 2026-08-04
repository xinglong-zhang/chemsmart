"""
GROMACS staged workflow CLI command.

This command creates the standard sequential workflow:

    EM -> NVT -> NPT -> MD

Example:
    chemsmart run gromacs workflow \
        --structure input.gro \
        --top topol.top \
        --output-folder workflow
"""

from __future__ import annotations

import logging
from pathlib import Path

import click

from chemsmart.cli.gromacs.gromacs import gromacs
from chemsmart.cli.job import click_job_options
from chemsmart.jobs.gromacs.workflow import GromacsWorkflow
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@gromacs.group(
    "workflow",
    cls=MyGroup,
    invoke_without_command=True,
)
@click_job_options
@click.option(
    "--structure",
    "-s",
    "structure",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Initial structure file used by the EM stage.",
)
@click.option(
    "--top",
    "top",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="GROMACS topology .top file shared by all stages.",
)
@click.option(
    "--index",
    "index",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Optional GROMACS index .ndx file.",
)
@click.option(
    "--itp",
    "itp",
    multiple=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    help="Optional GROMACS .itp file. Can be used multiple times.",
)
@click.option(
    "--output-folder",
    "folder",
    type=click.Path(
        file_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default="gromacs_workflow",
    show_default=True,
    help="Directory used for all workflow stage outputs.",
)
@click.option(
    "--em-mdp",
    "em_mdp",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Optional custom EM .mdp file.",
)
@click.option(
    "--nvt-mdp",
    "nvt_mdp",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Optional custom NVT .mdp file.",
)
@click.option(
    "--npt-mdp",
    "npt_mdp",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Optional custom NPT .mdp file.",
)
@click.option(
    "--md-mdp",
    "md_mdp",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Optional custom production MD .mdp file.",
)
@click.pass_context
def workflow(
    ctx,
    skip_completed,
    structure,
    top,
    index,
    itp,
    folder,
    em_mdp,
    nvt_mdp,
    npt_mdp,
    md_mdp,
    molecule=None,
    **kwargs,
):
    """
    Create the sequential EM, NVT, NPT, and MD workflow.
    """

    jobrunner = ctx.obj.get("jobrunner", None)
    project_settings = ctx.obj.get("project_settings", None)
    filename = ctx.obj.get("filename", None)

    if structure is None:
        structure = filename

    mdp_files = {
        "em": em_mdp,
        "nvt": nvt_mdp,
        "npt": npt_mdp,
        "md": md_mdp,
    }

    if project_settings is not None:
        if structure is None:
            structure = project_settings.structure_file

        if top is None:
            top = project_settings.top_file

        if index is None:
            index = project_settings.index_file

        if not itp and project_settings.itp_files:
            itp = tuple(project_settings.itp_files)

        stage = str(project_settings.job_type or "").lower()
        if (
            stage in mdp_files
            and mdp_files[stage] is None
            and project_settings.mdp_file is not None
        ):
            mdp_files[stage] = project_settings.mdp_file

    missing = []

    if structure is None:
        missing.append("--structure")

    if top is None:
        missing.append("--top")

    if missing:
        raise click.UsageError(
            "Missing required workflow options: " + ", ".join(missing)
        )

    logger.info("Creating GROMACS staged workflow: EM -> NVT -> NPT -> MD")
    logger.info("Initial structure file: %s", structure)
    logger.info("Topology file: %s", top)
    logger.info("Workflow output folder: %s", folder)
    logger.info("Skip completed stages: %s", skip_completed)

    if ctx.invoked_subcommand is None:
        return GromacsWorkflow(
            structure_file=Path(structure),
            top_file=Path(top),
            jobrunner=jobrunner,
            folder=Path(folder),
            itp_files=list(itp) if itp else [],
            index_file=Path(index) if index else None,
            mdp_files=mdp_files,
            skip_completed=skip_completed,
        )

    return None
