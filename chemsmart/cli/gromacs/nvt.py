from __future__ import annotations

"""
GROMACS NVT equilibration CLI command.

This command supports two creation modes:
1. YAML-driven mode: chemsmart run gromacs -p project.yaml nvt
2. Direct CLI mode: chemsmart run gromacs nvt --mdp nvt.mdp --structure em.gro --top topol.top
"""

import logging
from pathlib import Path

import click

from chemsmart.cli.gromacs.gromacs import gromacs
from chemsmart.cli.job import click_job_options
from chemsmart.jobs.gromacs.job import GromacsNVTJob
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


@gromacs.group("nvt", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click.option(
    "--mdp",
    "mdp_file",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="GROMACS .mdp file for NVT equilibration.",
)
@click.option(
    "--structure",
    "-s",
    "structure_file",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Input structure file for NVT equilibration, such as EM output .gro.",
)
@click.option(
    "--input-pdb",
    "input_pdb",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="Raw PDB file used by the full_setup workflow.",
)
@click.option(
    "--top",
    "top_file",
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    default=None,
    help="GROMACS topology .top file.",
)
@click.option(
    "--index",
    "index_file",
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
    "itp_files",
    multiple=True,
    type=click.Path(
        exists=True,
        dir_okay=False,
        resolve_path=True,
        path_type=Path,
    ),
    help="Optional GROMACS .itp include file. Can be used multiple times.",
)
@click.option(
    "--workflow",
    type=click.Choice(["prepared", "full_setup"]),
    default=None,
    help="GROMACS workflow mode.",
)
@click.pass_context
def nvt(
    ctx,
    skip_completed,
    mdp_file,
    structure_file,
    input_pdb,
    top_file,
    index_file,
    itp_files,
    workflow,
    molecule=None,
    **kwargs,
):
    """
    CLI subcommand for running GROMACS NVT equilibration.
    """

    jobrunner = ctx.obj.get("jobrunner", None)
    project = ctx.obj.get("project", None)
    project_yaml = ctx.obj.get("project_yaml", None)
    project_settings = ctx.obj.get("project_settings", None)
    filename = ctx.obj.get("filename", None)

    if structure_file is None:
        structure_file = filename

    if project_settings is not None:
        settings = project_settings.with_overrides(
            mdp_file=Path(mdp_file) if mdp_file else None,
            structure_file=Path(structure_file) if structure_file else None,
            input_pdb=Path(input_pdb) if input_pdb else None,
            top_file=Path(top_file) if top_file else None,
            index_file=Path(index_file) if index_file else None,
            workflow=workflow,
            itp_files=list(itp_files) if itp_files else None,
        )

        settings.validate()

        logger.info(
            "Creating GROMACS NVT job from project YAML: %s",
            project_yaml,
        )

        if ctx.invoked_subcommand is None:
            return GromacsNVTJob.from_project_settings(
                settings=settings,
                molecule=molecule,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
        return None

    # direct CLI mode
    label = "gromacs_nvt"
    if structure_file is not None:
        label = Path(structure_file).stem + "_nvt"

    workflow = workflow or "prepared"

    logger.info("Creating GROMACS NVT job from direct CLI options.")
    logger.info("GROMACS NVT structure file: %s", structure_file)
    logger.info("GROMACS NVT mdp file: %s", mdp_file)
    logger.info("GROMACS NVT topology file: %s", top_file)
    logger.info("GROMACS NVT itp files: %s", itp_files)
    logger.info("GROMACS NVT workflow: %s", workflow)

    if ctx.invoked_subcommand is None:
        return GromacsNVTJob(
            molecule=molecule,
            label=label,
            jobrunner=jobrunner,
            mdp_file=Path(mdp_file) if mdp_file else None,
            structure_file=Path(structure_file) if structure_file else None,
            input_pdb=Path(input_pdb) if input_pdb else None,
            top_file=Path(top_file) if top_file else None,
            itp_files=list(itp_files) if itp_files else [],
            index_file=Path(index_file) if index_file else None,
            workflow=workflow,
            skip_completed=skip_completed,
            **kwargs,
        )

