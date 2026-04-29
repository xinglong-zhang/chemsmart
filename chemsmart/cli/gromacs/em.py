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
@click.pass_context
def em(ctx, skip_completed, **kwargs):
    """CLI subcommand for running GROMACS energy minimization."""

    jobrunner = ctx.obj.get("jobrunner", None)
    project = ctx.obj.get("project", None)
    filename = ctx.obj.get("filename", None)

    logger.info(f"GROMACS EM project: {project}")
    logger.info(f"GROMACS EM filename: {filename}")

    label = "gromacs_em"
    if filename is not None:
        import os
        label = os.path.splitext(os.path.basename(filename))[0] + "_em"

    if ctx.invoked_subcommand is None:
        return GromacsEMJob(
            molecule=None,
            settings=None,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )