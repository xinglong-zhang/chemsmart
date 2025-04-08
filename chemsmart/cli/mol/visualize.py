import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import click_pymol_visualization_options, mol
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("visualize", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click.pass_context
def visualize(
    ctx, style_file, quiet, command_line_only, skip_completed, **kwargs
):

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]  # last molecule of the list for visualization
    logger.info(f"Visualizing molecule: {molecule}.")

    # get label for the job
    label = ctx.obj["label"]

    from chemsmart.jobs.mol.visualize import PyMOLVisualizationJob

    return PyMOLVisualizationJob(
        molecule=molecule,
        settings=None,
        label=label,
        pymol_script=style_file,
        quite_mode=quiet,
        command_line_only=command_line_only,
        skip_completed=skip_completed,
        **kwargs,
    )
