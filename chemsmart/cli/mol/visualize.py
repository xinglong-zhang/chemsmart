import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import mol
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("visualize", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization.",
)
@click.pass_context
def visualize(ctx, freeze_atoms, skip_completed, **kwargs):

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
        skip_completed=skip_completed,
        **kwargs,
    )
