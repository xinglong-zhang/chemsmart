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
    ctx,
    style_file,
    render_style,
    vdw,
    quiet,
    command_line_only,
    skip_completed,
    **kwargs,
):
    """CLI for running automatic PyMOL visualization and saving as pse file.
    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone.com visualize -v
    This visualizes phenyldioxazolone.com file and saves as phenyldioxazolone_visualize.pse
    with added Van der Waal's surface (-v) automatically."""

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
        render_style=render_style,
        vdw=vdw,
        quite_mode=quiet,
        command_line_only=command_line_only,
        skip_completed=skip_completed,
        **kwargs,
    )
