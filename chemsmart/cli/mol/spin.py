import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_pml_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("spin", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_pml_options
@click.pass_context
def spin(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    isosurface_value,
    transparency_value,
    surface_quality,
    antialias_value,
    ray_trace_mode,
    skip_completed,
    **kwargs,
):
    """Generate spin density visualization and save as PSE file.
    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone.log spin
    This visualizes phenyldioxazolone.log file and saves as
    phenyldioxazolone_spin.pse.
    Requires phenyldioxazolone.chk be present together with
    phenyldioxazolone.log
    """

    # get molecule
    molecules = ctx.obj["molecules"]
    logger.info(f"Visualizing spin density of molecule(s): {molecules}.")

    # get label for the job
    label = ctx.obj["label"]
    if coordinates is not None:
        logger.debug(f"Coordinates for visualization: {coordinates}")
        try:
            coordinates = ast.literal_eval(coordinates)
        except (ValueError, SyntaxError) as e:
            logger.error(
                f"Invalid coordinates input: {coordinates}. Error: {e}"
            )
            raise ValueError(
                "Invalid coordinates input. Please provide a valid Python "
                "literal."
            )
    from chemsmart.jobs.mol.spin import PyMOLSpinJob

    return PyMOLSpinJob(
        molecule=molecules,
        label=label,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        isosurface_value=isosurface_value,
        transparency_value=transparency_value,
        surface_quality=surface_quality,
        antialias_value=antialias_value,
        ray_trace_mode=ray_trace_mode,
        skip_completed=skip_completed,
        **kwargs,
    )
