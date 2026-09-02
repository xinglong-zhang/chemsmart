import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_esp_options,
    click_pymol_pml_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("esp", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_pml_options
@click_pymol_esp_options
@click.pass_context
def esp(
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
    npts,
    color_range,
    skip_completed,
    **kwargs,
):
    """Generate electrostatic potential surface visualization and save as PSE.

    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone.log esp

    Requires the corresponding .chk file together with the Gaussian output.
    """
    molecules = ctx.obj["molecules"]
    logger.info(f"Visualizing ESP of molecule(s): {molecules}.")

    label = ctx.obj["label"]
    source_basename = ctx.obj["source_basename"]
    esp_basename = label if ctx.obj["label_provided"] else None
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
    from chemsmart.jobs.mol.esp import PyMOLESPJob

    return PyMOLESPJob(
        molecule=molecules,
        label=label,
        esp_basename=esp_basename,
        color_range=color_range,
        npts=npts,
        isosurface_value=isosurface_value,
        source_basename=source_basename,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        transparency_value=transparency_value,
        surface_quality=surface_quality,
        antialias_value=antialias_value,
        ray_trace_mode=ray_trace_mode,
        skip_completed=skip_completed,
        **kwargs,
    )
