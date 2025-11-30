import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_nci_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("nci", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_nci_options
@click.pass_context
def nci(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    isosurface_value,
    color_range,
    binary,
    intermediate,
    skip_completed,
    **kwargs,
):
    """Generate automatic PyMOL NCI plot and save as PSE file.
    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone_nci.log nci
    This visualizes phenyldioxazolone_nci.log file and saves as
    phenyldioxazolone_nci.pse.
    Note that the phenyldioxazolone_nci-dens.cube and phenyldioxazolone_nci-grad.cube
    files should be present for the plotting of nci to work; otherwise,
    an error will be raised.
    """

    # get molecule
    molecules = ctx.obj["molecules"]
    logger.info(f"Visualizing molecule(s): {molecules}.")

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
    from chemsmart.jobs.mol.nci import PyMOLNCIJob

    return PyMOLNCIJob(
        molecule=molecules,
        label=label,
        isosurface_value=isosurface_value,
        color_range=color_range,
        binary=binary,
        intermediate=intermediate,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        skip_completed=skip_completed,
        **kwargs,
    )
