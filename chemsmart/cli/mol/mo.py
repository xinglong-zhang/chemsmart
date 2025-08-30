import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (click_pymol_mo_options,
                                   click_pymol_visualization_options, mol)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("mo", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_mo_options
@click.pass_context
def mo(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    number,
    homo,
    lumo,
    skip_completed,
    **kwargs,
):
    """CLI for generating molecular orbitals (MOs) and saving as pse file.
    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone.com mo --homo
    This visualizes the HOMO of phenyldioxazolone.com file and saves as phenyldioxazolone_HOMO.pse
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
                "Invalid coordinates input. Please provide a valid Python literal."
            )
    from chemsmart.jobs.mol.mo import PyMOLMOJob

    return PyMOLMOJob(
        molecule=molecules,
        label=label,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        number=number,
        homo=homo,
        lumo=lumo,
        skip_completed=skip_completed,
        **kwargs,
    )
