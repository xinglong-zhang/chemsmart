import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("align", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click.pass_context
def align(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    skip_completed,
    **kwargs,
):
    """CLI for PyMOL alignment of multiple molecule files.

    Example:
        chemsmart run mol -f a.log -f b.xyz -f c.gjf align
        chemsmart run mol -af '*.log(*.xyz/*.gjf)' align
    """

    # get molecule
    molecules = ctx.obj["molecules"]
    if not isinstance(molecules, list) or len(molecules) < 2:
        raise click.BadParameter("Need at least two molecules")

    # get label for the job
    label = ctx.obj["label"]
    if coordinates is not None:
        logger.debug(f"Coordinates for alignment: {coordinates}")
        try:
            coordinates = ast.literal_eval(coordinates)
        except (ValueError, SyntaxError) as e:
            logger.error(
                f"Invalid coordinates input: {coordinates}. Error: {e}"
            )
            raise ValueError(
                "Invalid coordinates input. Please provide a valid Python literal."
            )

    from chemsmart.jobs.mol.align import PyMOLAlignJob

    return PyMOLAlignJob(
        molecule=molecules,
        label=label,
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
