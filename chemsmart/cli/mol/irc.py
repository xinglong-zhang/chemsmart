import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_save_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("irc", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click_pymol_save_options
@click.option(
    "-r",
    "--reactant",
    type=str,
    required=True,
    help="IRC file leading to the reactant side.",
)
@click.option(
    "-p",
    "--product",
    type=str,
    required=True,
    help="IRC file leading to the product side.",
)
@click.pass_context
def irc(
    ctx,
    reactant,
    product,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    overwrite,
    skip_completed,
    **kwargs,
):
    """CLI for running automatic PyMOL visualization and saving as pse file.
    Example usage:
        chemsmart run --debug mol irc -r vhr_ox_modred_ts10_ircr.log -p vhr_ox_modred_ts10_ircf.log -c [1,12] -o
    This makes an IRC movie from vhr_ox_modred_ts10_ircr.log and vhr_ox_modred_ts10_ircf.log, with coordinate labels.
    If the movie mp4 file exists, it will not be overwritten unless -o is specified.
    """

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

    from chemsmart.jobs.mol.irc import PyMOLIRCMovieJob

    return PyMOLIRCMovieJob.from_files(
        reactant_file=reactant,
        product_file=product,
        label=label,
        pymol_script=file,
        style=style,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        overwrite=overwrite,
        skip_completed=skip_completed,
        **kwargs,
    )
