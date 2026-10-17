import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import click_pymol_visualization_options, mol
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("nbo", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click.option(
    "--threshold",
    type=float,
    default=10.0,
    show_default=True,
    help="Minimum E(2) value (kcal/mol) to include in NBO visualization.",
)
@click.option(
    "-N",
    "--max-interactions",
    type=int,
    default=20,
    show_default=True,
    help="Maximum number of NBO perturbation interactions to visualize.",
)
@click.pass_context
def nbo(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    threshold,
    max_interactions,
    skip_completed,
    **kwargs,
):
    """CLI subcommand for visualizing second-order NBO perturbation analysis."""

    molecules = ctx.obj["molecules"]
    label = ctx.obj["label"]
    source_basename = ctx.obj["source_basename"]
    nbo_basename = label if ctx.obj["label_provided"] else None
    filenames = ctx.obj["filenames"]

    if isinstance(filenames, (list, tuple)):
        if len(filenames) != 1:
            raise ValueError(
                "NBO visualization requires a single Gaussian output file."
            )
        analysis_filename = filenames[0]
    else:
        analysis_filename = filenames

    if analysis_filename is None:
        raise ValueError(
            "Could not determine Gaussian output file for NBO visualization."
        )

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

    from chemsmart.jobs.mol.nbo import PyMOLNBOJob

    return PyMOLNBOJob(
        molecule=molecules,
        label=label,
        source_basename=source_basename,
        analysis_filename=analysis_filename,
        threshold=threshold,
        max_interactions=max_interactions,
        nbo_basename=nbo_basename,
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
