import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import click_pymol_visualization_options, mol
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command("visualize", cls=MyCommand)
@click_job_options
@click_pymol_visualization_options
@click.option(
    "--hybrid",
    is_flag=True,
    default=False,
    help="Use hybrid visualization mode.",
)
@click.option(
    "-g1",
    "--group1",
    type=str,
    default=None,
    help="indexes of atoms to be selected as group 1.",
)
@click.option(
    "-g2",
    "--group2",
    type=str,
    default=None,
    help="indexes of atoms to be selected as group 2.",
)
@click.option(
    "-g3",
    "--group3",
    type=str,
    default=None,
    help="indexes of atoms to be selected as group 3.",
)
@click.option(
    "-g4",
    "--group4",
    type=str,
    default=None,
    help="indexes of atoms to be selected as group 4.",
)
@click.option(
    "-c1",
    "--color1",
    type=str,
    default=None,
    help="customized coloring style of group 1.",
)
@click.option(
    "-c2",
    "--color2",
    type=str,
    default=None,
    help="customized coloring style of group 2.",
)
@click.option(
    "-c3",
    "--color3",
    type=str,
    default=None,
    help="customized coloring style of group 3.",
)
@click.option(
    "-c4",
    "--color4",
    type=str,
    default=None,
    help="customized coloring style of group 4.",
)
@click.option(
    "-sc",
    "--surface-color",
    type=str,
    default=None,
    help="customized surface color.",
)
@click.option(
    "-st",
    "--surface-transparency",
    type=str,
    default=None,
    help="customized surface transparency.",
)
@click.option(
    "--hybrid",
    is_flag=True,
    default=False,
    help="Use hybrid visualization mode.",
)
# all other click options removed here to simplify CLI parsing via **kwargs
@click.pass_context
def visualize(
    ctx,
    file,
    style,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    skip_completed,
    hybrid,
    **kwargs,
):
    """CLI for running automatic PyMOL visualization and saving as pse file.
    Example usage:
        chemsmart run --debug mol -f phenyldioxazolone.com visualize -v
    This visualizes phenyldioxazolone.com file and saves as
    phenyldioxazolone_visualize.pse with added van der Waal's surface (-v)
    automatically.
        chemsmart run --debug mol -f vhr_ox_modred_ts10.log visualize -c
        [[1,2],[3,4,5],[1,3,4,5],[4,5],[4,6,9]]
    This visualizes vhr_ox_modred_ts10.log file and saves as
    vhr_ox_modred_ts10_visualize.pse and add in additional coordinates
    (bonds, angles and dihedrals) for labelling."""

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
    from chemsmart.jobs.mol.visualize import (
        PyMOLHybridVisualizationJob,
        PyMOLVisualizationJob,
    )

    visualizationjob = (
        PyMOLHybridVisualizationJob if hybrid else PyMOLVisualizationJob
    )

    # dynamically extract hybrid options from ctx.params
    hybrid_opts = {}
    if hybrid:
        for key in [
            "group1",
            "group2",
            "group3",
            "group4",
            "color1",
            "color2",
            "color3",
            "color4",
            "surface_color",
            "surface_transparency",
        ]:
            value = ctx.params.get(key)
            if value is not None:
                hybrid_opts[key] = value
            kwargs.pop(key, None)

    job = visualizationjob(
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
        **hybrid_opts,
        **kwargs,
    )

    return job
