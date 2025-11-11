import ast
import logging
import re
import sys

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import click_pymol_visualization_options, mol
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.repattern import pymol_hybrid_selection_pattern

logger = logging.getLogger(__name__)


@mol.command(
    "visualize",
    cls=MyCommand,
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
)
@click_job_options
@click_pymol_visualization_options
@click.option(
    "--hybrid",
    is_flag=True,
    default=False,
    help="Use hybrid visualization mode.",
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
    "-g",
    "--group",
    multiple=True,
    type=str,
    help="Indexes of atoms to select for a group. Repeatable for multiple groups, e.g., -g '1-5' -g '6,7,8'.",
)
@click.option(
    "-c",
    "--color",
    multiple=True,
    type=str,
    help="Color for each group. Repeatable to match -g options.",
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

    hybrid_opts = {}
    pattern = re.compile(pymol_hybrid_selection_pattern)
    # zip group indices and color style from command line args
    raw_args = sys.argv
    for i in range(len(raw_args) - 1):
        match = pattern.match(raw_args[i])
        if match:
            kind, idx = match.groups()
            value = raw_args[i + 1]
            key = f"{'group' if kind in ['g', 'group'] else 'color'}{idx}"
            hybrid_opts[key] = value

    groups = kwargs.pop("group", ())
    colors = kwargs.pop("color", ())
    for i, grp in enumerate(groups):
        hybrid_opts[f"group{i + 1}"] = grp
        color_i = colors[i] if i < len(colors) else None
        if color_i:
            hybrid_opts[f"color{i + 1}"] = color_i

    # Include surface options if specified
    if kwargs.get("surface_color"):
        hybrid_opts["surface_color"] = kwargs.pop("surface_color")
    if kwargs.get("surface_transparency"):
        hybrid_opts["surface_transparency"] = kwargs.pop(
            "surface_transparency"
        )

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
