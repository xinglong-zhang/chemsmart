import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_hybrid_visualization_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command(
    "visualize",
    cls=MyCommand,
    context_settings=dict(allow_extra_args=True),
)
@click_job_options
@click_pymol_visualization_options
@click_pymol_hybrid_visualization_options
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
    (bonds, angles and dihedrals) for labelling.

    When --hybrid flag is used, the hybrid visualization mode is enabled, which allows the user to draw different groups in different styles.
    Example usage:
    chemsmart run mol -f 'structure_file' visualize -g  '233,468-512' -g '308,397-414,416-423'
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
    from chemsmart.jobs.mol.visualize import (
        PyMOLHybridVisualizationJob,
        PyMOLVisualizationJob,
    )

    if hybrid:
        command_line_only = True

    visualizationjob = (
        PyMOLHybridVisualizationJob if hybrid else PyMOLVisualizationJob
    )

    hybrid_opts = {}

    groups = kwargs.pop("group", ())
    colors = kwargs.pop("color", ())
    # raise error if -g/-c/-sc/-st is provided when --hybrid is false
    hybrid_only_opts = [
        "groups",
        "colors",
        "surface_color",
        "surface_transparency",
    ]
    if any(kwargs.get(opt) for opt in hybrid_only_opts) and not hybrid:
        raise click.UsageError(
            "The options '-g/--group', '--colors', '--surface-color', and "
            "'--surface-transparency' can only be used with '--hybrid'. "
            "Please enable hybrid visualization mode with '--hybrid'."
        )
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
