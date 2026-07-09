import ast
import logging

import click
from click.core import ParameterSource

from chemsmart.cli.job import click_job_options
from chemsmart.cli.mol.mol import (
    click_pymol_hybrid_visualization_options,
    click_pymol_visualization_options,
    mol,
)
from chemsmart.jobs.mol.runner import normalize_pymol_style
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@mol.command(
    "visualize",
    cls=MyCommand,
    context_settings=dict(allow_extra_args=True),
)
@click_job_options
@click_pymol_visualization_options(include_visualize_styles=True)
@click_pymol_hybrid_visualization_options
@click.pass_context
def visualize(
    ctx,
    file,
    style,
    style_background,
    trace,
    vdw,
    quiet,
    command_line_only,
    coordinates,
    skip_completed,
    hybrid,
    **kwargs,
):
    """CLI subcommand for running automatic
    PyMOL visualization and save as PSE file.
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

    When --hybrid flag is used, the hybrid visualization mode is enabled,
    which allows the user to draw different groups in different styles.
    Example usage:
    chemsmart run mol -f 'structure_file' visualize
    -G '233,468-512' -G '308,397-414,416-423'

    When -s glossy is used, the glossy semi-metallic visualization mode is
    enabled. Example usage:
        chemsmart run mol -f complex.xyz visualize -s glossy --style-background dark

    When -s comic is used, the comic illustration mode is enabled. Example usage:
        chemsmart run mol -f complex.xyz visualize -s comic --style-background dark
    """

    molecules = ctx.obj["molecules"]
    logger.info(f"Visualizing molecule(s): {molecules}.")

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

    normalized_style = None
    if style is not None:
        normalized_style = normalize_pymol_style(style)

    is_glossy = normalized_style == "glossy"
    is_comic = normalized_style == "comic"
    is_special_style = is_glossy or is_comic

    style_background_explicit = (
        ctx.get_parameter_source("style_background") != ParameterSource.DEFAULT
    )

    if hybrid and is_special_style:
        raise click.UsageError(
            "Hybrid visualization (-H/--hybrid) cannot be combined with "
            f"-s {style}."
        )

    if style_background_explicit and not is_special_style:
        raise click.UsageError(
            "'--style-background' can only be used with '-s glossy', "
            "'-s comic', or '-s hybrid' (alias for comic)."
        )

    from chemsmart.jobs.mol.visualize import (
        PyMOLComicVisualizationJob,
        PyMOLGlossyVisualizationJob,
        PyMOLHybridVisualizationJob,
        PyMOLVisualizationJob,
    )

    if hybrid:
        visualization_job = PyMOLHybridVisualizationJob
    elif is_glossy:
        visualization_job = PyMOLGlossyVisualizationJob
    elif is_comic:
        visualization_job = PyMOLComicVisualizationJob
    else:
        visualization_job = PyMOLVisualizationJob

    hybrid_opts = {}

    groups = kwargs.pop("groups", ())
    hybrid_opts["groups"] = groups
    colors = kwargs.pop("colors", ())
    hybrid_opts["colors"] = colors

    hybrid_only_opts = [
        "groups",
        "colors",
        "surface_color",
        "surface_transparency",
        "new_color_carbon",
        "new_color_nitrogen",
        "new_color_oxygen",
        "new_color_sulfur",
        "new_color_phosphorus",
    ]
    if any(kwargs.get(opt) for opt in hybrid_only_opts) and not hybrid:
        raise click.UsageError(
            "The options '-G/--groups', '-C/--colors', '--surface-color', "
            "'--surface-transparency', and '--new-color-*' can only be used "
            "with '-H/--hybrid'. Please enable hybrid visualization mode "
            "with '-H/--hybrid'."
        )

    for hybrid_opt in [
        "surface_color",
        "surface_transparency",
        "new_color_carbon",
        "new_color_nitrogen",
        "new_color_oxygen",
        "new_color_sulfur",
        "new_color_phosphorus",
    ]:
        if kwargs.get(hybrid_opt):
            hybrid_opts[hybrid_opt] = kwargs.pop(hybrid_opt)

    logger.info(f"Hybrid visualization job options: {hybrid_opts}")

    job_kwargs = dict(
        molecule=molecules,
        label=label,
        pymol_script=file,
        trace=trace,
        vdw=vdw,
        quiet_mode=quiet,
        command_line_only=command_line_only,
        coordinates=coordinates,
        skip_completed=skip_completed,
        **hybrid_opts,
        **kwargs,
    )

    if is_special_style:
        job_kwargs["style_background"] = style_background
    else:
        job_kwargs["style"] = style

    job = visualization_job(**job_kwargs)

    return job
