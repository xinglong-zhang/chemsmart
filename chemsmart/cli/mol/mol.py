import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
    click_folder_options,
    click_pubchem_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label, select_items_by_index

logger = logging.getLogger(__name__)


def click_pymol_visualization_options(f):
    """Common click options for PyMOL visualization."""

    @click.option(
        "-f",
        "--file",
        type=str,
        default=None,
        help="PyMOL file script or style. If not specified, defaults to "
        "zhang_group_pymol_style.py.",
    )
    @click.option(
        "-s",
        "--style",
        type=click.Choice(["pymol", "cylview"], case_sensitive=False),
        default=None,
        help='PyMOL render style. Choices include "pymol" or "cylview" when '
        "using zhang_group_pymol_style.",
    )
    @click.option(
        "-t/",
        "--trace/--no-trace",
        type=bool,
        default=True,
        help="PyMOL option to ray trace or not. Defaults to True.",
    )
    @click.option(
        "-v",
        "--vdw",
        is_flag=True,
        default=False,
        help="Add Van der Waals surface. Defaults to False.",
    )
    @click.option(
        "-q",
        "--quiet",
        is_flag=True,
        default=False,
        help="Run PyMOL in quiet mode. Defaults to False.",
    )
    @click.option(
        "--command-line-only/--no-command-line-only",
        is_flag=True,
        default=True,
        help="Run PyMOL in command line only mode. Defaults to True.",
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="List of coordinates (bonds, angles, and dihedrals) for "
        "labelling. 1-indexed.",
    )
    @click.option(
        "--label-offset",
        type=str,
        default=None,
        help="Tuple for offsetting label position in mol jobs.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_hybrid_visualization_options(f):
    """Group all hybrid visualization options for reuse."""

    @click.option(
        "-H",
        "--hybrid",
        is_flag=True,
        default=False,
        help="Use hybrid visualization mode.",
    )
    @click.option(
        "-G",
        "--groups",
        multiple=True,
        type=str,
        help=(
            "Indexes of atoms to select for a group. Repeatable for multiple groups, "
            "e.g., -G '1-5' -G '6,7,8'."
        ),
    )
    @click.option(
        "-C",
        "--colors",
        multiple=True,
        type=str,
        help=(
            "Color for each group. Repeatable to match -G options. "
            "Example color scheme: ['cbap', 'cbac', 'cbay', 'cbag', ...]."
        ),
    )
    @click.option(
        "-SC",
        "--surface-color",
        type=str,
        default=None,
        help="Customized surface color.",
    )
    @click.option(
        "-ST",
        "--surface-transparency",
        type=str,
        default=None,
        help="Customized surface transparency.",
    )
    @click.option(
        "-NC",
        "--new-color-carbon",
        type=str,
        default=None,
        help="Color carbon atoms with user-specified color, "
        "e.g. -NC [0.8, 0.8, 0.9]",
    )
    @click.option(
        "-NN",
        "--new-color-nitrogen",
        type=str,
        default=None,
        help="Color nitrogen atoms with user-specified color, "
        "e.g. -NN [0.6, 0.8, 1.0]",
    )
    @click.option(
        "-NO",
        "--new-color-oxygen",
        type=str,
        default=None,
        help="Color oxygen atoms with user-specified color, "
        "e.g. -NO [1.0, 0.7, 0.7]",
    )
    @click.option(
        "-NS",
        "--new-color-sulfur",
        type=str,
        default=None,
        help="Color sulfur atoms with user-specified color, "
        "e.g. -NS [0.8, 0.8, 0.9]",
    )
    @click.option(
        "-NP",
        "--new-color-phosphorus",
        type=str,
        default=None,
        help="Color phosphorus atoms with user-specified color."
        "e.g. -NP [1.0, 0.85, 0.6]",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_nci_options(f):
    """Common click options for PyMOL NCI visualization."""

    @click.option(
        "-i",
        "--isosurface-value",
        type=float,
        default=None,
        help="Isosurface value for NCI plot. Defaults to 0.5.",
    )
    @click.option(
        "-r",
        "--color-range",
        type=float,
        default=1.0,
        help="Ramp range for NCI plot. Defaults to 1.0.",
    )
    @click.option(
        "-b",
        "--binary/--no-binary",
        is_flag=True,
        default=False,
        help="Plot NCI plots with two colors only. Defaults to False.",
    )
    @click.option(
        "--intermediate/--no-intermediate",
        is_flag=True,
        default=False,
        help="Plot NCI plots with intermediate range colors. Defaults to False.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_mo_options(f):
    """Common click options for PyMOL molecular orbital visualization."""

    @click.option(
        "-n",
        "--number",
        type=int,
        default=None,
        help="Molecular orbital number to be visualized (e.g., 31 will "
        "visualize MO #31). Defaults to None.",
    )
    @click.option(
        "-h",
        "--homo",
        is_flag=True,
        default=False,
        help="Plot the highest occupied molecular orbital (HOMO). "
        "Defaults to False.",
    )
    @click.option(
        "-l",
        "--lumo",
        is_flag=True,
        default=False,
        help="Plot the lowest unoccupied molecular orbital (LUMO). "
        "Defaults to False.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_pml_options(f):
    """Common click options for PyMOL .pml files."""

    @click.option(
        "-i",
        "--isosurface-value",
        type=float,
        default=None,
        help="Set isosurface value to be used in PyMOL .pml file.",
    )
    @click.option(
        "-tv",
        "--transparency-value",
        type=float,
        default=None,
        help="Set transparency value to be used in PyMOL .pml file. "
        "Value range: 0.0 - 1.0; 0.0 = fully opaque; 1.0 = fully transparent.",
    )
    @click.option(
        "-sq",
        "--surface-quality",
        type=int,
        default=None,
        help="Set surface quality in PyMOL .pml file. Controls the "
        "quality of molecular surfaces. Higher values yield smoother "
        "surfaces but may increase rendering time. 0 = Low quality "
        "(fast, faceted surfaces); 1 = Medium quality (balanced); "
        "2 = High quality (smooth surfaces, slower rendering); "
        "3 = Very high quality (very smooth, longest rendering time); "
        "4 = Ultra quality (maximum smoothness, may be very slow).",
    )
    @click.option(
        "-a",
        "--antialias-value",
        type=int,
        default=None,
        help="Set antialias value in PyMOL .pml file. Controls smoothing of edges "
        "in the 3D rendering (anti-aliasing). Helps remove jagged edges, "
        "especially useful for high-quality figures. 0 = Off (fast, jagged edges); "
        "1 = On (basic anti-aliasing); 2 = Higher quality anti-aliasing. "
        "Some builds allow up to 4.",
    )
    @click.option(
        "-m",
        "--ray-trace-mode",
        type=int,
        default=None,
        help="Set ray trace mode in PyMOL .pml file. Controls quality "
        "of ray-traced images. Higher values yield better quality "
        "but take longer to render. 0 = Normal shading (standard "
        "photorealistic render); 1 = Cartoon/line outlines (black "
        "outlines around objects, like cell-shading); 2 = Black "
        "outline only (no shading, wireframe-like appearance); "
        "3 = White outline mode (for figures on dark backgrounds).",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_pymol_save_options(f):
    """Common click options for PyMOL save options."""

    @click.option(
        "-o",
        "--overwrite",
        is_flag=True,
        default=False,
        help="Overwrite existing files. Defaults to False.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_filenames_options
@click_file_label_and_index_options
@click_folder_options
@click_pubchem_options
@click.pass_context
def mol(
    ctx,
    filenames,
    label,
    append_label,
    index,
    directory,
    filetype,
    pubchem,
):
    """
    CLI subcommand for running PyMOL visualization jobs using the chemsmart framework.

    Example usage:
        chemsmart run mol -f test.xyz visualize -c [[413,409],[413,412],[413,505],[413,507]]
    """
    # Initialize molecules variable
    molecules = None

    # obtain molecule structure
    if directory is not None and filetype is not None:
        ctx.obj["directory"] = directory
        ctx.obj["filetype"] = filetype
        ctx.obj["index"] = index
        ctx.obj["filenames"] = None
        ctx.obj["molecules"] = None
        ctx.obj["label"] = label
        return

    if filenames is None and pubchem is None:
        # this is fine for PyMOL IRC Movie Job and Align job
        logger.warning("[filename] or [pubchem] has not been specified!")
        ctx.obj["molecules"] = None
        ctx.obj["label"] = None
        return
    # if both filename and pubchem are specified, raise error
    if filenames and pubchem:
        raise ValueError(
            "Both [filename] and [pubchem] have been specified!\nPlease specify only one of them."
        )

    # if filename is specified, read the file and obtain molecule
    if filenames:
        # Check if this is an align task by looking for " align" in invokded subcommand
        is_align_task = ctx.invoked_subcommand == "align"

        if is_align_task:
            # align task can handle multiple files - pass to align command
            ctx.obj["filenames"] = filenames
            ctx.obj["index"] = index
            ctx.obj["directory"] = None
            ctx.obj["filetype"] = None
            ctx.obj["molecules"] = None
            ctx.obj["label"] = label
            return
        else:
            if len(filenames) == 1:
                filenames = filenames[0]
                molecules = Molecule.from_filepath(
                    filepath=filenames, index=":", return_list=True
                )
                assert (
                    molecules is not None
                ), f"Could not obtain molecule from {filenames}!"
                logger.debug(f"Obtained molecule {molecules} from {filenames}")
            else:
                raise ValueError(
                    f"This task can only process one file, but {len(filenames)} files were provided. "
                )

    # if pubchem is specified, obtain molecule from PubChem
    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give Gaussian input filename or name to be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filenames))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filenames))[0]

    label = clean_label(label)

    logger.info(
        f"Obtained molecules: {molecules} before applying indices,"
        f"with label: {label}"
    )

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        logger.debug(f"Using molecule with index: {index}")
        molecules = select_items_by_index(
            molecules,
            index,
            allow_duplicates=False,
            allow_out_of_range=False,
        )
    else:
        molecules = [molecules[-1]]  # Default: last molecule as list

    logger.debug(f"Obtained molecules: {molecules}")

    # store objects
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label


@mol.result_callback()
@click.pass_context
def mol_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
