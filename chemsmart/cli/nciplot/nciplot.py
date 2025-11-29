import ast
import functools
import logging

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filenames_options,
    click_job_options,
    click_pubchem_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.nciplot.job import NCIPLOTJob
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.io import clean_label

logger = logging.getLogger(__name__)


def click_nciplot_settings_options(f):
    """Common click options for NCIPLOT Job Settings."""

    @click.option(
        "-r",
        "--rthres",
        type=float,
        default=None,
        help="r distance along a cubic box. This extends the grid to include "
        "a larger region around the molecule for capturing NCIs that extend "
        "further out.",
    )
    @click.option(
        "--ligand-file-number",
        type=int,
        default=None,
        help="Ligand file number, corresponds to the nth file of the input.",
    )
    @click.option(
        "--ligand-radius",
        type=float,
        default=None,
        help="Radius of interaction from ligand file n.",
    )
    @click.option(
        "-rp",
        "--radius-positions",
        type=str,
        default=None,
        help="Positions around which interactions are represented.\n"
        "Accepts strings in the form of 'x,y,z' or as a tuple (x, y, z). ",
    )
    @click.option(
        "-rr",
        "--radius-r",
        type=float,
        default=None,
        help="Radius from which interactions are represented.",
    )
    @click.option(
        "-i1",
        "--intercut1",
        type=float,
        default=None,
        help="Cutoff 1, r1, for intermolecularity.",
    )
    @click.option(
        "-i2",
        "--intercut2",
        type=float,
        default=None,
        help="Cutoff 2, r2, for intermolecularity.",
    )
    @click.option(
        "--increments",
        type=str,
        default=None,
        help="Increments along the x, y, z directions of the cube in Å. "
        "The default is set to 0.1, 0.1, 0.1."
        "Accepts strings in the form of 'x,y,z' or as a tuple (x, y, z).",
    )
    @click.option(
        "--fragments",
        type=str,
        default=None,
        help="Fragments in a dictionary-type string.\n E.g., "
        "{1: [1, 2, 3], 2: [4, 5, 6]}",
    )
    @click.option(
        "-cdd",
        "--cutoff-density-dat",
        type=float,
        default=None,
        help="Cutoff for density used in creating the dat file.",
    )
    @click.option(
        "-crd",
        "--cutoff-rdg-dat",
        type=float,
        default=None,
        help="Cutoff for RDG (reduced density gradient) used in creating the "
        "dat file.",
    )
    @click.option(
        "-cdc",
        "--cutoff-density-cube",
        type=float,
        default=None,
        help="Cutoff for density used in creating the cube file.",
    )
    @click.option(
        "-crc",
        "--cutoff-rdg-cube",
        type=float,
        default=None,
        help="Cutoff for RDG (reduced density gradient) used in creating the "
        "cube file.",
    )
    @click.option(
        "--dgrid/--no-dgrid",
        type=bool,
        default=False,
        help="Turn on Radial grids for promolecular densities. Default is "
        "exponential fits.\nExponential fits are available up to Z=18 (Ar); "
        "radial grids are available up to Pu (Z=94).\nDefaults to exponential "
        "fits unless the molecule contains some atom with Z>18 or there are "
        "charged atoms (cations). Using DGRID, only radial grids are used to "
        "calculate promolecular densities. ",
    )
    @click.option(
        "--integrate/--no-integrate",
        type=bool,
        default=False,
        help="Trigger the integration of properties.",
    )
    @click.option(
        "--ranges",
        default=None,
        help="Ranges for computing properties. A lower and upper bounds "
        "are required for every interval (one per line).\n "
        "Example input: [[-0.1,-0.02],[-0.02,0.02],[0.02,0.1]]\n",
    )
    @click.option(
        "--grid-quality",
        type=click.Choice(
            ["coarse", "fine", "ultrafine"], case_sensitive=False
        ),
        default=None,
        help="Quality of the grid used for NCIPLOT calculations.\n "
        "Options are 'coarse', 'fine', or 'ultrafine'. ",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.command(cls=MyCommand)
@click_job_options
@click_filenames_options
@click_file_label_and_index_options
@click_nciplot_settings_options
@click_pubchem_options
@click.pass_context
def nciplot(
    ctx,
    filenames,
    label,
    rthres,
    ligand_file_number,
    ligand_radius,
    radius_positions,
    radius_r,
    intercut1,
    intercut2,
    increments,
    fragments,
    cutoff_density_dat,
    cutoff_rdg_dat,
    cutoff_density_cube,
    cutoff_rdg_cube,
    dgrid,
    integrate,
    ranges,
    grid_quality,
    pubchem,
    **kwargs,
):
    """CLI for running NCIPLOT jobs using the chemsmart framework.
    Example usage:
    chemsmart run nciplot -f test.xyz -f test2.xyz -l nci_test \\
    --fragments "{1: [1,4,5], 2: [3,4,5]}"
    """

    from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings

    # update NCIPLOT job settings using values from CLI
    if len(filenames) == 0:
        filenames = None

    if fragments is not None:
        try:
            fragments = ast.literal_eval(fragments)
        except (ValueError, SyntaxError) as e:
            raise click.BadParameter(f"Invalid format for --fragments: {e}")

    if ranges is not None:
        try:
            ranges = ast.literal_eval(ranges)
        except (ValueError, SyntaxError) as e:
            raise click.BadParameter(f"Invalid format for --ranges: {e}")

    job_settings = NCIPLOTJobSettings(
        rthres=rthres,
        ligand_file_number=ligand_file_number,
        ligand_radius=ligand_radius,
        radius_positions=radius_positions,
        radius_r=radius_r,
        intercut1=intercut1,
        intercut2=intercut2,
        increments=increments,
        fragments=fragments,
        cutoff_density_dat=cutoff_density_dat,
        cutoff_rdg_dat=cutoff_rdg_dat,
        cutoff_density_cube=cutoff_density_cube,
        cutoff_rdg_cube=cutoff_rdg_cube,
        dgrid=dgrid,
        integrate=integrate,
        ranges=ranges,
        grid_quality=grid_quality,
    )
    logger.info(f"Obtained NCIPLOT job settings: {job_settings.__dict__}")

    molecule = None
    if pubchem:
        molecule = Molecule.from_pubchem(identifier=pubchem, return_list=False)
        assert (
            molecule is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecule} from PubChem {pubchem}")

        assert (
            filenames is None
        ), "Cannot provide both filenames and PubChem identifier!"
        assert (
            label is not None
        ), "Label for file is required since creating molecule from PubChem!"
        if not label.endswith("promolecular"):
            label = f"{label}_promolecular"
    else:
        if filenames is None:
            raise ValueError(
                "Must provide at least one input file using -f option."
            )
        else:
            if not isinstance(filenames, (list, tuple)):
                raise TypeError(
                    f"Expected filenames to be a list or tuple, got "
                    f"{type(filenames).__name__}"
                )
            if len(filenames) == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide "
                    "at least one file."
                )
            elif len(filenames) == 1:
                label = filenames[0].split(".")[0] if label is None else label
                if not filenames[0].endswith((".wfn", ".wfx")):
                    if label is not None and not label.endswith(
                        "promolecular"
                    ):
                        label = f"{label}_promolecular"
            else:
                # add filenames together
                label = (
                    "_and_".join(
                        [filename.split(".")[0] for filename in filenames]
                    )
                    if label is None
                    else label
                )
    label = clean_label(label)

    return NCIPLOTJob(
        filenames=filenames,  # accepts multiple files
        molecule=molecule,
        settings=job_settings,
        label=label,
        jobrunner=None,
        **kwargs,
    )
