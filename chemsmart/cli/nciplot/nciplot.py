import functools
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import (
    return_objects_from_string_index,
)

logger = logging.getLogger(__name__)


def click_nciplot_settings_options(f):
    """Common click options for Gaussian Settings."""

    @click.option(
        "-f",
        "--filenames",
        type=str,
        multiple=True,
        default=None,
        help="filenames from which new NCIPLOT input is prepared.",
    )
    @click.option(
        "-l",
        "--label",
        type=str,
        default=None,
        help="write user input filenames for the job (without extension)",
    )
    @click.option(
        "-a",
        "--append-label",
        type=str,
        default=None,
        help="name to be appended to file for the job",
    )
    @click.option(
        "-r",
        "--rthres",
        type=float,
        default=None,
        help="r distance along a cubic box.",
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
        type=tuple,
        default=None,
        help="(x, y, z) positions around which interactions are represented.",
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
        type=tuple,
        default=None,
        help="Increments along the x, y, z directions of the cube in Å. "
        "The default is set to 0.1, 0.1, 0.1.",
    )
    @click.option(
        "--fragment-label",
        type=int,
        default=None,
        help="ifile label given by the input order.",
    )
    @click.option(
        "--fragment-atoms",
        type=tuple,
        default=None,
        help="Atoms for the fragment.",
    )
    @click.option(
        "-cdd",
        "--cutoff-density-dat",
        type=float,
        default=None,
        help="Cutoff for density used in creating the dat file.",
    )
    @click.option(
        "-cdr",
        "--cutoff-rdg-dat",
        type=float,
        default=None,
        help="Cutoff for RDG (reduced density gradient) used in creating the dat file.",
    )
    @click.option(
        "-cdc",
        "--cutoff-density-cube",
        type=float,
        default=None,
        help="Cutoff for density used in creating the cube file.",
    )
    @click.option(
        "-cdr",
        "--cutoff-rdg-cube",
        type=float,
        default=None,
        help="Cutoff for RDG (reduced density gradient) used in creating the cube file.",
    )
    @click.option(
        "--dgrid-promolecular",
        type=float,
        default=None,
        help="Grids for promolecular densities.",
    )
    @click.option(
        "--integrate/--no-integrate",
        type=bool,
        default=False,
        help="Trigger the integration of properties.",
    )
    @click.option(
        "--number-of-ranges",
        type=int,
        default=None,
        help="Number of ranges to compute properties.",
    )
    @click.option(
        "--ranges",
        type=tuple,
        default=None,
        help="Ranges for computing properties. A lower and upper bounds "
        "are required for every interval (one per line).\n "
        "Number of ranges given should match the number of ranges "
        "specified in --number-of-ranges option.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_nciplot_settings_options
@click.option(
    "-P",
    "--pubchem",
    type=str,
    default=None,
    help="Queries structure from PubChem using name, smiles, cid and conformer information.",
)
@click.pass_context
def gaussian(
    ctx,
    project,
    filename,
    label,
    append_label,
    title,
    charge,
    multiplicity,
    functional,
    basis,
    semiempirical,
    index,
    additional_opt_options,
    additional_route_parameters,
    append_additional_info,
    custom_solvent,
    dieze_tag,
    forces,
    pubchem,
):

    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.settings.gaussian import GaussianProjectSettings

    # get project settings
    project_settings = GaussianProjectSettings.from_project(project)

    # obtain Gaussian Settings from filenames, if supplied; otherwise return defaults

    if filename is None:
        # for cases where filenames is not supplied, eg, get structure from pubchem
        job_settings = GaussianJobSettings.default()
        logger.info(
            f"No filenames is supplied and Gaussian default settings are used:\n{job_settings.__dict__} "
        )
    elif filename.endswith((".com", "gjf", ".inp", ".out", ".log")):
        # filenames supplied - we would want to use the settings from here and do not use any defaults!
        job_settings = GaussianJobSettings.from_filepath(filename)
    # elif filenames.endswith((".xyz", ".pdb", ".mol", ".mol2", ".sdf", ".smi", ".cif", ".traj", ".gro", ".db")):
    else:
        job_settings = GaussianJobSettings.default()
    # else:
    #     raise ValueError(
    #         f"Unrecognised filetype {filenames} to obtain GaussianJobSettings"
    #     )

    # Update keywords
    keywords = (
        "charge",
        "multiplicity",
    )  # default keywords to merge filenames charge and multiplicity
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
    if functional is not None:
        job_settings.functional = functional
        keywords += ("functional",)  # update keywords
    if basis is not None:
        job_settings.basis = basis
        keywords += ("basis",)
    if semiempirical is not None:
        job_settings.semiempirical = semiempirical
        keywords += ("semiempirical",)
    if additional_opt_options is not None:
        job_settings.additional_opt_options_in_route = additional_opt_options
        keywords += ("additional_opt_options_in_route",)
    if additional_route_parameters is not None:
        job_settings.additional_route_parameters = additional_route_parameters
        keywords += ("additional_route_parameters",)
    if append_additional_info is not None:
        job_settings.append_additional_info = append_additional_info
        keywords += ("append_additional_info",)
    if custom_solvent is not None:
        job_settings.custom_solvent = custom_solvent
        keywords += ("custom_solvent",)
    if title is not None:
        job_settings.title = title
        keywords += ("title",)
    if dieze_tag is not None:
        job_settings.dieze_tag = dieze_tag
        keywords += ("dieze_tag",)
    if forces:
        job_settings.forces = forces
        keywords += ("forces",)

    # obtain molecule structure
    molecules = None
    if filename is None and pubchem is None:
        raise ValueError(
            "[filenames] or [pubchem] has not been specified!\nPlease specify one of them!"
        )
    if filename and pubchem:
        raise ValueError(
            "Both [filenames] and [pubchem] have been specified!\nPlease specify only one of them."
        )

    if filename:
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecule from {filename}!"
        logger.debug(
            f"Obtained {len(molecules)} molecule {molecules} from {filename}"
        )

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give Gaussian input filenames or name to be be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        molecules = return_objects_from_string_index(
            list_of_objects=molecules, index=index
        )

    if not isinstance(molecules, list):
        molecules = [molecules]

    logger.debug(f"Obtained molecules: {molecules}")

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label
    ctx.obj["filenames"] = filename


@gaussian.result_callback()
@click.pass_context
def gaussian_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
