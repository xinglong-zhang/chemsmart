"""
xTB Command Line Interface

This module provides the main CLI interface for xTB semiempirical quantum chemistry
calculations. It defines common options, settings configurations, and the
main xTB command group that serves as the entry point for all xTB-related
operations.
"""

import functools
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def click_xtb_options(f):
    """
    Common click options decorator for xTB jobs.

    This decorator adds common command-line options that are shared across
    different xTB job types, specifically project settings.
    """

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_xtb_settings_options(f):
    """
    Common click options decorator for xTB computational settings.

    This decorator adds comprehensive command-line options for configuring
    xTB calculations including file I/O, molecular properties, GFN methods,
    optimization levels, and solvation parameters specific to xTB.
    """

    @click.option(
        "-f",
        "--filename",
        type=str,
        default=None,
        help="filename from which new xTB input is prepared.",
    )
    @click.option(
        "-l",
        "--label",
        type=str,
        default=None,
        help="write user input filename for the job (without extension)",
    )
    @click.option(
        "-a",
        "--append-label",
        type=str,
        default=None,
        help="name to be appended to file for the job",
    )
    @click.option(
        "-t", "--title", type=str, default=None, help="xTB job title."
    )
    @click.option(
        "-c", "--charge", type=int, default=None, help="charge of the molecule"
    )
    @click.option(
        "-m",
        "--multiplicity",
        type=int,
        default=None,
        help="multiplicity of the molecule",
    )
    @click.option(
        "-g",
        "--gfn-version",
        type=click.Choice(
            ["gfn0", "gfn1", "gfn2", "gfnff"], case_sensitive=False
        ),
        default=None,
        help="GFN-xTB method version (gfn0, gfn1, gfn2, gfnff)",
    )
    @click.option(
        "-i",
        "--index",
        type=str,
        default=None,
        help="Index of molecules to use; 1-based indices. "
        "Default to the last molecule structure. 1-based index.",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="additional route parameters",
    )
    @click.option(
        "--grad/--no-grad",
        default=False,
        help="Whether to calculate gradient.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_xtb_jobtype_options(f):
    """Common click options for xTB link/crest jobs."""

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='xTB job type. Options: ["opt", "sp"]',
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="List of coordinates to be fixed for modred or scan job. 1-indexed.",
    )
    @click.option(
        "-s",
        "--step-size",
        default=None,
        help="Step size of coordinates to scan.",
    )
    @click.option(
        "-n",
        "--num-steps",
        default=None,
        help="Step size of coordinates to scan.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_xtb_solvent_options(f):
    """Common click solvent options for xTB jobs."""

    @click.option(
        "--remove-solvent/--no-remove-solvent",
        "-r/ ",
        type=bool,
        default=False,
        help="Whether to use solvent model in the job. Defaults to project settings.",
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="solvent model to be used for single point.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="solvent ID to be used for single point.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_xtb_options
@click_xtb_settings_options
@click.option(
    "-P",
    "--pubchem",
    type=str,
    default=None,
    help="Queries structure from PubChem using name, smiles, cid and conformer information.",
)
@click.pass_context
def xtb(
    ctx,
    project,
    filename,
    label,
    append_label,
    title,
    charge,
    multiplicity,
    gfn_version,
    index,
    additional_route_parameters,
    grad,
    pubchem,
):
    """CLI subcommand for running xTB semiempirical quantum chemistry jobs."""

    from chemsmart.jobs.xtb.settings import XTBJobSettings
    from chemsmart.settings.xtb import XTBProjectSettings

    # get project settings
    project_settings = XTBProjectSettings.from_project(project)

    # obtain xTB Settings from filename, if supplied; otherwise return defaults

    if filename is None:
        # for cases where filename is not supplied, eg, get structure from pubchem
        job_settings = XTBJobSettings.default()
        logger.info(
            f"No filename is supplied and xTB default settings are used:\n{job_settings.__dict__} "
        )
    elif filename.endswith((".com", ".gjf", ".inp", ".out", ".log")):
        # Read from Gaussian/ORCA files but only extract basic molecular info
        job_settings = XTBJobSettings.from_filepath(filename)
        logger.info(
            f"Extracted molecular information from {filename} for xTB calculation"
        )
    else:
        # all other input files such as .xyz, .pdb, .mol, .mol2, .sdf, .smi, .cif, .traj, .gro, .db
        job_settings = XTBJobSettings.default()

    # Update keywords
    keywords = (
        "charge",
        "multiplicity",
    )  # default keywords to merge filename charge and multiplicity
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
    if gfn_version is not None:
        job_settings.gfn_version = gfn_version.lower()
        keywords += ("gfn_version",)
    if additional_route_parameters is not None:
        job_settings.input_string = additional_route_parameters
        keywords += ("input_string",)
    if title is not None:
        job_settings.title = title
        keywords += ("title",)
    if grad is not None:
        job_settings.grad = grad
        keywords += ("grad",)

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give xTB input filename or name to be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = (
            f"{label}_{ctx.invoked_subcommand}"  # add in subcommand to label
        )

    # obtain molecule structure
    if filename is None and pubchem is None:
        raise ValueError(
            "[filename] or [pubchem] has not been specified!\nPlease specify one of them!"
        )
    if filename and pubchem:
        raise ValueError(
            "Both [filename] and [pubchem] have been specified!\nPlease specify only one of them."
        )

    if filename:
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        assert (
            molecules is not None
        ), f"Could not obtain molecule from {filename}!"
        logger.debug(f"Obtained molecule {molecules} from {filename}")

        # # create input file as .xyz for xTB jobs
        # if filename.endswith(".xyz"):
        #     # if xyz file, then copy this to label.xyz
        #     from shutil import copy
        #
        #     copy(filename, f"{label}.xyz")
        # else:
        #     # convert file to .xyz
        #     molecule = Molecule.from_filepath(
        #         filepath=filename, index="-1", return_list=False
        #     )
        #     logger.debug(f"Writting molecule {molecule} to: {label}.xyz")
        #     molecule.write_xyz(filename=f"{label}.xyz")

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    logger.debug(f"Obtained molecules: {molecules} before applying indices")

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        try:
            # try to get molecule using python style string indexing,
            # but in 1-based
            from chemsmart.utils.utils import string2index_1based

            index = string2index_1based(index)
            molecules = molecules[index]
            if not isinstance(molecules, list):
                molecules = [molecules]
        except ValueError:
            # except user defined indices such as s='[1-3,28-31,34-41]'
            # or s='1-3,28-31,34-41' which cannot be parsed by string2index_1based
            index = get_list_from_string_range(index)
            molecules = [molecules[i - 1] for i in index]

    logger.debug(f"Obtained molecules: {molecules}")

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@xtb.result_callback()
@click.pass_context
def xtb_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
