import functools
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def click_gaussian_options(f):
    """Common click options for Gaussian jobs."""

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_settings_options(f):
    """Common click options for Gaussian Settings."""

    @click.option(
        "-f",
        "--filename",
        type=str,
        default=None,
        help="filename from which new Gaussian input is prepared.",
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
        "-t", "--title", type=str, default=None, help="Gaussian job title."
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
        "-x",
        "--functional",
        type=str,
        default=None,
        help="New functional to run.",
    )
    @click.option(
        "-b", "--basis", type=str, default=None, help="New basis set to run."
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
        "-o",
        "--additional-opt-options",
        type=str,
        default=None,
        help="additional opt options",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="additional route parameters",
    )
    @click.option(
        "-A",
        "--append-additional-info",
        type=str,
        default=None,
        help="additional information to be appended at the end of the "
        "input file. E.g, scrf=read",
    )
    @click.option(
        "-C",
        "--custom-solvent",
        type=str,
        default=None,
        help="additional information to be appended at the end of the "
        "input file. E.g, scrf=read",
    )
    @click.option(
        "-d",
        "--dieze-tag",
        type=str,
        default=None,
        help="dieze tag for gaussian job; possible options include "
        '"n", "p", "t" to get "#n", "#p", "#t", respectively',
    )
    @click.option(
        "--forces/--no-forces",
        default=False,
        help="Whether to calculate forces.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_jobtype_options(f):
    """Common click options for Gaussian link/crest jobs."""

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp"]',
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


def click_gaussian_solvent_options(f):
    """Common click solvent options for Gaussian jobs."""

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


def click_gaussian_td_options(f):
    @click.option(
        "-s",
        "--states",
        type=click.Choice(
            ["singlets", "triplets", "50-50"], case_sensitive=False
        ),
        default="singlets",
        help="States for closed-shell singlet systems.\n"
        'Options choice =["Singlets", "Triplets", "50-50"]',
    )
    @click.option(
        "-r",
        "--root",
        type=int,
        default=1,
        help="Specifies the “state of interest”. The default is the first excited state (N=1).",
    )
    @click.option(
        "-n",
        "--nstates",
        type=int,
        default=3,
        help="Solve for M states (the default is 3). "
        "If 50-50, this gives the number of each type of state to solve "
        "(i.e., 3 singlets and 3 triplets).",
    )
    @click.option(
        "-e",
        "--eqsolv",
        type=str,
        default=None,
        help="Whether to perform equilibrium or non-equilibrium PCM solvation. "
        "NonEqSolv is the default except for excited state opt and when "
        "excited state density is requested (e.g., Density=Current or All).",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_gaussian_options
@click_gaussian_settings_options
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

    # obtain Gaussian Settings from filename, if supplied; otherwise return defaults

    if filename is None:
        # for cases where filename is not supplied, eg, get structure from pubchem
        job_settings = GaussianJobSettings.default()
        logger.info(
            f"No filename is supplied and Gaussian default settings are used:\n{job_settings.__dict__} "
        )
    elif filename.endswith((".com", ".inp", ".out", ".log")):
        # filename supplied - we would want to use the settings from here and do not use any defaults!
        job_settings = GaussianJobSettings.from_filepath(filename)
    elif filename.endswith(".xyz"):
        job_settings = GaussianJobSettings.default()
    else:
        raise ValueError(
            f"Unrecognised filetype {filename} to obtain GaussianJobSettings"
        )

    # Update keywords
    keywords = (
        "charge",
        "multiplicity",
    )  # default keywords to merge filename charge and multiplicity
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

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give Gaussian input filename or name to be be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"

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


@gaussian.result_callback()
@click.pass_context
def gaussian_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
