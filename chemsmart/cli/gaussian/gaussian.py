import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_pubchem_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import (
    return_objects_and_indices_from_string_index,
)

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
        "-t", "--title", type=str, default=None, help="Gaussian job title."
    )
    @click.option(
        "-c",
        "--charge",
        type=int,
        default=None,
        help="Charge of the molecule.",
    )
    @click.option(
        "-m",
        "--multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the molecule.",
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
        "-s",
        "--semiempirical",
        type=str,
        default=None,
        help="Semiempirical method to run.",
    )
    @click.option(
        "-o",
        "--additional-opt-options",
        type=str,
        default=None,
        help="Additional optimization options.",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="Additional route parameters.",
    )
    @click.option(
        "-A",
        "--append-additional-info",
        type=str,
        default=None,
        help="Additional information to be appended at the end of the "
        "input file (e.g., scrf=read).",
    )
    @click.option(
        "-C",
        "--custom-solvent",
        type=str,
        default=None,
        help="Custom solvent information to be appended at the end of the "
        "input file (e.g., scrf=read).",
    )
    @click.option(
        "-d",
        "--dieze-tag",
        type=str,
        default=None,
        help="Dieze tag for Gaussian job. Possible options include "
        '"n", "p", "t" to get "#n", "#p", "#t", respectively.',
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


def click_gaussian_irc_options(f):
    """Common click options for IRC-related jobs."""

    @click.option(
        "-fl/",
        "--flat-irc/--no-flat-irc",
        type=bool,
        default=False,
        help="Whether to run flat IRC or not.",
    )
    @click.option(
        "-pt",
        "--predictor",
        type=click.Choice(
            ["LQA", "HPC", "EulerPC", "DVV", "Euler"], case_sensitive=False
        ),
        default=None,
        help="Type of predictor used for IRC. Examples include [HPC, EulerPC, "
        "LQA, DVV, Euler].",
    )
    @click.option(
        "-rc",
        "--recorrect",
        type=click.Choice(["Never", "Always", "Test"], case_sensitive=False),
        default=None,
        help="Recorrection step of HPC and EulerPC IRCs. Options are: "
        '["Never", "Always", "Test"].',
    )
    @click.option(
        "-rs",
        "--recalc-step",
        type=int,
        default=6,
        help="Compute the Hessian analytically every N predictor steps or every "
        "|N| corrector steps if N < 0.",
    )
    @click.option(
        "-mp",
        "--maxpoints",
        type=int,
        default=512,
        help="Number of points along the reaction path to examine.",
    )
    @click.option(
        "-mc",
        "--maxcycles",
        type=int,
        default=128,
        help="Maximum number of steps along the IRC to run.",
    )
    @click.option(
        "-ss",
        "--stepsize",
        type=int,
        default=20,
        help="Step size along the reaction path, in units of 0.01 Bohr.",
    )
    @click.option(
        "-d",
        "--direction",
        type=click.Choice(["forward", "reverse"], case_sensitive=False),
        default=None,
        help="Only run the forward or reverse IRC. Defaults to run both directions.",
    )
    @functools.wraps(f)
    def wrapper_irc_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_irc_options


def click_gaussian_jobtype_options(f):
    """Common click options for Gaussian link/crest jobs."""

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='Gaussian job type. Options: ["opt", "ts", "modred", "scan", '
        '"sp", "irc"].',
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="List of coordinates to be fixed for modred or scan jobs. "
        "1-indexed.",
    )
    @click.option(
        "-s",
        "--step-size",
        default=None,
        help="Step size for coordinate scans.",
    )
    @click.option(
        "-n",
        "--num-steps",
        default=None,
        help="Number of steps for coordinate scans.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_grouper_options(f):
    """Common click options for Gaussian grouper jobs."""

    @click.option(
        "-g",
        "--grouping-strategy",
        type=click.Choice(
            ["rmsd", "tanimoto", "isomorphism", "formula", "connectivity"],
            case_sensitive=False,
        ),
        default=None,
        help="Grouping strategy to use for Gaussian jobs.",
    )
    @click.option(
        "-t",
        "--threshold",
        type=float,
        default=None,
        help="Threshold for grouping. If not specified, uses strategy-specific "
        "defaults: RMSD=0.5, Tanimoto=0.9, Connectivity=0.0.",
    )
    @click.option(
        "-i",
        "--ignore-hydrogens",
        is_flag=True,
        default=False,
        help="Whether to ignore hydrogens in grouping.",
    )
    @click.option(
        "-p",
        "--num-procs",
        type=int,
        default=1,
        help="Number of processors to use for the grouper.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_solvent_options(f):
    """Common click options for Gaussian solvent settings."""

    @click.option(
        "--remove-solvent/--no-remove-solvent",
        "-r/ ",
        type=bool,
        default=False,
        help="Whether to remove the solvent model from the job. Defaults to project "
        "settings.",
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="Solvent model to be used for single point calculations.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="Solvent ID to be used for single point calculations.",
    )
    @click.option(
        "-so",
        "--solvent-options",
        type=str,
        default=None,
        help="Additional solvent options in scrf=() route. "
        "E.g., `iterative` in scrf=(smd,water,iterative) via "
        "chemsmart sub -s xz gaussian -p dnam -f output.log -a scrf_iter sp "
        "-so iterative.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_gaussian_td_options(f):
    """Common click options for Gaussian TDDFT calculations."""

    @click.option(
        "-s",
        "--states",
        type=click.Choice(
            ["singlets", "triplets", "50-50"], case_sensitive=False
        ),
        default="singlets",
        help="States for closed-shell singlet systems. "
        'Options: ["Singlets", "Triplets", "50-50"].',
    )
    @click.option(
        "-r",
        "--root",
        type=int,
        default=1,
        help='Specifies the "state of interest". '
        "The default is the first excited state (N=1).",
    )
    @click.option(
        "-n",
        "--nstates",
        type=int,
        default=3,
        help="Number of states to solve for (default is 3). "
        "If 50-50, this gives the number of each type of state to solve "
        "(i.e., 3 singlets and 3 triplets).",
    )
    @click.option(
        "-e",
        "--eqsolv",
        type=str,
        default=None,
        help="Whether to perform equilibrium or non-equilibrium PCM solvation. "
        "NonEqSolv is the default except for excited state optimization and when "
        "excited state density is requested (e.g., Density=Current or All).",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_gaussian_options
@click_filename_options
@click_file_label_and_index_options
@click_gaussian_settings_options
@click_pubchem_options
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
    """CLI subcommand for running Gaussian jobs using the chemsmart framework."""

    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.settings.gaussian import GaussianProjectSettings

    # get project settings
    project_settings = GaussianProjectSettings.from_project(project)

    # obtain Gaussian Settings from filename, if supplied;
    #  otherwise return defaults

    if filename is None:
        # for cases where filename is not supplied, eg,
        #  get structure from pubchem
        job_settings = GaussianJobSettings.default()
        logger.info(
            f"No filename is supplied and Gaussian default settings are used:\n"
            f"{job_settings.__dict__} "
        )
    elif filename.endswith((".com", "gjf", ".inp", ".out", ".log")):
        # filename supplied - we would want to use the settings from here
        #  and do not use any defaults!
        job_settings = GaussianJobSettings.from_filepath(filename)
    # elif filename.endswith((".xyz", ".pdb", ".mol", ".mol2", ".sdf", ".smi",
    #  ".cif", ".traj", ".gro", ".db")):
    else:
        job_settings = GaussianJobSettings.default()
    # else:
    #     raise ValueError(
    #         f"Unrecognised filetype {filename} to obtain GaussianJobSettings"
    #     )

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
            "[filename] or [pubchem] has not been specified!\n"
            "Please specify one of them!"
        )
    if filename and pubchem:
        raise ValueError(
            "Both [filename] and [pubchem] have been specified!\n"
            "Please specify only one of them."
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
            "Only give Gaussian input filename or name to be be appended,"
            "but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"

    label = clean_label(label)

    # if user has specified an index to use to access particular structure
    # then return that structure as a list and track the original indices
    molecule_indices = None
    if index is not None:
        molecules, molecule_indices = (
            return_objects_and_indices_from_string_index(
                list_of_objects=molecules, index=index
            )
        )

    if not isinstance(molecules, list):
        molecules = [molecules]
        if molecule_indices is not None and not isinstance(
            molecule_indices, list
        ):
            molecule_indices = [molecule_indices]

    logger.debug(f"Obtained molecules: {molecules}")
    logger.debug(f"Molecule indices: {molecule_indices}")

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["molecule_indices"] = (
        molecule_indices  # Store original 1-based indices
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@gaussian.result_callback()
@click.pass_context
def gaussian_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
