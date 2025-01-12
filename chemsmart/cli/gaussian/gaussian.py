import click
import functools
import logging

from chemsmart.io.molecules.structure import Molecule
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.cli import MyGroup

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
        "-c", "--charge", type=int, default=None, help="charge of the atoms"
    )
    @click.option(
        "-m",
        "--multiplicity",
        type=int,
        default=None,
        help="multiplicity of the atoms",
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
        help="index of atom of ase db to be used",
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
        help="additional information to be appended at the end of the input file. E.g, scrf=read",
    )
    @click.option(
        "-C",
        "--custom-solvent",
        type=str,
        default=None,
        help="additional information to be appended at the end of the input file. E.g, scrf=read",
    )
    @click.option(
        "-d",
        "--dieze-tag",
        type=str,
        default=None,
        help="dieze tag for gaussian job; possible options include "
        '"n", "p", "t" to get "#n", "#p", "#t", respectively',
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
    help="Queries structure from PubChem using name, smiles, cid and conformer informaotion.",
)
@click.pass_context
def gaussian(  # noqa: PLR0912, PLR0915
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
    pubchem,
):
    import os
    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.settings.gaussian import GaussianProjectSettings

    from chemsmart.utils.utils import string2index

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

    # obtain atoms structure
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

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)

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

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        molecules = molecules[string2index(index)]
        molecules = [molecules]

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = (
        molecules  # molecules as a list, as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename
    jobrunner = ctx.obj["jobrunner"]
    logger.info(f"Scratch: {jobrunner.scratch}")


@gaussian.result_callback()
@click.pass_context
def gaussian_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]


@gaussian.command("opt", cls=MyCommand)
@click_job_options
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    help="Indices of atoms to freeze for constrained optimization.",
)
@click.pass_context
def opt(ctx, freeze_atoms, **kwargs):
    # folder = ctx.obj["folder"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `sub.py gaussian -c <user_charge> -m <user_multiplicity>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project opt settings with job settings from cli keywords from cli.gaussian.py subcommands
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    # get atoms
    molecules = ctx.obj["molecules"]
    molecule = molecules[
        -1
    ]  # get last atom from list of atoms from cli.gaussian.py subcommands
    # index = '-1' would access the right structure from the list of atoms returned from cli.gaussian.py subcommands
    # user specified index was used there to return the right atoms and store it as a list of single element/itself

    # get label for the job
    label = ctx.obj["label"]

    # # Set atoms to freeze
    # from ase.constraints import FixAtoms
    #
    # from chemsmart.utils.utils import get_list_from_string_range
    #
    # if freeze_atoms is not None:
    #     frozen_atoms_list = get_list_from_string_range(freeze_atoms)
    #     atoms.set_constraint(FixAtoms(frozen_atoms_list))
    #
    # logger.info(f"Opt settings from project: {opt_settings.__dict__}")
    #
    from chemsmart.jobs.gaussian import GaussianGeomOptJob

    return GaussianGeomOptJob(
        molecule=molecule,
        settings=opt_settings,
        label=label,
        **kwargs,
    )
