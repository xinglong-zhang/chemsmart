import functools
import logging
import os

import click

from pyatoms.utils.cli import MyGroup

logger = logging.getLogger(__name__)


def click_orca_options(f):
    """Common click options for ORCA jobs."""

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @click.option(
        "--scratch/--no-scratch",
        type=bool,
        default=True,
        help="To run in scratch or without scratch folder.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_orca_settings_options(f):
    """Common click options for ORCA Settings."""

    @click.option(
        "-f",
        "--filename",
        type=str,
        default=None,
        help="filename from which new ORCA input is prepared.",
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
        "-t", "--title", type=str, default=None, help="ORCA job title."
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
        default="-1",
        help="index of atom to use; default is last atom",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="additional route parameters",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_orca_options
@click_orca_settings_options
@click.option(
    "-P",
    "--pubchem",
    type=str,
    default=None,
    help="Queries structure from PubChem using name, smiles, cid and conformer informaotion.",
)
@click.pass_context
def orca(  # noqa: PLR0912, PLR0915
    ctx,
    project,
    scratch,
    filename,
    label,
    append_label,
    title,
    charge,
    multiplicity,
    functional,
    basis,
    index,
    additional_route_parameters,
    pubchem,
):
    from ase.io.formats import string2index

    from pyatoms.io.ase.atoms import AtomsWrapper
    from pyatoms.settings.projects.orca import ORCAProjectSettings

    # get project settings
    project_settings = ORCAProjectSettings.from_project_name(project)

    # obtain ORCA Settings from filename, if supplied; otherwise return defaults
    from pyatoms.jobs.orca.settings import ORCAJobSettings

    if filename is None:
        job_settings = ORCAJobSettings.default()
        logger.info(
            f"No filename is supplied and Gaussian default settings are used:\n{job_settings.__dict__} "
        )
    else:
        # filename supplied - we would want to use the settings from here and do not use any defaults!
        job_settings = ORCAJobSettings.from_filepath(filename)
        logger.info(f"Job settings from {filename}:\n{job_settings.__dict__}")

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
    if additional_route_parameters is not None:
        job_settings.additional_route_parameters = additional_route_parameters
        keywords += ("additional_route_parameters",)

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
        atoms = AtomsWrapper.from_filepath(
            filepath=filename, index=":", return_list=True
        )

    if pubchem:
        atoms = AtomsWrapper.from_pubchem(identifier=pubchem, return_list=True)

    # update labels
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give ORCA input filename or name to be be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"

    # return list of atoms
    atoms = atoms[string2index(index)]

    if not isinstance(atoms, list):
        # if somehow atoms is not a list, make it a list
        atoms = [atoms]

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["atoms"] = (
        atoms  # atoms as a list as some jobs requires all structures to be used
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename
    jobrunner = ctx.obj["jobrunner"]
    jobrunner.scratch = scratch


@orca.result_callback()
@click.pass_context
def orca_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
