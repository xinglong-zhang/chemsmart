"""
ORCA Command Line Interface

This module provides the main CLI interface for ORCA quantum chemistry
calculations. It defines common options, settings configurations, and the
main ORCA command group that serves as the entry point for all ORCA-related
operations.
"""

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


def click_orca_options(f):
    """
    Common click options decorator for ORCA jobs.

    This decorator adds common command-line options that are shared across
    different ORCA job types, specifically project settings.
    """

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_orca_settings_options(f):
    """
    Common click options decorator for ORCA computational settings.

    This decorator adds comprehensive command-line options for configuring
    ORCA calculations including file I/O, molecular properties, computational
    methods, basis sets, SCF settings, and various quantum chemistry
    parameters.
    """

    @click.option(
        "-t", "--title", type=str, default=None, help="ORCA job title."
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
        "-A",
        "--ab-initio",
        type=str,
        default=None,
        help="Ab initio method to be used.",
    )
    @click.option(
        "-x",
        "--functional",
        type=str,
        default=None,
        help="New functional to run.",
    )
    @click.option(
        "-D",
        "--dispersion",
        type=str,
        default=None,
        help="Dispersion for DFT functional.",
    )
    @click.option(
        "-b", "--basis", type=str, default=None, help="New basis set to run."
    )
    @click.option(
        "-a",
        "--aux-basis",
        type=str,
        default=None,
        help="Auxiliary basis set.",
    )
    @click.option(
        "-e",
        "--extrapolation-basis",
        type=str,
        default=None,
        help="Extrapolation basis set.",
    )
    @click.option(
        "-d",
        "--defgrid",
        type=click.Choice(
            ["defgrid1", "defgrid2", "defgrid3"], case_sensitive=False
        ),
        default="defgrid2",  # default used in ORCA is defgrid2
        help="Grid for numerical integration. Choices are "
        "['defgrid1', 'defgrid2', 'defgrid3'].",
    )
    @click.option(
        "--scf-tol",
        type=click.Choice(
            [
                "NormalSCF",
                "LooseSCF",
                "SloppySCF",
                "StrongSCF",
                "TightSCF",
                "VeryTightSCF",
                "ExtremeSCF",
            ]
        ),
        default=None,
        help="SCF convergence tolerance.",
    )
    @click.option(
        "--scf-algorithm",
        type=click.Choice(
            ["GDIIS", "DIIS", "SOSCF", "AutoTRAH"], case_sensitive=False
        ),  # SOSCF is an approximately quadratically convergent variant of
        # the SCF procedure.
        # In cases where conventional SCF procedures (DIIS/KDIIS/SOSCF) struggle,
        # we invoke TRAH-SCF automatically (AutoTRAH).
        default=None,
        help="SCF algorithm to use.",
    )
    @click.option(
        "--scf-maxiter",
        type=int,
        default=None,
        help="Maximum number of SCF iterations.",
    )
    @click.option(
        "--scf-convergence",
        type=float,
        default=None,
        help="SCF convergence criterion.",
    )
    @click.option(
        "--dipole/--no-dipole",
        default=None,
        type=bool,
        help="Enable dipole moment calculation.",
    )
    @click.option(
        "--quadrupole/--no-quadrupole",
        default=None,
        type=bool,
        help="Enable quadrupole moment calculation.",
    )
    @click.option(
        "--mdci-cutoff",
        type=click.Choice(["loose", "normal", "tight"], case_sensitive=False),
        default=None,
        help="MDCI cutoff. Choices are ['loose', 'normal', 'tight'].",
    )
    @click.option(
        "--mdci-density",
        type=click.Choice(
            ["none", "unrelaxed", "relaxed"], case_sensitive=False
        ),
        default=None,
        help="MDCI density. Choices are ['none', 'unrelaxed', 'relaxed'].",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="Additional route parameters.",
    )
    @click.option(
        "--forces/--no-forces",
        default=False,
        help="Enable forces calculation.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_orca_jobtype_options(f):
    """
    Common click options decorator for ORCA job type specifications.

    This decorator adds command-line options for specifying ORCA job types
    and related parameters for geometry optimizations, transition state
    searches, scans, and coordinate constraints.
    """

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='ORCA job type. Options: ["opt", "ts", "modred", "scan", "sp"].',
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="List of coordinates to be fixed for modred or scan jobs. "
        "1-indexed.",
    )
    @click.option(
        "-x",
        "--dist-start",
        default=None,
        help="Starting distance to scan, in Angstroms.",
    )
    @click.option(
        "-y",
        "--dist-end",
        default=None,
        help="Ending distance to scan, in Angstroms.",
    )
    @click.option(
        "-n",
        "--num-steps",
        default=None,
        help="Number of intermediate points for coordinate scans.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_orca_options
@click_filename_options
@click_file_label_and_index_options
@click_orca_settings_options
@click_pubchem_options
@click.pass_context
def orca(
    ctx,
    project,
    filename,
    label,
    append_label,
    title,
    charge,
    multiplicity,
    ab_initio,
    functional,
    dispersion,
    basis,
    aux_basis,
    extrapolation_basis,
    defgrid,
    scf_tol,
    scf_algorithm,
    scf_maxiter,
    scf_convergence,
    dipole,
    quadrupole,
    mdci_cutoff,
    mdci_density,
    index,
    additional_route_parameters,
    forces,
    pubchem,
):
    """
    Main CLI command group for running ORCA jobs using the chemsmart framework.

    This function serves as the primary entry point for all ORCA quantum
    chemistry calculations. It processes command-line arguments, configures
    job settings, loads molecular structures, and prepares the context for
    subcommands.
    """

    from chemsmart.jobs.orca.settings import ORCAJobSettings
    from chemsmart.settings.orca import ORCAProjectSettings

    # get project settings
    project_settings = ORCAProjectSettings.from_project(project)
    logger.debug(f"Loaded project settings: {project_settings}")

    # obtain ORCA Settings from filename, if supplied; otherwise return
    # defaults

    if filename is None:
        # for cases where filename is not supplied, eg, get structure from
        # pubchem
        job_settings = ORCAJobSettings.default()
        logger.info(
            f"No filename supplied, using ORCA default settings: "
            f"{job_settings.__dict__}"
        )
    elif filename.endswith((".com", ".inp", ".out", ".log")):
        # filename supplied - we would want to use the settings from here and
        # do not use any defaults!
        job_settings = ORCAJobSettings.from_filepath(filename)
        logger.info(f"Loaded ORCA settings from file: {filename}")
    elif filename.endswith(".xyz"):
        job_settings = ORCAJobSettings.default()
        logger.info(f"Using default ORCA settings for XYZ file: {filename}")
    else:
        raise ValueError(
            f"Unrecognised filetype {filename} to obtain ORCAJobSettings"
        )

    # Update keywords based on command-line arguments
    keywords = (
        "charge",
        "multiplicity",
    )  # default keywords to merge filename charge and multiplicity
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
    if title is not None:
        job_settings.title = title
        keywords += ("title",)
    if ab_initio is not None:
        job_settings.ab_initio = ab_initio
        keywords += ("ab_initio",)
    if functional is not None:
        job_settings.functional = functional
        keywords += ("functional",)  # update keywords
    if dispersion is not None:
        job_settings.dispersion = dispersion
        keywords += ("dispersion",)
    if basis is not None:
        job_settings.basis = basis
        keywords += ("basis",)
    if aux_basis is not None:
        job_settings.aux_basis = aux_basis
        keywords += ("aux_basis",)
    if extrapolation_basis is not None:
        job_settings.extrapolation_basis = extrapolation_basis
        keywords += ("extrapolation_basis",)
    if defgrid is not None:
        job_settings.defgrid = defgrid
        keywords += ("defgrid",)
    if scf_tol is not None:
        job_settings.scf_tol = scf_tol
        keywords += ("scf_tol",)
    if scf_algorithm is not None:
        job_settings.scf_algorithm = scf_algorithm
        keywords += ("scf_algorithm",)
    if scf_maxiter is not None:
        job_settings.scf_maxiter = scf_maxiter
        keywords += ("scf_maxiter",)
    if scf_convergence is not None:
        job_settings.scf_convergence = scf_convergence
        keywords += ("scf_convergence",)
    if dipole is not None:
        job_settings.dipole = dipole
        keywords += ("dipole",)
    if quadrupole is not None:
        job_settings.quadrupole = quadrupole
        keywords += ("quadrupole",)
    if mdci_cutoff is not None:
        job_settings.mdci_cutoff = mdci_cutoff
        keywords += ("mdci_cutoff",)
    if mdci_density is not None:
        job_settings.mdci_density = mdci_density
        keywords += ("mdci_density",)
    if additional_route_parameters is not None:
        job_settings.additional_route_parameters = additional_route_parameters
        keywords += ("additional_route_parameters",)
    if forces is not None:
        job_settings.forces = forces
        keywords += ("forces",)

    # obtain molecule structure from file or PubChem
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
        logger.debug(f"Obtained molecules {molecules} from {filename}")

    if pubchem:
        molecules = Molecule.from_pubchem(identifier=pubchem, return_list=True)
        assert (
            molecules is not None
        ), f"Could not obtain molecule from PubChem {pubchem}!"
        logger.debug(f"Obtained molecule {molecules} from PubChem {pubchem}")

    # update job labels for output file naming
    if label is not None and append_label is not None:
        raise ValueError(
            "Only give ORCA input filename or name to be appended, "
            "but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
        logger.debug(f"Created label with append: {label}")
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"
        logger.debug(f"Created default label: {label}")

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
    else:
        # If molecules is a list but molecule_indices is not set, create sequential indices
        if molecule_indices is None:
            molecule_indices = list(range(1, len(molecules) + 1))

    logger.debug(f"Final molecules list: {molecules}")
    logger.debug(f"Molecule indices: {molecule_indices}")
    logger.debug(f"Job settings keywords: {keywords}")

    # store objects in context for subcommands
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = molecules
    ctx.obj["molecule_indices"] = (
        molecule_indices  # Store original 1-based indices
    )
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@orca.result_callback()
@click.pass_context
def orca_process_pipeline(ctx, *args, **kwargs):
    """
    Result callback function for processing ORCA command pipeline.

    This function is executed after the ORCA subcommand completes and
    handles the final processing of results. It updates the context
    with subcommand information and returns the processed results.
    """
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    logger.debug(
        f"Pipeline completed for subcommand: {ctx.invoked_subcommand}"
    )
    return args[0]
