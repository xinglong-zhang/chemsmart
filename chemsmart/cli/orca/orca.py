import functools
import logging
import os

import click

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


def click_orca_options(f):
    """Common click options for ORCA jobs."""

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
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
        help="Grid for numerical integration. Choices are ['defgrid1', 'defgrid2', 'defgrid3']",
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
        ),  # SOSCF is an approximately quadratically convergent variant of the SCF procedure
        # In cases conventional SCF procedures (DIIS/KDIIS/SOSCF) struggle, we invoke TRAH-SCF
        # automatically (AutoTRAH).
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
        help="Dipole moment calculation.",
    )
    @click.option(
        "--quadrupole/--no-quadrupole",
        default=None,
        type=bool,
        help="Quadrupole moment calculation.",
    )
    @click.option(
        "--mdci-cutoff",
        type=click.Choice(["loose", "normal", "tight"], case_sensitive=False),
        default=None,
        help="MDCI cutoff. Choices are ['loose', 'normal', 'tight']",
    )
    @click.option(
        "--mdci-density",
        type=click.Choice(
            ["none", "unrelaxed", "relaxed"], case_sensitive=False
        ),
        default=None,
        help="MDCI density. Choices are ['none', 'unrelaxed', 'relaxed']",
    )
    @click.option(
        "-i",
        "--index",
        type=str,
        default=None,
        help="index of molecule to use; default to the last molecule structure.",
    )
    @click.option(
        "-r",
        "--additional-route-parameters",
        type=str,
        default=None,
        help="additional route parameters",
    )
    @click.option(
        "--forces/--no-forces", default=False, help="Forces calculation."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_orca_jobtype_options(f):
    """Common click options for ORCA link/crest jobs."""

    @click.option(
        "-j",
        "--jobtype",
        type=str,
        default=None,
        help='ORCA job type. Options: ["opt", "ts", "modred", "scan", "sp"]',
    )
    @click.option(
        "-c",
        "--coordinates",
        default=None,
        help="List of coordinates to be fixed for modred or scan job. 1-indexed.",
    )
    @click.option('-x', '--dist-start', default=None, help='starting distance to scan, in Angstroms.')
    @click.option('-y', '--dist-end', default=None, help='ending distance to scan, in Angstroms.')
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

    from chemsmart.jobs.orca.settings import ORCAJobSettings
    from chemsmart.settings.orca import ORCAProjectSettings

    # get project settings
    project_settings = ORCAProjectSettings.from_project(project)

    # obtain ORCA Settings from filename, if supplied; otherwise return defaults

    if filename is None:
        # for cases where filename is not supplied, eg, get structure from pubchem
        job_settings = ORCAJobSettings.default()
        logger.info(
            f"No filename is supplied and Gaussian default settings are used:\n{job_settings.__dict__} "
        )
    elif filename.endswith((".com", ".inp", ".out", ".log")):
        # filename supplied - we would want to use the settings from here and do not use any defaults!
        job_settings = ORCAJobSettings.from_filepath(filename)
    elif filename.endswith(".xyz"):
        job_settings = ORCAJobSettings.default()
    else:
        raise ValueError(
            f"Unrecognised filetype {filename} to obtain ORCAJobSettings"
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
            filepath=filename, index=";", return_list=True
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

    # if user has specified an index to use to access particular structure
    # then return that structure as a list
    if index is not None:
        # return list of molecules
        molecules = molecules[string2index_1based(index)]

    logger.debug(f"Obtained molecules: {molecules}")

    # store objects
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = molecules
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@orca.result_callback()
@click.pass_context
def orca_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
