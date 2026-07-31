import functools
import logging
import os

import click

from chemsmart.cli.database.database import click_database_id_options
from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_pubchem_options,
)
from chemsmart.database.utils import is_chemsmart_database
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)


def click_crest_options(f):
    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_crest_settings_options(f):
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
        "-g",
        "--gfn-version",
        type=click.Choice(
            ["gfn1", "gfn2", "gfnff", "gfn2//gfnff"], case_sensitive=False
        ),
        default=None,
        help="GFN-xTB method version.",
    )
    @click.option(
        "-O",
        "--optimization-level",
        type=click.Choice(
            [
                "crude",
                "sloppy",
                "loose",
                "lax",
                "normal",
                "tight",
                "vtight",
                "extreme",
            ],
            case_sensitive=False,
        ),
        default=None,
        help="Level for GFN-xTB optimizations.",
    )
    @click.option(
        "--ewin",
        "energy_window",
        type=float,
        default=None,
        help="Energy window in kcal/mol for conformer selection.",
    )
    @click.option(
        "--nci/--no-nci",
        default=None,
        help="Enable non-covalent interaction mode.",
    )
    @click.option(
        "-r",
        "--additional-flags",
        type=str,
        default=None,
        help="Additional CREST CLI flags.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_crest_solvent_options(f):
    """Common click options for CREST solvent settings."""

    @click.option(
        "--remove-solvent/--no-remove-solvent",
        default=False,
        help="Remove the solvent model from the job.",
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="Implicit solvent model.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="Implicit solvent identifier.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_crest_constrain_options(f):
    """Reusable CLI options for CREST geometry constraints."""

    @click.option(
        "-c",
        "--constraints",
        type=str,
        default=None,
        help="List of coordinates to constrain during conformational search. "
        "1-indexed. Distance: [i,j]; angle: [i,j,k]; dihedral: [i,j,k,l]. "
        "Example: [[1,2],[5,7,8],[3,4,1,7]].",
    )
    @click.option(
        "-f",
        "--force-constant",
        type=float,
        default=None,
        help="Force constant for constraints.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_crest_options
@click_filename_options
@click_file_label_and_index_options
@click_database_id_options
@click_crest_settings_options
@click_crest_solvent_options
@click_pubchem_options
@click.pass_context
def crest(
    ctx,
    project,
    filename,
    label,
    append_label,
    index,
    record_index,
    record_id,
    structure_id,
    structure_index,
    molecule_id,
    charge,
    multiplicity,
    gfn_version,
    optimization_level,
    energy_window,
    nci,
    additional_flags,
    remove_solvent,
    solvent_model,
    solvent_id,
    pubchem,
):
    """CLI subcommand for CREST conformational search jobs."""
    # --mid is not supported for job submission
    if molecule_id is not None:
        raise click.UsageError(
            "--mid/--molecule-id is not supported for CREST job submission. "
            "Use --sid/--structure-id or --ri/--rid with -i/--si instead."
        )
    # -i/--index and --si/--structure-index are equivalent aliases
    if index is not None and structure_index is not None:
        raise click.UsageError(
            "-i/--index and --si/--structure-index are mutually exclusive. "
            "Use only one to specify the structure index."
        )
    # If --si is given, treat it as -i so all downstream code uses index
    if structure_index is not None:
        index = structure_index

    is_chemsmart_db = is_chemsmart_database(filename)
    if is_chemsmart_db:
        record_selectors = [record_index is not None, record_id is not None]
        if sum(record_selectors) + (structure_id is not None) != 1:
            raise click.UsageError(
                "For chemsmart database input, select exactly one of "
                "--ri/--record-index, --rid/--record-id, or "
                "--sid/--structure-id."
            )
        if index is not None and not any(record_selectors):
            raise click.UsageError(
                "For chemsmart database input, -i/--index (or --si/--structure-index) "
                "can only be used together with --ri/--record-index or --rid/--record-id."
            )

    from chemsmart.jobs.crest.settings import CRESTJobSettings
    from chemsmart.settings.crest import CRESTProjectSettings

    project_settings = CRESTProjectSettings.from_project(project)

    # Obtain CREST settings from filename when possible.
    if filename is None:
        job_settings = CRESTJobSettings.default()
        logger.info(
            f"No filename is supplied and CREST default settings are used:\n"
            f"{job_settings.__dict__} "
        )
    elif filename.endswith((".com", ".gjf", ".inp", ".out", ".log")):
        job_settings = CRESTJobSettings.from_filepath(filename)
    elif filename.endswith(".db"):
        if is_chemsmart_db:
            job_settings = CRESTJobSettings.from_database(
                filepath=filename,
                record_index=record_index,
                record_id=record_id,
                structure_index=index or "-1",
                structure_id=structure_id,
            )
        else:
            logger.debug(
                f"File {filename} is not a valid chemsmart database file."
            )
            job_settings = CRESTJobSettings.default()
    else:
        logger.debug(
            f"Falling back to default CREST job settings for file {filename}."
        )
        job_settings = CRESTJobSettings.default()

    keywords = (
        "charge",
        "multiplicity",
    )
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
    if gfn_version is not None:
        job_settings.gfn_version = gfn_version.lower()
        keywords += ("gfn_version",)
    if optimization_level is not None:
        job_settings.optimization_level = optimization_level.lower()
        keywords += ("optimization_level",)
    if energy_window is not None:
        job_settings.energy_window = energy_window
        keywords += ("energy_window",)
    if nci is not None:
        job_settings.nci = nci
        keywords += ("nci",)
    if additional_flags is not None:
        job_settings.additional_flags = additional_flags
        keywords += ("additional_flags",)

    if remove_solvent:
        job_settings.solvent_model = None
        job_settings.solvent_id = None
        keywords += ("solvent_model", "solvent_id")
    else:
        if solvent_model is not None:
            job_settings.solvent_model = solvent_model.lower()
            keywords += ("solvent_model",)
        if solvent_id is not None:
            job_settings.solvent_id = solvent_id.lower()
            keywords += ("solvent_id",)

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
        if is_chemsmart_db:
            if structure_id is not None:
                molecules = Molecule.from_filepath(
                    filepath=filename,
                    return_list=True,
                    structure_id=structure_id,
                )
            else:
                molecules = Molecule.from_filepath(
                    filepath=filename,
                    index=index or "-1",
                    return_list=True,
                    record_index=record_index,
                    record_id=record_id,
                )

            assert (
                molecules is not None
            ), f"Could not obtain molecule from database {filename}!"
            logger.debug(
                f"Obtained database molecule {molecules} from {filename}"
            )
        else:
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

    if label is not None and append_label is not None:
        raise ValueError(
            "Only give CREST input filename or name to be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        if is_chemsmart_db:
            if structure_id is not None:
                label = f"{label}_SID-{structure_id}"
            elif record_id is not None:
                label = f"{label}_RID-{record_id}"
            elif record_index is not None:
                label = f"{label}_RI-{record_index}"
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        if filename:
            label = os.path.splitext(os.path.basename(filename))[0]
        else:
            label = "output"
        if is_chemsmart_db:
            if structure_id is not None:
                label = f"{label}_SID-{structure_id}"
            elif record_id is not None:
                label = f"{label}_RID-{record_id}"
            elif record_index is not None:
                label = f"{label}_RI-{record_index}"
        if ctx.invoked_subcommand:
            label = f"{label}_{ctx.invoked_subcommand}"
    label = clean_label(label)

    # if user has specified an index to use to access particular structure
    # then return that structure as a list and track the original indices
    molecule_indices = None
    if index is not None and not is_chemsmart_db:
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

    logger.debug(f"CREST project settings: {project_settings}")
    logger.debug(f"CREST job settings before merge: {job_settings.__dict__}")
    logger.debug(f"CREST merge keywords: {keywords}")
    logger.debug(f"CREST selected molecule count: {len(molecules)}")
    logger.debug(f"CREST selected molecule indices: {molecule_indices}")
    logger.debug(f"CREST job label: {label}")

    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = tuple(keywords)
    ctx.obj["molecules"] = molecules
    ctx.obj["molecule_indices"] = molecule_indices
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@crest.result_callback()
@click.pass_context
def crest_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
