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


def click_xtb_options(f):
    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_xtb_settings_options(f):
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
            ["gfn0", "gfn1", "gfn2", "gfnff"], case_sensitive=False
        ),
        default=None,
        help="GFN-xTB method version.",
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
        default=None,
        help="Enable gradient output.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_xtb_solvent_options(f):
    """Common click options for xTB solvent settings."""

    @click.option(
        "--remove-solvent/--no-remove-solvent",
        default=False,
        help="Remove the solvent model from the job. Defaults to project "
        "settings.",
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="xTB implicit solvent model.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="xTB implicit solvent identifier.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_xtb_options
@click_filename_options
@click_file_label_and_index_options
@click_database_id_options
@click_xtb_settings_options
@click_xtb_solvent_options
@click_pubchem_options
@click.pass_context
def xtb(
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
    additional_route_parameters,
    charge,
    multiplicity,
    gfn_version,
    remove_solvent,
    solvent_model,
    solvent_id,
    grad,
    pubchem,
):
    """CLI subcommand for xTB semiempirical quantum chemistry jobs."""
    # --mid is not supported for job submission
    if molecule_id is not None:
        raise click.UsageError(
            "--mid/--molecule-id is not supported for xTB job submission. "
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

    from chemsmart.jobs.xtb.settings import XTBJobSettings
    from chemsmart.settings.xtb import XTBProjectSettings

    # get project settings
    project_settings = XTBProjectSettings.from_project(project)

    # obtain xTB settings from filename, if supplied; otherwise return defaults
    if filename is None:
        # for cases where filename is not supplied, eg,
        #  get structure from pubchem
        job_settings = XTBJobSettings.default()
        logger.info(
            f"No filename is supplied and xTB default settings are used:\n"
            f"{job_settings.__dict__} "
        )
    elif filename.endswith((".com", ".gjf", ".inp", ".out", ".log")):
        job_settings = XTBJobSettings.from_filepath(filename)
    elif filename.endswith(".db"):
        if is_chemsmart_db:
            job_settings = XTBJobSettings.from_database(
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
            job_settings = XTBJobSettings.default()
    else:
        logger.debug(
            f"Falling back to default xTB job settings for file {filename}."
        )
        job_settings = XTBJobSettings.default()

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
    if additional_route_parameters is not None:
        job_settings.additional_route_parameters = additional_route_parameters
        keywords += ("additional_route_parameters",)
    if grad is not None:
        job_settings.grad = grad
        keywords += ("grad",)

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
            "Only give xTB input filename or name to be appended, but not both!"
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
        label = os.path.splitext(os.path.basename(filename))[0]
        if is_chemsmart_db:
            if structure_id is not None:
                label = f"{label}_SID-{structure_id}"
            elif record_id is not None:
                label = f"{label}_RID-{record_id}"
            elif record_index is not None:
                label = f"{label}_RI-{record_index}"
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

    logger.debug(f"xTB project settings: {project_settings}")
    logger.debug(f"xTB job settings before merge: {job_settings.__dict__}")
    logger.debug(f"xTB merge keywords: {keywords}")
    logger.debug(f"xTB selected molecule count: {len(molecules)}")
    logger.debug(f"xTB selected molecule indices: {molecule_indices}")
    logger.debug(f"xTB job label: {label}")

    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = tuple(keywords)
    ctx.obj["molecules"] = molecules
    ctx.obj["molecule_indices"] = molecule_indices
    ctx.obj["label"] = label
    ctx.obj["filename"] = filename


@xtb.result_callback()
@click.pass_context
def xtb_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
