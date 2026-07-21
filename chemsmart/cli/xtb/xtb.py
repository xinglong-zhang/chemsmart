import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)


def require_xtb_filename(ctx):
    """Reject real xTB executions started without -f/--filename.

    The group callback defers this check so that `chemsmart run xtb <job>
    --help` still renders; the leaf commands call this once to enforce it.
    """
    if ctx.obj.get("xtb_missing_filename"):
        raise ValueError("xTB jobs require -f/--filename.")


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
    @click.option(
        "--grad/--no-grad",
        default=None,
        help="Enable gradient output.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_xtb_options
@click_filename_options
@click_file_label_and_index_options
@click_xtb_settings_options
@click.pass_context
def xtb(
    ctx,
    project,
    filename,
    label,
    append_label,
    index,
    charge,
    multiplicity,
    gfn_version,
    solvent_model,
    solvent_id,
    grad,
):
    """CLI subcommand for xTB semiempirical quantum chemistry jobs."""
    from chemsmart.jobs.xtb.settings import XTBJobSettings
    from chemsmart.settings.xtb import XTBProjectSettings

    ctx.ensure_object(dict)
    if filename is None:
        # Let `chemsmart run xtb <job> --help` render without requiring
        # molecular input. Real executions are rejected once in the leaf
        # command via require_xtb_filename() so the check lives in one place.
        ctx.obj["xtb_missing_filename"] = True
        return

    project_settings = XTBProjectSettings.from_project(project)
    if filename.endswith((".com", ".gjf", ".inp", ".out", ".log")):
        job_settings = XTBJobSettings.from_filepath(filename)
    else:
        job_settings = XTBJobSettings.default()

    keywords = ["charge", "multiplicity"]
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
    if gfn_version is not None:
        job_settings.gfn_version = gfn_version.lower()
        keywords.append("gfn_version")
    if solvent_model is not None:
        job_settings.solvent_model = solvent_model.lower()
        keywords.append("solvent_model")
    if solvent_id is not None:
        job_settings.solvent_id = solvent_id.lower()
        keywords.append("solvent_id")
    if grad is not None:
        job_settings.grad = grad
        keywords.append("grad")

    # Project YAML provides the xTB defaults; CLI/file-derived settings are
    # merged in the leaf command so opt/sp/hess each keep their jobtype.
    molecules = Molecule.from_filepath(
        filepath=filename, index=":", return_list=True
    )
    assert molecules is not None, f"Could not obtain molecule from {filename}!"

    if label is not None and append_label is not None:
        raise ValueError(
            "Only give xTB input filename or name to be appended, but not both!"
        )
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None and append_label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{ctx.invoked_subcommand}"
    label = clean_label(label)

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
