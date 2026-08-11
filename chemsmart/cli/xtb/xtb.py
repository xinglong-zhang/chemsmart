"""Click group for ChemSmart-managed xTB calculations."""

import functools
import logging
import os

import click

from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xtb import XTB_ALL_METHODS, XTB_ALL_SOLVENT_MODELS
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)


def _root_num_gpus_was_explicit(ctx):
    """Return whether a parent CLI explicitly supplied ``--num-gpus``."""

    current = ctx
    while current is not None:
        try:
            source = current.get_parameter_source("num_gpus")
        except (KeyError, TypeError):
            source = None
        if source is not None:
            return source == click.core.ParameterSource.COMMANDLINE
        current = current.parent
    return False


def _bind_cpu_only_runner(ctx):
    """Use CPU xTB on GPU-capable hosts; reject an explicit GPU request."""

    jobrunner = ctx.obj.get("jobrunner")
    if jobrunner is None:
        return
    explicit = _root_num_gpus_was_explicit(ctx)
    requested = getattr(jobrunner, "num_gpus", None)
    if explicit and (
        isinstance(requested, bool)
        or not isinstance(requested, int)
        or requested != 0
    ):
        raise click.UsageError(
            "ChemSmart's xTB execution plane is CPU-only; an explicit "
            "--num-gpus request must be 0."
        )
    setattr(jobrunner, "xtb_gpu_request_explicit", explicit)
    if not explicit and isinstance(requested, int) and not isinstance(
        requested, bool
    ):
        # Generic runners inherit the server inventory.  That inventory is
        # not an xTB resource request, so normalize this program node to CPU.
        jobrunner.num_gpus = 0


def require_xtb_filename(ctx):
    """Reject real commands without input while preserving leaf help."""
    if ctx.obj.get("xtb_missing_filename"):
        raise ValueError("xTB jobs require -f/--filename.")


def click_xtb_options(function):
    @click.option(
        "--project",
        "-p",
        type=str,
        default=None,
        help="Named xTB project or explicit project-YAML path.",
    )
    @functools.wraps(function)
    def wrapper_common_options(*args, **kwargs):
        return function(*args, **kwargs)

    return wrapper_common_options


def click_xtb_settings_options(function):
    @click.option("-c", "--charge", type=int, default=None)
    @click.option("-m", "--multiplicity", type=int, default=None)
    @click.option(
        "-g",
        "--gfn-version",
        type=click.Choice(XTB_ALL_METHODS, case_sensitive=False),
        default=None,
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=click.Choice(XTB_ALL_SOLVENT_MODELS, case_sensitive=False),
        default=None,
    )
    @click.option("-si", "--solvent-id", type=str, default=None)
    @click.option("--grad/--no-grad", default=None)
    @functools.wraps(function)
    def wrapper_common_options(*args, **kwargs):
        return function(*args, **kwargs)

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
    """Prepare an xTB sp, opt, or hess job through ChemSmart."""
    from chemsmart.jobs.xtb.settings import XTBJobSettings
    from chemsmart.settings.xtb import XTBProjectSettings

    ctx.ensure_object(dict)
    if filename is None:
        ctx.obj["xtb_missing_filename"] = True
        return

    _bind_cpu_only_runner(ctx)

    if (solvent_model is None) != (solvent_id is None):
        raise click.UsageError(
            "--solvent-model and --solvent-id must be supplied together."
        )

    if project is not None and os.path.isfile(project):
        project_reference = os.path.abspath(project)
    else:
        project_reference = project
    filename = os.path.abspath(filename)
    project_settings = XTBProjectSettings.from_project(project_reference)
    project_explicit_fields = project_settings.explicit_fields(
        ctx.invoked_subcommand
    )
    lower_filename = filename.lower()
    if lower_filename.endswith((".com", ".gjf", ".inp", ".out", ".log")):
        job_settings = XTBJobSettings.from_filepath(filename)
    else:
        job_settings = XTBJobSettings.default()

    # Only explicit CLI/project fields may replace source electronic state.
    # build_xtb_jobs resolves omitted fields independently for every selected
    # molecule and re-runs the electron-count parity check.
    keywords = []
    explicit_state_fields = set(
        project_explicit_fields & {"charge", "multiplicity"}
    )
    if charge is not None:
        job_settings.charge = charge
        keywords.append("charge")
        explicit_state_fields.add("charge")
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity
        keywords.append("multiplicity")
        explicit_state_fields.add("multiplicity")
    if gfn_version is not None:
        job_settings.gfn_version = gfn_version.lower()
        keywords.append("gfn_version")
    if solvent_model is not None:
        job_settings.solvent_model = solvent_model.lower()
        job_settings.solvent_id = solvent_id.strip().lower()
        keywords.extend(("solvent_model", "solvent_id"))
    if grad is not None:
        job_settings.grad = grad
        keywords.append("grad")
    job_settings.validate()

    molecules = Molecule.from_filepath(
        filepath=filename, index=":", return_list=True
    )
    if not molecules:
        raise ValueError(f"Could not obtain a molecule from {filename!r}.")

    if label is not None and append_label is not None:
        raise ValueError("Give --label or --append-label, not both.")
    if append_label is not None:
        label = os.path.splitext(os.path.basename(filename))[0]
        label = f"{label}_{append_label}"
    if label is None:
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

    ctx.obj.update(
        {
            "project_settings": project_settings,
            "job_settings": job_settings,
            "keywords": tuple(keywords),
            "molecules": molecules,
            "molecule_indices": molecule_indices,
            "label": label,
            "filename": filename,
            "project_reference": project_reference,
            "project_source_file": project_settings.source_file,
            "explicit_state_fields": frozenset(explicit_state_fields),
        }
    )


xtb.semantic_required_options = (("filename", "file"),)


@xtb.result_callback()
@click.pass_context
def xtb_process_pipeline(ctx, *args, **kwargs):
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    return args[0]
