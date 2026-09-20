"""CLI group for YAML chain pipelines."""

import functools
import logging
import os

import click

from chemsmart.cli.chain.step_options import apply_step_option_tokens
from chemsmart.cli.chain.workflows import (
    CHAIN_WORKFLOW_COMMANDS,
    register_chain_workflows,
)
from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_job_options,
)
from chemsmart.cli.utils import (
    CHAIN_CLI_DEFAULTS_KEY,
    CHAIN_PROJECT_SETTINGS_KEY,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain_steps import build_chain_job, parse_chain_step
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)

_PIPELINE_STEPS_KEY = "chain_pipeline_steps"
_SKIP_COMPLETED_KEY = "chain_skip_completed"


def click_chain_options(f):
    """Project, charge, and multiplicity options for chain jobs."""

    @click.option(
        "--project",
        "-p",
        type=str,
        required=False,
        help=(
            "Chain project settings. Required for the -s/--steps pipeline "
            "and for pka/fukui/redox/reaction submit."
        ),
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
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _chain_label(filename, label, append_label):
    if label is not None and append_label is not None:
        raise click.UsageError(
            "Give only one of -l/--label and -a/--append-label."
        )
    if append_label is not None:
        stem = os.path.splitext(os.path.basename(filename))[0]
        label = f"{stem}_{append_label}"
    elif label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
    return clean_label(label)


def _molecule_from_filename(filename, index):
    molecules = Molecule.from_filepath(
        filepath=filename, index=":", return_list=True
    )
    if molecules is None:
        raise click.UsageError(f"Could not obtain molecule from {filename}.")
    if index is not None:
        molecules, _ = return_objects_and_indices_from_string_index(
            list_of_objects=molecules, index=index
        )
    if not isinstance(molecules, list):
        molecules = [molecules]
    if not molecules:
        raise click.UsageError(f"No molecule selected from {filename}.")
    return molecules[-1]


def _parse_cli_steps(steps, chain_settings):
    parsed = []
    for token in steps:
        try:
            step = parse_chain_step(token)
            step = apply_step_option_tokens(step)
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        if step.program not in ChainProjectSettings.PROGRAMS:
            allowed = ", ".join(ChainProjectSettings.PROGRAMS)
            raise click.UsageError(
                f"Unknown chain program {step.program!r}. "
                f"Allowed programs: {allowed}."
            )
        try:
            chain_settings.project_for(step.program)
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        parsed.append(step)
    return parsed


def _store_chain_invocation(
    ctx,
    chain_settings,
    filename,
    label,
    append_label,
    index,
    charge,
    multiplicity,
    skip_completed,
    steps,
):
    ctx.obj[CHAIN_PROJECT_SETTINGS_KEY] = chain_settings
    ctx.obj[CHAIN_CLI_DEFAULTS_KEY] = {
        "filename": filename,
        "label": label,
        "append_label": append_label,
        "index": index,
        "charge": charge,
        "multiplicity": multiplicity,
    }
    ctx.obj[_SKIP_COMPLETED_KEY] = skip_completed
    ctx.obj[_PIPELINE_STEPS_KEY] = steps


def build_chain_pipeline_from_ctx(ctx):
    """Build a ``ChainJob`` from chain-group ``-s/--steps`` stored on ``ctx``."""
    chain_settings = ctx.obj.get(CHAIN_PROJECT_SETTINGS_KEY)
    if chain_settings is None:
        raise click.UsageError(
            "-p/--project is required for the chain pipeline."
        )
    defaults = ctx.obj[CHAIN_CLI_DEFAULTS_KEY]
    steps = ctx.obj.get(_PIPELINE_STEPS_KEY) or ()
    if not steps:
        raise click.UsageError(
            "Custom pipeline requires at least one -s/--steps."
        )
    filename = defaults["filename"]
    if filename is None:
        raise click.UsageError(
            "Chain pipeline requires a structure file (-f/--filename)."
        )
    parsed_steps = _parse_cli_steps(steps, chain_settings)
    molecule = _molecule_from_filename(filename, defaults["index"])
    job_label = _chain_label(
        filename, defaults["label"], defaults["append_label"]
    )
    logger.info(
        "Building chain pipeline %s from project %s",
        job_label,
        chain_settings.PROJECT_NAME,
    )
    try:
        return build_chain_job(
            chain_settings,
            steps=parsed_steps,
            molecule=molecule,
            label=job_label,
            charge=defaults["charge"],
            multiplicity=defaults["multiplicity"],
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=ctx.obj[_SKIP_COMPLETED_KEY],
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


@click.group(cls=MyGroup, invoke_without_command=True)
@click_chain_options
@click_filename_options
@click_file_label_and_index_options
@click_job_options
@click.option(
    "-s",
    "--steps",
    "steps",
    multiple=True,
    metavar="STEP",
    help=(
        "Pipeline step (repeatable). Use PROGRAM:JOB or quoted "
        "PROGRAM: [OPTIONS] JOB, e.g. gaussian:opt or "
        '"gaussian: -o maxstep=8,maxsize=12 opt".'
    ),
)
@click.pass_context
def chain(
    ctx,
    project,
    filename,
    label,
    append_label,
    index,
    charge,
    multiplicity,
    skip_completed,
    steps,
):
    """Run a custom ``-s/--steps`` pipeline or a workflow subcommand.

    ``chemsmart run/sub chain -p NAME -f mol.xyz -s gaussian:opt`` runs a
    pipeline from repeatable ``-s/--steps``. The ``custom`` subcommand is
    the same pipeline. Predefined workflows are ``pka``, ``fukui``,
    ``redox``, and ``reaction`` and require ``--program`` for submit.
    """
    ctx.ensure_object(dict)
    if project is None:
        chain_settings = None
    else:
        try:
            chain_settings = ChainProjectSettings.from_project(project)
        except FileNotFoundError as exc:
            raise click.UsageError(str(exc)) from exc
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc

    if ctx.invoked_subcommand in CHAIN_WORKFLOW_COMMANDS and steps:
        names = ", ".join(CHAIN_WORKFLOW_COMMANDS)
        raise click.UsageError(
            "Cannot combine -s/--steps with chain workflow subcommands "
            f"({names})."
        )

    _store_chain_invocation(
        ctx,
        chain_settings,
        filename,
        label,
        append_label,
        index,
        charge,
        multiplicity,
        skip_completed,
        steps,
    )

    if ctx.invoked_subcommand is not None:
        return None

    return build_chain_pipeline_from_ctx(ctx)


@click.command("custom", cls=MyCommand)
@click.pass_context
def custom(ctx):
    """Run a custom ``-s/--steps`` pipeline from this chain project.

    Pass ``-s`` on the chain group before ``custom``, e.g.
    ``chain -p combined -f mol.xyz -s gaussian:opt custom``.
    """
    return build_chain_pipeline_from_ctx(ctx)


chain.add_command(custom)
register_chain_workflows(chain)
