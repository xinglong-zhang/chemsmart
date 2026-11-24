"""Build Gaussian fragment guesses with optional wavefunction stabilization."""

import re
from pathlib import Path

import click
from click.core import ParameterSource

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_irc_options,
    click_gaussian_jobtype_options,
    click_gaussian_stability_option,
    gaussian,
)
from chemsmart.cli.job import click_job_options
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.jobs.gaussian.frag import GaussianFragJob
from chemsmart.jobs.gaussian.settings import GaussianFragJobSettings
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
    update_irc_label,
)
from chemsmart.utils.utils import (
    get_list_from_string_range,
    get_range_from_list,
)


def _parse_fragment_option(value, name):
    """Convert a CLI fragment override to integers."""
    if value is None:
        return None
    value = value.strip()
    if value.startswith("[") and value.endswith("]"):
        value = value[1:-1]
    tokens = value.replace(",", " ").split()
    if not tokens:
        raise click.BadParameter(
            "Expected at least one integer.", param_hint=name
        )
    values = []
    for index, token in enumerate(tokens, start=1):
        try:
            values.append(int(token))
        except ValueError as exc:
            raise click.BadParameter(
                f"Value {index} ({token!r}) must be an integer.",
                param_hint=name,
            ) from exc
    return values


def _parse_fragment_ids(value, num_atoms):
    """Convert ordered atom groups to fragment IDs in the input atom order."""
    if not value.lstrip().startswith("["):
        return _parse_fragment_option(value, "--fragment-ids")

    value = value.strip()
    if not re.fullmatch(
        r"\[\s*\[[^\[\]]*\](?:\s*,\s*\[[^\[\]]*\])*\s*\]", value
    ):
        raise click.BadParameter(
            "Expected atom groups such as '[[1,4-10],[2,3]]'.",
            param_hint="--fragment-ids",
        )
    ids_by_atom = {}
    groups = re.findall(r"\[([^\[\]]*)\]", value)
    for fragment_id, group in enumerate(groups, start=1):
        try:
            atom_indices = get_list_from_string_range(
                group, max_index=num_atoms
            )
        except ValueError as exc:
            raise click.BadParameter(
                f"Fragment {fragment_id} has invalid atom indices. {exc}",
                param_hint="--fragment-ids",
            ) from exc
        for atom_index in atom_indices:
            if atom_index in ids_by_atom:
                raise click.BadParameter(
                    f"Atom {atom_index} is listed more than once.",
                    param_hint="--fragment-ids",
                )
            ids_by_atom[atom_index] = fragment_id
    missing = [
        index for index in range(1, num_atoms + 1) if index not in ids_by_atom
    ]
    if missing:
        missing_indices = ", ".join(get_range_from_list(missing))
        raise click.BadParameter(
            f"Missing atom indices {missing_indices}. "
            f"Assign all {num_atoms} atoms to fragment groups.",
            param_hint="--fragment-ids",
        )
    return [ids_by_atom[index] for index in range(1, num_atoms + 1)]


@gaussian.command("frag", cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click_gaussian_irc_options
@click.option(
    "-fi",
    "--fragment-ids",
    help="Groups of 1-based atom indices, e.g. '[[1,4-10],[2,3]]'.",
)
@click.option(
    "-fc",
    "--fragment-charges",
    help="Charges in fragment order, e.g. '[1,-1]'.",
)
@click.option(
    "-fm",
    "--fragment-multiplicities",
    help="Signed multiplicities in fragment order, e.g. '[2,-2]'.",
)
@click_gaussian_stability_option
@click.pass_context
def frag(
    ctx,
    jobtype,
    fragment_ids,
    fragment_charges,
    fragment_multiplicities,
    stable,
    coordinates,
    step_size,
    num_steps,
    predictor,
    recorrect,
    recalc_step,
    direction,
    maxpoints,
    maxcycles,
    stepsize,
    flat_irc,
    **kwargs,
):
    """Fragment guess -> stability step -> task.

    Read fragments from .com/.gjf/.inp, or supply all three fragment options.
    Every continuation uses the same checkpoint and Guess=Read.

    Supported continuations: opt, ts, modred, scan, sp, irc, td, wbi, nci.
    """
    if jobtype is None:
        raise click.BadParameter(
            "Provide -j/--jobtype, as for a link job.", param_hint="--jobtype"
        )
    jobtype = jobtype.lower()
    if jobtype not in GaussianFragJobSettings.CONTINUATION_TYPES:
        raise click.BadParameter(
            f"Unsupported job type: {jobtype}", param_hint="--jobtype"
        )
    try:
        settings = get_setting_from_jobtype_for_gaussian(
            ctx.obj["project_settings"],
            jobtype,
            coordinates,
            step_size,
            num_steps,
        )
    except (AssertionError, SyntaxError, ValueError) as exc:
        raise click.UsageError(str(exc)) from exc
    if settings is None:
        raise click.UsageError(
            f"The selected Gaussian project has no settings for '{jobtype}'."
        )
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    input_ids = None
    charges = _parse_fragment_option(fragment_charges, "--fragment-charges")
    multiplicities = _parse_fragment_option(
        fragment_multiplicities, "--fragment-multiplicities"
    )
    filename = ctx.obj.get("filename")
    source = None
    if filename and Path(filename).suffix.lower() in {".com", ".gjf", ".inp"}:
        source = Gaussian16Input(filename)
        try:
            if fragment_ids is None:
                input_ids = source.fragment_ids
            if charges is None or multiplicities is None:
                pairs = source.fragment_charge_multiplicities
                if pairs:
                    input_charges, input_multiplicities = zip(*pairs)
                    if charges is None:
                        charges = input_charges
                    if multiplicities is None:
                        multiplicities = input_multiplicities
        except ValueError as exc:
            raise click.ClickException(str(exc)) from exc
    settings_kwargs = settings.__dict__.copy()
    if jobtype == "irc":
        irc_options = dict(
            predictor=predictor,
            recorrect=recorrect,
            recalc_step=recalc_step,
            direction=direction,
            maxpoints=maxpoints,
            maxcycles=maxcycles,
            stepsize=stepsize,
            flat_irc=flat_irc,
        )
        # A replayed Click default would override the IRC project setting.
        recorded_options = ctx.obj["subcommand"][-1]["kwargs"]
        for name, value in irc_options.items():
            if ctx.get_parameter_source(name) == ParameterSource.DEFAULT:
                recorded_options.pop(name, None)
            else:
                settings_kwargs[name] = value
    molecules = ctx.obj["molecules"]
    if not molecules:
        raise click.UsageError("Select at least one structure for frag.")
    indices = ctx.obj.get("molecule_indices")
    selected = (
        list(zip(molecules, indices))
        if indices is not None
        else [(molecules[-1], None)]
    )
    jobs = []
    for molecule, index in selected:
        prefix = f"Structure {index}. " if index is not None else ""
        try:
            if (
                source is not None
                and molecule.frozen_atoms is not None
                and 0 < len(molecule.frozen_atoms) < molecule.num_atoms
            ):
                molecule = molecule.copy()
                molecule.frozen_atoms = (
                    source.coordinate_block.atom_aligned_constraints
                )
            ids = (
                input_ids
                if fragment_ids is None
                else _parse_fragment_ids(fragment_ids, molecule.num_atoms)
            )
            frag_settings = GaussianFragJobSettings(
                fragment_ids=ids,
                fragment_charges=charges,
                fragment_multiplicities=multiplicities,
                stable=stable,
                **settings_kwargs,
            )
            label = ctx.obj["label"]
            index_suffix = f"_idx{index}" if index is not None else ""
            if jobtype == "irc":
                label = update_irc_label(
                    label=f"{label}{index_suffix}_irc",
                    direction=frag_settings.direction,
                    flat_irc=frag_settings.flat_irc
                    and frag_settings.direction is not None,
                )
                label += "_frag"
            else:
                label = f"{label}_{jobtype}_frag{index_suffix}"
            jobs.append(
                GaussianFragJob(
                    molecule=molecule,
                    settings=frag_settings,
                    label=label,
                    jobrunner=ctx.obj["jobrunner"],
                    submission_index=index,
                    **kwargs,
                )
            )
        except click.BadParameter as exc:
            exc.message = prefix + exc.message
            raise
        except ValueError as exc:
            raise click.ClickException(f"{prefix}{exc}") from exc
    return jobs if len(jobs) > 1 else jobs[0]
