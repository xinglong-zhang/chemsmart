"""Shared redox CLI options, job builder, and backend-independent analysis.

Job submission lives under Gaussian, ORCA, or chain and uses
``register_redox_cli`` / ``build_redox_job``. Output analysis is
``chemsmart run redox analyze`` (also ``chemsmart run chain redox analyze``).
"""

from __future__ import annotations

import functools
import inspect
import logging
import os
from pathlib import Path

import click

from chemsmart.analysis.redox import (  # noqa: F401
    RedoxReference,
    compute_redox_potential,
    format_redox_summary,
    get_redox_reference,
    list_redox_references,
    register_redox_reference,
    resolve_redox_reference,
)
from chemsmart.cli.job import click_job_options
from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry_cutoff_options,
    thermochemistry_temp_pressure_conc_options,
)
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import get_program_type_from_file
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


def click_redox_reference_options(f):
    """Reference couple and electron-count options for redox submit."""

    @click.option(
        "-n",
        "--n-electrons",
        type=int,
        default=None,
        help=(
            "Electrons transferred. Defaults to the reference couple. "
            "Must match the reference couple when given."
        ),
    )
    @click.option(
        "-r",
        "--reference",
        type=str,
        default="fc_fc+",
        show_default=True,
        help="Registry name of the reference couple.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_redox_submit_structure_options(f):
    """Target/reference geometry and charge options for redox submit."""

    @click.option(
        "-rrm",
        "--ref-red-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reduced reference. Defaults to the "
            "registry value or the reduced-form multiplicity of Ref_ox."
        ),
    )
    @click.option(
        "-rrc",
        "--ref-red-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reduced reference. Defaults to the registry "
            "value or ox − n of the oxidized reference."
        ),
    )
    @click.option(
        "-rom",
        "--ref-ox-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the oxidized reference. Defaults to the "
            "registry value or the loaded geometry."
        ),
    )
    @click.option(
        "-roc",
        "--ref-ox-charge",
        type=int,
        default=None,
        help=(
            "Charge of the oxidized reference. Defaults to the registry "
            "value or the loaded geometry."
        ),
    )
    @click.option(
        "-rr",
        "--ref-red",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Reduced reference geometry. Defaults to the oxidized reference "
            "structure from --ref-ox (or the registry ox_file) with charge "
            "ox − n."
        ),
    )
    @click.option(
        "-ro",
        "--ref-ox",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Oxidized reference geometry. Required unless the registered "
            "couple provides ox_file. The reduced reference defaults to "
            "this geometry with charge ox − n."
        ),
    )
    @click.option(
        "-rdm",
        "--red-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reduced target. Defaults to the "
            "reduced-form multiplicity of the oxidized target."
        ),
    )
    @click.option(
        "-rdc",
        "--red-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reduced target. Defaults to ox − n from the "
            "parent oxidized charge."
        ),
    )
    @click.option(
        "-rd",
        "--red",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Reduced target geometry. Defaults to the oxidized structure "
            "from parent -f with charge ox − n."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def _click_thermochemistry_options(f):
    """T/P/c and quasi-RRHO cutoffs accepted by ``Thermochemistry``."""
    f = thermochemistry_temp_pressure_conc_options(
        f,
        temperature_required=False,
        temperature_default=298.15,
        concentration_default=1.0,
        pressure_default=1.0,
    )
    return thermochemistry_cutoff_options(f)


def click_redox_shared_options(f):
    """Submit options shared by Gaussian, ORCA, and chain redox commands."""
    f = click_redox_submit_structure_options(f)
    f = click_redox_reference_options(f)
    return _click_thermochemistry_options(f)


def click_redox_analyze_options(f):
    """Output-file options for exchange redox analysis."""
    f = click.option(
        "-o",
        "--output",
        default=None,
        type=str,
        help="Write redox results to this file instead of printing them.",
    )(f)
    f = click.option(
        "-n",
        "--n-electrons",
        type=int,
        default=1,
        show_default=True,
        help="Electrons transferred in the exchange.",
    )(f)
    f = click.option(
        "-er",
        "--e-ref",
        type=float,
        required=True,
        help=(
            "Reference reduction potential in volts. Used as E_ref in "
            "E_target = E_ref − ΔG_exchange / (n F)."
        ),
    )(f)
    f = click.option(
        "-rrs",
        "--ref-red-solv",
        "ref_red_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Solution-phase SP output for the reduced reference.",
    )(f)
    f = click.option(
        "-ros",
        "--ref-ox-solv",
        "ref_ox_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Solution-phase SP output for the oxidized reference.",
    )(f)
    f = click.option(
        "-rds",
        "--red-solv",
        "red_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Solution-phase SP output for the reduced target.",
    )(f)
    f = click.option(
        "-oxs",
        "--ox-solv",
        "ox_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Solution-phase SP output for the oxidized target.",
    )(f)
    f = click.option(
        "-rrg",
        "--ref-red-gas",
        "ref_red_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Gas-phase opt+freq output for the reduced reference.",
    )(f)
    f = click.option(
        "-rog",
        "--ref-ox-gas",
        "ref_ox_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Gas-phase opt+freq output for the oxidized reference.",
    )(f)
    f = click.option(
        "-rdg",
        "--red-gas",
        "red_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Gas-phase opt+freq output for the reduced target.",
    )(f)
    f = click.option(
        "-oxg",
        "--ox-gas",
        "ox_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help="Gas-phase opt+freq output for the oxidized target.",
    )(f)
    return f


def _auto_discover_redox_files(ox_gas_path, ref_ox_gas_path, program=None):
    """Infer companion output paths from Ox and Ref_ox gas-phase paths."""
    from chemsmart.utils.datasets import (
        REDOX_REFERENCE_SUFFIX_HELP,
        REDOX_TARGET_SUFFIX_HELP,
        discover_redox_reference_companion_outputs,
        discover_redox_target_companion_outputs,
    )

    if program is None:
        program = get_program_type_from_file(ox_gas_path)

    results = discover_redox_target_companion_outputs(
        ox_gas_path, program=program
    )
    results.update(discover_redox_reference_companion_outputs(ref_ox_gas_path))

    missing = [
        f"  {key}: {path}"
        for key, path in results.items()
        if not os.path.isfile(path)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            + REDOX_TARGET_SUFFIX_HELP
            + "\n"
            + REDOX_REFERENCE_SUFFIX_HELP
        )
    return results


def _thermochemistry_kwargs_from_shared(shared):
    keys = (
        "temperature",
        "concentration",
        "pressure",
        "cutoff_entropy_grimme",
        "cutoff_enthalpy",
        "entropy_method",
    )
    return {key: shared[key] for key in keys if shared.get(key) is not None}


def store_redox_shared(ctx, kwargs):
    """Record redox CLI options on ``ctx.obj`` for submit and analyze."""
    s_freq_cutoff, entropy_method = resolve_entropy_cutoff(
        kwargs.get("cutoff_entropy_grimme"),
        kwargs.get("cutoff_entropy_truhlar"),
    )
    ctx.ensure_object(dict)
    ctx.obj["redox_shared"] = dict(
        reference=kwargs.get("reference"),
        n_electrons=kwargs.get("n_electrons"),
        red_file=kwargs.get("red"),
        red_charge=kwargs.get("red_charge"),
        red_multiplicity=kwargs.get("red_multiplicity"),
        ref_ox_file=kwargs.get("ref_ox"),
        ref_red_file=kwargs.get("ref_red"),
        ref_ox_charge=kwargs.get("ref_ox_charge"),
        ref_ox_multiplicity=kwargs.get("ref_ox_multiplicity"),
        ref_red_charge=kwargs.get("ref_red_charge"),
        ref_red_multiplicity=kwargs.get("ref_red_multiplicity"),
        **_thermochemistry_kwargs_from_shared(
            {
                "temperature": kwargs.get("temperature"),
                "concentration": kwargs.get("concentration"),
                "pressure": kwargs.get("pressure"),
                "cutoff_entropy_grimme": s_freq_cutoff,
                "cutoff_enthalpy": kwargs.get("cutoff_enthalpy"),
                "entropy_method": entropy_method,
            }
        ),
    )


def _base_settings_kwargs(opt_settings, settings_cls):
    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.jobs.orca.settings import ORCAJobSettings

    if issubclass(settings_cls, GaussianJobSettings):
        base_cls = GaussianJobSettings
    elif issubclass(settings_cls, ORCAJobSettings):
        base_cls = ORCAJobSettings
    else:
        raise TypeError(
            "settings_cls must be a Gaussian or ORCA job settings class, "
            f"got {settings_cls!r}."
        )
    params = {
        name
        for name, param in inspect.signature(
            base_cls.__init__
        ).parameters.items()
        if name != "self"
        and param.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
    }
    return {
        key: value
        for key, value in vars(opt_settings).items()
        if key in params and value is not None
    }


def _first_non_none(*values):
    for value in values:
        if value is not None:
            return value
    return None


def build_redox_job(ctx, job_cls, settings_cls, skip_completed, **kwargs):
    """Create one redox job from parent ``-f`` and shared redox options."""
    if kwargs:
        store_redox_shared(ctx, kwargs)
    shared = ctx.obj["redox_shared"]

    molecules = ctx.obj.get("molecules")
    if not molecules:
        raise click.UsageError(
            "Redox submit requires a structure file via parent -f/--filename."
        )

    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = project_settings.opt_settings().merge(
        job_settings, keywords=keywords
    )
    sp_settings = project_settings.sp_settings()
    settings_kwargs = _base_settings_kwargs(opt_settings, settings_cls)
    settings_kwargs["solvent_model"] = _first_non_none(
        settings_kwargs.get("solvent_model"),
        opt_settings.solvent_model,
        sp_settings.solvent_model,
    )
    settings_kwargs["solvent_id"] = _first_non_none(
        settings_kwargs.get("solvent_id"),
        opt_settings.solvent_id,
        sp_settings.solvent_id,
    )
    settings_kwargs.update(
        n_electrons=shared["n_electrons"],
        red_file=shared["red_file"],
        red_charge=shared["red_charge"],
        red_multiplicity=shared["red_multiplicity"],
        reference=shared["reference"],
        ref_ox_file=shared["ref_ox_file"],
        ref_red_file=shared["ref_red_file"],
        ref_ox_charge=shared["ref_ox_charge"],
        ref_ox_multiplicity=shared["ref_ox_multiplicity"],
        ref_red_charge=shared["ref_red_charge"],
        ref_red_multiplicity=shared["ref_red_multiplicity"],
        **_thermochemistry_kwargs_from_shared(shared),
    )
    settings_kwargs = {
        key: value
        for key, value in settings_kwargs.items()
        if value is not None
    }

    label = ctx.obj["label"]
    logger.info("Creating %s job %s", job_cls.__name__, label)
    try:
        settings = settings_cls(**settings_kwargs)
        check_charge_and_multiplicity(settings)
        return job_cls(
            molecule=molecules[-1],
            settings=settings,
            label=label,
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=skip_completed,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def register_redox_cli(parent_group, job_cls, settings_cls):
    """Attach ``redox`` submit to a Gaussian or ORCA Click group."""

    @parent_group.command("redox", cls=MyCommand)
    @click_job_options
    @click_redox_shared_options
    @click.pass_context
    def redox(ctx, skip_completed, **kwargs):
        """Submit a dual-level exchange redox calculation.

        Oxidized target comes from parent ``-f``; reduced target uses the
        same geometry (or ``--red``) with charge ``ox − n``. The oxidized
        reference comes from ``--ref-ox`` (or registry ``ox_file``); the
        reduced reference uses the same geometry (or ``--ref-red``) with
        charge ``ox − n``.

        Analyze completed outputs with ``chemsmart run redox analyze``.
        """
        return build_redox_job(
            ctx, job_cls, settings_cls, skip_completed, **kwargs
        )

    return redox


@click.group(name="redox", cls=MyGroup)
@click.pass_context
def redox(ctx):
    """Backend-independent redox output analysis.

    Submit jobs with ``chemsmart run/sub gaussian|orca … redox`` or
    ``chemsmart run/sub chain … redox --program {gaussian,orca}``.
    Analyze requires ``--e-ref`` (reference potential in volts).
    """
    ctx.ensure_object(dict)


@redox.command("analyze", cls=MyCommand)
@click_redox_analyze_options
@_click_thermochemistry_options
@click.pass_context
def analyze(
    ctx,
    ox_gas_file,
    red_gas_file,
    ref_ox_gas_file,
    ref_red_gas_file,
    ox_solv_file,
    red_solv_file,
    ref_ox_solv_file,
    ref_red_solv_file,
    output,
    e_ref,
    n_electrons,
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
):
    """Compute the exchange redox potential from existing output files.

    Dual-level solution free energies: ``G_soln = E_solv + (qh-G − E_gas)``.

    Reaction: ``Ox + Ref_red → Red + Ref_ox``

    ``E_target = E_ref − ΔG_exchange / (n F)``

    ``--e-ref`` is the reference potential in volts.

    Only ``--ox-gas`` and ``--ref-ox-gas`` are required. The remaining six
    files are auto-discovered from the CHEMSMART redox naming convention:
      <basename>_redox_red_opt.<ext>   reduced target gas-phase
      <basename>_redox_ox_sp.<ext>      oxidized target solvent SP
      <basename>_redox_red_sp.<ext>     reduced target solvent SP
      (and the corresponding ``_redox_Ref*`` files for the reference couple)
    Override any auto-discovered path with the corresponding flag.
    """
    if ox_gas_file is None:
        raise click.UsageError("--ox-gas is required for redox analysis.")
    if ref_ox_gas_file is None:
        raise click.UsageError("--ref-ox-gas is required for redox analysis.")

    optional = {
        "red_gas_file": red_gas_file,
        "ox_solv_file": ox_solv_file,
        "red_solv_file": red_solv_file,
        "ref_red_gas_file": ref_red_gas_file,
        "ref_ox_solv_file": ref_ox_solv_file,
        "ref_red_solv_file": ref_red_solv_file,
    }
    if any(value is None for value in optional.values()):
        discovered = _auto_discover_redox_files(ox_gas_file, ref_ox_gas_file)
        discovery_map = {
            "red_gas_file": "red_gas",
            "ox_solv_file": "ox_solv",
            "red_solv_file": "red_solv",
            "ref_red_gas_file": "ref_red_gas",
            "ref_ox_solv_file": "ref_ox_solv",
            "ref_red_solv_file": "ref_red_solv",
        }
        for param, key in discovery_map.items():
            if optional[param] is None:
                optional[param] = discovered[key]
                logger.info("Auto-discovered %s: %s", param, discovered[key])
        red_gas_file = optional["red_gas_file"]
        ox_solv_file = optional["ox_solv_file"]
        red_solv_file = optional["red_solv_file"]
        ref_red_gas_file = optional["ref_red_gas_file"]
        ref_ox_solv_file = optional["ref_ox_solv_file"]
        ref_red_solv_file = optional["ref_red_solv_file"]

    s_freq_cutoff, entropy_method = resolve_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )
    try:
        result = compute_redox_potential(
            ox_gas_file=ox_gas_file,
            red_gas_file=red_gas_file,
            ref_ox_gas_file=ref_ox_gas_file,
            ref_red_gas_file=ref_red_gas_file,
            ox_solv_file=ox_solv_file,
            red_solv_file=red_solv_file,
            ref_ox_solv_file=ref_ox_solv_file,
            ref_red_solv_file=ref_red_solv_file,
            e_ref=e_ref,
            n_electrons=n_electrons,
            **_thermochemistry_kwargs_from_shared(
                {
                    "temperature": temperature,
                    "concentration": concentration,
                    "pressure": pressure,
                    "cutoff_entropy_grimme": s_freq_cutoff,
                    "cutoff_enthalpy": cutoff_enthalpy,
                    "entropy_method": entropy_method,
                }
            ),
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    summary = format_redox_summary(result)
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(summary + "\n")
    else:
        click.echo(summary)
    return None
