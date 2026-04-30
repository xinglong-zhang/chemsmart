"""
CLI for Gaussian pKa input generation (job submission).

Subcommands
-----------
submit         Submit single-molecule (or multi-fragment CDXML) pKa jobs.
batch          Table-driven batch job submission.

When ``pka`` is invoked without an explicit subcommand the ``submit``
path is executed automatically for backward compatibility.

Output analysis (analyze, batch-analyze, thermo) lives in the
backend-independent ``chemsmart run pka`` command.
"""

import logging
import os
from pathlib import Path

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.cli.pka import (
    click_pka_proton_options,
    click_pka_shared_options,
    click_pka_submit_options,
    resolve_proton_index,
    validate_reference_options,
)
from chemsmart.jobs.runner import get_serial_mode
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_shared_options
@click_pka_proton_options
@click_pka_submit_options
@click.pass_context
def pka(
    ctx,
    skip_completed,
    scheme,
    reference,
    reference_proton_index,
    reference_color_code,
    reference_charge,
    reference_multiplicity,
    reference_conjugate_base_charge,
    reference_conjugate_base_multiplicity,
    delta_g_proton,
    conjugate_base_charge,
    conjugate_base_multiplicity,
    solvent_model,
    solvent_id,
    temperature,
    concentration,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    proton_index,
    color_code,
    **kwargs,
):
    """Gaussian pKa job submission.

    For output analysis (backend-independent), use:
      chemsmart run pka analyze ...
      chemsmart run pka batch-analyze ...
    """
    shared = dict(
        scheme=scheme,
        reference=reference,
        reference_proton_index=reference_proton_index,
        reference_color_code=reference_color_code,
        reference_charge=reference_charge,
        reference_multiplicity=reference_multiplicity,
        reference_conjugate_base_charge=reference_conjugate_base_charge,
        reference_conjugate_base_multiplicity=reference_conjugate_base_multiplicity,
        delta_g_proton=delta_g_proton,
        conjugate_base_charge=conjugate_base_charge,
        conjugate_base_multiplicity=conjugate_base_multiplicity,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        temperature=temperature,
        concentration=concentration,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        skip_completed=skip_completed,
    )
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = shared
    ctx.obj["pka_proton_index"] = proton_index
    ctx.obj["pka_color_code"] = color_code

    if ctx.invoked_subcommand is None:
        ctx.invoke(submit, skip_completed=skip_completed)


@pka.command("submit", cls=MyCommand)
@click_job_options
@click_pka_submit_options
@click.pass_context
def submit(ctx, skip_completed, **kwargs):
    """Submit a single-molecule Gaussian pKa calculation."""
    shared = ctx.obj["pka_shared"]
    filename = ctx.obj.get("filename")
    proton_index = ctx.obj.get("pka_proton_index")
    color_code = ctx.obj.get("pka_color_code")
    jobrunner = ctx.obj["jobrunner"]
    serial_mode = get_serial_mode(jobrunner)

    # ── resolve proton index (CDXML auto-detect) ──
    proton_index, pka_molecules = resolve_proton_index(
        filename, proton_index, color_code
    )

    # ── multi-fragment CDXML -> one job per molecule ──
    if pka_molecules is not None:
        return _create_pka_jobs_from_molecules(
            ctx, pka_molecules, shared, skip_completed, **kwargs
        )

    # ── validate reference acid ──
    validate_reference_options(shared)

    # ── build settings ──
    from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob, GaussianpKaJob

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    pka_settings = _build_gaussian_pka_settings(
        proton_index, shared, opt_settings, project_settings.sp_settings()
    )

    if pka_settings.charge is None:
        raise click.UsageError(
            "Charge must be specified via -c/--charge option"
        )
    if pka_settings.multiplicity is None:
        raise click.UsageError(
            "Multiplicity must be specified via -m/--multiplicity option"
        )

    _log_pka_settings(pka_settings, proton_index, shared)

    # ── create job(s) ──
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj.get("molecule_indices")
    label = ctx.obj["label"]
    jobrunner = ctx.obj["jobrunner"]

    if len(molecules) > 1 and molecule_indices:
        logger.info(f"Creating {len(molecules)} pKa jobs")
        jobs = [
            GaussianpKaJob(
                molecule=mol,
                settings=pka_settings,
                label=f"{label}_idx{idx}",
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
            for mol, idx in zip(molecules, molecule_indices)
        ]
        # Batch execution policy is derived from the shared jobrunner mode.
        return GaussianpKaBatchJob(
            jobs=jobs,
            run_in_serial=serial_mode.run_in_serial,
            jobrunner=jobrunner,
        )

    job = GaussianpKaJob(
        molecule=molecules[-1],
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
    # Batch execution policy is derived from the shared jobrunner mode.
    return GaussianpKaBatchJob(
        jobs=[job],
        run_in_serial=serial_mode.run_in_serial,
        jobrunner=jobrunner,
    )


@pka.command("batch", cls=MyCommand)
@click_job_options
@click.pass_context
def batch(ctx, skip_completed, **kwargs):
    """Table-driven batch pKa job submission."""
    shared = ctx.obj["pka_shared"]
    jobrunner = ctx.obj["jobrunner"]
    serial_mode = get_serial_mode(jobrunner)

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent Gaussian -f/--filename to specify the table file path."
        )

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob, GaussianpKaJob
    from chemsmart.utils.utils import (
        parse_pka_table,
        validate_pka_table_entries,
    )

    try:
        entries = parse_pka_table(input_table_path)
        validate_pka_table_entries(entries, check_file_exists=True)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    if shared["scheme"] == "proton exchange":
        from chemsmart.cli.pka import resolve_reference_proton

        missing = []
        if shared["reference"] is None:
            missing.append("-r/--reference")
        elif shared["reference_proton_index"] is None:
            ref = shared["reference"]
            if ref.endswith((".cdx", ".cdxml")):
                try:
                    shared["reference_proton_index"] = (
                        resolve_reference_proton(
                            ref,
                            None,
                            shared["reference_color_code"],
                        )
                    )
                except click.UsageError:
                    missing.append("-rpi/--reference-proton-index")
            else:
                missing.append("-rpi/--reference-proton-index")
        if shared["reference_charge"] is None:
            missing.append("-rc/--reference-charge")
        if shared["reference_multiplicity"] is None:
            missing.append("-rm/--reference-multiplicity")
        if missing:
            raise click.UsageError(
                "For proton exchange cycle with batch input, these reference acid options are required:\n  "
                + "\n  ".join(missing)
            )

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj.get("job_settings")
    kw = ctx.obj.get("keywords", {})
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=kw)

    import copy

    jobs = []
    for entry in entries:
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        molecule = Molecule.from_filepath(filepath)
        label = Path(filepath).stem

        row_opt_settings = copy.copy(opt_settings)
        row_opt_settings.charge = int(entry.charge)
        row_opt_settings.multiplicity = int(entry.multiplicity)

        pka_settings = _build_gaussian_pka_settings(
            proton_index=int(entry.proton_index),
            shared=shared,
            opt_settings=row_opt_settings,
            sp_settings=project_settings.sp_settings(),
        )

        jobs.append(
            GaussianpKaJob(
                molecule=molecule,
                settings=pka_settings,
                label=label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
        )

    return GaussianpKaBatchJob(
        jobs=jobs,
        run_in_serial=serial_mode.run_in_serial,
        write_outcome_logs=True,
        jobrunner=jobrunner,
    )


def _build_gaussian_pka_settings(
    proton_index, shared, opt_settings, sp_settings=None
):
    """Build a GaussianpKaJobSettings from shared options and project."""
    import inspect

    from chemsmart.jobs.gaussian.settings import (
        GaussianJobSettings,
        GaussianpKaJobSettings,
    )

    cli_only = {"reference_color_code", "skip_completed"}
    rename = {
        "reference": "reference_file",
        "delta_g_proton": "delta_G_proton",
    }

    pka_kwargs = {}
    for k, v in shared.items():
        if k in cli_only:
            continue
        pka_kwargs[rename.get(k, k)] = v

    gs_params = {
        name
        for name, param in inspect.signature(
            GaussianJobSettings.__init__
        ).parameters.items()
        if name != "self"
        and param.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
    }
    opt_kwargs = {
        k: v
        for k, v in vars(opt_settings).items()
        if k in gs_params and v is not None and k not in pka_kwargs
    }

    def _first_non_none(*values):
        for val in values:
            if val is not None:
                return val
        return None

    solvent_model = _first_non_none(
        pka_kwargs.get("solvent_model"),
        getattr(opt_settings, "solvent_model", None),
        getattr(sp_settings, "solvent_model", None) if sp_settings else None,
        "SMD",
    )
    solvent_id = _first_non_none(
        pka_kwargs.get("solvent_id"),
        getattr(opt_settings, "solvent_id", None),
        getattr(sp_settings, "solvent_id", None) if sp_settings else None,
        "water",
    )

    pka_kwargs["solvent_model"] = solvent_model
    pka_kwargs["solvent_id"] = solvent_id

    return GaussianpKaJobSettings(
        proton_index=proton_index,
        **pka_kwargs,
        **opt_kwargs,
    )


def _create_pka_jobs_from_molecules(
    ctx, pka_molecules, shared, skip_completed, **kwargs
):
    """Create one GaussianpKaJob per PKaMolecule."""
    from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob, GaussianpKaJob

    validate_reference_options(shared)

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    jobrunner = ctx.obj["jobrunner"]
    serial_mode = get_serial_mode(jobrunner)

    filename = ctx.obj.get("filename", "")
    base_name = os.path.splitext(os.path.basename(filename))[0]

    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        mol_label = f"{base_name}_frag{idx}_pka"
        pka_settings = _build_gaussian_pka_settings(
            pka_mol.proton_index,
            shared,
            opt_settings,
            project_settings.sp_settings(),
        )

        jobs.append(
            GaussianpKaJob(
                molecule=pka_mol,
                settings=pka_settings,
                label=mol_label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
        )

    return GaussianpKaBatchJob(
        jobs=jobs,
        run_in_serial=serial_mode.run_in_serial,
        jobrunner=jobrunner,
    )


def _log_pka_settings(pka_settings, proton_index, shared):
    """Emit informational log messages for a pKa job."""
    logger.info(f"Proton index to remove: {proton_index}")
    logger.info(f"Thermodynamic cycle: {shared['scheme']}")
    logger.info(
        f"Protonated form (HA): charge={pka_settings.charge}, "
        f"mult={pka_settings.multiplicity}"
    )

    cb_charge = (
        shared["conjugate_base_charge"]
        if shared["conjugate_base_charge"] is not None
        else pka_settings.charge - 1
    )
    cb_mult = (
        shared["conjugate_base_multiplicity"]
        if shared["conjugate_base_multiplicity"] is not None
        else pka_settings.multiplicity
    )
    logger.info(f"Conjugate base (A-): charge={cb_charge}, mult={cb_mult}")
