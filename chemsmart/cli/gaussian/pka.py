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
from pathlib import Path

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.cli.pka import (
    apply_pka_molecule_charge_multiplicity,
    batch_pka_jobs_from_cdxml,
    click_pka_proton_options,
    click_pka_shared_options,
    is_pka_cdxml_input,
    require_pka_charge_multiplicity,
    resolve_pka_batch_row,
    resolve_pka_submit_proton_options,
    resolve_proton_index,
    validate_reference_options,
    wrap_pka_jobs_in_batch,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.jobs.gaussian.pka import GaussianpKaJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_shared_options
@click_pka_proton_options
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
    pressure,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
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
    from chemsmart.cli.pka import resolve_pka_entropy_cutoff

    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )
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
        pressure=pressure,
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
        skip_completed=skip_completed,
    )
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = shared
    ctx.obj["pka_proton_index"] = proton_index
    ctx.obj["pka_color_code"] = color_code

    if ctx.invoked_subcommand is None:
        from chemsmart.utils.datasets import PKaTableEntry

        if PKaTableEntry.is_submission_table(ctx.obj.get("filename")):
            return ctx.invoke(batch, skip_completed=skip_completed)
        return ctx.invoke(submit, skip_completed=skip_completed)


@pka.command("submit", cls=MyCommand)
@click_job_options
@click_pka_proton_options
@click.pass_context
def submit(ctx, skip_completed, proton_index, color_code, **kwargs):
    """Submit a single-molecule Gaussian pKa calculation.

    Builds Gaussian pKa settings from CLI/project context, validates required
    inputs (charge/multiplicity), and returns the created job(s). If the input
    is a multi-fragment CDXML, this will expand to one job per molecule.

    Example:
        chemsmart sub gaussian <gaussian_options> pka <pka_options> submit

    Args:
        ctx: Click context containing shared options and objects from the
            parent `pka` command (settings, molecules, jobrunner).
        skip_completed (bool): Skip execution for jobs already completed.
        **kwargs: Additional CLI options forwarded to job creation.

    Returns:
        GaussianpKaJob: Single job when one molecule is processed.
        GaussianBatchJob: Batch container when multiple molecules are detected.
    """
    shared = ctx.obj["pka_shared"]
    filename = ctx.obj.get("filename")
    jobrunner = ctx.obj["jobrunner"]

    from chemsmart.utils.datasets import PKaTableEntry

    if PKaTableEntry.is_submission_table(filename):
        return ctx.invoke(batch, skip_completed=skip_completed)

    proton_index, color_code = resolve_pka_submit_proton_options(
        ctx, proton_index=proton_index, color_code=color_code
    )

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

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    molecules = ctx.obj["molecules"]
    if molecules:
        opt_settings = apply_pka_molecule_charge_multiplicity(
            opt_settings, molecules[-1]
        )

    pka_settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
        proton_index, shared, opt_settings, project_settings.sp_settings()
    )
    require_pka_charge_multiplicity(
        pka_settings, source_hint=f"input file {filename}"
    )

    log_pka_settings(pka_settings, proton_index, shared)

    # ── create job(s) ──
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
        return wrap_pka_jobs_in_batch(
            jobs,
            GaussianBatchJob,
            jobrunner,
            label=f"{label}_pka_batch",
        )

    return GaussianpKaJob(
        molecule=molecules[-1],
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )


@pka.command("batch", cls=MyCommand)
@click_job_options
@click_pka_proton_options
@click.pass_context
def batch(ctx, skip_completed, proton_index, color_code, **kwargs):
    """Batch pKa job submission from a CSV table or multi-molecule CDXML.

    CSV tables provide filepath, proton_index, charge, and multiplicity per row.
    CDXML files create one job per fragment using coloured-proton detection.
    """
    shared = ctx.obj["pka_shared"]
    jobrunner = ctx.obj["jobrunner"]

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent Gaussian -f/--filename to specify the table file path."
        )

    if is_pka_cdxml_input(input_table_path):
        return batch_pka_jobs_from_cdxml(
            ctx,
            skip_completed,
            _create_pka_jobs_from_molecules,
            lambda ctx, **invoke_kwargs: ctx.invoke(submit, **invoke_kwargs),
            **kwargs,
        )
    from chemsmart.utils.datasets import PKaOutputTable, PKaTableEntry

    try:
        entries = PKaTableEntry.parse_pka_table(input_table_path)
        # Validate entries using the utility that supports optional file checks
        PKaOutputTable.validate_pka_table_entries(
            entries, check_file_exists=True
        )
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    if shared["scheme"] == "proton exchange":

        missing = []
        if shared["reference"] is None:
            missing.append("-r/--reference")
        elif shared["reference_proton_index"] is None:
            ref = shared["reference"]
            if ref.endswith((".cdx", ".cdxml")):
                try:
                    shared["reference_proton_index"] = (
                        PKaCDXFile.resolve_reference_proton(
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

    _, color_code = resolve_pka_submit_proton_options(
        ctx, proton_index=proton_index, color_code=color_code
    )

    jobs = []
    original_scheme = shared["scheme"]
    for index, entry in enumerate(entries):
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        try:
            row_proton_index, molecule = resolve_pka_batch_row(
                filepath,
                proton_index=entry.proton_index,
                color_code=color_code,
            )
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        label = Path(filepath).stem

        row_opt_settings = copy.copy(opt_settings)
        row_opt_settings.charge = int(entry.charge)
        row_opt_settings.multiplicity = int(entry.multiplicity)

        row_shared = copy.copy(shared)
        row_shared["scheme"] = (
            original_scheme
            if index == 0
            else (
                "direct"
                if original_scheme == "proton exchange"
                else original_scheme
            )
        )
        if row_shared["scheme"] != "proton exchange":
            row_shared["reference"] = None
            row_shared["reference_proton_index"] = None
            row_shared["reference_charge"] = None
            row_shared["reference_multiplicity"] = None
            row_shared["reference_conjugate_base_charge"] = None
            row_shared["reference_conjugate_base_multiplicity"] = None

        pka_settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            proton_index=row_proton_index,
            shared=row_shared,
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

    logger.info(f"Created {len(jobs)} pKa jobs from table")
    table_label = Path(input_table_path).stem or "pka_batch"
    return wrap_pka_jobs_in_batch(
        jobs,
        GaussianBatchJob,
        jobrunner,
        label=f"{table_label}_pka_batch",
    )


def _create_pka_jobs_from_molecules(
    ctx, pka_molecules, shared, skip_completed, **kwargs
):
    """Create one GaussianpKaJob per PKaMolecule."""

    validate_reference_options(shared)

    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    jobrunner = ctx.obj["jobrunner"]
    filename = ctx.obj.get("filename", "")

    opt_settings = project_settings.opt_settings()
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    basename = Path(filename).stem or "pka"

    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        label = f"{basename}_frag{idx}_pka"
        row_opt_settings = apply_pka_molecule_charge_multiplicity(
            opt_settings, pka_mol
        )
        require_pka_charge_multiplicity(
            row_opt_settings,
            source_hint=f"CDXML fragment {idx} in {filename}",
        )
        pka_settings = GaussianpKaJobSettings.build_gaussian_pka_settings(
            pka_mol.proton_index,
            shared,
            row_opt_settings,
            project_settings.sp_settings(),
        )

        logger.info(
            f"Creating pKa job for fragment {idx}: "
            f"proton_index={pka_mol.proton_index}, label={label}"
        )

        jobs.append(
            GaussianpKaJob(
                molecule=pka_mol,
                settings=pka_settings,
                label=label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                **kwargs,
            )
        )

    logger.info(f"Created {len(jobs)} pKa jobs from multi-fragment CDXML")
    basename = Path(filename).stem or "pka"
    return wrap_pka_jobs_in_batch(
        jobs,
        GaussianBatchJob,
        jobrunner,
        label=f"{basename}_pka_batch",
    )


def log_pka_settings(pka_settings, proton_index, shared):
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
