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
    validate_reference_options,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


def _build_gaussian_pka_settings(
    proton_index, shared, opt_settings, sp_settings=None
):
    """Build ``GaussianpKaJobSettings`` from CLI and project settings."""
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
    for key, value in shared.items():
        if key in cli_only:
            continue
        pka_kwargs[rename.get(key, key)] = value

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
        key: value
        for key, value in vars(opt_settings).items()
        if key in gs_params and value is not None and key not in pka_kwargs
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
        ctx.invoke(submit, skip_completed=skip_completed)


@pka.command("submit", cls=MyCommand)
@click_job_options
@click.pass_context
def submit(ctx, skip_completed, **kwargs):
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
        list[GaussianpKaJob]: List of jobs when multiple molecules detected.
    """
    shared = ctx.obj["pka_shared"]
    filename = ctx.obj.get("filename")
    proton_index = ctx.obj.get("pka_proton_index")
    color_code = ctx.obj.get("pka_color_code")
    jobrunner = ctx.obj["jobrunner"]

    # ── resolve proton index (CDXML auto-detect) ──
    proton_index, pka_molecules = PKaCDXFile.resolve_proton_index(
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
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob

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
        return [
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
@click.pass_context
def batch(ctx, skip_completed, **kwargs):
    """Table-driven batch pKa job submission.

    Returns a list of pKa jobs created from the input table, one per row.
    """
    shared = ctx.obj["pka_shared"]
    jobrunner = ctx.obj["jobrunner"]

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent Gaussian -f/--filename to specify the table file path."
        )

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob
    from chemsmart.utils.utils import (
        PKaTableEntry,
    )

    try:
        entries = PKaTableEntry.parse_pka_table(input_table_path)
        for entry in entries:
            entry.validate(check_file_exists=True)
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
    original_scheme = shared["scheme"]
    for index, entry in enumerate(entries):
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        molecule = Molecule.from_filepath(filepath)
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

        pka_settings = _build_gaussian_pka_settings(
            proton_index=int(entry.proton_index),
            shared=row_shared,
            opt_settings=row_opt_settings,
            sp_settings=project_settings.sp_settings(),
        )

        job = GaussianpKaJob(
            molecule=molecule,
            settings=pka_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        # Preserve row-level input so submit-script reconstruction can emit
        # one-entry commands instead of replaying the full table.
        job._batch_entry = {
            "filepath": str(filepath),
            "proton_index": int(entry.proton_index),
            "charge": int(entry.charge),
            "multiplicity": int(entry.multiplicity),
            "scheme": row_shared["scheme"],
        }
        jobs.append(job)

    logger.info(f"Created {len(jobs)} pKa jobs from table")
    return jobs


def _create_pka_jobs_from_molecules(
    ctx, pka_molecules, shared, skip_completed, **kwargs
):
    """Create one GaussianpKaJob per PKaMolecule."""
    from chemsmart.jobs.gaussian.pka import GaussianpKaJob

    validate_reference_options(shared)

    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    jobrunner = ctx.obj["jobrunner"]
    filename = ctx.obj.get("filename", "")

    opt_settings = project_settings.opt_settings()
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    base_name = os.path.splitext(os.path.basename(filename))[0] or "pka"

    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        label = f"{base_name}_frag{idx}_pka"
        pka_settings = _build_gaussian_pka_settings(
            pka_mol.proton_index,
            shared,
            opt_settings,
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
    return jobs


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
