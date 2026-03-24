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
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# pka group – holds shared options, dispatches to subcommands
# ═══════════════════════════════════════════════════════════════════════


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_shared_options
@click_pka_proton_options
@click_pka_submit_options
@click.pass_context
def pka(
    ctx,
    # job options
    skip_completed,
    # shared pka options
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
    # proton identification options (pka-layer only)
    proton_index,
    color_code,
    # submit options
    **kwargs,
):
    """Gaussian pKa job submission.

    \b
    Subcommands:
      submit         Single-molecule job submission (default).
      batch          Table-driven batch submission.

    When invoked without a subcommand the ``submit`` path runs
    automatically, so the following two invocations are equivalent:

    \b
      chemsmart run gaussian -f acid.xyz -c 0 -m 1 pka -pi 10 ...
      chemsmart run gaussian -f acid.xyz -c 0 -m 1 pka submit -pi 10 ...

    \b
    For output analysis (backend-independent):
      chemsmart run pka analyze ...
      chemsmart run pka batch-analyze ...
      chemsmart run pka thermo ...

    \b
    Thermodynamic cycles:
      proton exchange (default): HA + Ref- -> A- + HRef
      direct: uses absolute free energy of H+ in water
    """
    # Store shared options for subcommands
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

    # Default to ``submit`` when no subcommand is given
    if ctx.invoked_subcommand is None:
        ctx.invoke(
            submit,
            skip_completed=skip_completed,
        )


# ═══════════════════════════════════════════════════════════════════════
# pka submit
# ═══════════════════════════════════════════════════════════════════════


@pka.command("submit", cls=MyCommand)
@click_job_options
@click_pka_submit_options
@click.pass_context
def submit(ctx, skip_completed, **kwargs):
    """Submit a single-molecule Gaussian pKa calculation.

    \b
    Workflow:
      1. Gas-phase optimization + frequency for HA and A-.
      2. Solution-phase single-point for HA and A- at the SAME level.

    \b
    Examples:
      # Proton exchange cycle
      chemsmart run gaussian -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -s "proton exchange" -r ref.xyz -rpi 1 -rc 0 -rm 1 submit

      # Direct cycle
      chemsmart run gaussian -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -s direct submit

      # Auto-detect proton from ChemDraw file
      chemsmart run gaussian -f phenol.cdxml -c 0 -m 1 pka submit
    """
    shared = ctx.obj["pka_shared"]
    filename = ctx.obj.get("filename")
    proton_index = ctx.obj.get("pka_proton_index")
    color_code = ctx.obj.get("pka_color_code")
    jobrunner = ctx.obj["jobrunner"]

    # Align pKa parallel execution with global --run-in-serial flag
    parallel = not jobrunner.run_in_serial

    # ── resolve proton index (CDXML auto-detect) ──
    proton_index, pka_molecules = resolve_proton_index(
        filename, proton_index, color_code
    )

    # ── multi-fragment CDXML -> one job per molecule ──
    if pka_molecules is not None:
        return _create_pka_jobs_from_molecules(
            ctx, pka_molecules, shared, skip_completed, parallel, **kwargs
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
                parallel=parallel,
                **kwargs,
            )
            for mol, idx in zip(molecules, molecule_indices)
        ]
        return GaussianpKaBatchJob(
            jobs=jobs,
            run_in_serial=jobrunner.run_in_serial,
            jobrunner=jobrunner,
        )

    job = GaussianpKaJob(
        molecule=molecules[-1],
        settings=pka_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        parallel=parallel,
        **kwargs,
    )
    return GaussianpKaBatchJob(
        jobs=[job], run_in_serial=jobrunner.run_in_serial, jobrunner=jobrunner
    )


# ═══════════════════════════════════════════════════════════════════════
# pka batch
# ═══════════════════════════════════════════════════════════════════════


@pka.command("batch", cls=MyCommand)
@click_job_options
@click.pass_context
def batch(ctx, skip_completed, **kwargs):
    """Table-driven batch pKa job submission.

    The table file path is taken from the parent Gaussian ``-f/--filename``
    option.

    \b
    Table format (4 columns, whitespace or comma-delimited):
        filepath    proton_index    charge    multiplicity

    For proton exchange cycle the reference acid options
    (``-r``, ``-rpi``, ``-rc``, ``-rm``) on the parent ``pka`` group
    are still required.

    \b
    Examples:
      chemsmart run gaussian -p myproject -f molecules.txt pka \\
          -s "proton exchange" -r ref.xyz -rpi 5 -rc 0 -rm 1 batch

      chemsmart run gaussian -p myproject -f molecules.csv pka \\
          -s direct batch
    """
    shared = ctx.obj["pka_shared"]
    jobrunner = ctx.obj["jobrunner"]

    # Align pKa parallel execution with global --run-in-serial flag
    parallel = not jobrunner.run_in_serial

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent Gaussian -f/--filename to "
            "specify the table file path."
        )

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob, GaussianpKaJob
    from chemsmart.utils.utils import (
        parse_pka_table,
        validate_pka_table_entries,
    )

    logger.info(f"Reading pKa jobs from table: {input_table_path}")
    try:
        entries = parse_pka_table(input_table_path)
        validate_pka_table_entries(entries, check_file_exists=True)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(entries)} entries in table")

    # Validate reference acid for proton exchange
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
                "For proton exchange cycle with batch input, these "
                "reference acid options are required:\n  "
                + "\n  ".join(missing)
            )

    # Build jobs
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj.get("job_settings")
    kw = ctx.obj.get("keywords", {})
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=kw)
    jobrunner = ctx.obj[
        "jobrunner"
    ]  # Re-assigned but fine, can also remove duplicates if preferred, but existing code has it. Let's not remove it to minimize changes unless nearby.

    import copy

    jobs = []
    for entry in entries:
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        molecule = Molecule.from_filepath(filepath)
        label = Path(filepath).stem

        # Per-row charge and multiplicity override the opt_settings defaults.
        row_opt_settings = copy.copy(opt_settings)
        row_opt_settings.charge = int(entry.charge)
        row_opt_settings.multiplicity = int(entry.multiplicity)

        pka_settings = _build_gaussian_pka_settings(
            proton_index=int(entry.proton_index),
            shared=shared,
            opt_settings=row_opt_settings,
            sp_settings=project_settings.sp_settings(),
        )

        logger.info(
            f"Creating pKa job for {filepath}: "
            f"proton_index={entry.proton_index}, "
            f"charge={entry.charge}, mult={entry.multiplicity}"
        )

        job = GaussianpKaJob(
            molecule=molecule,
            settings=pka_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            parallel=parallel,
            **kwargs,
        )

        jobs.append(job)

    logger.info(f"Created {len(jobs)} pKa jobs from table")
    return GaussianpKaBatchJob(
        jobs=jobs, run_in_serial=jobrunner.run_in_serial, jobrunner=jobrunner
    )


# ═══════════════════════════════════════════════════════════════════════
# Private helpers
# ═══════════════════════════════════════════════════════════════════════


def _build_gaussian_pka_settings(
    proton_index, shared, opt_settings, sp_settings=None
):
    """Build a ``GaussianpKaJobSettings`` from shared options and project.

    All attributes of *opt_settings* that are accepted by
    ``GaussianJobSettings.__init__`` are forwarded automatically via
    ``**kwargs``, so new attributes on ``GaussianJobSettings`` are picked up
    without any changes here.
    """
    import inspect

    from chemsmart.jobs.gaussian.settings import (
        GaussianJobSettings,
        GaussianpKaJobSettings,
    )

    # ── 1. Translate CLI-side shared dict to GaussianpKaJobSettings names ──
    # Keys that are CLI-only (not constructor params) are dropped.
    # Two keys have different names between the CLI dict and the constructor.
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

    # ── 2. Forward all GaussianJobSettings attrs from opt_settings ──
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

    # ── 3. Prefer project solvent settings when CLI omitted ──
    # Solvent-related options: prefer CLI/shared, then opt settings, then sp settings
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
    additional_solvent_options = _first_non_none(
        pka_kwargs.get("additional_solvent_options"),
        getattr(opt_settings, "additional_solvent_options", None),
        (
            getattr(sp_settings, "additional_solvent_options", None)
            if sp_settings
            else None
        ),
    )
    custom_solvent = _first_non_none(
        pka_kwargs.get("custom_solvent"),
        getattr(opt_settings, "custom_solvent", None),
        getattr(sp_settings, "custom_solvent", None) if sp_settings else None,
    )

    pka_kwargs["solvent_model"] = solvent_model
    pka_kwargs["solvent_id"] = solvent_id
    if additional_solvent_options is not None:
        pka_kwargs["additional_solvent_options"] = additional_solvent_options
    if custom_solvent is not None:
        pka_kwargs["custom_solvent"] = custom_solvent

    return GaussianpKaJobSettings(
        proton_index=proton_index,
        **pka_kwargs,
        **opt_kwargs,
    )


def _create_pka_jobs_from_molecules(
    ctx, pka_molecules, shared, skip_completed, parallel, **kwargs
):
    """Create one ``GaussianpKaJob`` per ``PKaMolecule``."""
    from chemsmart.jobs.gaussian.pka import GaussianpKaBatchJob, GaussianpKaJob

    validate_reference_options(shared)

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    jobrunner = ctx.obj["jobrunner"]

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

        logger.info(
            f"Creating pKa job for fragment {idx}: "
            f"proton_index={pka_mol.proton_index}, label={mol_label}"
        )

        jobs.append(
            GaussianpKaJob(
                molecule=pka_mol,
                settings=pka_settings,
                label=mol_label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                parallel=parallel,
                **kwargs,
            )
        )

    logger.info(f"Created {len(jobs)} pKa jobs from multi-fragment CDXML file")
    return GaussianpKaBatchJob(
        jobs=jobs, run_in_serial=jobrunner.run_in_serial, jobrunner=jobrunner
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

    if shared["scheme"] == "proton exchange":
        if shared["reference"] is not None:
            logger.info(f"Reference acid file: {shared['reference']}")
            logger.info(
                f"Reference proton index: "
                f"{shared['reference_proton_index']}"
            )
            logger.info(
                f"Reference acid (HB): "
                f"charge={shared['reference_charge']}, "
                f"mult={shared['reference_multiplicity']}"
            )
        else:
            logger.info("No reference acid file provided")
    else:
        logger.info(f"delta_G(H+)_aq: {shared['delta_g_proton']} kcal/mol")

    logger.info(f"Gas phase: {pka_settings.functional}/{pka_settings.basis}")
    logger.info(
        f"Solution SP: {pka_settings.functional}/{pka_settings.basis} "
        f"with {shared['solvent_model']}({shared['solvent_id']})"
    )
    logger.info(
        f"Thermochem: T={shared['temperature']}K, "
        f"c={shared['concentration']}mol/L, "
        f"csg={shared['cutoff_entropy_grimme']}cm^-1, "
        f"ch={shared['cutoff_enthalpy']}cm^-1"
    )
