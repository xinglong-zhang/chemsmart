"""
CLI for ORCA pKa input generation (job submission).

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

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.pka import (
    build_per_entry_subcommands,
    click_pka_proton_options,
    click_pka_shared_options,
    click_pka_submit_options,
    resolve_proton_index,
    validate_reference_options,
)
from chemsmart.utils.cli import MyCommand, MyGroup

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# pka group
# ═══════════════════════════════════════════════════════════════════════


@orca.group("pka", cls=MyGroup, invoke_without_command=True)
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
    """ORCA pKa job submission.

    \b
    Subcommands:
      submit         Single-molecule job submission (default).
      batch          Table-driven batch submission.

    When invoked without a subcommand the ``submit`` path runs
    automatically.

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
    """Submit a single-molecule ORCA pKa calculation.

    \b
    Examples:
      chemsmart run orca -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -r ref.xyz -rpi 1 -rc 0 -rm 1 submit

      chemsmart run orca -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -s direct submit
    """
    shared = ctx.obj["pka_shared"]
    filename = ctx.obj.get("filename")
    proton_index = ctx.obj.get("pka_proton_index")
    color_code = ctx.obj.get("pka_color_code")
    jobrunner = ctx.obj["jobrunner"]

    # Align pKa parallel execution with global --run-in-serial flag
    parallel = not jobrunner.run_in_serial

    proton_index, pka_molecules = resolve_proton_index(
        filename, proton_index, color_code
    )

    if pka_molecules is not None:
        return _create_orca_pka_jobs_from_molecules(
            ctx, pka_molecules, shared, skip_completed, parallel, **kwargs
        )

    validate_reference_options(shared)

    from chemsmart.jobs.orca.pka import ORCApKaJob

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = opt_settings.merge(job_settings, keywords=keywords)

    pka_settings = _build_orca_pka_settings(proton_index, shared, opt_settings)

    if pka_settings.charge is None:
        raise click.UsageError(
            "Charge must be specified via -c/--charge option"
        )
    if pka_settings.multiplicity is None:
        raise click.UsageError(
            "Multiplicity must be specified via -m/--multiplicity option"
        )

    logger.info(f"ORCA pKa job settings: {pka_settings.__dict__}")
    logger.info(f"Proton index to remove: {proton_index}")
    logger.info(f"Thermodynamic cycle: {shared['scheme']}")

    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj.get("molecule_indices")
    label = ctx.obj["label"]
    base_label = label if label.endswith("_pka") else f"{label}_pka"

    if len(molecules) > 1 and molecule_indices:
        logger.info(f"Creating {len(molecules)} ORCA pKa jobs")
        return [
            ORCApKaJob(
                molecule=mol,
                settings=pka_settings,
                label=f"{base_label}_idx{idx}",
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                parallel=parallel,
                **kwargs,
            )
            for mol, idx in zip(molecules, molecule_indices)
        ]

    return ORCApKaJob(
        molecule=molecules[-1],
        settings=pka_settings,
        label=base_label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        parallel=parallel,
        **kwargs,
    )


# ═══════════════════════════════════════════════════════════════════════
# pka batch
# ═══════════════════════════════════════════════════════════════════════


@pka.command("batch", cls=MyCommand)
@click_job_options
@click.pass_context
def batch(ctx, skip_completed, **kwargs):
    """Table-driven batch ORCA pKa job submission.

    \b
    Table format (4 columns, whitespace or comma-delimited):
        filepath    proton_index    charge    multiplicity

    \b
    Examples:
      chemsmart run orca -p myproject -f molecules.txt pka \\
          -s "proton exchange" -r ref.xyz -rpi 5 -rc 0 -rm 1 batch
    """
    shared = ctx.obj["pka_shared"]
    jobrunner = ctx.obj["jobrunner"]

    # Align pKa parallel execution with global --run-in-serial flag
    parallel = not jobrunner.run_in_serial

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent ORCA -f/--filename to "
            "specify the table file path."
        )

    from chemsmart.io.molecules.structure import Molecule
    from chemsmart.jobs.orca.pka import ORCApKaJob
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings
    from chemsmart.utils.utils import (
        parse_pka_table,
        validate_pka_table_entries,
    )

    logger.info(f"Reading ORCA pKa jobs from table: {input_table_path}")
    try:
        entries = parse_pka_table(input_table_path)
        validate_pka_table_entries(entries, check_file_exists=True)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(entries)} entries in table")

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
                            ref, None, shared["reference_color_code"]
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

    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj.get("job_settings")
    kw = ctx.obj.get("keywords", {})
    if job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=kw)

    jobs = []
    for entry in entries:
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        molecule = Molecule.from_filepath(filepath)
        label = Path(filepath).stem
        base_label = label if label.endswith("_pka") else f"{label}_pka"

        solvent_model = shared["solvent_model"]
        if solvent_model is None:
            try:
                solvent_model = opt_settings.solvent_model
            except AttributeError:
                solvent_model = None
        solvent_id = shared["solvent_id"]
        if solvent_id is None:
            try:
                solvent_id = opt_settings.solvent_id
            except AttributeError:
                solvent_id = None
        if solvent_model is None:
            solvent_model = "CPCM"
        if solvent_id is None:
            solvent_id = "water"

        pka_settings = ORCApKaJobSettings(
            proton_index=int(entry.proton_index),
            scheme=shared["scheme"],
            reference_file=shared["reference"],
            reference_proton_index=shared["reference_proton_index"],
            reference_charge=shared["reference_charge"],
            reference_multiplicity=shared["reference_multiplicity"],
            reference_conjugate_base_charge=shared[
                "reference_conjugate_base_charge"
            ],
            reference_conjugate_base_multiplicity=shared[
                "reference_conjugate_base_multiplicity"
            ],
            delta_G_proton=shared["delta_g_proton"],
            conjugate_base_charge=None,
            conjugate_base_multiplicity=None,
            solvent_model=solvent_model,
            solvent_id=solvent_id,
            temperature=shared["temperature"],
            concentration=shared["concentration"],
            cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
            cutoff_enthalpy=shared["cutoff_enthalpy"],
            charge=int(entry.charge),
            multiplicity=int(entry.multiplicity),
            functional=opt_settings.functional,
            basis=opt_settings.basis,
            ab_initio=opt_settings.ab_initio,
            dispersion=opt_settings.dispersion,
            aux_basis=opt_settings.aux_basis,
            defgrid=opt_settings.defgrid,
            semiempirical=opt_settings.semiempirical,
            additional_route_parameters=(
                opt_settings.additional_route_parameters
            ),
            gen_genecp_file=opt_settings.gen_genecp_file,
            heavy_elements=opt_settings.heavy_elements,
            heavy_elements_basis=opt_settings.heavy_elements_basis,
            light_elements_basis=opt_settings.light_elements_basis,
        )

        jobs.append(
            ORCApKaJob(
                molecule=molecule,
                settings=pka_settings,
                label=base_label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                parallel=parallel,
                **kwargs,
            )
        )

        # Build per-entry subcommand override so that `chemsmart sub`
        # generates a run script that processes ONLY this entry instead
        # of re-running the entire batch CSV.
        jobs[-1]._batch_subcommands_override = build_per_entry_subcommands(
            ctx,
            filepath=filepath,
            charge=int(entry.charge),
            multiplicity=int(entry.multiplicity),
            proton_index=int(entry.proton_index),
        )

    logger.info(f"Created {len(jobs)} ORCA pKa jobs from table")
    return jobs


# ═══════════════════════════════════════════════════════════════════════
# Private helpers
# ═══════════════════════════════════════════════════════════════════════


def _build_orca_pka_settings(proton_index, shared, opt_settings):
    """Build an ``ORCApKaJobSettings`` from shared options and project."""
    from chemsmart.jobs.orca.settings import ORCApKaJobSettings

    solvent_model = shared["solvent_model"]
    if solvent_model is None:
        try:
            solvent_model = opt_settings.solvent_model
        except AttributeError:
            solvent_model = None
    solvent_id = shared["solvent_id"]
    if solvent_id is None:
        try:
            solvent_id = opt_settings.solvent_id
        except AttributeError:
            solvent_id = None
    if solvent_model is None:
        solvent_model = "CPCM"
    if solvent_id is None:
        solvent_id = "water"

    return ORCApKaJobSettings(
        proton_index=proton_index,
        scheme=shared["scheme"],
        reference_file=shared["reference"],
        reference_proton_index=shared["reference_proton_index"],
        reference_charge=shared["reference_charge"],
        reference_multiplicity=shared["reference_multiplicity"],
        reference_conjugate_base_charge=shared[
            "reference_conjugate_base_charge"
        ],
        reference_conjugate_base_multiplicity=shared[
            "reference_conjugate_base_multiplicity"
        ],
        delta_G_proton=shared["delta_g_proton"],
        conjugate_base_charge=shared["conjugate_base_charge"],
        conjugate_base_multiplicity=shared["conjugate_base_multiplicity"],
        solvent_model=solvent_model,
        solvent_id=solvent_id,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        charge=opt_settings.charge,
        multiplicity=opt_settings.multiplicity,
        functional=opt_settings.functional,
        basis=opt_settings.basis,
        ab_initio=opt_settings.ab_initio,
        dispersion=opt_settings.dispersion,
        aux_basis=opt_settings.aux_basis,
        defgrid=opt_settings.defgrid,
        semiempirical=opt_settings.semiempirical,
        additional_route_parameters=opt_settings.additional_route_parameters,
        gen_genecp_file=opt_settings.gen_genecp_file,
        heavy_elements=opt_settings.heavy_elements,
        heavy_elements_basis=opt_settings.heavy_elements_basis,
        light_elements_basis=opt_settings.light_elements_basis,
    )


def _create_orca_pka_jobs_from_molecules(
    ctx, pka_molecules, shared, skip_completed, parallel, **kwargs
):
    """Create one ``ORCApKaJob`` per ``PKaMolecule``."""
    from chemsmart.jobs.orca.pka import ORCApKaJob

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
        pka_settings = _build_orca_pka_settings(
            pka_mol.proton_index, shared, opt_settings
        )

        logger.info(
            f"Creating ORCA pKa job for fragment {idx}: "
            f"proton_index={pka_mol.proton_index}, label={mol_label}"
        )

        jobs.append(
            ORCApKaJob(
                molecule=pka_mol,
                settings=pka_settings,
                label=mol_label,
                jobrunner=jobrunner,
                skip_completed=skip_completed,
                parallel=parallel,
                **kwargs,
            )
        )

    logger.info(f"Created {len(jobs)} ORCA pKa jobs from multi-fragment CDXML")
    return jobs
