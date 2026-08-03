"""Shared job construction for the PySCF CLI leaves."""

import logging

from chemsmart.utils.cli import create_sp_label
logger = logging.getLogger(__name__)


def resolve_engine(ctx, settings):
    """Resolve the execution engine and record it on ``settings``.

    GPU execution is never inferred from the server inventory. An explicit
    ``--gpu`` plus a positive resolved allocation is required; otherwise the
    fail-closed default is CPU. The resolved value is written back onto the
    settings so it reaches the generated script and results provenance -- GPU
    and CPU results are scientifically distinct execution requests.
    """
    jobrunner = ctx.obj.get("jobrunner")
    num_gpus = getattr(jobrunner, "num_gpus", None) if jobrunner else None
    if isinstance(num_gpus, bool) or not isinstance(num_gpus, int):
        raise ValueError(
            "Resolved server NUM_GPUS must be a non-negative integer, got "
            f"{num_gpus!r}."
        )
    if num_gpus < 0:
        raise ValueError(
            "Resolved server NUM_GPUS must be a non-negative integer, got "
            f"{num_gpus!r}."
        )

    if "engine" in ctx.obj.get("keywords", ()):
        if settings.engine == "gpu" and num_gpus <= 0:
            raise ValueError(
                "Explicit --gpu requires a positive resolved NUM_GPUS; "
                "ChemSmart will not invent device 0."
            )
        logger.debug(f"Engine set explicitly: {settings.engine}")
        return settings.engine

    settings.engine = "cpu"
    logger.debug(
        "No explicit --gpu request; resolved fail-closed CPU engine with "
        f"server NUM_GPUS={num_gpus}"
    )
    return settings.engine


def build_pyscf_jobs(ctx, job_class, settings, skip_completed, kwargs):
    """Build one PySCF job per selected molecule.

    Validation runs here rather than at execution time so that an
    unsupported request fails before any SCF cycle is spent.
    """
    resolve_engine(ctx, settings)

    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]
    explicit_state_fields = ctx.obj.get(
        "explicit_state_fields", frozenset()
    )

    def settings_for_molecule(molecule):
        """Resolve source state unless project or CLI explicitly replaced it."""

        resolved = settings.copy()
        if (
            "charge" not in explicit_state_fields
            and molecule.charge is not None
        ):
            resolved.charge = molecule.charge
        if (
            "multiplicity" not in explicit_state_fields
            and molecule.multiplicity is not None
        ):
            resolved.multiplicity = molecule.multiplicity
        if resolved.charge is None or resolved.multiplicity is None:
            raise ValueError(
                "PySCF charge and multiplicity must come from the molecular "
                "source or be supplied explicitly in project YAML/CLI."
            )
        resolved.validate()
        return resolved

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_settings = settings_for_molecule(molecule)
            molecule_label = create_sp_label(
                f"{label}_idx{idx}", molecule_settings
            )
            jobs.append(
                job_class(
                    molecule=molecule,
                    settings=molecule_settings,
                    label=molecule_label,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        logger.debug(f"Created {len(jobs)} PySCF jobs")
        return jobs

    molecule = molecules[-1]
    molecule_settings = settings_for_molecule(molecule)
    final_label = create_sp_label(label, molecule_settings)
    logger.info(
        f"Running PySCF {molecule_settings.jobtype} on {molecule} "
        f"with label {final_label}"
    )
    return job_class(
        molecule=molecule,
        settings=molecule_settings,
        label=final_label,
        skip_completed=skip_completed,
        **kwargs,
    )
