"""Shared job construction for the PySCF CLI leaves."""

import logging

from chemsmart.utils.cli import create_sp_label

logger = logging.getLogger(__name__)


def resolve_engine(ctx, settings):
    """Resolve the execution engine and record it on ``settings``.

    The merged project/CLI settings are authoritative: omission of the CLI
    flag preserves the project engine, while an explicit ``--gpu`` or
    ``--no-gpu`` has already replaced it through Click's merge whitelist.
    GPU execution is never inferred from server inventory and every GPU
    request used for execution requires a positive resolved allocation. A
    ``--fake`` invocation may still render the exact GPU script on a CPU host;
    it remains non-executable and never changes the request to CPU.
    """
    explicit_engine = "engine" in ctx.obj.get("keywords", ())
    requested_engine = str(settings.engine or "").strip().lower()
    if requested_engine not in {"cpu", "gpu"}:
        raise ValueError(
            f"Unsupported PySCF engine {settings.engine!r}; expected cpu or "
            "gpu."
        )

    if requested_engine == "cpu":
        # A CPU request does not depend on scheduler GPU inventory.  This is
        # particularly important for ``sub --test --no-gpu`` previews, where
        # no server allocation has been resolved yet.  Requiring NUM_GPUS in
        # that path made a valid CPU command fail for an unrelated reason.
        logger.debug(
            "Resolved CPU engine from %s; GPU inventory is unused",
            "explicit CLI override" if explicit_engine else "project YAML",
        )
        return settings.engine

    jobrunner = ctx.obj.get("jobrunner")
    num_gpus = getattr(jobrunner, "num_gpus", None) if jobrunner else None
    fake_preview = getattr(jobrunner, "fake", False) is True

    if (
        isinstance(num_gpus, bool)
        or not isinstance(num_gpus, int)
        or num_gpus <= 0
    ):
        if fake_preview:
            logger.debug(
                "Preserving GPU engine for fake preview without a device "
                "allocation; execution remains unavailable"
            )
            return settings.engine
        origin = "explicit --gpu" if explicit_engine else "project YAML"
        raise ValueError(
            f"PySCF GPU engine requested by {origin}, but execution requires "
            "a positive resolved NUM_GPUS; ChemSmart will not fall back to "
            "CPU or invent device 0."
        )
    logger.debug(
        "Resolved GPU engine from %s with NUM_GPUS=%s",
        "explicit CLI override" if explicit_engine else "project YAML",
        num_gpus,
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
    explicit_state_fields = ctx.obj.get("explicit_state_fields", frozenset())

    def settings_for_molecule(molecule):
        """Resolve source state unless project or CLI explicitly replaced it."""

        resolved = settings.copy()
        if (resolved.charge is None) != (resolved.multiplicity is None):
            raise ValueError(
                "PySCF charge and multiplicity overrides must be supplied "
                "together, or both inherited from the molecular source."
            )
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
