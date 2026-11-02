"""Shared job construction for the PySCF CLI leaves."""

import logging

from chemsmart.utils.cli import create_sp_label
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


def resolve_engine(ctx, settings):
    """Resolve the execution engine and record it on ``settings``.

    Mirrors the scratch ladder: an explicit ``--gpu/--no-gpu`` wins, then
    ``engine`` in the project YAML, then the server's GPU count, then CPU.
    The resolved value is written back onto the
    settings so it reaches the generated script and the results provenance --
    GPU and CPU results are not bit-identical, so every number has to state
    which engine produced it.
    """
    if settings.engine is not None:
        logger.debug(f"Engine set explicitly: {settings.engine}")
        return settings.engine

    jobrunner = ctx.obj.get("jobrunner")
    num_gpus = jobrunner.num_gpus if jobrunner else None
    settings.engine = "gpu" if num_gpus else "cpu"
    logger.debug(
        f"Engine resolved from server NUM_GPUS={num_gpus}: {settings.engine}"
    )
    return settings.engine


def settings_for_molecule(settings, molecule):
    """Return settings with charge and multiplicity resolved for a molecule.

    An explicit ``-c``/``-m`` or project value wins; otherwise the input file's
    own values apply, as when reading a calculation output. Formats that carry
    neither, such as plain XYZ, still have to be given them explicitly.
    """
    settings = settings.copy()
    if settings.charge is None:
        settings.charge = molecule.charge
    if settings.multiplicity is None:
        settings.multiplicity = molecule.multiplicity
    check_charge_and_multiplicity(settings)
    return settings


def build_pyscf_jobs(ctx, job_class, settings, skip_completed, kwargs):
    """Build one PySCF job per selected molecule.

    Validation runs here rather than at execution time so that an
    unsupported request fails before any SCF cycle is spent.
    """
    resolve_engine(ctx, settings)
    settings.validate()

    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = create_sp_label(f"{label}_idx{idx}", settings)
            jobs.append(
                job_class(
                    molecule=molecule,
                    settings=settings_for_molecule(settings, molecule),
                    label=molecule_label,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        logger.debug(f"Created {len(jobs)} PySCF jobs")
        return jobs

    molecule = molecules[-1]
    final_label = create_sp_label(label, settings)
    logger.info(
        f"Running PySCF {settings.jobtype} on {molecule} "
        f"with label {final_label}"
    )
    return job_class(
        molecule=molecule,
        settings=settings_for_molecule(settings, molecule),
        label=final_label,
        skip_completed=skip_completed,
        **kwargs,
    )
