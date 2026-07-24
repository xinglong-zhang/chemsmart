import logging

logger = logging.getLogger(__name__)


def build_xtb_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    """Build one or more xTB jobs shared across the opt/sp/hess subcommands."""
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            # Preserve one output directory per selected structure.
            molecule_label = f"{label}_idx{idx}"
            logger.info(f"Creating xTB job {molecule_label}")
            jobs.append(
                job_cls(
                    molecule=molecule,
                    settings=settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        return jobs

    molecule = molecules[-1]
    logger.info(f"Creating xTB job {label}")
    return job_cls(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
