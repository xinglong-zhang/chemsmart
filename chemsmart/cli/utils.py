import logging

logger = logging.getLogger(__name__)

CHAIN_PROJECT_SETTINGS_KEY = "chain_project_settings"
CHAIN_CLI_DEFAULTS_KEY = "chain_cli_defaults"


def build_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    """Build one or more jobs from molecules stored on the Click context."""
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            # Preserve one output directory per selected structure.
            molecule_label = f"{label}_idx{idx}"
            logger.info(f"Creating job {molecule_label}")
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
    logger.info(f"Creating job {label}")
    return job_cls(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
