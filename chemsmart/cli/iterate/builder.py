"""Input-neutral Iterate job construction."""

import logging
import sys

import click

from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.settings import (
    IterateJobSettings,
    resolve_algorithm_config,
)

logger = logging.getLogger(__name__)


def build_iterate_job(
    ctx, cli_algorithm_name=None, cli_options=None
) -> IterateJob:
    """Build an Iterate job from an already normalized input config."""
    data = ctx.obj["iterate"]
    config = data.get("config")
    if config is None:
        raise click.ClickException(
            "Iterate input adapter did not provide a normalized config."
        )

    nprocs = ctx.obj["jobrunner"].num_cores
    if nprocs < 1:
        raise click.BadParameter(
            "must be at least 1",
            param_hint="'-n' / '--num-cores'",
        )

    try:
        algorithm_config = resolve_algorithm_config(
            yaml_algorithm=config.get("algorithm"),
            cli_algorithm_name=cli_algorithm_name,
            cli_options=cli_options,
        )
    except ValueError as exc:
        raise click.BadParameter(str(exc))

    logger.debug(f"  Skeletons: {len(config['skeletons'])}")
    logger.debug(f"  Substituents: {len(config['substituents'])}")
    logger.debug(
        f"Using algorithm '{algorithm_config.name}' "
        f"with options {algorithm_config.options}"
    )
    logger.debug(f"  Combination mode: {data['combination_mode']}")

    job_settings = IterateJobSettings(
        config_file=data.get("config_file"),
        algorithm_config=algorithm_config,
        combination_mode=data["combination_mode"],
        max_substituted_sites=data["max_substituted_sites"],
    )
    job_settings.skeleton_list = config["skeletons"]
    job_settings.substituent_list = config["substituents"]

    job = IterateJob(
        settings=job_settings,
        jobrunner=ctx.obj.get("jobrunner"),
        nprocs=nprocs,
        timeout=data["timeout"],
        outputfile=data["outputfile"],
        separate_outputs=data["separate_outputs"],
        output_directory=data["directory"],
        command_line=" ".join(sys.argv),
        show_worker_logs=logger.isEnabledFor(logging.DEBUG),
        skip_completed=data["skip_completed"],
    )

    logger.debug(f"Created IterateJob with {nprocs} process(es)")
    return job
