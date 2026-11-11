import logging
import os
import tempfile

import click

logger = logging.getLogger(__name__)


def build_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    """Build one or more jobs from molecules stored on the Click context."""
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]
    batch_label = ctx.obj.get("batch_label")
    if batch_label is not None:
        if len(molecules) != 1:
            raise ValueError(
                "A batch subtask must resolve to exactly one molecule."
            )
        # Establish the final label before constructors derive their job folder.
        label = batch_label

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


def prepare_structure_inputs(jobs, source, commands, prepare_single=False):
    """Yield ``(job, input_path)`` pairs for geometry-only batch inputs.

    The files must remain accessible until the submitted jobs have finished.
    Paths are None when the input is not suitable for preparation.
    """
    eligible = (
        (len(jobs) > 1 or prepare_single)
        and source
        and any(c.get("name") in {"xtb", "gaussian", "orca"} for c in commands)
        and all(not getattr(job, "_batch_entry", None) for job in jobs)
    )
    if eligible and not str(source).lower().endswith((".npz", ".xyz")):
        raise click.UsageError(
            "Multi-structure submission currently supports NPZ and XYZ "
            "inputs only; convert the input to XYZ before submitting."
        )
    prepare = eligible
    input_dir = (
        tempfile.mkdtemp(prefix="chemsmart_inputs_", dir=os.getcwd())
        if prepare
        else None
    )
    for index, job in enumerate(jobs, start=1):
        path = None
        if input_dir:
            path = os.path.join(input_dir, f"structure_{index:05d}.xyz")
            job.molecule.write_xyz(path, mode="w")
        yield job, path


def prepared_structure_arguments(cli_args, input_path, label):
    """Replace the source with one prepared structure and preserve its label."""
    args = list(cli_args)
    for option in ("--filename", "-f"):
        if option in args:
            args[args.index(option) + 1] = input_path
            break
    else:
        raise ValueError("Cannot locate the structure input filename.")
    for option in ("--index", "-i", "--structure-index", "--si"):
        if option in args:
            position = args.index(option)
            del args[position : position + 2]
    return ["--batch-label", str(label), *args]
