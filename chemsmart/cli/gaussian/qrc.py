import logging

import click

from chemsmart.cli.gaussian.gaussian import (
    click_gaussian_jobtype_options,
    gaussian,
)
from chemsmart.cli.job import (
    click_job_options,
    click_molecule_vibrational_displacement_options,
)
from chemsmart.utils.cli import (
    MyCommand,
    get_setting_from_jobtype_for_gaussian,
)
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@gaussian.command(cls=MyCommand)
@click_job_options
@click_gaussian_jobtype_options
@click_molecule_vibrational_displacement_options
@click.pass_context
def qrc(
    ctx,
    jobtype,
    coordinates,
    step_size,
    num_steps,
    mode_idx,
    amp,
    nframes,
    phase,
    normalize,
    return_xyz,
    skip_completed,
    **kwargs,
):
    """
    CLI subcommand for running Gaussian QRC jobs.

    Examples:
        `chemsmart sub gaussian -p proj -f ts.log qrc` runs QRC optimization jobs
        after displacing along vibrational mode 1 (default) by +/- 0.5 Angstroms (default).

        `chemsmart sub gaussian -p proj -f ts.log qrc -m 2 -a 1.2` runs QRC optimization jobs
        after displacing along vibrational mode 2 (-m 2) by +/- 1.2 Angstroms (-a 1.2).

        `chemsmart sub gaussian -p proj -f ts.log qrc -j ts -m 2 -a 1.5` runs QRC TS jobs
        after displacing along vibrational mode 2 (-m 2) by +/- 1.5 Angstroms (-a 1.5).

        `chemsmart sub gaussian -p proj -f ts.log qrc -j modred -c [1,2] -m 2 -a 1.5`
        runs QRC modred jobs after displacing along vibrational mode 2 by +/- 1.5 Angstroms.
    """

    # get jobrunner for running Gaussian crest jobs
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project
    project_settings = ctx.obj["project_settings"]

    # get job settings from job type
    if jobtype is None:
        jobtype = "opt"  # default using opt job for QRC, if jobtype not given

    qrc_settings = get_setting_from_jobtype_for_gaussian(
        project_settings, jobtype, coordinates, step_size, num_steps
    )

    # job setting from cli
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    qrc_settings = qrc_settings.merge(job_settings, keywords=keywords)

    if ctx.invoked_subcommand is not None:
        return

    check_charge_and_multiplicity(qrc_settings)

    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]  # use last structure

    # get label
    label = ctx.obj["label"]

    logger.info(
        f"QRC {jobtype} settings from project: {qrc_settings.__dict__}"
    )

    from chemsmart.jobs.gaussian.qrc import GaussianQRCJob

    logger.debug(f"Creating GaussianQRCJob using molecule {molecule}")

    return GaussianQRCJob(
        molecule=molecule,
        settings=qrc_settings,
        label=label,
        jobrunner=jobrunner,
        mode_idx=mode_idx,
        amp=amp,
        nframes=nframes,
        phase=phase,
        normalize=normalize,
        return_xyz=return_xyz,
        skip_completed=skip_completed,
        **kwargs,
    )
