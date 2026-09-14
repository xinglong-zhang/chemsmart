"""
ORCA Minimum Energy Cross Point (MECP) CLI Module.

This module provides the command-line interface for ORCA MECP calculations
using the native ``SurfCrossOpt`` feature.  Unlike the Gaussian MECP driver
(which runs two SP+forces sub-jobs per Python-driven iteration), ORCA
SurfCrossOpt performs the entire crossing-seam optimisation internally in
a single ORCA invocation.

The two spin states must share the same charge and level of theory. Their
multiplicities are specified with the numbered state options.
"""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import (
    click_orca_solvent_options,
    orca,
)
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("mecp", cls=MyCommand)
@click_job_options
@click_orca_solvent_options
@click.option(
    "--multiplicity1",
    "-m1",
    type=click.IntRange(min=1),
    default=None,
    required=True,
    help="PES1 spin multiplicity; written to the * xyz line (required).",
)
@click.option(
    "--multiplicity2",
    "-m2",
    type=click.IntRange(min=1),
    default=None,
    required=True,
    help="PES2 spin multiplicity; written to %mecp Mult (required).",
)
@click.option(
    "--mode",
    type=click.Choice(["opt", "numfreq"], case_sensitive=False),
    default="opt",
    show_default=True,
    help=(
        "Optimisation mode. 'opt' (default) runs straight SurfCrossOpt; "
        "'numfreq' adds SurfCrossOpt NumFreq for a numerical-frequency "
        "verification step after convergence."
    ),
)
@click.option(
    "--maxiter",
    type=int,
    default=200,
    show_default=True,
    help="Maximum number of SurfCrossOpt iterations.",
)
@click.option(
    "--broken-sym",
    default=None,
    help=(
        "PES2 broken-symmetry pair NA,NB: unpaired electrons on two "
        "antiferromagnetically coupled centres. For example, 1,1 produces "
        "an open-shell singlet and therefore requires --multiplicity2 1."
    ),
)
@click.option(
    "--moinp",
    default=None,
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help="PES2 GBW guess file.",
)
@click.option(
    "-f",
    "--freeze-atoms",
    type=str,
    default=None,
    help="1-based atom indices to freeze, for example 1,3-5.",
)
@click.option(
    "-i",
    "--invert-constraints/--no-invert-constraints",
    default=False,
    help="Invert the frozen-atom selection.",
)
@click.option("--casscf-nel", type=int, default=None)
@click.option("--casscf-norb", type=int, default=None)
@click.option(
    "--casscf-mult",
    default=None,
    help="Comma-separated CASSCF multiplicities.",
)
@click.option(
    "--casscf-nroots", default=None, help="Comma-separated CASSCF root counts."
)
@click.option(
    "--casscf-bweight", default=None, help="Comma-separated CASSCF weights."
)
@click.pass_context
def mecp(
    ctx,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    solventfilename,
    multiplicity1,
    multiplicity2,
    mode,
    maxiter,
    broken_sym,
    moinp,
    freeze_atoms,
    invert_constraints,
    casscf_nel,
    casscf_norb,
    casscf_mult,
    casscf_nroots,
    casscf_bweight,
    skip_completed,
    **kwargs,
):
    """
    Run ORCA Minimum Energy Cross Point (MECP) calculations.

    Performs a geometry optimisation on the crossing seam between two spin
    states of the same charge and using the same level of theory, using
    ORCA's native SurfCrossOpt feature.

    The multiplicities of the two states are required.  Charge and
    computational method are inherited from the project configuration and
    can be overridden via the standard ``-c`` / ``-m`` / ``--method`` /
    ``--basis`` flags of the parent ``orca`` command.
    """
    from chemsmart.jobs.orca.settings import ORCAMECPJobSettings

    # get jobrunner from context
    jobrunner = ctx.obj["jobrunner"]

    # get settings from project — MECP inherits from the project's opt
    # settings (method, basis, solvent, SCF, etc.)
    project_settings = ctx.obj["project_settings"]
    mecp_project_settings = project_settings.opt_settings()

    # job setting from filename or default, with updates from user in cli
    # specified in keywords (e.g., `chemsmart orca -c <charge> -m <mult>`)
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job-level settings from cli keywords
    mecp_project_settings = mecp_project_settings.merge(
        job_settings, keywords=keywords
    )

    # cli-supplied solvent model, solvent id, and additional solvent options
    mecp_project_settings.modify_solvent(
        remove_solvent=remove_solvent,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
    )
    if solvent_options is not None:
        mecp_project_settings.additional_solvent_options = solvent_options
    if solventfilename is not None:
        mecp_project_settings.solventfilename = solventfilename

    # convert to MECP-specific settings — this is where we add the
    # SurfCrossOpt-specific fields on top of the inherited project settings
    mecp_settings = ORCAMECPJobSettings.from_settings(mecp_project_settings)

    # The two surfaces have independent multiplicities but ORCA SurfCrossOpt
    # requires them to share one charge and one level of theory.
    mecp_settings.multiplicity1 = multiplicity1
    mecp_settings.multiplicity2 = multiplicity2
    mecp_settings.multiplicity = multiplicity1

    # SurfCrossOpt optimisation mode (opt or numfreq)
    mecp_settings.mode = mode.lower()

    mecp_settings.maxiter = maxiter
    logger.debug(f"Set SurfCrossOpt MaxIter: {maxiter}")

    def _comma_values(value):
        if value is None:
            return None
        return [int(item) for item in value.split(",")]

    mecp_settings.broken_sym = _comma_values(broken_sym)
    mecp_settings.moinp = moinp
    mecp_settings.invert_constraints = invert_constraints
    mecp_settings.casscf_nel = casscf_nel
    mecp_settings.casscf_norb = casscf_norb
    mecp_settings.casscf_mult = _comma_values(casscf_mult)
    mecp_settings.casscf_nroots = _comma_values(casscf_nroots)
    mecp_settings.casscf_bweight = _comma_values(casscf_bweight)

    # validate charge and multiplicity consistency (state-A multiplicity
    # is inherited from project/job settings)
    check_charge_and_multiplicity(mecp_settings)
    mecp_settings.validate()

    # get molecule from context
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1].copy()  # get last molecule from list
    if freeze_atoms is not None:
        from chemsmart.utils.utils import (
            convert_list_to_gaussian_frozen_list,
            get_list_from_string_range,
        )

        frozen_atoms = get_list_from_string_range(freeze_atoms)
        molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
            frozen_atoms, molecule
        )
    logger.info(f"Running ORCA MECP on molecule: {molecule}")

    # get label for the job output files
    label = ctx.obj["label"]

    logger.info(f"Final MECP job settings: {mecp_settings.__dict__}")

    from chemsmart.jobs.orca.mecp import ORCAMECPJob

    job = ORCAMECPJob(
        molecule=molecule,
        settings=mecp_settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
    logger.debug(f"Created ORCA MECP job: {job}")
    return job
