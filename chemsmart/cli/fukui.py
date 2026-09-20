"""Fukui submit options, job builder, and backend-independent analysis.

Job submission lives under ``chemsmart sub … gaussian … fukui`` or
``chemsmart sub … orca … fukui`` via ``register_fukui_cli``. Analysis is
``chemsmart run fukui`` (also ``chemsmart run chain fukui analyze``).
"""

import functools
import logging

import click

from chemsmart.analysis.fukui import (
    FUKUI_MODES,
    analyze_fukui,
    discover_fukui_companion_outputs,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


def click_fukui_submit_options(f=None, *, modes=None):
    """Fukui submit options shared by Gaussian and ORCA program commands."""
    if modes is None:
        modes = FUKUI_MODES

    def decorator(func):
        @click.option(
            "--mode",
            default="mulliken",
            show_default=True,
            type=click.Choice(list(modes), case_sensitive=False),
            help="Charges to be used for Fukui Indices calculations.",
        )
        @click.option(
            "-rcc",
            "--radical-cation-charge",
            type=int,
            default=None,
            help=(
                "Override charge for the radical-cation job. "
                "Default is derived from the neutral charge."
            ),
        )
        @click.option(
            "-rcm",
            "--radical-cation-multiplicity",
            type=int,
            default=None,
            help=(
                "Override multiplicity for the radical-cation job. "
                "Default is derived from the neutral multiplicity."
            ),
        )
        @click.option(
            "-rac",
            "--radical-anion-charge",
            type=int,
            default=None,
            help=(
                "Override charge for the radical-anion job. "
                "Default is derived from the neutral charge."
            ),
        )
        @click.option(
            "-ram",
            "--radical-anion-multiplicity",
            type=int,
            default=None,
            help=(
                "Override multiplicity for the radical-anion job. "
                "Default is derived from the neutral multiplicity."
            ),
        )
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

    if f is None:
        return decorator
    return decorator(f)


def build_fukui_job(
    ctx,
    job_cls,
    skip_completed,
    mode,
    radical_cation_charge,
    radical_cation_multiplicity,
    radical_anion_charge,
    radical_anion_multiplicity,
    nbo_uses_wbi=False,
    **kwargs,
):
    """Create one Fukui job from parent ``-f`` and shared Fukui options."""
    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    if nbo_uses_wbi and mode.lower() == "nbo":
        pop_settings = project_settings.wbi_settings()
    else:
        pop_settings = project_settings.sp_settings()
    pop_settings = pop_settings.merge(job_settings, keywords=keywords)
    check_charge_and_multiplicity(pop_settings)

    molecule = ctx.obj["molecules"][-1]
    label = ctx.obj["label"]
    logger.info("Creating %s job: mode=%s", job_cls.__name__, mode)
    return job_cls(
        molecule=molecule,
        settings=pop_settings,
        label=label,
        jobrunner=jobrunner,
        mode=mode,
        radical_cation_charge=radical_cation_charge,
        radical_cation_multiplicity=radical_cation_multiplicity,
        radical_anion_charge=radical_anion_charge,
        radical_anion_multiplicity=radical_anion_multiplicity,
        skip_completed=skip_completed,
        **kwargs,
    )


def register_fukui_cli(
    parent_group, job_cls, *, modes=FUKUI_MODES, nbo_uses_wbi=False
):
    """Attach ``fukui`` submit to a Gaussian or ORCA Click group."""

    @parent_group.command("fukui", cls=MyCommand)
    @click_job_options
    @click_fukui_submit_options(modes=modes)
    @click.pass_context
    def fukui(
        ctx,
        skip_completed,
        mode,
        radical_cation_charge,
        radical_cation_multiplicity,
        radical_anion_charge,
        radical_anion_multiplicity,
        **kwargs,
    ):
        """Submit Fukui charge-state calculations.

        Always runs neutral, radical-cation, and radical-anion population
        jobs from the parent ``-f`` structure. Analyze completed outputs
        with ``chemsmart run fukui``.
        """
        return build_fukui_job(
            ctx,
            job_cls,
            skip_completed,
            mode,
            radical_cation_charge,
            radical_cation_multiplicity,
            radical_anion_charge,
            radical_anion_multiplicity,
            nbo_uses_wbi=nbo_uses_wbi,
            **kwargs,
        )

    return fukui


@click.command(name="fukui", cls=MyCommand)
@click.option(
    "-n",
    "--neutral-filename",
    required=True,
    type=str,
    help="Gaussian or ORCA output file for the neutral system.",
)
@click.option(
    "-c",
    "--radical-cation-filename",
    default=None,
    type=str,
    help="Gaussian or ORCA output file for the radical cationic system.",
)
@click.option(
    "-a",
    "--radical-anion-filename",
    default=None,
    type=str,
    help="Gaussian or ORCA output file for the radical anionic system.",
)
@click.option(
    "-m",
    "--mode",
    default="mulliken",
    show_default=True,
    type=click.Choice(list(FUKUI_MODES), case_sensitive=False),
    help="Charges to be used for Fukui Indices calculations.",
)
@click.option(
    "-o",
    "--output",
    default=None,
    type=str,
    help="Write Fukui results to this file instead of logging them.",
)
def fukui(
    neutral_filename,
    radical_cation_filename=None,
    radical_anion_filename=None,
    mode="mulliken",
    output=None,
):
    """Compute Fukui reactivity indices from existing output files.

    Companion ``_rc`` / ``_ra`` files are auto-discovered from ``-n`` when
    omitted (same labels as Gaussian Fukui job submission).

    \b
    Examples:
      chemsmart run fukui -n mol_n.log -c mol_rc.log -a mol_ra.log

      chemsmart run fukui -n mol_n.log -m nbo

      chemsmart run fukui -n mol_n.log -o fukui.dat
    """
    if radical_cation_filename is None or radical_anion_filename is None:
        discovered = discover_fukui_companion_outputs(neutral_filename)
        if radical_cation_filename is None and discovered["radical_cation"]:
            radical_cation_filename = discovered["radical_cation"]
            if output is None:
                logger.info(
                    f"Auto-discovered radical cation: {radical_cation_filename}"
                )
        if radical_anion_filename is None and discovered["radical_anion"]:
            radical_anion_filename = discovered["radical_anion"]
            if output is None:
                logger.info(
                    f"Auto-discovered radical anion: {radical_anion_filename}"
                )

    if radical_cation_filename is None and radical_anion_filename is None:
        raise click.UsageError(
            "At least one of -c/--radical-cation-filename or "
            "-a/--radical-anion-filename must be provided (or auto-discoverable "
            "as <base>_rc / <base>_ra beside -n)."
        )

    analyze_fukui(
        neutral_filename=neutral_filename,
        radical_cation_filename=radical_cation_filename,
        radical_anion_filename=radical_anion_filename,
        mode=mode,
        output=output,
    )
    return None
