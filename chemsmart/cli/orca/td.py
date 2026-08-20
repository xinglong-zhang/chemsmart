"""Canonical ORCA vertical CIS/TD-DFT command."""

import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_solvent_options, orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)


@orca.command("td", cls=MyCommand)
@click_job_options
@click_orca_solvent_options
@click.option(
    "--response-method",
    type=click.Choice(("tda", "tddft"), case_sensitive=False),
    default=None,
    help="Excited-state response approximation; defaults to project YAML.",
)
@click.option(
    "--nstates",
    type=click.IntRange(min=1),
    default=None,
    help="Number of excited-state roots; defaults to project YAML.",
)
@click.option(
    "--state-manifold",
    type=click.Choice(("singlet", "singlet_triplet"), case_sensitive=False),
    default=None,
    help=(
        "Closed-shell singlets, or singlets together with spin-adapted "
        "triplets; defaults to project YAML."
    ),
)
@click.pass_context
def td(
    ctx,
    response_method,
    nstates,
    state_manifold,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    solventfilename,
    skip_completed,
    **kwargs,
):
    """Run a fixed-geometry ORCA electronic excitation calculation."""

    project_settings = ctx.obj["project_settings"]
    settings = project_settings.td_settings()
    settings = settings.merge(
        ctx.obj["job_settings"], keywords=ctx.obj["keywords"]
    )
    settings.jobtype = "td"
    settings.freq = False
    if response_method is not None:
        settings.response_method = response_method.lower()
    if nstates is not None:
        settings.nstates = nstates
    if state_manifold is not None:
        settings.state_manifold = state_manifold.lower()
    settings.modify_solvent(
        remove_solvent=remove_solvent,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
    )
    if solvent_options is not None:
        settings.additional_solvent_options = solvent_options
    if solventfilename is not None:
        settings.solventfilename = solventfilename

    missing = [
        name
        for name in ("response_method", "nstates", "state_manifold")
        if getattr(settings, name, None) is None
    ]
    if missing:
        raise click.UsageError(
            "ORCA td requires project or CLI values for " + ", ".join(missing)
        )
    if settings.state_manifold in {"singlet", "singlet_triplet"} and (
        settings.multiplicity is not None and int(settings.multiplicity) != 1
    ):
        raise click.UsageError(
            "ORCA closed-shell singlet TD roots require multiplicity 1"
        )
    check_charge_and_multiplicity(settings)

    from chemsmart.jobs.orca.td import ORCATDDFTJob
    from chemsmart.utils.cli import create_sp_label

    molecules = ctx.obj["molecules"]
    indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]
    if len(molecules) > 1 and indices is not None:
        return [
            ORCATDDFTJob(
                molecule=molecule,
                settings=settings,
                label=create_sp_label(f"{label}_idx{index}", settings),
                skip_completed=skip_completed,
                **kwargs,
            )
            for molecule, index in zip(molecules, indices)
        ]
    return ORCATDDFTJob(
        molecule=molecules[-1],
        settings=settings,
        label=create_sp_label(label, settings),
        skip_completed=skip_completed,
        **kwargs,
    )
