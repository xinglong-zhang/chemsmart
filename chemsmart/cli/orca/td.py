"""ORCA TDDFT / TDA CLI command."""

import functools
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import click_orca_solvent_options, orca
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import (
    check_charge_and_multiplicity,
    get_list_from_string_range,
)

logger = logging.getLogger(__name__)


# Task keywords that switch a TD job from vertical to excited-state
# optimization / frequencies.  Matched case-insensitively as whole tokens.
_TD_TASK_TOKENS = {"opt", "freq", "numfreq"}


def click_orca_td_options(f):
    """Common TD-specific CLI options for the ORCA ``td`` command."""

    @click.option(
        "-n",
        "--nroots",
        type=int,
        default=3,
        show_default=True,
        help="Number of excited states to solve for (NRoots).",
    )
    @click.option(
        "--triplets/--no-triplets",
        default=None,
        help="Compute triplet excitations in addition to singlets. "
        "Implicitly enabled by --dosoc and by --root-mult triplet when "
        "not set.",
    )
    @click.option(
        "--tda/--no-tda",
        default=False,
        show_default=True,
        help="Enable the Tamm-Dancoff approximation. Default writes "
        "`TDA false`, requesting full TDDFT.",
    )
    @click.option(
        "--dosoc/--no-dosoc",
        default=None,
        help="Compute singlet-triplet SOC matrix elements (DoSOC). "
        "Only writes `%rel` if --soc-type is also given.",
    )
    @click.option(
        "--printlevel",
        type=int,
        default=None,
        help="TDDFT print level.",
    )
    @click.option(
        "--cpcmeq/--no-cpcmeq",
        default=None,
        help="Use equilibrium CPCM response for excited states (CPCMEQ). "
        "Left to ORCA's default when not set.",
    )
    @click.option(
        "--soc-type",
        type=int,
        default=None,
        help="SOC operator type written as `SOCType` in `%rel`. If omitted, "
        "ORCA uses its default operator type (3 in ORCA 6.1). This option "
        "does not enable SOC; use --dosoc to request SOC.",
    )
    @click.option(
        "--nto/--no-nto",
        default=None,
        help="Perform natural transition orbital analysis (DoNTO).",
    )
    @click.option(
        "--nto-states",
        type=str,
        default=None,
        help='State list for NTOStates, e.g. "1,2,3" or "1-3". Implicitly '
        "enables --nto if not set.",
    )
    @click.option(
        "--nto-thresh",
        type=float,
        default=None,
        help="Weight threshold for NTOs (NTOThresh).",
    )
    @click.option(
        "--td-maxiter",
        type=int,
        default=None,
        help="Maximum TD iterations (MaxIter inside %%tddft). Independent "
        "of the SCF MaxIter.",
    )
    @click.option(
        "--td-maxdim",
        type=int,
        default=None,
        help="Davidson subspace dimension multiplier (MaxDim).",
    )
    @click.option(
        "--td-etol",
        type=float,
        default=None,
        help="TD excitation-energy convergence tolerance in Hartree (Eh), "
        "written as ETol.",
    )
    @click.option(
        "--td-rtol",
        type=float,
        default=None,
        help="Residual convergence for TD (RTol).",
    )
    @click.option(
        "--tprint",
        type=float,
        default=None,
        help="Print threshold for transition amplitudes (TPrint).",
    )
    @click.option(
        "--root",
        type=int,
        default=None,
        help="Target excited-state index (IRoot). Defaults to 1 when an "
        "excited-state Opt/Freq task is requested via -r on the orca group. "
        "Must satisfy 1 <= root <= nroots.",
    )
    @click.option(
        "--root-mult",
        type=click.Choice(["singlet", "triplet"], case_sensitive=False),
        default=None,
        help="Multiplicity of the target excited state (IRootMult). Defaults "
        "to `singlet` for excited-state tasks; `triplet` implies --triplets.",
    )
    @click.option(
        "--follow-root/--no-follow-root",
        default=None,
        help="Follow the target root during excited-state optimization "
        "(FollowIRoot).",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def _parse_nto_states(value):
    """Parse comma/range list of positive state indices."""
    try:
        states = get_list_from_string_range(value)
    except Exception as exc:
        raise click.BadParameter(
            f"Invalid --nto-states value {value!r}: {exc}"
        ) from exc
    if not states or any(s < 1 for s in states):
        raise click.BadParameter(
            "--nto-states must contain positive 1-based state indices."
        )
    return states


def _split_task_tokens(additional_route_parameters):
    """Return ``(task_flags, remaining_extras)``.

    ``task_flags`` maps lower-cased keys ``{"opt", "freq", "numfreq"}`` to
    bool.  ``remaining_extras`` is the leftover ``-r`` string with the task
    tokens removed, or ``None`` if empty.
    """
    task_flags = {"opt": False, "freq": False, "numfreq": False}
    if additional_route_parameters is None:
        return task_flags, None
    raw = str(additional_route_parameters).strip()
    if not raw:
        return task_flags, None
    remaining_tokens = []
    for tok in raw.split():
        key = tok.lower()
        if key in _TD_TASK_TOKENS:
            task_flags[key] = True
        else:
            remaining_tokens.append(tok)
    remaining = " ".join(remaining_tokens) if remaining_tokens else None
    return task_flags, remaining


@orca.command("td", cls=MyCommand)
@click_job_options
@click_orca_solvent_options
@click_orca_td_options
@click.pass_context
def td(
    ctx,
    remove_solvent,
    solvent_model,
    solvent_id,
    solvent_options,
    solventfilename,
    nroots,
    triplets,
    tda,
    dosoc,
    printlevel,
    cpcmeq,
    soc_type,
    nto,
    nto_states,
    nto_thresh,
    td_maxiter,
    td_maxdim,
    td_etol,
    td_rtol,
    tprint,
    root,
    root_mult,
    follow_root,
    skip_completed,
    **kwargs,
):
    """Run TDDFT / TDA excited-state calculations with ORCA.

    Uses the ``td:`` section of the ORCA project YAML for method, basis,
    grid and solvent, plus the CLI options for the TD-specific block. By
    default this command writes ``NRoots 3`` and ``TDA false`` inside
    ``%tddft``; other TD options are omitted unless explicitly supplied or
    implicitly enabled.

    Excited-state optimization and frequencies are opt-in via the
    ``additional_route_parameters`` field, which can be set either on the
    project YAML (``td:`` section) or on the ``orca`` group via
    ``-r`` / ``--additional-route-parameters``. When neither the YAML nor the
    CLI requests ``Opt`` / ``Freq`` / ``NumFreq`` there, the CLI runs a
    vertical TD calculation and clears any Opt/Freq/NumFreq inherited from
    the ``ORCAJobSettings`` defaults or from parsing a ``.log`` / ``.inp``
    file. ``-r`` **replaces** the YAML value in full — passing ``-r ""``
    clears the extras, and passing ``-r TightSCF`` drops the YAML's Opt/Freq
    entirely; there is no per-token merge. ``-r`` does not replace any other
    structured setting.

    ``Opt`` / ``Freq`` / ``NumFreq`` act on the target excited state
    selected by ``--root`` / ``--root-mult`` (default: first singlet), not
    on the ground state. ``--root-mult triplet`` selects the target's
    multiplicity and implies ``--triplets``; that is distinct from
    ``--triplets`` alone, which only controls whether triplet excitations
    are solved for. ORCA's ``DoSOC`` computes singlet-triplet couplings
    and is not a SOC-based gradient method; ``Opt`` combined with
    ``--dosoc`` is refused because CHEMSMART does not currently drive
    ``SOCGrad`` — use ``chemsmart run orca inp`` for that workflow.

    \b
    Examples:
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 td
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 -r Opt        td --root 1
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 -r Freq       td --root 1
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 -r 'Opt Freq' td --root 1
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 -r 'Opt NumFreq' td --root 1
      chemsmart run orca -p td -f mol.xyz -c 0 -m 1 -r '' td   # force vertical
    """
    from chemsmart.jobs.orca.settings import ORCATDDFTJobSettings
    from chemsmart.jobs.orca.tddft import ORCATDDFTJob

    project_settings = ctx.obj["project_settings"]
    td_settings = project_settings.td_settings()
    if td_settings is None:
        raise click.UsageError(
            "The ORCA project YAML does not define a `td:` section. "
            "Add functional/basis/solvent for TD there before running "
            "`orca td`."
        )
    logger.debug(f"Loaded TDDFT settings from project: {td_settings}")

    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    td_settings = td_settings.merge(job_settings, keywords=keywords)

    td_settings.modify_solvent(
        remove_solvent=remove_solvent,
        solvent_model=solvent_model,
        solvent_id=solvent_id,
    )
    if solvent_options is not None:
        td_settings.additional_solvent_options = solvent_options
    if solventfilename is not None:
        td_settings.solventfilename = solventfilename

    # Split any Opt/Freq/NumFreq tokens out of the effective
    # additional_route_parameters (from YAML or from -r; -r fully replaces
    # the YAML value, no per-token merge) so we can validate + write them
    # once via ORCATDDFTJobSettings.
    task_flags, remaining_extras = _split_task_tokens(
        td_settings.additional_route_parameters
    )
    td_settings.additional_route_parameters = remaining_extras

    # Promote to TDDFT settings and clear any Opt/Freq/NumFreq inherited
    # from ORCAJobSettings defaults or from parsing a .log/.inp file — only
    # tasks requested via the effective additional_route_parameters apply.
    td_settings = ORCATDDFTJobSettings(**td_settings.__dict__)
    td_settings.jobtype = "td"
    td_settings.opt_excited = task_flags["opt"]
    td_settings.freq = task_flags["freq"]
    td_settings.numfreq = task_flags["numfreq"]

    if td_settings.freq and td_settings.numfreq:
        raise click.UsageError(
            "additional_route_parameters cannot request both Freq and "
            "NumFreq for a TD job."
        )

    # ---- validation & implicit-enable rules -------------------------------
    if nroots is not None and nroots < 1:
        raise click.BadParameter(
            "--nroots must be a positive integer.",
            param_hint="'-n' / '--nroots'",
        )

    if dosoc is True and triplets is False:
        raise click.UsageError(
            "--dosoc requires triplet excitations; remove --no-triplets or "
            "drop --dosoc."
        )
    effective_triplets = triplets
    if dosoc is True and effective_triplets is None:
        effective_triplets = True

    parsed_nto_states = None
    if nto_states is not None:
        parsed_nto_states = _parse_nto_states(nto_states)

    donto = nto
    nto_extra_given = parsed_nto_states is not None or nto_thresh is not None
    if nto is False and nto_extra_given:
        raise click.UsageError(
            "--no-nto conflicts with --nto-states/--nto-thresh; drop the "
            "extras or omit --no-nto."
        )
    if donto is None and nto_extra_given:
        donto = True

    # ---- excited-state target handling -----------------------------------
    excited_task = (
        td_settings.opt_excited or td_settings.freq or td_settings.numfreq
    )

    if td_settings.opt_excited and dosoc is True:
        raise click.UsageError(
            "-r Opt combined with --dosoc is not supported: ORCA's DoSOC "
            "computes SOC matrix elements, not gradients. SOC-based "
            "excited-state optimization requires SOCGrad; use "
            "`chemsmart run orca inp` for that workflow."
        )

    effective_root_mult = root_mult.lower() if root_mult is not None else None
    if excited_task and effective_root_mult is None:
        effective_root_mult = "singlet"

    if effective_root_mult == "triplet":
        if effective_triplets is False:
            raise click.UsageError(
                "--root-mult triplet requires triplet excitations; drop "
                "--no-triplets."
            )
        if effective_triplets is None:
            effective_triplets = True

    effective_root = root
    if excited_task and effective_root is None:
        effective_root = 1

    if effective_root is not None:
        if effective_root < 1:
            raise click.BadParameter(
                "--root must be a positive integer.",
                param_hint="'--root'",
            )
        if nroots is not None and effective_root > nroots:
            raise click.UsageError(
                f"--root ({effective_root}) exceeds --nroots ({nroots}); "
                "increase --nroots or lower --root."
            )

    if follow_root is not None and not excited_task:
        raise click.UsageError(
            "--follow-root only applies to excited-state Opt/Freq tasks; "
            "request one via `-r Opt` / `-r Freq` on the orca group."
        )

    td_settings.nroots = nroots
    td_settings.triplets = effective_triplets
    td_settings.tda = tda
    td_settings.dosoc = dosoc
    td_settings.printlevel = printlevel
    td_settings.cpcmeq = cpcmeq
    td_settings.soc_type = soc_type
    td_settings.donto = donto
    td_settings.ntostates = parsed_nto_states
    td_settings.ntothresh = nto_thresh
    td_settings.td_maxiter = td_maxiter
    td_settings.td_maxdim = td_maxdim
    td_settings.td_etol = td_etol
    td_settings.td_rtol = td_rtol
    td_settings.tprint = tprint
    td_settings.iroot = effective_root
    td_settings.iroot_mult = effective_root_mult
    td_settings.follow_iroot = follow_root

    logger.info(f"Final ORCA TDDFT settings: {td_settings.__dict__}")

    check_charge_and_multiplicity(td_settings)

    molecules = ctx.obj["molecules"]
    label = ctx.obj["label"]
    molecule_indices = ctx.obj["molecule_indices"]

    if len(molecules) > 1 and molecule_indices is not None:
        logger.info(f"Creating {len(molecules)} ORCA TDDFT jobs")
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            molecule_label = f"{label}_idx{idx}"
            jobs.append(
                ORCATDDFTJob(
                    molecule=molecule,
                    settings=td_settings,
                    label=molecule_label,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        return jobs

    molecule = molecules[-1]
    return ORCATDDFTJob(
        molecule=molecule,
        settings=td_settings,
        label=label,
        skip_completed=skip_completed,
        **kwargs,
    )
