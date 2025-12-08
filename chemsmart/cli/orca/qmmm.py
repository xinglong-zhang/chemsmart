import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.utils.cli import (
    MyCommand,
)
from chemsmart.utils.utils import convert_string_to_slices

logger = logging.getLogger(__name__)


def _populate_charge_and_multiplicity_on_settings(qs):
    """
    Populate top-level charge/multiplicity from QMMM-specific fields.

    Sets the parent class charge and multiplicity attributes from
    layer-specific values following the preference order:
    1. medium (charge_medium, mult_medium)
    2. qm (charge_qm, mult_qm)
    3. total (charge_total, mult_total)

    This ensures proper charge/multiplicity values are available
    for downstream ORCA input file writers.

    Args:
        qs: ORCAQMMMJobSettings instance to update
    """
    charge = getattr(qs, "charge", None)
    mult = getattr(qs, "multiplicity", None)

    # Preference: medium -> qm -> total
    if (
        getattr(qs, "charge_medium", None) is not None
        and getattr(qs, "mult_medium", None) is not None
    ):
        charge = qs.charge_medium
        mult = qs.mult_medium
    elif (
        getattr(qs, "charge_qm", None) is not None
        and getattr(qs, "mult_qm", None) is not None
    ):
        charge = qs.charge_qm
        mult = qs.mult_qm
    elif (
        getattr(qs, "charge_total", None) is not None
        and getattr(qs, "mult_total", None) is not None
    ):
        charge = qs.charge_total
        mult = qs.mult_total

    if charge is not None:
        qs.charge = charge
    if mult is not None:
        qs.multiplicity = mult


@orca.command("qmmm", cls=MyCommand)
@click_job_options
@click.option(
    "-j",
    "--job-type",
    type=click.Choice(
        [
            "QMMM",
            "QM/QM2",
            "QM/QM2/MM",
            "MOL-CRYSTAL-QMMM",
            "IONIC-CRYSTAL-QMMM",
        ],
        case_sensitive=False,
    ),
    help="Multiscale calculation type",
)
@click.option(
    "-qx",
    "--qm-functional",
    type=str,
    help="DFT functional for QM region",
)
@click.option(
    "-qb",
    "--qm-basis",
    type=str,
    help="Basis set for QM region",
)
@click.option(
    "-qx2",
    "--qm2-functional",
    type=str,
    help="DFT functional for QM2 region",
)
@click.option(
    "-qb2",
    "--qm2-basis",
    type=str,
    help="Basis set for QM2 region",
)
@click.option(
    "-qm2m",
    "--qm2-methods",
    type=str,
    help="Built-in method for QM2 region",
)
@click.option(
    "-mf",
    "--mm-force-field",
    type=str,
    help="Force field for MM region",
)
@click.option(
    "-qa",
    "--qm-atoms",
    type=str,
    help="QM atom indices (e.g., '1-15,20')",
)
@click.option(
    "-qa2",
    "--qm2-atoms",
    type=str,
    help="QM2 atom indices (e.g., '16-30')",
)
@click.option(
    "-ct",
    "--charge-total",
    type=int,
    help="Total system charge",
)
@click.option(
    "-mt",
    "--mult-total",
    type=str,
    help="Total system multiplicity",
)
@click.option(
    "-cm",
    "--charge-medium",
    type=int,
    help="Medium layer charge",
)
@click.option(
    "-mm",
    "--mult-medium",
    type=str,
    help="Medium layer multiplicity",
)
@click.option(
    "-cq",
    "--charge-qm",
    type=int,
    help="QM region charge",
)
@click.option(
    "-mq",
    "--mult-qm",
    type=str,
    help="QM region multiplicity",
)
@click.option(
    "-s",
    "--qm2-solvation",
    type=str,
    help="Solvation model for QM2 region",
)
@click.option(
    "-a",
    "--active-atoms",
    type=str,
    help="Active atom indices for optimization",
)
@click.option(
    "-ua",
    "--use-active-info-from-pbc",
    type=str,
    help="Use active atom info from PDB file",
)
@click.option(
    "-o",
    "--optregion-fixed-atoms",
    type=str,
    help="Fixed atom indices in optimization",
)
@click.option(
    "-h",
    "--qm-h-bond-length",
    type=dict,
    help="Custom QM-H bond lengths",
)
@click.option(
    "-d",
    "--delete-la-double-counting",
    type=bool,
    help="Whether to neglect bonds (QM1-MM1)",
)
@click.option(
    "-db",
    "--delete-la-bond-double-counting_atoms",
    type=bool,
    help="Whether to neglect bonds (QM1-MM1)",
)
@click.option(
    "-e",
    "--embedding-type",
    type=str,
    help="whether to use electronic (default)/mechanical embedding between QM and MM region",
)
@click.option(
    "-cc",
    "--conv-charges",
    type=bool,
    help="Use converged charges for crystal QM/MM",
)
@click.option(
    "-xn",
    "--conv-charges-max-n-cycles",
    type=int,
    help="Max cycles for charge convergence",
)
@click.option(
    "-t",
    "--conv-charges-conv-thresh",
    type=float,
    help="Charge convergence threshold",
)
@click.option(
    "-sc",
    "--scale-formal-charge-mm-atom",
    type=float,
    help="MM atom charge scaling factor",
)
@click.option(
    "-nc",
    "--n-unit-cell-atoms",
    type=int,
    help="Atoms per unit cell (MOL-CRYSTAL-QMMM)",
)
@click.option(
    "-ecp",
    "--ecp-layer-ecp",
    type=str,
    help="ECP type for boundary region",
)
@click.option(
    "-ecpn",
    "--ecp-layer",
    type=int,
    help="Number of ECP layers around QM region",
)
@click.option(
    "-sc2",
    "--scale-formal-charge-ecp-atom",
    type=float,
    help="ECP atom charge scaling factor",
)
@click.pass_context
def qmmm(
    ctx,
    job_type,
    qm_functional,
    qm_basis,
    qm2_functional,
    qm2_basis,
    qm2_methods,
    mm_force_field,
    qm_atoms,
    qm2_atoms,
    charge_total,
    mult_total,
    charge_medium,
    mult_medium,
    charge_qm,
    mult_qm,
    qm2_solvation,
    active_atoms,
    use_active_info_from_pbc,
    optregion_fixed_atoms,
    qm_h_bond_length,
    delete_la_double_counting,
    delete_la_bond_double_counting_atoms,
    embedding_type,
    conv_charges,
    conv_charges_max_n_cycles,
    conv_charges_conv_thresh,
    scale_formal_charge_mm_atom,
    n_unit_cell_atoms,
    ecp_layer_ecp,
    ecp_layer,
    scale_formal_charge_ecp_atom,
    skip_completed,
    **kwargs,
):
    """
    Set up and submit ORCA multiscale QM/MM calculations.

    Supports additive QM/MM, subtractive ONIOM schemes, and crystal QM/MM
    methods. Configures calculation parameters, validates settings, and
    creates ORCAQMMMJob instance for execution.
    """
    from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

    # Get execution context and initialize QM/MM flag
    jobrunner = ctx.obj["jobrunner"]
    ctx.obj["qmmm"] = True
    project_settings = ctx.obj["project_settings"]
    logger.debug("Project settings: %s", ctx.obj["project_settings"].__dict__)

    # Get base job settings and CLI keywords
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # Initialize QMMM settings from project or defaults
    qmmm_settings = project_settings.qmmm_settings()
    if qmmm_settings is None:
        logger.warning(
            "Project qmmm settings not found; using ORCAQMMMJobSettings defaults."
        )
        qmmm_settings = ORCAQMMMJobSettings()

    # Merge settings with error handling
    try:
        qmmm_merged = qmmm_settings.merge(job_settings, keywords=keywords)
    except Exception as exc:
        logger.debug("qmmm_settings.merge failed or is unavailable: %s", exc)
        # Fallback to job_settings or keep defaults
        if job_settings is not None:
            try:
                qmmm_merged = ORCAQMMMJobSettings(
                    **getattr(job_settings, "__dict__", job_settings)
                )
            except Exception:
                qmmm_merged = qmmm_settings
        else:
            qmmm_merged = qmmm_settings

    # Ensure final settings is ORCAQMMMJobSettings instance
    if isinstance(qmmm_merged, ORCAQMMMJobSettings):
        qmmm_settings = qmmm_merged
    else:
        try:
            qmmm_settings = ORCAQMMMJobSettings(
                **getattr(qmmm_merged, "__dict__", {})
            )
        except Exception:
            qmmm_settings = ORCAQMMMJobSettings()

    label = ctx.obj.get("label")
    logger.debug("Label for job: %s", label)

    # Convert settings to ORCAQMMMJobSettings if needed
    qmmm_settings = ORCAQMMMJobSettings(**qmmm_settings.__dict__)

    # Update settings with CLI-provided options (only if specified)
    if job_type is not None:
        qmmm_settings.jobtype = job_type
    if qm_functional is not None:
        qmmm_settings.qm_functional = qm_functional
    if qm_basis is not None:
        qmmm_settings.qm_basis = qm_basis
    if qm2_functional is not None:
        qmmm_settings.qm2_functional = qm2_functional
    if qm2_basis is not None:
        qmmm_settings.qm2_basis = qm2_basis
    if qm2_methods is not None:
        qmmm_settings.qm2_methods = qm2_methods
    if mm_force_field is not None:
        qmmm_settings.mm_force_field = mm_force_field
    if qm_atoms is not None:
        qmmm_settings.qm_atoms = qm_atoms
    if qm2_atoms is not None:
        qmmm_settings.qm2_atoms = qm2_atoms
    if charge_total is not None:
        qmmm_settings.charge_total = charge_total
    if mult_total is not None:
        qmmm_settings.mult_total = mult_total
    if charge_medium is not None:
        qmmm_settings.charge_medium = charge_medium
    if mult_medium is not None:
        qmmm_settings.mult_medium = mult_medium
    if charge_qm is not None:
        qmmm_settings.charge_qm = charge_qm
    if mult_qm is not None:
        qmmm_settings.mult_qm = mult_qm
    if qm2_solvation is not None:
        qmmm_settings.qm2_solvation = qm2_solvation
    if active_atoms is not None:
        qmmm_settings.active_atoms = active_atoms
    if use_active_info_from_pbc is not None:
        qmmm_settings.use_active_info_from_pbc = use_active_info_from_pbc
    if optregion_fixed_atoms is not None:
        qmmm_settings.optregion_fixed_atoms = optregion_fixed_atoms
    if qm_h_bond_length is not None:
        qmmm_settings.qm_h_bond_length = qm_h_bond_length
    if delete_la_double_counting is not None:
        qmmm_settings.delete_la_double_counting = delete_la_double_counting
    if delete_la_bond_double_counting_atoms is not None:
        qmmm_settings.delete_la_bond_double_counting_atoms = (
            delete_la_bond_double_counting_atoms
        )
    if embedding_type is not None:
        qmmm_settings.embedding_type = embedding_type
    if conv_charges is not None:
        qmmm_settings.conv_charges = conv_charges
    if conv_charges_max_n_cycles is not None:
        qmmm_settings.conv_charges_max_n_cycles = conv_charges_max_n_cycles
    if conv_charges_conv_thresh is not None:
        qmmm_settings.conv_charges_conv_thresh = conv_charges_conv_thresh
    if scale_formal_charge_mm_atom is not None:
        qmmm_settings.scale_formal_charge_mm_atom = scale_formal_charge_mm_atom
    if n_unit_cell_atoms is not None:
        qmmm_settings.n_unit_cell_atoms = n_unit_cell_atoms
    if ecp_layer_ecp is not None:
        qmmm_settings.ecp_layer_ecp = ecp_layer_ecp
    if ecp_layer is not None:
        qmmm_settings.ecp_layer = ecp_layer
    if scale_formal_charge_ecp_atom is not None:
        qmmm_settings.scale_formal_charge_ecp_atom = (
            scale_formal_charge_ecp_atom
        )

    # Set top-level charge/multiplicity from layer-specific values
    _populate_charge_and_multiplicity_on_settings(qmmm_settings)

    # Configure molecule with QM/MM partition parameters
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"QMMM job settings from project: {qmmm_settings.__dict__}")

    if qm_atoms is not None:
        qm_atoms = convert_string_to_slices(qm_atoms)
        molecule.high_level_atoms = qm_atoms
    if qm2_atoms is not None:
        qm2_atoms = convert_string_to_slices(qm2_atoms)
        molecule.medium_level_atoms = qm2_atoms
    if qm_h_bond_length is not None:
        qm_h_bond_length = ast.literal_eval(qm_h_bond_length)
        molecule.scale_factors = qm_h_bond_length

    from chemsmart.jobs.orca.qmmm import ORCAQMMMJob

    return ORCAQMMMJob(
        molecule=molecule,
        settings=qmmm_settings,
        label=label,
        skip_completed=skip_completed,
        jobrunner=jobrunner,
        # **ctx.obj["kwargs"],
        **kwargs,
    )
