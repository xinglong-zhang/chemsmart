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
    help="Please specify the job type you want to run.",
)
@click.option(
    "-qx",
    "--qm-functional",
    type=str,
    help="Functional of QM region",
)
@click.option(
    "-qb",
    "--qm-basis",
    type=str,
    help="Basis set of QM region",
)
@click.option(
    "-qx2",
    "--qm2-functional",
    type=str,
    help="Functional of QM2 region",
)
@click.option(
    "-qb2",
    "--qm2-basis",
    type=str,
    help="Basis set of QM2 region",
)
@click.option(
    "-qm2m",
    "--qm2-methods",
    type=str,
    help="Methods of MM region",
)
@click.option(
    "-mf",
    "--mm-force-field",
    type=str,
    help="Force field of MM region",
)
@click.option(
    "-qa",
    "--qm-atoms",
    type=str,
    help="Indices of QM atoms",
)
@click.option(
    "-qa2",
    "--qm2-atoms",
    type=str,
    help="Indices of QM2 atoms",
)
@click.option(
    "-ct",
    "--charge-total",
    type=int,
    help="Charge of the total system",
)
@click.option(
    "-mt",
    "--mult-total",
    type=str,
    help="Multiplicity of the total system",
)
@click.option(
    "-cm",
    "--charge-medium",
    type=int,
    help="Charge of the medium system",
)
@click.option(
    "-mm",
    "--mult-medium",
    type=str,
    help="Multiplicity of the medium system",
)
@click.option(
    "-cq",
    "--charge-qm",
    type=int,
    help="Charge of the QM system",
)
@click.option(
    "-mq",
    "--mult-qm",
    type=str,
    help="Multiplicity of the QM system",
)
@click.option(
    "-s",
    "--qm2-solvation",
    type=str,
    help="Solvation model for QM2 level of theory",
)
@click.option(
    "-a",
    "--active-atoms",
    type=str,
    help="Indices of acitve atoms",
)
@click.option(
    "-ua",
    "--use-active-info-from-pbc",
    type=str,
    help="Get the definition of active atoms from the pdb file, default False",
)
@click.option(
    "-o",
    "--optregion-fixed-atoms",
    type=str,
    help="Indices of fixed atoms in the optimization region",
)
@click.option(
    "-h",
    "--qm-h-bond-length",
    type=dict,
    help="Customized bond length for two bonding atoms",
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
    help="Whether to use the converged charges",
)
@click.option(
    "-xn",
    "--conv-charges-max-n-cycles",
    type=int,
    help="Maximum number of cycles for charge convergence",
)
@click.option(
    "-t",
    "--conv-charges-conv-thresh",
    type=float,
    help="Convergence threshold for charge convergence",
)
@click.option(
    "-sc",
    "--scale-formal-charge-mm-atom",
    type=float,
    help="Scale the formal charge of MM atoms",
)
@click.option(
    "-nc",
    "--n-unit-cell-atoms",
    type=int,
    help="Number of atoms in the unit cell",
)
@click.option(
    "-ecp",
    "--ecp-layer-ecp",
    type=str,
    help="ECP layer for ECP atoms",
)
@click.option(
    "-ecpn",
    "--ecp-layer",
    type=int,
    help="number of cECP layers around the QM region",
)
@click.option(
    "-sc2",
    "--scale-formal-charge-ecp-atom",
    type=float,
    help="Scale the formal charge of ECP atoms",
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
    QM/MM calculation using ORCA.
    """
    from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

    jobrunner = ctx.obj["jobrunner"]
    project_settings = ctx.obj["project_settings"]
    qmmm_settings = project_settings.qmmm_settings()

    # job setting from filename or default, with updates from user in cli specified in keywords
    # e.g., `chemsmart sub orca qmmm -qf <qm_functional> -qb <qm_basis>`
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # merge project settings with job settings from cli keywords from cli.orca.py subcommands
    qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)

    # get label for the job
    label = ctx.obj["label"]
    logger.debug(f"Label for job: {label}")

    # convert from ORCAJobSettings instance to ORCAQMMMJobSettings instance
    qmmm_settings = ORCAQMMMJobSettings(**qmmm_settings.__dict__)
    qmmm_settings.jobtype = job_type
    qmmm_settings.qm_functional = qm_functional
    qmmm_settings.qm_basis = qm_basis
    qmmm_settings.qm2_functional = qm2_functional
    qmmm_settings.qm2_basis = qm2_basis
    qmmm_settings.qm2_methods = qm2_methods
    qmmm_settings.mm_force_field = mm_force_field
    qmmm_settings.qm_atoms = qm_atoms
    qmmm_settings.qm2_atoms = qm2_atoms
    qmmm_settings.charge_total = charge_total
    qmmm_settings.mult_total = mult_total
    qmmm_settings.charge_medium = charge_medium
    qmmm_settings.mult_medium = mult_medium
    qmmm_settings.charge_qm = charge_qm
    qmmm_settings.mult_qm = mult_qm
    qmmm_settings.qm2_solvation = qm2_solvation
    qmmm_settings.active_atoms = active_atoms
    qmmm_settings.use_active_info_from_pbc = use_active_info_from_pbc
    qmmm_settings.optregion_fixed_atoms = optregion_fixed_atoms
    qmmm_settings.qm_h_bond_length = qm_h_bond_length
    qmmm_settings.delete_la_double_counting = delete_la_double_counting
    qmmm_settings.delete_la_bond_double_counting_atoms = (
        delete_la_bond_double_counting_atoms
    )
    qmmm_settings.embedding_type = embedding_type
    qmmm_settings.conv_charges = conv_charges
    qmmm_settings.conv_charges_max_n_cycles = conv_charges_max_n_cycles
    qmmm_settings.conv_charges_conv_thresh = conv_charges_conv_thresh
    qmmm_settings.scale_formal_charge_mm_atom = scale_formal_charge_mm_atom
    qmmm_settings.n_unit_cell_atoms = n_unit_cell_atoms
    qmmm_settings.ecp_layer_ecp = ecp_layer_ecp
    qmmm_settings.ecp_layer = ecp_layer
    qmmm_settings.scale_formal_charge_ecp_atom = scale_formal_charge_ecp_atom
    # get molecule
    molecules = ctx.obj["molecules"]
    molecule = molecules[-1]
    logger.info(f"QMMM job settings from project: {qmmm_settings.__dict__}")
    # populate cli options by attaching QMMM parameters to the molecule
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
