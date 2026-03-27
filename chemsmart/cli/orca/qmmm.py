"""Helper utilities for attaching ORCA QMMM subcommands to Click commands."""

import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import convert_string_to_slices

logger = logging.getLogger(__name__)


def create_orca_qmmm_subcommand(parent_command):
    """Attach the ORCA QMMM CLI (with full option set) to *parent_command*."""

    @parent_command.command("qmmm", cls=MyCommand)
    @click_job_options
    @click.option(
        "-j",
        "--jobtype",
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
        "-hx",
        "--high-level-functional",
        type=str,
        help="DFT functional for high-level (QM) region",
    )
    @click.option(
        "-hb",
        "--high-level-basis",
        type=str,
        help="Basis set for high-level (QM) region",
    )
    @click.option(
        "-ix",
        "--intermediate-level-functional",
        type=str,
        help="DFT functional for intermediate-level (QM2) region",
    )
    @click.option(
        "-ib",
        "--intermediate-level-basis",
        type=str,
        help="Basis set for intermediate-level (QM2) region",
    )
    @click.option(
        "-im",
        "--intermediate-level-method",
        type=str,
        help="Built-in method for intermediate-level region",
    )
    @click.option(
        "-lm",
        "--low-level-method",
        "--low-level-force-field",  # legacy alias
        type=str,
        help="Method/force field for low-level (MM) region",
    )
    @click.option(
        "-ha",
        "--high-level-atoms",
        type=str,
        help="High-level atom indices (e.g., '1-15,20')",
    )
    @click.option(
        "-ia",
        "--intermediate-level-atoms",
        type=str,
        help="Intermediate-level atom indices (e.g., '16-30')",
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
        type=int,
        help="Total system multiplicity",
    )
    @click.option(
        "-ci",
        "--charge-intermediate",
        type=int,
        help="Intermediate layer charge",
    )
    @click.option(
        "-mi",
        "--mult-intermediate",
        type=int,
        help="Intermediate layer multiplicity",
    )
    @click.option(
        "-ch",
        "--charge-high",
        type=int,
        help="High-level region charge",
    )
    @click.option(
        "-mh",
        "--mult-high",
        type=int,
        help="High-level region multiplicity",
    )
    @click.option(
        "-s",
        "--intermediate-level-solvation",
        type=str,
        help="Solvation model for intermediate-level region",
    )
    @click.option(
        "-ss",
        "--intermediate-solv-sheme",
        type=click.Choice(["CMCM_B", "CPCM_C"], case_sensitive=False),
        default="CMCM_B",
        show_default=True,
        help="Solvation model for intermediate-level region",
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
        "--high-level-h-bond-length",
        type=dict,
        help="Custom high-level-H bond lengths",
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
        jobtype,
        high_level_functional,
        high_level_basis,
        intermediate_level_functional,
        intermediate_level_basis,
        intermediate_level_method,
        low_level_method,
        high_level_atoms,
        intermediate_level_atoms,
        charge_total,
        mult_total,
        charge_intermediate,
        mult_intermediate,
        charge_high,
        mult_high,
        intermediate_level_solvation,
        active_atoms,
        use_active_info_from_pbc,
        optregion_fixed_atoms,
        high_level_h_bond_length,
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
        from chemsmart.jobs.orca.settings import ORCAQMMMJobSettings

        parent_settings = ctx.obj.get("parent_settings")
        parent_kwargs = ctx.obj.get("parent_kwargs", {})
        parent_skip_completed = ctx.obj.get("parent_skip_completed")
        parent_jobtype = ctx.obj.get("parent_jobtype")

        jobrunner = ctx.obj["jobrunner"]
        ctx.obj["qmmm"] = True
        project_settings = ctx.obj["project_settings"]
        logger.debug(
            "Project settings: %s", ctx.obj["project_settings"].__dict__
        )

        job_settings = ctx.obj["job_settings"]
        keywords = ctx.obj["keywords"]

        qmmm_settings = project_settings.qmmm_settings()
        if qmmm_settings is None:
            logger.warning(
                "Project qmmm settings not found; using ORCAQMMMJobSettings defaults."
            )
            qmmm_settings = ORCAQMMMJobSettings()

        try:
            qmmm_merged = qmmm_settings.merge(job_settings, keywords=keywords)
        except Exception as exc:
            logger.debug(
                "qmmm_settings.merge failed or is unavailable: %s", exc
            )
            if job_settings is not None:
                try:
                    qmmm_merged = ORCAQMMMJobSettings(
                        **getattr(job_settings, "__dict__", job_settings)
                    )
                except Exception:
                    qmmm_merged = qmmm_settings
            else:
                qmmm_merged = qmmm_settings

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

        # Distinguish QMMM subcommand outputs from parent job outputs
        if label and "qmmm" not in label.lower():
            label = f"{label}_qmmm"

        qmmm_settings = ORCAQMMMJobSettings(**qmmm_settings.__dict__)
        if parent_jobtype is not None:
            qmmm_settings.parent_jobtype = parent_jobtype

        if jobtype is not None:
            qmmm_settings.jobtype = jobtype
        if high_level_functional is not None:
            qmmm_settings.high_level_functional = high_level_functional
        if high_level_basis is not None:
            qmmm_settings.high_level_basis = high_level_basis
        if intermediate_level_functional is not None:
            qmmm_settings.intermediate_level_functional = (
                intermediate_level_functional
            )
        if intermediate_level_basis is not None:
            qmmm_settings.intermediate_level_basis = intermediate_level_basis
        if intermediate_level_method is not None:
            qmmm_settings.intermediate_level_method = intermediate_level_method
        if low_level_method is not None:
            # accept both new and legacy flag
            qmmm_settings.low_level_method = low_level_method
        if high_level_atoms is not None:
            qmmm_settings.high_level_atoms = high_level_atoms
        if intermediate_level_atoms is not None:
            qmmm_settings.intermediate_level_atoms = intermediate_level_atoms
        if charge_total is not None:
            qmmm_settings.charge_total = charge_total
        if mult_total is not None:
            qmmm_settings.mult_total = mult_total
        if charge_intermediate is not None:
            qmmm_settings.charge_intermediate = charge_intermediate
        if mult_intermediate is not None:
            qmmm_settings.mult_intermediate = mult_intermediate
        if charge_high is not None:
            qmmm_settings.charge_high = charge_high
        if mult_high is not None:
            qmmm_settings.mult_high = mult_high
        if intermediate_level_solvation is not None:
            qmmm_settings.intermediate_level_solvation = (
                intermediate_level_solvation
            )
        if active_atoms is not None:
            qmmm_settings.active_atoms = active_atoms
        if use_active_info_from_pbc is not None:
            qmmm_settings.use_active_info_from_pbc = use_active_info_from_pbc
        if optregion_fixed_atoms is not None:
            qmmm_settings.optregion_fixed_atoms = optregion_fixed_atoms
        if high_level_h_bond_length is not None:
            qmmm_settings.high_level_h_bond_length = high_level_h_bond_length
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
            qmmm_settings.scale_formal_charge_mm_atom = (
                scale_formal_charge_mm_atom
            )
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

        if parent_settings is not None:
            inherited_keywords = [
                "modred",
                "custom_solvent",
                "append_additional_info",
                "additional_route_parameters",
            ]
            try:
                qmmm_settings = qmmm_settings.merge(
                    parent_settings, keywords=inherited_keywords
                )
            except Exception as exc:
                logger.debug(
                    "Failed to merge parent settings into QMMM: %s", exc
                )

        effective_skip_completed = (
            skip_completed
            if skip_completed is not None
            else parent_skip_completed
        )

        qmmm_settings.re_init_and_validate()

        _populate_charge_and_multiplicity_on_settings(qmmm_settings)

        molecules = ctx.obj["molecules"]
        molecule = molecules[-1]
        logger.info(
            f"QMMM job settings from project: {qmmm_settings.__dict__}"
        )

        if high_level_atoms is not None:
            high_level_atoms_converted = convert_string_to_slices(
                high_level_atoms
            )
            molecule.high_level_atoms = high_level_atoms_converted
        if intermediate_level_atoms is not None:
            intermediate_level_atoms_converted = convert_string_to_slices(
                intermediate_level_atoms
            )
            molecule.intermediate_level_atoms = (
                intermediate_level_atoms_converted
            )
        if high_level_h_bond_length is not None:
            high_level_h_bond_length_dict = ast.literal_eval(
                high_level_h_bond_length
            )
            molecule.scale_factors = high_level_h_bond_length_dict

        combined_kwargs = {**parent_kwargs, **kwargs}

        from chemsmart.jobs.orca.qmmm import ORCAQMMMJob

        return ORCAQMMMJob(
            molecule=molecule,
            settings=qmmm_settings,
            label=label,
            skip_completed=effective_skip_completed,
            jobrunner=jobrunner,
            **combined_kwargs,
        )

    return qmmm


def _populate_charge_and_multiplicity_on_settings(qs):
    charge = getattr(qs, "charge", None)
    mult = getattr(qs, "multiplicity", None)

    if (
        getattr(qs, "charge_intermediate", None) is not None
        and getattr(qs, "mult_intermediate", None) is not None
    ):
        charge = qs.charge_intermediate
        mult = qs.mult_intermediate
    elif (
        getattr(qs, "charge_high", None) is not None
        and getattr(qs, "mult_high", None) is not None
    ):
        charge = qs.charge_high
        mult = qs.mult_high
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
