"""Helper functions and decorators for QMMM job support across Gaussian commands."""

import ast
import logging

import click

from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)


def create_qmmm_subcommand(parent_command):
    """
    Create a QMMM subcommand for a given parent Gaussian job command.

    This function creates a 'qmmm' subcommand that can be attached to any
    Gaussian job type (opt, ts, modred, sp, etc.) to enable QM/MM calculations.

    Args:
        parent_command: The Click group command to attach the qmmm subcommand to

    Returns:
        The qmmm command function
    """

    @parent_command.command("qmmm", cls=MyCommand)
    @click.option(
        "-hf",
        "--high-level-functional",
        type=str,
        default=None,
        help="Functional for high-level QM region.",
    )
    @click.option(
        "-hb",
        "--high-level-basis",
        type=str,
        default=None,
        help="Basis set for high-level QM region.",
    )
    @click.option(
        "-hff",
        "--high-level-force-field",
        type=str,
        default=None,
        help="Force field for high-level region.",
    )
    @click.option(
        "-mf",
        "--medium-level-functional",
        type=str,
        default=None,
        help="Functional for medium-level QM region.",
    )
    @click.option(
        "-mb",
        "--medium-level-basis",
        type=str,
        default=None,
        help="Basis set for medium-level QM region.",
    )
    @click.option(
        "-mff",
        "--medium-level-force-field",
        type=str,
        default=None,
        help="Force field for medium-level region.",
    )
    @click.option(
        "-lf",
        "--low-level-functional",
        type=str,
        default=None,
        help="Functional for low-level MM region.",
    )
    @click.option(
        "-lb",
        "--low-level-basis",
        type=str,
        default=None,
        help="Basis set for low-level MM region.",
    )
    @click.option(
        "-lff",
        "--low-level-force-field",
        type=str,
        default=None,
        help="Force field for low-level MM region.",
    )
    @click.option(
        "-ct",
        "--charge-total",
        type=int,
        default=None,
        help="Total system charge (1-indexed atoms)",
    )
    @click.option(
        "-mt",
        "--mult-total",
        type=int,
        default=None,
        help="Total system multiplicity (2S+1)",
    )
    @click.option(
        "-ci",
        "--charge-intermediate",
        type=int,
        default=None,
        help="Charge of intermediate (high+medium) system",
    )
    @click.option(
        "-mi",
        "--mult-intermediate",
        type=int,
        default=None,
        help="Multiplicity of intermediate (high+medium) system",
    )
    @click.option(
        "-ch",
        "--charge-high",
        type=int,
        default=None,
        help="Charge of high-level (model) system",
    )
    @click.option(
        "-mh",
        "--mult-high",
        type=int,
        default=None,
        help="Multiplicity of high-level (model) system",
    )
    @click.option(
        "-ha",
        "--high-level-atoms",
        type=str,
        default=None,
        help="Atom indices for high-level QM region (1-indexed, e.g., '1-10,15,20').",
    )
    @click.option(
        "-ma",
        "--medium-level-atoms",
        type=str,
        default=None,
        help="Atom indices for medium-level region (1-indexed).",
    )
    @click.option(
        "-la",
        "--low-level-atoms",
        type=str,
        default=None,
        help="Atom indices for low-level MM region (1-indexed).",
    )
    @click.option(
        "-ba",
        "--bonded-atoms",
        type=str,
        default=None,
        help="Bonded atom pairs at QM/MM boundary as string representation of dict.",
    )
    @click.option(
        "-sf",
        "--scale-factors",
        type=str,
        default=None,
        help="Scale factors for QM/MM calculations as string representation of dict.",
    )
    @click_job_options
    @click.pass_context
    def qmmm(
        ctx,
        high_level_functional,
        high_level_basis,
        high_level_force_field,
        medium_level_functional,
        medium_level_basis,
        medium_level_force_field,
        low_level_functional,
        low_level_basis,
        low_level_force_field,
        charge_total,
        mult_total,
        charge_intermediate,
        mult_intermediate,
        charge_high,
        mult_high,
        high_level_atoms,
        medium_level_atoms,
        low_level_atoms,
        bonded_atoms,
        scale_factors,
        **kwargs,
    ):
        """Run a QM/MM calculation for this job type."""

        from chemsmart.jobs.gaussian.qmmm import GaussianQMMMJob
        from chemsmart.jobs.gaussian.settings import GaussianQMMMJobSettings
        from chemsmart.utils.utils import convert_list_to_gaussian_frozen_list

        # Get parent command context
        jobrunner = ctx.obj["jobrunner"]
        project_settings = ctx.obj["project_settings"]
        label = ctx.obj.get("label")

        # Distinguish QMMM subcommand outputs from parent job outputs
        if label and "qmmm" not in label.lower():
            label = f"{label}_qmmm"

        # Infer jobtype from parent command name
        parent_cmd_name = ctx.parent.command.name if ctx.parent else None
        jobtype = parent_cmd_name  # The parent command name (opt, ts, sp, modred, scan, etc.)

        # Get parent command options
        # Subcommand value wins over parent; fall back to parent then default False
        skip_cli = kwargs.get("skip_completed", None)
        skip_completed = (
            skip_cli
            if skip_cli is not None
            else ctx.obj.get("parent_skip_completed", False)
        )
        freeze_atoms = ctx.obj.get("parent_freeze_atoms", None)
        parent_kwargs = ctx.obj.get("parent_kwargs", {})
        parent_settings = ctx.obj.get("parent_settings", None)
        parent_jobtype = ctx.obj.get("parent_jobtype", None)

        # Build QMMM settings in proper order:
        # 1. Start with project QMMM settings (from YAML)
        # 2. Merge parent job settings (to inherit modred, custom_solvent, etc.)
        # 3. Apply CLI options (done later in the code)

        # Step 1: Get project QMMM settings from YAML
        qmmm_settings = project_settings.qmmm_settings()
        if qmmm_settings is None:
            logger.debug(
                "Project qmmm settings not found; starting with GaussianQMMMJobSettings defaults."
            )
            qmmm_settings = GaussianQMMMJobSettings()
        else:
            logger.debug("Loaded project QMMM settings from YAML")

        # Step 2: Merge parent job settings to inherit job-specific attributes
        if parent_settings is not None:
            logger.debug(
                f"Merging parent settings ({parent_settings.__class__.__name__}) into QMMM settings"
            )

            # Define which attributes should ALWAYS be inherited from parent (if not None in parent)
            # These are job-specific attributes that come from the parent command or CLI
            always_inherit_from_parent = [
                "modred",  # From modred/scan parent command
                "custom_solvent",  # From parent settings or -A option
                "append_additional_info",  # From parent settings or -A option
                "additional_route_parameters",  # From parent settings or -r option
            ]
            qmmm_settings = qmmm_settings.merge(
                parent_settings, keywords=always_inherit_from_parent
            )

        #     for attr, value in parent_settings.__dict__.items():
        #         # Only copy if:
        #         # 1. QMMM settings has this attribute
        #         # 2. It's not a QMMM-specific attribute (those will be set from CLI or project)
        #         # 3. The value is not None (don't overwrite with None)
        #         if hasattr(qmmm_settings, attr) and value is not None:
        #             # Skip QMMM-specific attributes - these come from project or CLI
        #             qmmm_specific_attrs = {
        #                 'jobtype', 'high_level_functional', 'high_level_basis',
        #                 'high_level_force_field', 'medium_level_functional',
        #                 'medium_level_basis', 'medium_level_force_field',
        #                 'low_level_functional', 'low_level_basis',
        #                 'low_level_force_field', 'charge_total', 'real_multiplicity',
        #                 'int_charge', 'int_multiplicity', 'model_charge',
        #                 'model_multiplicity', 'high_level_atoms', 'medium_level_atoms',
        #                 'low_level_atoms', 'bonded_atoms', 'scale_factors'
        #             }
        #
        #             if attr not in qmmm_specific_attrs:
        #                 # Check if qmmm_settings already has a non-None value from project
        #                 current_value = getattr(qmmm_settings, attr, None)
        #                 # If project didn't set it (None or not set), inherit from parent
        #                 if current_value is None:
        #                     setattr(qmmm_settings, attr, value)
        #                     logger.debug(f"  Inherited {attr} from parent: {value}")
        #
        #     logger.debug(f"Inherited attributes from parent settings: modred={qmmm_settings.modred}, "
        #                 f"custom_solvent={qmmm_settings.custom_solvent}")
        # else:
        #     # No parent settings - try to merge with job_settings from file if available
        #     if job_settings is not None:
        #         logger.debug("No parent settings; merging job_settings from file")
        #         try:
        #             qmmm_settings = qmmm_settings.merge(job_settings, keywords=keywords)
        #         except Exception as exc:
        #             logger.debug(f"qmmm_settings.merge with job_settings failed: {exc}")

        # Set jobtype inferred from parent command
        if jobtype is not None:
            qmmm_settings.jobtype = jobtype
            logger.debug(f"Inferred jobtype from parent command: {jobtype}")

        # Step 3: Apply CLI options to qmmm_settings (these override project and parent)
        logger.debug("Applying QMMM-specific CLI options")
        if high_level_functional is not None:
            qmmm_settings.high_level_functional = high_level_functional
        if high_level_basis is not None:
            qmmm_settings.high_level_basis = high_level_basis
        if high_level_force_field is not None:
            qmmm_settings.high_level_force_field = high_level_force_field
        if medium_level_functional is not None:
            qmmm_settings.medium_level_functional = medium_level_functional
        if medium_level_basis is not None:
            qmmm_settings.medium_level_basis = medium_level_basis
        if medium_level_force_field is not None:
            qmmm_settings.medium_level_force_field = medium_level_force_field
        if low_level_functional is not None:
            qmmm_settings.low_level_functional = low_level_functional
        if low_level_basis is not None:
            qmmm_settings.low_level_basis = low_level_basis
        if low_level_force_field is not None:
            qmmm_settings.low_level_force_field = low_level_force_field
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
        if high_level_atoms is not None:
            qmmm_settings.high_level_atoms = high_level_atoms
        if medium_level_atoms is not None:
            qmmm_settings.medium_level_atoms = medium_level_atoms
        if low_level_atoms is not None:
            qmmm_settings.low_level_atoms = low_level_atoms
        if bonded_atoms is not None:
            qmmm_settings.bonded_atoms = bonded_atoms
        if scale_factors is not None:
            qmmm_settings.scale_factors = scale_factors

        # Get molecule
        molecules = ctx.obj["molecules"]
        molecule = molecules[-1]

        # Handle freeze_atoms from parent command
        if freeze_atoms is not None:
            frozen_atoms_list = get_list_from_string_range(freeze_atoms)
            logger.debug(f"Freezing atoms: {frozen_atoms_list}")
            molecule.frozen_atoms = convert_list_to_gaussian_frozen_list(
                frozen_atoms_list, molecule
            )

        # Attach QMMM parameters to molecule
        if high_level_atoms is not None:
            molecule.high_level_atoms = get_list_from_string_range(
                high_level_atoms
            )
        if medium_level_atoms is not None:
            molecule.medium_level_atoms = get_list_from_string_range(
                medium_level_atoms
            )
        if low_level_atoms is not None:
            molecule.low_level_atoms = get_list_from_string_range(
                low_level_atoms
            )
        if bonded_atoms is not None:
            molecule.bonded_atoms = ast.literal_eval(bonded_atoms)
        if scale_factors is not None:
            molecule.scale_factors = ast.literal_eval(scale_factors)

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

        if parent_jobtype is not None:
            qmmm_settings.parent_jobtype = parent_jobtype

        logger.info(f"QMMM job settings: {qmmm_settings.__dict__}")

        # Append "_qmmm" suffix to label for QMMM jobs
        if label and not label.endswith("_qmmm"):
            label = f"{label}_qmmm"

        return GaussianQMMMJob(
            molecule=molecule,
            settings=qmmm_settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **parent_kwargs,
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
        qs.charge_total = charge
    if mult is not None:
        qs.multiplicity = mult
        qs.mult_total = mult
