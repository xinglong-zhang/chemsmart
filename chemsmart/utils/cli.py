"""
CLI utilities for ChemSmart.

Provides helpers and Click integrations used by the ChemSmart command line
tools, including:

- Custom Click `Group`/`Command` that persist subcommand metadata on the
  context for later reconstruction.
- Reconstruction of command-line arguments from a Click context object into a
  flat `list[str]` suitable for logging or re-execution.
- Job settings selection helpers for Gaussian and ORCA (opt, ts, modred, irc,
  scan, sp, td, wbi, nci).
- Coordinate/parameter validation helpers for scan and modred workflows with
  clear error messages.

Public API:
- `MyGroup`, `MyCommand`
- `CtxObjArguments`
- `get_setting_from_jobtype_for_gaussian`, `check_scan_coordinates_gaussian`
- `get_setting_from_jobtype_for_orca`, `check_scan_coordinates_orca`
"""

import ast
import copy
import logging
import pprint
import typing as t

import click
from click import Context

logger = logging.getLogger(__name__)


def _add_subcommand_info_to_ctx(ctx):
    """
    Add subcommand information to the Click context object.

    Extracts parameter information from the current command context and
    stores it in a structured format for later reconstruction of command
    line arguments.

    Args:
        ctx (Context): Click context object containing command information.
    """
    if "subcommand" not in ctx.obj:
        ctx.obj["subcommand"] = []

    subcommand = {}
    subcommand["name"] = ctx.info_name

    # Order of ctx.params may not be the same as order of ctx.command.params!
    nargs = {param.name: param.nargs for param in ctx.command.params}
    is_multiple = {param.name: param.multiple for param in ctx.command.params}
    types = {param.name: param.type for param in ctx.command.params}
    is_flag = {param.name: param.is_flag for param in ctx.command.params}
    secondary_opts = {
        param.name: param.secondary_opts for param in ctx.command.params
    }

    subcommand["kwargs"] = {}
    for param, value in ctx.params.items():
        subcommand["kwargs"][param] = {}
        subcommand["kwargs"][param]["value"] = value
        subcommand["kwargs"][param]["nargs"] = nargs[param]
        subcommand["kwargs"][param]["is_multiple"] = is_multiple[param]
        subcommand["kwargs"][param]["type"] = types[param]
        subcommand["kwargs"][param]["is_flag"] = is_flag[param]
        subcommand["kwargs"][param]["secondary_opts"] = secondary_opts[param]

    parent = ctx.parent.info_name if ctx.parent is not None else None

    subcommand["parent"] = parent

    ctx.obj["subcommand"] += [subcommand]


class MyGroup(click.Group):
    """
    Click Group that records invocation metadata in the context.

    On each invocation it appends this group's name and parameter values to
    `ctx.obj["subcommand"]`, so the full CLI can be reconstructed later
    (see `CtxObjArguments`). Works in tandem with `MyCommand` to capture the
    entire subcommand chain. Expects `ctx.obj` to be a mutable mapping.
    """

    def invoke(self, ctx: Context) -> t.Any:
        """
        Record group invocation metadata, then delegate to Click.

        Args:
            ctx (Context): Click context object.

        Returns:
            Any: Result of `click.Group.invoke`.
        """
        _add_subcommand_info_to_ctx(ctx)
        return super().invoke(ctx)


class MyCommand(click.Command):
    """
    Click Command that records invocation metadata in the context.

    On invocation it appends this command's name and parameter values to
    `ctx.obj["subcommand"]` so the full CLI can be reconstructed later
    (see `CtxObjArguments`). Typically used together with `MyGroup` to
    capture the entire command chain. Expects `ctx.obj` to be a mutable
    mapping.
    """

    def invoke(self, ctx):
        """
        Record command invocation metadata, then delegate to Click.

        Args:
            ctx (Context): Click context object.

        Returns:
            Any: Result of `click.Command.invoke`.
        """
        _add_subcommand_info_to_ctx(ctx)
        return super().invoke(ctx)


class CtxObjArguments:
    """
    Reconstruct command-line arguments from Click context metadata.

    Consumes the subcommand records produced by `MyGroup`/`MyCommand` and
    rebuilds a flat list of CLI arguments suitable for logging or
    re-execution.

    Attributes:
        commands (list[dict]): Shallow copy of the context subcommand list;
            each dict contains keys like "name", "parent", and "kwargs".
        entry_point (str | None): Optional explicit entry command name; if
            None, the command with no parent is used as the root.
    """

    def __init__(self, commands, entry_point=None):
        """
        Initialize the argument reconstruction manager.

        Args:
            commands (list): List of command dictionaries from context.
            entry_point (str, optional): Specific entry point command name.
        """
        self.commands = copy.copy(commands)
        self.entry_point = entry_point

    @staticmethod
    def _value(v):
        """
        Convert a parameter value to its CLI string representation.

        Booleans return an empty string because flags do not take values.

        Args:
            v: Parameter value of any type.

        Returns:
            str: String representation (empty for booleans).
        """
        if v is True or v is False:
            return ""
        return str(v)

    @staticmethod
    def argument_keyword(arg, value, is_flag, secondary_opts):
        """
        Generate the appropriate CLI option keyword for an argument.

        Handles flag options, boolean negation via secondary options, and
        short vs. long option formatting. Returns an empty string when a
        False-valued simple flag should be omitted entirely.

        Args:
            arg (str): Argument name.
            value: Argument value.
            is_flag (bool): Whether the argument is a flag.
            secondary_opts (list[str]): Secondary option names from Click.

        Returns:
            str: Formatted option keyword (may be empty string).
        """
        arg = arg.replace("_", "-")

        if value is False:
            if is_flag and len(secondary_opts) == 0:
                # Click "is_flag" option. Click param class has no way to check if
                # flag is True or not only can determine by checking arg_secondary_opts
                return ""

            # If variable is bool, then there will be non-empty secondary_opts
            if len(secondary_opts) == 1:
                arg = secondary_opts[0].strip("-")
            elif len(secondary_opts) > 1:
                arg = secondary_opts[-1].strip(
                    "--"
                )  # Instead of assert, use last option

        return "-" + arg if len(arg) == 1 else "--" + arg

    @property
    def _subcommands(self):
        """
        Get the list of subcommands.

        Returns:
            list[dict]: List of subcommand dictionaries.
        """
        return self.commands

    def _entry_point(self):
        """
        Find and return the entry point command.

        Locates the main entry point command either by specified name
        or by finding the command with no parent.

        Returns:
            dict: Entry point command dictionary.

        Raises:
            ValueError: If specified entry point is not found.
            AssertionError: If no entry point can be determined.
        """
        if self.entry_point is not None:
            subcommand = [
                s for s in self._subcommands if s["name"] == self.entry_point
            ]
            if len(subcommand) == 0:
                raise ValueError(
                    f"Specified entry point: {self.entry_point} could not be found"
                )
            assert len(subcommand) == 1
            return subcommand[0]

        for subcommand in self._subcommands:
            if subcommand["parent"] is None:
                return subcommand
        raise AssertionError("Could not find entry point")

    def _reconstruct_command(self, command):
        """
        Reconstruct command line arguments for a single command.

        Processes a command dictionary to generate the corresponding
        command line string with all arguments and options properly
        formatted.

        Args:
            command (dict): Command dictionary with name and kwargs.

        Returns:
            list[str]: Command line arguments for this command.
        """
        subcommand_args = command["kwargs"]
        subcommand_name = command["name"]
        logger.debug(f"Subcommand name {subcommand_name}")
        logger.debug(f"Subcommand args {pprint.pformat(subcommand_args)}")

        command_line_string = [subcommand_name]
        for k, subdict in subcommand_args.items():
            v = subdict["value"]
            multiple_arg = subdict["is_multiple"]
            arg_type = subdict["type"]
            arg_nargs = subdict["nargs"]
            arg_is_flag = subdict["is_flag"]
            arg_secondary_opts = subdict["secondary_opts"]

            if v is None or (isinstance(v, tuple) and len(v) == 0):
                continue

            keyword = self.argument_keyword(
                k, v, arg_is_flag, arg_secondary_opts
            )
            if multiple_arg:
                logger.debug(f"Keyword {keyword}, value {v}")
                for i in v:
                    command_line_string += [keyword]
                    if arg_nargs == 1:
                        command_line_string += [self._value(i)]
                    else:
                        command_line_string += [self._value(j) for j in i]
            else:
                command_line_string += [keyword]
                if arg_type.name == "literal_eval":
                    command_line_string += [self._value(v)]
                elif isinstance(v, (list, tuple)):
                    command_line_string += [self._value(i) for i in v]
                else:
                    command_line_string += [self._value(v)]
        logger.debug(f"Command line string: {command_line_string}")
        return command_line_string

    def _reconstruct_family(self, parent):
        """
        Reconstruct command line for a command family (parent and children).

        Recursively processes a parent command and all its child commands
        to generate a complete command line sequence.

        Args:
            parent (dict): Parent command dictionary.

        Returns:
            list[str]: Complete command line arguments for parent and children.
        """
        children = self._subcommands_of(parent)
        logger.debug(
            f'Parent: {parent["name"]};  children: {[c["name"] for c in children]}'
        )

        parent_command_line = self._reconstruct_command(parent)
        children_command_line = []
        for child in children:
            child_command_line = self._reconstruct_family(child)
            children_command_line.extend(child_command_line)
        return parent_command_line + children_command_line

    def _subcommands_of(self, parent):
        """
        Get all direct child subcommands of a parent command.

        Args:
            parent (dict): Parent command dictionary.

        Returns:
            list[dict]: List of child command dictionaries.
        """
        return [
            subcommand
            for subcommand in self._subcommands
            if subcommand["parent"] == parent["name"]
        ]

    def reconstruct_command_line(self):
        """
        Reconstruct the complete command line from stored context data.

        Processes all stored command information to generate a complete
        command line that can be used for job execution or logging.

        Returns:
            list[str]: Complete command line arguments as list of strings.
        """
        parent = self._entry_point()
        command_line = self._reconstruct_family(parent)
        # logger.info(f'Cmd: {command_line}')
        return [i for i in command_line if len(i) != 0]


def get_setting_from_jobtype_for_gaussian(
    project_settings, jobtype, coordinates, step_size, num_steps
):
    """
    Get Gaussian job settings based on job type and parameters.

    Retrieves appropriate settings configuration for different types of
    Gaussian calculations including optimization, transition state search,
    IRC calculations, and coordinate scans.

    Args:
        project_settings: Gaussian project settings object.
        jobtype (str): Type of calculation (opt, ts, modred, irc, scan, etc.).
        coordinates: Coordinate specification for modred/scan jobs.
        step_size: Step size for scan calculations.
        num_steps: Number of steps for scan calculations.

    Returns:
        GaussianJobSettings or None: Configured settings object for the job
        type (or subclass), or None if the jobtype is unsupported.

    Raises:
        ValueError: If jobtype is None.
        AssertionError: If required parameters are missing for specific job types.
    """
    if jobtype is None:
        raise ValueError("Jobtype must be provided for Crest and Link job.")

    settings = None

    if jobtype.lower() == "opt":
        settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for modred job."
        settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        check_scan_coordinates_gaussian(coordinates, step_size, num_steps)
        settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        settings = project_settings.nci_settings()

    if coordinates is not None:
        modred_info = ast.literal_eval(coordinates)
        if jobtype == "modred":
            settings.modred = modred_info
        elif jobtype == "scan":
            num_steps_info = ast.literal_eval(num_steps)
            step_size_info = ast.literal_eval(step_size)
            if not isinstance(num_steps_info, list):
                num_steps_info = [num_steps_info]
            if not isinstance(step_size_info, list):
                step_size_info = [step_size_info]
            check_scan_parameters_consistency(
                modred_info, step_size_info, num_steps_info
            )
            scan_info = {
                "coords": modred_info,
                "num_steps": num_steps_info,
                "step_size": step_size_info,
            }
            settings.modred = scan_info

    return settings


def check_scan_coordinates_gaussian(coordinates, step_size, num_steps):
    """
    Validate scan coordinate parameters for Gaussian calculations.

    Ensures all required parameters for coordinate scan calculations are
    provided and non-None. Provides detailed error messages with usage
    examples for missing parameters.

    Args:
        coordinates: Coordinate specification for the scan.
        step_size: Step size for the scan.
        num_steps: Number of scan steps.

    Raises:
        AssertionError: If any required parameter is None, with detailed
            usage instructions.
    """
    assert all(v is not None for v in [coordinates, step_size, num_steps]), (
        "Scanning coordinates, step size and number of steps of scan required!\n"
        f"But is coordinates: {coordinates}, step_size: {step_size}, num_steps: {num_steps}\n"
        "Use the flags `-c -s -n` for coordinates, step-size and num-steps respectively.\n"
        "Example usage: `-c [[2,3],[6,7]] -s 0.1 -n 15`"
    )


def check_scan_parameters_consistency(
    coords: list, step_size: list, num_steps: list
):
    """
    Validate consistency between scan coordinates and their parameters.

    Ensures that the number of coordinates matches the number of step_size
    and num_steps values provided. Handles both single coordinate and
    multi-coordinate scan scenarios.

    Args:
        coords: Coordinate specification - list[int] for single coordinate
                or list[list[int]] for multiple coordinates.
        step_size: List of step sizes - must match number of coordinates.
        num_steps: List of step counts - must match number of coordinates.

    Raises:
        ValueError: If the number of coordinates doesn't match the number
                   of step_size or num_steps values, with detailed error
                   message showing the mismatch.

    Examples:
        # Valid single coordinate
        check_scan_parameters_consistency([1,2], [0.1], [10])

        # Valid multiple coordinates
        check_scan_parameters_consistency([[1,2],[3,4]], [0.1,0.2], [10,15])

        # Invalid - mismatched counts
        check_scan_parameters_consistency([[1,2],[3,4]], [0.1], [10,15])
        # Raises ValueError
    """
    if isinstance(coords[0], list):
        coords_num = len(coords)
    elif isinstance(coords[0], int):
        coords_num = 1
    else:
        raise ValueError("Invalid format for coordinates.")
    if len(step_size) != coords_num or len(num_steps) != coords_num:
        raise ValueError(
            f"Mismatch in number of coordinates and step parameters: "
            f"coords={coords_num}, step_size={len(step_size)}, num_steps={len(num_steps)}"
        )


def get_setting_from_jobtype_for_orca(
    project_settings, jobtype, coordinates, dist_start, dist_end, num_steps
):
    """
    Get ORCA job settings based on job type and parameters.

    Retrieves appropriate settings configuration for different types of
    ORCA calculations including optimization, transition state search,
    IRC calculations, and coordinate scans with distance specifications.

    Args:
        project_settings: ORCA project settings object.
        jobtype (str): Type of calculation (opt, ts, modred, irc, scan, etc.).
        coordinates: Coordinate specification for modred/scan jobs.
        dist_start (float): Starting distance for scan calculations.
        dist_end (float): Ending distance for scan calculations.
        num_steps (int): Number of steps for scan calculations.

    Returns:
        ORCAJobSettings or None: Configured settings object for the job type
        (or subclass), or None if the jobtype is unsupported.

    Raises:
        ValueError: If jobtype is None.
        AssertionError: If required parameters are missing for specific job types.
    """
    if jobtype is None:
        raise ValueError("Jobtype must be provided for Crest and Link job.")

    settings = None

    if jobtype.lower() == "opt":
        settings = project_settings.opt_settings()
    elif jobtype.lower() == "ts":
        settings = project_settings.ts_settings()
    elif jobtype.lower() == "modred":
        assert (
            coordinates is not None
        ), "Coordinates must be provided for modred job."
        settings = project_settings.modred_settings()
    elif jobtype.lower() == "irc":
        settings = project_settings.irc_settings()
    elif jobtype.lower() == "scan":
        check_scan_coordinates_orca(
            coordinates, dist_start, dist_end, num_steps
        )
        settings = project_settings.scan_settings()
    elif jobtype.lower() == "sp":
        settings = project_settings.sp_settings()
    elif jobtype.lower() == "td":
        settings = project_settings.td_settings()
    elif jobtype.lower() == "wbi":
        settings = project_settings.wbi_settings()
    elif jobtype.lower() == "nci":
        settings = project_settings.nci_settings()

    if coordinates is not None:
        modred_info = ast.literal_eval(coordinates)
        if jobtype == "modred":
            settings.modred = modred_info
        elif jobtype == "scan":
            scan_info = {
                "coords": modred_info,
                "dist_start": float(dist_start),
                "dist_end": float(dist_end),
                "num_steps": int(num_steps),
            }
            settings.modred = scan_info

    return settings


def check_scan_coordinates_orca(coordinates, dist_start, dist_end, num_steps):
    """
    Validate scan coordinate parameters for ORCA calculations.

    Ensures all required parameters for ORCA coordinate scan calculations
    are provided and non-None. Provides detailed error messages with usage
    examples and index conversion notes.

    Args:
        coordinates: Coordinate specification for the scan.
        dist_start (float): Starting distance for the scan.
        dist_end (float): Ending distance for the scan.
        num_steps (int): Number of scan steps.

    Raises:
        AssertionError: If any required parameter is None, with detailed
            usage instructions and indexing information.
    """
    assert all(
        v is not None for v in [coordinates, dist_start, dist_end, num_steps]
    ), (
        "Scanning coordinates, starting distance, ending distance and number of steps of scan required!\n"
        "Use flags `-c -a -b -n` for coordinates, starting distance, ending distance and num-steps respectively.\n"
        "Example usage: `-c [[2,3],[6,7]] -x 3.0 -y 1.2 -n 15` to scan the distance between atom 2 and atom 3 "
        "and distance between atom 6 and 7 from distance 3.0 Angstrom to 1.2 Angstrom in 15 points.\n "
        "Note: all indices should be 1-indexed. Chemsmart has already taken care of converting 0-indexed "
        "(used in ORCA) to 1-indexed (used in visualization software such as PyMOL, Gaussview, etc)."
    )


def update_irc_label(label, direction, flat_irc):
    """
    Update the job label based on IRC direction and flat IRC flag.

    Appends 'f' for forward direction, 'r' for reverse direction,
    and '_flat' if flat_irc is True.

    Args:
        label (str): Original job label.
        direction (str | None): IRC direction ('forward' or 'reverse').
        flat_irc (bool): Whether the IRC is flat.

    Returns:
        str: Updated job label with appropriate suffixes.
    """
    if direction is not None:
        if direction.lower() == "forward":
            label += "f"
        elif direction.lower() == "reverse":
            label += "r"
        else:
            raise ValueError(
                "Invalid direction for IRC job. Must be 'forward' or 'reverse'."
            )
    if flat_irc:
        label += "_flat"
    return label
