"""
CLI interface for chemsmart project.

This module provides the main entry point for the chemsmart command-line
interface, organizing various subcommands and providing the ASCII art
banner display.
"""

import json
import sys
from pathlib import Path

import click

from chemsmart import __version__
from chemsmart.utils.cli import MyGroup

from .config import config
from .run import run
from .sub import sub
from .update import update


@click.group(cls=MyGroup)
@click.pass_context
@click.version_option(version=__version__, prog_name="CHEMSMART")
@click.option("--verbose", is_flag=True, default=True)
def entry_point(ctx, verbose):
    """
    Main entry point for the chemsmart CLI.
    """
    if verbose:
        debug = True
        stream = True
    else:
        debug = False
        stream = False
    # Set up logging
    from chemsmart.utils.logger import create_logger

    logger = create_logger(debug=debug, stream=stream)

    # ASCII Arts for CHEMSMART
    logger.info("\n")
    logger.info(
        "   "
        + " " * 25
        + "  ____ _   _ _____ __  __ ____  __  __    _    ____ _____ "
    )
    logger.info(
        "   "
        + " " * 25
        + " / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|"
    )
    logger.info(
        "   "
        + " " * 25
        + "| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || |  "
    )
    logger.info(
        "   "
        + " " * 25
        + "| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | |  "
    )
    logger.info(
        "   "
        + " " * 25
        + " \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_|  \n"
    )


entry_point.add_command(run)
entry_point.add_command(sub)
entry_point.add_command(config)
entry_point.add_command(update)


def _option_value(arguments, *names):
    """Return the last value supplied for one of ``names``."""
    value = None
    for position, argument in enumerate(arguments):
        if argument in names:
            if position + 1 < len(arguments):
                value = arguments[position + 1]
            continue
        for name in names:
            if name.startswith("--") and argument.startswith(f"{name}="):
                value = argument.split("=", 1)[1]
    return value


def _has_option(arguments, *names):
    return any(argument in names for argument in arguments)


def maybe_submit_deferred_npz(arguments):
    """Submit NPZ node batches before Click reads the archive.

    The submitted allocation streams structures from the archive into a
    bounded worker queue, so expansion happens on the compute system rather
    than in the user's login terminal.
    """
    if "sub" not in arguments:
        return False
    sub_position = arguments.index("sub")
    sub_arguments = list(arguments[sub_position + 1 :])
    program_position = next(
        (
            position
            for position, value in enumerate(sub_arguments)
            if value in {"xtb", "gaussian", "orca"}
        ),
        len(sub_arguments),
    )
    if program_position == len(sub_arguments):
        return False
    global_arguments = sub_arguments[:program_position]
    program_arguments = sub_arguments[program_position:]
    program = program_arguments[0]

    program_command = sub.commands[program]
    program_context = program_command.make_context(
        program, program_arguments[1:], resilient_parsing=False
    )
    program_params = program_context.params
    invoked_subcommand = (
        program_context.protected_args[0]
        if program_context.protected_args
        else None
    )
    if invoked_subcommand is None:
        return False

    nodes = _option_value(global_arguments, "-N", "--num-nodes")
    filename = program_params.get("filename")
    if (
        nodes is None
        or filename is None
        or not filename.lower().endswith(".npz")
    ):
        return False

    from chemsmart.settings.server import Server

    server_name = _option_value(global_arguments, "-s", "--server")
    server = (
        Server.current()
        if server_name is None
        else Server.from_servername(server_name)
    )
    if (server.scheduler or "").upper() != "SLURM":
        raise click.UsageError(
            "Deferred NPZ node batches currently require a SLURM server."
        )

    queue = (
        _option_value(global_arguments, "-q", "--queue") or server.queue_name
    )
    hours = _option_value(global_arguments, "-t", "--time-hours")
    hours = hours or server.num_hours
    label = program_params.get("label")
    append_label = program_params.get("append_label")
    if label is not None and append_label is not None:
        raise click.UsageError(
            "--label and --append-label cannot be used together."
        )
    if append_label is not None:
        label = f"{Path(filename).stem}_{append_label}"
    elif label is None:
        label = f"{Path(filename).stem}_{invoked_subcommand}"
    from chemsmart.utils.io import clean_label

    label = clean_label(label)[:100]
    test = _has_option(global_arguments, "--test")
    script = Path(f"{label}_deferred_npz.sh")

    if not test:
        server.ensure_batch_name_available(label)
    options_with_values = {
        "-N",
        "--num-nodes",
        "-t",
        "--time-hours",
        "-q",
        "--queue",
    }
    replay_global_arguments = list(global_arguments)
    position = 0
    while position < len(replay_global_arguments):
        if replay_global_arguments[position] in options_with_values:
            del replay_global_arguments[position : position + 2]
            continue
        if any(
            option.startswith("--")
            and replay_global_arguments[position].startswith(f"{option}=")
            for option in options_with_values
        ):
            del replay_global_arguments[position]
            continue
        if replay_global_arguments[position] in {
            "--test",
            "--no-test",
            "--print-command",
            "--no-print-command",
            "-v",
            "--verbose",
            "--no-verbose",
        }:
            del replay_global_arguments[position]
            continue
        position += 1
    task_arguments = [*replay_global_arguments, *program_arguments]
    manifest = Path(f"{label}_deferred_npz.json")
    manifest.write_text(
        json.dumps(
            {
                "filename": str(Path(filename).resolve()),
                "index": (
                    program_params.get("index")
                    or program_params.get("structure_index")
                    or ":"
                ),
                "base_label": label,
                "queue_capacity": 100,
                "num_nodes": int(nodes),
                "cores_per_node": server.num_cores,
                "memory_per_node_gb": server.mem_gb,
                "gpus_per_node": server.num_gpus,
                "task_cores": int(
                    _option_value(global_arguments, "-n", "--num-cores")
                    or server.num_cores
                ),
                "task_memory_gb": int(
                    _option_value(global_arguments, "-m", "--mem-gb")
                    or server.mem_gb
                ),
                "task_gpus": int(
                    _option_value(global_arguments, "-g", "--num-gpus")
                    or server.num_gpus
                ),
                "task_arguments": task_arguments,
            },
            indent=2,
        )
        + "\n"
    )
    command = [
        "python",
        "-c",
        "import sys; from chemsmart.jobs.batch import run_streaming_npz_batch; "
        "raise SystemExit(run_streaming_npz_batch(sys.argv[1]))",
        str(manifest),
    ]
    try:
        server.submit_streaming_npz_batch(
            script_path=script,
            label=label,
            num_nodes=int(nodes),
            program=program,
            command=command,
            queue=queue,
            hours=hours,
            test=test,
        )
    except RuntimeError as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"Deferred NPZ batch script: {script}")
    return True


def main():  # pragma: no cover
    """
    The main function executes on commands:
    ``python -m chemsmart`` and ``$ chemsmart``.

    This is your program's entry point.

    You can change this function to do whatever you want.
    Examples:
        * Run a test suite
        * Run a server
        * Do some other stuff
        * Run a command line application (Click, Typer, ArgParse)
        * List all available tasks
        * Run an application (Flask, FastAPI, Django, etc.)
    """
    if maybe_submit_deferred_npz(sys.argv[1:]):
        return
    obj = {}
    entry_point(obj=obj)


if __name__ == "__main__":
    main()
