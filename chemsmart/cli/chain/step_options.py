"""Parse program CLI options embedded in quoted chain ``-s`` steps."""

from dataclasses import replace

import click
from click.core import ParameterSource

_CHAIN_OWNED_PARAMS = frozenset(
    {
        "project",
        "filename",
        "label",
        "append_label",
        "index",
        "charge",
        "multiplicity",
        "pubchem",
        "record_index",
        "record_id",
        "structure_id",
        "structure_index",
        "molecule_id",
    }
)

_PARAM_TO_SETTING = {
    "additional_opt_options": "additional_opt_options_in_route",
}


def apply_step_option_tokens(step):
    """Attach settings overlay parsed from ``step.extra_option_tokens``."""
    if not step.extra_option_tokens:
        return step
    overlay, keywords = overlay_from_program_cli(
        step.program, step.extra_option_tokens
    )
    return replace(
        step,
        extra_settings=tuple(overlay.items()),
        extra_keywords=tuple(keywords),
    )


def overlay_from_program_cli(program, tokens):
    """Parse ``tokens`` with the program group CLI into a settings overlay."""
    command = _program_command(program)
    try:
        ctx = command.make_context(
            command.name, list(tokens), allow_extra_args=True
        )
    except click.ClickException as exc:
        raise ValueError(str(exc)) from exc
    if ctx.args:
        raise ValueError(
            f"Unrecognized arguments in chain step: {' '.join(ctx.args)}"
        )

    overlay = {}
    keywords = []
    remove_solvent = False
    for name, value in ctx.params.items():
        source = ctx.get_parameter_source(name)
        if source is not ParameterSource.COMMANDLINE:
            continue
        if name in _CHAIN_OWNED_PARAMS:
            raise ValueError(
                f"Set {name} on the chain command, not inside -s/--steps."
            )
        if name == "remove_solvent":
            remove_solvent = bool(value)
            continue
        setting_name = _PARAM_TO_SETTING.get(name, name)
        overlay[setting_name] = value
        keywords.append(setting_name)

    if remove_solvent:
        for setting_name in (
            "solvent_model",
            "solvent_id",
            "custom_solvent",
            "solventfilename",
        ):
            overlay[setting_name] = None
            if setting_name not in keywords:
                keywords.append(setting_name)
    return overlay, tuple(keywords)


def _program_command(program):
    from chemsmart.cli.crest.crest import crest
    from chemsmart.cli.gaussian.gaussian import gaussian
    from chemsmart.cli.orca.orca import orca
    from chemsmart.cli.xtb.xtb import xtb

    commands = {
        "crest": crest,
        "xtb": xtb,
        "gaussian": gaussian,
        "orca": orca,
    }
    command = commands.get(program)
    if command is None:
        raise ValueError(f"Unknown chain program {program!r}.")
    return command
