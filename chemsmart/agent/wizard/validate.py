"""Validation helpers for rendered wizard server YAML."""

from __future__ import annotations

from dataclasses import dataclass

import yaml

from chemsmart.settings.server import Server

_REQUIRED_SERVER_KEYS = [
    "SCHEDULER",
    "QUEUE_NAME",
    "NUM_HOURS",
    "MEM_GB",
    "NUM_CORES",
    "SUBMIT_COMMAND",
]


@dataclass(frozen=True)
class ValidationResult:
    ok: bool
    errors: list[str]
    parsed: dict | None


def validate_server_yaml(
    yaml_text: str,
    server_name: str = "__wizard_candidate__",
) -> ValidationResult:
    """Load and validate a rendered server YAML candidate."""

    try:
        parsed = yaml.safe_load(yaml_text)
    except yaml.YAMLError as exc:
        return ValidationResult(
            ok=False,
            errors=[f"Invalid YAML: {exc}"],
            parsed=None,
        )

    if not isinstance(parsed, dict):
        return ValidationResult(
            ok=False,
            errors=["YAML root must be a mapping."],
            parsed=None,
        )

    errors: list[str] = []
    server_block = parsed.get("SERVER")
    if not isinstance(server_block, dict):
        errors.append("Missing SERVER mapping.")
    else:
        for key in _REQUIRED_SERVER_KEYS:
            if server_block.get(key) is None:
                errors.append(f"Missing required SERVER.{key}.")

    program_blocks = {
        key: value
        for key, value in parsed.items()
        if key != "SERVER" and isinstance(value, dict)
    }
    if not program_blocks:
        errors.append("At least one program block is required.")

    if isinstance(server_block, dict):
        try:
            Server(name=server_name, **server_block)
        except (TypeError, KeyError) as exc:
            errors.append(f"SERVER block could not initialize Server: {exc}")

    return ValidationResult(ok=not errors, errors=errors, parsed=parsed)
