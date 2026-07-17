"""Value normalization shared by project-YAML services."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

KNOWN_PROGRAMS = {"gaussian", "orca"}
_DEF2_BASIS_PATTERN = (
    r"def2[- ]?(?:svp|svpd|tzvp|tzvpd|tzvpp|tzvppd|qzvp|qzvpd)"
)
_SOLVENT_ALIASES = {
    "acetonitrile": "acetonitrile",
    "chloroform": "chloroform",
    "cyclohexane": "cyclohexane",
    "dichloroethane": "dichloroethane",
    "dichloromethane": "dichloromethane",
    "dmso": "dmso",
    "ethanol": "ethanol",
    "methanol": "methanol",
    "thf": "thf",
    "tetrahydrofuran": "thf",
    "toluene": "toluene",
    "water": "water",
}


def normalize_program(program: str) -> str:
    normalized = (program or "").strip().lower()
    if normalized not in KNOWN_PROGRAMS:
        raise ValueError(f"unsupported project YAML program: {program!r}")
    return normalized


def normalize_project_name(project_name: str) -> str:
    stem = Path(project_name or "project").stem.lower()
    stem = re.sub(r"[^a-z0-9_.-]+", "_", stem).strip("._-")
    return stem or "project"


def normalize_basis_name(value: str) -> str:
    return value.lower().replace("-", "").replace(" ", "")


def normalize_basis_if_known(value: str | None) -> str | None:
    if value is None:
        return None
    if re.fullmatch(rf"(?i){_DEF2_BASIS_PATTERN}", value.strip()):
        return normalize_basis_name(value)
    return value.strip()


def normalize_solvent_id(value: str) -> str:
    key = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    return _SOLVENT_ALIASES.get(key, key.replace(" ", ""))


def normalize_yaml_text(yaml_text: str) -> str:
    return yaml_text.rstrip() + "\n"


def string_or_none(value: Any) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    return normalized or None


def string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


def first_or_none(values: Any) -> Any:
    for value in values:
        return value
    return None


__all__ = [
    "KNOWN_PROGRAMS",
    "first_or_none",
    "normalize_basis_if_known",
    "normalize_basis_name",
    "normalize_program",
    "normalize_project_name",
    "normalize_solvent_id",
    "normalize_yaml_text",
    "string_list",
    "string_or_none",
]
