"""Conservative backstops for recurring compact-SPEC kind confusions."""

from __future__ import annotations

import re
from typing import Any

_WIBERG = re.compile(r"\bwiberg\b|bond[\s-]*index|bond[\s-]*order", re.I)
_INTERNAL = re.compile(
    r"modredundant"
    r"|\b(freez\w*|constrain\w*|fix\w*|hold\w*|keep\w*|lock\w*)\b"
    r"[^.]*\b(bond|distance|angle|dihedral|length)\b"
    r"|\b(bond|distance|angle|dihedral|length)\b"
    r"[^.]*\b(fixed|frozen|constrained|held|kept)\b",
    re.I,
)
_PAIR = re.compile(r"(\d+)\s*[-–to&,]+\s*(\d+)")
_ORCA_TS = re.compile(
    r"\borca\b(?=[^.]{0,120}\b("
    r"opt[-\s]*ts|transition[-\s]*state|first[-\s]*order\s+saddle|"
    r"saddle\s+(?:point\s+)?(?:search|optimization|opt|hunt)"
    r")\b)"
    r"|\b("
    r"opt[-\s]*ts|transition[-\s]*state|first[-\s]*order\s+saddle|"
    r"saddle\s+(?:point\s+)?(?:search|optimization|opt|hunt)"
    r")\b(?=[^.]{0,120}\borca\b)",
    re.I,
)


def _atoms_from_freeze(value: Any) -> list[int]:
    if isinstance(value, str):
        return [
            int(item)
            for item in re.split(r"[,\s]+", value.strip())
            if item.lstrip("-").isdigit()
        ]
    if isinstance(value, list):
        return [int(item) for item in value if isinstance(item, int)]
    return []


def disambiguate(query: str, spec: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    """Mutate and return ``spec`` when a high-confidence rule fires."""
    if not isinstance(spec, dict) or spec.get("intent") != "workflow":
        return spec, False
    changed = False
    for job in spec.get("jobs", []):
        if not isinstance(job, dict):
            continue
        kind = str(job.get("kind") or "")
        settings = job.get("settings") or {}
        program = kind.split(".", 1)[0] if "." in kind else "gaussian"

        if kind.endswith(".dias") and (
            _WIBERG.search(query or "") or "fragment_indices" not in settings
        ):
            job["kind"] = f"{program}.wbi"
            job.pop("settings", None)
            changed = True
            continue
        if kind.endswith(".wbi") and _WIBERG.search(query or ""):
            job.pop("settings", None)

        if program == "orca" and not kind.endswith(".ts") and _ORCA_TS.search(
            query or ""
        ):
            job["kind"] = "orca.ts"
            allowed = {"recalc_hess", "trust_radius", "tssearch_type"}
            job["settings"] = {
                key: value
                for key, value in settings.items()
                if key in allowed
            }
            if not job["settings"]:
                job.pop("settings", None)
            changed = True
            continue

        if (
            kind.endswith(".opt")
            and settings.get("freeze_atoms")
            and _INTERNAL.search(query or "")
        ):
            atoms = _atoms_from_freeze(settings["freeze_atoms"])
            if len(atoms) >= 2:
                job["kind"] = f"{program}.modred"
                job["settings"] = {"modred": [atoms]}
                changed = True
                continue

        if kind.endswith(".modred") and not settings.get("modred"):
            pairs = [
                [int(left), int(right)]
                for left, right in _PAIR.findall(query or "")
            ]
            if pairs:
                updated = dict(settings)
                updated["modred"] = pairs
                job["settings"] = updated
                changed = True
    return spec, changed


__all__ = ["disambiguate"]
