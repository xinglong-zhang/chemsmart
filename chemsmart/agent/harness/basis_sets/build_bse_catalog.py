"""Build the chemsmart harness basis-set catalog from Basis Set Exchange.

Run from the repository root:

    python -m chemsmart.agent.harness.basis_sets.build_bse_catalog

The generated JSON is intentionally compact. It stores BSE names, family/role
metadata, element coverage, auxiliaries, and whether the basis can be rendered
through the BSE Gaussian/ORCA writers.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
from typing import Any

import basis_set_exchange as bse

from chemsmart.agent.harness.basis_sets.catalog import normalize_basis_name

PROGRAM_FORMATS = {
    "gaussian": "gaussian94",
    "orca": "orca",
}


def build_catalog() -> dict[str, Any]:
    metadata_by_key = bse.get_metadata()
    basis_sets: dict[str, Any] = {}
    aliases: dict[str, str] = {}
    display_name_to_key: dict[str, str] = {}
    programs: dict[str, dict[str, Any]] = {
        program: {"format": fmt, "basis_names": []}
        for program, fmt in PROGRAM_FORMATS.items()
    }

    for key in sorted(metadata_by_key):
        item = metadata_by_key[key]
        latest = str(item.get("latest_version", "0"))
        version = item.get("versions", {}).get(latest, {})
        elements = sorted(int(z) for z in version.get("elements", []))
        display_name = item.get("display_name") or item.get("basename") or key
        entry_programs = [
            program
            for program, fmt in PROGRAM_FORMATS.items()
            if _renderable(display_name, elements, fmt)
        ]
        entry = {
            "display_name": display_name,
            "bse_key": key,
            "family": item.get("family"),
            "role": item.get("role"),
            "function_types": sorted(item.get("function_types", [])),
            "latest_version": latest,
            "elements": elements,
            "auxiliaries": item.get("auxiliaries", {}),
            "programs": entry_programs,
        }
        basis_sets[key] = entry
        display_name_to_key[display_name] = key
        for alias in [display_name, item.get("basename"), key, *item.get("other_names", [])]:
            if not alias:
                continue
            aliases.setdefault(normalize_basis_name(str(alias)), key)
        for program in entry_programs:
            programs[program]["basis_names"].append(display_name)

    for program_data in programs.values():
        program_data["basis_names"].sort()
        program_data["count"] = len(program_data["basis_names"])

    return {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "source": "Basis Set Exchange Python package",
            "source_package": "basis_set_exchange",
            "source_version": metadata.version("basis_set_exchange"),
            "basis_set_count": len(basis_sets),
            "schema_version": 1,
        },
        "programs": programs,
        "aliases": aliases,
        "display_name_to_key": display_name_to_key,
        "basis_sets": basis_sets,
    }


def write_catalog(path: Path | None = None) -> Path:
    target = path or Path(__file__).with_name("bse_basis_catalog.json")
    catalog = build_catalog()
    target.write_text(
        json.dumps(catalog, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return target


def _renderable(name: str, elements: list[int], fmt: str) -> bool:
    if not elements:
        return False
    try:
        bse.get_basis(name, elements=[elements[0]], fmt=fmt, header=False)
    except Exception:
        return False
    return True


if __name__ == "__main__":
    print(write_catalog())
