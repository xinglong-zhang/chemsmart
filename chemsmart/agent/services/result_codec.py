"""Serialize agent tool results and resolve plan references."""

from __future__ import annotations

import importlib
import re
from functools import lru_cache
from typing import Any

from pydantic import BaseModel

from chemsmart.agent.handles import HandleStore, is_handle_id
from chemsmart.agent.serialization import generic_json_safe

REFERENCE_RE = re.compile(
    r"^\$step(?P<index>\d+)(?P<path>(?:\.[A-Za-z_][A-Za-z0-9_]*)*)$"
)


def resolve_refs(
    value: Any,
    prior_results: list[Any],
    handle_store: HandleStore | None = None,
) -> Any:
    """Resolve nested step references and persisted handle ids."""
    if isinstance(value, str):
        return _resolve_ref_string(
            value,
            prior_results,
            handle_store=handle_store,
        )
    if isinstance(value, list):
        return [
            resolve_refs(item, prior_results, handle_store=handle_store)
            for item in value
        ]
    if isinstance(value, dict):
        return {
            key: resolve_refs(item, prior_results, handle_store=handle_store)
            for key, item in value.items()
        }
    return value


def json_safe(value: Any) -> Any:
    """Project runtime chemistry objects into stable JSON-safe values."""
    if _is_instance_of(
        value,
        module_name="chemsmart.io.molecules.structure",
        class_names=("Molecule",),
    ):
        return _serialize_molecule(value)
    if _is_settings(value):
        return _serialize_settings(value)
    if _is_job(value):
        return _serialize_job(value)
    if isinstance(value, BaseModel):
        return json_safe(value.model_dump())
    return generic_json_safe(value, recurse=json_safe)


def preview_value(value: Any) -> Any:
    return json_safe(value)


def restore_json_result(value: Any) -> Any:
    """Restore serialized molecules, settings, and jobs recursively."""
    if isinstance(value, list):
        return [restore_json_result(item) for item in value]
    if not isinstance(value, dict):
        return value

    marker = value.get("__chemsmart_type__")
    if marker == "molecule":
        return _restore_molecule(value)
    if marker == "settings":
        return _restore_settings(value)
    if marker == "job":
        return _restore_job(value)
    return {key: restore_json_result(item) for key, item in value.items()}


def _resolve_ref_string(
    value: str,
    prior_results: list[Any],
    handle_store: HandleStore | None = None,
) -> Any:
    match = REFERENCE_RE.match(value)
    if match is not None:
        index = int(match.group("index")) - 1
        if index < 0 or index >= len(prior_results):
            raise IndexError(f"Step reference {value!r} is out of range")
        resolved = prior_results[index]
        for part in match.group("path").split("."):
            if not part:
                continue
            resolved = (
                resolved[part]
                if isinstance(resolved, dict)
                else getattr(resolved, part)
            )
        return restore_json_result(resolved)

    if handle_store is not None and is_handle_id(value):
        try:
            return handle_store.get(value)
        except KeyError:
            return restore_json_result(handle_store.get_summary(value))
    return value


def _serialize_molecule(value: Any) -> dict[str, Any]:
    positions = value.positions
    if hasattr(positions, "tolist"):
        positions = positions.tolist()
    return {
        "__chemsmart_type__": "molecule",
        "symbols": json_safe(list(value.symbols)),
        "positions": json_safe(positions),
        "charge": value.charge,
        "multiplicity": value.multiplicity,
        "frozen_atoms": json_safe(value.frozen_atoms),
        "pbc_conditions": json_safe(value.pbc_conditions),
        "translation_vectors": json_safe(value.translation_vectors),
        "energy": value.energy,
        "forces": json_safe(value.forces),
        "velocities": json_safe(value.velocities),
        "info": json_safe(value.info),
    }


def _serialize_settings(value: Any) -> dict[str, Any]:
    return {
        "__chemsmart_type__": "settings",
        "module": value.__class__.__module__,
        "class": value.__class__.__name__,
        **{
            key: json_safe(item)
            for key, item in vars(value).items()
            if not key.startswith("_")
        },
    }


def _serialize_job(value: Any) -> dict[str, Any]:
    return {
        "__chemsmart_type__": "job",
        "module": value.__class__.__module__,
        "class": value.__class__.__name__,
        "molecule": json_safe(value.molecule),
        "settings": json_safe(value.settings),
        "label": value.label,
        "folder": value.folder,
    }


def _restore_molecule(value: dict[str, Any]) -> Any:
    return _molecule_class()(
        symbols=value.get("symbols"),
        positions=value.get("positions"),
        charge=value.get("charge"),
        multiplicity=value.get("multiplicity"),
        frozen_atoms=value.get("frozen_atoms"),
        pbc_conditions=value.get("pbc_conditions"),
        translation_vectors=value.get("translation_vectors"),
        energy=value.get("energy"),
        forces=value.get("forces"),
        velocities=value.get("velocities"),
        info=value.get("info"),
    )


def _restore_settings(value: dict[str, Any]) -> Any:
    settings_cls = _load_class(value["module"], value["class"])
    kwargs = {
        key: restore_json_result(item)
        for key, item in value.items()
        if key not in {"__chemsmart_type__", "module", "class"}
    }
    return settings_cls(**kwargs)


def _restore_job(value: dict[str, Any]) -> Any:
    job_cls = _load_class(value["module"], value["class"])
    job = job_cls(
        molecule=restore_json_result(value["molecule"]),
        settings=restore_json_result(value["settings"]),
        label=value["label"],
        jobrunner=None,
    )
    if value.get("folder"):
        job.set_folder(value["folder"])
    return job


def _is_settings(value: Any) -> bool:
    return _is_instance_of(
        value,
        module_name="chemsmart.jobs.gaussian.settings",
        class_names=("GaussianJobSettings",),
    ) or _is_instance_of(
        value,
        module_name="chemsmart.jobs.orca.settings",
        class_names=("ORCAJobSettings",),
    )


def _is_job(value: Any) -> bool:
    return (
        value.__class__.__module__.startswith("chemsmart.jobs.")
        and hasattr(value, "molecule")
        and hasattr(value, "settings")
    )


def _is_instance_of(
    value: Any,
    *,
    module_name: str,
    class_names: tuple[str, ...],
) -> bool:
    value_type = type(value)
    return (
        value_type.__module__ == module_name
        and value_type.__name__ in class_names
    )


@lru_cache(maxsize=1)
def _molecule_class():
    from chemsmart.io.molecules.structure import Molecule

    return Molecule


def _load_class(module_name: str, class_name: str) -> type[Any]:
    module = importlib.import_module(module_name)
    return getattr(module, class_name)


__all__ = [
    "REFERENCE_RE",
    "json_safe",
    "preview_value",
    "resolve_refs",
    "restore_json_result",
]
