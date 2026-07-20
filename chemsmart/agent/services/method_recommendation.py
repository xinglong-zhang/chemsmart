"""Conservative project-based method recommendation service."""

from __future__ import annotations

import os
from typing import Any

import yaml

from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.periodictable import PeriodicTable

_TASK_PROJECT_MAP = {
    "opt": ["dft_default", "organics"],
    "opt+freq": ["dft_default", "organics"],
    "ts": ["tspaths"],
}
_TASK_JOBTYPE_MAP = {
    "opt": "opt",
    "opt+freq": "opt",
    "ts": "ts",
}
_PERIODIC_TABLE = PeriodicTable()


def recommend_method(
    task: str,
    charge: int = 0,
    multiplicity: int = 1,
    atomic_numbers: list[int] | None = None,
    project_hint: str | None = None,
) -> dict[str, Any]:
    """Return a conservative project-based method recommendation or no-match."""

    user_settings = ChemsmartUserSettings()
    available_project_paths = _get_available_project_paths(user_settings)
    available_projects = sorted(available_project_paths)
    normalized_task = (task or "").strip().lower()
    heavy_symbols = _get_heavy_symbols(atomic_numbers)

    if project_hint in available_project_paths:
        return _matched_recommendation(
            project_name=project_hint,
            project_path=available_project_paths[project_hint],
            task=normalized_task,
            rationale=(f"matched project_hint rule: {project_hint}.yaml"),
            available_projects=available_projects,
        )
    if charge != 0:
        return _no_match_response(available_projects, f"charge={charge}")
    if multiplicity > 1:
        return _no_match_response(
            available_projects,
            f"multiplicity={multiplicity}",
        )

    heavy_element_projects = _get_heavy_element_projects(
        available_project_paths=available_project_paths,
        task=normalized_task,
        heavy_symbols=heavy_symbols,
    )
    if heavy_symbols and not heavy_element_projects:
        return _no_match_response(
            available_projects,
            f"no project lists heavy_elements for {', '.join(heavy_symbols)}",
        )

    task_candidates = [
        name
        for name in _TASK_PROJECT_MAP.get(normalized_task, [])
        if name in available_project_paths
    ]
    if heavy_symbols:
        task_candidates = [
            name for name in task_candidates if name in heavy_element_projects
        ]
    if task_candidates:
        name = task_candidates[0]
        return _matched_recommendation(
            project_name=name,
            project_path=available_project_paths[name],
            task=normalized_task,
            rationale=(
                f"matched task rule for '{normalized_task}': {name}.yaml"
            ),
            available_projects=available_projects,
        )
    if heavy_element_projects:
        name = heavy_element_projects[0]
        return _matched_recommendation(
            project_name=name,
            project_path=available_project_paths[name],
            task=normalized_task,
            rationale=(
                "matched heavy_elements rule for "
                f"{', '.join(heavy_symbols)}: {name}.yaml"
            ),
            available_projects=available_projects,
        )
    return _no_match_response(
        available_projects,
        f"no task rule matched for '{normalized_task}'",
    )


def _get_available_project_paths(
    user_settings: ChemsmartUserSettings,
) -> dict[str, str]:
    paths = {}
    for filepath in user_settings.gaussian_project_yaml_files:
        name = os.path.basename(filepath).removesuffix(".yaml")
        if name != "defaults":
            paths[name] = filepath
    return paths


def _get_heavy_symbols(atomic_numbers: list[int] | None) -> list[str]:
    if not atomic_numbers:
        return []
    symbols = {
        _PERIODIC_TABLE.to_symbol(number)
        for number in atomic_numbers
        if number >= 54
    }
    return _PERIODIC_TABLE.sorted_periodic_table_list(list(symbols))


def _get_heavy_element_projects(
    available_project_paths: dict[str, str],
    task: str,
    heavy_symbols: list[str],
) -> list[str]:
    if not heavy_symbols:
        return []
    matches = []
    for name, path in available_project_paths.items():
        settings = _get_project_settings(path, task)
        configured = set(
            _normalize_heavy_elements(settings.get("heavy_elements")) or []
        )
        if set(heavy_symbols).issubset(configured):
            matches.append(name)
    return sorted(matches)


def _matched_recommendation(
    project_name: str,
    project_path: str,
    task: str,
    rationale: str,
    available_projects: list[str],
) -> dict[str, Any]:
    settings = _get_project_settings(project_path, task)
    return {
        "match": project_name,
        "functional": settings.get("functional"),
        "basis": settings.get("basis"),
        "solvent_model": settings.get("solvent_model"),
        "solvent_id": settings.get("solvent_id"),
        "heavy_elements": _normalize_heavy_elements(
            settings.get("heavy_elements")
        ),
        "heavy_elements_basis": settings.get("heavy_elements_basis"),
        "rationale": rationale,
        "available_projects": available_projects,
    }


def _no_match_response(
    available_projects: list[str],
    reason: str,
) -> dict[str, Any]:
    return {
        "match": None,
        "functional": None,
        "basis": None,
        "solvent_model": None,
        "solvent_id": None,
        "heavy_elements": None,
        "heavy_elements_basis": None,
        "rationale": (
            f"no rule matched: {reason}; pick from available_projects"
        ),
        "available_projects": available_projects,
    }


def _get_project_settings(project_path: str, task: str) -> dict[str, Any]:
    project_settings = _load_yaml(project_path)
    gas_settings = project_settings.get("gas")
    solv_settings = project_settings.get("solv")
    jobtype = _TASK_JOBTYPE_MAP.get(task, "opt")
    if jobtype == "sp":
        settings = solv_settings or gas_settings or project_settings
    else:
        settings = gas_settings or solv_settings or project_settings
    return settings if isinstance(settings, dict) else {}


def _normalize_heavy_elements(
    heavy_elements: list[str] | str | None,
) -> list[str] | None:
    if heavy_elements is None:
        return None
    if isinstance(heavy_elements, str):
        heavy_elements = heavy_elements.replace(",", " ").split()
    return [_PERIODIC_TABLE.to_element(item) for item in heavy_elements]


def _load_yaml(filepath: str) -> dict[str, Any]:
    with open(filepath) as file:
        return yaml.safe_load(file) or {}


__all__ = ["recommend_method"]
