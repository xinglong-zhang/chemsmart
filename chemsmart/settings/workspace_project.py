"""Workspace-local project YAML discovery."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PROJECT_PROGRAMS = ("gaussian", "orca")


@dataclass(frozen=True)
class WorkspaceProjectStatus:
    loaded: bool
    project: str
    program: str
    path: Path | None
    candidates: tuple[Path, ...]
    message: str


def workspace_project_dir(
    program: str,
    *,
    cwd: str | Path | None = None,
) -> Path:
    return Path(cwd or Path.cwd()).resolve() / ".chemsmart" / program


def workspace_project_path(
    project: str,
    program: str,
    *,
    cwd: str | Path | None = None,
) -> Path:
    name = _normalize_project_name(project)
    return workspace_project_dir(program, cwd=cwd) / f"{name}.yaml"


def iter_workspace_project_yaml(
    *,
    cwd: str | Path | None = None,
) -> tuple[Path, ...]:
    root = Path(cwd or Path.cwd()).resolve() / ".chemsmart"
    paths: list[Path] = []
    for program in PROJECT_PROGRAMS:
        folder = root / program
        if not folder.is_dir():
            continue
        paths.extend(
            path
            for path in folder.iterdir()
            if path.is_file() and path.suffix.lower() in {".yaml", ".yml"}
        )
    return tuple(sorted(paths))


def workspace_project_names(
    program: str,
    *,
    cwd: str | Path | None = None,
) -> list[str]:
    folder = workspace_project_dir(program, cwd=cwd)
    if not folder.is_dir():
        return []
    return sorted(
        path.stem
        for path in folder.iterdir()
        if path.is_file() and path.suffix.lower() in {".yaml", ".yml"}
    )


def resolve_workspace_project(
    *,
    cwd: str | Path | None = None,
) -> WorkspaceProjectStatus:
    candidates = iter_workspace_project_yaml(cwd=cwd)
    if not candidates:
        return WorkspaceProjectStatus(
            loaded=False,
            project="",
            program="",
            path=None,
            candidates=(),
            message=(
                "No workspace project YAML found. Create one with /init, then "
                "write it with /write-project."
            ),
        )
    if len(candidates) > 1:
        names = ", ".join(
            f"{path.parent.name}:{path.stem}" for path in candidates[:4]
        )
        if len(candidates) > 4:
            names += ", ..."
        return WorkspaceProjectStatus(
            loaded=False,
            project="",
            program="",
            path=None,
            candidates=candidates,
            message=(
                "Multiple workspace project YAML files found; specify one "
                f"with -p/--project. Candidates: {names}"
            ),
        )
    path = candidates[0]
    return WorkspaceProjectStatus(
        loaded=True,
        project=path.stem,
        program=path.parent.name,
        path=path,
        candidates=candidates,
        message=f"Loaded workspace project YAML: {path}",
    )


def _normalize_project_name(project: str) -> str:
    name = str(project or "").strip()
    if name.endswith((".yaml", ".yml")):
        name = Path(name).stem
    return name
