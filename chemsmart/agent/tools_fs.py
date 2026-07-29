from __future__ import annotations

import os
from pathlib import Path
from typing import Any

_BINARY_SNIFF_BYTES = 8192

_GEOMETRY_SUFFIXES = (".xyz", ".com", ".gjf", ".inp", ".sdf", ".pdb")
_OUTPUT_SUFFIXES = (".out", ".log")
_LIST_SKIP_DIRS = frozenset(
    {"__pycache__", ".git", ".venv", "node_modules", "scratch"}
)
_LIST_MAX_DEPTH = 3


def read(
    path: str,
    start_line: int = 1,
    limit: int = 200,
) -> dict[str, Any]:
    """Read text from a local file using 1-based numbered lines."""

    cwd = Path.cwd().resolve(strict=False)
    resolved_path = Path(path).resolve(strict=False)

    if not resolved_path.is_relative_to(cwd):
        return {
            "error": "path_outside_cwd",
            "path": str(resolved_path),
            "cwd": str(cwd),
        }

    if not resolved_path.exists() or not resolved_path.is_file():
        return {
            "error": "file_not_found",
            "path": str(resolved_path),
        }

    with resolved_path.open("rb") as handle:
        if b"\x00" in handle.read(_BINARY_SNIFF_BYTES):
            return {
                "error": "binary_file_refused",
                "path": str(resolved_path),
            }

    start_index = max(start_line - 1, 0)
    line_limit = max(limit, 0)
    lines = resolved_path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines()
    total_lines = len(lines)
    selected_lines = lines[start_index : start_index + line_limit]
    if selected_lines:
        end_line = start_line + len(selected_lines) - 1
    else:
        end_line = start_line - 1
    truncated = start_index + len(selected_lines) < total_lines
    content = "\n".join(
        f"{line_number:>6}\t{line}"
        for line_number, line in enumerate(
            selected_lines,
            start=start_index + 1,
        )
    )

    return {
        "path": str(resolved_path),
        "start_line": start_line,
        "end_line": end_line,
        "total_lines": total_lines,
        "truncated": truncated,
        "content": content,
    }


def list_workspace(
    subdir: str = ".",
    max_entries: int = 200,
) -> dict[str, Any]:
    """Summarize workspace files: geometries, project YAMLs, job outputs."""

    from chemsmart.settings.workspace_project import (
        iter_workspace_project_yaml,
    )

    cwd = Path.cwd().resolve(strict=False)
    root = (cwd / subdir).resolve(strict=False)
    if not root.is_relative_to(cwd):
        return {
            "error": "path_outside_cwd",
            "path": str(root),
            "cwd": str(cwd),
        }
    if not root.is_dir():
        return {"error": "directory_not_found", "path": str(root)}

    entry_cap = max(int(max_entries), 1)
    geometry_files: list[str] = []
    outputs: list[dict[str, Any]] = []
    other_file_count = 0
    truncated = False
    for current, dirnames, filenames in os.walk(root):
        current_path = Path(current)
        depth = len(current_path.relative_to(root).parts)
        # .chemsmart project YAMLs are reported separately below; every
        # other hidden tree is runtime noise for a chemistry workspace.
        dirnames[:] = sorted(
            name
            for name in dirnames
            if name not in _LIST_SKIP_DIRS
            and (not name.startswith(".") or name == ".chemsmart")
            and depth < _LIST_MAX_DEPTH
        )
        if current_path.name == ".chemsmart" or ".chemsmart" in (
            current_path.relative_to(root).parts
        ):
            continue
        for filename in sorted(filenames):
            if filename.startswith("."):
                continue
            if len(geometry_files) + len(outputs) >= entry_cap:
                truncated = True
                break
            path = current_path / filename
            relative = str(path.relative_to(cwd))
            suffix = path.suffix.lower()
            if suffix in _GEOMETRY_SUFFIXES:
                geometry_files.append(relative)
            elif suffix in _OUTPUT_SUFFIXES:
                try:
                    mtime = path.stat().st_mtime
                except OSError:
                    mtime = 0.0
                outputs.append({"path": relative, "mtime": mtime})
            else:
                other_file_count += 1
        if truncated:
            break

    outputs.sort(key=lambda item: item["mtime"], reverse=True)
    for item in outputs:
        item.pop("mtime", None)

    project_yamls = [
        {"program": path.parent.name, "project": path.stem}
        for path in iter_workspace_project_yaml(cwd=cwd)
    ]
    return {
        "root": str(root),
        "geometry_files": geometry_files,
        "project_yamls": project_yamls,
        "outputs": outputs,
        "other_file_count": other_file_count,
        "truncated": truncated,
        "workspace_rules_file": (cwd / "CHEMSMART.md").is_file(),
        "user_rules_file": (
            Path.home() / ".chemsmart" / "CHEMSMART.md"
        ).is_file(),
    }


def save_geometry(
    molecule: Any,
    filename: str,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write a molecule's coordinates to a workspace .xyz file.

    Persists an in-memory geometry — typically the result of
    ``extract_optimized_geometry`` after a cheap xTB pre-optimization — so a
    follow-up Gaussian/ORCA job can consume it as its ``-f`` input.
    """

    cwd = Path.cwd().resolve(strict=False)
    resolved_path = (cwd / filename).resolve(strict=False)
    if not resolved_path.is_relative_to(cwd):
        return {
            "error": "path_outside_cwd",
            "path": str(resolved_path),
            "cwd": str(cwd),
        }
    if resolved_path.suffix.lower() != ".xyz":
        return {
            "error": "unsupported_suffix",
            "path": str(resolved_path),
            "message": "save_geometry writes .xyz files only",
        }
    if resolved_path.exists() and not overwrite:
        return {
            "error": "file_exists",
            "path": str(resolved_path),
            "message": "pass overwrite=true to replace the existing file",
        }
    symbols = getattr(molecule, "symbols", None)
    positions = getattr(molecule, "positions", None)
    if symbols is None or positions is None:
        return {
            "error": "not_a_molecule",
            "message": "molecule must be a build_molecule or "
            "extract_optimized_geometry result handle",
        }
    resolved_path.parent.mkdir(parents=True, exist_ok=True)
    molecule.write(str(resolved_path), format="xyz", mode="w")
    # Ground the saved file as this molecule's source so any later
    # `chemsmart run/sub ... -f` command renders with a real path.
    setattr(molecule, "_agent_source_filepath", str(resolved_path))
    setattr(molecule, "_agent_source_index", "-1")
    return {
        "path": str(resolved_path),
        "atom_count": len(list(symbols)),
        "charge": getattr(molecule, "charge", None),
        "multiplicity": getattr(molecule, "multiplicity", None),
    }
