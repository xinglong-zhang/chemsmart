"""File indexing helpers for composer popups."""

from __future__ import annotations

from pathlib import Path

_BIASED_EXTENSIONS = [".xyz", ".log", ".com", ".inp", ".gjf", ".out"]


def iter_candidate_files(root: str | Path, *, depth: int = 3) -> list[Path]:
    root_path = Path(root)
    if not root_path.exists():
        return []

    candidates: list[Path] = []
    for path in root_path.rglob("*"):
        if not path.is_file():
            continue
        try:
            relative = path.relative_to(root_path)
        except ValueError:
            continue
        if len(relative.parts) > depth:
            continue
        candidates.append(path)

    def sort_key(path: Path):
        suffix = path.suffix.lower()
        bias = (
            _BIASED_EXTENSIONS.index(suffix)
            if suffix in _BIASED_EXTENSIONS
            else len(_BIASED_EXTENSIONS)
        )
        return (bias, len(path.parts), str(path).lower())

    return sorted(candidates, key=sort_key)
