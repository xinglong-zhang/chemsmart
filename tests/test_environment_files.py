"""The Windows CI environment may differ from the canonical one only by
the documented ffmpeg exclusion.

``environment-windows.yml`` exists solely to dodge an upstream conda-forge
defect (gdk-pixbuf's win-64 post-link script). Without a guard the two
files drift silently and Windows CI stops testing what the other platforms
test, so this asserts the difference stays exactly one known line.
"""

from __future__ import annotations

from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[1]
_CANONICAL = _ROOT / "environment.yml"
_WINDOWS = _ROOT / "environment-windows.yml"

#: The only dependency the Windows environment is allowed to omit.
_WINDOWS_EXCLUDED_PREFIXES = ("ffmpeg",)


def _load(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _requirements(document: dict) -> tuple[list[str], list[str]]:
    """Split a dependency list into conda specs and pip specs."""
    conda: list[str] = []
    pip: list[str] = []
    for entry in document["dependencies"]:
        if isinstance(entry, dict):
            pip.extend(entry.get("pip", []))
        else:
            conda.append(str(entry))
    return conda, pip


def test_windows_environment_exists() -> None:
    assert (
        _WINDOWS.is_file()
    ), "the Windows CI job references environment-windows.yml"


def test_environment_name_and_channels_match() -> None:
    canonical = _load(_CANONICAL)
    windows = _load(_WINDOWS)
    assert windows["name"] == canonical["name"]
    assert windows["channels"] == canonical["channels"]


def test_windows_environment_omits_only_the_documented_dependency() -> None:
    canonical_conda, canonical_pip = _requirements(_load(_CANONICAL))
    windows_conda, windows_pip = _requirements(_load(_WINDOWS))

    assert (
        windows_pip == canonical_pip
    ), "pip requirements must be identical across environments"

    expected = [
        spec
        for spec in canonical_conda
        if not spec.startswith(_WINDOWS_EXCLUDED_PREFIXES)
    ]
    assert windows_conda == expected, (
        "environment-windows.yml may differ from environment.yml only by "
        f"omitting {_WINDOWS_EXCLUDED_PREFIXES}; fix the drift or update "
        "both files together"
    )


def test_the_exclusion_is_actually_present_upstream() -> None:
    """A stale guard would pass vacuously once ffmpeg leaves the canonical
    environment; require the excluded dependency to still exist there."""
    canonical_conda, _pip = _requirements(_load(_CANONICAL))
    assert any(
        spec.startswith(_WINDOWS_EXCLUDED_PREFIXES) for spec in canonical_conda
    ), (
        "environment.yml no longer pins ffmpeg; drop environment-windows.yml "
        "and point the Windows CI job back at environment.yml"
    )
