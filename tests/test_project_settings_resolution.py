from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.user import ChemsmartUserSettings


def _write_project(path: Path, functional: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "gas:\n"
        f"  functional: {functional}\n"
        "  basis: def2svp\n"
        "solv:\n"
        f"  functional: {functional}\n"
        "  basis: def2svp\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("program", "settings_cls"),
    [
        ("gaussian", GaussianProjectSettings),
        ("orca", ORCAProjectSettings),
    ],
)
def test_workspace_project_precedes_legacy_cli_project(
    monkeypatch, tmp_path, program, settings_cls
):
    home = tmp_path / "home"
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    monkeypatch.chdir(workspace)
    monkeypatch.setattr(
        ChemsmartUserSettings,
        "USER_CONFIG_DIR",
        str(home / ".chemsmart"),
    )
    legacy = home / ".chemsmart" / program / "shared.yaml"
    _write_project(legacy, "pbe0")

    legacy_settings = settings_cls.from_project("shared")
    assert legacy_settings.opt_settings().functional == "pbe0"

    local = workspace / ".chemsmart" / program / "shared.yaml"
    _write_project(local, "b3lyp")

    workspace_settings = settings_cls.from_project("shared")
    assert workspace_settings.opt_settings().functional == "b3lyp"
