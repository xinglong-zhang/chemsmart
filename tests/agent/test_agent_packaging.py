"""Distribution boundaries for the optional ChemSmart agent runtimes."""

from __future__ import annotations

from pathlib import Path

import tomlkit


REPO_ROOT = Path(__file__).parents[2]


def _project_metadata() -> dict:
    return tomlkit.parse(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )


def test_optional_agent_dependencies_are_not_core_requirements():
    project = _project_metadata()["project"]
    core = "\n".join(project["dependencies"]).lower()
    extras = project["optional-dependencies"]

    for dependency in (
        "anthropic",
        "openai",
        "python-dotenv",
        "textual",
        "transformers",
        "mlx-lm",
        "cclib",
        "fastapi",
    ):
        assert dependency not in core

    assert set(
        (
            "agent",
            "agent-tui",
            "local-transformers",
            "local-mlx",
            "result-parsers",
            "local-server",
        )
    ).issubset(extras)


def test_distribution_includes_offline_basis_catalog():
    package_data = _project_metadata()["tool"]["setuptools"]["package-data"]

    assert "agent/harness/basis_sets/*.json" in package_data["chemsmart"]


def test_distribution_discovers_only_chemsmart_packages():
    package_find = _project_metadata()["tool"]["setuptools"]["packages"][
        "find"
    ]

    assert package_find["include"] == ["chemsmart*"]
    assert {"tests*", "scripts*", "docs*"}.issubset(package_find["exclude"])
