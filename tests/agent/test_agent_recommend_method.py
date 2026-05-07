from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from chemsmart.agent.tools import recommend_method
from chemsmart.settings.user import ChemsmartUserSettings


@pytest.fixture()
def user_gaussian_settings_dir(monkeypatch, tmp_path: Path) -> Path:
    config_dir = tmp_path / ".chemsmart"
    gaussian_dir = config_dir / "gaussian"
    gaussian_dir.mkdir(parents=True)
    monkeypatch.setattr(
        ChemsmartUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    return gaussian_dir


def write_project(
    gaussian_dir: Path,
    name: str,
    project_data: dict,
) -> None:
    with open(gaussian_dir / f"{name}.yaml", "w") as file:
        yaml.safe_dump(project_data, file, sort_keys=False)


class TestRecommendMethod:
    def test_project_hint_matches_yaml(self, user_gaussian_settings_dir):
        write_project(
            user_gaussian_settings_dir,
            "hinted",
            {
                "gas": {
                    "functional": "wb97xd",
                    "basis": "def2svp",
                },
                "solv": {
                    "functional": "wb97xd",
                    "basis": "def2tzvp",
                    "solvent_model": "smd",
                    "solvent_id": "water",
                },
            },
        )

        recommendation = recommend_method(task="opt", project_hint="hinted")

        assert recommendation["match"] == "hinted"
        assert recommendation["functional"] == "wb97xd"
        assert recommendation["basis"] == "def2svp"
        assert recommendation["solvent_model"] is None
        assert recommendation["solvent_id"] is None
        assert "project_hint" in recommendation["rationale"]
        assert "hinted.yaml" in recommendation["rationale"]

    def test_unknown_task_returns_no_match(self, user_gaussian_settings_dir):
        write_project(
            user_gaussian_settings_dir,
            "alpha",
            {"gas": {"functional": "pbe0", "basis": "def2svp"}},
        )
        write_project(
            user_gaussian_settings_dir,
            "beta",
            {"gas": {"functional": "m062x", "basis": "def2tzvp"}},
        )

        recommendation = recommend_method(task="scan-ts")

        assert recommendation["match"] is None
        assert recommendation["available_projects"] == ["alpha", "beta"]
        assert (
            "no task rule matched for 'scan-ts'" in recommendation["rationale"]
        )

    def test_charge_returns_no_match(self, user_gaussian_settings_dir):
        write_project(
            user_gaussian_settings_dir,
            "dft_default",
            {"gas": {"functional": "pbe0", "basis": "def2svp"}},
        )

        recommendation = recommend_method(task="opt", charge=-1)

        assert recommendation["match"] is None
        assert "charge=-1" in recommendation["rationale"]

    def test_multiplicity_returns_no_match(self, user_gaussian_settings_dir):
        write_project(
            user_gaussian_settings_dir,
            "dft_default",
            {"gas": {"functional": "pbe0", "basis": "def2svp"}},
        )

        recommendation = recommend_method(task="opt", multiplicity=2)

        assert recommendation["match"] is None
        assert "multiplicity=2" in recommendation["rationale"]

    def test_heavy_elements_require_explicit_project(
        self,
        user_gaussian_settings_dir,
    ):
        write_project(
            user_gaussian_settings_dir,
            "dft_default",
            {"gas": {"functional": "pbe0", "basis": "def2svp"}},
        )

        no_match = recommend_method(task="opt", atomic_numbers=[79])

        assert no_match["match"] is None
        assert "Au" in no_match["rationale"]

        write_project(
            user_gaussian_settings_dir,
            "gold_ecp",
            {
                "gas": {
                    "functional": "pbe0",
                    "basis": "genecp",
                    "heavy_elements": ["Au"],
                    "heavy_elements_basis": "def2tzvp",
                }
            },
        )

        match = recommend_method(task="opt", atomic_numbers=[79])

        assert match["match"] == "gold_ecp"
        assert match["functional"] == "pbe0"
        assert match["basis"] == "genecp"
        assert match["heavy_elements"] == ["Au"]
        assert match["heavy_elements_basis"] == "def2tzvp"

    def test_functional_is_only_yaml_value_or_none(
        self,
        user_gaussian_settings_dir,
    ):
        write_project(
            user_gaussian_settings_dir,
            "hinted",
            {"gas": {"functional": "wb97xd", "basis": "def2svp"}},
        )
        write_project(
            user_gaussian_settings_dir,
            "gold_ecp",
            {
                "gas": {
                    "functional": "pbe0",
                    "basis": "genecp",
                    "heavy_elements": ["Au"],
                    "heavy_elements_basis": "def2tzvp",
                }
            },
        )

        allowed_functionals = {"wb97xd", "pbe0"}
        recommendations = [
            recommend_method(task="opt", project_hint="hinted"),
            recommend_method(task="unknown-task"),
            recommend_method(task="opt", atomic_numbers=[79]),
        ]

        for recommendation in recommendations:
            assert recommendation["functional"] in allowed_functionals | {None}
