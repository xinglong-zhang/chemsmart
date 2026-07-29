"""list_workspace visibility and save_geometry handoff tools."""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.tools import build_molecule
from chemsmart.agent.tools_fs import list_workspace, save_geometry

_WATER_XYZ = """3
water
O 0.0 0.0 0.119
H 0.0 0.763 -0.477
H 0.0 -0.763 -0.477
"""


def _seed_workspace(root: Path) -> None:
    (root / "geoms").mkdir()
    (root / "geoms" / "water.xyz").write_text(_WATER_XYZ)
    (root / "ts_guess.com").write_text("#p opt\n")
    (root / "old_run.out").write_text("* finished run\n")
    (root / "notes.txt").write_text("misc\n")
    (root / ".chemsmart" / "xtb").mkdir(parents=True)
    (root / ".chemsmart" / "xtb" / "fast.yaml").write_text("opt: {}\n")
    (root / ".hidden").mkdir()
    (root / ".hidden" / "secret.xyz").write_text(_WATER_XYZ)


class TestListWorkspace:
    def test_classifies_geometries_projects_and_outputs(
        self, monkeypatch, tmp_path
    ) -> None:
        _seed_workspace(tmp_path)
        monkeypatch.chdir(tmp_path)

        result = list_workspace()

        assert sorted(result["geometry_files"]) == [
            "geoms/water.xyz",
            "ts_guess.com",
        ]
        assert result["project_yamls"] == [
            {"program": "xtb", "project": "fast"}
        ]
        assert result["outputs"] == [{"path": "old_run.out"}]
        assert result["other_file_count"] == 1
        assert result["truncated"] is False
        # Hidden trees other than .chemsmart stay invisible.
        assert not any("secret" in path for path in result["geometry_files"])

    def test_rejects_paths_outside_cwd(self, monkeypatch, tmp_path) -> None:
        monkeypatch.chdir(tmp_path)
        assert list_workspace(subdir="../..")["error"] == "path_outside_cwd"
        assert (
            list_workspace(subdir="missing")["error"] == "directory_not_found"
        )

    def test_caps_entries_and_reports_truncation(
        self, monkeypatch, tmp_path
    ) -> None:
        monkeypatch.chdir(tmp_path)
        for index in range(30):
            (tmp_path / f"mol_{index:03d}.xyz").write_text(_WATER_XYZ)

        result = list_workspace(max_entries=10)

        assert len(result["geometry_files"]) == 10
        assert result["truncated"] is True


class TestSaveGeometry:
    def test_saves_and_grounds_the_molecule(
        self, monkeypatch, tmp_path
    ) -> None:
        monkeypatch.chdir(tmp_path)
        (tmp_path / "water.xyz").write_text(_WATER_XYZ)
        molecule = build_molecule("water.xyz")

        result = save_geometry(molecule, "refined/water_refined.xyz")

        saved = Path(result["path"])
        assert saved == tmp_path / "refined" / "water_refined.xyz"
        assert result["atom_count"] == 3
        assert saved.read_text().strip().startswith("3")
        # Grounding now points at the saved file so a follow-up job's
        # `-f` renders with a real path.
        assert getattr(molecule, "_agent_source_filepath") == str(saved)

    def test_refuses_overwrite_bad_suffix_and_escape(
        self, monkeypatch, tmp_path
    ) -> None:
        monkeypatch.chdir(tmp_path)
        (tmp_path / "water.xyz").write_text(_WATER_XYZ)
        molecule = build_molecule("water.xyz")

        assert (
            save_geometry(molecule, "../escape.xyz")["error"]
            == "path_outside_cwd"
        )
        assert (
            save_geometry(molecule, "geom.com")["error"]
            == "unsupported_suffix"
        )
        first = save_geometry(molecule, "once.xyz")
        assert "error" not in first
        assert save_geometry(molecule, "once.xyz")["error"] == "file_exists"
        replaced = save_geometry(molecule, "once.xyz", overwrite=True)
        assert "error" not in replaced

    def test_rejects_non_molecule_payloads(
        self, monkeypatch, tmp_path
    ) -> None:
        monkeypatch.chdir(tmp_path)
        assert (
            save_geometry({"not": "a molecule"}, "x.xyz")["error"]
            == "not_a_molecule"
        )
