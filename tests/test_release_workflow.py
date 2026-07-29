from pathlib import Path


def test_package_release_accepts_only_version_tags() -> None:
    workflow = (
        Path(__file__).resolve().parents[1]
        / ".github"
        / "workflows"
        / "release.yml"
    ).read_text(encoding="utf-8")

    assert "- 'v[0-9]*'" in workflow
    assert "- '*'" not in workflow
    assert "agent-hpc-rc1" not in workflow
