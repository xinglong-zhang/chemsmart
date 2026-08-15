from pathlib import Path

from click.testing import CliRunner

import chemsmart.cli.update as update_module
from chemsmart.cli.update import ProjectsUpdater, update


def _write(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def _invoke_update_projects(config_root: Path):
    return CliRunner().invoke(
        update,
        ["projects"],
        env={"CHEMSMART_CONFIG_DIR": str(config_root)},
    )


def test_update_projects_adds_missing_template_content_without_overwriting(
    tmp_path, monkeypatch
):
    assert ProjectsUpdater._template_config_dir().is_dir()

    template_root = tmp_path / "templates" / ".chemsmart"
    config_root = tmp_path / "user" / ".chemsmart"

    _write(template_root / "fake_a" / "existing.yaml", "template existing\n")
    _write(template_root / "fake_a" / "missing.yaml", "template missing\n")
    _write(
        template_root / "fake_a" / "nested" / "nested.yaml",
        "template nested\n",
    )
    _write(template_root / "fake_b" / "new.yaml", "template new\n")
    (template_root / "fake_empty").mkdir(parents=True)
    _write(template_root / "server" / "ignored.yaml", "server template\n")
    _write(template_root / "usersettings.yaml", "user settings template\n")

    existing = _write(
        config_root / "fake_a" / "existing.yaml", "user existing\n"
    )
    custom = _write(config_root / "fake_a" / "custom.yaml", "user custom\n")
    personal = _write(
        config_root / "custom_project" / "personal.yaml", "user personal\n"
    )

    monkeypatch.setattr(
        ProjectsUpdater,
        "_template_config_dir",
        staticmethod(lambda: template_root),
    )

    result = _invoke_update_projects(config_root)

    assert result.exit_code == 0, result.output
    assert "added missing content: fake_a" in result.output
    assert "added missing content: fake_b" in result.output
    assert "added missing content: fake_empty" in result.output
    assert existing.read_text(encoding="utf-8") == "user existing\n"
    assert custom.read_text(encoding="utf-8") == "user custom\n"
    assert personal.read_text(encoding="utf-8") == "user personal\n"
    assert (config_root / "fake_a" / "missing.yaml").read_text(
        encoding="utf-8"
    ) == "template missing\n"
    assert (config_root / "fake_a" / "nested" / "nested.yaml").read_text(
        encoding="utf-8"
    ) == "template nested\n"
    assert (config_root / "fake_b" / "new.yaml").read_text(
        encoding="utf-8"
    ) == "template new\n"
    assert (config_root / "fake_empty").is_dir()
    assert not (config_root / "server").exists()
    assert not (config_root / "usersettings.yaml").exists()

    files_before = {
        path.relative_to(config_root): path.read_bytes()
        for path in config_root.rglob("*")
        if path.is_file()
    }
    repeated = _invoke_update_projects(config_root)
    files_after = {
        path.relative_to(config_root): path.read_bytes()
        for path in config_root.rglob("*")
        if path.is_file()
    }

    assert repeated.exit_code == 0, repeated.output
    assert "already up to date" in repeated.output
    assert files_after == files_before


def test_update_projects_reports_errors_without_damaging_user_content(
    tmp_path, monkeypatch
):
    template_root = tmp_path / "templates" / ".chemsmart"
    config_root = tmp_path / "user" / ".chemsmart"
    _write(template_root / "fake_a" / "settings.yaml", "template value\n")
    conflicting_path = _write(config_root / "fake_a", "user value\n")

    monkeypatch.setattr(
        ProjectsUpdater,
        "_template_config_dir",
        staticmethod(lambda: template_root),
    )

    result = _invoke_update_projects(config_root)

    assert result.exit_code != 0
    assert "Expected a directory" in result.output
    assert "Traceback" not in result.output
    assert conflicting_path.read_text(encoding="utf-8") == "user value\n"

    copy_config_root = tmp_path / "copy-user" / ".chemsmart"
    custom = _write(copy_config_root / "custom.yaml", "user custom\n")

    def fail_during_copy(source_file, destination_file):
        destination_file.write(b"partial content")
        raise OSError("copy failed")

    monkeypatch.setattr(update_module.shutil, "copyfileobj", fail_during_copy)

    copy_result = _invoke_update_projects(copy_config_root)

    assert copy_result.exit_code != 0
    assert "copy failed" in copy_result.output
    assert not (copy_config_root / "fake_a" / "settings.yaml").exists()
    assert custom.read_text(encoding="utf-8") == "user custom\n"
