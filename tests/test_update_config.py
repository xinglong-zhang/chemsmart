import builtins
import errno
from importlib import resources
from pathlib import Path
from textwrap import dedent

import pytest
import yaml
from click.testing import CliRunner

from chemsmart.cli.main import entry_point
from chemsmart.cli.update import update


def _server_dir(config_root: Path) -> Path:
    server_dir = config_root / "server"
    server_dir.mkdir(parents=True, exist_ok=True)
    return server_dir


def _write_server_yaml(config_root: Path, name: str, content: str) -> Path:
    path = _server_dir(config_root) / name
    path.write_text(dedent(content).lstrip(), encoding="utf-8")
    return path


def _invoke_update_config(config_root: Path, args=None, env_value=None):
    return CliRunner().invoke(
        update,
        ["config"] + (args or []),
        env={"CHEMSMART_CONFIG_DIR": str(env_value or config_root)},
    )


def _read_yaml(path: Path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _template_text(template_name="SLURM.yaml") -> str:
    template = (
        resources.files("chemsmart.settings")
        / "templates"
        / ".chemsmart"
        / "server"
        / template_name
    )
    return template.read_text(encoding="utf-8")


def _template_yaml(template_name="SLURM.yaml"):
    return yaml.safe_load(_template_text(template_name))


MINIMAL_SLURM = """
SERVER:
    SCHEDULER: SLURM
"""


OLD_SERVER_WITH_CUSTOM_GAUSSIAN = """
SERVER:
    SCHEDULER: SLURM
    QUEUE_NAME: custom

GAUSSIAN:
    EXEFOLDER: /custom/g16

CUSTOM_FIELD:
    VALUE: keep
"""


PRESERVED_SECTIONS_YAML = """
SERVER:
    SCHEDULER: SLURM
    QUEUE_NAME: custom

GAUSSIAN:
    EXEFOLDER: /custom/g16

XTB:
    EXEFOLDER: /old/xtb

CUSTOM_FIELD:
    VALUE: keep
"""


def test_update_config_registered_without_field_or_interactive(tmp_path):
    result = CliRunner().invoke(
        entry_point,
        ["update", "config", "--help"],
        env={"CHEMSMART_CONFIG_DIR": str(tmp_path)},
    )

    assert result.exit_code == 0
    assert "--server" in result.output
    assert "--field" not in result.output
    assert "--interactive" not in result.output


@pytest.mark.parametrize("args", [["-f", "XTB"], ["-i"]])
def test_update_config_removed_options_fail_without_writing(tmp_path, args):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path, args)

    assert result.exit_code != 0
    assert "No such option" in result.output
    assert slurm.read_text(encoding="utf-8") == before


def test_update_config_default_updates_existing_yaml_only(tmp_path):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_SLURM.replace("SCHEDULER: SLURM", "SCHEDULER: PBS"),
    )

    result = _invoke_update_config(tmp_path)
    server_dir = tmp_path / "server"

    assert result.exit_code == 0, result.output
    assert "SLURM.yaml:" in result.output
    assert "PBS.yaml:" in result.output
    assert "XTB" in _read_yaml(slurm)
    assert "XTB" in _read_yaml(pbs)
    assert not (server_dir / "local.yaml").exists()
    assert not (server_dir / "small.yaml").exists()
    assert not (server_dir / "my_cluster.yaml").exists()


@pytest.mark.parametrize("server_arg", ["SLURM", "SLURM.yaml"])
def test_update_config_server_option_accepts_name_or_yaml(
    tmp_path, server_arg
):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_SLURM.replace("SCHEDULER: SLURM", "SCHEDULER: PBS"),
    )

    result = _invoke_update_config(tmp_path, ["-s", server_arg])

    assert result.exit_code == 0, result.output
    assert "XTB" in _read_yaml(slurm)
    assert "XTB" not in _read_yaml(pbs)


def test_update_config_missing_server_fails_without_creating_file(tmp_path):
    _server_dir(tmp_path)

    result = _invoke_update_config(tmp_path, ["-s", "missing"])

    assert result.exit_code != 0
    assert "Server YAML not found" in result.output
    assert not (tmp_path / "server" / "missing.yaml").exists()


@pytest.mark.parametrize(
    "server_arg, expected",
    [
        ("../SLURM", "not a path"),
        ("/tmp/SLURM.yaml", "not a path"),
        ("", "--server must not be empty"),
        ("SLURM.json", "--server must name a .yaml file"),
    ],
)
def test_update_config_server_option_rejects_invalid_values(
    tmp_path, server_arg, expected
):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path, ["-s", server_arg])

    assert result.exit_code != 0
    assert expected in result.output
    assert slurm.read_text(encoding="utf-8") == before
    assert not (tmp_path / "server" / "SLURM.json").exists()


def test_update_config_prefers_same_name_template(tmp_path):
    local = _write_server_yaml(
        tmp_path,
        "local.yaml",
        """
        SERVER:
            SCHEDULER: UNKNOWN
        """,
    )

    result = _invoke_update_config(tmp_path, ["-s", "local"])

    assert result.exit_code == 0, result.output
    assert _read_yaml(local)["XTB"] == _template_yaml("local.yaml")["XTB"]


@pytest.mark.parametrize(
    "scheduler, template_name",
    [
        ("SLURM", "SLURM.yaml"),
        ("pbs", "PBS.yaml"),
        ("local", "local.yaml"),
        ("none", "local.yaml"),
        ("null", "local.yaml"),
        ("Null", "local.yaml"),
    ],
)
def test_update_config_custom_file_matches_scheduler(
    tmp_path, scheduler, template_name
):
    custom = _write_server_yaml(
        tmp_path,
        "my_cluster.yaml",
        f"""
        SERVER:
            SCHEDULER: {scheduler}
        """,
    )

    result = _invoke_update_config(tmp_path, ["-s", "my_cluster"])

    assert result.exit_code == 0, result.output
    assert _read_yaml(custom)["XTB"] == _template_yaml(template_name)["XTB"]


def test_update_config_copies_multiple_missing_programs(tmp_path):
    slurm = _write_server_yaml(
        tmp_path, "SLURM.yaml", OLD_SERVER_WITH_CUSTOM_GAUSSIAN
    )

    result = _invoke_update_config(tmp_path)
    data = _read_yaml(slurm)
    template = _template_yaml()

    assert result.exit_code == 0, result.output
    for program in ["ORCA", "XTB", "NCIPLOT"]:
        assert data[program] == template[program]
        assert f"added program: {program}" in result.output


def test_update_config_future_template_program_is_discovered(
    tmp_path, monkeypatch
):
    from chemsmart.cli import update_config as update_config_module

    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    templates = update_config_module._load_server_templates()
    for template_data in templates.values():
        template_data["CREST"] = {"EXEFOLDER": None, "LOCAL_RUN": True}
    monkeypatch.setattr(
        update_config_module, "_load_server_templates", lambda: templates
    )

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert _read_yaml(slurm)["CREST"] == {
        "EXEFOLDER": None,
        "LOCAL_RUN": True,
    }


def test_update_config_preserves_existing_sections_and_custom_fields(tmp_path):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", PRESERVED_SECTIONS_YAML)
    before = _read_yaml(slurm)

    result = _invoke_update_config(tmp_path)
    after = _read_yaml(slurm)

    assert result.exit_code == 0, result.output
    assert after["SERVER"] == before["SERVER"]
    assert after["GAUSSIAN"] == {"EXEFOLDER": "/custom/g16"}
    assert after["XTB"] == {"EXEFOLDER": "/old/xtb"}
    assert after["CUSTOM_FIELD"] == {"VALUE": "keep"}
    assert after["ORCA"] == _template_yaml()["ORCA"]
    assert after["NCIPLOT"] == _template_yaml()["NCIPLOT"]


def test_update_config_repeated_execution_is_idempotent(tmp_path):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)

    first = _invoke_update_config(tmp_path)
    after_first = slurm.read_text(encoding="utf-8")
    second = _invoke_update_config(tmp_path)

    assert first.exit_code == 0, first.output
    assert second.exit_code == 0, second.output
    assert slurm.read_text(encoding="utf-8") == after_first
    assert "already up to date" in second.output


def test_update_config_complete_yaml_does_not_rewrite_file(tmp_path):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", _template_text())
    before_mtime = slurm.stat().st_mtime_ns

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert slurm.stat().st_mtime_ns == before_mtime
    assert "already up to date" in result.output


@pytest.mark.parametrize("content", ["SERVER: [", "", "- not-a-mapping\n"])
def test_update_config_invalid_yaml_roots_are_not_overwritten(
    tmp_path, content
):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", content)
    before = slurm.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert slurm.read_text(encoding="utf-8") == before


@pytest.mark.parametrize(
    "content",
    [
        """
        SERVER:
            SCHEDULER: UNKNOWN
        """,
        """
        SERVER: []
        """,
        """
        SERVER: {}
        """,
    ],
)
def test_update_config_unmatched_scheduler_is_skipped_by_default(
    tmp_path, content
):
    custom = _write_server_yaml(
        tmp_path,
        "custom.yaml",
        content,
    )
    before = custom.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert (
        "skipped: could not match a bundled server template" in result.output
    )
    assert custom.read_text(encoding="utf-8") == before


def test_update_config_missing_server_directory_is_clean_error(tmp_path):
    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "Server config directory not found" in result.output
    assert not (tmp_path / "server").exists()
    assert not list(tmp_path.glob("*.yaml"))


def test_update_config_missing_ruamel_shows_install_hint(
    tmp_path, monkeypatch
):
    original_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "chemsmart.cli.update_config":
            raise ModuleNotFoundError(name="ruamel")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "ruamel.yaml is required" in result.output
    assert "pip install" in result.output
    assert "make env" in result.output
    assert "Traceback" not in result.output


def test_update_config_unrelated_import_error_is_not_mislabeled(
    tmp_path, monkeypatch
):
    original_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "chemsmart.cli.update_config":
            raise ModuleNotFoundError(name="unexpected_dependency")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    result = _invoke_update_config(tmp_path)

    assert isinstance(result.exception, ModuleNotFoundError)
    assert result.exception.name == "unexpected_dependency"
    assert "ruamel.yaml is required" not in result.output


def test_update_config_user_yaml_read_error_is_clean(tmp_path, monkeypatch):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")
    original_read_text = Path.read_text

    def fake_read_text(path, *args, **kwargs):
        if path == slurm:
            raise PermissionError(errno.EACCES, "Permission denied")
        return original_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", fake_read_text)

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "Could not read SLURM.yaml" in result.output
    assert "Permission denied" in result.output
    assert "Traceback" not in result.output
    assert original_read_text(slurm, encoding="utf-8") == before


@pytest.mark.parametrize(
    "template_content, expected",
    [
        ("SERVER: [", "is not valid YAML"),
        ("- invalid\n- template\n", "must contain a mapping"),
    ],
)
def test_load_server_templates_rejects_invalid_template(
    tmp_path, monkeypatch, template_content, expected
):
    from chemsmart.cli import update_config as update_config_module

    template_dir = tmp_path / "templates" / ".chemsmart" / "server"
    template_dir.mkdir(parents=True)
    (template_dir / "INVALID.yaml").write_text(
        template_content, encoding="utf-8"
    )
    monkeypatch.setattr(
        update_config_module.resources, "files", lambda package: tmp_path
    )

    with pytest.raises(
        update_config_module.ConfigUpdateError, match=expected
    ):
        update_config_module._load_server_templates()


def test_update_config_unknown_scheduler_errors_when_selected(tmp_path):
    custom = _write_server_yaml(
        tmp_path,
        "custom.yaml",
        """
        SERVER:
            SCHEDULER: UNKNOWN
        """,
    )
    before = custom.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path, ["-s", "custom"])

    assert result.exit_code != 0
    assert "Could not match a bundled template" in result.output
    assert custom.read_text(encoding="utf-8") == before


def test_update_config_respects_relative_config_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config_root = tmp_path / "cfg"
    slurm = _write_server_yaml(config_root, "SLURM.yaml", MINIMAL_SLURM)

    result = _invoke_update_config(config_root, env_value="cfg")
    assert result.exit_code == 0, result.output
    assert "XTB" in _read_yaml(slurm)


def test_update_config_preserves_comments_order_and_block_scalar(tmp_path):
    slurm = _write_server_yaml(
        tmp_path,
        "SLURM.yaml",
        """
        # keep file comment
        SERVER:
            SCHEDULER: SLURM
            EXTRA_COMMANDS: |
                echo keep

        CUSTOM_FIELD:
            VALUE: keep

        GAUSSIAN:
            EXEFOLDER: /custom/g16
        """,
    )

    result = _invoke_update_config(tmp_path)
    text = slurm.read_text(encoding="utf-8")

    assert result.exit_code == 0, result.output
    assert "# keep file comment" in text
    assert "EXTRA_COMMANDS: |" in text
    assert "echo keep" in text
    assert text.index("SERVER:") < text.index("CUSTOM_FIELD:")
    assert text.index("CUSTOM_FIELD:") < text.index("GAUSSIAN:")
    assert text.index("GAUSSIAN:") < text.index("ORCA:")


def test_update_config_refuses_symlink_without_partial_writes(tmp_path):
    server_dir = _server_dir(tmp_path)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_SLURM.replace("SCHEDULER: SLURM", "SCHEDULER: PBS"),
    )
    target = tmp_path / "target_SLURM.yaml"
    target.write_text(dedent(MINIMAL_SLURM).lstrip(), encoding="utf-8")
    symlink = server_dir / "SLURM.yaml"
    try:
        symlink.symlink_to(target)
    except (OSError, NotImplementedError) as exc:
        pytest.skip(f"Symlinks are not available in this environment: {exc}")
    pbs_before = pbs.read_text(encoding="utf-8")
    target_before = target.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "symlink" in result.output.lower()
    assert symlink.is_symlink()
    assert pbs.read_text(encoding="utf-8") == pbs_before
    assert target.read_text(encoding="utf-8") == target_before
    assert not list(server_dir.glob(".*.tmp"))


def test_update_config_mkstemp_error_is_clean_click_error(
    tmp_path, monkeypatch
):
    from chemsmart.cli import update_config as update_config_module

    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    def raise_permission_error(*args, **kwargs):
        raise PermissionError("Permission denied")

    monkeypatch.setattr(
        update_config_module.tempfile, "mkstemp", raise_permission_error
    )

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "Could not update SLURM.yaml: Permission denied" in result.output
    assert "Traceback" not in result.output
    assert slurm.read_text(encoding="utf-8") == before


def test_update_config_replace_error_cleans_temp_file(tmp_path, monkeypatch):
    from chemsmart.cli import update_config as update_config_module

    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    def raise_permission_error(*args, **kwargs):
        raise PermissionError("Permission denied")

    monkeypatch.setattr(
        update_config_module.os, "replace", raise_permission_error
    )

    result = _invoke_update_config(tmp_path)
    assert result.exit_code != 0
    assert "Could not update SLURM.yaml: Permission denied" in result.output
    assert "Traceback" not in result.output
    assert slurm.read_text(encoding="utf-8") == before
    assert not list((tmp_path / "server").glob(".SLURM.yaml.*.tmp"))


def test_updated_yaml_can_be_read_by_existing_settings_classes(
    tmp_path, monkeypatch
):
    from chemsmart.io.yaml import YAMLFile
    from chemsmart.settings import executable as executable_module
    from chemsmart.settings.executable import XTBExecutable
    from chemsmart.settings.server import Server
    from chemsmart.settings.user import CHEMSMARTUserSettings

    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    result = _invoke_update_config(tmp_path)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
    monkeypatch.setattr(
        executable_module, "user_settings", CHEMSMARTUserSettings()
    )

    assert result.exit_code == 0, result.output
    assert YAMLFile(filename=str(slurm)).yaml_contents_dict["XTB"]
    assert Server.from_yaml(str(slurm)).scheduler == "SLURM"
    assert XTBExecutable.from_servername("SLURM").get_executable() == "xtb"
