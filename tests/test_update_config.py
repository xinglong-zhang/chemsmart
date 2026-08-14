import builtins
import errno
from importlib import resources
from pathlib import Path
from textwrap import dedent

import pytest
import yaml
from click.testing import CliRunner

import chemsmart.cli.update as update_module
from chemsmart.cli.main import entry_point
from chemsmart.cli.update import ConfigurationUpdater, update


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


def _template_yaml():
    template = (
        resources.files("chemsmart.settings") / "templates" / "server.yaml"
    )
    return yaml.safe_load(template.read_text(encoding="utf-8"))


MINIMAL_SLURM = """
SERVER:
    SCHEDULER: SLURM
"""


MINIMAL_PBS = """
SERVER:
    SCHEDULER: PBS
"""


SERVER_WITH_CUSTOM_FIELDS = """
SERVER:
    SCHEDULER: SLURM
    QUEUE_NAME: custom

GAUSSIAN:
    EXEFOLDER: /custom/g16

CUSTOM_FIELD:
    VALUE: keep
"""


def test_update_config_registered_on_main_cli(tmp_path):
    result = CliRunner().invoke(
        entry_point,
        ["update", "config", "--help"],
        env={"CHEMSMART_CONFIG_DIR": str(tmp_path)},
    )

    assert result.exit_code == 0
    assert "--server" in result.output


def test_update_config_default_updates_existing_yaml_only(tmp_path):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_PBS,
    )

    result = _invoke_update_config(tmp_path)
    server_dir = tmp_path / "server"

    assert result.exit_code == 0, result.output
    assert "SLURM.yaml:" in result.output
    assert "PBS.yaml:" in result.output
    template = _template_yaml()
    for path in (slurm, pbs):
        data = _read_yaml(path)
        for program in ConfigurationUpdater._template_programs(template):
            assert data[program] == template[program]
            assert f"added program: {program}" in result.output
    assert {path.name for path in server_dir.glob("*.yaml")} == {
        "SLURM.yaml",
        "PBS.yaml",
    }


def test_update_config_server_option_accepts_multiple_files(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FAKE_PROGRAM_A": {"EXEFOLDER": None, "LOCAL_RUN": True},
        "FAKE_PROGRAM_B": {"EXEFOLDER": None, "LOCAL_RUN": True},
        "FAKE_PROGRAM_C": {"EXEFOLDER": None, "LOCAL_RUN": True},
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    a_yaml = _write_server_yaml(
        tmp_path,
        "A.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        FAKE_PROGRAM_B:
            EXEFOLDER: /old/a-fake-b
        FAKE_PROGRAM_C:
            EXEFOLDER: /old/a-fake-c
        """,
    )
    b_yaml = _write_server_yaml(
        tmp_path,
        "B.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        FAKE_PROGRAM_A:
            EXEFOLDER: /old/b-fake-a
        FAKE_PROGRAM_C:
            EXEFOLDER: /old/b-fake-c
        """,
    )
    c_yaml = _write_server_yaml(tmp_path, "C.yaml", MINIMAL_SLURM)
    c_before = c_yaml.read_text(encoding="utf-8")
    prompts = []
    answers = iter(["/new/fake-a", "/new/fake-b"])

    def prompt_for_program(text, **kwargs):
        prompts.append(text)
        return next(answers)

    monkeypatch.setattr(update_module.click, "prompt", prompt_for_program)

    result = _invoke_update_config(tmp_path, ["-s", "A", "-s", "B.yaml"])

    assert result.exit_code == 0, result.output
    assert "A.yaml:" in result.output
    assert "B.yaml:" in result.output
    assert "C.yaml:" not in result.output
    assert len(prompts) == 2
    assert prompts[0].startswith("FAKE_PROGRAM_A EXEFOLDER")
    assert prompts[1].startswith("FAKE_PROGRAM_B EXEFOLDER")

    a_data = _read_yaml(a_yaml)
    b_data = _read_yaml(b_yaml)
    assert a_data["FAKE_PROGRAM_A"] == {
        "EXEFOLDER": "/new/fake-a",
        "LOCAL_RUN": True,
    }
    assert a_data["FAKE_PROGRAM_B"] == {"EXEFOLDER": "/old/a-fake-b"}
    assert a_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/a-fake-c"}
    assert b_data["FAKE_PROGRAM_A"] == {"EXEFOLDER": "/old/b-fake-a"}
    assert b_data["FAKE_PROGRAM_B"] == {
        "EXEFOLDER": "/new/fake-b",
        "LOCAL_RUN": True,
    }
    assert b_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/b-fake-c"}
    assert c_yaml.read_text(encoding="utf-8") == c_before


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


def test_update_config_uses_server_template_regardless_of_scheduler(tmp_path):
    custom = _write_server_yaml(
        tmp_path,
        "my_cluster.yaml",
        """
        SERVER:
            SCHEDULER: UNKNOWN
        """,
    )

    result = _invoke_update_config(tmp_path, ["-s", "my_cluster.yaml"])

    assert result.exit_code == 0, result.output
    data = _read_yaml(custom)
    template = _template_yaml()
    for program in ConfigurationUpdater._template_programs(template):
        assert data[program] == template[program]


def test_update_config_prompts_once_per_program_and_targets_missing_files(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FAKE_PROGRAM_A": {"EXEFOLDER": None, "LOCAL_RUN": True},
        "FAKE_PROGRAM_B": {"EXEFOLDER": None, "LOCAL_RUN": True},
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )

    a_yaml = _write_server_yaml(
        tmp_path,
        "A.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        FAKE_PROGRAM_B:
            EXEFOLDER: /old/a-fake-b
        """,
    )
    b_yaml = _write_server_yaml(
        tmp_path,
        "B.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        FAKE_PROGRAM_A:
            EXEFOLDER: /old/b-fake-a
        """,
    )
    c_yaml = _write_server_yaml(
        tmp_path,
        "C.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        """,
    )
    prompts = []
    answers = iter(["/new/fake-a", "/new/fake-b"])

    def prompt_for_program(text, **kwargs):
        prompts.append(text)
        return next(answers)

    monkeypatch.setattr(update_module.click, "prompt", prompt_for_program)

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert len(prompts) == 2
    assert prompts[0].startswith("FAKE_PROGRAM_A EXEFOLDER")
    assert prompts[1].startswith("FAKE_PROGRAM_B EXEFOLDER")
    assert not any(
        filename in prompt
        for prompt in prompts
        for filename in ("A.yaml", "B.yaml", "C.yaml")
    )

    a_data = _read_yaml(a_yaml)
    b_data = _read_yaml(b_yaml)
    c_data = _read_yaml(c_yaml)
    assert a_data["FAKE_PROGRAM_A"] == {
        "EXEFOLDER": "/new/fake-a",
        "LOCAL_RUN": True,
    }
    assert a_data["FAKE_PROGRAM_B"] == {"EXEFOLDER": "/old/a-fake-b"}
    assert b_data["FAKE_PROGRAM_A"] == {"EXEFOLDER": "/old/b-fake-a"}
    assert b_data["FAKE_PROGRAM_B"] == {
        "EXEFOLDER": "/new/fake-b",
        "LOCAL_RUN": True,
    }
    assert c_data["FAKE_PROGRAM_A"] == {
        "EXEFOLDER": "/new/fake-a",
        "LOCAL_RUN": True,
    }
    assert c_data["FAKE_PROGRAM_B"] == {
        "EXEFOLDER": "/new/fake-b",
        "LOCAL_RUN": True,
    }


def test_update_config_enter_keeps_server_template_exefolder(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FAKE_PROGRAM": {"EXEFOLDER": "/template/fake"},
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    prompt_count = 0

    def use_template_value(text, **kwargs):
        nonlocal prompt_count
        prompt_count += 1
        return ""

    monkeypatch.setattr(update_module.click, "prompt", use_template_value)
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_PBS,
    )

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert prompt_count == 1
    assert _read_yaml(slurm)["FAKE_PROGRAM"]["EXEFOLDER"] == "/template/fake"
    assert _read_yaml(pbs)["FAKE_PROGRAM"]["EXEFOLDER"] == "/template/fake"


def test_update_config_existing_program_missing_exefolder_is_untouched(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FAKE_PROGRAM": {
            "EXEFOLDER": "/template/fake",
            "LOCAL_RUN": True,
        },
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    monkeypatch.setattr(
        update_module.click,
        "prompt",
        lambda *args, **kwargs: pytest.fail("prompt should not be called"),
    )
    slurm = _write_server_yaml(
        tmp_path,
        "SLURM.yaml",
        """
        SERVER:
            SCHEDULER: SLURM
        FAKE_PROGRAM:
            LOCAL_RUN: false
        """,
    )
    before = slurm.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert slurm.read_text(encoding="utf-8") == before
    assert _read_yaml(slurm)["FAKE_PROGRAM"] == {"LOCAL_RUN": False}


def test_update_config_program_without_exefolder_is_copied_without_prompt(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FUTURE_PROGRAM": {"LOCAL_RUN": True},
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    monkeypatch.setattr(
        update_module.click,
        "prompt",
        lambda *args, **kwargs: pytest.fail("prompt should not be called"),
    )
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert _read_yaml(slurm)["FUTURE_PROGRAM"] == {"LOCAL_RUN": True}


def test_update_config_non_tty_uses_template_without_prompt(
    tmp_path, monkeypatch
):
    template = {
        "SERVER": {"SCHEDULER": "PBS"},
        "FAKE_PROGRAM": {"EXEFOLDER": "/template/fake"},
    }
    monkeypatch.setattr(
        ConfigurationUpdater, "_load_server_template", lambda self: template
    )
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: False)
    )
    monkeypatch.setattr(
        update_module.click,
        "prompt",
        lambda *args, **kwargs: pytest.fail("prompt should not be called"),
    )
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)

    result = _invoke_update_config(tmp_path)

    assert result.exit_code == 0, result.output
    assert _read_yaml(slurm)["FAKE_PROGRAM"] == {"EXEFOLDER": "/template/fake"}


def test_update_config_abort_during_prompt_writes_nothing(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    monkeypatch.setattr(
        update_module.click,
        "prompt",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            update_module.click.Abort()
        ),
    )
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    pbs = _write_server_yaml(
        tmp_path,
        "PBS.yaml",
        MINIMAL_PBS,
    )
    slurm_before = slurm.read_text(encoding="utf-8")
    pbs_before = pbs.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert slurm.read_text(encoding="utf-8") == slurm_before
    assert pbs.read_text(encoding="utf-8") == pbs_before


def test_update_config_preserves_existing_sections_and_custom_fields(tmp_path):
    slurm = _write_server_yaml(
        tmp_path, "SLURM.yaml", SERVER_WITH_CUSTOM_FIELDS
    )
    before = _read_yaml(slurm)

    result = _invoke_update_config(tmp_path)
    after = _read_yaml(slurm)

    assert result.exit_code == 0, result.output
    assert after["SERVER"] == before["SERVER"]
    assert after["GAUSSIAN"] == {"EXEFOLDER": "/custom/g16"}
    assert after["CUSTOM_FIELD"] == {"VALUE": "keep"}


def test_update_config_repeated_execution_is_idempotent(tmp_path, monkeypatch):
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    prompt_count = 0

    def use_template_value(text, **kwargs):
        nonlocal prompt_count
        prompt_count += 1
        return ""

    monkeypatch.setattr(update_module.click, "prompt", use_template_value)
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)

    first = _invoke_update_config(tmp_path)
    after_first = slurm.read_text(encoding="utf-8")
    after_first_mtime = slurm.stat().st_mtime_ns
    first_prompt_count = prompt_count
    second = _invoke_update_config(tmp_path)

    assert first.exit_code == 0, first.output
    assert second.exit_code == 0, second.output
    assert first_prompt_count > 0
    assert prompt_count == first_prompt_count
    assert slurm.read_text(encoding="utf-8") == after_first
    assert slurm.stat().st_mtime_ns == after_first_mtime
    assert "already up to date" in second.output


@pytest.mark.parametrize("content", ["SERVER: [", "", "- not-a-mapping\n"])
def test_update_config_invalid_yaml_roots_are_not_overwritten(
    tmp_path, monkeypatch, content
):
    monkeypatch.setattr(
        ConfigurationUpdater, "_is_interactive", staticmethod(lambda: True)
    )
    monkeypatch.setattr(
        update_module.click,
        "prompt",
        lambda *args, **kwargs: pytest.fail("prompt should not be called"),
    )
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", content)
    before = slurm.read_text(encoding="utf-8")

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert slurm.read_text(encoding="utf-8") == before


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
        if name == "ruamel.yaml":
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
        if name == "ruamel.yaml":
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
def test_load_server_template_rejects_invalid_template(
    tmp_path, monkeypatch, template_content, expected
):
    template_dir = tmp_path / "templates"
    template_dir.mkdir()
    (template_dir / "server.yaml").write_text(
        template_content, encoding="utf-8"
    )
    monkeypatch.setattr(
        update_module.resources, "files", lambda package: tmp_path
    )

    with pytest.raises(update_module.ConfigUpdateError, match=expected):
        ConfigurationUpdater()._load_server_template()


def test_update_config_respects_relative_config_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config_root = tmp_path / "cfg"
    slurm = _write_server_yaml(config_root, "SLURM.yaml", MINIMAL_SLURM)

    result = _invoke_update_config(config_root, env_value="cfg")
    assert result.exit_code == 0, result.output
    assert set(_template_yaml()) <= set(_read_yaml(slurm))


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
        MINIMAL_PBS,
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
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    def raise_permission_error(*args, **kwargs):
        raise PermissionError("Permission denied")

    monkeypatch.setattr(
        update_module.tempfile, "mkstemp", raise_permission_error
    )

    result = _invoke_update_config(tmp_path)

    assert result.exit_code != 0
    assert "Could not update SLURM.yaml: Permission denied" in result.output
    assert "Traceback" not in result.output
    assert slurm.read_text(encoding="utf-8") == before


def test_update_config_replace_error_cleans_temp_file(tmp_path, monkeypatch):
    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    before = slurm.read_text(encoding="utf-8")

    def raise_permission_error(*args, **kwargs):
        raise PermissionError("Permission denied")

    monkeypatch.setattr(update_module.os, "replace", raise_permission_error)

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
    from chemsmart.settings.executable import GaussianExecutable
    from chemsmart.settings.server import Server
    from chemsmart.settings.user import CHEMSMARTUserSettings

    slurm = _write_server_yaml(tmp_path, "SLURM.yaml", MINIMAL_SLURM)
    result = _invoke_update_config(tmp_path)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
    monkeypatch.setattr(
        executable_module, "user_settings", CHEMSMARTUserSettings()
    )

    assert result.exit_code == 0, result.output
    assert YAMLFile(filename=str(slurm)).yaml_contents_dict["GAUSSIAN"]
    assert Server.from_yaml(str(slurm)).scheduler == "SLURM"
    gaussian_folder = _template_yaml()["GAUSSIAN"]["EXEFOLDER"]
    assert GaussianExecutable.from_servername("SLURM").get_executable() == str(
        (Path(gaussian_folder) / "g16").expanduser()
    )
