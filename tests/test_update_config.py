import builtins
from pathlib import Path
from textwrap import dedent

import yaml
from click.testing import CliRunner

import chemsmart.cli.update as update_module
from chemsmart.cli.update import ConfigurationUpdater, update


def _write_server_yaml(config_root: Path, name: str, content: str) -> Path:
    server_dir = config_root / "server"
    server_dir.mkdir(parents=True, exist_ok=True)
    path = server_dir / name
    path.write_text(dedent(content).lstrip(), encoding="utf-8")
    return path


def _invoke_update_config(config_root: Path, args=None):
    return CliRunner().invoke(
        update,
        ["config"] + (args or []),
        env={"CHEMSMART_CONFIG_DIR": str(config_root)},
    )


def _read_yaml(path: Path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def test_update_config_updates_selected_files_by_missing_program(
    tmp_path, monkeypatch
):
    fake_home = tmp_path / "home"
    monkeypatch.setenv("HOME", str(fake_home))
    template = {
        "SERVER": {"TEMPLATE_VALUE": "should-not-be-copied"},
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
            USER_VALUE: preserved
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
            USER_VALUE: preserved
        FAKE_PROGRAM_A:
            EXEFOLDER: /old/b-fake-a
        FAKE_PROGRAM_C:
            EXEFOLDER: /old/b-fake-c
        """,
    )
    c_yaml = _write_server_yaml(
        tmp_path,
        "C.yaml",
        """
        SERVER:
            USER_VALUE: preserved
        FAKE_PROGRAM_C:
            EXEFOLDER: /old/c-fake-c
        """,
    )
    unselected = _write_server_yaml(
        tmp_path,
        "unselected.yaml",
        """
        SERVER:
            USER_VALUE: preserved
        """,
    )
    unselected_before = unselected.read_text(encoding="utf-8")

    prompts = []
    answers = iter(["~/path/to/fake-a", "nUlL"])

    def prompt_for_program(text, **kwargs):
        prompts.append(text)
        return next(answers)

    monkeypatch.setattr(update_module.click, "prompt", prompt_for_program)

    result = _invoke_update_config(
        tmp_path, ["-s", "A", "-s", "B.yaml", "-s", "C"]
    )

    assert result.exit_code == 0, result.output
    assert "A.yaml:" in result.output
    assert "B.yaml:" in result.output
    assert "C.yaml:" in result.output
    assert "unselected.yaml:" not in result.output
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
    expanded_fake_a = str(fake_home / "path" / "to" / "fake-a")
    assert a_data["SERVER"] == {"USER_VALUE": "preserved"}
    assert a_data["FAKE_PROGRAM_A"] == {
        "EXEFOLDER": expanded_fake_a,
        "LOCAL_RUN": True,
    }
    assert a_data["FAKE_PROGRAM_B"] == {"EXEFOLDER": "/old/a-fake-b"}
    assert a_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/a-fake-c"}
    assert b_data["FAKE_PROGRAM_A"] == {"EXEFOLDER": "/old/b-fake-a"}
    assert b_data["FAKE_PROGRAM_B"] == {
        "EXEFOLDER": None,
        "LOCAL_RUN": True,
    }
    assert b_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/b-fake-c"}
    assert c_data["FAKE_PROGRAM_A"] == {
        "EXEFOLDER": expanded_fake_a,
        "LOCAL_RUN": True,
    }
    assert c_data["FAKE_PROGRAM_B"] == {
        "EXEFOLDER": None,
        "LOCAL_RUN": True,
    }
    assert c_data["FAKE_PROGRAM_C"] == {"EXEFOLDER": "/old/c-fake-c"}
    assert "EXEFOLDER: Null" in b_yaml.read_text(encoding="utf-8")
    assert "EXEFOLDER: Null" in c_yaml.read_text(encoding="utf-8")
    assert unselected.read_text(encoding="utf-8") == unselected_before


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
