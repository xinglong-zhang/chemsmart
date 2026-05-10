from __future__ import annotations

import pytest

from chemsmart.agent.wizard.write import write_server_yaml

YAML_TEXT = """SERVER:
  SCHEDULER: SLURM
"""

UPDATED_TEXT = """SERVER:
  SCHEDULER: PBS
"""


def test_write_server_yaml_writes_to_user_server_dir(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))

    written = write_server_yaml("perlmutter", YAML_TEXT)

    expected = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert written == str(expected)
    assert expected.read_text(encoding="utf-8") == YAML_TEXT


def test_write_server_yaml_rejects_existing_file_without_overwrite(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("HOME", str(tmp_path))
    write_server_yaml("perlmutter", YAML_TEXT)

    with pytest.raises(FileExistsError):
        write_server_yaml("perlmutter", UPDATED_TEXT)


def test_write_server_yaml_overwrites_when_requested(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    write_server_yaml("perlmutter", YAML_TEXT)

    write_server_yaml("perlmutter", UPDATED_TEXT, overwrite=True)

    expected = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert expected.read_text(encoding="utf-8") == UPDATED_TEXT
