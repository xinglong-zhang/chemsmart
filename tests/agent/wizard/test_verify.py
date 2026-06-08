from __future__ import annotations

from pathlib import Path

from chemsmart.agent.wizard.verify import verify_server_yaml

_BASE_YAML = """SERVER:
  SCHEDULER: SLURM
  QUEUE_NAME: debug
  NUM_HOURS: 8
  MEM_GB: 64
  NUM_CORES: 16
  SUBMIT_COMMAND: sbatch
  USE_HOSTS: true
GAUSSIAN:
  EXEFOLDER: /apps/gaussian
  LOCAL_RUN: true
  SCRATCH: true
  CONDA_ENV: ""
  MODULES: ""
  ENVARS: ""
"""


def _write_server_yaml(tmp_path: Path, name: str, contents: str) -> Path:
    target = tmp_path / ".chemsmart" / "server" / f"{name}.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(contents, encoding="utf-8")
    return target


def test_verify_server_yaml_mode_a(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    _write_server_yaml(
        tmp_path,
        "local-cluster",
        _BASE_YAML.replace(
            "  USE_HOSTS: true\n", "  HOST: localhost\n  USE_HOSTS: true\n"
        ),
    )

    result = verify_server_yaml("local-cluster")

    assert result.host == "localhost"
    assert result.mode == "local"
    assert result.would_submit_via == "local"
    assert result.transport_invocation == [
        "sbatch",
        "/tmp/chemsmart-wizard-verify.sh",
    ]
    assert result.warnings == []
    assert result.errors == []


def test_verify_server_yaml_mode_b(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    _write_server_yaml(
        tmp_path,
        "remote-cluster",
        _BASE_YAML.replace(
            "  USE_HOSTS: true\n", "  HOST: cluster\n  USE_HOSTS: true\n"
        ),
    )

    result = verify_server_yaml("remote-cluster")

    assert result.host == "cluster"
    assert result.mode == "ssh"
    assert result.would_submit_via == "ssh"
    assert result.transport_invocation is not None
    assert result.transport_invocation[0] == "ssh"
    assert result.transport_invocation[1] == "cluster"
    assert result.errors == []


def test_verify_server_yaml_warns_for_legacy_yaml_without_host(
    monkeypatch, tmp_path
):
    monkeypatch.setenv("HOME", str(tmp_path))
    _write_server_yaml(tmp_path, "legacy-cluster", _BASE_YAML)

    result = verify_server_yaml("legacy-cluster")

    assert result.host is None
    assert result.mode == "local"
    assert result.would_submit_via == "ssh"
    assert result.transport_invocation is not None
    assert result.transport_invocation[0] == "ssh"
    assert result.warnings
    assert any(
        "SERVER.HOST is missing" in warning for warning in result.warnings
    )
    assert any(
        "filename-stem fallback" in warning for warning in result.warnings
    )
    assert result.errors == []


def test_verify_server_yaml_reports_missing_server(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))

    result = verify_server_yaml("does-not-exist")

    assert result.host is None
    assert result.mode == "unknown"
    assert result.would_submit_via == "unknown"
    assert result.transport_invocation is None
    assert result.warnings == []
    assert result.errors
    assert "Server YAML not found" in result.errors[0]
