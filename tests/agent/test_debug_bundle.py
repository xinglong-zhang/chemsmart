from __future__ import annotations

import json
import tarfile
from pathlib import Path

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.services.debug_bundle import create_debug_bundle


def _write_session(root: Path, session_id: str, workspace: Path) -> Path:
    session = root / session_id
    session.mkdir(parents=True)
    (session / "decision_log.jsonl").write_text(
        json.dumps(
            {
                "kind": "provider_turn_error",
                "payload": {
                    "message": "Bearer sensitive-token",
                    "workspace": str(workspace),
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (session / "runtime_state.json").write_text(
        json.dumps(
            {
                "phase": "blocked",
                "api_key": "sk-super-secret",
                "home_file": str(Path.home() / "private.txt"),
            }
        ),
        encoding="utf-8",
    )
    (session / "session_metadata.json").write_text(
        json.dumps({"blocked": True, "token": "secret-token"}),
        encoding="utf-8",
    )
    (session / "agent.yaml").write_text(
        "api_key: must-never-enter-bundle\n",
        encoding="utf-8",
    )
    (session / "api.env").write_text(
        "ai_api_key=must-never-enter-bundle\n",
        encoding="utf-8",
    )
    return session


def test_debug_bundle_is_explicit_minimal_and_redacted(tmp_path):
    root = tmp_path / "sessions"
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    _write_session(root, "session-123", workspace)
    output = tmp_path / "bundle.tar.gz"

    receipt = create_debug_bundle(
        root,
        "session-123",
        output,
        workspace_path=workspace,
    )

    assert receipt["session_id"] == "session-123"
    with tarfile.open(output, "r:gz") as archive:
        names = set(archive.getnames())
        assert names == {
            "session/decision_log.jsonl",
            "session/runtime_state.json",
            "session/session_metadata.json",
            "environment_manifest.json",
        }
        payload = "\n".join(
            archive.extractfile(name).read().decode("utf-8")
            for name in sorted(names)
        )
    assert "sensitive-token" not in payload
    assert "sk-super-secret" not in payload
    assert "secret-token" not in payload
    assert "must-never-enter-bundle" not in payload
    assert str(Path.home()) not in payload
    assert str(workspace) not in payload
    assert "<HOME>" in payload
    assert "<WORKSPACE>" in payload
    assert "<REDACTED>" in payload


def test_debug_bundle_cli_rejects_latest(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "chemsmart.agent.cli_commands._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    result = CliRunner().invoke(
        agent,
        [
            "debug-bundle",
            "latest",
            "--output",
            str(tmp_path / "bundle.tar.gz"),
        ],
    )

    assert result.exit_code != 0
    assert "explicit session ID" in result.output


def test_debug_bundle_cli_packages_requested_session(monkeypatch, tmp_path):
    root = tmp_path / "sessions"
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    _write_session(root, "chosen-session", workspace)
    _write_session(root, "newer-unrelated-session", workspace)
    output = tmp_path / "bundle.tar.gz"
    monkeypatch.setattr(
        "chemsmart.agent.cli_commands._default_session_root",
        lambda: str(root),
    )
    monkeypatch.chdir(workspace)

    result = CliRunner().invoke(
        agent,
        [
            "debug-bundle",
            "chosen-session",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.output
    assert "session: chosen-session" in result.output
    with tarfile.open(output, "r:gz") as archive:
        manifest = json.loads(
            archive.extractfile("environment_manifest.json")
            .read()
            .decode("utf-8")
        )
    assert manifest["session_id"] == "chosen-session"
