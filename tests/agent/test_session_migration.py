from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.models import SessionState
from chemsmart.agent.services.session_store import (
    LegacySessionFormatError,
    SessionMigrationError,
    load_current_session_state,
    migrate_legacy_session,
)
from chemsmart.agent.tui.services.job_poller import JobStateReader


def _write_legacy_session(root: Path, session_id: str = "legacy-001") -> Path:
    session_dir = root / session_id
    session_dir.mkdir(parents=True)
    SessionState(
        session_id=session_id,
        cwd=str(root),
        request="optimize water",
    ).save(session_dir / "state.json")
    (session_dir / "decision_log.jsonl").write_text(
        '{"kind":"request","payload":{"request":"optimize water"}}\n',
        encoding="utf-8",
    )
    return session_dir


def test_legacy_session_requires_explicit_migration(tmp_path: Path):
    source = _write_legacy_session(tmp_path)

    with pytest.raises(LegacySessionFormatError, match="migrate-session"):
        load_current_session_state(source, required=True)

    with pytest.raises(LegacySessionFormatError, match="migrate-session"):
        JobStateReader.load(source)


def test_migration_preserves_source_and_creates_canonical_copy(tmp_path: Path):
    source = _write_legacy_session(tmp_path)
    source_state = (source / "state.json").read_bytes()

    result = migrate_legacy_session(source)

    destination = Path(str(result["destination"]))
    assert destination.name == "legacy-001-migrated-v2"
    assert not (source / "session.json").exists()
    assert (source / "state.json").read_bytes() == source_state
    migrated = load_current_session_state(destination, required=True)
    assert migrated is not None
    assert migrated.session_id == destination.name
    manifest = json.loads(
        (destination / "migration_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == 2
    assert manifest["source_session_id"] == "legacy-001"
    assert "state.json" in manifest["source_artifact_sha256"]


def test_migration_refuses_existing_destination(tmp_path: Path):
    source = _write_legacy_session(tmp_path)
    destination = tmp_path / "existing"
    destination.mkdir()

    with pytest.raises(SessionMigrationError, match="already exists"):
        migrate_legacy_session(source, destination)


def test_cli_migrates_session_id_under_configured_root(
    monkeypatch,
    tmp_path: Path,
):
    session_root = tmp_path / "sessions"
    source = _write_legacy_session(session_root)
    monkeypatch.setattr(
        "chemsmart.agent.cli_commands._default_session_root",
        lambda: str(session_root),
    )

    result = CliRunner().invoke(agent, ["migrate-session", source.name])

    assert result.exit_code == 0, result.output
    assert "source preserved:" in result.output
    assert "legacy-001-migrated-v2" in result.output
    assert (session_root / "legacy-001-migrated-v2" / "session.json").is_file()
