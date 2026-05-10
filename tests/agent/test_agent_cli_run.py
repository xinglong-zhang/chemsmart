from __future__ import annotations

import json
from datetime import UTC, datetime, timedelta
from pathlib import Path

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.core import Plan

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import openai_final_response
from .tui._helpers import write_session_fixture


def test_agent_cli_run_prints_plan_and_writes_decision_log(
    monkeypatch,
    tmp_path: Path,
):
    def fake_run_loop(self, request, **kwargs):
        session_dir = tmp_path / "sessions" / "session-001"
        session_dir.mkdir(parents=True, exist_ok=True)
        (session_dir / "decision_log.jsonl").write_text(
            '{"kind":"request","payload":{"request":"optimize water"}}\n',
            encoding="utf-8",
        )
        return {
            "session_id": "session-001",
            "session_dir": str(session_dir),
            "plan": Plan(steps=[], rationale="Answer directly."),
            "plan_text": "Plan:",
            "critic_verdict": None,
            "completed_steps": 0,
            "blocked": False,
            "dry_run_result": None,
            "dry_run_results": [],
            "runtime_result": None,
            "preview_submit": None,
            "results": [],
            "assistant_output": "Done.",
            "tool_requests": [],
            "tool_outcomes": [],
            "loop_state": {},
            "final_message": "Done.",
            "limit_reason": None,
            "advisory_only": True,
            "is_chitchat": False,
            "approval_mode": "driving",
            "driving_mode": True,
            "yolo": False,
            "denials_count": 0,
            "approvals_count": 0,
        }

    monkeypatch.setattr(
        "chemsmart.agent.cli.AgentSession.run_loop", fake_run_loop
    )
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        ["run", "--dry-submit", "optimize water"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Advice:" in result.output
    assert "assistant: Done." in result.output
    assert "decision log:" in result.output

    session_dir = tmp_path / "sessions"
    decision_logs = list(session_dir.rglob("decision_log.jsonl"))
    assert len(decision_logs) == 1
    assert Path(decision_logs[0]).read_text()


def test_agent_tools_shows_descriptions(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["tools"], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert "Registered tools" in result.output
    assert "build_molecule" in result.output
    assert "Load one molecule from a" in result.output
    assert "chemsmart parsing." in result.output
    assert "filepath, index" in result.output


def test_resume_missing_emits_friendly_error(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["resume", "does-not-exist"])

    assert result.exit_code == 1
    assert (
        "Error: session 'does-not-exist' not found. List recent sessions "
        "with: chemsmart agent sessions"
    ) in result.output
    assert "Traceback" not in result.output


def test_agent_run_resume_missing_emits_friendly_error(
    monkeypatch, tmp_path: Path
):
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(agent, ["run", "--resume", "does-not-exist"])

    assert result.exit_code == 1
    assert (
        "Error: session 'does-not-exist' not found. List recent sessions "
        "with: chemsmart agent sessions"
    ) in result.output
    assert "Traceback" not in result.output


def test_agent_sessions_lists_recent(monkeypatch, tmp_path: Path):
    session_root = tmp_path / "sessions"
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(session_root),
    )

    now = datetime.now(UTC)
    for index in range(12):
        session_id = f"session-{index:02d}"
        session_dir = write_session_fixture(
            session_root, session_id=session_id
        )
        started_at = (now - timedelta(minutes=index + 5)).isoformat()
        ended_at = (now - timedelta(minutes=index)).isoformat()
        session_json = json.loads(
            (session_dir / "session.json").read_text(encoding="utf-8")
        )
        session_json["request"] = (
            f"Prepare session {index:02d} with a very long request string "
            f"that should be trimmed in the sessions table output."
        )
        session_json["started_at"] = started_at
        (session_dir / "session.json").write_text(
            json.dumps(session_json, indent=2), encoding="utf-8"
        )
        (session_dir / "state.json").write_text(
            json.dumps(session_json, indent=2), encoding="utf-8"
        )

        metadata = json.loads(
            (session_dir / "session_metadata.json").read_text(encoding="utf-8")
        )
        metadata["request"] = session_json["request"]
        metadata["started_at"] = started_at
        metadata["ended_at"] = ended_at
        metadata["blocked"] = index == 1
        metadata["critic_verdict"] = "ok"
        (session_dir / "session_metadata.json").write_text(
            json.dumps(metadata, indent=2), encoding="utf-8"
        )

    runner = CliRunner()
    result = runner.invoke(agent, ["sessions"], catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert "Agent sessions" in result.output
    assert "session-00" in result.output
    assert "session-09" in result.output
    assert "session-10" not in result.output
    assert "session-11" not in result.output
    assert "blocked" in result.output
    assert "Prepare session 00" in result.output


def test_agent_cli_run_wraps_runtime_errors_as_click_exceptions(
    monkeypatch,
):
    def fake_run_loop(self, request, **kwargs):
        raise RuntimeError(
            "run_local failed with returncode 1; see /tmp/run.stderr"
        )

    monkeypatch.setattr(
        "chemsmart.agent.cli.AgentSession.run_loop", fake_run_loop
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        ["run", "--dry-submit", "optimize examples/h2o.xyz"],
        catch_exceptions=False,
    )

    assert result.exit_code == 1
    assert "Error: run_local failed with returncode 1" in result.output
    assert "Traceback" not in result.output


def test_agent_cli_run_surfaces_advisory_only_answers(
    monkeypatch,
    tmp_path: Path,
):
    def fake_run_loop(self, request, **kwargs):
        return {
            "session_id": "session-002",
            "session_dir": str(tmp_path / "sessions" / "session-002"),
            "plan": Plan(
                steps=[],
                rationale=(
                    "Use M06-2X/def2-SVP for the TS search and refine with "
                    "def2-TZVP after one imaginary frequency check."
                ),
            ),
            "plan_text": "Plan:",
            "critic_verdict": None,
            "completed_steps": 0,
            "blocked": False,
            "dry_run_result": None,
            "dry_run_results": [],
            "runtime_result": None,
            "preview_submit": None,
            "results": [],
            "assistant_output": (
                "Use M06-2X/def2-SVP for the TS search and refine with "
                "def2-TZVP after one imaginary frequency check."
            ),
            "tool_requests": [],
            "tool_outcomes": [],
            "loop_state": {},
            "final_message": "Done.",
            "limit_reason": None,
            "advisory_only": True,
            "is_chitchat": False,
            "approval_mode": "driving",
            "driving_mode": True,
            "yolo": False,
            "denials_count": 0,
            "approvals_count": 0,
        }

    monkeypatch.setattr(
        "chemsmart.agent.cli.AgentSession.run_loop", fake_run_loop
    )
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        [
            "run",
            "--dry-submit",
            "Recommend a TS method and basis set for a Cope rearrangement.",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Advice:" in result.output
    assert "M06-2X/def2-SVP" in result.output
    assert "critic verdict:" not in result.output


def test_ask_renders_advisory_plan(
    monkeypatch,
    tmp_path: Path,
):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_final_response(
                    "Advisory text: start with wB97X-D/def2-SVP, confirm "
                    "one imaginary frequency, then refine with def2-TZVP."
                )
            }
        ]
    )

    def fake_get_provider():
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )
    monkeypatch.setattr(
        "chemsmart.agent.cli._default_session_root",
        lambda: str(tmp_path / "sessions"),
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        [
            "ask",
            "Recommend method/basis for a Cope rearrangement TS and explain trade-offs.",
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert "Assistant" in result.output
    assert "Advisory text" in result.output
    assert "Plan" not in result.output
