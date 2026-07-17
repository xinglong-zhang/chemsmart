from __future__ import annotations

import asyncio
import json

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.screens.chat import (
    _find_project_yaml_candidate_for_write,
    _latest_project_yaml_candidate,
)
from chemsmart.agent.tui.widgets.composer import Composer

from .._agent_session_helpers import FakeProvider
from ._helpers import (
    assert_matches_snapshot,
    normalize_svg,
    write_session_fixture,
)


def _set_composer_text(app: ChemsmartTuiApp, text: str) -> None:
    composer = app.query_one(Composer)
    composer.load_text(text)


def test_phase1_slash_commands_match_snapshots(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)

    def fake_get_provider():
        return FakeProvider([])

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("ai_api_key", "test-key")
    monkeypatch.setenv("CHEMSMART_AGENT_TUI_MODE", "run")

    async def snapshot_for(
        command: str, name: str, *, patch_exit: bool = False
    ):
        app = ChemsmartTuiApp(session_root=session_root)
        exit_calls: list[bool] = []
        if patch_exit:
            app.exit = lambda *args, **kwargs: exit_calls.append(True)  # type: ignore[method-assign]
        async with app.run_test() as pilot:
            await pilot.pause()
            _set_composer_text(app, command)
            await pilot.press("enter")
            await pilot.pause(0.2)
            if command == "/sessions":
                await pilot.pause()
            if command.startswith("/resume"):
                await pilot.pause(0.2)
            if patch_exit:
                assert exit_calls == [True]
            assert_matches_snapshot(
                name, normalize_svg(app.export_screenshot())
            )

    async def scenario() -> None:
        await snapshot_for("/help", "slash_help")
        await snapshot_for("/sessions", "slash_sessions")
        await snapshot_for("/resume session-001", "slash_resume")
        await snapshot_for("/clear", "slash_clear")
        await snapshot_for("/quit", "slash_quit", patch_exit=True)

    asyncio.run(scenario())


def test_init_alias_uses_unified_project_yaml_request(monkeypatch, tmp_path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)

    def fake_get_provider():
        return FakeProvider([])

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setenv("AI_PROVIDER", "openai")
    monkeypatch.setenv("ai_api_key", "test-key")
    monkeypatch.setenv("CHEMSMART_AGENT_TUI_MODE", "run")

    class _FakeWorker:
        is_finished = True

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            captured: list[str] = []
            monkeypatch.setattr(
                screen,
                "run_unified_session",
                lambda text: captured.append(text) or _FakeWorker(),
            )

            # Bare /init is now a helper/alias, not a persistent mode.
            _set_composer_text(app, "/init")
            await pilot.press("enter")
            await pilot.pause()
            assert screen._build_mode is False

            # A subsequent natural message routes through the unified session.
            _set_composer_text(
                app, "Optimize in water with B3LYP-D3BJ/def2-SVP, call it h2o."
            )
            await pilot.press("enter")
            await pilot.pause()
            assert captured == [
                "Optimize in water with B3LYP-D3BJ/def2-SVP, call it h2o."
            ]
            assert screen._build_mode is False

            # /init with an argument injects project-YAML intent.
            _set_composer_text(app, "/init Use Gaussian B3LYP/def2-SVP, call it co2.")
            await pilot.press("enter")
            await pilot.pause()
            assert captured[-1].startswith(
                "Build a workspace project YAML from this reported method."
            )
            assert "call it co2" in captured[-1]

    asyncio.run(scenario())


def test_latest_project_yaml_candidate_reads_tool_loop_log(tmp_path):
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    yaml_text = (
        "gas:\n"
        "  functional: b3lyp empiricaldispersion=gd3bj\n"
        "  basis: gen\n"
        "  freq: true\n"
    )
    entries = [
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "render_project_yaml",
                "status": "ok",
                "payload": {
                    "ok": True,
                    "project_name": "co2",
                    "program": "gaussian",
                    "yaml_text": yaml_text,
                },
            },
        },
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "validate_project_yaml",
                "status": "ok",
                "payload": {
                    "ok": True,
                    "project_name": "co2",
                    "program": "gaussian",
                    "verdict": "ok",
                    "issues": [],
                },
            },
        }
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )

    candidate = _latest_project_yaml_candidate(session_dir)

    assert candidate == {
        "project_name": "co2",
        "program": "gaussian",
        "yaml_text": yaml_text,
    }


def _render_entry(yaml_text: str, project="co2", program="gaussian"):
    return {
        "kind": "tool_use_result",
        "payload": {
            "tool": "render_project_yaml",
            "status": "ok",
            "payload": {
                "ok": True,
                "project_name": project,
                "program": program,
                "yaml_text": yaml_text,
            },
        },
    }


def _validate_request_entry(
    yaml_text: str,
    *,
    call_id: str = "validate-1",
    project: str = "co2",
    program: str = "gaussian",
):
    return {
        "kind": "tool_use_request",
        "payload": {
            "tool": "validate_project_yaml",
            "provider_call_id": call_id,
            "args": {
                "project_name": project,
                "program": program,
                "yaml_text": yaml_text,
            },
        },
    }


def _validate_entry(verdict="ok", *, call_id: str | None = None):
    payload = {
        "tool": "validate_project_yaml",
        "status": "ok",
        "payload": {"ok": verdict == "ok", "verdict": verdict},
    }
    if call_id is not None:
        payload["provider_call_id"] = call_id
    return {
        "kind": "tool_use_result",
        "payload": payload,
    }


def test_candidate_survives_identical_rerender_after_validate(tmp_path):
    # Build-mode over-iteration: the model re-renders the SAME candidate after a
    # successful validate. That must NOT make the candidate "unvalidated" (the
    # /write-project "Project write unavailable" chain bug).
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    yaml_text = "gas:\n  functional: b3lyp\n  basis: gen\n  freq: true\n"
    entries = [
        _render_entry(yaml_text),
        _validate_entry("ok"),
        _render_entry(yaml_text),  # re-render, last tool call
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(e) for e in entries) + "\n", encoding="utf-8"
    )

    candidate = _latest_project_yaml_candidate(session_dir)
    assert candidate is not None
    assert candidate["yaml_text"] == yaml_text


def test_changed_candidate_after_validate_requires_revalidation(tmp_path):
    # A genuinely NEW/changed candidate rendered after a validate must NOT be
    # treated as validated until it is itself validated.
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    first = "gas:\n  functional: b3lyp\n  basis: def2svp\n"
    changed = "gas:\n  functional: m062x\n  basis: def2tzvp\n"
    entries = [
        _render_entry(first),
        _validate_entry("ok"),
        _render_entry(changed),  # different candidate, not yet validated
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(e) for e in entries) + "\n", encoding="utf-8"
    )

    assert _latest_project_yaml_candidate(session_dir) is None


def test_candidate_uses_yaml_from_matching_validation_request(tmp_path):
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    rejected_yaml = (
        "gas:\n"
        "  functional: null\n"
        "  basis: def2svp\n"
        "  freq: true\n"
    )
    validated_yaml = (
        "gas:\n"
        "  functional: pbe0\n"
        "  basis: def2-SVP\n"
        "  freq: true\n"
    )
    rejected_render = _render_entry(rejected_yaml, project="water")
    rejected_render["payload"]["payload"].update(
        {
            "ok": False,
            "validation": {"verdict": "reject"},
        }
    )
    entries = [
        rejected_render,
        _validate_request_entry(
            validated_yaml,
            call_id="validate-pbe0",
            project="water",
        ),
        _validate_entry("ok", call_id="validate-pbe0"),
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )

    candidate = _latest_project_yaml_candidate(session_dir)

    assert candidate == {
        "project_name": "water",
        "program": "gaussian",
        "yaml_text": validated_yaml,
    }


def test_rejected_render_cannot_be_promoted_by_unbound_validation(tmp_path):
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    rejected_yaml = (
        "gas:\n"
        "  functional: null\n"
        "  basis: def2svp\n"
    )
    rejected_render = _render_entry(rejected_yaml, project="water")
    rejected_render["payload"]["payload"].update(
        {
            "ok": False,
            "validation": {"verdict": "reject"},
        }
    )
    entries = [rejected_render, _validate_entry("ok")]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )

    assert _latest_project_yaml_candidate(session_dir) is None


def test_validation_result_does_not_borrow_different_request(tmp_path):
    session_dir = tmp_path / "session"
    session_dir.mkdir()
    yaml_text = "gas:\n  functional: pbe0\n  basis: def2-SVP\n"
    entries = [
        _validate_request_entry(yaml_text, call_id="request-a"),
        _validate_entry("ok", call_id="result-b"),
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )

    assert _latest_project_yaml_candidate(session_dir) is None


def test_project_yaml_write_candidate_falls_back_to_latest_session(tmp_path):
    stale = tmp_path / "stale-session"
    stale.mkdir()
    (stale / "decision_log.jsonl").write_text("", encoding="utf-8")
    latest = tmp_path / "latest-session"
    latest.mkdir()
    yaml_text = "gas:\n  functional: b3lyp\n  basis: 6-31g(d)\n"
    entries = [
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "render_project_yaml",
                "status": "ok",
                "payload": {
                    "ok": True,
                    "project_name": "co2",
                    "program": "gaussian",
                    "yaml_text": yaml_text,
                },
            },
        },
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "validate_project_yaml",
                "status": "ok",
                "payload": {"verdict": "ok"},
            },
        },
    ]
    (latest / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )

    candidate = _find_project_yaml_candidate_for_write(
        tmp_path,
        preferred_session_dir=None,
    )

    assert candidate == {
        "project_name": "co2",
        "program": "gaussian",
        "yaml_text": yaml_text,
    }


def test_init_slash_command_is_registered_in_palette():
    from chemsmart.agent.tui.screens.chat import _SLASH_PALETTE_COMMANDS

    commands = {name for name, _desc in _SLASH_PALETTE_COMMANDS}
    assert "/init" in commands
