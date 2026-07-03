from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.provider_config import AgentProviderConfig
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    ErrorCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.slash_palette import SlashCommandPalette
from chemsmart.agent.tui.widgets.transcript import Transcript


def _local_config() -> AgentProviderConfig:
    return AgentProviderConfig(
        name="local_chemsmart_v13_1_mlx4",
        type="local",
        api_key="",
        model="chemsmart-qwen2.5-coder-3b-instruct-v13_1-mlx-4bit",
        base_url="",
        extra_headers={},
        runtime="mlx",
        project="test",
    )


class _SemanticOk:
    def to_dict(self) -> dict:
        return {
            "verdict": "ok",
            "failed_rule_ids": [],
            "issues": [],
            "generated_inputs": [
                {
                    "path": "/tmp/h2o.com",
                    "route": "#p b3lyp/6-31g(d) sp",
                }
            ],
        }


class _SemanticReject:
    def to_dict(self) -> dict:
        return {
            "verdict": "reject",
            "failed_rule_ids": ["cmd.semantic.generated_input_missing"],
            "issues": [
                {
                    "rule_id": "cmd.semantic.generated_input_missing",
                    "severity": "reject",
                    "message": "safe dry-run did not produce an input file",
                    "evidence": {},
                }
            ],
            "generated_inputs": [],
        }


class _FakeSynthesisSession:
    requests: list[str] = []

    def __init__(self) -> None:
        self._last_semantic_result = _SemanticOk()
        self._last_raw_response = '{"intent":"workflow","jobs":[]}'

    def prepare_command(self, request: str) -> dict:
        self.requests.append(request)
        return {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p test -f examples/h2o.xyz "
                "-c 0 -m 1 -l h2o_01_001 sp"
            ),
            "explanation": "Prepared chemsmart command from compact SPEC.",
            "confidence": "high",
            "missing_info": [],
            "alternatives": [],
        }


class _RejectingSynthesisSession:
    requests: list[str] = []

    def __init__(self) -> None:
        self._last_semantic_result = _SemanticReject()
        self._last_raw_response = '{"intent":"workflow","jobs":[]}'

    def prepare_command(self, request: str) -> dict:
        self.requests.append(request)
        return {
            "status": "infeasible",
            "command": "",
            "explanation": (
                "Synthesized command failed validation and could not be repaired."
            ),
            "confidence": "low",
            "missing_info": [],
            "alternatives": [],
        }


class _BrokenMLXSynthesisSession:
    def prepare_command(self, request: str) -> dict:
        raise RuntimeError(
            "MLX runtime requires Apple Silicon/Metal and mlx-lm. Install in "
            "an MLX-compatible environment with: pip install 'mlx-lm==0.31.3'"
        )


def test_slash_palette_lists_and_filters_commands(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            palette = app.query_one(SlashCommandPalette)

            composer.load_text("/")
            await pilot.pause()
            assert "hidden" not in palette.classes
            commands = [command for command, _ in palette.matches]
            assert "/help" in commands
            assert "/doctor" in commands

            composer.load_text("/d")
            await pilot.pause()
            commands = [command for command, _ in palette.matches]
            assert "/doctor" in commands
            assert "/dryrun" in commands
            assert "/deny" in commands
            assert "/help" not in commands

            composer.load_text("/doctor ")
            await pilot.pause()
            assert "hidden" in palette.classes

    asyncio.run(scenario())


def test_local_provider_can_use_tui_synthesis_mode(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.load_active_provider_config",
        lambda: _local_config(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.SynthesisSession",
        _FakeSynthesisSession,
    )
    _FakeSynthesisSession.requests = []

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            # Local providers default to ask mode (run/harness is incompatible).
            assert app.chat_screen._interaction_mode == "ask"
            assert "ask:local" in str(app.query_one(FooterWidget).renderable)

            composer = app.query_one(Composer)
            composer.load_text("single-point on examples/h2o.xyz with Gaussian")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if app.query_one(FooterWidget).phase == Phase.FINISHED:
                    break

            assert _FakeSynthesisSession.requests == [
                "single-point on examples/h2o.xyz with Gaussian"
            ]
            cells = list(app.query_one(Transcript).query_one("#cells").children)
            assert [type(cell) for cell in cells] == [
                UserMessageCell,
                AgentMessageCell,
                CommandInterpretationCell,
            ]
            synthesis = cells[1]
            assert synthesis.border_title == "Synthesis"
            assert "chemsmart run gaussian" in synthesis.source_text
            assert "runtime semantic gate:" in synthesis.source_text
            assert "verdict: `ok`" in synthesis.source_text
            interpretation = cells[2]
            assert interpretation.border_title == "Command Interpretation"
            assert "deterministic command parser:" in interpretation.source_text
            assert "- program: `gaussian`" in interpretation.source_text
            assert "- job: `sp` (single-point energy)" in interpretation.source_text
            assert "- `-p` meaning: program-level -p/--project" in interpretation.source_text
            assert interpretation.source_text.splitlines()[-1].startswith("Summary:")

    asyncio.run(scenario())


def test_semantic_reject_is_rendered_for_infeasible_synthesis(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.load_active_provider_config",
        lambda: _local_config(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.SynthesisSession",
        _RejectingSynthesisSession,
    )
    _RejectingSynthesisSession.requests = []

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("make a bad command")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if app.query_one(FooterWidget).phase == Phase.FINISHED:
                    break

            cells = list(app.query_one(Transcript).query_one("#cells").children)
            synthesis = next(
                cell
                for cell in cells
                if isinstance(cell, AgentMessageCell)
                and cell.border_title == "Synthesis"
            )
            assert "Synthesized command failed validation" in synthesis.source_text
            assert "runtime semantic gate:" in synthesis.source_text
            assert "verdict: `reject`" in synthesis.source_text
            assert "cmd.semantic.generated_input_missing" in synthesis.source_text

    asyncio.run(scenario())


def test_mlx_runtime_error_is_reported_inside_synthesis_mode(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.load_active_provider_config",
        lambda: _local_config(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.SynthesisSession",
        _BrokenMLXSynthesisSession,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._handle_slash_command("/mode ask")
            await pilot.pause()

            composer = app.query_one(Composer)
            composer.load_text("single-point on examples/h2o.xyz with Gaussian")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if app.query_one(FooterWidget).phase == Phase.FINISHED:
                    break

            cells = list(app.query_one(Transcript).query_one("#cells").children)
            assert not any(isinstance(cell, ErrorCell) for cell in cells)
            synthesis = next(
                cell
                for cell in cells
                if isinstance(cell, AgentMessageCell)
                and cell.border_title == "Synthesis"
            )
            assert "MLX-enabled interpreter" in synthesis.source_text
            assert "mlx-lm==0.31.3" in synthesis.source_text

    asyncio.run(scenario())


def test_local_provider_refuses_run_mode(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.load_active_provider_config",
        lambda: _local_config(),
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            # Local defaults to ask and cannot switch into the run/harness path
            # (the tool-calling loop does not support the local provider).
            assert app.chat_screen._interaction_mode == "ask"
            assert "ask:local" in str(app.query_one(FooterWidget).renderable)

            app.chat_screen._handle_slash_command("/mode run")
            await pilot.pause()
            assert app.chat_screen._interaction_mode == "ask"
            assert "ask:local" in str(app.query_one(FooterWidget).renderable)

    asyncio.run(scenario())
