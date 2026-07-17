from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.provider_config import AgentProviderConfig
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.events import ToolUseEvent
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    ErrorCell,
    FinalAnswerCell,
    SynthesisTraceCell,
    ToolChainToggleCell,
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


def _openai_config() -> AgentProviderConfig:
    return AgentProviderConfig(
        name="frontier_openai",
        type="openai",
        api_key="test-key",
        model="gpt-test",
        base_url="https://api.example.test/v1",
        extra_headers={},
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
            "project": "test",
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


class _InformationalSynthesisSession:
    requests: list[str] = []

    def __init__(self) -> None:
        self._last_semantic_result = None
        self._last_raw_response = (
            '{"action":"explain_command","decision_summary":"explain previous"}'
        )

    def prepare_command(self, request: str) -> dict:
        self.requests.append(request)
        return {
            "status": "informational",
            "command": (
                "chemsmart run gaussian -p test -f examples/h2o.xyz "
                "-c 0 -m 1 opt"
            ),
            "explanation": "deterministic command parser:\n- program: gaussian",
            "confidence": "high",
            "project": "test",
            "missing_info": [],
            "alternatives": [],
            "action": "explain_command",
            "decision_trace": {
                "router": "api_frontier_intent_router",
                "action": "explain_command",
                "confidence": "high",
                "decision_summary": "The user asked what the command does.",
                "target_command": (
                    "chemsmart run gaussian -p test -f examples/h2o.xyz "
                    "-c 0 -m 1 opt"
                ),
                "default_project": "test",
                "last_command_available": True,
                "request_excerpt": request,
                "evidence": ["The request asks for command meaning."],
                "rejected_actions": {
                    "synthesize_command": "No new job setup was requested."
                },
                "note": "public trace",
            },
        }


class _BrokenMLXSynthesisSession:
    def prepare_command(self, request: str) -> dict:
        raise RuntimeError(
            "MLX runtime requires Apple Silicon/Metal and mlx-lm. Install in "
            "an MLX-compatible environment with: pip install 'mlx-lm==0.31.3'"
        )


class _FakeUnifiedAgentSession:
    requests: list[str] = []
    prompts: list[str] = []

    def __init__(self, *args, **kwargs) -> None:
        self.stage_prompt = kwargs.get("stage_prompt", "tool_loop.md")
        root = Path(kwargs.get("session_root") or ".")
        self.session_dir = root / "unified_fake"
        self.session_dir.mkdir(parents=True, exist_ok=True)
        self.prompts.append(self.stage_prompt)

    def run_loop(self, request: str, **kwargs) -> dict:
        self.requests.append(request)
        policy = kwargs["policy"]
        return {
            "session_id": "unified_fake",
            "session_dir": "",
            "plan": None,
            "plan_text": None,
            "critic_verdict": None,
            "completed_steps": 0,
            "blocked": False,
            "dry_run_result": None,
            "dry_run_results": [],
            "runtime_result": None,
            "preview_submit": None,
            "results": [],
            "assistant_output": "Unified loop handled request.",
            "tool_requests": [],
            "tool_outcomes": [],
            "final_message": "Unified loop handled request.",
            "advisory_only": True,
            "is_chitchat": False,
            "session_allow_tools": sorted(policy.session_allow),
        }


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
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
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
            assert app.chat_screen._interaction_mode == "synthesis"
            assert "synthesis:local" in str(
                app.query_one(FooterWidget).renderable
            )

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
                SynthesisTraceCell,
                CommandInterpretationCell,
                AgentMessageCell,
                ToolChainToggleCell,
                FinalAnswerCell,
            ]
            trace = cells[1]
            assert trace.border_title == "Agent Trace"
            assert "lane: local synthesis" in trace.source_text
            assert "semantic gate: ok; failed rules: none" in trace.source_text
            assert (tmp_path / "sessions" / "local" / "synthesis").is_dir()
            interpretation = cells[2]
            assert (
                interpretation.border_title
                == "Deterministic Harness: CLI Command Synthesis ▾"
            )
            assert "deterministic command parser:" in interpretation.source_text
            assert "- program: `gaussian`" in interpretation.source_text
            assert "- job: `sp` (single-point energy)" in interpretation.source_text
            assert "- `-p` meaning: program-level -p/--project" in interpretation.source_text
            assert interpretation.source_text.splitlines()[-1].startswith("Summary:")
            assert not cells[1].display
            assert not cells[2].display
            assert not cells[3].display
            final = cells[5]
            assert final.border_title == "Final Command"
            assert "chemsmart run gaussian" in final.source_text
            assert "confidence:" not in final.source_text

    asyncio.run(scenario())


def test_frontier_provider_defaults_to_unified_tool_loop(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("CHEMSMART_AGENT_TUI_MODE", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
        lambda: _openai_config(),
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.AgentSession",
        _FakeUnifiedAgentSession,
    )
    _FakeUnifiedAgentSession.requests = []
    _FakeUnifiedAgentSession.prompts = []

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            assert app.chat_screen._interaction_mode == "unified"
            footer_text = str(app.query_one(FooterWidget).renderable)
            assert "unified:openai" in footer_text
            assert "project test" in footer_text

            composer = app.query_one(Composer)
            composer.load_text("single-point on examples/h2o.xyz with Gaussian")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if app.query_one(FooterWidget).phase == Phase.FINISHED:
                    break

            assert _FakeUnifiedAgentSession.requests == [
                "single-point on examples/h2o.xyz with Gaussian"
            ]
            assert _FakeUnifiedAgentSession.prompts == ["unified_agent.md"]
            cells = list(app.query_one(Transcript).query_one("#cells").children)
            assert [type(cell) for cell in cells] == [UserMessageCell]

    asyncio.run(scenario())


def test_command_tool_result_renders_synthesis_trace(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("CHEMSMART_AGENT_TUI_MODE", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
        lambda: _openai_config(),
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_tool_use_event(
                ToolUseEvent(
                    kind="tool_use_result",
                    step_index=1,
                    tool="synthesize_command",
                    status="ok",
                    args={},
                    description="Synthesize chemsmart command",
                    provider_call_id="call_1",
                    payload={
                        "status": "ready",
                        "command": (
                            "chemsmart run gaussian -p test "
                            "-f examples/h2o.xyz -c 0 -m 1 opt"
                        ),
                        "explanation": "Prepared by unified command tool.",
                        "confidence": "high",
                        "project": "test",
                        "semantic": _SemanticOk().to_dict(),
                        "decision_trace": {
                            "router": "unified_agent",
                            "action": "synthesize_command",
                            "confidence": "high",
                            "decision_summary": "The user requested a job.",
                            "evidence": ["The request asks for command setup."],
                        },
                    },
                )
            )
            await pilot.pause()

            cells = list(app.query_one(Transcript).query_one("#cells").children)
            assert any(
                cell.__class__.__name__ == "ToolCallCell" for cell in cells
            )
            rendered = [
                cell
                for cell in cells
                if isinstance(
                    cell,
                    (
                        SynthesisTraceCell,
                        CommandInterpretationCell,
                        FinalAnswerCell,
                    ),
                )
            ]
            assert [type(cell) for cell in rendered] == [
                SynthesisTraceCell,
                CommandInterpretationCell,
                FinalAnswerCell,
            ]
            trace_cell = rendered[0]
            assert trace_cell.border_title == "Agent Trace"
            assert "lane: api synthesis" in trace_cell.source_text
            assert "The request asks for command setup." in trace_cell.source_text
            final = rendered[2]
            assert final.border_title == "Final Command"
            assert "chemsmart run gaussian" in final.source_text

    asyncio.run(scenario())


def test_semantic_reject_is_rendered_for_infeasible_synthesis(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
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
                if isinstance(cell, FinalAnswerCell)
                and cell.border_title == "Final Status"
            )
            assert "Synthesized command failed validation" in synthesis.source_text
            trace = next(cell for cell in cells if isinstance(cell, SynthesisTraceCell))
            assert "semantic gate: reject" in trace.source_text
            assert "cmd.semantic.generated_input_missing" in synthesis.source_text

    asyncio.run(scenario())


def test_mlx_runtime_error_is_reported_inside_synthesis_mode(
    monkeypatch, tmp_path: Path
):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
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
                if isinstance(cell, FinalAnswerCell)
                and cell.border_title == "Final Status"
            )
            assert "MLX-enabled interpreter" in synthesis.source_text
            assert "mlx-lm==0.31.3" in synthesis.source_text

    asyncio.run(scenario())


def test_removed_mode_command_is_rejected(monkeypatch, tmp_path: Path):
    monkeypatch.delenv("AI_PROVIDER", raising=False)
    monkeypatch.setattr(
        "chemsmart.agent.tui.chat_helpers.load_active_provider_config",
        lambda: _local_config(),
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            assert app.chat_screen._interaction_mode == "synthesis"
            assert "synthesis:local" in str(
                app.query_one(FooterWidget).renderable
            )

            app.chat_screen._handle_slash_command("/mode")
            await pilot.pause()
            cells = list(app.query_one(Transcript).query_one("#cells").children)
            error = next(cell for cell in cells if isinstance(cell, ErrorCell))
            assert error.error_title == "Unknown command"
            assert error.message == "/mode"
            assert app.chat_screen._interaction_mode == "synthesis"
            assert "synthesis:local" in str(
                app.query_one(FooterWidget).renderable
            )

    asyncio.run(scenario())
