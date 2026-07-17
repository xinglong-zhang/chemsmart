"""Tests for natural-language ChemSmart CLI command synthesis."""

from __future__ import annotations

import json
from typing import Any

import pytest
from click.testing import CliRunner

from chemsmart.agent.harness.command_semantics import (
    CommandSemanticIssue,
    CommandSemanticResult,
)
from chemsmart.agent.harness.workflow_state import current_workflow_state
from chemsmart.agent.cli import agent
from chemsmart.agent.synthesis import SynthesisSession, resolve_default_project

READY = {
    "status": "ready",
    "command": "chemsmart sub gaussian ts -p oxetane -b def2-svp",
    "explanation": "Prepare a Gaussian transition-state job.",
    "confidence": "high",
    "missing_info": [],
    "alternatives": [],
}


class DummyProvider:
    """Deterministic provider double returning queued responses."""

    def __init__(self, responses: list[Any], *, name: str = "local") -> None:
        self.name = name
        self.responses = list(responses)
        self.messages: list[list[dict[str, str]]] = []

    def chat(self, messages: list[dict[str, str]], **_kwargs: Any) -> Any:
        self.messages.append(messages)
        if not self.responses:
            raise AssertionError("provider response queue exhausted")
        return self.responses.pop(0)


def _json_response(payload: dict[str, Any]) -> str:
    return json.dumps(payload)


def test_synthesize_parses_happy_path_json() -> None:
    provider = DummyProvider([_json_response(READY)])
    session = SynthesisSession(provider=provider, schema={"subcommands": {}})

    result = session.synthesize("make a ts job")

    assert result == READY
    assert provider.messages[0][-1] == {
        "role": "user",
        "content": "make a ts job",
    }


def test_synthesize_injects_default_project_for_api_command() -> None:
    ready_without_project = {
        "status": "ready",
        "command": (
            "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 "
            "-x B3LYP -b '6-31G(d)' opt"
        ),
        "explanation": "Prepare a Gaussian optimization.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    provider = DummyProvider([_json_response(ready_without_project)])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        default_project="test",
    )

    result = session.synthesize("optimize water")

    assert result["project"] == "test"
    assert result["command"].startswith("chemsmart run gaussian -p test ")


def test_synthesize_preserves_explicit_api_project() -> None:
    provider = DummyProvider([_json_response(READY)])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        default_project="test",
    )

    result = session.synthesize("make a ts job")

    assert " -p oxetane " in result["command"]
    assert " -p test " not in result["command"]


def test_frontier_router_explains_previous_command_without_resynthesis() -> (
    None
):
    ready_without_project = {
        "status": "ready",
        "command": (
            "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 "
            "-x B3LYP -b '6-31G(d)' opt"
        ),
        "explanation": "Prepare a Gaussian optimization.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    synthesize_decision = {
        "action": "synthesize_command",
        "target_command": "",
        "confidence": "high",
        "missing_info": [],
        "decision_summary": "The user asked to prepare a job command.",
        "evidence": ["The request says prepare a Gaussian optimization."],
        "rejected_actions": {
            "explain_command": "No existing command explanation was requested."
        },
    }
    explain_decision = {
        "action": "explain_command",
        "target_command": "",
        "confidence": "high",
        "missing_info": [],
        "decision_summary": "The user asked what the previous command does.",
        "evidence": ["The request asks what this command means."],
        "rejected_actions": {
            "synthesize_command": "No new calculation setup was requested."
        },
    }
    # After the router picks explain_command, the provider is called again to
    # COMPOSE the user-facing answer from grounded facts (tool-grounded design).
    composed_answer = {
        "answer": (
            "这个命令让 ChemSmart 调用 gaussian，对 examples/h2o.xyz "
            "在本地运行一次几何优化 (opt)。电荷/多重度是 0/1，项目是 test。"
        ),
        "facts_used": {
            "program": "gaussian",
            "job": "opt",
            "action": "run",
            "project": "test",
            "filename": "examples/h2o.xyz",
            "charge": "0",
            "multiplicity": "1",
        },
        "reasoning": [
            "用户要求用中文解释上一条命令。",
            "命令是 gaussian opt。",
        ],
        "caveats": [],
    }
    provider = DummyProvider(
        [
            _json_response(synthesize_decision),
            _json_response(ready_without_project),
            _json_response(explain_decision),
            _json_response(composed_answer),
        ],
        name="openai",
    )
    session = SynthesisSession(
        provider=provider,
        default_project="test",
        semantic_gate=False,
    )

    first = session.prepare_command(
        "Prepare a Gaussian optimization for examples/h2o.xyz."
    )
    second = session.prepare_command(
        "이 명령이 수행하는 작업은 뭐야? 중국어로 설명해줘."
    )

    assert first["status"] == "ready"
    assert first["command"].startswith("chemsmart run gaussian -p test ")
    assert first["decision_trace"]["action"] == "synthesize_command"
    assert second["status"] == "informational"
    assert second["action"] == "explain_command"
    assert second["command"] == first["command"]
    # The composed (model-authored) Chinese answer is the user-facing text —
    # not a deterministic English field dump.
    assert second["explanation"] == composed_answer["answer"]
    assert "deterministic command parser:" not in second["explanation"]
    assert second["grounded"] is True
    assert second["reasoning"] == composed_answer["reasoning"]
    assert second["decision_trace"]["action"] == "explain_command"
    assert (
        second["decision_trace"]["reasoning"] == composed_answer["reasoning"]
    )
    assert second["decision_trace"]["target_command"] == first["command"]
    assert len(provider.messages) == 4


def test_local_provider_skips_frontier_intent_router() -> None:
    provider = DummyProvider([_json_response(READY)], name="local")
    session = SynthesisSession(
        provider=provider,
        semantic_gate=False,
        enable_intent_router=None,
    )

    result = session.prepare_command("이 명령이 수행하는 작업은 뭐야?")

    assert result["status"] == "ready"
    assert "decision_trace" not in result
    assert len(provider.messages) == 1


def test_resolve_default_project_uses_single_workspace_project(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "test.yaml").write_text("gas: {}\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("CHEMSMART_AGENT_PROJECT", raising=False)

    assert resolve_default_project() == "test"


def test_prepare_command_falls_back_when_workspace_yaml_missing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    monkeypatch.chdir(tmp_path)
    provider = DummyProvider([_json_response(READY)], name="local")
    session = SynthesisSession(
        provider=provider,
        semantic_gate=False,
        enable_intent_router=False,
    )

    result = session.prepare_command("Prepare a Gaussian optimization.")

    assert result["status"] == "needs_clarification"
    assert "No workspace project YAML" in result["explanation"]
    assert result["yaml_check"]["loaded"] is False
    assert provider.messages == []


def test_prepare_command_uses_unique_workspace_project_and_records_state(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "demo.yaml").write_text("gas: {}\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    ready_without_project = {
        "status": "ready",
        "command": (
            "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 opt"
        ),
        "explanation": "Prepare a Gaussian optimization.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    provider = DummyProvider(
        [_json_response(ready_without_project)], name="local"
    )
    session = SynthesisSession(
        provider=provider,
        semantic_gate=False,
        enable_intent_router=False,
    )

    result = session.prepare_command("Prepare a Gaussian optimization.")

    assert result["status"] == "ready"
    assert result["command"].startswith("chemsmart run gaussian -p demo ")
    assert session.default_project == "demo"
    assert current_workflow_state().project is not None
    assert current_workflow_state().project.name == "demo"
    assert current_workflow_state().project.program == "gaussian"


def test_synthesize_retries_once_after_invalid_json() -> None:
    provider = DummyProvider(["not json", _json_response(READY)])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        max_retries=1,
    )

    result = session.synthesize("make a ts job")

    assert result["status"] == "ready"
    assert len(provider.messages) == 2
    assert (
        "previous response was invalid" in provider.messages[1][-1]["content"]
    )


def test_synthesize_raises_after_two_invalid_json_responses() -> None:
    provider = DummyProvider(["not json", "still not json"])
    session = SynthesisSession(
        provider=provider,
        schema={"subcommands": {}},
        max_retries=1,
    )

    with pytest.raises(ValueError, match="failed to return valid"):
        session.synthesize("make a ts job")

    assert len(provider.messages) == 2


def test_validate_command_accepts_legal_sub_gaussian_ts_command() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command(
        "chemsmart sub gaussian ts -p oxetane -b def2-svp"
    )

    assert ok is True
    assert error == ""


def test_validate_command_rejects_invented_option() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command("chemsmart sub gaussian ts -z bogus")

    assert ok is False
    assert "-z" in error


def test_validate_command_rejects_non_chemsmart_command() -> None:
    session = SynthesisSession(provider=DummyProvider([]))

    ok, error = session.validate_command("python -m chemsmart")

    assert ok is False
    assert "chemsmart" in error


def test_run_interactive_retries_after_runtime_semantic_reject(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    bad = {
        "status": "ready",
        "command": "chemsmart run orca opt -f water.xyz -c 0 -m 1",
        "explanation": "Bad order.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    repaired = {
        "status": "ready",
        "command": "chemsmart run orca -f water.xyz -c 0 -m 1 opt",
        "explanation": "Repaired.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    provider = DummyProvider([_json_response(bad), _json_response(repaired)])
    session = SynthesisSession(provider=provider, max_retries=1)
    semantic_calls: list[str] = []
    confirmed: list[str] = []

    def fake_semantic(command: str, **_kwargs: Any) -> CommandSemanticResult:
        semantic_calls.append(command)
        if len(semantic_calls) == 1:
            return CommandSemanticResult(
                verdict="reject",
                command=command,
                checked_argv=("chemsmart", "run", "orca", "opt"),
                issues=(
                    CommandSemanticIssue(
                        rule_id="cmd.semantic.safe_execution_failed",
                        severity="reject",
                        message="No server implemented for local.yaml.",
                        missing_info=(
                            "valid chemsmart server configuration",
                            "available servers: ['colab_orca']",
                        ),
                    ),
                ),
                stderr_tail=(
                    "No server implemented for local.yaml.\n"
                    "Currently available servers: ['colab_orca']"
                ),
            )
        return CommandSemanticResult(verdict="ok", command=command)

    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        fake_semantic,
    )
    monkeypatch.setattr(
        session,
        "confirm_and_execute",
        lambda command, *_args: confirmed.append(command),
    )

    session.run_interactive("run ORCA water opt")

    assert confirmed == ["chemsmart run orca -f water.xyz -c 0 -m 1 opt"]
    assert (
        "Notice: synthesized command failed runtime semantic validation"
        in (capsys.readouterr().out)
    )
    retry_prompt = provider.messages[1][-1]["content"]
    assert "runtime semantic validation" in retry_prompt
    assert "valid chemsmart server configuration" in retry_prompt
    assert "available servers: ['colab_orca']" in retry_prompt


def test_prepare_command_repairs_missing_job_before_intent_gate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    bad = {
        "status": "ready",
        "command": (
            "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 "
            "-r freq"
        ),
        "explanation": "Optimization and frequency.",
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }
    repaired = {
        **bad,
        "command": bad["command"] + " opt",
        "explanation": "Repaired optimization and frequency command.",
    }
    provider = DummyProvider(
        [_json_response(bad), _json_response(repaired)],
        name="openai",
    )
    session = SynthesisSession(
        provider=provider,
        default_project="demo",
        max_retries=1,
        enable_intent_router=False,
    )
    semantic_calls: list[str] = []

    def fake_semantic(command: str, **_kwargs: Any) -> CommandSemanticResult:
        semantic_calls.append(command)
        if not command.endswith(" opt"):
            return CommandSemanticResult(
                verdict="reject",
                command=command,
                issues=(
                    CommandSemanticIssue(
                        rule_id="cmd.contract.job_subcommand_required",
                        severity="reject",
                        message="explicit Gaussian job subcommand required",
                    ),
                ),
            )
        return CommandSemanticResult(verdict="ok", command=command)

    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        fake_semantic,
    )

    result = session.prepare_command(
        "Prepare a Gaussian optimization and frequency calculation for "
        "water.xyz, charge 0 multiplicity 1."
    )

    assert semantic_calls == [bad["command"], repaired["command"]]
    assert result["status"] == "ready"
    assert result["command"].endswith(" opt")
    assert result["intent_assertion"]["verdict"] == "ok"


def test_prepare_command_repairs_explicit_intent_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    command_12 = (
        "chemsmart run orca -p orca_demo -f complex.xyz -c 0 -m 1 scan "
        "--coordinates '[2,7]' --dist-start 1.8 --dist-end 3.0 "
        "--num-steps 12"
    )
    command_13 = command_12.replace("--num-steps 12", "--num-steps 13")

    def response(command: str) -> str:
        return _json_response(
            {
                "status": "ready",
                "command": command,
                "explanation": "Prepare an ORCA scan.",
                "confidence": "high",
                "missing_info": [],
                "alternatives": [],
            }
        )

    provider = DummyProvider(
        [response(command_12), response(command_13)],
        name="openai",
    )
    session = SynthesisSession(
        provider=provider,
        default_project="orca_demo",
        max_retries=1,
        enable_intent_router=False,
    )
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        lambda command, **_kwargs: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )

    result = session.prepare_command(
        "Using ORCA project orca_demo, scan bond 2-7 in complex.xyz from "
        "1.8 to 3.0 angstrom with 13 points."
    )

    assert result["status"] == "ready"
    assert result["command"] == command_13
    assert result["intent_assertion"]["verdict"] == "ok"
    retry_prompt = provider.messages[1][-1]["content"]
    assert "intent.chemistry.num_steps" in retry_prompt
    assert "expected 13" in retry_prompt
    assert "observed '12'" in retry_prompt


def test_run_interactive_reports_final_runtime_semantic_reject(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    provider = DummyProvider([_json_response(READY)])
    session = SynthesisSession(provider=provider, max_retries=0)

    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        lambda command, **_kwargs: CommandSemanticResult(
            verdict="reject",
            command=command,
            issues=(
                CommandSemanticIssue(
                    rule_id="cmd.semantic.generated_input_missing",
                    severity="reject",
                    message="no input generated",
                ),
            ),
        ),
    )

    session.run_interactive("make a ts job")

    output = capsys.readouterr().out
    assert (
        "Notice: synthesized command failed runtime semantic validation"
        in output
    )
    assert "no input generated" in output


def test_confirm_and_execute_inserts_test_only_for_sub(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    calls: list[list[str]] = []
    choices = iter(["T", "T"])

    monkeypatch.setattr(
        "click.prompt", lambda *_args, **_kwargs: next(choices)
    )
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian ts -p oxetane",
        "explain",
        "medium",
    )
    session.confirm_and_execute("chemsmart config show", "explain", "medium")

    assert calls == [
        [
            "chemsmart",
            "--no-verbose",
            "sub",
            "--test",
            "gaussian",
            "ts",
            "-p",
            "oxetane",
        ]
    ]


def test_confirm_and_execute_injects_no_verbose_on_yes_by_default(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    calls: list[list[str]] = []

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "Y")
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian opt -p water -b 6-31g*",
        "explain",
        "high",
    )

    assert calls == [
        [
            "chemsmart",
            "--no-verbose",
            "sub",
            "gaussian",
            "opt",
            "-p",
            "water",
            "-b",
            "6-31g*",
        ]
    ]


def test_confirm_and_execute_skips_no_verbose_when_debug(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    session = SynthesisSession(provider=DummyProvider([]), debug=True)
    calls: list[list[str]] = []

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "Y")
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )

    session.confirm_and_execute(
        "chemsmart sub gaussian opt -p water -b 6-31g*",
        "explain",
        "high",
    )

    assert calls == [
        [
            "chemsmart",
            "sub",
            "gaussian",
            "opt",
            "-p",
            "water",
            "-b",
            "6-31g*",
        ]
    ]


def test_quiet_argv_respects_existing_verbose_flags() -> None:
    session = SynthesisSession(provider=DummyProvider([]))
    argv_with_verbose = session._quiet_argv(
        ["chemsmart", "--verbose", "sub", "gaussian", "opt"]
    )
    argv_with_no_verbose = session._quiet_argv(
        ["chemsmart", "--no-verbose", "sub", "gaussian", "opt"]
    )

    assert argv_with_verbose == [
        "chemsmart",
        "--verbose",
        "sub",
        "gaussian",
        "opt",
    ]
    assert argv_with_no_verbose == [
        "chemsmart",
        "--no-verbose",
        "sub",
        "gaussian",
        "opt",
    ]


def test_multiturn_clarification_is_carried_in_memory(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clarification = {
        "status": "needs_clarification",
        "command": "",
        "explanation": "Need method.",
        "confidence": "medium",
        "missing_info": ["What basis set?"],
        "alternatives": [],
    }
    provider = DummyProvider(
        [_json_response(clarification), _json_response(READY)]
    )
    session = SynthesisSession(provider=provider, semantic_gate=False)

    monkeypatch.setattr("click.prompt", lambda *_args, **_kwargs: "def2-svp")
    monkeypatch.setattr(session, "confirm_and_execute", lambda *_args: None)

    session.run_interactive("optimize oxetane")

    assert len(provider.messages) == 2
    second_call = provider.messages[1]
    assert [message["role"] for message in second_call] == ["system", "user"]
    assert "What basis set?: def2-svp" in second_call[-1]["content"]


def test_agent_ask_e2e_uses_synthesis_session(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    provider = DummyProvider([_json_response(READY)])
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider",
        lambda: provider,
    )
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        lambda command, **_kwargs: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )
    monkeypatch.setattr(
        "chemsmart.agent.cli_commands._agent_log_root",
        lambda: tmp_path,
    )

    result = CliRunner().invoke(agent, ["ask", "make a ts job"], input="N\n")

    assert result.exit_code == 0
    assert "chemsmart sub gaussian ts -p oxetane -b def2-svp" in result.output
    assert "Run? [Y]es/[E]dit/[N]o/[T]est" in result.output
    assert calls == []


def test_agent_ask_debug_flag_plumbs_to_synthesis_session(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    provider = DummyProvider([_json_response(READY)])
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider",
        lambda: provider,
    )
    captured: dict[str, bool] = {}
    real_init = SynthesisSession.__init__

    def spy_init(self, *args: Any, **kwargs: Any) -> None:
        captured["debug"] = bool(kwargs.get("debug", False))
        real_init(self, *args, **kwargs)

    monkeypatch.setattr(SynthesisSession, "__init__", spy_init)
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.subprocess.run",
        lambda args, check=False: calls.append(args),
    )
    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        lambda command, **_kwargs: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )

    result = CliRunner().invoke(
        agent, ["ask", "--debug", "make a ts job"], input="Y\n"
    )

    assert result.exit_code == 0, result.output
    assert captured.get("debug") is True
    assert calls == [
        [
            "chemsmart",
            "sub",
            "gaussian",
            "ts",
            "-p",
            "oxetane",
            "-b",
            "def2-svp",
        ]
    ]


def test_clarification_from_semantic_asks_for_missing_charge_multiplicity():
    from chemsmart.agent.harness.command_semantics import (
        CommandSemanticResult,
    )
    from chemsmart.agent.synthesis import _clarification_from_semantic

    assert _clarification_from_semantic(None) == []

    reject = CommandSemanticResult(
        verdict="reject",
        command="chemsmart run gaussian -p co2 -f examples/h2o.xyz opt",
        stderr_tail=(
            "ValueError: Charge and multiplicity must be set for "
            "Gaussian/ORCA jobs."
        ),
    )
    questions = _clarification_from_semantic(reject)
    assert len(questions) == 1
    assert "charge" in questions[0].lower()
    assert "multiplicity" in questions[0].lower()

    # An unrelated runtime failure yields no synthetic clarification.
    other = CommandSemanticResult(
        verdict="reject",
        command="chemsmart run gaussian -p co2 -f examples/h2o.xyz opt",
        stderr_tail="RuntimeError: some other failure",
    )
    assert _clarification_from_semantic(other) == []


def test_synthesize_tolerates_think_preamble_and_captures_reasoning() -> None:
    # Reasoning models (Qwen3-Thinking / DeepSeek) emit <think>…</think>
    # before the JSON; synthesis must parse the JSON AND keep the reasoning.
    reasoning = "The user wants a TS job; gaussian ts subcommand, def2-svp."
    content = f"<think>{reasoning}</think>\n" + json.dumps(READY)
    provider = DummyProvider([content])
    session = SynthesisSession(provider=provider, schema={"subcommands": {}})

    result = session.synthesize("make a ts job")

    assert result == READY
    assert session._last_reasoning == reasoning


def test_synthesize_captures_reasoning_content_field() -> None:
    # DashScope/OpenAI-reasoning convention: reasoning in a separate field,
    # content is the clean JSON.
    reasoning = "Single point, gaussian sp, neutral singlet."
    response = {
        "choices": [
            {
                "message": {
                    "role": "assistant",
                    "content": json.dumps(READY),
                    "reasoning_content": reasoning,
                }
            }
        ]
    }
    provider = DummyProvider([response])
    session = SynthesisSession(provider=provider, schema={"subcommands": {}})

    result = session.synthesize("sp please")

    assert result == READY
    assert session._last_reasoning == reasoning


def test_synthesize_reasoning_empty_when_absent() -> None:
    provider = DummyProvider([_json_response(READY)])
    session = SynthesisSession(provider=provider, schema={"subcommands": {}})
    session.synthesize("make a ts job")
    assert session._last_reasoning == ""
