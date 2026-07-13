"""Regression tests for the live-harness hardening (H1-H4).

Each test pins one verified live-harness defect that let a *silent-wrong*
command (syntactically valid, but quietly changing the intended job) reach the
user:

* H1 - the synthesis prompt falsely denied the real ``gaussian scan``
  subcommand and reused ``-n`` for cores/steps/states.
* H2 - the schema pruner hid the correct sibling subcommand, forcing a
  collapse to the ``{opt, sp}`` floor.
* H3 - the direct ``SynthesisSession.synthesize`` path never applied the intent
  gate, so a coordinate scan collapsed to a plain ``opt`` surfaced as ready.
* H4 - the intent oracle returned ``None`` for common scan/freeze phrasings, so
  the kind assertion never fired.
"""

from __future__ import annotations

import json
from typing import Any

from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
from chemsmart.agent.harness.intent import (
    IntentSpec,
    _kind_from_request,
    evaluate_intent,
)
from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt
from chemsmart.agent.schema_prune import prune_schema_for_request
from chemsmart.agent.synthesis import SynthesisSession


# --------------------------------------------------------------------------- #
# H1 - prompt is grounded in the true CLI
# --------------------------------------------------------------------------- #
def test_prompt_teaches_gaussian_scan_and_drops_false_claim() -> None:
    prompt = build_synthesis_system_prompt(build_chemsmart_cli_schema())
    assert "gaussian scan" in prompt
    # the old, factually-wrong claim must be gone
    assert "There is NO" not in prompt
    # the -n triple-overload example must use the long form for state count
    assert "td --nstates 5" in prompt
    assert "td -n 5" not in prompt


# --------------------------------------------------------------------------- #
# H2 - the pruner keeps the confusable sibling so the model can pick it
# --------------------------------------------------------------------------- #
def _gaussian_kinds(request: str) -> set[str]:
    pruned = prune_schema_for_request(build_chemsmart_cli_schema(), request)
    run = (pruned.get("subcommands") or {}).get("run") or {}
    gaussian = (run.get("subcommands") or {}).get("gaussian") or {}
    return set((gaussian.get("subcommands") or {}).keys())


def test_pruner_keeps_scan_and_modred_together() -> None:
    # A coordinate-scan request must expose BOTH scan and modred; either is a
    # legal target and the disambiguator/model may switch between them.
    kinds = _gaussian_kinds(
        "scan the bond between atoms 1 and 2 of oxetane.xyz with gaussian"
    )
    assert {"scan", "modred"} <= kinds
    # A freeze request likewise keeps both.
    kinds = _gaussian_kinds(
        "freeze bond 1 2 and optimize oxetane.xyz with gaussian"
    )
    assert {"scan", "modred"} <= kinds


def test_pruner_keeps_dias_and_wbi_together() -> None:
    kinds = _gaussian_kinds("distortion interaction analysis on gaussian")
    assert {"dias", "wbi"} <= kinds


# --------------------------------------------------------------------------- #
# H4 - the intent oracle catches the collapse without false-rejecting
# --------------------------------------------------------------------------- #
def _verdict(request: str, command: str) -> str:
    return evaluate_intent(command, IntentSpec.from_request(request)).verdict


def test_intent_rejects_scan_collapsed_to_opt() -> None:
    assert (
        _verdict(
            "scan the bond between atoms 1 and 2 of oxetane.xyz with gaussian",
            "chemsmart run gaussian -f oxetane.xyz opt",
        )
        == "reject"
    )
    assert (
        _verdict(
            "relaxed PES scan of the C-C bond in ethane with gaussian",
            "chemsmart run gaussian -p ethane opt",
        )
        == "reject"
    )


def test_intent_rejects_freeze_collapsed_to_opt() -> None:
    assert (
        _verdict(
            "hold the 1-2 bond fixed and optimize oxetane.xyz gaussian",
            "chemsmart run gaussian -f oxetane.xyz opt",
        )
        == "reject"
    )


def test_intent_accepts_correct_scan_and_rejects_modred_substitution() -> None:
    # correct scan
    assert (
        _verdict(
            "relaxed PES scan of the C-C bond in ethane with gaussian",
            'chemsmart run gaussian -p ethane scan -c "1 2" '
            "--step-size 0.1 --num-steps 10",
        )
        == "ok"
    )
    # Candidate pruning exposes both siblings, but the final intent gate must
    # preserve the user's choice to vary rather than freeze a coordinate.
    assert (
        _verdict(
            "scan the bond 1 2 gaussian",
            'chemsmart run gaussian -f a.xyz modred -c "1 2"',
        )
        == "reject"
    )


def test_intent_does_not_false_reject_tddft_alias() -> None:
    # expected kind is gaussian.tddft; the CLI subcommand is `td`. The oracle
    # must canonicalize so a correct command is not flagged.
    assert (
        _verdict(
            "TDDFT 5 excited states of pyridine with gaussian",
            "chemsmart run gaussian -p pyridine td --nstates 5",
        )
        == "ok"
    )


def test_intent_freeze_does_not_steal_transition_state() -> None:
    # "...with the bond frozen" inside a TS request must stay a TS, not modred.
    assert (
        _verdict(
            "transition state search with the forming bond frozen, gaussian",
            "chemsmart run gaussian -f a.xyz ts",
        )
        == "ok"
    )


def test_kind_from_request_recovers_scan_and_freeze_phrasings() -> None:
    assert (
        _kind_from_request(
            "scan the bond between atoms 1 and 2 with gaussian", "gaussian"
        )
        == "gaussian.scan"
    )
    assert (
        _kind_from_request(
            "relaxed pes scan of the c-c bond in ethane", "gaussian"
        )
        == "gaussian.scan"
    )
    assert (
        _kind_from_request(
            "hold the 1-2 bond fixed and optimize", "gaussian"
        )
        == "gaussian.modred"
    )


# --------------------------------------------------------------------------- #
# H3 - the direct synthesize path enforces intent end-to-end
# --------------------------------------------------------------------------- #
class _MockProvider:
    name = "openai"

    def __init__(self, command: str) -> None:
        self._command = command

    def chat(self, messages: list[dict[str, Any]], **_kwargs: Any) -> dict:
        return {
            "content": json.dumps(
                {
                    "status": "ready",
                    "command": self._command,
                    "explanation": "draft",
                    "confidence": "high",
                    "missing_info": [],
                    "alternatives": [],
                }
            )
        }


def _synthesize(request: str, command: str) -> dict:
    session = SynthesisSession(
        provider=_MockProvider(command), default_project="demo"
    )
    return session.synthesize(request)


def test_synthesize_downgrades_silent_wrong_command() -> None:
    result = _synthesize(
        "scan the bond 1 2 of oxetane.xyz with gaussian",
        "chemsmart run gaussian -f oxetane.xyz opt",
    )
    assert result["status"] == "needs_clarification"
    assert result.get("intent_reject") == ["intent.kind"]


def test_synthesize_keeps_correct_command_ready() -> None:
    result = _synthesize(
        "scan the bond 1 2 of oxetane.xyz with gaussian",
        'chemsmart run gaussian -f oxetane.xyz scan -c "1 2" '
        "--step-size 0.1 --num-steps 10",
    )
    assert result["status"] == "ready"
    assert "intent_reject" not in result
