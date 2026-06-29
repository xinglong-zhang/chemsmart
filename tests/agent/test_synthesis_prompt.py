"""Locks in the strengthened synthesis system prompt.

These assertions are deliberate guardrails so future refactors do not
accidentally drop the routing rules, engine guidance, subcommand picks,
option-spelling clarifications, or the few-shot example set that drive
synthesis quality.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.prompts.synthesis import build_synthesis_system_prompt

_EMPTY_SCHEMA = {"name": "chemsmart", "subcommands": {}}


@pytest.fixture
def prompt() -> str:
    return build_synthesis_system_prompt(_EMPTY_SCHEMA)


def test_prompt_declares_routing_rule_for_remote_server(prompt: str) -> None:
    assert "chemsmart sub" in prompt
    assert "chemnode1" in prompt
    assert "qsub" in prompt or "sbatch" in prompt


def test_prompt_declares_routing_rule_for_local(prompt: str) -> None:
    assert "chemsmart run" in prompt
    assert "locally" in prompt or "without submission" in prompt


def test_prompt_distinguishes_project_p_from_pubchem_P(prompt: str) -> None:
    assert "-p/--project" in prompt
    assert "-P/--pubchem" in prompt
    assert "from PubChem" in prompt


def test_prompt_prohibits_orca_td(prompt: str) -> None:
    assert "ORCA does NOT have td" in prompt or (
        "ORCA" in prompt and "does NOT have td" in prompt
    )


def test_prompt_prohibits_gaussian_scan_subcommand(prompt: str) -> None:
    assert "NO `gaussian scan`" in prompt or "NO gaussian scan" in prompt


def test_prompt_includes_nci_family_guidance(prompt: str) -> None:
    assert "nciplot" in prompt
    assert "gaussian nci" in prompt
    assert ".wfn" in prompt


def test_prompt_includes_modred_canonical_form(prompt: str) -> None:
    assert "gaussian modred" in prompt
    assert "-c/--coordinates" in prompt


def test_prompt_carries_few_shot_examples_for_each_axis(prompt: str) -> None:
    assert "chemsmart run gaussian -p h2o" in prompt
    assert "chemsmart sub -s chemnode1" in prompt
    assert "chemsmart run nciplot -f dimer.wfn" in prompt
    assert "chemsmart run thermochemistry boltzmann" in prompt
    assert "chemsmart run gaussian -f oxetane.xyz modred" in prompt


def test_prompt_keeps_no_debug_flag_rule(prompt: str) -> None:
    assert "--debug" in prompt and "--verbose" in prompt
    assert "unless the user explicitly asked" in prompt
