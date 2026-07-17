"""Generated-input route invariants beyond gaussian.ts + evaluate_harness dispatch."""

from __future__ import annotations

from chemsmart.agent.harness.invariants.route_checks import (
    check_freq_route,
    check_irc_route,
    check_orca_ts_route,
)
from chemsmart.agent.harness.runner import evaluate_harness


def _plan(kind: str) -> dict:
    return {
        "steps": [
            {"tool": "build_molecule", "args": {}},
            {"tool": "build_job", "args": {"kind": kind}},
            {"tool": "dry_run_input", "args": {"job": "$step2"}},
        ]
    }


def _run(kind: str, content: str):
    result = evaluate_harness(
        _plan(kind),
        [
            {
                "content": content,
                "inputfile": "x",
                "command": "chemsmart run orca -f x.xyz -c 0 -m 1 opt",
                "cli_grounded": True,
            }
        ],
    )
    return result.verdict, result.failed_rule_ids


# --- unit: individual invariants -------------------------------------------


def test_orca_freq_missing_freq_rejects():
    r = check_freq_route("!  hf sto-3g defgrid2", software="orca")
    assert r.verdict == "reject"


def test_orca_freq_present_ok():
    assert (
        check_freq_route("! Opt Freq hf sto-3g", software="orca").verdict
        == "ok"
    )


def test_freq_duplicate_rejects():
    assert (
        check_freq_route("! Opt Freq Freq hf", software="orca").verdict
        == "reject"
    )


def test_orca_ts_missing_optts_rejects():
    assert check_orca_ts_route("! Opt hf sto-3g").verdict == "reject"


def test_orca_ts_present_ok():
    assert check_orca_ts_route("! OptTS Freq hf").verdict == "ok"


def test_gaussian_irc_missing_keyword_rejects():
    assert (
        check_irc_route("#p opt b3lyp 6-31g*", software="gaussian").verdict
        == "reject"
    )


def test_gaussian_irc_present_ok():
    assert (
        check_irc_route(
            "#p irc=(calcfc,forward) b3lyp", software="gaussian"
        ).verdict
        == "ok"
    )


# --- integration: evaluate_harness now dispatches per kind ------------------


def test_evaluate_harness_dispatches_orca_freq():
    # regression for the silent ORCA freq bug: a freq job whose route lacks Freq
    verdict, failed = _run(
        "orca.freq", "!  hf sto-3g defgrid2\n* xyz 0 1\n*\n"
    )
    assert verdict == "reject"
    assert "orca.freq.route" in failed


def test_evaluate_harness_orca_freq_ok_when_route_has_freq():
    verdict, failed = _run("orca.freq", "! Opt Freq hf sto-3g\n")
    assert verdict == "ok"
    assert failed == []


def test_evaluate_harness_orca_ts_missing_optts():
    verdict, failed = _run("orca.ts", "! Opt hf sto-3g\n")
    assert verdict == "reject"
    assert "orca.ts.route" in failed


def test_evaluate_harness_gaussian_irc_missing_keyword():
    verdict, failed = _run("gaussian.irc", "#p opt b3lyp 6-31g*\n")
    assert verdict == "reject"
    assert "gaussian.irc.route" in failed


def test_evaluate_harness_no_invariant_kind_is_ok():
    # a kind without a generated-input invariant must not spuriously fail
    verdict, failed = _run("orca.opt", "! Opt hf sto-3g\n")
    assert verdict == "ok"
    assert failed == []
