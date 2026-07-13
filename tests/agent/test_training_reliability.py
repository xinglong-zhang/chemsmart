from __future__ import annotations

from scripts.training.evaluate_agent_reliability import evaluate


def _terminal(*, passed: bool):
    status = "passed" if passed else "failed"
    assertions = [
        {
            "id": "command.returncode",
            "expected": 0,
            "observed": 0 if passed else 1,
            "status": "pass" if passed else "fail",
        }
    ]
    return {
        "schema_version": 1,
        "action": "submit_job",
        "command": "chemsmart sub -s mock-pbs gaussian -p mock opt",
        "status": status,
        "all_passed": passed,
        "returncode": 0 if passed else 1,
        "expected_returncode": None,
        "assertions": assertions,
        "artifacts": [],
    }


def test_sub_terminal_requires_valid_terminal_receipt_and_reports_pass_power():
    rows = [
        {
            "batch_id": "repeat-a",
            "harness_id": "fixed-harness",
            "scenario": "gaussian_opt",
            "session_id": "s1",
            "grade": "PASS_SUB_TERMINAL",
            "terminal_state": _terminal(passed=True),
        },
        {
            "batch_id": "repeat-b",
            "harness_id": "fixed-harness",
            "scenario": "gaussian_opt",
            "session_id": "s2",
            "grade": "PASS_SUB_TERMINAL",
            "terminal_state": _terminal(passed=True),
        },
        {
            "batch_id": "repeat-c",
            "harness_id": "fixed-harness",
            "scenario": "gaussian_opt",
            "session_id": "s3",
            "grade": "PASS_SUB_TERMINAL",
            "terminal_state": _terminal(passed=False),
        },
    ]

    report = evaluate(rows, source="sub_terminal", k=3)
    group = report["groups"][0]
    assert group["trials"] == 3
    assert group["positive_trials"] == 2
    assert group["pass_at_1"] == round(2 / 3, 6)
    assert group["pass_power_k"] == round((2 / 3) ** 3, 6)
    assert group["all_trials_pass"] is False


def test_qmmm_requires_command_and_generated_input_evidence():
    rows = [
        {
            "batch_id": "repeat-a",
            "harness_id": "fixed-harness",
            "scenario": "orca_qmmm",
            "session_id": "s1",
            "grade": "PASS_QMMM",
            "semantic_verdict": "ok",
            "command": "chemsmart run orca opt qmmm",
            "generated_input_evidence": [{"route": "! Opt QMMM"}],
        },
        {
            "batch_id": "repeat-b",
            "harness_id": "fixed-harness",
            "scenario": "orca_qmmm",
            "session_id": "s2",
            "grade": "PASS_QMMM",
            "semantic_verdict": "ok",
            "command": "chemsmart run orca opt qmmm",
            "generated_input_evidence": [],
        },
    ]

    report = evaluate(rows, source="qmmm", k=3)
    group = report["groups"][0]
    assert group["positive_trials"] == 1
    assert group["status"] == "insufficient_trials"
    assert group["pass_power_k"] is None


def test_sub_can_require_request_intent_assertions():
    row = {
        "batch_id": "repeat-a",
        "harness_id": "intent-harness",
        "scenario": "gaussian_opt",
        "session_id": "s1",
        "grade": "PASS_SUB_TERMINAL",
        "terminal_state": _terminal(passed=True),
        "intent_assertions": [
            {"id": "sub.intent_filename", "status": "pass"},
        ],
    }
    assert evaluate(
        [row], source="sub_terminal", k=1, require_intent=True
    )["groups"][0]["positive_trials"] == 1
    row["intent_assertions"][0]["status"] = "fail"
    assert evaluate(
        [row], source="sub_terminal", k=1, require_intent=True
    )["groups"][0]["positive_trials"] == 0
