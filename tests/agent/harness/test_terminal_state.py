from chemsmart.agent.harness.terminal_state import (
    assertion,
    build_terminal_state,
    terminal_state_is_positive,
    validate_terminal_state,
)


def test_terminal_state_requires_all_observable_assertions_to_pass():
    state = build_terminal_state(
        action="submit_job",
        command="chemsmart sub -s mock-pbs gaussian ...",
        returncode=0,
        assertions=[
            assertion("server.yaml_exists", expected=True, observed=True),
            assertion("scheduler.returncode", expected=0, observed=0),
        ],
    )

    assert terminal_state_is_positive(state)
    assert state["status"] == "passed"
    assert validate_terminal_state(state) == []


def test_terminal_state_rejects_nonzero_returncode_and_failed_assertion():
    state = build_terminal_state(
        action="submit_job",
        command="chemsmart sub -s mock-pbs gaussian ...",
        returncode=1,
        assertions=[
            assertion("server.yaml_exists", expected=True, observed=True),
            assertion("scheduler.returncode", expected=0, observed=1),
        ],
    )

    assert not terminal_state_is_positive(state)
    issues = validate_terminal_state(state)
    assert "terminal_state.failed:scheduler.returncode" in issues
    assert "terminal_state.returncode_nonzero" in issues


def test_terminal_state_accepts_expected_negative_terminal_outcome():
    state = build_terminal_state(
        action="repair_submit_job",
        command="chemsmart sub -s missing-pbs gaussian ...",
        returncode=1,
        expected_returncode=1,
        assertions=[
            assertion("server.yaml_exists", expected=True, observed=True),
            assertion("scheduler.rejection_returncode", expected=1, observed=1),
        ],
    )

    assert validate_terminal_state(state) == []
    assert state["status"] == "passed"
    assert state["all_passed"] is True
    assert not terminal_state_is_positive(state)
