"""The Agent can only plan in a model it has been told about.

The wave contract is part of the natural-language policy the host puts in
front of the model, and that policy is a registry rather than free text:
`POLICY_RULES` renders the stem prompt, the wake context and the tool
descriptions, and `test_every_rule_renders_once` holds it. So teaching
the Agent to plan in waves is registered capability work in one place.

Three things it must learn, because none is derivable from the tools:
that execution proceeds in waves it selects and it is woken once when
every member is terminal; that scientific width is its own while machine
concurrency is the host's, so asking for seven costs one turn and not
two; and that it never sizes the machine.

A rule with no reader is the defect class this whole round exists to
remove, so the rendered surface is asserted, not just the registry.
"""

from __future__ import annotations

from chemsmart.agent.rules import POLICY_RULES, render_rules


def _ids():
    return {rule.rule_id for rule in POLICY_RULES}


def test_the_wave_contract_is_registered():
    for rule_id in (
        "stem.wave_execution",
        "stem.width_is_yours_concurrency_is_the_hosts",
        "stem.hardware_is_the_hosts",
        "wake.cohort_evidence",
    ):
        assert rule_id in _ids(), f"{rule_id} is not registered"


def test_the_rendered_stem_actually_carries_it():
    stem = render_rules("stem").lower()
    assert "wave" in stem
    # Woken once, and on terminality rather than success.
    assert "once" in stem
    assert "terminal" in stem
    # The two quantities that must not be collapsed.
    assert "at the same time" in stem or "concurrency" in stem


def test_the_wake_says_a_whole_wave_is_being_read():
    wake = render_rules("wake").lower()
    assert "wave" in wake
    assert "fail" in wake or "cancel" in wake, (
        "the wake does not tell the model that a failed member is part of "
        "the evidence it asked for, so it may read one as an error to "
        "route around"
    )


def test_the_model_is_told_it_does_not_size_the_machine():
    stem = render_rules("stem").lower()
    for forbidden in ("core", "memory"):
        assert forbidden in stem, (
            "the model is not told which knobs are the host's, so it may "
            f"try to choose {forbidden}"
        )


def test_no_tool_lets_the_model_write_a_hardware_setting():
    """Inspect the rooms; never rewrite them.

    Asserted over the whole tool surface rather than one tool, because a
    single writable path anywhere makes the sentence false.
    """

    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    specs = build_command_compiled_tool_surface().tool_definitions
    names = {str(spec.get("function", {}).get("name") or "") for spec in specs}
    for writer in (
        "set_execution_resources",
        "configure_server",
        "set_num_cores",
        "set_max_concurrent_tasks",
    ):
        assert (
            writer not in names
        ), f"{writer} lets the model rewrite host hardware policy"
    blob = repr(specs).lower()
    for key in ("num_cores", "num_threads", "mem_gb", "max_concurrent_tasks"):
        assert f'"{key}"' not in blob or "inspect" in blob, (
            f"{key} appears in the tool surface as something the model "
            "supplies rather than reads"
        )
