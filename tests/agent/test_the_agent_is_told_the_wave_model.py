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

    # Walk the declared input properties of every tool rather than
    # matching quoted substrings in a repr. The first version searched
    # `repr(specs)` for double-quoted keys while a repr generally yields
    # single quotes, and it passed whenever the word "inspect" appeared
    # anywhere on the surface -- so a writable hardware parameter under
    # any new name would have kept it green.
    # A denylist of exact names is a list of the mistakes already made:
    # a writable `nprocs`, `threads`, `partition`, `gpus` or `wall_hours`
    # walks past it, which is how this test's own predecessor failed. The
    # rule is about *what the field is*, so it is matched by the parts a
    # hardware field is spelled from, at every depth -- nested objects
    # and array items included, because a wave's members arrive as an
    # array and a resource block would arrive as an object.
    hardware_words = {
        "core",
        "cores",
        "cpu",
        "cpus",
        "nproc",
        "nprocs",
        "thread",
        "threads",
        "omp",
        "mem",
        "memory",
        "ram",
        "gpu",
        "gpus",
        "queue",
        "partition",
        "qos",
        "account",
        "ntasks",
        "walltime",
        "wall",
        "concurrency",
        "concurrent",
        "parallel",
        "slurm",
        "sbatch",
        "scheduler",
    }
    # Fields that name hardware in order to *read* or *record* it. Each
    # is listed because the owner's ruling is "inspect, never set", and
    # each is a value the host wrote that the model quotes back.
    readers = {
        "inspect_program_environment.engine",
        "inspect_program.engine",
        "synthesize_command.engine",
        "prepare_program_node.engine",
    }
    # "node" is this codebase's word for a step of a scientific DAG, so
    # it cannot be matched as a word. A *compute* node count has its own
    # spellings, and those are matched exactly.
    node_counts = {"num_nodes", "nnodes", "node_count", "nodes_per_job"}

    def _words(key: str) -> set:
        return {
            part
            for part in str(key).replace("-", "_").lower().split("_")
            if part
        }

    def _walk(schema, path, offenders):
        if not isinstance(schema, dict):
            return
        for key, value in (schema.get("properties") or {}).items():
            here = f"{path}.{key}"
            spelled = str(key).replace("-", "_").lower()
            if (
                _words(key) & hardware_words or spelled in node_counts
            ) and here not in readers:
                offenders.append(here)
            _walk(value, here, offenders)
        items = schema.get("items")
        if isinstance(items, dict):
            _walk(items, f"{path}[]", offenders)

    offenders: list[str] = []
    for spec in specs:
        function = spec.get("function") or {}
        _walk(
            function.get("parameters") or {},
            str(function.get("name") or ""),
            offenders,
        )
    assert not offenders, (
        "these tool inputs let the model supply host hardware policy "
        f"rather than read it: {sorted(offenders)}"
    )


def test_the_hardware_guard_can_fail():
    """A guard that cannot go red guards nothing.

    Its predecessor searched a repr for double-quoted keys while a repr
    yields single ones, and passed whenever the word "inspect" appeared
    anywhere on the surface.
    """

    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    specs = list(build_command_compiled_tool_surface().tool_definitions)
    specs.append(
        {
            "function": {
                "name": "invented",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "members": {
                            "type": "array",
                            "items": {
                                "type": "object",
                                "properties": {"nprocs": {"type": "integer"}},
                            },
                        }
                    },
                },
            }
        }
    )

    hardware_words = {"nproc", "nprocs"}

    def _words(key):
        return {p for p in str(key).lower().split("_") if p}

    found = []

    def _walk(schema, path):
        if not isinstance(schema, dict):
            return
        for key, value in (schema.get("properties") or {}).items():
            if _words(key) & hardware_words:
                found.append(f"{path}.{key}")
            _walk(value, f"{path}.{key}")
        if isinstance(schema.get("items"), dict):
            _walk(schema["items"], f"{path}[]")

    for spec in specs:
        function = spec.get("function") or {}
        _walk(function.get("parameters") or {}, function.get("name") or "")
    assert found == ["invented.members[].nprocs"], found
