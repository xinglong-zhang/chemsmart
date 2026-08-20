"""A tool the harness cannot let succeed must not be advertised.

Four live sessions across three cases called ``repair_command`` with no
counterexample bound.  Three separate attempts to discourage it by wording --
describing the argument, stating the precondition on that tool, then deriving
the precondition for every tool in its class -- all failed.

The reason is that ``CommandCounterexampleV1`` is never constructed anywhere in
the package and ``register_counterexample`` has no callers, so the registry is
empty for the whole lifetime of every session.  The tool cannot succeed.  It
was not a wording problem, and no wording could have fixed it.

The general rule these tests pin: the surface must advertise only what the
harness can actually deliver.  That is the same principle the capability
registry already applies when it declares PySCF ``td`` as
``execution_supported: false`` instead of offering it as runnable.
"""

import inspect
from pathlib import Path

import pytest

from chemsmart.agent import tool_runtime
from chemsmart.agent.tool_specs import (
    LATE_BOUND_ARGUMENTS,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)

_SURFACES = (
    build_command_compiled_tool_surface,
    build_approved_execution_tool_surface,
)

#: Registry attribute on the host for each late-bound argument, so "can this
#: ever be filled?" is answered from code rather than from belief.
#:
#: The last three were missing, and their absence is why an unfillable tool
#: stayed advertised in the preview surface for so long: the rule below was
#: stated generally but only ever checked against the registries someone had
#: remembered to list here.  Every label in ``LATE_BOUND_ARGUMENTS`` that names
#: a host registry now appears, and the test below asserts that too, so the
#: next omission fails rather than passing silently.
_REGISTRY_ATTRIBUTES = {
    "counterexample": "counterexamples",
    "canonical invocation": "invocations",
    "command inspection receipt": "command_inspections",
    "run receipt": "run_receipts",
    "settings object": "settings_objects",
    "scientific claim evidence": "scientific_claim_evidence",
}

#: Entry points that construct the live host.  A registry with no in-session
#: writer could still be filled at construction, so both are checked before a
#: tool is called unreachable.
_LIVE_ENTRY_POINTS = (
    "chemsmart/agent/live_session.py",
    "chemsmart/agent/bootstrap.py",
    "chemsmart/cli/agent.py",
)


def _has_a_producer(attribute: str) -> bool:
    """Whether anything can put an entry in this registry during a session."""

    source = inspect.getsource(tool_runtime)
    if f"self.{attribute}[" in source:
        return True
    # Constructor seeding counts too: the V1 surface was filled that way.
    root = Path(inspect.getsourcefile(tool_runtime)).parents[2]
    return any(
        attribute in (root / entry).read_text(encoding="utf-8")
        for entry in _LIVE_ENTRY_POINTS
        if (root / entry).exists()
    )


@pytest.mark.parametrize("build", _SURFACES)
def test_no_exposed_tool_depends_on_a_registry_nothing_fills(build):
    unreachable = []
    for item in build().tool_definitions:
        function = item["function"]
        properties = (function.get("parameters") or {}).get("properties") or {}
        for name in properties:
            label = LATE_BOUND_ARGUMENTS.get(name)
            attribute = _REGISTRY_ATTRIBUTES.get(label or "")
            if attribute and not _has_a_producer(attribute):
                unreachable.append(f"{function['name']} needs {label}")
    assert not unreachable, (
        "these tools are advertised but cannot succeed, because nothing in "
        f"the runtime fills the registry they read: {sorted(set(unreachable))}"
    )


def test_the_counterexample_registry_still_has_no_producer():
    """Pins the reason repair_command is withheld.

    When a producer is wired this test fails, which is the signal to re-expose
    the tool in the same change.
    """

    source = inspect.getsource(tool_runtime)
    assert "self.counterexamples[" in source, "the registry write exists"
    assert source.count("register_counterexample") == 1, (
        "register_counterexample is defined once and called nowhere; if that "
        "changed, re-expose repair_command"
    )


@pytest.mark.parametrize("build", _SURFACES)
def test_repair_command_is_not_offered_while_it_cannot_work(build):
    exposed = {item["function"]["name"] for item in build().tool_definitions}
    assert "repair_command" not in exposed


def test_the_legacy_seeded_verifier_is_offered_by_neither_surface():
    """Runtime V2 has no legacy settings/run-ID producers -- in either surface.

    This test previously asserted the verifier belonged in the preview surface
    and only had to go from the execution one.  That was belief, not evidence:
    ``settings_objects`` and ``run_receipts`` are constructor-injected only, no
    live entry point seeds them, and nothing writes them during a session, so
    the tool could never succeed in *either* surface.  Advertising it could
    only induce guesses at identifiers the live host never binds -- the same
    defect as ``repair_command``, missed because the general rule above was
    checked against an incomplete registry list.

    ``assess_program_candidate`` is here for the same reason: it reads claim
    evidence, which is equally unbound.
    """

    preview_tools = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface().tool_definitions
    }
    execution_tools = {
        item["function"]["name"]
        for item in build_approved_execution_tool_surface().tool_definitions
    }
    for withheld in (
        "inspect_calculation_artifact",
        "assess_program_candidate",
    ):
        assert withheld not in preview_tools
        assert withheld not in execution_tools
    assert "execute_approved_program_node" in execution_tools


def test_every_late_bound_registry_is_actually_checked():
    """The omission that let an unfillable tool stay advertised.

    ``LATE_BOUND_ARGUMENTS`` names the registry behind each deferred argument.
    Any of those that is a host registry must appear in the map above, or the
    reachability rule silently stops applying to it.
    """

    host_registries = {
        label
        for label in LATE_BOUND_ARGUMENTS.values()
        if hasattr(tool_runtime.CommandCompiledToolHostV1, "__init__")
        and label not in {"project render receipt"}
    }
    missing = sorted(host_registries.difference(_REGISTRY_ATTRIBUTES))
    assert not missing, (
        "these late-bound registries are never checked for a producer, so a "
        f"tool needing one could be advertised unreachably: {missing}"
    )


def test_tools_needing_an_unbound_registry_keep_their_handlers():
    """Withholding a tool must not mean deleting the feature."""

    source = inspect.getsource(tool_runtime)
    for handler in (
        "_inspect_calculation_artifact",
        "_assess_program_candidate",
    ):
        assert (
            handler in source
        ), f"{handler} must survive so re-exposing it is a surface change"


def test_the_handler_and_contract_are_kept_so_re_enabling_is_one_change():
    """Withholding the tool must not mean deleting the feature."""

    from chemsmart.agent.commands import CommandCounterexampleV1

    assert CommandCounterexampleV1 is not None
    assert "repair_command" in inspect.getsource(
        tool_runtime
    ), "the runtime handler must survive so re-exposing is a surface change"
