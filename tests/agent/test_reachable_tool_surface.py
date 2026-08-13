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
_REGISTRY_ATTRIBUTES = {
    "counterexample": "counterexamples",
    "canonical invocation": "invocations",
    "command inspection receipt": "command_inspections",
}


def _has_a_producer(attribute: str) -> bool:
    source = inspect.getsource(tool_runtime)
    return f"self.{attribute}[" in source


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


def test_approved_execution_does_not_offer_the_legacy_seeded_result_verifier():
    """Runtime V2 execution has no legacy settings/run-ID producers.

    Newly executed results are validated by the execution receipt and then
    exposed through typed quantity extraction.  Offering the older verifier
    here can only induce guesses at IDs that the live host never binds.
    """

    preview_tools = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface().tool_definitions
    }
    execution_tools = {
        item["function"]["name"]
        for item in build_approved_execution_tool_surface().tool_definitions
    }
    assert "inspect_calculation_artifact" in preview_tools
    assert "inspect_calculation_artifact" not in execution_tools
    assert "execute_approved_program_node" in execution_tools


def test_the_handler_and_contract_are_kept_so_re_enabling_is_one_change():
    """Withholding the tool must not mean deleting the feature."""

    from chemsmart.agent.commands import CommandCounterexampleV1

    assert CommandCounterexampleV1 is not None
    assert "repair_command" in inspect.getsource(tool_runtime), (
        "the runtime handler must survive so re-exposing is a surface change"
    )
