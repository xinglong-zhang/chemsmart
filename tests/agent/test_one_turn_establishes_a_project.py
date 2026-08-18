"""Three calls that are never used apart cost three turns per node.

`render_project_yaml`, `promote_project_yaml` and `validate_project_yaml` are
always used together, always in that order, and every calculation node needs its
own. The ceremony therefore scales with the size of the graph, and it is not a
small tax: measured across three live sessions it took 34%, 56% and 61% of the
entire tool budget. The paper-reproduction task ran out of turns inside it,
having reasoned about the chemistry correctly and never reached a reviewable
workflow.

`establish_project` runs the same three handlers, in the same order, and returns
all three receipts. The point of these checks is that it is *equivalent* -- if
it ever became a shortcut that skipped a step, it would be worth less than the
turns it saves.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_SECTIONS = {"gas": {"functional": "b3lyp", "basis": "def2-tzvp"}}


def _field(item, name):
    """Dispatch serialises receipts; the handler returns them as objects."""

    if isinstance(item, dict):
        return item[name]
    return getattr(item, name)


def _digest(item):
    return _field(item, "receipt_sha256")


def _host(tmp_path, name):
    workspace = tmp_path / name
    workspace.mkdir()
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / f"{name}.jsonl", session_id=name
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=workspace,
    )


def _call(host, tool, **arguments):
    return host.dispatch(turn_id="t1", tool_name=tool, arguments=arguments)[
        "result"
    ]


def _capability(host):
    return _call(
        host,
        "inspect_program_capability",
        program="orca",
        jobtype="opt",
        engine="cpu",
    )["receipt_sha256"]


def _separately(host, artifact_id):
    capability = _capability(host)
    rendered = _call(
        host, "render_project_yaml", program="orca", sections=_SECTIONS
    )
    promoted = _call(
        host,
        "promote_project_yaml",
        render_receipt_sha256=_digest(rendered),
        artifact_id=artifact_id,
    )
    validated = _call(
        host,
        "validate_project_yaml",
        project_artifact_id=artifact_id,
        capability_receipt_sha256=capability,
    )
    return rendered, promoted, validated


def _together(host, artifact_id):
    capability = _capability(host)
    result = _call(
        host,
        "establish_project",
        program="orca",
        sections=_SECTIONS,
        artifact_id=artifact_id,
        capability_receipt_sha256=capability,
    )
    return result["rendered"], result["promoted"], result["validated"]


def test_the_same_receipts_come_back_either_way(tmp_path):
    """Equivalence, digest for digest."""

    apart = _separately(_host(tmp_path, "apart"), "project-a")
    together = _together(_host(tmp_path, "together"), "project-a")

    assert _digest(apart[0]) == _digest(together[0])
    assert _field(apart[1]["artifact"], "sha256") == _field(
        together[1]["artifact"], "sha256"
    )
    assert _field(apart[2], "receipt_sha256") == _field(
        together[2], "receipt_sha256"
    )


def test_the_project_is_bound_and_usable_afterwards(tmp_path):
    """A node must be able to name it, which is the whole purpose."""

    host = _host(tmp_path, "bound")
    _together(host, "project-b")

    assert "project-b" in host.artifacts
    assert "project-b" in host.project_promotions


def test_the_validation_result_is_reported_not_swallowed(tmp_path):
    """A combined call must not hide a project that failed to validate."""

    host = _host(tmp_path, "reported")
    _, _, validated = _together(host, "project-c")

    assert _field(validated, "status")


def test_the_separate_tools_remain_available(tmp_path):
    """Combining is an option, not a replacement; replay needs the old path."""

    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    names = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface().tool_definitions
    }

    assert {
        "render_project_yaml",
        "promote_project_yaml",
        "validate_project_yaml",
        "establish_project",
    } <= names


def test_a_taken_artifact_id_is_still_refused(tmp_path):
    """The promotion guard must survive being called from inside."""

    from chemsmart.agent._contracts import ContractError

    host = _host(tmp_path, "taken")
    _together(host, "project-d")

    with pytest.raises(ContractError):
        _together(host, "project-d")
