"""Changing which project a node uses must not change what it calculates.

`amend_scientific_workflow` exists so a session can repair one node instead of
re-authoring the graph. When the repair is a project replacement, the
calculation side has to be recompiled, and the node is handed to
`_plan_command_workflow` as a dict whose fields are listed by hand. Anything
missing from that list is dropped in silence.

`internal_coordinates` was missing. A scan whose project role was amended came
back out of the rebuild still typed `scan`, still displayed as a scan, and
compiled with no `--coordinates`, `--dist-start`, `--dist-end` or `--num-steps`
at all -- a plain optimisation wearing a scan's name. That was visible in the
qualification session's second attempt, where one compiled scan invocation
carried only server, project, filename, charge and multiplicity.

The existing amend checks could not see it: their draft holds no calculation
nodes, so the reconstruction they cover is the empty one.

This drives the real `_rebind_scientific_workflow_projects` and asserts on what
it hands downstream, because the handoff is exactly where the loss happened.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)

_COORDINATES = {
    "scan": {
        "kind": "dihedral",
        "atoms": (3, 1, 2, 4),
        "start": 0.0,
        "stop": 360.0,
        "points": 25,
    }
}


def _scan_node(project_role="scan-project-v1"):
    return CommandNodeIntentV1(
        node_id="orca-scan-hooh",
        program="orca",
        jobtype="scan",
        project_role=project_role,
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="filename",
                artifact_id="geometry-hooh",
                artifact_class="geometry_xyz",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="scan-reader", artifact_class="orca_output"
            ),
        ),
        unresolved_fields=(),
        internal_coordinates=_COORDINATES,
    )


class _Draft:
    workflow_id = "hooh-torsion"
    task_spec_id = "a" * 64
    draft_sha256 = "d" * 64
    nodes = (_scan_node(),)


class _Validation:
    project_artifact_id = "scan-project-v2"
    project_sha256 = "b" * 64
    program = "orca"
    jobtype = "scan"
    status = "valid"


@pytest.fixture()
def host_and_capture():
    """A host holding one planned scan, plus the downstream call captured."""

    host = CommandCompiledToolHostV1.__new__(CommandCompiledToolHostV1)
    replacement = TrustedArtifactRefV1(
        artifact_id="scan-project-v2",
        kind="project_yaml",
        sha256="b" * 64,
        size_bytes=44,
        path="/tmp/scan-project-v2.yaml",
        cli_value="/tmp/scan-project-v2.yaml",
    )
    host.artifacts = {replacement.artifact_id: replacement}
    host.project_validations = {"v": _Validation()}
    host.frozen_workflow_approval = None
    host.scientific_plans = {}
    host._bind_program_toolchain = lambda *args, **kwargs: None

    captured = {}

    def _capture(turn_id, values, node_annotations=None):
        captured["nodes"] = values["nodes"]
        return {
            "workflow_draft": _Draft(),
            "actionable_node_ids": (),
            "unresolved_node_ids": (),
        }

    host._plan_command_workflow = _capture
    return host, captured


def _rebind(host):
    return host._rebind_scientific_workflow_projects(
        "turn-1",
        {
            "workflow_id": "hooh-torsion",
            "replacements": (
                {
                    "node_id": "orca-scan-hooh",
                    "project_role": "scan-project-v2",
                },
            ),
            "analysis_repairs": (),
        },
    )


def _prepare(host, captured, plan_state):
    host.scientific_toolchain_plans = {plan_state.plan_sha256: plan_state}
    host._scientific_toolchain_command_results = {
        plan_state.plan_sha256: {
            "workflow_draft": _Draft(),
            "actionable_node_ids": (),
            "unresolved_node_ids": (),
        }
    }


def test_the_rebuilt_node_still_drives_its_coordinate(host_and_capture):
    host, captured = host_and_capture
    from chemsmart.agent.scientific_toolchain import (
        build_scientific_toolchain_plan,
    )

    plan = build_scientific_toolchain_plan(
        plan_id="hooh-torsion",
        workflow_id="hooh-torsion",
        command_workflow_draft_sha256="d" * 64,
        calculation_nodes=_Draft.nodes,
        calculation_observables={"orca-scan-hooh": ("scan_energies",)},
        analysis_nodes=(),
        required_output_ids=(),
    )
    _prepare(host, captured, plan)

    _rebind(host)

    rebuilt = {node["node_id"]: node for node in captured["nodes"]}
    node = rebuilt["orca-scan-hooh"]
    assert node["project_role"] == "scan-project-v2"
    assert node["internal_coordinates"] == _COORDINATES


def test_a_node_that_drives_nothing_stays_absent(host_and_capture):
    """The field must not appear as None on an ordinary optimisation."""

    host, captured = host_and_capture
    from dataclasses import replace as _replace

    from chemsmart.agent.scientific_toolchain import (
        build_scientific_toolchain_plan,
    )

    plain = _replace(_scan_node(), jobtype="opt", internal_coordinates=None)

    class _PlainDraft(_Draft):
        nodes = (plain,)

    plan = build_scientific_toolchain_plan(
        plan_id="hooh-torsion",
        workflow_id="hooh-torsion",
        command_workflow_draft_sha256="d" * 64,
        calculation_nodes=_PlainDraft.nodes,
        calculation_observables={"orca-scan-hooh": ("energy",)},
        analysis_nodes=(),
        required_output_ids=(),
    )
    host.scientific_toolchain_plans = {plan.plan_sha256: plan}
    host._scientific_toolchain_command_results = {
        plan.plan_sha256: {
            "workflow_draft": _PlainDraft(),
            "actionable_node_ids": (),
            "unresolved_node_ids": (),
        }
    }
    _Validation.jobtype = "opt"
    try:
        _rebind(host)
    finally:
        _Validation.jobtype = "scan"

    node = {item["node_id"]: item for item in captured["nodes"]}[
        "orca-scan-hooh"
    ]
    assert "internal_coordinates" not in node
