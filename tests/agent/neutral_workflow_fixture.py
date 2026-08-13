"""Small provider-neutral scientific workflow fixture for Runtime tests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from chemsmart.agent._contracts import (
    TrustedArtifactRefV1,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.commands import build_scientific_identity_binding
from chemsmart.agent.capabilities import (
    CapabilityQueryV1,
    build_command_compiled_preview_overlay,
    load_program_capabilities,
    query_capability,
)
from chemsmart.agent.cli_schema import build_live_click_schema


@dataclass(frozen=True)
class NeutralWorkflowFixtureV1:
    public_context: Any
    host_inputs: dict[str, Any]


def build_neutral_workflow_fixture(root: Path) -> NeutralWorkflowFixtureV1:
    """Bind one exact geometry and an OPT to HESS data dependency."""

    root.mkdir(parents=True, exist_ok=True)
    geometry_path = root / "water.xyz"
    geometry_path.write_text(
        "3\nneutral workflow fixture\n"
        "O 0.000000 0.000000 0.000000\n"
        "H 0.000000 0.757000 0.586000\n"
        "H 0.000000 -0.757000 0.586000\n",
        encoding="utf-8",
    )
    geometry = TrustedArtifactRefV1(
        artifact_id="geometry.initial",
        kind="xyz",
        sha256=file_sha256(geometry_path),
        size_bytes=geometry_path.stat().st_size,
        path=str(geometry_path.resolve()),
        cli_value=str(geometry_path.resolve()),
    )
    task_spec_sha256 = canonical_sha256(
        {
            "geometry_sha256": geometry.sha256,
            "program": "pyscf",
            "stages": ("opt", "hess"),
        }
    )
    identity = build_scientific_identity_binding(
        task_spec_sha256=task_spec_sha256,
        geometry_artifact=geometry,
        charge=0,
        multiplicity=1,
    )
    project_path = root / "pyscf.yaml"
    project_path.write_text(
        "opt:\n  ab_initio: hf\n  basis: def2-svp\n"
        "hess:\n  ab_initio: hf\n  basis: def2-svp\n",
        encoding="utf-8",
    )
    project = TrustedArtifactRefV1(
        artifact_id="project.pyscf.cpu-hf",
        kind="project_yaml",
        sha256=file_sha256(project_path),
        size_bytes=project_path.stat().st_size,
        path=str(project_path.resolve()),
        cli_value=str(project_path.resolve()),
    )
    registry = load_program_capabilities()
    live_schema = build_live_click_schema()
    overlay = build_command_compiled_preview_overlay(
        registry, live_schema=live_schema
    )
    capabilities = tuple(
        query_capability(
            CapabilityQueryV1(
                program="pyscf", jobtype=jobtype, engine="cpu"
            ),
            registry=registry,
            overlay=overlay,
            live_schema=live_schema,
        )
        for jobtype in ("opt", "hess")
    )
    nodes = [
        {
            "node_id": "node.opt",
            "program": "pyscf",
            "jobtype": "opt",
            "project_role": "project.primary",
            "dependencies": [],
            "inputs": [
                {
                    "binding_id": "geometry.initial",
                    "artifact_id": geometry.artifact_id,
                    "artifact_class": "xyz",
                    "producer_node_id": "",
                    "producer_output_id": "",
                }
            ],
            "expected_outputs": [
                {
                    "output_id": "optimized_geometry",
                    "artifact_class": "xyz",
                },
                {
                    "output_id": "structured_result",
                    "artifact_class": "pyscf_hdf5",
                },
            ],
            "unresolved_fields": [],
        },
        {
            "node_id": "node.hess",
            "program": "pyscf",
            "jobtype": "hess",
            "project_role": "project.primary",
            "dependencies": ["node.opt"],
            "inputs": [
                {
                    "binding_id": "geometry.optimized",
                    "artifact_class": "xyz",
                    "producer_node_id": "node.opt",
                    "producer_output_id": "optimized_geometry",
                }
            ],
            "expected_outputs": [
                {
                    "output_id": "structured_result",
                    "artifact_class": "pyscf_hdf5",
                }
            ],
            "unresolved_fields": ["input.optimized_geometry_hash"],
        },
    ]
    action = SimpleNamespace(
        tool_name="plan_command_workflow",
        fields=tuple(
            sorted(
                {
                    "workflow_id": "workflow-neutral-opt-hess",
                    "task_spec_id": task_spec_sha256,
                    "nodes": nodes,
                }.items()
            )
        ),
    )
    return NeutralWorkflowFixtureV1(
        public_context=SimpleNamespace(
            task_spec_sha256=task_spec_sha256,
            next_actions=(action,),
        ),
        host_inputs={
            "artifacts": {
                geometry.artifact_id: geometry,
                project.artifact_id: project,
            },
            "scientific_identities": {identity.binding_sha256: identity},
            "capability_receipts": {
                item.receipt_sha256: item for item in capabilities
            },
            "registry": registry,
            "live_schema": live_schema,
        },
    )


__all__ = ["NeutralWorkflowFixtureV1", "build_neutral_workflow_fixture"]
