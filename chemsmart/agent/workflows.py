"""Immutable multi-command workflow contracts; planning and preview only."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
    require_sha256,
)


@dataclass(frozen=True)
class ArtifactInputIntentV1:
    """A semantic input slot; it may remain unresolved during planning."""

    binding_id: str
    artifact_class: str
    artifact_id: str
    producer_node_id: str
    producer_output_id: str

    def __post_init__(self) -> None:
        require_identifier(self.binding_id, "binding_id")
        require_identifier(self.artifact_class, "artifact_class")
        produced = bool(self.producer_node_id or self.producer_output_id)
        if produced and not (self.producer_node_id and self.producer_output_id):
            raise ContractError("producer input needs node and output IDs")
        if produced:
            require_identifier(self.producer_node_id, "producer_node_id")
            require_identifier(self.producer_output_id, "producer_output_id")
            if self.artifact_id:
                raise ContractError(
                    "future producer input cannot claim a materialized artifact"
                )


@dataclass(frozen=True)
class ArtifactOutputIntentV1:
    output_id: str
    artifact_class: str

    def __post_init__(self) -> None:
        require_identifier(self.output_id, "output_id")
        require_identifier(self.artifact_class, "artifact_class")


@dataclass(frozen=True)
class CommandNodeIntentV1:
    """Broad scientific command intent before execution-grade grounding."""

    node_id: str
    program: str
    jobtype: str
    project_role: str
    dependencies: tuple[str, ...]
    inputs: tuple[ArtifactInputIntentV1, ...]
    expected_outputs: tuple[ArtifactOutputIntentV1, ...]
    unresolved_fields: tuple[str, ...]

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        require_identifier(self.program, "program")
        require_identifier(self.jobtype, "jobtype")
        if self.dependencies != tuple(sorted(set(self.dependencies))):
            raise ContractError("draft dependencies must be sorted and unique")
        binding_ids = tuple(item.binding_id for item in self.inputs)
        if binding_ids != tuple(sorted(set(binding_ids))):
            raise ContractError("draft input bindings must be sorted and unique")
        output_ids = tuple(item.output_id for item in self.expected_outputs)
        if output_ids != tuple(sorted(set(output_ids))):
            raise ContractError("draft outputs must be sorted and unique")
        if not self.expected_outputs:
            raise ContractError("draft command node requires an expected output")
        if self.unresolved_fields != tuple(sorted(set(self.unresolved_fields))):
            raise ContractError("unresolved fields must be sorted and unique")


@dataclass(frozen=True)
class CommandWorkflowDraftV1:
    """Planning-level DAG that permits honest missing evidence."""

    schema_version: str
    workflow_id: str
    task_spec_id: str
    nodes: tuple[CommandNodeIntentV1, ...]
    status: str
    draft_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.command-workflow-draft.v1":
            raise ContractError("unsupported command workflow draft schema")
        require_identifier(self.workflow_id, "workflow_id")
        if not str(self.task_spec_id).strip():
            raise ContractError("task_spec_id must not be empty")
        if self.status != "planned":
            raise ContractError("workflow draft can only be planned")
        _validate_draft_dag(self.nodes)
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "task_spec_id": self.task_spec_id,
            "nodes": self.nodes,
            "status": self.status,
        }
        if self.draft_sha256 != canonical_sha256(body):
            raise ContractError("command workflow draft digest mismatch")


def build_command_workflow_draft(
    *,
    workflow_id: str,
    task_spec_id: str,
    nodes: tuple[CommandNodeIntentV1, ...],
) -> CommandWorkflowDraftV1:
    body = {
        "schema_version": "chemsmart.command-workflow-draft.v1",
        "workflow_id": workflow_id,
        "task_spec_id": task_spec_id,
        "nodes": nodes,
        "status": "planned",
    }
    return CommandWorkflowDraftV1(
        **body, draft_sha256=canonical_sha256(body)
    )


def _validate_draft_dag(nodes: tuple[CommandNodeIntentV1, ...]) -> None:
    if not nodes:
        raise ContractError("command workflow draft must contain a node")
    seen: dict[str, set[str]] = {}
    for node in nodes:
        if node.node_id in seen:
            raise ContractError("draft workflow node IDs must be unique")
        if set(node.dependencies).difference(seen):
            raise ContractError("draft workflow must be topologically ordered")
        for item in node.inputs:
            if not item.producer_node_id:
                continue
            if item.producer_node_id not in node.dependencies:
                raise ContractError("producer must be a direct dependency")
            if item.producer_output_id not in seen[item.producer_node_id]:
                raise ContractError("draft input references an unknown output")
        seen[node.node_id] = {item.output_id for item in node.expected_outputs}


@dataclass(frozen=True)
class ArtifactOutputV1:
    output_id: str
    artifact_id: str
    artifact_class: str

    def __post_init__(self) -> None:
        require_identifier(self.output_id, "output_id")
        require_identifier(self.artifact_id, "artifact_id")
        require_identifier(self.artifact_class, "artifact_class")


@dataclass(frozen=True)
class ArtifactBindingV1:
    """Exact external artifact or declared producer-output edge."""

    binding_id: str
    artifact_id: str
    artifact_class: str
    artifact_sha256: str
    producer_node_id: str
    producer_output_id: str
    producer_receipt_sha256: str

    def __post_init__(self) -> None:
        require_identifier(self.binding_id, "binding_id")
        require_identifier(self.artifact_id, "artifact_id")
        require_identifier(self.artifact_class, "artifact_class")
        produced = bool(self.producer_node_id or self.producer_output_id)
        if produced and not (self.producer_node_id and self.producer_output_id):
            raise ContractError("produced artifact needs node and output IDs")
        if produced:
            require_identifier(self.producer_node_id, "producer_node_id")
            require_identifier(self.producer_output_id, "producer_output_id")
            if bool(self.artifact_sha256) != bool(self.producer_receipt_sha256):
                raise ContractError(
                    "resolved producer artifact needs hash and producer receipt"
                )
        elif not self.artifact_sha256:
            raise ContractError("external artifact binding requires an exact hash")
        if self.artifact_sha256:
            require_sha256(self.artifact_sha256, "artifact_sha256")
        if self.producer_receipt_sha256:
            require_sha256(
                self.producer_receipt_sha256, "producer_receipt_sha256"
            )

    @property
    def resolved(self) -> bool:
        return bool(self.artifact_sha256)


@dataclass(frozen=True)
class CommandNodeV1:
    node_id: str
    program: str
    jobtype: str
    capability_receipt_sha256: str
    program_binding_sha256: str
    engine_binding_sha256: str
    scientific_identity_sha256: str
    project_artifact_id: str
    project_sha256: str
    dependencies: tuple[str, ...]
    inputs: tuple[ArtifactBindingV1, ...]
    expected_outputs: tuple[ArtifactOutputV1, ...]
    invocation_sha256: str
    preflight_receipt_sha256: str
    state: str

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        require_identifier(self.program, "program")
        require_identifier(self.jobtype, "jobtype")
        for field_name, digest in (
            ("capability_receipt_sha256", self.capability_receipt_sha256),
            ("program_binding_sha256", self.program_binding_sha256),
            ("engine_binding_sha256", self.engine_binding_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
        ):
            require_sha256(digest, field_name)
        if self.project_sha256:
            require_sha256(self.project_sha256, "project_sha256")
        if bool(self.project_artifact_id) != bool(self.project_sha256):
            raise ContractError("project artifact ID and hash must be paired")
        if self.dependencies != tuple(sorted(set(self.dependencies))):
            raise ContractError("node dependencies must be sorted and unique")
        binding_ids = tuple(item.binding_id for item in self.inputs)
        if binding_ids != tuple(sorted(set(binding_ids))):
            raise ContractError("input bindings must be sorted and unique")
        output_ids = tuple(item.output_id for item in self.expected_outputs)
        if output_ids != tuple(sorted(set(output_ids))):
            raise ContractError("expected outputs must be sorted and unique")
        if not self.inputs or not self.expected_outputs:
            raise ContractError("command node requires inputs and outputs")
        if self.state not in {"planned", "compiled", "previewed"}:
            raise ContractError("command workflow cannot claim execution")
        if self.state in {"compiled", "previewed"}:
            require_sha256(self.invocation_sha256, "invocation_sha256")
            if not all(item.resolved for item in self.inputs):
                raise ContractError("unresolved inputs cannot be compiled")
        elif self.invocation_sha256:
            raise ContractError("planned node cannot carry an invocation")
        if self.state == "previewed":
            require_sha256(
                self.preflight_receipt_sha256, "preflight_receipt_sha256"
            )
        elif self.preflight_receipt_sha256:
            raise ContractError("only previewed node carries preflight evidence")


@dataclass(frozen=True)
class CommandWorkflowSpecV1:
    schema_version: str
    workflow_id: str
    task_spec_sha256: str
    live_cli_schema_sha256: str
    nodes: tuple[CommandNodeV1, ...]
    workflow_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.command-workflow-spec.v1":
            raise ContractError("unsupported command workflow schema")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        require_sha256(
            self.live_cli_schema_sha256, "live_cli_schema_sha256"
        )
        _validate_dag(self.nodes)
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "task_spec_sha256": self.task_spec_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "nodes": self.nodes,
        }
        if self.workflow_sha256 != canonical_sha256(body):
            raise ContractError("command workflow digest mismatch")


def build_command_workflow_spec(
    *,
    workflow_id: str,
    task_spec_sha256: str,
    live_cli_schema_sha256: str,
    nodes: tuple[CommandNodeV1, ...],
) -> CommandWorkflowSpecV1:
    body = {
        "schema_version": "chemsmart.command-workflow-spec.v1",
        "workflow_id": require_identifier(workflow_id, "workflow_id"),
        "task_spec_sha256": task_spec_sha256,
        "live_cli_schema_sha256": live_cli_schema_sha256,
        "nodes": nodes,
    }
    return CommandWorkflowSpecV1(
        **body, workflow_sha256=canonical_sha256(body)
    )


def _validate_dag(nodes: tuple[CommandNodeV1, ...]) -> None:
    if not nodes:
        raise ContractError("command workflow must contain a node")
    node_ids = tuple(item.node_id for item in nodes)
    if len(node_ids) != len(set(node_ids)):
        raise ContractError("command workflow node IDs must be unique")
    seen: dict[str, CommandNodeV1] = {}
    global_artifact_ids: dict[str, tuple[str, str]] = {}
    for node in nodes:
        unknown_or_forward = set(node.dependencies).difference(seen)
        if unknown_or_forward:
            raise ContractError(
                "workflow order must be topological and dependencies known"
            )
        for output in node.expected_outputs:
            if output.artifact_id in global_artifact_ids:
                raise ContractError("two nodes declare the same output artifact")
            global_artifact_ids[output.artifact_id] = (
                node.node_id,
                output.output_id,
            )
        for binding in node.inputs:
            if not binding.producer_node_id:
                continue
            if binding.producer_node_id not in node.dependencies:
                raise ContractError(
                    "producer node must be an explicit direct dependency"
                )
            producer = seen[binding.producer_node_id]
            output = next(
                (
                    item
                    for item in producer.expected_outputs
                    if item.output_id == binding.producer_output_id
                ),
                None,
            )
            if output is None:
                raise ContractError("artifact binds an undeclared producer output")
            if (
                binding.artifact_id != output.artifact_id
                or binding.artifact_class != output.artifact_class
            ):
                raise ContractError(
                    "artifact binding differs from exact producer output"
                )
        seen[node.node_id] = node


__all__ = [
    "ArtifactInputIntentV1",
    "ArtifactOutputIntentV1",
    "ArtifactBindingV1",
    "ArtifactOutputV1",
    "CommandNodeIntentV1",
    "CommandNodeV1",
    "CommandWorkflowDraftV1",
    "CommandWorkflowSpecV1",
    "build_command_workflow_draft",
    "build_command_workflow_spec",
]
