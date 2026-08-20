"""Immutable multi-command workflow contracts; planning and preview only."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

from chemsmart.agent._contracts import (
    AuxiliaryArtifactBindingV1,
    ContractError,
    canonical_data,
    canonical_sha256,
    require_auxiliary_artifact_bindings,
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
        if produced and not (
            self.producer_node_id and self.producer_output_id
        ):
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


def _require_sorted_unique_in_node(
    values: tuple[str, ...], field: str, *, node_id: str
) -> None:
    """Refuse a malformed list while naming the node and the exact value.

    A rejection the caller cannot act on is barely better than a crash.  The
    old message said only "must be sorted and unique" -- no node, no value, and
    misleading besides, because the runtime sorts these lists before building
    the node, so duplication is the only reachable cause.  A live run hit it,
    retried unchanged, and produced the identical rejection.
    """

    duplicates = sorted(
        {value for value in values if list(values).count(value) > 1}
    )
    if duplicates:
        raise ContractError(
            f"node {node_id!r} repeats {field} {duplicates}; each must appear "
            "once, so give every entry a distinct name"
        )
    if tuple(values) != tuple(sorted(values)):
        raise ContractError(
            f"node {node_id!r} lists {field} out of order: {list(values)}; "
            f"expected {sorted(values)}"
        )


#: A DAG whose only node kind is "invoke a program" cannot express the stage
#: that turns finished results into the number the question asked for.  Without
#: the second kind, such an observable has no declared producer, so repairing a
#: plan to clear its findings also removes the answer.
WORKFLOW_NODE_KINDS = ("aggregate", "program_call")

#: The engine of an aggregate node.  ChemSmart performs the arithmetic itself,
#: so this is deliberately not a member of the executable program registry.
AGGREGATE_NODE_PROGRAM = "chemsmart"

#: The single stage name an aggregate node carries.  A post-processing step is
#: generally a small DAG of operations -- a complete-basis-set limit is several
#: -- so naming one operation here would misrepresent it.  The operations
#: themselves live in the quantity-expression vocabulary, which stays the one
#: source of truth; this node only declares that such a stage exists and what
#: it produces.
AGGREGATE_NODE_STAGE = "quantity_expression"


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
    #: ``program_call`` invokes a program; ``aggregate`` combines finished
    #: results into the requested number.  Defaulted, and omitted from the
    #: draft digest when default, so every already-recorded workflow keeps the
    #: identity it has and stays replayable.
    node_kind: str = "program_call"
    #: Optional node-local electronic state.  Most workflows inherit the
    #: task-bound state of their molecular input, so the fields stay absent.
    #: A producer geometry may, however, be deliberately reused on another
    #: charge/multiplicity surface (for example a vertical four-point energy).
    #: The pair is explicit and indivisible: geometry transfer never implies
    #: an electronic-state transfer.
    charge: int | None = None
    multiplicity: int | None = None
    #: Which internal coordinates this node drives or holds fixed.  Same class
    #: of fact as the state above -- a property of this molecule in this
    #: calculation rather than reusable method rationale -- so it lives here
    #: and not in the project.  Carried as canonical data so the draft digest
    #: covers it; the host renders it into each program's own idiom.
    internal_coordinates: Mapping[str, Any] | None = None

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        require_identifier(self.program, "program")
        self._validate_node_kind()
        self._validate_internal_coordinates()
        _validate_optional_electronic_state(
            self.charge,
            self.multiplicity,
            node_id=self.node_id,
        )
        require_identifier(self.jobtype, "jobtype")
        _require_sorted_unique_in_node(
            self.dependencies, "dependencies", node_id=self.node_id
        )
        _require_sorted_unique_in_node(
            tuple(item.binding_id for item in self.inputs),
            "input binding_id",
            node_id=self.node_id,
        )
        _require_sorted_unique_in_node(
            tuple(item.output_id for item in self.expected_outputs),
            "expected output_id",
            node_id=self.node_id,
        )
        if not self.expected_outputs:
            raise ContractError(
                f"node {self.node_id!r} declares no expected output; every "
                "node must state at least one artifact it produces"
            )
        _require_sorted_unique_in_node(
            self.unresolved_fields, "unresolved field", node_id=self.node_id
        )
        for field_name in self.unresolved_fields:
            try:
                require_identifier(field_name, "unresolved field")
            except ContractError as exc:
                raise ContractError(
                    f"node {self.node_id!r} has invalid unresolved field "
                    f"{field_name!r}; use a lower-case public identifier such "
                    "as 'input.geometry'"
                ) from exc

    def _validate_internal_coordinates(self) -> None:
        """Refuse a coordinate specification the jobtype cannot mean.

        A scan without a driven coordinate, or a driven coordinate on a job
        type that does not drive one, is a plan that reads as science but
        cannot compile.  Say so where it is written rather than at argv.
        """

        spec = self.internal_coordinates
        if spec is None:
            return
        if not isinstance(spec, Mapping):
            raise ContractError(
                f"node {self.node_id!r} internal_coordinates must be an "
                "object with 'scan' and/or 'constrained'"
            )
        unknown = sorted(set(spec).difference({"scan", "constrained"}))
        if unknown:
            raise ContractError(
                f"node {self.node_id!r} internal_coordinates does not accept "
                f"{unknown}; it accepts ['constrained', 'scan']"
            )
        has_scan = bool(spec.get("scan"))
        has_constrained = bool(spec.get("constrained"))
        if not has_scan and not has_constrained:
            raise ContractError(
                f"node {self.node_id!r} declares internal_coordinates with "
                "neither a scan nor a constraint; omit the field instead"
            )
        if self.jobtype == "scan" and not has_scan:
            raise ContractError(
                f"node {self.node_id!r} is a scan but names no driven "
                "coordinate; a scan is defined by what it drives"
            )
        if self.jobtype == "modred" and not has_constrained:
            raise ContractError(
                f"node {self.node_id!r} is a constrained optimisation but "
                "names no constrained coordinate"
            )
        if has_scan and self.jobtype not in {"scan", "ts"}:
            raise ContractError(
                f"node {self.node_id!r} names a driven coordinate but its "
                f"jobtype {self.jobtype!r} does not drive one; use 'scan'"
            )

    def _validate_node_kind(self) -> None:
        if self.node_kind not in WORKFLOW_NODE_KINDS:
            raise ContractError(
                f"node {self.node_id!r} has node_kind {self.node_kind!r}; "
                f"expected one of {sorted(WORKFLOW_NODE_KINDS)}"
            )
        if self.node_kind != "aggregate":
            if self.program == AGGREGATE_NODE_PROGRAM:
                raise ContractError(
                    f"node {self.node_id!r} names program "
                    f"{AGGREGATE_NODE_PROGRAM!r}, which is not an executable "
                    "program; use node_kind 'aggregate' for a stage that "
                    "combines finished results"
                )
            return
        # An aggregation is host arithmetic over finished results: it runs no
        # program, carries one canonical stage name, and must consume the
        # results it combines.
        if self.program != AGGREGATE_NODE_PROGRAM:
            raise ContractError(
                f"node {self.node_id!r} is an aggregate stage and runs no "
                f"program; its program must be {AGGREGATE_NODE_PROGRAM!r}, "
                f"not {self.program!r}"
            )
        if self.jobtype != AGGREGATE_NODE_STAGE:
            raise ContractError(
                f"node {self.node_id!r} is an aggregate stage; its jobtype "
                f"must be {AGGREGATE_NODE_STAGE!r}, not {self.jobtype!r}. The "
                "operations themselves belong to the quantity expression, not "
                "to the DAG node"
            )
        if not self.inputs:
            raise ContractError(
                f"node {self.node_id!r} is an aggregate stage but consumes "
                "nothing; it must declare the results it combines"
            )


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
            "nodes": tuple(
                _draft_node_digest_body(node) for node in self.nodes
            ),
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
    # ``node_kind`` joins the identity only when it is not the default, so a
    # workflow of ordinary program calls keeps the digest it has always had.
    digest_body = {
        **body,
        "nodes": tuple(_draft_node_digest_body(node) for node in nodes),
    }
    return CommandWorkflowDraftV1(
        **body, draft_sha256=canonical_sha256(digest_body)
    )


def _draft_node_digest_body(node: CommandNodeIntentV1) -> dict:
    body = canonical_data(node)
    if body.get("node_kind") == "program_call":
        body.pop("node_kind")
    if body.get("charge") is None and body.get("multiplicity") is None:
        # Additive state declarations must not invalidate existing workflow
        # identities when they are unused.
        body.pop("charge")
        body.pop("multiplicity")
    return body


def _validate_optional_electronic_state(
    charge: int | None,
    multiplicity: int | None,
    *,
    node_id: str,
) -> None:
    """Validate an explicit state without inventing one from geometry."""

    if (charge is None) != (multiplicity is None):
        raise ContractError(
            f"node {node_id!r} must declare charge and multiplicity together"
        )
    if charge is None:
        return
    if isinstance(charge, bool) or not isinstance(charge, int):
        raise ContractError(f"node {node_id!r} charge must be an integer")
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int):
        raise ContractError(
            f"node {node_id!r} multiplicity must be an integer"
        )
    if multiplicity < 1:
        raise ContractError(f"node {node_id!r} multiplicity must be positive")


def _validate_draft_dag(nodes: tuple[CommandNodeIntentV1, ...]) -> None:
    # A scientific workflow may be analysis-only when every root is an existing
    # registered result.  The model-facing calculation-only tool still requires
    # at least one node; permitting an empty internal draft avoids inventing a
    # fake program call solely to carry a real analysis DAG.
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
        if produced and not (
            self.producer_node_id and self.producer_output_id
        ):
            raise ContractError("produced artifact needs node and output IDs")
        if produced:
            require_identifier(self.producer_node_id, "producer_node_id")
            require_identifier(self.producer_output_id, "producer_output_id")
            if bool(self.artifact_sha256) != bool(
                self.producer_receipt_sha256
            ):
                raise ContractError(
                    "resolved producer artifact needs hash and producer receipt"
                )
        elif not self.artifact_sha256:
            raise ContractError(
                "external artifact binding requires an exact hash"
            )
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
            raise ContractError(
                "only previewed node carries preflight evidence"
            )


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
        require_sha256(self.live_cli_schema_sha256, "live_cli_schema_sha256")
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
                raise ContractError(
                    "two nodes declare the same output artifact"
                )
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
                raise ContractError(
                    "artifact binds an undeclared producer output"
                )
            if (
                binding.artifact_id != output.artifact_id
                or binding.artifact_class != output.artifact_class
            ):
                raise ContractError(
                    "artifact binding differs from exact producer output"
                )
        seen[node.node_id] = node


_COMPLEXITY_FACTORS = frozenset(
    {
        "conflicting_evidence",
        "excited_state",
        "gpu_engine",
        "multiple_stages",
        "producer_artifact_edge",
        "program_substitution",
        "solvent_semantics",
        "unresolved_scientific_setting",
    }
)


#: Support states a scientific workflow node may declare.
WORKFLOW_NODE_SUPPORT_STATES = (
    "blocked_unsupported",
    "resolvable",
    "unresolved_future",
)


@dataclass(frozen=True)
class ScientificWorkflowNodeV2:
    """Scientific stage before command materialization.

    A node names scientific intent and trusted project ownership only.  It
    deliberately carries neither argv nor a model-selected filesystem path.
    """

    node_id: str
    stage: str
    requested_program: str
    program: str
    engine: str
    project_role: str
    unresolved_fields: tuple[str, ...]
    #: Observables this node produces.  A node that carries the answer must say
    #: so, otherwise nothing can tell that dropping it loses the answer.
    produces_observables: tuple[str, ...] = ()
    #: ``resolvable`` -- materializable now; ``unresolved_future`` -- waiting on
    #: a producer artifact; ``blocked_unsupported`` -- required, expressible as
    #: intent, but outside the current program surface.  The third state exists
    #: so an honest plan can be finding-free without silently dropping a stage.
    support_state: str = "resolvable"
    blocked_reason: str = ""
    #: Explicit target state for this scientific stage.  When omitted the
    #: state is inherited from the task-bound input, preserving the historical
    #: plan meaning and digest.  When present on a data-edge consumer it
    #: authorizes a state rebind on the exact producer geometry.
    charge: int | None = None
    multiplicity: int | None = None

    def __post_init__(self) -> None:
        self._validate_identity()
        _validate_optional_electronic_state(
            self.charge,
            self.multiplicity,
            node_id=self.node_id,
        )
        if self.support_state not in WORKFLOW_NODE_SUPPORT_STATES:
            raise ContractError(
                "support_state must be one of "
                f"{sorted(WORKFLOW_NODE_SUPPORT_STATES)}"
            )
        if self.support_state == "blocked_unsupported":
            if not str(self.blocked_reason).strip():
                raise ContractError(
                    "a blocked node must state why it is unsupported"
                )
        elif str(self.blocked_reason).strip():
            raise ContractError(
                "blocked_reason applies only to a blocked_unsupported node"
            )
        observables = tuple(self.produces_observables)
        # Canonical JSON renders a tuple as a list, so a replayed record must be
        # normalized back or an identical plan would compare unequal.
        object.__setattr__(self, "produces_observables", observables)
        if observables != tuple(sorted(set(observables))):
            raise ContractError(
                "produces_observables must be sorted and unique"
            )
        for observable in self.produces_observables:
            require_identifier(observable, "observable")
        return

    def _validate_identity(self) -> None:
        for name, value in (
            ("node_id", self.node_id),
            ("stage", self.stage),
            ("requested_program", self.requested_program),
            ("program", self.program),
            ("engine", self.engine),
            ("project_role", self.project_role),
        ):
            require_identifier(value, name)
        if self.unresolved_fields != tuple(
            sorted(set(self.unresolved_fields))
        ):
            raise ContractError("unresolved_fields must be sorted and unique")
        for field_name in self.unresolved_fields:
            require_identifier(field_name, "unresolved field")


@dataclass(frozen=True)
class ScientificWorkflowEdgeV2:
    """A scheduling edge or an exact producer-artifact edge."""

    edge_id: str
    source_node_id: str
    target_node_id: str
    edge_kind: str
    artifact_class: str = ""
    producer_output_id: str = ""
    consumer_input_id: str = ""

    def __post_init__(self) -> None:
        require_identifier(self.edge_id, "edge_id")
        require_identifier(self.source_node_id, "source_node_id")
        require_identifier(self.target_node_id, "target_node_id")
        if self.source_node_id == self.target_node_id:
            raise ContractError("workflow edge cannot be a self edge")
        if self.edge_kind not in {"control", "data"}:
            raise ContractError("edge_kind must be control or data")
        artifact_fields = (
            self.artifact_class,
            self.producer_output_id,
            self.consumer_input_id,
        )
        if self.edge_kind == "data":
            if not all(artifact_fields):
                raise ContractError("data edge requires exact artifact roles")
            require_identifier(self.artifact_class, "artifact_class")
            require_identifier(self.producer_output_id, "producer_output_id")
            require_identifier(self.consumer_input_id, "consumer_input_id")
        elif any(artifact_fields):
            raise ContractError("control edge cannot claim an artifact")


@dataclass(frozen=True)
class ScientificWorkflowPlanV2:
    """Task-bound scientific DAG with explicit control and data edges."""

    schema_version: str
    workflow_id: str
    task_spec_sha256: str
    scientific_identity_sha256: str
    nodes: tuple[ScientificWorkflowNodeV2, ...]
    edges: tuple[ScientificWorkflowEdgeV2, ...]
    complexity_factors: tuple[str, ...]
    status: str
    plan_sha256: str
    #: Observables the task asked for.  Recording them on the plan is what lets
    #: a later repair be checked against the question rather than only against
    #: its own findings.
    required_observables: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-workflow-plan.v2":
            raise ContractError("unsupported scientific workflow plan schema")
        require_identifier(self.workflow_id, "workflow_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        require_sha256(
            self.scientific_identity_sha256, "scientific_identity_sha256"
        )
        _validate_scientific_dag(self.nodes, self.edges)
        derived = _derive_complexity_factors(
            self.nodes,
            self.edges,
            explicit=tuple(
                factor
                for factor in self.complexity_factors
                if factor in {"conflicting_evidence", "solvent_semantics"}
            ),
        )
        if self.complexity_factors != derived:
            raise ContractError(
                "complexity_factors differ from deterministic plan factors"
            )
        if self.status != "planned":
            raise ContractError("scientific workflow plan can only be planned")
        required = tuple(self.required_observables)
        object.__setattr__(self, "required_observables", required)
        if required != tuple(sorted(set(required))):
            raise ContractError(
                "required_observables must be sorted and unique"
            )
        produced = {
            observable
            for node in self.nodes
            for observable in node.produces_observables
        }
        missing = sorted(set(required) - produced)
        if missing:
            raise ContractError(
                "no node produces the required observable(s) "
                f"{missing}; a stage that cannot be materialized must stay in "
                "the plan as blocked_unsupported rather than be removed"
            )
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "task_spec_sha256": self.task_spec_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "nodes": tuple(
                _scientific_node_digest_body(node) for node in self.nodes
            ),
            "edges": self.edges,
            "complexity_factors": self.complexity_factors,
            "status": self.status,
        }
        if self.required_observables:
            # An absent field does not participate in identity, so plans that
            # declare no observable keep their existing digest.
            body["required_observables"] = self.required_observables
        if self.plan_sha256 != canonical_sha256(body):
            raise ContractError("scientific workflow plan digest mismatch")


def build_scientific_workflow_plan(
    *,
    workflow_id: str,
    task_spec_sha256: str,
    scientific_identity_sha256: str,
    nodes: Sequence[ScientificWorkflowNodeV2],
    edges: Sequence[ScientificWorkflowEdgeV2] = (),
    explicit_complexity_factors: Sequence[str] = (),
    required_observables: Sequence[str] = (),
) -> ScientificWorkflowPlanV2:
    normalized_nodes = tuple(nodes)
    normalized_edges = tuple(edges)
    _validate_scientific_dag(normalized_nodes, normalized_edges)
    explicit = tuple(sorted(set(explicit_complexity_factors)))
    if set(explicit).difference({"conflicting_evidence", "solvent_semantics"}):
        raise ContractError("unsupported explicit complexity factor")
    body = {
        "schema_version": "chemsmart.scientific-workflow-plan.v2",
        "workflow_id": require_identifier(workflow_id, "workflow_id"),
        "task_spec_sha256": require_sha256(
            task_spec_sha256, "task_spec_sha256"
        ),
        "scientific_identity_sha256": require_sha256(
            scientific_identity_sha256, "scientific_identity_sha256"
        ),
        "nodes": normalized_nodes,
        "edges": normalized_edges,
        "complexity_factors": _derive_complexity_factors(
            normalized_nodes, normalized_edges, explicit=explicit
        ),
        "status": "planned",
    }
    observables = tuple(
        sorted(set(str(item) for item in required_observables))
    )
    digest_body = dict(body)
    digest_body["nodes"] = tuple(
        _scientific_node_digest_body(node) for node in normalized_nodes
    )
    if observables:
        digest_body["required_observables"] = observables
    return ScientificWorkflowPlanV2(
        **body,
        required_observables=observables,
        plan_sha256=canonical_sha256(digest_body),
    )


def _scientific_node_digest_body(node: ScientificWorkflowNodeV2) -> dict:
    """Canonical node projection with legacy-compatible optional state."""

    body = canonical_data(node)
    if body.get("charge") is None and body.get("multiplicity") is None:
        body.pop("charge")
        body.pop("multiplicity")
    return body


def _validate_scientific_dag(
    nodes: tuple[ScientificWorkflowNodeV2, ...],
    edges: tuple[ScientificWorkflowEdgeV2, ...],
) -> None:
    if not nodes:
        raise ContractError("scientific workflow must contain a node")
    node_ids = tuple(node.node_id for node in nodes)
    if len(node_ids) != len(set(node_ids)):
        raise ContractError("scientific workflow node IDs must be unique")
    if tuple(edge.edge_id for edge in edges) != tuple(
        sorted({edge.edge_id for edge in edges})
    ):
        raise ContractError(
            "scientific workflow edges must be sorted and unique"
        )
    positions = {node_id: index for index, node_id in enumerate(node_ids)}
    data_inputs: set[tuple[str, str]] = set()
    for edge in edges:
        if (
            edge.source_node_id not in positions
            or edge.target_node_id not in positions
        ):
            raise ContractError(
                "scientific workflow edge names an unknown node"
            )
        if positions[edge.source_node_id] >= positions[edge.target_node_id]:
            raise ContractError(
                "scientific workflow must be topologically ordered"
            )
        if edge.edge_kind == "data":
            target_input = (edge.target_node_id, edge.consumer_input_id)
            if target_input in data_inputs:
                raise ContractError(
                    "consumer input has multiple artifact producers"
                )
            data_inputs.add(target_input)


def _derive_complexity_factors(
    nodes: tuple[ScientificWorkflowNodeV2, ...],
    edges: tuple[ScientificWorkflowEdgeV2, ...],
    *,
    explicit: tuple[str, ...],
) -> tuple[str, ...]:
    factors = set(explicit)
    if len(nodes) > 1:
        factors.add("multiple_stages")
    if any(edge.edge_kind == "data" for edge in edges):
        factors.add("producer_artifact_edge")
    if any(node.engine == "gpu" for node in nodes):
        factors.add("gpu_engine")
    if any(node.stage in {"td", "tda", "tddft"} for node in nodes):
        factors.add("excited_state")
    if any(node.requested_program != node.program for node in nodes):
        factors.add("program_substitution")
    if any(node.unresolved_fields for node in nodes):
        factors.add("unresolved_scientific_setting")
    if factors.difference(_COMPLEXITY_FACTORS):
        raise ContractError("unsupported workflow complexity factor")
    return tuple(sorted(factors))


@dataclass(frozen=True)
class MaterializedNodeV1:
    """Host-grounded receipts for one scientific stage."""

    node_id: str
    input_artifact_sha256: str
    project_artifact_sha256: str
    project_validation_receipt_sha256: str
    environment_receipt_sha256: str
    invocation_sha256: str
    preflight_receipt_sha256: str
    state: str
    #: Path-independent identity of the compiled command.  ``invocation_sha256``
    #: covers argv, which carries absolute paths into this session's own
    #: timestamped run directory, so no later process can reproduce it.
    #: Omitted from the workflow digest when empty, so a materialization
    #: recorded before this field existed keeps the digest it was stored with.
    invocation_identity_sha256: str = ""
    auxiliary_input_bindings: tuple[AuxiliaryArtifactBindingV1, ...] = ()

    def __post_init__(self) -> None:
        require_identifier(self.node_id, "node_id")
        if self.invocation_identity_sha256:
            require_sha256(
                self.invocation_identity_sha256, "invocation_identity_sha256"
            )
        require_auxiliary_artifact_bindings(self.auxiliary_input_bindings)
        for name, digest in (
            ("input_artifact_sha256", self.input_artifact_sha256),
            ("project_artifact_sha256", self.project_artifact_sha256),
            (
                "project_validation_receipt_sha256",
                self.project_validation_receipt_sha256,
            ),
            ("environment_receipt_sha256", self.environment_receipt_sha256),
        ):
            require_sha256(digest, name)
        if self.state not in {"grounded", "compiled", "previewed"}:
            raise ContractError("unsupported node materialization state")
        if self.state in {"compiled", "previewed"}:
            require_sha256(self.invocation_sha256, "invocation_sha256")
        elif self.invocation_sha256:
            raise ContractError("grounded node cannot carry an invocation")
        if self.state == "previewed":
            require_sha256(
                self.preflight_receipt_sha256, "preflight_receipt_sha256"
            )
        elif self.preflight_receipt_sha256:
            raise ContractError(
                "only previewed node carries preflight evidence"
            )


def _materialized_node_body(node: MaterializedNodeV1) -> dict[str, Any]:
    """Canonical node projection that omits unset optional evidence.

    The workflow digest folds its nodes, so serialising the dataclass whole
    would have made every stored ``materialized_sha256`` invalid the moment an
    optional field was added.  Omitting an empty optional keeps a workflow
    recorded before the field readable, and a workflow that carries the field
    digests it.
    """

    body = {
        "node_id": node.node_id,
        "input_artifact_sha256": node.input_artifact_sha256,
        "project_artifact_sha256": node.project_artifact_sha256,
        "project_validation_receipt_sha256": (
            node.project_validation_receipt_sha256
        ),
        "environment_receipt_sha256": node.environment_receipt_sha256,
        "invocation_sha256": node.invocation_sha256,
        "preflight_receipt_sha256": node.preflight_receipt_sha256,
        "state": node.state,
    }
    if node.invocation_identity_sha256:
        body["invocation_identity_sha256"] = node.invocation_identity_sha256
    if node.auxiliary_input_bindings:
        body["auxiliary_input_bindings"] = node.auxiliary_input_bindings
    return body


@dataclass(frozen=True)
class MaterializedWorkflowV1:
    """Partially or fully grounded view of a scientific workflow plan."""

    schema_version: str
    workflow_id: str
    plan_sha256: str
    task_spec_sha256: str
    scientific_identity_sha256: str
    live_cli_schema_sha256: str
    resource_sha256: str
    nodes: tuple[MaterializedNodeV1, ...]
    unresolved_node_ids: tuple[str, ...]
    status: str
    materialized_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.materialized-workflow.v1":
            raise ContractError("unsupported materialized workflow schema")
        require_identifier(self.workflow_id, "workflow_id")
        for name, digest in (
            ("plan_sha256", self.plan_sha256),
            ("task_spec_sha256", self.task_spec_sha256),
            ("scientific_identity_sha256", self.scientific_identity_sha256),
            ("live_cli_schema_sha256", self.live_cli_schema_sha256),
            ("resource_sha256", self.resource_sha256),
        ):
            require_sha256(digest, name)
        node_ids = tuple(node.node_id for node in self.nodes)
        if node_ids != tuple(sorted(set(node_ids))):
            raise ContractError("materialized nodes must be sorted and unique")
        if self.unresolved_node_ids != tuple(
            sorted(set(self.unresolved_node_ids))
        ):
            raise ContractError(
                "unresolved node IDs must be sorted and unique"
            )
        if set(node_ids).intersection(self.unresolved_node_ids):
            raise ContractError("node cannot be materialized and unresolved")
        for node_id in self.unresolved_node_ids:
            require_identifier(node_id, "unresolved_node_id")
        if self.status not in {
            "partial",
            "materialized",
            "previewed",
            "ready_for_approval",
        }:
            raise ContractError("unsupported materialized workflow status")
        if self.status == "partial" and not self.unresolved_node_ids:
            raise ContractError("partial workflow requires unresolved nodes")
        if self.status != "partial" and self.unresolved_node_ids:
            raise ContractError(
                "non-partial workflow cannot have unresolved nodes"
            )
        if self.status in {"previewed", "ready_for_approval"} and not all(
            node.state == "previewed" for node in self.nodes
        ):
            raise ContractError(
                "preview-ready workflow requires previewed nodes"
            )
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "plan_sha256": self.plan_sha256,
            "task_spec_sha256": self.task_spec_sha256,
            "scientific_identity_sha256": self.scientific_identity_sha256,
            "live_cli_schema_sha256": self.live_cli_schema_sha256,
            "resource_sha256": self.resource_sha256,
            "nodes": tuple(
                _materialized_node_body(node) for node in self.nodes
            ),
            "unresolved_node_ids": self.unresolved_node_ids,
            "status": self.status,
        }
        if self.materialized_sha256 != canonical_sha256(body):
            raise ContractError("materialized workflow digest mismatch")


#: The ``resource_sha256`` a workflow carries when it was materialized without
#: execution resources -- that is, by ``agent plan``, which makes no engine
#: call.  It is a sentinel meaning "not resource-bound yet", not a profile, and
#: the freeze step is what binds real resources to a reviewed plan.
PREVIEW_RESOURCE_SHA256 = canonical_sha256(
    {
        "schema_version": "chemsmart.preview-resource.v1",
        "chemistry_engine_calls": 0,
    }
)


def build_materialized_workflow(
    *,
    plan: ScientificWorkflowPlanV2,
    live_cli_schema_sha256: str,
    resource_sha256: str,
    nodes: Sequence[MaterializedNodeV1],
    unresolved_node_ids: Sequence[str],
    status: str,
) -> MaterializedWorkflowV1:
    normalized_nodes = tuple(sorted(nodes, key=lambda node: node.node_id))
    unresolved = tuple(sorted(set(unresolved_node_ids)))
    planned_ids = {node.node_id for node in plan.nodes}
    supplied_ids = {node.node_id for node in normalized_nodes}
    if supplied_ids.union(unresolved) != planned_ids:
        raise ContractError("materialization must account for every plan node")
    body = {
        "schema_version": "chemsmart.materialized-workflow.v1",
        "workflow_id": plan.workflow_id,
        "plan_sha256": plan.plan_sha256,
        "task_spec_sha256": plan.task_spec_sha256,
        "scientific_identity_sha256": plan.scientific_identity_sha256,
        "live_cli_schema_sha256": require_sha256(
            live_cli_schema_sha256, "live_cli_schema_sha256"
        ),
        "resource_sha256": require_sha256(resource_sha256, "resource_sha256"),
        "nodes": normalized_nodes,
        "unresolved_node_ids": unresolved,
        "status": str(status),
    }
    # The digest is taken over the same node projection the validator uses.
    # Two spellings of one body is how a record becomes unconstructible.
    digest_body = dict(body)
    digest_body["nodes"] = tuple(
        _materialized_node_body(node) for node in normalized_nodes
    )
    return MaterializedWorkflowV1(
        **body, materialized_sha256=canonical_sha256(digest_body)
    )


@dataclass(frozen=True)
class StationaryPointValidationPolicyV1:
    schema_version: str
    policy_id: str
    task_spec_sha256: str
    hessian_node_id: str
    stationary_point_kind: str
    expected_imaginary_mode_count: int
    imaginary_mode_cutoff_cm1: float
    require_finite_modes: bool
    require_symmetric_hessian: bool
    policy_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.stationary-point-policy.v1":
            raise ContractError("unsupported stationary point policy schema")
        require_identifier(self.policy_id, "policy_id")
        require_identifier(self.hessian_node_id, "hessian_node_id")
        require_sha256(self.task_spec_sha256, "task_spec_sha256")
        if self.stationary_point_kind not in {"minimum", "transition_state"}:
            raise ContractError("unsupported stationary point kind")
        expected = 0 if self.stationary_point_kind == "minimum" else 1
        if self.expected_imaginary_mode_count != expected:
            raise ContractError(
                "imaginary-mode expectation differs from point kind"
            )
        if self.imaginary_mode_cutoff_cm1 <= 0:
            raise ContractError("imaginary-mode cutoff must be positive")
        body = {
            "schema_version": self.schema_version,
            "policy_id": self.policy_id,
            "task_spec_sha256": self.task_spec_sha256,
            "hessian_node_id": self.hessian_node_id,
            "stationary_point_kind": self.stationary_point_kind,
            "expected_imaginary_mode_count": self.expected_imaginary_mode_count,
            "imaginary_mode_cutoff_cm1": self.imaginary_mode_cutoff_cm1,
            "require_finite_modes": self.require_finite_modes,
            "require_symmetric_hessian": self.require_symmetric_hessian,
        }
        if self.policy_sha256 != canonical_sha256(body):
            raise ContractError("stationary point policy digest mismatch")


__all__ = [
    "ArtifactInputIntentV1",
    "ArtifactOutputIntentV1",
    "ArtifactBindingV1",
    "ArtifactOutputV1",
    "CommandNodeIntentV1",
    "CommandNodeV1",
    "CommandWorkflowDraftV1",
    "CommandWorkflowSpecV1",
    "MaterializedNodeV1",
    "MaterializedWorkflowV1",
    "ScientificWorkflowEdgeV2",
    "ScientificWorkflowNodeV2",
    "ScientificWorkflowPlanV2",
    "StationaryPointValidationPolicyV1",
    "build_command_workflow_draft",
    "build_command_workflow_spec",
    "build_materialized_workflow",
    "build_scientific_workflow_plan",
]
