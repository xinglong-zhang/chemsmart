"""Execute an approved scientific workflow without consulting a model.

Why this exists
---------------

``agent run --approval-file`` is ``run_live_agent_session`` with one extra
tool, so a model holding a digest-frozen plan is asked to re-derive it turn by
turn.  Measured over the campaign's live sessions, that is where the harness
stops working: seventeen ``agent plan`` sessions rejected 42 of 868 tool calls
(4.8%), while the one ``agent run`` session rejected 20 of 69 (29%) -- twelve
of them attempts to execute, eight of them the model re-planning inside an
execution session -- and completed no nodes at all.

Nothing in that work needs judgement.  Which project, which program, which
node graph, which charge and multiplicity: all of it was decided by the model
at plan time and frozen under ``plan_sha256``.  What remains is bookkeeping
with exactly one right answer -- which node is ready, which digest satisfies
which contract, which receipt the host itself just produced -- and the harness
already owns a deterministic answer to each.

So this module supplies the caller, not the science.  It composes the same
host ``agent run`` composes, and calls the same ``host.dispatch`` a model
calls, with host-computed arguments.  The event log, the receipt chain and the
YAML->CLI DAG are the same records either way; only the author of the
bookkeeping changes.  That strengthens the provenance claim rather than
weakening it, because unlike a model this executor *cannot* silently re-plan:
``derive_ready_node_ids`` and the frozen approval both forbid deviation.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.execution import (
    build_workflow_run_state,
    derive_ready_node_ids,
    transition_workflow_node,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import build_approved_execution_tool_surface

#: The exact host calls one program node needs, in order.  This sequence is
#: not a heuristic: it is the order the tool contracts already require, and it
#: is what the model in the ``agent run`` session attempted and was refused
#: twelve times.
PROGRAM_NODE_SEQUENCE = (
    "inspect_program_capability",
    "inspect_program_environment",
    "bind_scientific_identity",
    "validate_project_yaml",
    "synthesize_command",
    "preview_command",
    "preflight_program_node",
    "execute_approved_program_node",
)


@dataclass(frozen=True)
class ExecutedNodeV1:
    """What one node did, reported without interpretation."""

    node_id: str
    program: str
    jobtype: str
    state: str
    invocation_identity_sha256: str
    execution_receipt_sha256: str
    rule_ids: tuple[str, ...]
    failure: str
    validated: bool = False
    result_validation_receipt_sha256: str = ""
    invocation_sha256: str = ""


@dataclass(frozen=True)
class WorkflowExecutionResultV1:
    """Outcome of one host-driven walk over an approved plan."""

    workflow_id: str
    plan_sha256: str
    approval_sha256: str
    run_directory: str
    nodes: tuple[ExecutedNodeV1, ...]
    status: str
    provider_calls: int = 0

    @property
    def executed_node_ids(self) -> tuple[str, ...]:
        return tuple(
            node.node_id for node in self.nodes if node.state == "executed"
        )


def _result_of(payload: Any) -> Any:
    """Unwrap a dispatch envelope; handlers return typed objects or dicts."""

    if isinstance(payload, Mapping) and "result" in payload:
        return payload["result"]
    return payload


def _field(value: Any, *names: str) -> Any:
    """Read a field whether the host handed back a record or a mapping."""

    for name in names:
        if isinstance(value, Mapping) and name in value:
            return value[name]
        if hasattr(value, name):
            return getattr(value, name)
    raise ContractError(
        f"host result does not carry any of {names!r}; the executor and the "
        "tool contract have drifted apart"
    )


def _locate_by_digest(workspace: Path, sha256: str) -> Path:
    """Find an approved artifact by content, never by remembered path.

    An approval names bytes.  Where those bytes sat during planning is a
    property of that session's directory layout, so the file is resolved here
    the way the approval identifies it.
    """

    for candidate in sorted(workspace.rglob("*")):
        if not candidate.is_file() or candidate.is_symlink():
            continue
        try:
            if file_sha256(candidate) == sha256:
                return candidate
        except OSError:
            continue
    raise ContractError(
        f"no file under {workspace} has the approved digest {sha256[:16]}...; "
        "the approved input is not present in this workspace"
    )


def _approved_initial_artifacts(
    workspace: Path, approval: Any
) -> dict[str, TrustedArtifactRefV1]:
    """Resolve every independently approved workflow entry geometry.

    A comparison workflow can have several sibling roots: two charge states,
    reactant and product, or multiple conformers.  Each node must receive the
    artifact named by its own approval binding rather than whichever initial
    geometry happened to appear first.
    """

    artifacts: dict[str, TrustedArtifactRefV1] = {}
    for binding in approval.node_bindings:
        if binding.input_mode != "initial":
            continue
        existing = artifacts.get(binding.initial_artifact_id)
        if existing is not None:
            if existing.sha256 != binding.initial_artifact_sha256:
                raise ContractError(
                    "one approved initial artifact ID names different bytes"
                )
            continue
        geometry_path = _locate_by_digest(
            workspace, binding.initial_artifact_sha256
        )
        artifacts[binding.initial_artifact_id] = TrustedArtifactRefV1(
            artifact_id=binding.initial_artifact_id,
            kind="geometry_xyz",
            sha256=binding.initial_artifact_sha256,
            size_bytes=geometry_path.stat().st_size,
            path=str(geometry_path),
            cli_value=str(geometry_path),
        )
    if not artifacts:
        raise ContractError(
            "no approved node takes an initial geometry, so the workflow has "
            "no entry point"
        )
    return artifacts


class ApprovedWorkflowExecutor:
    """Walk an approved DAG, dispatching host tools with host-computed args."""

    def __init__(
        self,
        *,
        host: CommandCompiledToolHostV1,
        plan: Any,
        approval: Any,
        frozen_approval: Any,
        initial_artifacts: Mapping[str, TrustedArtifactRefV1],
        project_artifacts: Mapping[str, Any],
        task_spec_sha256: str,
        run_directory: Path,
    ) -> None:
        self.host = host
        self.plan = plan
        self.approval = approval
        self.frozen_approval = frozen_approval
        self.initial_artifacts = dict(initial_artifacts)
        self.project_by_digest = {
            item.sha256: item.artifact_id for item in project_artifacts
        }
        self.task_spec_sha256 = task_spec_sha256
        self.run_directory = run_directory
        self._turn = 0
        self._handoff_inputs: dict[str, str] = {}

    def _call(self, tool_name: str, **arguments: Any) -> Any:
        self._turn += 1
        return _result_of(
            self.host.dispatch(
                turn_id=f"exec-{self._turn:04d}",
                tool_name=tool_name,
                arguments=arguments,
            )
        )

    def _binding(self, node_id: str) -> Any:
        for binding in self.approval.node_bindings:
            if binding.node_id == node_id:
                return binding
        raise ContractError(
            f"node {node_id!r} is in the approved plan but has no approved "
            "binding; the approval and the plan disagree"
        )

    def _input_artifact_id(self, binding: Any) -> str:
        if binding.input_mode == "initial":
            artifact = self.initial_artifacts.get(
                binding.initial_artifact_id
            )
            if artifact is None:
                raise ContractError(
                    f"approved initial artifact {binding.initial_artifact_id!r} "
                    "is unavailable"
                )
            if artifact.sha256 != binding.initial_artifact_sha256:
                raise ContractError(
                    "approved initial artifact bytes differ from node binding"
                )
            return artifact.artifact_id
        handoff = self._handoff_inputs.get(binding.node_id)
        if not handoff:
            raise ContractError(
                f"node {binding.node_id!r} consumes a producer edge, but no "
                "upstream node has handed off a validated artifact yet"
            )
        return handoff

    def run_node(self, node_id: str) -> ExecutedNodeV1:
        """Drive one approved node through the fixed host sequence."""

        binding = self._binding(node_id)
        identity_sha256 = ""
        invocation_sha256 = ""
        try:
            capability = self._call(
                "inspect_program_capability",
                program=binding.program,
                jobtype=binding.jobtype,
                engine=binding.engine,
            )
            capability_sha256 = _field(capability, "receipt_sha256")
            environment = self._call(
                "inspect_program_environment",
                capability_receipt_sha256=capability_sha256,
            )
            program_binding_sha256 = _field(
                _field(environment, "program_binding"), "binding_sha256"
            )
            engine_binding_sha256 = _field(
                _field(environment, "engine_binding"), "binding_sha256"
            )

            input_artifact_id = self._input_artifact_id(binding)
            input_artifact = self.host.artifacts.get(input_artifact_id)
            if input_artifact is None:
                raise ContractError(
                    f"input artifact {input_artifact_id!r} is unavailable"
                )
            scientific_identity = self._call(
                "bind_scientific_identity",
                input_artifact_id=input_artifact_id,
                task_spec_sha256=self.task_spec_sha256,
                charge=binding.charge,
                multiplicity=binding.multiplicity,
            )
            scientific_identity_sha256 = _field(
                scientific_identity, "binding_sha256"
            )

            project_artifact_id = self.project_by_digest[
                binding.project_artifact_sha256
            ]
            self._call(
                "validate_project_yaml",
                project_artifact_id=project_artifact_id,
                capability_receipt_sha256=capability_sha256,
            )

            synthesized = self._call(
                "synthesize_command",
                node_id=node_id,
                program=binding.program,
                jobtype=binding.jobtype,
                project_artifact_id=project_artifact_id,
                input_artifact_id=input_artifact_id,
                scientific_identity_sha256=scientific_identity_sha256,
                charge=binding.charge,
                multiplicity=binding.multiplicity,
                capability_receipt_sha256=capability_sha256,
                engine_binding_sha256=engine_binding_sha256,
            )
            invocation_sha256 = _field(
                _field(synthesized, "invocation"), "invocation_sha256"
            )
            inspection_sha256 = _field(
                _field(synthesized, "inspection"), "receipt_sha256"
            )
            identity_sha256 = self.host._invocation_identity(node_id)

            self._call("preview_command", invocation_sha256=invocation_sha256)
            self._call(
                "preflight_program_node",
                node_id=node_id,
                capability_receipt_sha256=capability_sha256,
                program_binding_sha256=program_binding_sha256,
                engine_binding_sha256=engine_binding_sha256,
                geometry_artifact_sha256=(
                    input_artifact.sha256
                ),
                scientific_identity_sha256=scientific_identity_sha256,
                charge=binding.charge,
                multiplicity=binding.multiplicity,
                invocation_sha256=invocation_sha256,
                command_inspection_receipt_sha256=inspection_sha256,
            )

            executed = self._call(
                "execute_approved_program_node", node_id=node_id
            )
            receipt = _field(executed, "execution")
            for item in _field(executed, "produced_handoffs") or ():
                handoff = (
                    item["handoff"] if isinstance(item, Mapping) else item
                )
                artifact = _field(handoff, "geometry_artifact_id")
                self._handoff_inputs[
                    _field(handoff, "consumer_node_id")
                ] = artifact
            return ExecutedNodeV1(
                node_id=node_id,
                program=binding.program,
                jobtype=binding.jobtype,
                state=_field(receipt, "execution_state"),
                invocation_identity_sha256=identity_sha256,
                execution_receipt_sha256=_field(receipt, "receipt_sha256"),
                rule_ids=tuple(_field(receipt, "findings") or ()),
                failure="",
                validated=bool(_field(receipt, "validated")),
                result_validation_receipt_sha256=_field(
                    receipt, "result_validation_receipt_sha256"
                ),
                invocation_sha256=_field(receipt, "invocation_sha256"),
            )
        except ContractError as error:
            return ExecutedNodeV1(
                node_id=node_id,
                program=binding.program,
                jobtype=binding.jobtype,
                state="blocked",
                invocation_identity_sha256=identity_sha256,
                execution_receipt_sha256="",
                rule_ids=(),
                failure=str(error),
                invocation_sha256=invocation_sha256,
            )

    def _settle(self, run_state: Any, outcome: ExecutedNodeV1) -> Any:
        """Move a node to the state its own execution receipt reports.

        The host already transitioned its event-store copy; this keeps the
        local frontier in step, and refuses to invent a state the receipt did
        not claim -- a consumer becomes ready only after its producer is
        genuinely ``validated``.
        """

        stamp = datetime.now(timezone.utc).isoformat()
        if not outcome.invocation_sha256:
            # Nothing was compiled, so the node never left the frontier in a
            # runnable state; record why rather than inventing a transition.
            return transition_workflow_node(
                run_state,
                node_id=outcome.node_id,
                new_state="blocked",
                timestamp=stamp,
            )
        run_state = transition_workflow_node(
            run_state,
            node_id=outcome.node_id,
            new_state="running",
            invocation_sha256=outcome.invocation_sha256,
            timestamp=stamp,
        )
        if not outcome.validated:
            return transition_workflow_node(
                run_state,
                node_id=outcome.node_id,
                new_state="failed",
                execution_receipt_sha256=outcome.execution_receipt_sha256,
                failure_rule_ids=outcome.rule_ids
                or ("execution.state." + (outcome.state or "unknown"),),
                timestamp=stamp,
            )
        run_state = transition_workflow_node(
            run_state,
            node_id=outcome.node_id,
            new_state="engine_complete",
            execution_receipt_sha256=outcome.execution_receipt_sha256,
            timestamp=stamp,
        )
        receipt = self.host.result_validation_receipts.get(
            outcome.result_validation_receipt_sha256
        )
        if receipt is None:
            raise ContractError(
                f"node {outcome.node_id!r} reported a validated execution but "
                "the host holds no result validation receipt for it"
            )
        return transition_workflow_node(
            run_state,
            node_id=outcome.node_id,
            new_state="validated",
            execution_receipt_sha256=outcome.execution_receipt_sha256,
            result_validation_receipt=receipt,
            timestamp=stamp,
        )

    def run(self) -> WorkflowExecutionResultV1:
        """Execute every ready node until the frontier stops advancing."""

        run_id = "run." + self.approval.approval_id
        run_state = build_workflow_run_state(
            run_id=run_id,
            plan=self.plan,
            approval=self.frozen_approval,
            # This walk is the act of consuming the approval. An unconsumed
            # approval has an empty frontier by construction, so leaving this
            # false would report "nothing to do" for a fully approved plan.
            approval_consumed=True,
        )
        executed: list[ExecutedNodeV1] = []
        seen: set[str] = set()
        data_edge_bindings = ()
        while True:
            ready = tuple(
                node_id
                for node_id in derive_ready_node_ids(
                    self.plan, run_state, data_edge_bindings
                )
                if node_id not in seen
            )
            if not ready:
                break
            progressed = False
            for node_id in ready:
                seen.add(node_id)
                outcome = self.run_node(node_id)
                executed.append(outcome)
                if outcome.execution_receipt_sha256:
                    frontier = self.host.event_store.workflow_frontier(
                        workflow_id=self.plan.workflow_id,
                        run_id=run_id,
                    )
                    if frontier.run_state is None:
                        raise ContractError(
                            "executed node has no durable workflow run state"
                        )
                    run_state = frontier.run_state
                    data_edge_bindings = frontier.data_edge_bindings
                else:
                    run_state = self._settle(run_state, outcome)
                progressed = progressed or outcome.validated
            if not progressed:
                break
        planned = {node.node_id for node in self.plan.nodes}
        done = {item.node_id for item in executed if item.validated}
        status = "completed" if done == planned else "partial"
        return WorkflowExecutionResultV1(
            workflow_id=self.plan.workflow_id,
            plan_sha256=self.plan.plan_sha256,
            approval_sha256=self.frozen_approval.approval_sha256,
            run_directory=str(self.run_directory),
            nodes=tuple(executed),
            status=status,
            provider_calls=0,
        )


def execute_approved_workflow(
    *,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    task_spec_sha256: str = "",
) -> WorkflowExecutionResultV1:
    """Execute an approved workflow end to end, contacting no provider.

    Every scientific decision is read from the approval; this function chooses
    nothing except the order in which already-approved nodes become ready,
    which ``derive_ready_node_ids`` determines from the plan's own edges.
    """

    from chemsmart.agent.live_session import (
        _bootstrap_conformance,
        _execution_composition_inputs,
        _observe_environments,
    )

    workspace = Path(workspace).resolve()
    run_directory = Path(run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)

    inputs = _execution_composition_inputs(
        host_type=CommandCompiledToolHostV1,
        workspace=workspace,
        run_directory=run_directory,
        approval_file=Path(approval_file).resolve(),
        task_spec_sha256=task_spec_sha256,
    )
    plan = inputs.pop("approved_scientific_plan")
    if plan is None:
        raise ContractError(
            "approval file carries no approved scientific plan, so there is "
            "no DAG to execute; freeze the plan session's plan into it"
        )
    materialized = inputs.pop("approved_materialized_workflow", None)
    environment_identities = inputs.pop("approved_environment_identities", ())
    project_artifacts = inputs.pop("approved_project_artifacts")
    approval = inputs["workflow_execution_approval"]
    frozen_approval = inputs["frozen_workflow_approval"]

    initial_artifacts = _approved_initial_artifacts(workspace, approval)
    bootstrap_artifact = next(iter(initial_artifacts.values()))

    registry = load_program_capabilities()
    live_schema = build_live_click_schema()
    conformance, _records = _bootstrap_conformance(
        run_directory=run_directory,
        input_artifact=bootstrap_artifact,
        registry_sha256=registry.registry_sha256,
        live_schema=live_schema,
    )
    environment_targets, compute_receipts, _observed = _observe_environments()

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            run_directory / "events.jsonl",
            session_id=f"execute-{plan.workflow_id}",
        ),
        artifacts={
            **{item.artifact_id: item for item in project_artifacts},
            **initial_artifacts,
        },
        component_conformance_receipts=conformance,
        environment_targets=environment_targets,
        compute_environment_receipts=compute_receipts,
        tool_surface=build_approved_execution_tool_surface(registry),
        registry=registry,
        live_schema=live_schema,
        task_spec_sha256s=(task_spec_sha256,) if task_spec_sha256 else (),
        # ``_execution_composition_inputs`` already writes and supplies the
        # execution server profile; naming it again here would be a second
        # source of truth for where real jobs are launched.
        scientific_workflow_plan=plan,
        materialized_workflow=materialized,
        approved_environment_identities=environment_identities,
        **inputs,
    )
    return ApprovedWorkflowExecutor(
        host=host,
        plan=plan,
        approval=approval,
        frozen_approval=frozen_approval,
        initial_artifacts=initial_artifacts,
        project_artifacts=project_artifacts,
        task_spec_sha256=task_spec_sha256,
        run_directory=run_directory,
    ).run()


__all__ = [
    "ApprovedWorkflowExecutor",
    "ExecutedNodeV1",
    "PROGRAM_NODE_SEQUENCE",
    "WorkflowExecutionResultV1",
    "execute_approved_workflow",
]
