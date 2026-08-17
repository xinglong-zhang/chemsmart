"""Execute one exact human-approved scientific workflow without a provider.

The executor follows the frozen DAG and normal ChemSmart host tools.  It does
not choose science, re-plan, or widen approval.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.agent.api_access import (
    DEFAULT_KEY_LABELS,
    PROVIDER_KEY_LABEL_TOKENS,
    normalize_key_label,
)
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.execution import (
    WorkflowExecutionApprovalBundleV1,
    WorkflowExecutionReviewV1,
    approve_workflow_execution_review,
    build_workflow_run_state,
    derive_ready_node_ids,
    transition_workflow_node,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.tool_specs import build_approved_execution_tool_surface

#: The present-tense host contract for one approved program node.
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


def _provider_secret_environment_labels() -> tuple[str, ...]:
    """Return known and provider-shaped credential labels in this process.

    Profiles may use lab- or project-specific ``api_key_env`` names.  The
    provider token is deliberately required in those labels, so engine
    isolation can discover and remove them without retaining provider config
    in the approved executor.
    """

    return tuple(
        sorted(
            {
                label
                for labels in DEFAULT_KEY_LABELS.values()
                for label in labels
            }
            | {
                label
                for label in os.environ
                if any(
                    token in normalize_key_label(label)
                    for token in PROVIDER_KEY_LABEL_TOKENS.values()
                )
            }
        )
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
    non_executable_node_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "non_executable_node_ids",
            tuple(self.non_executable_node_ids),
        )
        if self.non_executable_node_ids != tuple(
            sorted(set(self.non_executable_node_ids))
        ):
            raise ContractError(
                "non-executable result node ids must be sorted and unique"
            )
        if set(self.non_executable_node_ids).intersection(
            node.node_id for node in self.nodes
        ):
            raise ContractError(
                "a non-executable workflow node cannot have an execution result"
            )

    @property
    def executed_node_ids(self) -> tuple[str, ...]:
        return tuple(node.node_id for node in self.nodes if node.validated)


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


def _engine_lines_for(receipt: Any, host: Any) -> tuple[str, ...]:
    """Recover the program's own words about this failure, if it left any.

    The parser already quotes them into the result-validation receipt's
    observations; they simply never travelled any further.  Everything here is
    read defensively, because a reader who has already lost the run must not
    also lose the reason to a missing field.
    """

    if host is None:
        return ()
    try:
        digest = _field(receipt, "result_validation_receipt_sha256")
        validation = host.result_validation_receipts.get(digest)
        observations = getattr(validation, "observations", None) or {}
    except (ContractError, AttributeError):
        return ()
    for observation in observations.values():
        if not isinstance(observation, Mapping):
            continue
        failure = observation.get("native_failure")
        if isinstance(failure, Mapping):
            lines = failure.get("engine_lines") or ()
            if isinstance(lines, (list, tuple)):
                return tuple(str(item) for item in lines if str(item).strip())
    return ()


def _execution_failure_summary(receipt: Any, host: Any = None) -> str:
    """Describe why a node did not succeed, from its own receipt.

    Empty for a node that succeeded.  Otherwise name the terminal state and
    whatever the receipt actually observed -- the wrapper and child exit
    statuses, any findings, and the program's own account of the failure.

    That last part matters more than the rest.  The host's own vocabulary of
    failure classes is closed, so an unanticipated failure classifies as
    ``native_runtime`` and its canonical template says only that an error
    occurred.  Four coupled-cluster nodes once died because the requested
    approximation is not size-consistent -- fatal for exactly the quantity
    being computed -- and the engine explained that plainly while every layer
    above it reported the empty string.  Quoting the engine is not
    interpreting it, and it is what lets the reason reach whoever must re-plan.
    """

    state = str(_field(receipt, "execution_state") or "")
    if state in {"validated", "engine_complete"}:
        return ""
    parts = [f"execution_state={state or 'unknown'}"]
    for label in ("wrapper_exit_status", "child_exit_status"):
        # These are optional on the receipt, and a summary must never be the
        # thing that fails: a reader who has already lost the run should not
        # also lose the reason to a missing field.
        try:
            value = _field(receipt, label)
        except ContractError:
            continue
        if value not in (None, ""):
            parts.append(f"{label}={value}")
    try:
        findings = tuple(_field(receipt, "findings") or ())
    except ContractError:
        findings = ()
    if findings:
        parts.append("findings=" + ", ".join(str(item) for item in findings))
    engine_lines = _engine_lines_for(receipt, host)
    if engine_lines:
        # Attributed, so nothing downstream mistakes the engine's words for a
        # host claim about readiness, validity, or what to do next.
        parts.append(
            "engine reported (verbatim): " + " | ".join(engine_lines)
        )
    return "; ".join(parts)

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
        execution_bundle: Any,
        approval_workspace: Path,
        claim_workspace_bundle: bool = True,
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
        self.execution_bundle = execution_bundle
        self.approval_workspace = approval_workspace
        self.claim_workspace_bundle = bool(claim_workspace_bundle)
        self._bundle_claimed = False
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

    def _verify_launch_and_claim_once(
        self, *, node_id: str, invocation_sha256: str
    ) -> None:
        """Compare the real launch first, then claim the one-shot bundle."""

        review = self.execution_bundle.node_review(node_id)
        self.host.verify_reviewed_real_execution_argv(
            node_id=node_id,
            invocation_sha256=invocation_sha256,
            review=review,
        )
        if self._bundle_claimed:
            return
        if self.claim_workspace_bundle:
            from chemsmart.agent.live_session import (
                claim_workflow_execution_approval_bundle,
            )

            claim_workflow_execution_approval_bundle(
                self.execution_bundle,
                workspace=self.approval_workspace,
            )
        self._bundle_claimed = True

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

            self._verify_launch_and_claim_once(
                node_id=node_id,
                invocation_sha256=invocation_sha256,
            )
            executed = self._call(
                "execute_approved_program_node", node_id=node_id
            )
            receipt = _field(executed, "execution")
            for item in _field(executed, "produced_handoffs") or ():
                handoff = (
                    item["handoff"] if isinstance(item, Mapping) else item
                )
                artifact = (
                    handoff.get("geometry_artifact_id", "")
                    if isinstance(handoff, Mapping)
                    else getattr(handoff, "geometry_artifact_id", "")
                )
                if artifact:
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
                # A node whose engine run failed reported nothing but the word
                # "failed": the reason sat in its own receipt and never reached
                # the operator, who had to open the raw program output to learn
                # that the engine had refused the route.  Say what the receipt
                # knows.
                failure=_execution_failure_summary(receipt, self.host),
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
                plan=self.plan,
                failure_rule_ids=("execution.prepare.blocked",),
                timestamp=stamp,
            )
        run_state = transition_workflow_node(
            run_state,
            node_id=outcome.node_id,
            new_state="running",
            plan=self.plan,
            invocation_sha256=outcome.invocation_sha256,
            timestamp=stamp,
        )
        if not outcome.validated:
            return transition_workflow_node(
                run_state,
                node_id=outcome.node_id,
                new_state="failed",
                plan=self.plan,
                execution_receipt_sha256=outcome.execution_receipt_sha256,
                failure_rule_ids=outcome.rule_ids
                or ("execution.state." + (outcome.state or "unknown"),),
                timestamp=stamp,
            )
        run_state = transition_workflow_node(
            run_state,
            node_id=outcome.node_id,
            new_state="engine_complete",
            plan=self.plan,
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
            plan=self.plan,
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
        # The approval, not the plan, names what may run.  A scientific plan
        # keeps a stage this release cannot execute rather than dropping it,
        # and such a stage carries no approved binding.
        approved_node_ids = {
            binding.node_id for binding in self.approval.node_bindings
        }
        while True:
            ready = tuple(
                node_id
                for node_id in derive_ready_node_ids(
                    self.plan, run_state, data_edge_bindings
                )
                if node_id not in seen and node_id in approved_node_ids
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
        done = {item.node_id for item in executed if item.validated}
        status = "completed" if done == approved_node_ids else "partial"
        return WorkflowExecutionResultV1(
            workflow_id=self.plan.workflow_id,
            plan_sha256=self.plan.plan_sha256,
            approval_sha256=self.frozen_approval.approval_sha256,
            run_directory=str(self.run_directory),
            nodes=tuple(executed),
            status=status,
            provider_calls=0,
            non_executable_node_ids=(
                self.execution_bundle.non_executable_node_ids
            ),
        )


def _execution_inputs_from_bundle(
    *,
    bundle: WorkflowExecutionApprovalBundleV1,
    workspace: Path,
    run_directory: Path,
) -> dict[str, Any]:
    """Compose the normal ChemSmart execution host from one typed review.

    This is deliberately object based.  The TUI has already displayed the
    project settings, molecular state, CLI operations, DAG and resource
    bounds held by ``bundle``; no approval file or user-entered digest is a
    second authority over those ChemSmart objects.
    """

    from chemsmart.agent.live_session import (
        _approved_project_artifacts,
        _parse_bounded_execution_envelope_record,
        _write_execution_server_profile,
    )

    approval = bundle.workflow_approval
    resources = bundle.execution_resources
    if Path(approval.workspace).resolve() != workspace:
        raise ContractError("prepared workflow targets another workspace")
    if approval.task_spec_sha256 != (
        bundle.approved_scientific_plan.task_spec_sha256
    ):
        raise ContractError("prepared workflow targets another task")
    envelope = _parse_bounded_execution_envelope_record(
        bundle.execution_envelope,
        resources=resources,
    )
    requested_scratch_root = Path(envelope.scratch_root)
    if requested_scratch_root.is_symlink():
        raise ContractError("execution scratch root cannot be a symlink")
    scratch_root = requested_scratch_root.resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    server_profile = _write_execution_server_profile(
        run_directory,
        resources,
        scratch_root=scratch_root,
    )
    path_value = os.environ.get("PATH", "")
    xtb_executable = os.environ.get("CHEMSMART_XTB_EXECUTABLE") or shutil.which(
        "xtb"
    )
    executable_directory = (
        str(Path(xtb_executable).expanduser().parent)
        if xtb_executable
        else ""
    )
    environment = {
        "PATH": (
            path_value
            if not executable_directory
            else (
                executable_directory
                if not path_value
                else executable_directory + os.pathsep + path_value
            )
        ),
        "PYTHONNOUSERSITE": "1",
    }
    return {
        "approved_workspace": workspace,
        "execution_resources": resources,
        "workflow_execution_approval": approval,
        "frozen_workflow_approval": bundle.frozen_workflow_approval,
        "execution_server": str(server_profile),
        "execution_environment": environment,
        "approved_project_artifacts": _approved_project_artifacts(
            workspace, approval
        ),
        "approved_scientific_plan": bundle.approved_scientific_plan,
        "approved_materialized_workflow": (
            bundle.approved_materialized_workflow
        ),
        "approved_environment_identities": (
            bundle.approved_environment_identities
        ),
        "stationary_point_policy": bundle.stationary_point_policy,
    }


def _execute_workflow_bundle(
    *,
    bundle: WorkflowExecutionApprovalBundleV1,
    workspace: Path,
    run_directory: Path,
    claim_workspace_bundle: bool,
) -> WorkflowExecutionResultV1:
    """Execute the typed ChemSmart DAG represented by an approved review."""

    from chemsmart.agent.live_session import (
        _bootstrap_conformance,
        _observe_environments,
    )

    workspace = Path(workspace).resolve()
    run_directory = Path(run_directory).resolve()
    run_directory.mkdir(parents=True, exist_ok=True)
    effective_task_spec_sha256 = (
        bundle.approved_scientific_plan.task_spec_sha256
    )
    inputs = _execution_inputs_from_bundle(
        bundle=bundle,
        workspace=workspace,
        run_directory=run_directory,
    )
    plan = inputs.pop("approved_scientific_plan")
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
        task_spec_sha256s=(effective_task_spec_sha256,),
        # ``_execution_composition_inputs`` already writes and supplies the
        # execution server profile; naming it again here would be a second
        # source of truth for where real jobs are launched.
        scientific_workflow_plan=plan,
        materialized_workflow=materialized,
        approved_environment_identities=environment_identities,
        execution_server_file_sha256=file_sha256(inputs["execution_server"]),
        execution_environment_remove=_provider_secret_environment_labels(),
        **inputs,
    )
    return ApprovedWorkflowExecutor(
        host=host,
        plan=plan,
        approval=approval,
        frozen_approval=frozen_approval,
        initial_artifacts=initial_artifacts,
        project_artifacts=project_artifacts,
        task_spec_sha256=effective_task_spec_sha256,
        run_directory=run_directory,
        execution_bundle=bundle,
        approval_workspace=workspace,
        claim_workspace_bundle=claim_workspace_bundle,
    ).run()


def execute_prepared_workflow(
    *,
    review: WorkflowExecutionReviewV1,
    actor: str,
    execution_id: str,
    workspace: Path,
    run_directory: Path,
) -> WorkflowExecutionResultV1:
    """Approve the displayed in-memory ChemSmart workflow and run it once.

    The explicit TUI action is the authority.  Internal content digests remain
    provenance and mutation evidence, but the user does not create, reload or
    retype an approval artifact.  Atomic node-launch reservation prevents an
    accidental duplicate within this run.
    """

    bundle = approve_workflow_execution_review(
        review,
        approval_id=str(execution_id),
        approved_review_sha256=review.review_sha256,
        actor=str(actor).strip(),
        resolution_id=str(execution_id) + "-approval",
    )
    return _execute_workflow_bundle(
        bundle=bundle,
        workspace=workspace,
        run_directory=run_directory,
        claim_workspace_bundle=False,
    )


def execute_approved_workflow(
    *,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    task_spec_sha256: str = "",
    expected_approval_file_sha256: str = "",
) -> WorkflowExecutionResultV1:
    """Compatibility adapter for a previously persisted v1 approval file."""

    from chemsmart.agent.live_session import (
        load_workflow_execution_approval_bundle,
    )

    bundle = load_workflow_execution_approval_bundle(
        Path(approval_file).resolve(),
        expected_file_sha256=expected_approval_file_sha256,
    )
    effective_task = (
        str(task_spec_sha256).strip()
        or bundle.approved_scientific_plan.task_spec_sha256
    )
    if effective_task != bundle.approved_scientific_plan.task_spec_sha256:
        raise ContractError("task specification differs from approval bundle")
    return _execute_workflow_bundle(
        bundle=bundle,
        workspace=workspace,
        run_directory=run_directory,
        claim_workspace_bundle=True,
    )


__all__ = [
    "ApprovedWorkflowExecutor",
    "ExecutedNodeV1",
    "PROGRAM_NODE_SEQUENCE",
    "WorkflowExecutionResultV1",
    "execute_approved_workflow",
    "execute_prepared_workflow",
]
