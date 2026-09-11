"""Execute one exact human-approved scientific workflow without a provider.

The executor follows the frozen DAG and normal ChemSmart host tools.  It does
not choose science, re-plan, or widen approval.
"""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    canonical_json,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.api_access import (
    DEFAULT_KEY_LABELS,
    PROVIDER_KEY_LABEL_TOKENS,
    normalize_key_label,
)
from chemsmart.agent.capabilities import load_program_capabilities
from chemsmart.agent.cli_schema import build_live_click_schema
from chemsmart.agent.execution import (
    WorkflowExecutionApprovalBundleV1,
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
class ExecutedAnalysisNodeV1:
    """One approved analysis node's fate in the provider-free walk."""

    node_id: str
    analysis_kind: str
    #: ``executed`` | ``failed`` | ``skipped`` | ``blocked_unsupported``
    state: str
    receipt_sha256s: tuple[str, ...] = ()
    reason: str = ""
    #: Outputs this node was asked for that its result does not carry.  A
    #: node that ran and delivered some of what was requested still counts
    #: as executed: the state records what the walk did with the node, and
    #: the receipt's own ``partial`` status records what the evidence
    #: contains.  Consumers naming one of these skip; consumers naming only
    #: the delivered siblings do not.
    absent_output_ids: tuple[str, ...] = ()


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
    #: The approved analysis chain's fate, when the bundle carried one.
    analysis_nodes: tuple[ExecutedAnalysisNodeV1, ...] = ()
    #: "" (no chain) | "completed" | "partial" | "not_run"
    analysis_status: str = ""
    analysis_completion_receipt_sha256s: tuple[str, ...] = ()
    analysis_report_path: str = ""
    #: Per-record delivery verdicts, present only when the plan carried
    #: more than one record (disconnected sub-DAG).  Derived from node
    #: states by a fixed rule, never stored beside them; reached states
    #: and verdicts stay separate fields, and never-attempted is not
    #: failed.  A batch of N is N observations -- there is deliberately
    #: no aggregate quantity here.
    record_delivery: tuple[dict[str, Any], ...] = ()

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
    # The goal's original geometry behind a reached or displaced input,
    # located the same way, so the sensors can walk from the root.  Bytes
    # the workspace no longer holds leave the node without a root
    # comparison rather than without a run.
    for binding in approval.node_bindings:
        root_sha256 = str(getattr(binding, "root_artifact_sha256", "") or "")
        root_id = str(getattr(binding, "root_artifact_id", "") or "")
        if not root_sha256 or not root_id:
            continue
        existing = artifacts.get(root_id)
        if existing is not None:
            if existing.sha256 != root_sha256:
                raise ContractError(
                    "one approved root artifact ID names different bytes"
                )
            continue
        try:
            root_path = _locate_by_digest(workspace, root_sha256)
        except ContractError:
            continue
        artifacts[root_id] = TrustedArtifactRefV1(
            artifact_id=root_id,
            kind="geometry_xyz",
            sha256=root_sha256,
            size_bytes=root_path.stat().st_size,
            path=str(root_path),
            cli_value=str(root_path),
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

    # Programs record the summary at different depths: ORCA puts it on the
    # per-program observation, Gaussian on each per-artifact row beneath it.
    # Search one level into nested collections rather than requiring every
    # branch to hoist it, so a program whose shape differs loses its
    # diagnostic silently no longer.
    # observations -> program -> "outputs" tuple -> row is already four levels,
    # so the bound is set above the deepest shape any program uses today
    # rather than at the shallowest one that happens to work.
    def _summaries(value: Any, depth: int = 0):
        if depth > 5:
            return
        if isinstance(value, Mapping):
            failure = value.get("native_failure")
            if isinstance(failure, Mapping):
                yield failure
            for item in value.values():
                yield from _summaries(item, depth + 1)
        elif isinstance(value, (list, tuple)):
            for item in value:
                yield from _summaries(item, depth + 1)

    for failure in _summaries(observations):
        lines = failure.get("engine_lines") or ()
        if isinstance(lines, (list, tuple)):
            quoted = tuple(str(item) for item in lines if str(item).strip())
            if quoted:
                return quoted
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
        parts.append("engine reported (verbatim): " + " | ".join(engine_lines))
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
        should_stop: Any = None,
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
        self.should_stop = should_stop
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
        self._refuse_launch_the_program_already_refused(node_id)
        if self._bundle_claimed:
            return
        if self.claim_workspace_bundle:
            from chemsmart.agent.live_session import (
                claim_workflow_execution_approval_bundle,
                continue_workflow_execution_approval_bundle,
            )
            from chemsmart.agent.runtime.event_store import (
                ExecutionBundleAlreadyConsumedError,
            )

            try:
                claim_workflow_execution_approval_bundle(
                    self.execution_bundle,
                    workspace=self.approval_workspace,
                )
            except ExecutionBundleAlreadyConsumedError:
                # A spent approval admits exactly one more shape: the same
                # recorded run continuing in its original run directory,
                # still incomplete.  Anything else re-raises the refusal
                # inside the continuation check itself.
                continue_workflow_execution_approval_bundle(
                    self.execution_bundle,
                    workspace=self.approval_workspace,
                    run_event_store=self.host.event_store,
                )
        self._bundle_claimed = True

    def _refuse_launch_the_program_already_refused(self, node_id: str) -> None:
        """Do not spend an engine call on bytes the program rejected.

        The input-check probe is an observation and never a verdict,
        because the charter reserves the decision for the human reading
        the review. Under a goal's standing approval there is no human at
        this moment, so that category of host fact had no consumer at
        all: ino3-r17 launched six nodes over six aborted checks and all
        six died; po3-r18 launched two and both died with the error the
        probe had printed, verbatim, 0.107 s and no charge earlier.
        Eight engine calls across two windows.

        This costs no scientific freedom. An aborted input check is not
        a judgement about chemistry -- it is the program refusing its own
        input, in its own words, on these exact bytes. What it costs is
        the ability to spend the budget the human granted on a run that
        cannot start.

        The override is to re-probe: repair the field the program named,
        recompile, and the new input gets its own check. That keeps the
        model sovereign over the repair -- ORCA offers three legal
        keywords and choosing among them is a method decision -- and it
        is environment-safe, because a probe aborted on an executable
        that has since changed is superseded by the new one rather than
        standing forever.
        """

        observations = tuple(
            getattr(self.execution_bundle, "node_observations", {}).get(
                node_id, ()
            )
            or ()
        )
        aborted = [
            str(line)
            for line in observations
            if str(line).startswith("input-check probe: aborted")
        ]
        if not aborted:
            return
        engine_lines = [
            str(line).strip()
            for line in observations
            if str(line).lstrip().startswith("orca:")
        ]
        detail = "; ".join(engine_lines[:4]) or aborted[0]
        raise ContractError(
            f"node {node_id!r} is not launched: the program's own input "
            f"check refused these exact bytes -- {detail}. This is the "
            "program's word on its own input, not a judgement about the "
            "chemistry. Repair the field it names and compile the node "
            "again; the new input is probed on its own, and a passing "
            "check admits the launch."
        )

    def _input_artifact_id(self, binding: Any) -> str:
        if binding.input_mode == "initial":
            artifact = self.initial_artifacts.get(binding.initial_artifact_id)
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
                **(
                    {"excursion": binding.excursion}
                    if binding.excursion
                    else {}
                ),
                **(
                    {"internal_coordinates": binding.internal_coordinates}
                    if binding.internal_coordinates
                    else {}
                ),
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
                geometry_artifact_sha256=(input_artifact.sha256),
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
            # A refusal that lives only in this process is
            # indistinguishable, afterwards, from a node that was never
            # approved: a revision that reused a failed node's id met the
            # node-workspace guard, the outcome stayed in memory, and the
            # goal settled naming nothing (live, 2026-09-02). The word goes
            # into the run stream before anything else happens.
            from chemsmart.agent.runtime.events import EventKind

            self.host.event_store.append(
                turn_id=f"exec-refused-{node_id}",
                kind=EventKind.WORKFLOW_NODE_LAUNCH_REFUSED.value,
                payload={
                    "node_id": node_id,
                    "program": binding.program,
                    "jobtype": binding.jobtype,
                    "reason": str(error),
                },
            )
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

    def _run_analysis_phase(
        self, toolchain: Any
    ) -> tuple[
        tuple["ExecutedAnalysisNodeV1", ...], str, tuple[str, ...], str
    ]:
        """Walk the approved analysis chain with host-computed arguments.

        No provider exists in this process.  Every argument is derived from
        the digest-bound plan the human approved; the dispatched tools emit
        their own typed receipts and events, exactly as they do in a session.
        A failed validation VERDICT is a completed determination; only a
        refused kernel call fails a node, and its dependents are skipped
        with the failure named.
        """

        # The approved plan opens the guides its own operations need, so
        # the chain never meets the stem's leaf gate; the activation is
        # recorded like any other, with its signal.
        from chemsmart.agent.guides import guides_from_plan
        from chemsmart.agent.runtime.events import EventKind
        from chemsmart.agent.scientific_toolchain import (
            RegisteredResultInputIntentV1,
        )
        from chemsmart.analysis.quantity_expressions import (
            QuantityExpressionError,
        )
        from chemsmart.analysis.result_quantities import (
            QuantityContractError,
            canonical_thermochemistry_quantity,
        )

        expression_items = [
            item
            for node in toolchain.analysis_nodes
            for item in getattr(node, "expression_nodes", ())
        ]
        self.host.activate_guides(
            f"exec-analysis-{toolchain.plan_sha256[:8]}",
            guides_from_plan(
                operations=[
                    str(item.get("operation", "")) for item in expression_items
                ],
                constants=[
                    str(item.get("constant_name", ""))
                    for item in expression_items
                    if str(item.get("operation", "")) == "constant"
                ],
            ),
            signal="approved_plan",
        )
        program_by_kind = {
            "pyscf_hdf5": "pyscf",
            "orca_output": "orca",
            "gaussian_output": "gaussian",
            "xtb_output": "xtb",
            "geometry_xyz": "xyz",
        }
        analysis_by_id = {
            node.node_id: node for node in toolchain.analysis_nodes
        }
        calculation_ids = set(toolchain.calculation_node_ids)
        # The producer's typed result kind, from the approved V2 plan: an
        # engine run registers several artifact kinds (an xTB node leaves a
        # log AND a geometry), and only the program's reader kind is the
        # result an extraction reads. Observed live: an xtb single point
        # offered both xtb_output and geometry_xyz and the walker refused
        # the ambiguity instead of resolving it.
        program_by_node = {
            node.node_id: node.program
            for node in getattr(self.plan, "nodes", ())
        }

        def _producer_result_kind(producer: str) -> str:
            program = program_by_node.get(producer, "")
            from chemsmart.agent.postprocessing import (
                typed_result_artifact_kind,
            )

            return typed_result_artifact_kind(program)

        settled: dict[str, ExecutedAnalysisNodeV1] = {}
        #: (producer node, output id) -> (receipt digest, receipt quantity id)
        outputs: dict[tuple[str, str], tuple[str, str]] = {}
        ledger: list[str] = []
        # Settlement rounds are first-class: a continuation re-enters the
        # same run directory and re-walks the chain, and a node that
        # settled `failed` before the suspension may settle `executed`
        # after its producer validated on resume.  Each round's
        # settlements and report therefore carry a round tag derived from
        # the stream length at phase start -- unique across invocations,
        # shared within one -- so the append-only stream records every
        # round and the latest one supersedes by order, never by edit.
        # The first live batch resume conflicted here on the un-tagged
        # keys.
        round_tag = f"r{len(self.host.event_store.read_events())}"

        def _emit(record: ExecutedAnalysisNodeV1) -> None:
            settled[record.node_id] = record
            self.host.event_store.append(
                turn_id=f"exec-analysis-{len(settled):03d}",
                kind=EventKind.WORKFLOW_ANALYSIS_NODE_SETTLED.value,
                payload={
                    "node_id": record.node_id,
                    "analysis_kind": record.analysis_kind,
                    "state": record.state,
                    "reason": record.reason,
                    "receipt_sha256s": record.receipt_sha256s,
                    "toolchain_plan_sha256": toolchain.plan_sha256,
                },
                idempotency_key=(
                    "analysis-node:"
                    + toolchain.plan_sha256
                    + ":"
                    + record.node_id
                    + ":"
                    + round_tag
                ),
            )

        requested_registered_ids = {
            item.artifact_id
            for node in toolchain.analysis_nodes
            for item in node.inputs
            if isinstance(item, RegisteredResultInputIntentV1)
        }
        if any(
            artifact_id not in self.host.artifacts
            for artifact_id in requested_registered_ids
        ):
            # A registered-result id is content-derived, so the approved
            # workspace can rebuild exactly the artifact the planning
            # session registered. Observed live: the first composed chain
            # to reuse a registered acid result crashed here with a bare
            # KeyError, because the fresh executor host held only the
            # artifacts this run produced.
            from chemsmart.agent.live_session import (
                discover_registered_result_artifacts,
            )

            for artifact in discover_registered_result_artifacts(
                Path(self.approval_workspace)
            ):
                if artifact.artifact_id in requested_registered_ids:
                    self.host.artifacts.setdefault(
                        artifact.artifact_id, artifact
                    )

        def _registered_artifact(artifact_id: str) -> Any:
            artifact = self.host.artifacts.get(artifact_id)
            if artifact is None:
                raise ContractError(
                    f"registered result {artifact_id!r} is not present in "
                    "the approved workspace"
                )
            return artifact

        def _producer_artifact(producer: str) -> Any:
            prefix = f"result.{producer}."
            result_kind = _producer_result_kind(producer)
            candidates = [
                artifact
                for artifact_id, artifact in self.host.artifacts.items()
                if artifact_id.startswith(prefix)
                and artifact.kind == result_kind
            ]
            if len(candidates) != 1:
                raise ContractError(
                    f"expected exactly one registered {result_kind!r} result "
                    f"artifact for producer {producer!r}; found "
                    f"{len(candidates)}"
                )
            return candidates[0]

        def _extraction_quantity_id(node: Any, output: Any) -> str:
            selector_ids = {
                selector.quantity_id for selector in node.selectors
            }
            if output.output_id in selector_ids:
                return output.output_id
            if len(selector_ids) == 1:
                return next(iter(selector_ids))
            raise ContractError(
                f"extraction output {output.output_id!r} names no selector "
                "quantity and the selector set is not a singleton"
            )

        def _resolve_source(item: Any) -> tuple[str, str]:
            key = (item.producer_node_id, item.producer_output_id)
            if key not in outputs:
                raise ContractError(
                    "analysis input references an output the walk has not "
                    f"produced: {key!r}"
                )
            return outputs[key]

        def _run_node(node: Any) -> ExecutedAnalysisNodeV1:
            kind = node.analysis_kind
            if kind == "result_extraction":
                sources: list[Any] = []
                for item in node.inputs:
                    if isinstance(item, RegisteredResultInputIntentV1):
                        sources.append(_registered_artifact(item.artifact_id))
                    elif item.producer_node_id in calculation_ids:
                        sources.append(
                            _producer_artifact(item.producer_node_id)
                        )
                if len(sources) != 1:
                    raise ContractError(
                        "result extraction requires exactly one result "
                        f"artifact; resolved {len(sources)}"
                    )
                artifact = sources[0]
                receipt = self._call(
                    "extract_result_quantities",
                    program=program_by_kind[artifact.kind],
                    artifact_id=artifact.artifact_id,
                    selectors=[
                        {
                            "quantity_id": selector.quantity_id,
                            "selector": selector.selector,
                        }
                        for selector in node.selectors
                    ],
                )
                digest = _field(receipt, "receipt_sha256")
                absent_quantity_ids = {
                    str(item[0]) for item in (_field(receipt, "absent") or ())
                }
                missing: list[str] = []
                for output in node.outputs:
                    quantity_id = _extraction_quantity_id(node, output)
                    if quantity_id in absent_quantity_ids:
                        # Publishing nothing under this key is what makes a
                        # consumer that names it skip.  A consumer naming a
                        # sibling finds its value exactly where it always
                        # was.
                        missing.append(output.output_id)
                        continue
                    outputs[(node.node_id, output.output_id)] = (
                        digest,
                        quantity_id,
                    )
                return ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=kind,
                    state="executed",
                    receipt_sha256s=(digest,),
                    absent_output_ids=tuple(sorted(missing)),
                )
            if kind == "thermochemistry":
                if node.temperature_k is None or node.pressure_atm is None:
                    raise ContractError(
                        "a thermochemistry stage requires explicit "
                        "temperature and pressure"
                    )
                sources = []
                for item in node.inputs:
                    if isinstance(item, RegisteredResultInputIntentV1):
                        sources.append(_registered_artifact(item.artifact_id))
                    elif item.producer_node_id in calculation_ids:
                        sources.append(
                            _producer_artifact(item.producer_node_id)
                        )
                if len(sources) != 1:
                    raise ContractError(
                        "thermochemistry requires exactly one result "
                        f"artifact; resolved {len(sources)}"
                    )
                artifact = sources[0]
                arguments: dict[str, Any] = {
                    "program": program_by_kind[artifact.kind],
                    "artifact_id": artifact.artifact_id,
                    "temperature_k": float(node.temperature_k),
                    "pressure_atm": float(node.pressure_atm),
                    "entropy_method": node.entropy_method,
                    "alpha": node.alpha,
                    "use_weighted_mass": node.use_weighted_mass,
                    "frequency_scale_factor": node.frequency_scale_factor,
                }
                if node.concentration_mol_l is not None:
                    arguments["concentration_mol_l"] = node.concentration_mol_l
                if node.entropy_cutoff_cm1 is not None:
                    arguments["entropy_cutoff_cm1"] = node.entropy_cutoff_cm1
                if node.enthalpy_cutoff_cm1 is not None:
                    arguments["enthalpy_cutoff_cm1"] = node.enthalpy_cutoff_cm1
                receipt = self._call("derive_thermochemistry", **arguments)
                digest = _field(receipt, "receipt_sha256")
                for output in node.outputs:
                    outputs[(node.node_id, output.output_id)] = (
                        digest,
                        canonical_thermochemistry_quantity(
                            output.quantity_kind
                        ),
                    )
                return ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=kind,
                    state="executed",
                    receipt_sha256s=(digest,),
                )
            if kind == "quantity_expression":
                expression_inputs = []
                for item in node.inputs:
                    if isinstance(item, RegisteredResultInputIntentV1):
                        raise ContractError(
                            "a quantity expression consumes typed receipts, "
                            "not raw registered results"
                        )
                    digest, quantity_id = _resolve_source(item)
                    expression_inputs.append(
                        {
                            "input_id": item.input_id,
                            "receipt_sha256": digest,
                            "quantity_id": quantity_id,
                        }
                    )
                receipt = self._call(
                    "evaluate_quantity_expression",
                    expression_id=node.node_id,
                    inputs=expression_inputs,
                    nodes=json.loads(
                        canonical_json(tuple(node.expression_nodes))
                    ),
                    output_node_ids=list(node.expression_output_node_ids),
                )
                digest = _field(receipt, "receipt_sha256")
                for output in node.outputs:
                    outputs[(node.node_id, output.output_id)] = (
                        digest,
                        output.output_id,
                    )
                return ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=kind,
                    state="executed",
                    receipt_sha256s=(digest,),
                )
            if kind == "scientific_validation":
                validation_inputs = []
                for item in node.inputs:
                    digest, quantity_id = _resolve_source(item)
                    validation_inputs.append(
                        {
                            "input_id": item.input_id,
                            "receipt_sha256": digest,
                            "quantity_id": quantity_id,
                        }
                    )
                receipt = self._call(
                    "evaluate_scientific_validation",
                    workflow_id=toolchain.workflow_id,
                    node_id=node.node_id,
                    inputs=validation_inputs,
                )
                digest = _field(receipt, "receipt_sha256")
                raw_outputs = (
                    receipt.get("outputs", ())
                    if isinstance(receipt, Mapping)
                    else getattr(receipt, "outputs", ())
                )
                verdicts = tuple(
                    "{}={}".format(
                        _field(item, "quantity_id"), _field(item, "value")
                    )
                    for item in raw_outputs
                )
                for output in node.outputs:
                    outputs[(node.node_id, output.output_id)] = (
                        digest,
                        output.output_id,
                    )
                return ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=kind,
                    state="executed",
                    receipt_sha256s=(digest,),
                    reason="verdicts: " + "; ".join(verdicts),
                )
            if kind == "claim_rendering":
                producer_units = {
                    (producer.node_id, output.output_id): output.unit
                    for producer in toolchain.analysis_nodes
                    for output in producer.outputs
                }
                claims = []
                for item in node.inputs:
                    digest, quantity_id = _resolve_source(item)
                    unit = producer_units.get(
                        (item.producer_node_id, item.producer_output_id)
                    )
                    if unit is None:
                        raise ContractError(
                            "claim input names an output with no declared "
                            f"unit: {item.producer_output_id!r}"
                        )
                    claim = {
                        "claim_id": item.input_id,
                        "receipt_sha256": digest,
                        "quantity_id": quantity_id,
                        "display_unit": unit,
                    }
                    # The estimator the plan named, evaluated in this
                    # same provider-free walk. An estimate cannot
                    # predate its results and its estimator can, so a
                    # precision-bearing delivery no longer arrives
                    # `unstated` by construction and no longer costs a
                    # further cycle to say what it already computed.
                    # The number is the host's own, cited to the output
                    # that produced it, which is what `measured` means.
                    estimator_node = getattr(
                        item, "uncertainty_producer_node_id", ""
                    )
                    estimator_output = getattr(
                        item, "uncertainty_producer_output_id", ""
                    )
                    if estimator_node and estimator_output:
                        key = (estimator_node, estimator_output)
                        if key not in outputs:
                            raise ContractError(
                                "a planned uncertainty names an output the "
                                f"walk has not produced: {key!r}"
                            )
                        spread_digest, spread_quantity = outputs[key]
                        # The executor names the number and never
                        # writes it, exactly as it does for the claimed
                        # value: the host reads the cited quantity and
                        # copies it in.
                        claim["uncertainty_basis"] = "measured"
                        claim["uncertainty_reference"] = (
                            f"{spread_digest}:{spread_quantity}"
                        )
                    claims.append(claim)
                receipt = self._call(
                    "record_analysis_claims",
                    task_spec_sha256=self.task_spec_sha256,
                    claims=claims,
                )
                digest = _field(receipt, "receipt_sha256")
                return ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=kind,
                    state="executed",
                    receipt_sha256s=(digest,),
                )
            raise ContractError(
                f"the executor cannot run analysis kind {kind!r}"
            )

        for node_id in toolchain.node_order:
            node = analysis_by_id.get(node_id)
            if node is None:
                continue
            if node.support_state == "blocked_unsupported":
                _emit(
                    ExecutedAnalysisNodeV1(
                        node_id=node.node_id,
                        analysis_kind=node.analysis_kind,
                        state="blocked_unsupported",
                        reason=node.blocked_reason
                        or "declared non-executable intent",
                    )
                )
                continue
            broken = tuple(
                dependency
                for dependency in node.dependencies
                if dependency in settled
                and settled[dependency].state
                in {"failed", "skipped", "blocked_unsupported"}
            )
            # A dependency that ran and delivered some of what it was asked
            # for breaks only the consumers that named something it did not
            # deliver.  This used to key on the node alone, so one refused
            # selector settled every consumer of every sibling quantity --
            # in the run that prompted this, twelve nodes that named no
            # absent value at all.  The output-level edge was already here,
            # validated three lines before the plan collapsed it into a node
            # id, and used at both plan time and consume time.
            starved = tuple(
                sorted(
                    {
                        f"{item.producer_node_id}.{item.producer_output_id}"
                        for item in node.inputs
                        if getattr(item, "producer_node_id", None) in settled
                        and item.producer_output_id
                        in settled[item.producer_node_id].absent_output_ids
                    }
                )
            )
            if broken or starved:
                if broken:
                    reason = "upstream analysis did not execute: " + ", ".join(
                        broken
                    )
                else:
                    reason = "the result carries no value for: " + ", ".join(
                        starved
                    )
                _emit(
                    ExecutedAnalysisNodeV1(
                        node_id=node.node_id,
                        analysis_kind=node.analysis_kind,
                        state="skipped",
                        reason=reason,
                    )
                )
                continue
            try:
                record = _run_node(node)
            except (
                ContractError,
                QuantityContractError,
                QuantityExpressionError,
            ) as exc:
                # The analysis package raises its own typed errors -- an
                # extraction asking a closed-shell result for <S^2>, an
                # expression called with the wrong arity.  Six benchmark
                # runs died on the first and two on the second because
                # only ContractError was settled here; the escape killed
                # the executor and erased every refusal message with it.
                record = ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=node.analysis_kind,
                    state="failed",
                    reason=str(exc),
                )
            except Exception as exc:  # noqa: BLE001
                # This overturns a deliberate decision, so it states the
                # argument. The rule was "typed errors settle the node,
                # a bare exception stays what it should be -- a defect
                # that crashes" (tool_runtime._derive_thermochemistry).
                # The intent was to keep host defects loud. Measured
                # against a live run, crashing does not make a defect
                # loud: a kernel reached into a result's absent
                # thermochemistry, `float + None` raised TypeError, and
                # the escape killed the goal. Three hours twenty minutes
                # of engine time and five converged spin states were
                # lost; the ledger ended three entries long with no
                # run_recorded and no settlement, so a resume would have
                # believed the budget unspent. The traceback went to a
                # driver.log nothing reads (NOVEL-2 ino1, 2026-09-04).
                #
                # Loud and settled are not in conflict, and this is the
                # louder of the two: the failure is recorded with its
                # exception type in the durable stream the driver and
                # the human both read, marked a host defect rather than
                # a scientific finding, while every sibling node's
                # receipt survives and the goal settles.
                record = ExecutedAnalysisNodeV1(
                    node_id=node.node_id,
                    analysis_kind=node.analysis_kind,
                    state="failed",
                    reason=(
                        "host defect, not a scientific finding: "
                        f"{type(exc).__name__}: {exc}"
                    ),
                )
            _emit(record)
            ledger.extend(record.receipt_sha256s)

        executed_all = all(
            record.state in {"executed", "blocked_unsupported"}
            for record in settled.values()
        )
        # A node that ran and was refused one quantity still counts as
        # executed, so a chain can now finish with an absence in it.  The
        # run reads completed when it delivered what the session declared
        # the chain was for, and partial when an absence took one of those
        # required outputs -- which can happen with nothing skipped at all,
        # if the missing output was a deliverable no later node consumed.
        starved_requirements = tuple(
            sorted(
                {
                    output_id
                    for record in settled.values()
                    for output_id in record.absent_output_ids
                    if output_id in set(toolchain.required_output_ids)
                }
            )
        )
        executed_all = executed_all and not starved_requirements
        analysis_status = "completed" if executed_all else "partial"
        completion_receipts: tuple[str, ...] = ()
        report_path = ""

        def _write_report(report: str, filename: str) -> str:
            report_directory = self.run_directory / "analysis"
            report_directory.mkdir(parents=True, exist_ok=True)
            target = report_directory / filename
            target.write_text(report + "\n", encoding="utf-8")
            return str(target)

        partial_findings: tuple[str, ...] = ()
        if executed_all and ledger:
            try:
                completion_receipts = (
                    self.host.evaluate_approved_toolchain_completion(
                        toolchain,
                        source_receipt_sha256s=tuple(ledger),
                    )
                )
                report = self.host.render_completed_analysis_report(
                    completion_receipts[0]
                )
            except ContractError as exc:
                # Every node executed and then the completion binding or the
                # report renderer refused.  A crash here discards the whole
                # chain's receipts and says nothing; the refusal is instead
                # recorded as durable evidence, the status drops to partial,
                # and the partial envelope below still delivers what the
                # chain validated.
                completion_receipts = ()
                analysis_status = "partial"
                partial_findings = (
                    "completion or report rendering refused: " + str(exc),
                )
                self.host.event_store.append(
                    turn_id="exec-analysis-report",
                    kind=(
                        EventKind.WORKFLOW_ANALYSIS_COMPLETION_REFUSED.value
                    ),
                    payload={
                        "reason": str(exc),
                        "toolchain_plan_sha256": toolchain.plan_sha256,
                    },
                    idempotency_key=(
                        "analysis-completion-refused:"
                        + toolchain.plan_sha256
                        + ":"
                        + round_tag
                    ),
                )
            else:
                report_path = _write_report(
                    report, "completed-analysis-report.md"
                )
                self.host.event_store.append(
                    turn_id="exec-analysis-report",
                    kind=EventKind.WORKFLOW_ANALYSIS_REPORT_RENDERED.value,
                    payload={
                        "completion_receipt_sha256": completion_receipts[0],
                        "report_sha256": canonical_sha256(report),
                        "toolchain_plan_sha256": toolchain.plan_sha256,
                    },
                    idempotency_key=(
                        "analysis-report:"
                        + toolchain.plan_sha256
                        + ":"
                        + round_tag
                    ),
                )
        if analysis_status == "partial" and not report_path:
            # The partial-failure envelope: bind the receipts that DID land
            # plus findings naming every node that did not execute, and
            # render them as evidence at their rung.  A reader gets what was
            # validated with its limitation stated instead of nothing.  What
            # this enforces is disclosure, not exclusion: a number read off
            # a node that did not validate is in the claim record either
            # way, and the comment used to say it was kept out, which the
            # code has never done (2026-09-04).
            # A calculation node that never validated is disclosed in the
            # fixed shape the host's partial-completion gate demands, so a
            # batch record's failed engine run appears in the delivered
            # report by name rather than silently narrowing the chain.
            calculation_findings = tuple(
                f"{node_id} (calculation): did not validate"
                for node_id in toolchain.calculation_node_ids
                if not bool(
                    getattr(
                        self.host.execution_receipts.get(node_id),
                        "validated",
                        False,
                    )
                )
            )
            findings = calculation_findings + (
                partial_findings
                or tuple(
                    f"{record.node_id} ({record.analysis_kind}): "
                    f"{record.state}"
                    + (f" -- {record.reason}" if record.reason else "")
                    for record in settled.values()
                    if record.state in {"failed", "skipped"}
                )
            )
            try:
                completion_receipts = (
                    self.host.evaluate_partial_toolchain_completion(
                        toolchain,
                        source_receipt_sha256s=tuple(ledger),
                        findings=findings,
                    )
                )
                report = self.host.render_completed_analysis_report(
                    completion_receipts[0]
                )
            except ContractError as exc:
                completion_receipts = ()
                self.host.event_store.append(
                    turn_id="exec-analysis-report",
                    kind=(
                        EventKind.WORKFLOW_ANALYSIS_COMPLETION_REFUSED.value
                    ),
                    payload={
                        "reason": str(exc),
                        "toolchain_plan_sha256": toolchain.plan_sha256,
                    },
                    idempotency_key=(
                        "partial-analysis-refused:"
                        + toolchain.plan_sha256
                        + ":"
                        + round_tag
                    ),
                )
            else:
                report_path = _write_report(
                    report, "partial-analysis-report.md"
                )
                self.host.event_store.append(
                    turn_id="exec-analysis-report",
                    kind=EventKind.WORKFLOW_ANALYSIS_REPORT_RENDERED.value,
                    payload={
                        "completion_receipt_sha256": completion_receipts[0],
                        "report_sha256": canonical_sha256(report),
                        "toolchain_plan_sha256": toolchain.plan_sha256,
                        "report_kind": "partial",
                    },
                    idempotency_key=(
                        "partial-analysis-report:"
                        + toolchain.plan_sha256
                        + ":"
                        + round_tag
                    ),
                )
        return (
            tuple(settled[node_id] for node_id in sorted(settled)),
            analysis_status,
            completion_receipts,
            report_path,
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

    def _record_component_index(self) -> dict[str, int]:
        """Map every plan node to its record via the shared derivation.

        The record boundary is a connected component of the plan's own
        edges (``workflow_record_components``), never a stored field, so
        the walk, the review's per-record rows, and the delivery summary
        can never disagree about where one record ends.
        """

        from chemsmart.agent.workflows import workflow_record_components

        return {
            node_id: index
            for index, group in enumerate(
                workflow_record_components(
                    self.plan.nodes, getattr(self.plan, "edges", ())
                )
            )
            for node_id in group
        }

    def _analysis_calc_roots(self, toolchain: Any) -> dict[str, frozenset]:
        """Map each analysis node to the calculation nodes it rests on."""

        calc_ids = set(toolchain.calculation_node_ids)
        by_id = {node.node_id: node for node in toolchain.analysis_nodes}
        roots: dict[str, frozenset] = {}
        for node_id in toolchain.node_order:
            node = by_id.get(node_id)
            if node is None:
                continue
            collected: set[str] = set()
            for item in node.inputs:
                producer = getattr(item, "producer_node_id", "")
                if producer in calc_ids:
                    collected.add(producer)
                elif producer in roots:
                    collected |= set(roots[producer])
            roots[node.node_id] = frozenset(collected)
        return roots

    def _record_delivery_summary(
        self,
        executed: tuple[ExecutedNodeV1, ...],
        analysis_nodes: tuple[ExecutedAnalysisNodeV1, ...],
        toolchain: Any,
    ) -> tuple[dict[str, Any], ...]:
        """Per-record delivery verdicts for a plan that carries several.

        Verdict rules are fixed and small: a record's calculation is
        ``validated`` only when every node validated, ``not_delivered``
        when any attempted node did not validate, ``not_attempted`` when
        nothing ran (a suspended run's later records), and ``incomplete``
        otherwise; its analysis is judged over the chain nodes resting
        entirely on that record's calculations.  Analysis resting on
        several records or none (registered results) reports under
        ``shared``.
        """

        record_of = self._record_component_index()
        record_count = (max(record_of.values()) + 1) if record_of else 0
        if record_count <= 1:
            return ()
        outcome_by_id = {item.node_id: item for item in executed}
        roots = (
            self._analysis_calc_roots(toolchain)
            if toolchain is not None
            else {}
        )
        entries: list[dict[str, Any]] = []
        shared: list[ExecutedAnalysisNodeV1] = []
        claimed: set[str] = set()
        for index in range(record_count):
            node_ids = tuple(
                node.node_id
                for node in self.plan.nodes
                if record_of[node.node_id] == index
            )
            group = set(node_ids)
            states = {
                node_id: (
                    outcome_by_id[node_id].state
                    if node_id in outcome_by_id
                    else "not_attempted"
                )
                for node_id in node_ids
            }
            attempted = [
                node_id for node_id in node_ids if node_id in outcome_by_id
            ]
            if any(
                not outcome_by_id[node_id].validated for node_id in attempted
            ):
                calculation = "not_delivered"
            elif len(attempted) == len(node_ids) and attempted:
                calculation = "validated"
            elif not attempted:
                calculation = "not_attempted"
            else:
                calculation = "incomplete"
            mine = [
                record
                for record in analysis_nodes
                if roots.get(record.node_id)
                and set(roots[record.node_id]) <= group
            ]
            claimed.update(record.node_id for record in mine)
            if not mine:
                analysis = "none"
            elif all(record.state == "executed" for record in mine):
                analysis = "executed"
            elif all(
                record.state in {"skipped", "blocked_unsupported"}
                for record in mine
            ):
                analysis = "skipped"
            else:
                analysis = "partial"
            entries.append(
                {
                    "record": index + 1,
                    "node_ids": node_ids,
                    "node_states": states,
                    "calculation": calculation,
                    "analysis": analysis,
                }
            )
        shared = [
            record
            for record in analysis_nodes
            if record.node_id not in claimed
        ]
        if shared:
            entries.append(
                {
                    "record": "shared",
                    "node_ids": tuple(record.node_id for record in shared),
                    "node_states": {
                        record.node_id: record.state for record in shared
                    },
                    "calculation": "",
                    "analysis": (
                        "executed"
                        if all(record.state == "executed" for record in shared)
                        else "partial"
                    ),
                }
            )
        return tuple(entries)

    def _append_record_delivery_section(
        self,
        report_path: str,
        record_delivery: tuple[dict[str, Any], ...],
    ) -> None:
        """Render the per-record verdicts into the delivered report."""

        from chemsmart.agent.report_format import RECORD_DELIVERY_HEADING

        lines = [
            "",
            RECORD_DELIVERY_HEADING,
            "",
            "| Record | Nodes (reached state) | Calculation | Analysis |",
            "| --- | --- | --- | --- |",
        ]
        for entry in record_delivery:
            nodes_text = (
                ", ".join(
                    f"{node_id} ({entry['node_states'][node_id]})"
                    for node_id in entry["node_ids"]
                )
                or "-"
            )
            lines.append(
                f"| {entry['record']} | {nodes_text} | "
                f"{entry['calculation'] or '-'} | {entry['analysis']} |"
            )
        lines.append("")
        lines.append(
            "A batch of N is N observations; each verdict above is one "
            "record's, and no aggregate quantity is rendered."
        )
        with open(report_path, "a", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")

    def _replayed_outcome(
        self, node_state: Any, frontier: Any
    ) -> ExecutedNodeV1:
        """Report a durably settled node from its own receipts, unmoved.

        A continuation never re-drives the host sequence for a node the
        stream already settled: the preflight gate rightly refuses to
        re-prepare a non-pending node (the first live batch resume
        surfaced exactly that as six blocked rows), and the launch fence
        already guarantees no re-execution.  The durable receipt is the
        truth; the node's result artifacts are rehydrated from its own
        node workspace so the approved analysis chain can read them in
        this process -- the same scan that bound them when the engine
        finished.
        """

        binding = self._binding(node_state.node_id)
        receipt = frontier.receipt_for_node(node_state.node_id)
        if receipt is None:
            raise ContractError(
                f"durably settled node {node_state.node_id!r} holds no "
                "execution receipt; the stream is not continuable"
            )
        self.host.execution_receipts[node_state.node_id] = receipt
        node_workspace = (
            Path(self.approval_workspace) / "nodes" / node_state.node_id
        )
        if node_workspace.is_dir():
            self.host._execution_output_artifacts(
                node_state.node_id,
                node_workspace,
                program=binding.program,
            )
        validated = bool(getattr(receipt, "validated", False))
        return ExecutedNodeV1(
            node_id=node_state.node_id,
            program=binding.program,
            jobtype=binding.jobtype,
            state=str(getattr(receipt, "execution_state", node_state.state)),
            invocation_identity_sha256="",
            execution_receipt_sha256=str(
                getattr(receipt, "receipt_sha256", "")
            ),
            rule_ids=tuple(getattr(receipt, "findings", ()) or ()),
            failure=(
                ""
                if validated
                else _execution_failure_summary(receipt, self.host)
            ),
            validated=validated,
            result_validation_receipt_sha256=str(
                getattr(receipt, "result_validation_receipt_sha256", "")
            ),
            invocation_sha256=str(getattr(receipt, "invocation_sha256", "")),
        )

    def run(self) -> WorkflowExecutionResultV1:
        """Execute every ready node until the frontier stops advancing."""

        run_id = "run." + self.approval.approval_id
        executed: list[ExecutedNodeV1] = []
        seen: set[str] = set()
        # The approval, not the plan, names what may run.  A scientific plan
        # keeps a stage this release cannot execute rather than dropping it,
        # and such a stage carries no approved binding.
        approved_node_ids = {
            binding.node_id for binding in self.approval.node_bindings
        }
        durable = self.host.event_store.workflow_frontier(
            workflow_id=self.plan.workflow_id,
            run_id=run_id,
        )
        if durable.run_state is not None:
            # A continuation: the durable stream is the starting truth.
            # Admission happens here, at entry, not lazily at the first
            # launch -- a re-invocation whose remaining work is zero
            # launches must still be recorded as resumed, and a completed
            # approval must refuse before it re-delivers anything.
            if self.claim_workspace_bundle:
                from chemsmart.agent.live_session import (
                    continue_workflow_execution_approval_bundle,
                )

                continue_workflow_execution_approval_bundle(
                    self.execution_bundle,
                    workspace=self.approval_workspace,
                    run_event_store=self.host.event_store,
                )
                self._bundle_claimed = True
            run_state = durable.run_state
            data_edge_bindings = durable.data_edge_bindings
            for node_state in run_state.nodes:
                if node_state.node_id not in approved_node_ids:
                    continue
                if node_state.state in {
                    "validated",
                    "failed",
                    "engine_complete",
                }:
                    executed.append(
                        self._replayed_outcome(node_state, durable)
                    )
                    seen.add(node_state.node_id)
                elif node_state.state in {"blocked", "ambiguous"}:
                    binding = self._binding(node_state.node_id)
                    executed.append(
                        ExecutedNodeV1(
                            node_id=node_state.node_id,
                            program=binding.program,
                            jobtype=binding.jobtype,
                            state=node_state.state,
                            invocation_identity_sha256="",
                            execution_receipt_sha256="",
                            rule_ids=node_state.failure_rule_ids,
                            failure=(
                                "settled in a prior invocation of this "
                                "approval"
                            ),
                        )
                    )
                    seen.add(node_state.node_id)
                elif node_state.state == "running":
                    # A durable reservation with no receipt: the prior
                    # invocation died mid-engine.  Interrupted is a third
                    # state, and the default is refusal -- the fence
                    # forbids relaunching an engine whose original
                    # process and outputs nobody reconciled -- but the
                    # record must say so rather than vanish from the
                    # delivery table as never-attempted.  Observed live
                    # on the first mid-engine SIGTERM.
                    binding = self._binding(node_state.node_id)
                    executed.append(
                        ExecutedNodeV1(
                            node_id=node_state.node_id,
                            program=binding.program,
                            jobtype=binding.jobtype,
                            state="ambiguous",
                            invocation_identity_sha256="",
                            execution_receipt_sha256="",
                            # The typed id that distinguishes this from a
                            # timeout-ambiguous or a launch-ambiguous node;
                            # the English sentence below carried the whole
                            # distinction and was never persisted.
                            rule_ids=("execution.interrupted.mid_engine",),
                            failure=(
                                "launch reservation remains unresolved; "
                                "relaunch is forbidden -- a prior "
                                "invocation was interrupted mid-engine "
                                "and the original process and output "
                                "state need human reconciliation"
                            ),
                        )
                    )
                    seen.add(node_state.node_id)
        else:
            run_state = build_workflow_run_state(
                run_id=run_id,
                plan=self.plan,
                approval=self.frozen_approval,
                # This walk is the act of consuming the approval. An
                # unconsumed approval has an empty frontier by
                # construction, so leaving this false would report
                # "nothing to do" for a fully approved plan.
                approval_consumed=True,
            )
            data_edge_bindings = ()
        record_of = self._record_component_index()
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
            # Record-major order: a batch plans N records as N disconnected
            # sub-DAGs, and the bare frontier walk is breadth-first -- every
            # record's root before any record's second stage -- so a
            # suspended or partially failed run held N half-done records and
            # no finished one.  Completing the earliest incomplete record's
            # subgraph before opening the next root delivers whole records
            # early; a single-record plan is one component and unchanged.
            # Every attempted node enters ``seen``, so the walk terminates
            # without a progress check, and one record's failure leaves
            # later records' roots ready on the next pass.
            earliest = min(record_of[node_id] for node_id in ready)
            ready = tuple(
                node_id for node_id in ready if record_of[node_id] == earliest
            )
            cancelled = False
            for node_id in ready:
                # The human may cancel at a node boundary: an engine
                # already launched is never killed here -- that is what
                # timeouts and signals are for -- but no further node
                # launches under a grant the human has withdrawn.
                if self.should_stop is not None and self.should_stop():
                    binding = self._binding(node_id)
                    executed.append(
                        ExecutedNodeV1(
                            node_id=node_id,
                            program=binding.program,
                            jobtype=binding.jobtype,
                            state="cancelled",
                            invocation_identity_sha256="",
                            execution_receipt_sha256="",
                            rule_ids=("execution.cancelled.human",),
                            failure=(
                                "cancelled by the human at a node "
                                "boundary; the node was never launched"
                            ),
                        )
                    )
                    # The word must survive the process: an in-memory
                    # cancelled node derived as not_launched, so a
                    # withdrawn grant and a node that never came up
                    # were indistinguishable afterwards. A stop before
                    # the first launch leaves no durable run at all --
                    # there is nothing to mark in a stream that does
                    # not exist.
                    frontier = self.host.event_store.workflow_frontier(
                        workflow_id=self.plan.workflow_id,
                        run_id=run_id,
                    )
                    if frontier.run_state is not None:
                        _event, run_state = (
                            self.host.event_store.transition_workflow_run_node(
                                turn_id="execution-cancel",
                                run_id=run_id,
                                node_id=node_id,
                                new_state="cancelled",
                                plan=self.plan,
                                failure_rule_ids=(
                                    "execution.cancelled.human",
                                ),
                                timestamp=datetime.now(
                                    timezone.utc
                                ).isoformat(),
                            )
                        )
                    seen.add(node_id)
                    cancelled = True
                    continue
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
            if cancelled:
                break
        done = {item.node_id for item in executed if item.validated}
        status = "completed" if done == approved_node_ids else "partial"
        analysis_nodes: tuple[ExecutedAnalysisNodeV1, ...] = ()
        analysis_status = ""
        completion_receipts: tuple[str, ...] = ()
        report_path = ""
        toolchain = getattr(
            self.execution_bundle, "scientific_toolchain_plan", None
        )
        if toolchain is not None:
            # The chain walks whether the calculation partition completed or
            # not.  Under the old all-or-nothing gate one record's failed
            # SCF suppressed every other record's approved analysis; inside
            # the phase a node whose producer never validated settles as a
            # typed finding with the producer named, dependents skip, and
            # the partial envelope still renders every receipt that
            # survived -- the same uniform rule for one molecule or N.
            (
                analysis_nodes,
                analysis_status,
                completion_receipts,
                report_path,
            ) = self._run_analysis_phase(toolchain)
        record_delivery = self._record_delivery_summary(
            tuple(executed), analysis_nodes, toolchain
        )
        if record_delivery and report_path:
            self._append_record_delivery_section(report_path, record_delivery)
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
            analysis_nodes=analysis_nodes,
            analysis_status=analysis_status,
            analysis_completion_receipt_sha256s=completion_receipts,
            analysis_report_path=report_path,
            record_delivery=record_delivery,
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
    xtb_executable = os.environ.get(
        "CHEMSMART_XTB_EXECUTABLE"
    ) or shutil.which("xtb")
    executable_directory = (
        str(Path(xtb_executable).expanduser().parent) if xtb_executable else ""
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
        "run_evidence_root": workspace,
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
        "approved_scientific_toolchain_plan": getattr(
            bundle, "scientific_toolchain_plan", None
        ),
        # Without these the completion gate and the expectation rows read
        # an empty declaration set: they are recorded on the planning
        # session's tool host and evaluated on this one.
        "approved_requested_observable_declarations": getattr(
            bundle, "requested_observable_declarations", ()
        ),
        # The displayed budget must be the executing budget.  This
        # envelope was parsed above for the scratch root and then
        # dropped, so the provider-free executor ran with no episode
        # clock and no engine-call count -- only per-node timeouts.
        # Invisible in every one- and two-node workflow; the first
        # batch that outlived an episode launched all eight nodes
        # across thirty-one minutes of a displayed 1200 s episode.
        "bounded_execution_envelope": envelope,
    }


def _prior_anomalies(run_directory: Path) -> tuple[dict[str, Any], ...]:
    """Anomalies earlier cycles recorded, handed down through the run
    directory so this run's completion receipt can carry them."""

    from chemsmart.agent.terminal_states import PRIOR_ANOMALIES_FILE

    path = Path(run_directory) / PRIOR_ANOMALIES_FILE
    try:
        records = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ()
    return tuple(dict(item) for item in records if isinstance(item, Mapping))


def _execute_workflow_bundle(
    *,
    bundle: WorkflowExecutionApprovalBundleV1,
    workspace: Path,
    run_directory: Path,
    claim_workspace_bundle: bool,
    stop_file: Path | None = None,
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
        prior_anomaly_observations=_prior_anomalies(run_directory),
        # ``inputs`` already carries the execution server profile from the
        # bundle; naming it again here would be a second source of truth for
        # where real jobs are launched.
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
        should_stop=(
            (lambda: stop_file.exists()) if stop_file is not None else None
        ),
    ).run()


def execute_analysis_only_toolchain(
    *,
    host: Any,
    toolchain: Any,
    run_directory: Path,
    task_spec_sha256: str,
    workspace: Path,
) -> dict[str, Any]:
    """Walk a toolchain that plans no calculation, on the host that holds
    the results it reads.

    A woken session that had every number in hand planned a correct
    seventeen-node analysis-only chain naming the declared ids verbatim
    and then stopped, writing that approval was required before the
    analysis operations could execute. No review is built for a plan with
    zero calculation nodes, so the session ended ``planned`` and the goal
    returned to the human with sixteen engine calls unspent (NOVEL-2 po2,
    2026-09-04). Such a plan changes no identity, no state, no condition,
    and launches no engine: the goal's one decision already covers the
    results it reads. The owner ruled it is admitted and executed with no
    displayed review. The walker is the executor's own analysis phase --
    the same kernels, the same receipts, the same completion gate -- run
    without an approval bundle, because there is nothing to approve.
    """

    from types import SimpleNamespace

    from chemsmart.agent.runtime.events import EventKind

    if tuple(getattr(toolchain, "calculation_node_ids", ())):
        raise ContractError(
            "an analysis-only walk takes a toolchain with no calculation "
            "node; this plan names "
            f"{list(toolchain.calculation_node_ids)}"
        )
    if not tuple(getattr(toolchain, "analysis_nodes", ())):
        raise ContractError("an analysis-only walk needs an analysis node")
    # The completion renderer keys on the plan the host holds; the plan
    # tool stores it under its digest, and a walker reached another way
    # registers the same plan under the same digest.
    host.scientific_toolchain_plans.setdefault(
        toolchain.plan_sha256, toolchain
    )
    walker = ApprovedWorkflowExecutor(
        host=host,
        plan=SimpleNamespace(
            workflow_id=toolchain.workflow_id,
            plan_sha256=toolchain.plan_sha256,
            nodes=(),
        ),
        approval=SimpleNamespace(node_bindings=()),
        frozen_approval=SimpleNamespace(approval_sha256=""),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256=task_spec_sha256,
        run_directory=Path(run_directory),
        execution_bundle=SimpleNamespace(non_executable_node_ids=()),
        approval_workspace=Path(workspace),
        claim_workspace_bundle=False,
    )
    nodes, status, completions, report_path = walker._run_analysis_phase(
        toolchain
    )
    record = {
        "toolchain_plan_sha256": toolchain.plan_sha256,
        "workflow_id": toolchain.workflow_id,
        "analysis_status": status,
        "executed_nodes": tuple(
            {
                "node_id": node.node_id,
                "analysis_kind": node.analysis_kind,
                "state": node.state,
                "reason": node.reason,
            }
            for node in nodes
        ),
        "completion_receipt_sha256s": tuple(completions),
        "report_path": report_path,
        "engine_calls_consumed": 0,
    }
    host.event_store.append(
        turn_id="analysis-only-plan",
        kind=EventKind.ANALYSIS_ONLY_PLAN_EXECUTED.value,
        payload=record,
        idempotency_key="analysis-only-plan:" + toolchain.plan_sha256,
    )
    return record


def execute_approved_workflow(
    *,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    task_spec_sha256: str = "",
    expected_approval_file_sha256: str = "",
    stop_file: Path | None = None,
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
        stop_file=stop_file,
    )


__all__ = [
    "ApprovedWorkflowExecutor",
    "ExecutedNodeV1",
    "PROGRAM_NODE_SEQUENCE",
    "WorkflowExecutionResultV1",
    "execute_approved_workflow",
]
