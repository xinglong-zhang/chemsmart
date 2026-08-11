from __future__ import annotations

from copy import deepcopy

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.specialists import (
    BoundedSpecialistOrchestratorV1,
    READ_ONLY_CRITIC,
    SpecialistAdvisoryValidationError,
    SpecialistBudgetV1,
    build_specialist_session_response,
)
from chemsmart.agent.tool_specs import AgentToolSurfaceV1
from chemsmart.agent.workflows import (
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)


def _tool(name: str) -> dict:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": name,
            "parameters": {
                "type": "object",
                "properties": {},
                "additionalProperties": False,
            },
        },
    }


def _surface() -> AgentToolSurfaceV1:
    tools = tuple(
        _tool(name)
        for name in (
            "inspect_program_capability",
            "inspect_program_environment",
            "read_project_yaml",
            "validate_project_yaml",
            "inspect_calculation_artifact",
            "render_project_yaml",
            "repair_command",
            "execute_approved_program_node",
        )
    )
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="test-full-surface",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def _node(node_id: str, stage: str):
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="pyscf",
        program="pyscf",
        engine="cpu",
        project_role="water-project",
        unresolved_fields=(),
    )


def _plan():
    return build_scientific_workflow_plan(
        workflow_id="water-workflow",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            _node("opt-initial", "opt"),
            _node("hess-optimized", "hess"),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="opt-to-hess",
                source_node_id="opt-initial",
                target_node_id="hess-optimized",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
        ),
    )


class _FakeSession:
    def __init__(self, request, factory):
        self.session_id = request.session_id
        self.request = request
        self.factory = factory
        self.closed = False

    def run(self, request):
        assert request is self.request
        output = deepcopy(self.factory.outputs[request.role])
        tool_calls = request.context_manifest.allowed_tools[:1]
        return build_specialist_session_response(
            session_id=request.session_id,
            public_output=output,
            tool_calls=tool_calls,
            input_tokens=7,
            output_tokens=5,
            wall_time_millis=9,
        )

    def close(self):
        self.closed = True


class _FakeFactory:
    def __init__(self, *, critic_disposition="open"):
        self.requests = []
        self.sessions = []
        self.outputs = {
            "scientific-specialist": {
                "status": "complete",
                "public_summary": "B3LYP is the requested DFT method.",
                "field_proposals": [
                    {
                        "field_path": "project.sp.method",
                        "value": "b3lyp",
                        "evidence_sha256s": ["c" * 64],
                    }
                ],
                "unresolved_fields": [],
            },
            "pyscf-specialist": {
                "status": "complete",
                "public_summary": "PySCF uses the same typed method literal.",
                "field_proposals": [
                    {
                        "field_path": "project.sp.method",
                        "value": "b3lyp",
                        "evidence_sha256s": ["d" * 64],
                    },
                    {
                        "field_path": "project.sp.basis",
                        "value": "def2-svp",
                    },
                ],
                "unresolved_fields": [],
            },
            "dag-specialist": {
                "status": "complete",
                "public_summary": "The Hessian consumes the optimized geometry.",
                "field_proposals": [
                    {
                        "field_path": "workflow.edges",
                        "value": [
                            {
                                "source_node_id": "opt-initial",
                                "target_node_id": "hess-optimized",
                                "edge_kind": "data",
                            }
                        ],
                    }
                ],
                "unresolved_fields": [],
            },
            READ_ONLY_CRITIC: {
                "status": "complete",
                "public_summary": "The candidate preserves the requested method.",
                "findings": [
                    {
                        "rule_id": "method-preservation",
                        "severity": "info",
                        "expected": "B3LYP remains explicit.",
                        "observed": "B3LYP remains explicit.",
                        "disposition": critic_disposition,
                        "evidence_sha256s": ["e" * 64],
                    }
                ],
            },
        }

    def __call__(self, request):
        self.requests.append(request)
        session = _FakeSession(request, self)
        self.sessions.append(session)
        return session


def _orchestrator(factory):
    return BoundedSpecialistOrchestratorV1(
        base_tool_surface=_surface(),
        session_factory=factory,
        session_id_factory=(
            lambda role, serial: f"worker.{serial}.{role}"
        ),
    )


def _budget():
    return SpecialistBudgetV1(256, 2, 10)


def test_decomposition_uses_fresh_sessions_contexts_and_narrow_surfaces():
    factory = _FakeFactory()
    orchestrator = _orchestrator(factory)

    outcomes = orchestrator.run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Plan OPT then HESS for water."},
        source_sha256s=("c" * 64,),
        artifact_sha256s=("f" * 64,),
        roles=(
            "scientific_specialist",
            "pyscf-specialist",
            "dag-specialist",
        ),
        budget=_budget(),
    )

    assert len(outcomes) == 3
    assert len({request.session_id for request in factory.requests}) == 3
    assert "coordinator-session" not in {
        request.session_id for request in factory.requests
    }
    assert len(
        {
            request.context_manifest.manifest_sha256
            for request in factory.requests
        }
    ) == 3
    assert all(session.closed for session in factory.sessions)
    for request in factory.requests:
        assert set(request.context_manifest.allowed_tools).issubset(
            {
                "inspect_program_capability",
                "inspect_program_environment",
                "read_project_yaml",
                "validate_project_yaml",
                "inspect_calculation_artifact",
            }
        )
        assert "render_project_yaml" not in request.context_manifest.allowed_tools
        assert "execute_approved_program_node" not in (
            request.context_manifest.allowed_tools
        )
        with pytest.raises(ContractError, match="outside"):
            request.require_tool_allowed("repair_command")


def test_coordinator_merge_is_order_independent_and_retains_no_authority():
    factory = _FakeFactory()
    orchestrator = _orchestrator(factory)
    outcomes = orchestrator.run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Plan OPT then HESS for water."},
        roles=(
            "scientific-specialist",
            "pyscf-specialist",
            "dag-specialist",
        ),
        budget=_budget(),
    )

    forward = orchestrator.merge_before_coordinator(outcomes)
    reverse = orchestrator.merge_before_coordinator(tuple(reversed(outcomes)))

    assert forward.receipt_sha256 == reverse.receipt_sha256
    assert forward.status == "merged"
    assert tuple(item.field_path for item in forward.merged_fields) == (
        "project.sp.basis",
        "project.sp.method",
        "workflow.edges",
    )
    rendered = repr(forward)
    assert "execute_approved_program_node" not in rendered
    assert "terminal_state" not in rendered


def test_coordinator_merge_reports_conflict_instead_of_silently_selecting():
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"][0][
        "value"
    ] = "pbe0"
    orchestrator = _orchestrator(factory)
    outcomes = orchestrator.run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Plan OPT then HESS for water."},
        roles=("scientific-specialist", "pyscf-specialist"),
        budget=_budget(),
    )

    receipt = orchestrator.merge_before_coordinator(outcomes)

    assert receipt.status == "needs_clarification"
    assert tuple(item.field_path for item in receipt.conflicts) == (
        "project.sp.method",
    )
    assert "project.sp.method" in receipt.unresolved_fields


def test_specialist_cannot_transfer_shell_or_host_authority():
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"] = [
        {
            "field_path": "project.sp.command",
            "value": "python run.py",
        }
    ]
    orchestrator = _orchestrator(factory)

    with pytest.raises(ContractError, match="host authority"):
        orchestrator.run_before_coordinator(
            plan=_plan(),
            coordinator_session_id="coordinator-session",
            public_context={"task": "Plan OPT then HESS for water."},
            roles=("scientific-specialist",),
            budget=_budget(),
        )


def test_specialist_scientific_prose_allows_semicolon_punctuation():
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"] = [
        {
            "field_path": "scientific.suitability_assessment",
            "value": (
                "Closed-shell DFT is suitable for this state; exact numerical "
                "accuracy remains unverified."
            ),
        }
    ]
    outcomes = _orchestrator(factory).run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Assess method suitability."},
        roles=("scientific-specialist",),
        budget=_budget(),
    )

    assert outcomes[0].candidate.proposals[0].value.endswith("unverified.")


def test_specialist_accepts_observation_vocabulary_without_authority_transfer():
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"] = [
        {
            "field_path": "scientific.environment.pyscf_cpu.environment_observation",
            "value": {
                "observed_status": "available",
                "evidence_kind": "typed_host_compute_receipt",
            },
            "evidence_sha256s": ["c" * 64],
        }
    ]

    outcomes = _orchestrator(factory).run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Report the observed PySCF environment."},
        roles=("scientific-specialist",),
        budget=_budget(),
    )

    proposal = outcomes[0].candidate.proposals[0]
    assert proposal.field_path.endswith("environment_observation")
    assert proposal.value["observed_status"] == "available"


@pytest.mark.parametrize("authority_segment", ("readiness", "execution_ready"))
def test_specialist_still_rejects_readiness_authority_segment(
    authority_segment,
):
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"] = [
        {
            "field_path": (
                "scientific.environment.pyscf_cpu." + authority_segment
            ),
            "value": "environment-ready",
        }
    ]

    with pytest.raises(
        SpecialistAdvisoryValidationError,
        match="host authority",
    ):
        _orchestrator(factory).run_before_coordinator(
            plan=_plan(),
            coordinator_session_id="coordinator-session",
            public_context={"task": "Report the observed PySCF environment."},
            roles=("scientific-specialist",),
            budget=_budget(),
        )


def test_specialist_prompt_names_safe_observation_vocabulary():
    factory = _FakeFactory()
    _orchestrator(factory).run_before_coordinator(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        public_context={"task": "Assess method suitability."},
        roles=("scientific-specialist",),
        budget=_budget(),
    )

    prompt = factory.requests[0].system_instruction
    assert "environment_observation" in prompt
    assert "suitability_assessment" in prompt
    assert "Host-returned authority values" in prompt
    assert "must not be recast as a specialist decision" in prompt


@pytest.mark.parametrize(
    ("value", "rule_id"),
    (
        (
            "python run.py; echo finished",
            "specialist.advisory.shell-authority.v1",
        ),
        (
            "Review the generated file at /tmp/pyscf.py",
            "specialist.advisory.filesystem-authority.v1",
        ),
        (
            "from pyscf import gto\nmol = gto.M(atom='H 0 0 0')",
            "specialist.advisory.native-input.v1",
        ),
    ),
)
def test_specialist_rejects_executable_authority_with_role_local_finding(
    value,
    rule_id,
):
    factory = _FakeFactory()
    factory.outputs["scientific-specialist"]["field_proposals"] = [
        {
            "field_path": "scientific.suitability_assessment",
            "value": value,
        }
    ]

    with pytest.raises(SpecialistAdvisoryValidationError) as observed:
        _orchestrator(factory).run_before_coordinator(
            plan=_plan(),
            coordinator_session_id="coordinator-session",
            public_context={"task": "Assess method suitability."},
            roles=("scientific-specialist",),
            budget=_budget(),
        )

    finding = observed.value.public_finding()
    assert finding["role"] == "scientific-specialist"
    assert finding["rule_id"] == rule_id
    assert finding["disposition"] == "rejected"
    assert finding["finding_sha256"] == canonical_sha256(
        {key: value for key, value in finding.items() if key != "finding_sha256"}
    )


def test_specialist_unresolved_fields_reject_objects_with_local_feedback():
    factory = _FakeFactory()
    factory.outputs["dag-specialist"]["unresolved_fields"] = [
        {
            "field_path": "workflow.nodes.hess.geometry_input",
            "reason": "The optimized geometry does not exist yet.",
        }
    ]
    orchestrator = _orchestrator(factory)

    with pytest.raises(
        ContractError,
        match=r"unresolved_fields\[0\] must be a lower-case field-path string",
    ):
        orchestrator.run_before_coordinator(
            plan=_plan(),
            coordinator_session_id="coordinator-session",
            public_context={"task": "Plan OPT then HESS for water."},
            roles=("dag-specialist",),
            budget=_budget(),
        )


def test_fresh_critic_is_read_only_and_candidate_is_immutable():
    factory = _FakeFactory()
    orchestrator = _orchestrator(factory)
    candidate = {
        "project": {"sp": {"method": "b3lyp", "basis": "def2-svp"}},
        "workflow": ["opt-initial", "hess-optimized"],
    }
    before = deepcopy(candidate)

    outcome = orchestrator.run_after_coordinator_critic(
        plan=_plan(),
        coordinator_session_id="coordinator-session",
        candidate_id="water-candidate",
        candidate_record=candidate,
        public_context={"question": "Does the plan preserve the method?"},
        budget=_budget(),
    )

    request = factory.requests[-1]
    assert request.role == READ_ONLY_CRITIC
    assert request.session_id != "coordinator-session"
    assert set(request.context_manifest.allowed_tools).issubset(
        {
            "inspect_program_capability",
            "inspect_program_environment",
            "read_project_yaml",
            "validate_project_yaml",
            "inspect_calculation_artifact",
        }
    )
    assert "repair_command" not in request.context_manifest.allowed_tools
    assert "execute_approved_program_node" not in (
        request.context_manifest.allowed_tools
    )
    assert candidate == before
    assert outcome.candidate_sha256_before == outcome.candidate_sha256_after
    assert len(outcome.review.findings) == 1
    assert outcome.review.findings[0].disposition == "open"


def test_critic_cannot_resolve_its_own_finding():
    factory = _FakeFactory(critic_disposition="resolved")
    orchestrator = _orchestrator(factory)

    with pytest.raises(ContractError, match="cannot dispose"):
        orchestrator.run_after_coordinator_critic(
            plan=_plan(),
            coordinator_session_id="coordinator-session",
            candidate_id="water-candidate",
            candidate_record={"project": {"method": "b3lyp"}},
            public_context={"question": "Review the candidate."},
            budget=_budget(),
        )
