"""Bounded specialist and critic orchestration for scientific workflows.

This module is deliberately provider neutral.  A caller supplies a factory
that creates one disposable provider session per request.  The dispatcher
owns immutable context packets, role-specific read-only tool surfaces, budget
checks, typed result parsing, and deterministic coordinator merge semantics.

Specialists can propose scientific fields; they cannot transfer paths, shell
syntax, approval, readiness, execution, or terminal-state authority.  A critic
is always a fresh read-only session and can only return open findings.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
from typing import Any, Callable, Mapping, Protocol, Sequence
from uuid import uuid4

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
    require_identifier,
    require_sha256,
)
from chemsmart.agent.tool_specs import AgentToolSurfaceV1
from chemsmart.agent.workflows import (
    ContextManifestV1,
    ScientificReviewFindingV1,
    ScientificWorkflowPlanV2,
    SpecialistResultPacketV1,
    SpecialistTaskPacketV1,
    eligible_specialist_roles,
)


SCIENTIFIC_SPECIALIST = "scientific-specialist"
PYSCF_SPECIALIST = "pyscf-specialist"
DAG_SPECIALIST = "dag-specialist"
READ_ONLY_CRITIC = "critic"

_ROLE_ALIASES = {
    "scientific-specialist": SCIENTIFIC_SPECIALIST,
    "scientific_specialist": SCIENTIFIC_SPECIALIST,
    "pyscf-specialist": PYSCF_SPECIALIST,
    "pyscf_specialist": PYSCF_SPECIALIST,
    "dag-specialist": DAG_SPECIALIST,
    "dag_specialist": DAG_SPECIALIST,
    "critic": READ_ONLY_CRITIC,
    "fresh-read-only-critic": READ_ONLY_CRITIC,
}

_READ_ONLY_TOOLS = frozenset(
    {
        "inspect_program_capability",
        "inspect_program_environment",
        "read_project_yaml",
        "validate_project_yaml",
        "inspect_calculation_artifact",
    }
)
_ROLE_TOOLS = {
    SCIENTIFIC_SPECIALIST: frozenset(
        {
            "inspect_program_capability",
            "inspect_program_environment",
            "read_project_yaml",
            "validate_project_yaml",
            "inspect_calculation_artifact",
        }
    ),
    PYSCF_SPECIALIST: frozenset(
        {
            "inspect_program_capability",
            "inspect_program_environment",
            "read_project_yaml",
            "validate_project_yaml",
            "inspect_calculation_artifact",
        }
    ),
    DAG_SPECIALIST: frozenset(
        {
            "inspect_program_capability",
            "read_project_yaml",
        }
    ),
    READ_ONLY_CRITIC: _READ_ONLY_TOOLS,
}
_ROLE_FIELD_PREFIXES = {
    SCIENTIFIC_SPECIALIST: ("identity.", "project.", "scientific."),
    PYSCF_SPECIALIST: ("project.", "scientific."),
    DAG_SPECIALIST: ("workflow.",),
}
_FORBIDDEN_AUTHORITY_SEGMENTS = frozenset(
    {
        "approval",
        "approved",
        "argv",
        "command",
        "cwd",
        "execute",
        "executable",
        "execution",
        "execution_ready",
        "execution_state",
        "path",
        "readiness",
        "script",
        "shell",
        "success",
        "terminal",
        "terminal_state",
    }
)
_PRIVATE_REASONING_KEYS = frozenset(
    {"reasoning_content", "thinking", "chain_of_thought", "hidden_reasoning"}
)
_ABSOLUTE_PATH_VALUE = re.compile(r"(?:^|[\s\"'])(?:~/|/[A-Za-z0-9._-])")
_WINDOWS_ABSOLUTE_PATH = re.compile(r"(?:^|[\s\"'])[a-zA-Z]:[\\/]")
_CLEAR_EXECUTABLE_PAYLOAD = re.compile(
    r"""(?ix)
    (?:^|[;\r\n])\s*(?:
        python(?:3(?:\.\d+)*)?\s+(?:-[cm]\b|[^\s;]+\.py\b)
        |(?:bash|sh|zsh)\s+(?:-[a-z]+\b|[^\s;]+)
        |chemsmart\s+(?:run|sub|agent)\b
        |(?:xtb|orca|g16)\s+(?:--?[a-z]|[^\s;]+\.(?:xyz|inp|com|gjf)\b)
        |(?:rm|mv|cp|curl|wget|chmod|chown|sbatch|srun|qsub|echo|touch)\s+\S
        |export\s+[A-Za-z_][A-Za-z0-9_]*=
        |(?:cd|source|exec|sudo)\s+\S
        |env\s+(?:-[A-Za-z]+\s+)*\S
    )
    """
)
_CLEAR_NATIVE_INPUT_PAYLOAD = re.compile(
    r"(?im)^\s*(?:#!\s*/|from\s+pyscf\b|import\s+pyscf\b|%pal\b|"
    r"\*\s+xyz(?:file)?\b|\$coord\b|#p\s+\S+|!\s+\S+)"
)


class SpecialistAdvisoryValidationError(ContractError):
    """Publicly classifiable rejection of one role's untrusted advisory."""

    def __init__(self, *, role: str, rule_id: str, detail: str = "") -> None:
        self.role = _normalized_role(role)
        self.rule_id = str(rule_id)
        super().__init__(
            detail
            or f"{self.role} output was rejected by the typed advisory contract"
        )

    def public_finding(self) -> dict[str, Any]:
        body = {
            "schema_version": (
                "chemsmart.specialist-output-validation-finding.v1"
            ),
            "role": self.role,
            "rule_id": self.rule_id,
            "severity": "error",
            "disposition": "rejected",
            "public_summary": (
                "The role-local advisory was rejected before coordinator merge."
            ),
        }
        return {**body, "finding_sha256": canonical_sha256(body)}


def _normalized_role(role: str) -> str:
    normalized = str(role or "").strip().lower()
    try:
        return _ROLE_ALIASES[normalized]
    except KeyError as exc:
        raise ContractError(f"unsupported specialist role: {role!r}") from exc


def _tool_name(definition: Mapping[str, Any]) -> str:
    function = definition.get("function")
    if not isinstance(function, Mapping):
        raise ContractError("tool definition requires a function object")
    return require_identifier(str(function.get("name") or ""), "tool_name")


def build_specialist_tool_surface(
    *, role: str, base_surface: AgentToolSurfaceV1
) -> AgentToolSurfaceV1:
    """Project a normal tool surface to the role's least-privilege subset."""

    normalized_role = _normalized_role(role)
    permitted = _ROLE_TOOLS[normalized_role]
    tools = tuple(
        canonical_data(definition)
        for definition in base_surface.tool_definitions
        if _tool_name(definition) in permitted
    )
    observed_names = tuple(sorted(_tool_name(item) for item in tools))
    if set(observed_names).difference(_READ_ONLY_TOOLS):
        raise ContractError("specialist surface contains a mutating tool")
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile=f"{normalized_role}-read-only",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


@dataclass(frozen=True)
class SpecialistBudgetV1:
    token_budget: int
    tool_call_budget: int
    wall_time_seconds: int

    def __post_init__(self) -> None:
        if min(
            self.token_budget, self.tool_call_budget, self.wall_time_seconds
        ) < 1:
            raise ContractError("specialist budgets must be positive")


@dataclass(frozen=True)
class SpecialistSessionRequestV1:
    """Immutable input to one fresh provider session."""

    schema_version: str
    session_id: str
    coordinator_session_id: str
    role: str
    context_manifest: ContextManifestV1
    task_packet: SpecialistTaskPacketV1
    system_instruction: str
    public_context_json: str
    tool_definitions_json: str
    request_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-session-request.v1":
            raise ContractError("unsupported specialist session request")
        require_identifier(self.session_id, "session_id")
        if not str(self.coordinator_session_id).strip():
            raise ContractError("coordinator session ID is required")
        if self.session_id == self.coordinator_session_id:
            raise ContractError("specialist must use a fresh provider session")
        normalized_role = _normalized_role(self.role)
        if normalized_role != self.role:
            raise ContractError("specialist request role is not canonical")
        if self.task_packet.role != self.role:
            raise ContractError("specialist request and task packet roles differ")
        if (
            self.task_packet.context_manifest_sha256
            != self.context_manifest.manifest_sha256
        ):
            raise ContractError("task packet uses another context manifest")
        if not self.system_instruction.strip():
            raise ContractError("specialist system instruction is required")
        public_context = _canonical_json_record(
            self.public_context_json, "public context"
        )
        if not isinstance(public_context, Mapping):
            raise ContractError("public specialist context must be an object")
        _reject_private_reasoning(public_context)
        tool_definitions = _canonical_json_record(
            self.tool_definitions_json, "tool definitions"
        )
        if not isinstance(tool_definitions, list):
            raise ContractError("specialist tools must be an array")
        tool_names = tuple(sorted(_tool_name(item) for item in tool_definitions))
        if tool_names != self.context_manifest.allowed_tools:
            raise ContractError("context tools differ from provider tool schema")
        if canonical_sha256(tool_definitions) != (
            self.context_manifest.tool_schema_sha256
        ):
            raise ContractError("context tool digest differs from provider schema")
        body = {
            "schema_version": self.schema_version,
            "session_id": self.session_id,
            "coordinator_session_id": self.coordinator_session_id,
            "role": self.role,
            "context_manifest": self.context_manifest,
            "task_packet": self.task_packet,
            "system_instruction": self.system_instruction,
            "public_context_json": self.public_context_json,
            "tool_definitions_json": self.tool_definitions_json,
        }
        if self.request_sha256 != canonical_sha256(body):
            raise ContractError("specialist session request digest mismatch")

    @property
    def tool_definitions(self) -> tuple[dict[str, Any], ...]:
        decoded = json.loads(self.tool_definitions_json)
        return tuple(decoded)

    def require_tool_allowed(self, tool_name: str) -> str:
        normalized = require_identifier(tool_name, "tool_name")
        if normalized not in self.context_manifest.allowed_tools:
            raise ContractError(
                f"tool {normalized!r} is outside the specialist surface"
            )
        return normalized


@dataclass(frozen=True)
class SpecialistSessionResponseV1:
    """Public, reasoning-free result returned by a disposable session."""

    schema_version: str
    session_id: str
    public_output_json: str
    tool_calls: tuple[str, ...]
    input_tokens: int
    output_tokens: int
    wall_time_millis: int
    response_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-session-response.v1":
            raise ContractError("unsupported specialist session response")
        require_identifier(self.session_id, "session_id")
        output = _canonical_json_record(self.public_output_json, "public output")
        if not isinstance(output, Mapping):
            raise ContractError("specialist public output must be an object")
        _reject_private_reasoning(output)
        for tool_name in self.tool_calls:
            require_identifier(tool_name, "tool_call")
        if min(self.input_tokens, self.output_tokens, self.wall_time_millis) < 0:
            raise ContractError("specialist usage cannot be negative")
        body = {
            "schema_version": self.schema_version,
            "session_id": self.session_id,
            "public_output_json": self.public_output_json,
            "tool_calls": self.tool_calls,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "wall_time_millis": self.wall_time_millis,
        }
        if self.response_sha256 != canonical_sha256(body):
            raise ContractError("specialist session response digest mismatch")

    @property
    def public_output(self) -> dict[str, Any]:
        return json.loads(self.public_output_json)


def build_specialist_session_response(
    *,
    session_id: str,
    public_output: Mapping[str, Any],
    tool_calls: Sequence[str] = (),
    input_tokens: int = 0,
    output_tokens: int = 0,
    wall_time_millis: int = 0,
) -> SpecialistSessionResponseV1:
    body = {
        "schema_version": "chemsmart.specialist-session-response.v1",
        "session_id": require_identifier(session_id, "session_id"),
        "public_output_json": canonical_json(public_output),
        "tool_calls": tuple(str(item) for item in tool_calls),
        "input_tokens": int(input_tokens),
        "output_tokens": int(output_tokens),
        "wall_time_millis": int(wall_time_millis),
    }
    return SpecialistSessionResponseV1(
        **body, response_sha256=canonical_sha256(body)
    )


class SpecialistProviderSessionV1(Protocol):
    session_id: str

    def run(
        self, request: SpecialistSessionRequestV1
    ) -> SpecialistSessionResponseV1: ...

    def close(self) -> None: ...


class SpecialistSessionFactoryV1(Protocol):
    def __call__(
        self, request: SpecialistSessionRequestV1
    ) -> SpecialistProviderSessionV1: ...


@dataclass(frozen=True)
class SpecialistFieldProposalV1:
    schema_version: str
    proposal_id: str
    role: str
    field_path: str
    value_json: str
    evidence_sha256s: tuple[str, ...]
    uncertainty: str
    proposal_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-field-proposal.v1":
            raise ContractError("unsupported specialist field proposal")
        require_identifier(self.proposal_id, "proposal_id")
        normalized_role = _normalized_role(self.role)
        if normalized_role != self.role or self.role == READ_ONLY_CRITIC:
            raise ContractError("critic cannot emit field proposals")
        require_identifier(self.field_path, "field_path")
        if not self.field_path.startswith(_ROLE_FIELD_PREFIXES[self.role]):
            raise ContractError("field is outside specialist role ownership")
        _reject_authority_path(self.field_path)
        value = _canonical_json_record(self.value_json, "proposal value")
        _reject_authority_value(value)
        if self.evidence_sha256s != tuple(sorted(set(self.evidence_sha256s))):
            raise ContractError("proposal evidence must be sorted and unique")
        for digest in self.evidence_sha256s:
            require_sha256(digest, "evidence_sha256")
        body = {
            "schema_version": self.schema_version,
            "proposal_id": self.proposal_id,
            "role": self.role,
            "field_path": self.field_path,
            "value_json": self.value_json,
            "evidence_sha256s": self.evidence_sha256s,
            "uncertainty": self.uncertainty,
        }
        if self.proposal_sha256 != canonical_sha256(body):
            raise ContractError("specialist field proposal digest mismatch")

    @property
    def value(self) -> Any:
        return json.loads(self.value_json)


@dataclass(frozen=True)
class SpecialistCandidateV1:
    schema_version: str
    packet_id: str
    workflow_id: str
    role: str
    proposals: tuple[SpecialistFieldProposalV1, ...]
    unresolved_fields: tuple[str, ...]
    status: str
    candidate_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-candidate.v1":
            raise ContractError("unsupported specialist candidate")
        require_identifier(self.packet_id, "packet_id")
        require_identifier(self.workflow_id, "workflow_id")
        normalized_role = _normalized_role(self.role)
        if normalized_role != self.role or self.role == READ_ONLY_CRITIC:
            raise ContractError("critic cannot emit a coordinator candidate")
        order = tuple((item.field_path, item.proposal_sha256) for item in self.proposals)
        if order != tuple(sorted(order)):
            raise ContractError("specialist proposals are not canonical")
        paths = tuple(item.field_path for item in self.proposals)
        if len(paths) != len(set(paths)):
            raise ContractError("specialist candidate repeats a field")
        if any(item.role != self.role for item in self.proposals):
            raise ContractError("candidate contains another role's proposal")
        if self.unresolved_fields != tuple(sorted(set(self.unresolved_fields))):
            raise ContractError("candidate unresolved fields are not canonical")
        for field_path in self.unresolved_fields:
            require_identifier(field_path, "unresolved_field")
            _reject_authority_path(field_path)
        if self.status not in {"complete", "blocked", "failed"}:
            raise ContractError("unsupported specialist candidate status")
        if self.status == "failed" and self.proposals:
            raise ContractError("failed specialist cannot transfer proposals")
        body = {
            "schema_version": self.schema_version,
            "packet_id": self.packet_id,
            "workflow_id": self.workflow_id,
            "role": self.role,
            "proposals": self.proposals,
            "unresolved_fields": self.unresolved_fields,
            "status": self.status,
        }
        if self.candidate_sha256 != canonical_sha256(body):
            raise ContractError("specialist candidate digest mismatch")


@dataclass(frozen=True)
class SpecialistSessionOutcomeV1:
    schema_version: str
    request_sha256: str
    response_sha256: str
    result_packet: SpecialistResultPacketV1
    candidate: SpecialistCandidateV1
    outcome_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.specialist-session-outcome.v1":
            raise ContractError("unsupported specialist outcome")
        require_sha256(self.request_sha256, "request_sha256")
        require_sha256(self.response_sha256, "response_sha256")
        if self.result_packet.output_record_sha256 != self.candidate.candidate_sha256:
            raise ContractError("result packet does not bind specialist candidate")
        if self.result_packet.packet_id != self.candidate.packet_id:
            raise ContractError("result and candidate packets differ")
        body = {
            "schema_version": self.schema_version,
            "request_sha256": self.request_sha256,
            "response_sha256": self.response_sha256,
            "result_packet": self.result_packet,
            "candidate": self.candidate,
        }
        if self.outcome_sha256 != canonical_sha256(body):
            raise ContractError("specialist outcome digest mismatch")


@dataclass(frozen=True)
class MergedFieldProposalV1:
    schema_version: str
    field_path: str
    value_json: str
    source_proposal_sha256s: tuple[str, ...]
    source_roles: tuple[str, ...]
    evidence_sha256s: tuple[str, ...]
    merged_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.merged-field-proposal.v1":
            raise ContractError("unsupported merged field proposal")
        require_identifier(self.field_path, "field_path")
        _reject_authority_path(self.field_path)
        _reject_authority_value(
            _canonical_json_record(self.value_json, "merged value")
        )
        for values, name in (
            (self.source_proposal_sha256s, "source proposal"),
            (self.evidence_sha256s, "evidence"),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} digests are not canonical")
            for digest in values:
                require_sha256(digest, name)
        if self.source_roles != tuple(sorted(set(self.source_roles))):
            raise ContractError("merged source roles are not canonical")
        body = {
            "schema_version": self.schema_version,
            "field_path": self.field_path,
            "value_json": self.value_json,
            "source_proposal_sha256s": self.source_proposal_sha256s,
            "source_roles": self.source_roles,
            "evidence_sha256s": self.evidence_sha256s,
        }
        if self.merged_sha256 != canonical_sha256(body):
            raise ContractError("merged field proposal digest mismatch")


@dataclass(frozen=True)
class CoordinatorMergeConflictV1:
    schema_version: str
    field_path: str
    proposal_sha256s: tuple[str, ...]
    value_sha256s: tuple[str, ...]
    conflict_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.coordinator-merge-conflict.v1":
            raise ContractError("unsupported coordinator merge conflict")
        require_identifier(self.field_path, "field_path")
        for values, name in (
            (self.proposal_sha256s, "proposal"),
            (self.value_sha256s, "value"),
        ):
            if values != tuple(sorted(set(values))) or len(values) < 2:
                raise ContractError(f"merge conflict requires distinct {name}s")
            for digest in values:
                require_sha256(digest, name)
        body = {
            "schema_version": self.schema_version,
            "field_path": self.field_path,
            "proposal_sha256s": self.proposal_sha256s,
            "value_sha256s": self.value_sha256s,
        }
        if self.conflict_sha256 != canonical_sha256(body):
            raise ContractError("coordinator merge conflict digest mismatch")


@dataclass(frozen=True)
class CoordinatorMergeReceiptV1:
    schema_version: str
    workflow_id: str
    coordinator_owner: str
    source_result_sha256s: tuple[str, ...]
    merged_fields: tuple[MergedFieldProposalV1, ...]
    conflicts: tuple[CoordinatorMergeConflictV1, ...]
    unresolved_fields: tuple[str, ...]
    noncomplete_result_sha256s: tuple[str, ...]
    status: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.coordinator-merge-receipt.v1":
            raise ContractError("unsupported coordinator merge receipt")
        require_identifier(self.workflow_id, "workflow_id")
        if self.coordinator_owner != "coordinator":
            raise ContractError("only coordinator can own a specialist merge")
        for values, name in (
            (self.source_result_sha256s, "source result"),
            (self.noncomplete_result_sha256s, "noncomplete result"),
        ):
            if values != tuple(sorted(set(values))):
                raise ContractError(f"{name} digests are not canonical")
            for digest in values:
                require_sha256(digest, name)
        if tuple(item.field_path for item in self.merged_fields) != tuple(
            sorted(item.field_path for item in self.merged_fields)
        ):
            raise ContractError("merged fields are not canonical")
        if tuple(item.field_path for item in self.conflicts) != tuple(
            sorted(item.field_path for item in self.conflicts)
        ):
            raise ContractError("merge conflicts are not canonical")
        if self.unresolved_fields != tuple(sorted(set(self.unresolved_fields))):
            raise ContractError("merge unresolved fields are not canonical")
        if self.status not in {"merged", "needs_clarification", "no_candidate"}:
            raise ContractError("unsupported coordinator merge status")
        if self.conflicts and self.status != "needs_clarification":
            raise ContractError("conflicted merge must need clarification")
        body = {
            "schema_version": self.schema_version,
            "workflow_id": self.workflow_id,
            "coordinator_owner": self.coordinator_owner,
            "source_result_sha256s": self.source_result_sha256s,
            "merged_fields": self.merged_fields,
            "conflicts": self.conflicts,
            "unresolved_fields": self.unresolved_fields,
            "noncomplete_result_sha256s": self.noncomplete_result_sha256s,
            "status": self.status,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("coordinator merge receipt digest mismatch")


@dataclass(frozen=True)
class CriticReviewRecordV1:
    schema_version: str
    candidate_id: str
    candidate_sha256: str
    findings: tuple[ScientificReviewFindingV1, ...]
    status: str
    review_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.critic-review-record.v1":
            raise ContractError("unsupported critic review record")
        require_identifier(self.candidate_id, "candidate_id")
        require_sha256(self.candidate_sha256, "candidate_sha256")
        order = tuple(item.finding_id for item in self.findings)
        if order != tuple(sorted(set(order))):
            raise ContractError("critic findings are not canonical")
        if any(item.disposition != "open" for item in self.findings):
            raise ContractError("read-only critic can only emit open findings")
        if self.status not in {"complete", "blocked", "failed"}:
            raise ContractError("unsupported critic review status")
        if self.status == "failed" and self.findings:
            raise ContractError("failed critic cannot transfer findings")
        body = {
            "schema_version": self.schema_version,
            "candidate_id": self.candidate_id,
            "candidate_sha256": self.candidate_sha256,
            "findings": self.findings,
            "status": self.status,
        }
        if self.review_sha256 != canonical_sha256(body):
            raise ContractError("critic review digest mismatch")


@dataclass(frozen=True)
class CriticSessionOutcomeV1:
    schema_version: str
    request_sha256: str
    response_sha256: str
    result_packet: SpecialistResultPacketV1
    review: CriticReviewRecordV1
    candidate_sha256_before: str
    candidate_sha256_after: str
    outcome_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.critic-session-outcome.v1":
            raise ContractError("unsupported critic outcome")
        for digest, name in (
            (self.request_sha256, "request_sha256"),
            (self.response_sha256, "response_sha256"),
            (self.candidate_sha256_before, "candidate_sha256_before"),
            (self.candidate_sha256_after, "candidate_sha256_after"),
        ):
            require_sha256(digest, name)
        if self.candidate_sha256_before != self.candidate_sha256_after:
            raise ContractError("critic mutated the reviewed candidate")
        if self.review.candidate_sha256 != self.candidate_sha256_before:
            raise ContractError("critic review targets another candidate")
        if self.result_packet.role != READ_ONLY_CRITIC:
            raise ContractError("critic result packet has another role")
        if self.result_packet.output_record_sha256 != self.review.review_sha256:
            raise ContractError("critic result packet does not bind the review")
        body = {
            "schema_version": self.schema_version,
            "request_sha256": self.request_sha256,
            "response_sha256": self.response_sha256,
            "result_packet": self.result_packet,
            "review": self.review,
            "candidate_sha256_before": self.candidate_sha256_before,
            "candidate_sha256_after": self.candidate_sha256_after,
        }
        if self.outcome_sha256 != canonical_sha256(body):
            raise ContractError("critic outcome digest mismatch")


class BoundedSpecialistOrchestratorV1:
    """Composition boundary called before and after the coordinator session."""

    def __init__(
        self,
        *,
        base_tool_surface: AgentToolSurfaceV1,
        session_factory: SpecialistSessionFactoryV1,
        session_id_factory: Callable[[str, int], str] | None = None,
    ) -> None:
        self._base_tool_surface = base_tool_surface
        self._session_factory = session_factory
        self._session_id_factory = session_id_factory or _default_session_id
        self._serial = 0
        self._used_session_ids: set[str] = set()
        self._used_context_sha256s: set[str] = set()

    def run_before_coordinator(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        coordinator_session_id: str,
        public_context: Mapping[str, Any],
        source_sha256s: Sequence[str] = (),
        artifact_sha256s: Sequence[str] = (),
        input_record_sha256s: Sequence[str] = (),
        roles: Sequence[str] | None = None,
        budget: SpecialistBudgetV1,
    ) -> tuple[SpecialistSessionOutcomeV1, ...]:
        """Run independent workers before the coordinator integrates a plan."""

        selected = roles
        if selected is None:
            selected = eligible_specialist_roles(plan)
        normalized_roles = tuple(sorted({_normalized_role(role) for role in selected}))
        if READ_ONLY_CRITIC in normalized_roles:
            raise ContractError("critic runs only after the coordinator candidate")
        outcomes = []
        for role in normalized_roles:
            request = self._build_request(
                plan=plan,
                coordinator_session_id=coordinator_session_id,
                role=role,
                public_context={
                    "workflow_plan": canonical_data(plan),
                    "task_context": canonical_data(public_context),
                },
                source_sha256s=source_sha256s,
                artifact_sha256s=artifact_sha256s,
                input_record_sha256s=(plan.plan_sha256, *input_record_sha256s),
                budget=budget,
                expected_output_schema="specialist-candidate-v1",
            )
            response = self._run_fresh_session(request)
            try:
                candidate, summary = _parse_specialist_candidate(
                    request=request, response=response
                )
            except ContractError as exc:
                raise SpecialistAdvisoryValidationError(
                    role=role,
                    rule_id=_specialist_validation_rule_id(exc),
                    detail=str(exc),
                ) from exc
            packet = _result_packet(
                request=request,
                output_schema="specialist-candidate-v1",
                output_record_sha256=candidate.candidate_sha256,
                public_summary=summary,
                status=candidate.status,
            )
            body = {
                "schema_version": "chemsmart.specialist-session-outcome.v1",
                "request_sha256": request.request_sha256,
                "response_sha256": response.response_sha256,
                "result_packet": packet,
                "candidate": candidate,
            }
            outcomes.append(
                SpecialistSessionOutcomeV1(
                    **body, outcome_sha256=canonical_sha256(body)
                )
            )
        return tuple(outcomes)

    def merge_before_coordinator(
        self, outcomes: Sequence[SpecialistSessionOutcomeV1]
    ) -> CoordinatorMergeReceiptV1:
        """Deterministically merge proposals; conflicts remain unresolved."""

        return merge_specialist_outcomes(outcomes)

    def run_after_coordinator_critic(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        coordinator_session_id: str,
        candidate_id: str,
        candidate_record: Any,
        public_context: Mapping[str, Any],
        source_sha256s: Sequence[str] = (),
        artifact_sha256s: Sequence[str] = (),
        input_record_sha256s: Sequence[str] = (),
        budget: SpecialistBudgetV1,
    ) -> CriticSessionOutcomeV1:
        """Cross-examine a detached candidate in a fresh read-only session."""

        candidate_id = require_identifier(candidate_id, "candidate_id")
        detached_candidate = canonical_data(candidate_record)
        candidate_sha256_before = canonical_sha256(detached_candidate)
        request = self._build_request(
            plan=plan,
            coordinator_session_id=coordinator_session_id,
            role=READ_ONLY_CRITIC,
            public_context={
                "candidate_id": candidate_id,
                "candidate_record": detached_candidate,
                "review_context": canonical_data(public_context),
            },
            source_sha256s=source_sha256s,
            artifact_sha256s=artifact_sha256s,
            input_record_sha256s=(
                plan.plan_sha256,
                candidate_sha256_before,
                *input_record_sha256s,
            ),
            budget=budget,
            expected_output_schema="critic-review-record-v1",
        )
        response = self._run_fresh_session(request)
        review, summary = _parse_critic_review(
            request=request,
            response=response,
            candidate_id=candidate_id,
            candidate_sha256=candidate_sha256_before,
        )
        candidate_sha256_after = canonical_sha256(detached_candidate)
        packet = _result_packet(
            request=request,
            output_schema="critic-review-record-v1",
            output_record_sha256=review.review_sha256,
            public_summary=summary,
            status=review.status,
        )
        body = {
            "schema_version": "chemsmart.critic-session-outcome.v1",
            "request_sha256": request.request_sha256,
            "response_sha256": response.response_sha256,
            "result_packet": packet,
            "review": review,
            "candidate_sha256_before": candidate_sha256_before,
            "candidate_sha256_after": candidate_sha256_after,
        }
        return CriticSessionOutcomeV1(
            **body, outcome_sha256=canonical_sha256(body)
        )

    def _build_request(
        self,
        *,
        plan: ScientificWorkflowPlanV2,
        coordinator_session_id: str,
        role: str,
        public_context: Mapping[str, Any],
        source_sha256s: Sequence[str],
        artifact_sha256s: Sequence[str],
        input_record_sha256s: Sequence[str],
        budget: SpecialistBudgetV1,
        expected_output_schema: str,
    ) -> SpecialistSessionRequestV1:
        role = _normalized_role(role)
        self._serial += 1
        session_id = require_identifier(
            self._session_id_factory(role, self._serial), "session_id"
        )
        if (
            session_id == coordinator_session_id
            or session_id in self._used_session_ids
        ):
            raise ContractError("provider session factory reused a session ID")
        self._used_session_ids.add(session_id)
        suffix = canonical_sha256(
            {"session_id": session_id, "serial": self._serial}
        )[:16]
        surface = build_specialist_tool_surface(
            role=role, base_surface=self._base_tool_surface
        )
        tool_names = tuple(
            sorted(_tool_name(item) for item in surface.tool_definitions)
        )
        context_body = {
            "schema_version": "chemsmart.context-manifest.v1",
            "manifest_id": f"context.{role}.{suffix}",
            "workflow_id": plan.workflow_id,
            "task_spec_sha256": plan.task_spec_sha256,
            "scientific_identity_sha256": plan.scientific_identity_sha256,
            "source_sha256s": tuple(sorted(set(source_sha256s))),
            "artifact_sha256s": tuple(sorted(set(artifact_sha256s))),
            "tool_schema_sha256": surface.tool_schema_sha256,
            "allowed_tools": tool_names,
            "token_budget": budget.token_budget,
            "tool_call_budget": budget.tool_call_budget,
            "wall_time_seconds": budget.wall_time_seconds,
        }
        context = ContextManifestV1(
            **context_body, manifest_sha256=canonical_sha256(context_body)
        )
        if context.manifest_sha256 in self._used_context_sha256s:
            raise ContractError("specialist context manifest was reused")
        self._used_context_sha256s.add(context.manifest_sha256)
        packet_body = {
            "schema_version": "chemsmart.specialist-task-packet.v1",
            "packet_id": f"packet.{role}.{suffix}",
            "workflow_id": plan.workflow_id,
            "role": role,
            "context_manifest_sha256": context.manifest_sha256,
            "input_record_sha256s": tuple(
                sorted(set(input_record_sha256s))
            ),
            "expected_output_schema": expected_output_schema,
            "owner": "coordinator",
            "merge_key": f"{plan.workflow_id}.{role}",
        }
        packet = SpecialistTaskPacketV1(
            **packet_body, packet_sha256=canonical_sha256(packet_body)
        )
        body = {
            "schema_version": "chemsmart.specialist-session-request.v1",
            "session_id": session_id,
            "coordinator_session_id": str(coordinator_session_id),
            "role": role,
            "context_manifest": context,
            "task_packet": packet,
            "system_instruction": _role_instruction(role),
            "public_context_json": canonical_json(public_context),
            "tool_definitions_json": canonical_json(surface.tool_definitions),
        }
        return SpecialistSessionRequestV1(
            **body, request_sha256=canonical_sha256(body)
        )

    def _run_fresh_session(
        self, request: SpecialistSessionRequestV1
    ) -> SpecialistSessionResponseV1:
        session = self._session_factory(request)
        if session.session_id != request.session_id:
            raise ContractError("provider session identity differs from request")
        try:
            response = session.run(request)
        finally:
            session.close()
        if not isinstance(response, SpecialistSessionResponseV1):
            raise ContractError("specialist session returned an untyped response")
        if response.session_id != request.session_id:
            raise ContractError("specialist response came from another session")
        for tool_name in response.tool_calls:
            request.require_tool_allowed(tool_name)
        context = request.context_manifest
        if response.input_tokens + response.output_tokens > context.token_budget:
            raise ContractError("specialist exceeded its token budget")
        if len(response.tool_calls) > context.tool_call_budget:
            raise ContractError("specialist exceeded its tool-call budget")
        if response.wall_time_millis > context.wall_time_seconds * 1000:
            raise ContractError("specialist exceeded its wall-time budget")
        return response


def merge_specialist_outcomes(
    outcomes: Sequence[SpecialistSessionOutcomeV1],
) -> CoordinatorMergeReceiptV1:
    if not outcomes:
        raise ContractError("coordinator merge requires specialist outcomes")
    workflow_ids = {item.candidate.workflow_id for item in outcomes}
    if len(workflow_ids) != 1:
        raise ContractError("cannot merge specialists from different workflows")
    by_field: dict[str, list[SpecialistFieldProposalV1]] = {}
    unresolved: set[str] = set()
    noncomplete: set[str] = set()
    for outcome in outcomes:
        if outcome.result_packet.status != "complete":
            noncomplete.add(outcome.result_packet.result_sha256)
        unresolved.update(outcome.candidate.unresolved_fields)
        if outcome.candidate.status == "failed":
            continue
        for proposal in outcome.candidate.proposals:
            by_field.setdefault(proposal.field_path, []).append(proposal)

    merged: list[MergedFieldProposalV1] = []
    conflicts: list[CoordinatorMergeConflictV1] = []
    for field_path, proposals in sorted(by_field.items()):
        value_groups: dict[str, list[SpecialistFieldProposalV1]] = {}
        for proposal in proposals:
            value_groups.setdefault(proposal.value_json, []).append(proposal)
        if len(value_groups) > 1:
            conflict_body = {
                "schema_version": "chemsmart.coordinator-merge-conflict.v1",
                "field_path": field_path,
                "proposal_sha256s": tuple(
                    sorted(item.proposal_sha256 for item in proposals)
                ),
                "value_sha256s": tuple(
                    sorted(canonical_sha256(value) for value in value_groups)
                ),
            }
            conflicts.append(
                CoordinatorMergeConflictV1(
                    **conflict_body,
                    conflict_sha256=canonical_sha256(conflict_body),
                )
            )
            unresolved.add(field_path)
            continue
        value_json = next(iter(value_groups))
        merged_body = {
            "schema_version": "chemsmart.merged-field-proposal.v1",
            "field_path": field_path,
            "value_json": value_json,
            "source_proposal_sha256s": tuple(
                sorted(item.proposal_sha256 for item in proposals)
            ),
            "source_roles": tuple(sorted({item.role for item in proposals})),
            "evidence_sha256s": tuple(
                sorted(
                    {
                        digest
                        for item in proposals
                        for digest in item.evidence_sha256s
                    }
                )
            ),
        }
        merged.append(
            MergedFieldProposalV1(
                **merged_body, merged_sha256=canonical_sha256(merged_body)
            )
        )

    if conflicts or noncomplete or unresolved:
        status = "needs_clarification"
    elif merged:
        status = "merged"
    else:
        status = "no_candidate"
    body = {
        "schema_version": "chemsmart.coordinator-merge-receipt.v1",
        "workflow_id": next(iter(workflow_ids)),
        "coordinator_owner": "coordinator",
        "source_result_sha256s": tuple(
            sorted(item.result_packet.result_sha256 for item in outcomes)
        ),
        "merged_fields": tuple(merged),
        "conflicts": tuple(conflicts),
        "unresolved_fields": tuple(sorted(unresolved)),
        "noncomplete_result_sha256s": tuple(sorted(noncomplete)),
        "status": status,
    }
    return CoordinatorMergeReceiptV1(
        **body, receipt_sha256=canonical_sha256(body)
    )


def _parse_specialist_candidate(
    *,
    request: SpecialistSessionRequestV1,
    response: SpecialistSessionResponseV1,
) -> tuple[SpecialistCandidateV1, str]:
    output = response.public_output
    _require_exact_output_keys(
        output,
        required={"status", "public_summary", "field_proposals", "unresolved_fields"},
    )
    status = str(output["status"])
    summary = str(output["public_summary"]).strip()
    if not summary:
        raise ContractError("specialist public summary is required")
    raw_proposals = output["field_proposals"]
    if not isinstance(raw_proposals, list):
        raise ContractError("field_proposals must be an array")
    raw_unresolved = output["unresolved_fields"]
    if not isinstance(raw_unresolved, list):
        raise ContractError("unresolved_fields must be an array")
    unresolved_fields = []
    for index, item in enumerate(raw_unresolved):
        if not isinstance(item, str):
            raise ContractError(
                f"unresolved_fields[{index}] must be a lower-case field-path "
                "string, not an object; put explanation in public_summary or "
                "a proposal uncertainty"
            )
        normalized = item.strip().lower()
        try:
            require_identifier(normalized, "unresolved_field")
        except ContractError as exc:
            raise ContractError(
                f"unresolved_fields[{index}] has invalid field path {item!r}"
            ) from exc
        unresolved_fields.append(normalized)
    proposals = []
    for index, raw in enumerate(raw_proposals):
        if not isinstance(raw, Mapping):
            raise ContractError("field proposal must be an object")
        _require_exact_output_keys(
            raw,
            required={"field_path", "value"},
            optional={"evidence_sha256s", "uncertainty"},
        )
        field_path = str(raw["field_path"]).strip().lower()
        proposal_body = {
            "schema_version": "chemsmart.specialist-field-proposal.v1",
            "proposal_id": (
                f"proposal.{canonical_sha256({'packet': request.task_packet.packet_id, 'index': index, 'field': field_path})[:16]}"
            ),
            "role": request.role,
            "field_path": field_path,
            "value_json": canonical_json(raw["value"]),
            "evidence_sha256s": tuple(
                sorted(set(raw.get("evidence_sha256s") or ()))
            ),
            "uncertainty": str(raw.get("uncertainty") or ""),
        }
        proposals.append(
            SpecialistFieldProposalV1(
                **proposal_body,
                proposal_sha256=canonical_sha256(proposal_body),
            )
        )
    proposals = sorted(
        proposals, key=lambda item: (item.field_path, item.proposal_sha256)
    )
    candidate_body = {
        "schema_version": "chemsmart.specialist-candidate.v1",
        "packet_id": request.task_packet.packet_id,
        "workflow_id": request.task_packet.workflow_id,
        "role": request.role,
        "proposals": tuple(proposals),
        "unresolved_fields": tuple(sorted(set(unresolved_fields))),
        "status": status,
    }
    return (
        SpecialistCandidateV1(
            **candidate_body,
            candidate_sha256=canonical_sha256(candidate_body),
        ),
        summary,
    )


def _parse_critic_review(
    *,
    request: SpecialistSessionRequestV1,
    response: SpecialistSessionResponseV1,
    candidate_id: str,
    candidate_sha256: str,
) -> tuple[CriticReviewRecordV1, str]:
    output = response.public_output
    _require_exact_output_keys(
        output, required={"status", "public_summary", "findings"}
    )
    status = str(output["status"])
    summary = str(output["public_summary"]).strip()
    if not summary:
        raise ContractError("critic public summary is required")
    raw_findings = output["findings"]
    if not isinstance(raw_findings, list):
        raise ContractError("critic findings must be an array")
    findings = []
    for index, raw in enumerate(raw_findings):
        if not isinstance(raw, Mapping):
            raise ContractError("critic finding must be an object")
        _require_exact_output_keys(
            raw,
            required={"rule_id", "severity", "expected", "observed"},
            optional={"evidence_sha256s", "disposition", "target_id"},
        )
        if str(raw.get("disposition") or "open") != "open":
            raise ContractError("read-only critic cannot dispose its own finding")
        target_id = str(raw.get("target_id") or candidate_id).strip().lower()
        if target_id != candidate_id:
            raise ContractError("critic finding targets another candidate")
        seed = {
            "packet": request.task_packet.packet_id,
            "index": index,
            "rule_id": raw["rule_id"],
            "target_id": target_id,
        }
        finding_body = {
            "schema_version": "chemsmart.scientific-review-finding.v1",
            "finding_id": f"finding.{canonical_sha256(seed)[:16]}",
            "reviewer_role": READ_ONLY_CRITIC,
            "rule_id": str(raw["rule_id"]).strip().lower(),
            "severity": str(raw["severity"]).strip().lower(),
            "target_id": target_id,
            "evidence_sha256s": tuple(
                sorted(set(raw.get("evidence_sha256s") or ()))
            ),
            "expected": str(raw["expected"]),
            "observed": str(raw["observed"]),
            "disposition": "open",
        }
        findings.append(
            ScientificReviewFindingV1(
                **finding_body,
                finding_sha256=canonical_sha256(finding_body),
            )
        )
    findings = sorted(findings, key=lambda item: item.finding_id)
    review_body = {
        "schema_version": "chemsmart.critic-review-record.v1",
        "candidate_id": candidate_id,
        "candidate_sha256": candidate_sha256,
        "findings": tuple(findings),
        "status": status,
    }
    return (
        CriticReviewRecordV1(
            **review_body, review_sha256=canonical_sha256(review_body)
        ),
        summary,
    )


def _result_packet(
    *,
    request: SpecialistSessionRequestV1,
    output_schema: str,
    output_record_sha256: str,
    public_summary: str,
    status: str,
) -> SpecialistResultPacketV1:
    body = {
        "schema_version": "chemsmart.specialist-result-packet.v1",
        "packet_id": request.task_packet.packet_id,
        "workflow_id": request.task_packet.workflow_id,
        "role": request.role,
        "context_manifest_sha256": request.context_manifest.manifest_sha256,
        "output_schema": output_schema,
        "output_record_sha256": output_record_sha256,
        "public_summary": public_summary,
        "status": status,
    }
    return SpecialistResultPacketV1(
        **body, result_sha256=canonical_sha256(body)
    )


def _role_instruction(role: str) -> str:
    common = (
        "Return exactly one JSON object without Markdown or surrounding prose. "
        "Never provide native input, shell, paths, approval, readiness, "
        "execution, or terminal state. "
        "Never use authority-bearing segments such as approval, approved, "
        "readiness, execution, execute, executable, success, terminal, command, "
        "argv, cwd, path, shell, script, or execution_ready. Distinguish loader observations, "
        "preview conformance, environment observations, and scientific suitability; "
        "none implies another. Host-returned authority values may be discussed in "
        "public_summary but must not be recast as a specialist decision. Do not claim quantitative "
        "accuracy, cost, or density-fitting effects, and do not call an alternative "
        "runnable without typed host evidence. Use only the supplied tools and "
        "budget. Do not claim RI or density fitting unless the exact project "
        "explicitly enables density_fit. Do not state an implementation-specific "
        "functional alias or correlation convention unless the immutable packet "
        "contains a host functional-resolution receipt; otherwise report the "
        "requested literal and leave applied semantics unresolved. Functional "
        "resolution is not environment-observation or scientific-suitability "
        "evidence."
    )
    proposal_common = (
        " Propose typed scientific fields. In proposal field paths and nested "
        "value-object keys, describe evidence only with advisory names such as "
        "loader_observation, preview_conformance, environment_observation, "
        "suitability_assessment, and unresolved_fields. "
        "unresolved_fields must be a JSON array of lower-case field-path strings "
        "only, for example [\"workflow.nodes.hess.geometry_input\"]; never put "
        "objects there. Put explanations in public_summary or proposal uncertainty. "
    )
    if role == READ_ONLY_CRITIC:
        return (
            "You are a fresh read-only computational-chemistry critic. "
            "Cross-examine the detached candidate and return only open findings. "
            "Do not repair, approve, execute, or set readiness. Return exactly "
            "these keys: status, public_summary, findings. status is complete, "
            "blocked, or failed and describes only this review packet, not workflow "
            "readiness. Each finding has exactly rule_id, severity, expected, "
            "observed, optional evidence_sha256s, optional target_id, and "
            "optional disposition=open. severity must be exactly info, warning, "
            "or critical. Omit target_id unless it is exactly the detached "
            "candidate_id; describe a narrower subtarget in observed instead. "
            "Never add proposal, repair, recommendation, or any other finding "
            "key. Before returning, verify that the entire response is one valid "
            "JSON object. " + common
        )
    ownership = {
        SCIENTIFIC_SPECIALIST: "identity., project., or scientific.",
        PYSCF_SPECIALIST: "project. or scientific.",
        DAG_SPECIALIST: "workflow.",
    }[role]
    return (
        f"You are the bounded {role}. Work independently from other workers. "
        "Return exactly these keys: status, public_summary, field_proposals, "
        "and unresolved_fields. status is complete, blocked, or failed and "
        "describes only this advisory packet, not workflow readiness. Each "
        "field_proposal has field_path, value, optional evidence_sha256s, and "
        "optional uncertainty, with no other keys. Consolidate related facts "
        "instead of repeating one proposal per receipt. Before returning, parse "
        "the entire response mentally as one valid JSON object and check that "
        "every array and object is closed exactly once. Field paths for this role start with "
        f"{ownership} "
        + proposal_common
        + common
    )


def _default_session_id(role: str, serial: int) -> str:
    return f"specialist.{role}.{serial}.{uuid4().hex}"


def _canonical_json_record(value: str, name: str) -> Any:
    try:
        parsed = json.loads(value)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ContractError(f"{name} must be canonical JSON") from exc
    if canonical_json(parsed) != value:
        raise ContractError(f"{name} is not canonical JSON")
    return parsed


def _require_exact_output_keys(
    value: Mapping[str, Any],
    *,
    required: set[str],
    optional: set[str] | None = None,
) -> None:
    optional = optional or set()
    observed = set(value)
    missing = required.difference(observed)
    extra = observed.difference(required.union(optional))
    if missing:
        raise ContractError(
            "typed specialist output is missing: " + ", ".join(sorted(missing))
        )
    if extra:
        raise ContractError(
            "typed specialist output has unsupported fields: "
            + ", ".join(sorted(extra))
        )


def _reject_private_reasoning(value: Any) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if str(key).strip().lower() in _PRIVATE_REASONING_KEYS:
                raise ContractError("provider-private reasoning cannot be persisted")
            _reject_private_reasoning(item)
    elif isinstance(value, list):
        for item in value:
            _reject_private_reasoning(item)


def _reject_authority_path(field_path: str) -> None:
    segments = set(field_path.lower().split("."))
    forbidden = segments.intersection(_FORBIDDEN_AUTHORITY_SEGMENTS)
    if forbidden:
        raise ContractError(
            "specialist proposal requests host authority: "
            + ", ".join(sorted(forbidden))
        )


def _reject_authority_value(value: Any) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized = str(key).strip().lower()
            if normalized in _FORBIDDEN_AUTHORITY_SEGMENTS:
                raise ContractError(
                    f"specialist value contains host authority field {normalized!r}"
                )
            _reject_authority_value(item)
        return
    if isinstance(value, list):
        for item in value:
            _reject_authority_value(item)
        return
    if not isinstance(value, str):
        return
    if _ABSOLUTE_PATH_VALUE.search(value) or _WINDOWS_ABSOLUTE_PATH.search(value):
        raise ContractError("specialist value cannot transfer a filesystem path")
    if any(token in value for token in ("\x00", "`", "$(", "&&", "||")):
        raise ContractError("specialist value cannot transfer shell syntax")
    if _CLEAR_EXECUTABLE_PAYLOAD.search(value):
        raise ContractError("specialist value cannot transfer shell syntax")
    if _CLEAR_NATIVE_INPUT_PAYLOAD.search(value):
        raise ContractError("specialist value cannot transfer native input")


def validate_specialist_advisory_value(value: Any) -> None:
    """Apply the canonical path, shell, and native-input advisory boundary."""

    _reject_authority_value(value)


def _specialist_validation_rule_id(exc: ContractError) -> str:
    message = str(exc)
    if "filesystem path" in message:
        return "specialist.advisory.filesystem-authority.v1"
    if "shell syntax" in message:
        return "specialist.advisory.shell-authority.v1"
    if "native input" in message:
        return "specialist.advisory.native-input.v1"
    if "host authority" in message:
        return "specialist.advisory.host-authority.v1"
    return "specialist.advisory.typed-contract.v1"


__all__ = [
    "BoundedSpecialistOrchestratorV1",
    "CoordinatorMergeConflictV1",
    "CoordinatorMergeReceiptV1",
    "CriticReviewRecordV1",
    "CriticSessionOutcomeV1",
    "DAG_SPECIALIST",
    "MergedFieldProposalV1",
    "PYSCF_SPECIALIST",
    "READ_ONLY_CRITIC",
    "SCIENTIFIC_SPECIALIST",
    "SpecialistAdvisoryValidationError",
    "SpecialistBudgetV1",
    "SpecialistCandidateV1",
    "SpecialistFieldProposalV1",
    "SpecialistProviderSessionV1",
    "SpecialistSessionFactoryV1",
    "SpecialistSessionOutcomeV1",
    "SpecialistSessionRequestV1",
    "SpecialistSessionResponseV1",
    "validate_specialist_advisory_value",
    "build_specialist_session_response",
    "build_specialist_tool_surface",
    "merge_specialist_outcomes",
]
