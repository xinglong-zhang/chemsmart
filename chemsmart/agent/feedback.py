"""Versioned provider-visible projections of canonical tool results.

The event stream always records the complete public tool result.  ``full-v1``
passes those bytes to the provider unchanged.  ``causal-v1`` is deliberately
narrower: it exposes only the values registered as necessary for a subsequent
typed action, scientific interpretation, or an honest stop/repair decision.

The causal contract does *not* claim that a projection is semantically
equivalent to the complete result.  It proves only that a versioned causal
action signature, independently recomputed from the canonical result and the
provider projection, is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_json,
    canonical_sha256,
)
from chemsmart.agent.workflows import CausalToolFeedbackV1


FULL_FEEDBACK_V1 = "full-v1"
CAUSAL_FEEDBACK_V1 = "causal-v1"
FEEDBACK_MODES = frozenset({FULL_FEEDBACK_V1, CAUSAL_FEEDBACK_V1})

CAUSAL_ACTION_PROJECTION_V1 = "chemsmart.causal-action-projection.v1"
CAUSAL_ACTION_SIGNATURE_V1 = "chemsmart.causal-action-signature.v1"

# These views are either executable/rendered material, provider-irrelevant
# local paths, or broad inventories.  Their canonical bytes and digest remain
# in the host event record; causal feedback never sends them to a provider.
_FORBIDDEN_CAUSAL_KEYS = frozenset(
    {
        "allowed_engines",
        "allowed_jobtypes",
        "argv",
        "canonical_command",
        "cli_value",
        "covered_engines",
        "covered_jobtypes",
        "engine_job_capabilities",
        "engines",
        "generated_input",
        "jobtypes",
        "native_input",
        "options",
        "path",
        "programs",
        "project_owned_parameters",
        "project_parameter_domains",
        "public_transcript",
        "raw_output",
        "rendered_command",
        "rendered_yaml",
        "script",
        "script_text",
        "stderr",
        "stdout",
        "tool_definitions",
        "yaml_text",
    }
)

# Scalar names with causal meaning across typed ChemSmart contracts.  Digests,
# stable IDs, explicit statuses, counts, and unit-bearing fields are also kept
# by :func:`_is_action_scalar_key`.
_COMMON_ACTION_KEYS = frozenset(
    {
        "actionable_node_ids",
        "affected_field",
        "affected_fields",
        "approval_required",
        "artifact_class",
        "basis",
        "basis_mode",
        "binding_id",
        "charge",
        "command_path",
        "decision",
        "dependencies",
        "diagnostic",
        "dispersion",
        "engine",
        "error_class",
        "evidence_refs",
        "execution_ready",
        "engine_complete",
        "execution_target",
        "exit_status",
        "expected",
        "fake_mode",
        "functional",
        "jobtype",
        "kind",
        "loader_id",
        "message",
        "method",
        "method_family",
        "method_name",
        "multiplicity",
        "no_scratch_mode",
        "observation_method",
        "observed_version",
        "observed",
        "output_id",
        "plan_state",
        "program",
        "project_role",
        "reason",
        "request_id",
        "requested_engine",
        "requested_program",
        "response_method",
        "rule_id",
        "rule_ids",
        "scheduler_test_mode",
        "schema_version",
        "selected_engine",
        "selected_program",
        "severity",
        "sha256",
        "size_bytes",
        "solvent",
        "solvent_model",
        "spin",
        "stage_order",
        "standard_state",
        "state",
        "status",
        "temperature",
        "target_kind",
        "transport_status",
        "unit",
        "units",
        "unresolved_fields",
        "unresolved_node_ids",
        "validation_status",
        "validated",
        "scientifically_validated",
        "value",
    }
)

# Branches whose arbitrary child names are scientific data rather than an
# inventory.  For example, project setting names are program-specific and must
# survive without teaching this module a second settings registry.
_COMMON_WHOLE_BRANCHES = frozenset(
    {
        "findings",
        "sections",
        "settings",
    }
)


@dataclass(frozen=True)
class _ToolProjectionContractV1:
    tool: str
    whole_branches: frozenset[str] = frozenset()
    extra_action_keys: frozenset[str] = frozenset()
    dropped_branches: frozenset[str] = frozenset()


# Causal mode is fail-closed: adding a model-visible tool requires consciously
# registering its action surface here.  Most contracts use the shared typed
# key policy; per-tool branches retain only scientific values which cannot be
# named generically (project settings, explicit findings, and small queries).
_PROJECTION_CONTRACTS = {
    item.tool: item
    for item in (
        _ToolProjectionContractV1(
            "inspect_program_capability", whole_branches=frozenset({"query"})
        ),
        _ToolProjectionContractV1(
            "inspect_program_environment",
        ),
        _ToolProjectionContractV1(
            "assess_program_candidate",
            whole_branches=frozenset({"findings"}),
        ),
        _ToolProjectionContractV1("render_project_yaml"),
        _ToolProjectionContractV1("promote_project_yaml"),
        _ToolProjectionContractV1("bind_scientific_identity"),
        _ToolProjectionContractV1(
            "read_project_yaml", whole_branches=frozenset({"sections"})
        ),
        _ToolProjectionContractV1(
            "validate_project_yaml",
            whole_branches=frozenset(
                {
                    "settings",
                    "findings",
                    "scientific_materializations",
                    "decision_binding",
                }
            ),
        ),
        _ToolProjectionContractV1(
            "plan_command_workflow",
            # ``workflow_context`` is retained whole: it is host-derived, and
            # every field of it -- what a node waits for, which upstream output
            # feeds it, what depends on it -- is exactly what the next typed
            # action needs. Projecting it field-by-field would deliver the
            # readiness verdict without the reason for it.
            whole_branches=frozenset({"findings", "workflow_context"}),
            dropped_branches=frozenset({"nodes"}),
        ),
        _ToolProjectionContractV1("synthesize_command"),
        _ToolProjectionContractV1("repair_command"),
        _ToolProjectionContractV1(
            "preview_command", dropped_branches=frozenset({"artifacts"})
        ),
        _ToolProjectionContractV1(
            "preflight_program_node",
            whole_branches=frozenset({"findings"}),
        ),
        _ToolProjectionContractV1("record_scientific_decision"),
        _ToolProjectionContractV1(
            "inspect_calculation_artifact",
            whole_branches=frozenset({"findings", "properties"}),
        ),
        _ToolProjectionContractV1(
            "extract_result_quantities",
            whole_branches=frozenset({"quantities"}),
        ),
        _ToolProjectionContractV1(
            "derive_thermochemistry",
            whole_branches=frozenset({"quantities", "assumptions"}),
        ),
        _ToolProjectionContractV1(
            "evaluate_quantity_expression",
            whole_branches=frozenset(
                {"outputs", "node_values", "output_dependencies"}
            ),
        ),
        _ToolProjectionContractV1(
            "record_analysis_claims",
            whole_branches=frozenset({"claims"}),
        ),
        _ToolProjectionContractV1(
            "execute_approved_program_node",
            whole_branches=frozenset({"findings", "native_failure"}),
        ),
    )
}


@dataclass(frozen=True)
class ToolFeedbackEquivalenceReceiptV1:
    """Historical receipt reader for already-recorded v1 experiment data.

    New projections never create this receipt because its verdict overstated
    causal projection as complete semantic equivalence.  It remains available
    so downstream readers can recognize historical event payloads.
    """

    schema_version: str
    tool: str
    mode: str
    canonical_result_sha256: str
    provider_feedback_sha256: str
    causal_feedback_sha256: str
    semantic_signature_sha256: str
    omitted_paths: tuple[str, ...]
    canonical_bytes: int
    provider_feedback_bytes: int
    verdict: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.tool-feedback-equivalence.v1":
            raise ContractError("unsupported historical tool-feedback receipt")
        if self.mode not in FEEDBACK_MODES:
            raise ContractError("unsupported tool-feedback projection")
        if self.verdict != "causal_semantics_preserved":
            raise ContractError("historical tool-feedback verdict is invalid")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("historical tool-feedback digest mismatch")


@dataclass(frozen=True)
class CausalActionSignatureV1:
    """The exact subset promised by the causal projection contract."""

    schema_version: str
    projection_contract: str
    tool: str
    outcome_status: str
    action_values_sha256: str
    rule_ids: tuple[str, ...]
    affected_fields: tuple[str, ...]
    signature_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != CAUSAL_ACTION_SIGNATURE_V1:
            raise ContractError("unsupported causal action signature")
        if self.projection_contract != CAUSAL_ACTION_PROJECTION_V1:
            raise ContractError("unsupported causal action projection contract")
        if self.tool not in _PROJECTION_CONTRACTS:
            raise ContractError("causal action signature targets an unregistered tool")
        if self.outcome_status not in {"valid", "invalid", "blocked", "failed"}:
            raise ContractError("unsupported causal action outcome")
        for values in (self.rule_ids, self.affected_fields):
            if values != tuple(sorted(set(values))):
                raise ContractError("causal action identifiers are not canonical")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "signature_sha256"
        }
        if self.signature_sha256 != canonical_sha256(body):
            raise ContractError("causal action signature digest mismatch")


@dataclass(frozen=True)
class ToolFeedbackProjectionReceiptV2:
    """Bind a provider view to a canonical result and a narrow action contract."""

    schema_version: str
    tool: str
    mode: str
    projection_contract: str
    canonical_result_sha256: str
    provider_feedback_sha256: str
    causal_feedback_sha256: str
    causal_action_signature_sha256: str
    omitted_paths: tuple[str, ...]
    canonical_bytes: int
    provider_feedback_bytes: int
    bytes_reduced: int
    reduction_basis_points: int
    verdict: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.tool-feedback-projection-receipt.v2":
            raise ContractError("unsupported tool-feedback projection receipt")
        if self.mode not in FEEDBACK_MODES:
            raise ContractError("unsupported tool-feedback projection")
        expected_contract = (
            "chemsmart.full-tool-result.v1"
            if self.mode == FULL_FEEDBACK_V1
            else CAUSAL_ACTION_PROJECTION_V1
        )
        if self.projection_contract != expected_contract:
            raise ContractError("tool-feedback projection contract mismatch")
        expected_verdict = (
            "full_result_preserved"
            if self.mode == FULL_FEEDBACK_V1
            else "causal_action_signature_preserved"
        )
        if self.verdict != expected_verdict:
            raise ContractError("tool-feedback projection verdict mismatch")
        if self.canonical_bytes <= 0 or self.provider_feedback_bytes <= 0:
            raise ContractError("tool-feedback byte counts must be positive")
        if self.bytes_reduced != self.canonical_bytes - self.provider_feedback_bytes:
            raise ContractError("tool-feedback byte reduction is inconsistent")
        expected_bps = round(self.bytes_reduced * 10_000 / self.canonical_bytes)
        if self.reduction_basis_points != expected_bps:
            raise ContractError("tool-feedback reduction ratio is inconsistent")
        if self.omitted_paths != tuple(sorted(set(self.omitted_paths))):
            raise ContractError("omitted tool-feedback paths are not canonical")
        body = {
            key: value
            for key, value in self.__dict__.items()
            if key != "receipt_sha256"
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("tool-feedback projection receipt digest mismatch")


@dataclass(frozen=True)
class ProjectedToolFeedbackV1:
    """One provider-visible result and its host-owned projection evidence."""

    content: dict[str, Any]
    receipt: ToolFeedbackProjectionReceiptV2
    action_signature: CausalActionSignatureV1 | None = None


def project_tool_feedback(
    *, tool: str, result: Mapping[str, Any], mode: str
) -> ProjectedToolFeedbackV1:
    """Return full bytes or a registered causal-action projection.

    Full feedback is intentionally byte-for-byte compatible with ``full-v1``.
    Causal feedback rejects unregistered tools before any information can be
    silently discarded.
    """

    normalized_mode = str(mode).strip().lower()
    if normalized_mode not in FEEDBACK_MODES:
        raise ContractError("unsupported tool-feedback projection")
    normalized_tool = str(tool).strip().lower()
    canonical = canonical_data(dict(result))
    if not isinstance(canonical, dict):  # pragma: no cover - defensive
        raise ContractError("tool result must be an object")
    observed_tool = str(canonical.get("tool") or "").strip().lower()
    if (
        normalized_mode == CAUSAL_FEEDBACK_V1
        and observed_tool
        and observed_tool != normalized_tool
    ):
        raise ContractError("causal feedback tool identity mismatch")

    canonical_bytes = len(canonical_json(canonical).encode("utf-8"))
    if normalized_mode == FULL_FEEDBACK_V1:
        provider_feedback = canonical
        signature = None
        omitted: list[str] = []
        projection_contract = "chemsmart.full-tool-result.v1"
        causal_feedback_sha256 = canonical_sha256(canonical)
        causal_action_signature_sha256 = canonical_sha256(
            {
                "schema_version": "chemsmart.full-tool-result-signature.v1",
                "tool": normalized_tool,
                "canonical_result_sha256": canonical_sha256(canonical),
            }
        )
        verdict = "full_result_preserved"
    else:
        contract = _PROJECTION_CONTRACTS.get(normalized_tool)
        if contract is None:
            raise ContractError(
                "causal feedback requires a registered tool projection"
            )
        omitted = []
        payload = _project_causal_payload(
            canonical,
            contract=contract,
            omitted=omitted,
        )
        canonical_signature = _build_action_signature(
            tool=normalized_tool,
            canonical=canonical,
            action_values=payload,
        )
        causal_record = _causal_record(
            tool=normalized_tool, canonical=canonical
        )
        # One non-redundant provider contract.  The former envelope duplicated
        # status/rules/full-result identity in both ``feedback`` and ``result``
        # and could exceed the canonical response for ordinary tools.
        provider_feedback = {
            "schema_version": CAUSAL_ACTION_PROJECTION_V1,
            "tool": normalized_tool,
            "status": canonical_signature.outcome_status,
            "transport_status": str(canonical.get("status") or "unknown"),
            "result": payload,
        }
        observed_signature = _action_signature_from_provider_projection(
            provider_feedback
        )
        if observed_signature != canonical_signature:
            raise ContractError(
                "provider feedback changed the causal action signature"
            )
        signature = canonical_signature
        projection_contract = CAUSAL_ACTION_PROJECTION_V1
        causal_feedback_sha256 = causal_record.feedback_sha256
        causal_action_signature_sha256 = canonical_signature.signature_sha256
        verdict = "causal_action_signature_preserved"

    provider_bytes = len(canonical_json(provider_feedback).encode("utf-8"))
    bytes_reduced = canonical_bytes - provider_bytes
    body = {
        "schema_version": "chemsmart.tool-feedback-projection-receipt.v2",
        "tool": normalized_tool,
        "mode": normalized_mode,
        "projection_contract": projection_contract,
        "canonical_result_sha256": canonical_sha256(canonical),
        "provider_feedback_sha256": canonical_sha256(provider_feedback),
        "causal_feedback_sha256": causal_feedback_sha256,
        "causal_action_signature_sha256": causal_action_signature_sha256,
        "omitted_paths": tuple(sorted(set(omitted))),
        "canonical_bytes": canonical_bytes,
        "provider_feedback_bytes": provider_bytes,
        "bytes_reduced": bytes_reduced,
        "reduction_basis_points": round(bytes_reduced * 10_000 / canonical_bytes),
        "verdict": verdict,
    }
    receipt = ToolFeedbackProjectionReceiptV2(
        **body, receipt_sha256=canonical_sha256(body)
    )
    return ProjectedToolFeedbackV1(
        content=provider_feedback,
        receipt=receipt,
        action_signature=signature,
    )


def _project_causal_payload(
    canonical: Mapping[str, Any],
    *,
    contract: _ToolProjectionContractV1,
    omitted: list[str],
) -> dict[str, Any]:
    """Project only the canonical tool payload or its public error envelope."""

    result = canonical.get("result")
    if isinstance(result, Mapping):
        projected = _project_mapping(
            result,
            path="$.result",
            contract=contract,
            omitted=omitted,
            preserve_all=False,
        )
    else:
        projected = {}
        if result is not None:
            omitted.append("$.result")

    # Transport-level error fields live outside ``result``.  They remain
    # available for repair/clarification but do not make the tool authoritative.
    for key in ("error_class", "message", "rule_ids", "affected_fields"):
        if key in canonical:
            projected[key] = canonical_data(canonical[key])
    return canonical_data(projected)


def _project_mapping(
    value: Mapping[str, Any],
    *,
    path: str,
    contract: _ToolProjectionContractV1,
    omitted: list[str],
    preserve_all: bool,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    whole_branches = _COMMON_WHOLE_BRANCHES | contract.whole_branches
    for raw_key, item in sorted(value.items(), key=lambda pair: str(pair[0])):
        key = str(raw_key)
        normalized_key = key.lower()
        child_path = f"{path}.{key}"
        # A forbidden key is unsafe to send at any depth, and stays so.  A
        # dropped branch only suppresses a redundant echo, and it does not
        # apply inside a branch the contract consciously retains whole:
        # otherwise a common name such as ``nodes`` would strip an unrelated
        # nested structure the model has never seen.
        if normalized_key in _FORBIDDEN_CAUSAL_KEYS or (
            not preserve_all and normalized_key in contract.dropped_branches
        ):
            omitted.append(child_path)
            continue
        keep_whole = preserve_all or normalized_key in whole_branches
        if isinstance(item, Mapping):
            projected = _project_mapping(
                item,
                path=child_path,
                contract=contract,
                omitted=omitted,
                preserve_all=keep_whole,
            )
            if projected:
                result[key] = projected
            elif item:
                omitted.append(child_path)
            elif _is_action_scalar_key(normalized_key, contract):
                result[key] = {}
            continue
        if isinstance(item, (tuple, list)):
            if keep_whole or _is_action_scalar_key(normalized_key, contract):
                result[key] = _sanitize_retained_sequence(
                    item,
                    path=child_path,
                    contract=contract,
                    omitted=omitted,
                    preserve_all=keep_whole,
                )
            else:
                projected_rows = _project_sequence(
                    item,
                    path=child_path,
                    contract=contract,
                    omitted=omitted,
                )
                if projected_rows:
                    result[key] = projected_rows
                elif item:
                    omitted.append(child_path)
            continue
        if keep_whole or _is_action_scalar_key(normalized_key, contract):
            result[key] = canonical_data(item)
        else:
            omitted.append(child_path)
    return result


def _project_sequence(
    value: tuple[Any, ...] | list[Any],
    *,
    path: str,
    contract: _ToolProjectionContractV1,
    omitted: list[str],
) -> list[Any]:
    result: list[Any] = []
    for index, item in enumerate(value):
        child_path = f"{path}[{index}]"
        if isinstance(item, Mapping):
            projected = _project_mapping(
                item,
                path=child_path,
                contract=contract,
                omitted=omitted,
                preserve_all=False,
            )
            if projected:
                result.append(projected)
            elif item:
                omitted.append(child_path)
        elif isinstance(item, (tuple, list)):
            nested = _project_sequence(
                item,
                path=child_path,
                contract=contract,
                omitted=omitted,
            )
            if nested:
                result.append(nested)
            elif item:
                omitted.append(child_path)
        else:
            omitted.append(child_path)
    return result


def _sanitize_retained_sequence(
    value: tuple[Any, ...] | list[Any],
    *,
    path: str,
    contract: _ToolProjectionContractV1,
    omitted: list[str],
    preserve_all: bool,
) -> list[Any]:
    result: list[Any] = []
    for index, item in enumerate(value):
        child_path = f"{path}[{index}]"
        if isinstance(item, Mapping):
            result.append(
                _project_mapping(
                    item,
                    path=child_path,
                    contract=contract,
                    omitted=omitted,
                    preserve_all=preserve_all,
                )
            )
        elif isinstance(item, (tuple, list)):
            result.append(
                _sanitize_retained_sequence(
                    item,
                    path=child_path,
                    contract=contract,
                    omitted=omitted,
                    preserve_all=preserve_all,
                )
            )
        else:
            result.append(canonical_data(item))
    return result


def _is_action_scalar_key(
    key: str, contract: _ToolProjectionContractV1
) -> bool:
    if key in _COMMON_ACTION_KEYS or key in contract.extra_action_keys:
        return True
    return key.endswith(
        (
            "_sha256",
            "_digest",
            "_id",
            "_ids",
            "_status",
            "_count",
            "_unit",
            "_units",
            "_cm1",
            "_kelvin",
            "_angstrom",
            "_hartree",
        )
    )


def _build_action_signature(
    *, tool: str, canonical: Mapping[str, Any], action_values: Mapping[str, Any]
) -> CausalActionSignatureV1:
    body = {
        "schema_version": CAUSAL_ACTION_SIGNATURE_V1,
        "projection_contract": CAUSAL_ACTION_PROJECTION_V1,
        "tool": tool,
        "outcome_status": _causal_outcome_status(canonical),
        "action_values_sha256": canonical_sha256(action_values),
        "rule_ids": tuple(sorted(set(_collect_identifiers(canonical, "rule_id")))),
        "affected_fields": tuple(
            sorted(set(_collect_identifiers(canonical, "affected_field")))
        ),
    }
    return CausalActionSignatureV1(
        **body, signature_sha256=canonical_sha256(body)
    )


def _action_signature_from_provider_projection(
    provider_feedback: Mapping[str, Any],
) -> CausalActionSignatureV1:
    action_values = provider_feedback.get("result")
    if not isinstance(action_values, Mapping):
        raise ContractError("causal provider feedback lacks action values")
    body = {
        "schema_version": CAUSAL_ACTION_SIGNATURE_V1,
        "projection_contract": str(provider_feedback.get("schema_version") or ""),
        "tool": str(provider_feedback.get("tool") or ""),
        "outcome_status": str(provider_feedback.get("status") or ""),
        "action_values_sha256": canonical_sha256(action_values),
        "rule_ids": tuple(
            sorted(set(_collect_identifiers(action_values, "rule_id")))
        ),
        "affected_fields": tuple(
            sorted(
                set(_collect_identifiers(action_values, "affected_field"))
            )
        ),
    }
    return CausalActionSignatureV1(
        **body, signature_sha256=canonical_sha256(body)
    )


def _causal_record(
    *, tool: str, canonical: Mapping[str, Any]
) -> CausalToolFeedbackV1:
    status = _causal_outcome_status(canonical)
    rule_ids = tuple(sorted(set(_collect_identifiers(canonical, "rule_id"))))
    affected_fields = tuple(
        sorted(set(_collect_identifiers(canonical, "affected_field")))
    )
    body = {
        "schema_version": "chemsmart.causal-tool-feedback.v1",
        "tool_name": tool,
        "full_receipt_sha256": canonical_sha256(canonical),
        "status": status,
        "rule_ids": rule_ids,
        "affected_fields": affected_fields,
        "public_summary": (
            f"Host tool {tool} produced a {status} action outcome; "
            "the complete public result remains authoritative in the host "
            "event record."
        ),
    }
    return CausalToolFeedbackV1(
        **body, feedback_sha256=canonical_sha256(body)
    )


def _causal_outcome_status(canonical: Mapping[str, Any]) -> str:
    """Return the most severe typed outcome, not transport success."""

    mapped = {
        "ok": "valid",
        "supported": "valid",
        "available": "valid",
        "candidate_rendered": "valid",
        "materialized": "valid",
        "planned": "valid",
        "previewed": "valid",
        "recorded": "valid",
        "resolved": "valid",
        "compiled": "valid",
        "validated": "valid",
        "preview_only": "valid",
        "valid": "valid",
        "invalid": "invalid",
        "rejected": "invalid",
        "loader_unavailable": "blocked",
        "missing": "blocked",
        "not_declared": "blocked",
        "not_queried": "blocked",
        "reference_only": "blocked",
        "disabled": "blocked",
        "unknown_program": "blocked",
        "unsupported_jobtype": "blocked",
        "unsupported_engine": "blocked",
        "unsupported_engine_job_combination": "blocked",
        "cli_schema_mismatch": "blocked",
        "blocked": "blocked",
        "failed": "failed",
    }
    rank = {"valid": 0, "invalid": 1, "blocked": 2, "failed": 3}
    observed = [
        mapped[item]
        for item in _collect_statuses(canonical)
        if item in mapped
    ]
    # A completed wrapper or engine is not a scientific success when its
    # explicit result validator is red.  This check is intentionally limited
    # to validation flags; capability booleans such as execution_supported are
    # observations, not a tool failure.
    if _contains_explicit_false(
        canonical, {"validated", "scientifically_validated"}
    ):
        observed.append("invalid")
    if not observed:
        return "failed"
    return max(observed, key=rank.__getitem__)


def _contains_explicit_false(value: Any, targets: set[str]) -> bool:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if str(key).lower() in targets and item is False:
                return True
            if _contains_explicit_false(item, targets):
                return True
    elif isinstance(value, (tuple, list)):
        return any(_contains_explicit_false(item, targets) for item in value)
    return False


def _collect_statuses(value: Any) -> list[str]:
    result: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized_key = str(key).lower()
            if normalized_key in {
                "status",
                "state",
                "plan_state",
                "execution_state",
                "validation_status",
                "program_validation_status",
            } and isinstance(item, str):
                result.append(item.strip().lower())
            result.extend(_collect_statuses(item))
    elif isinstance(value, (tuple, list)):
        for item in value:
            result.extend(_collect_statuses(item))
    return result


def _collect_identifiers(value: Any, target: str) -> list[str]:
    result: list[str] = []
    if isinstance(value, Mapping):
        for key, item in value.items():
            normalized_key = str(key).lower()
            if normalized_key in {target, target + "s"}:
                rows = item if isinstance(item, (tuple, list)) else (item,)
                for row in rows:
                    normalized = str(row or "").strip().lower()
                    if normalized and all(
                        char.isalnum() or char in "_.-" for char in normalized
                    ) and normalized[0].isalpha():
                        result.append(normalized)
            result.extend(_collect_identifiers(item, target))
    elif isinstance(value, (tuple, list)):
        for item in value:
            result.extend(_collect_identifiers(item, target))
    return result


__all__ = [
    "CAUSAL_ACTION_PROJECTION_V1",
    "CAUSAL_ACTION_SIGNATURE_V1",
    "CAUSAL_FEEDBACK_V1",
    "CausalActionSignatureV1",
    "FEEDBACK_MODES",
    "FULL_FEEDBACK_V1",
    "ProjectedToolFeedbackV1",
    "ToolFeedbackEquivalenceReceiptV1",
    "ToolFeedbackProjectionReceiptV2",
    "project_tool_feedback",
]
