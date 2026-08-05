"""Approved model-visible command-compiled tool surface."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.skills import skills_enabled
from chemsmart.agent.workflows import (
    AGGREGATE_NODE_PROGRAM,
    AGGREGATE_NODE_STAGE,
    WORKFLOW_NODE_KINDS,
)
from chemsmart.agent.capabilities import (
    ProgramCapabilityRegistryV1,
    load_program_capabilities,
)


@dataclass(frozen=True)
class AgentToolSurfaceV1:
    schema_version: str
    profile: str
    tool_definitions: tuple[dict, ...]
    tool_schema_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.agent-tool-surface.v1":
            raise ContractError("unsupported agent tool surface schema")
        if self.tool_schema_sha256 != canonical_sha256(self.tool_definitions):
            raise ContractError("agent tool schema digest mismatch")


def build_command_compiled_tool_surface(
    registry: ProgramCapabilityRegistryV1 | None = None,
) -> AgentToolSurfaceV1:
    registry = registry or load_program_capabilities()
    programs = [item.program for item in registry.programs]
    program = {"type": "string", "enum": programs}
    structured_result_program = {"type": "string", "enum": ["pyscf"]}
    digest = {"type": "string", "pattern": "^[0-9a-f]{64}$"}
    tools = (
        _tool(
            "inspect_program_capability",
            "Join core capability, active support overlay, and live Click schema.",
            {"program": program, "jobtype": _string(), "engine": _string()},
            ("program", "jobtype", "engine"),
        ),
        _tool(
            "inspect_program_environment",
            "Inspect a host-provided target-interpreter or executable receipt.",
            {"capability_receipt_sha256": digest},
            ("capability_receipt_sha256",),
        ),
        _tool(
            "assess_program_candidate",
            "Assess the closed, evidence-bound program substitution matrix.",
            {
                "request_id": _string(),
                "requested_program": program,
                "selected_program": program,
                "requested_engine": _string(),
                "selected_engine": _string(),
                "job_families": {"type": "array", "items": _string()},
                "method_family": _string(),
                "method_name": _string(),
                "basis_mode": {
                    "type": "string",
                    "enum": ["uniform", "mixed", "ecp", "mixed_ecp"],
                },
                "constraint_kinds": {"type": "array", "items": _string()},
                "requires_post_hf": {"type": "boolean"},
                "requires_double_hybrid": {"type": "boolean"},
                "functional_equivalence_receipt_sha256": digest,
                "source_claim_sha256s": {"type": "array", "items": digest},
                "capability_receipt_sha256": digest,
            },
            (
                "request_id",
                "requested_program",
                "selected_program",
                "requested_engine",
                "selected_engine",
                "job_families",
                "method_family",
                "method_name",
                "basis_mode",
                "constraint_kinds",
                "requires_post_hf",
                "requires_double_hybrid",
                "source_claim_sha256s",
                "capability_receipt_sha256",
            ),
        ),
        _tool(
            "render_project_yaml",
            (
                "Render a typed project candidate; this does not validate or "
                "write it. Sections must be an object whose stage names map "
                "to setting-name/value objects, never a list."
            ),
            {
                "program": program,
                "sections": {
                    "type": "object",
                    "additionalProperties": {"type": "object"},
                },
            },
            ("program", "sections"),
        ),
        _tool(
            "promote_project_yaml",
            "Promote one rendered candidate into the host-owned task workspace.",
            {
                "render_receipt_sha256": digest,
                "artifact_id": _string(),
            },
            ("render_receipt_sha256", "artifact_id"),
        ),
        _tool(
            "bind_scientific_identity",
            (
                "Bind an explicit charge and multiplicity to the exact host "
                "geometry artifact. input_artifact_id must identify a "
                "geometry_xyz artifact, never a project YAML or result."
            ),
            {
                "input_artifact_id": _string(),
                "task_spec_sha256": digest,
                "charge": {"type": "integer"},
                "multiplicity": {"type": "integer", "minimum": 1},
            },
            (
                "input_artifact_id",
                "task_spec_sha256",
                "charge",
                "multiplicity",
            ),
        ),
        _tool(
            "read_project_yaml",
            "Read an already host-bound project artifact by stable ID.",
            {"program": program, "project_artifact_id": _string()},
            ("program", "project_artifact_id"),
        ),
        _tool(
            "validate_project_yaml",
            "Validate a bound project through the checked-out settings loader.",
            {
                "project_artifact_id": _string(),
                "capability_receipt_sha256": digest,
            },
            ("project_artifact_id", "capability_receipt_sha256"),
        ),
        _tool(
            "consult_domain_skill",
            (
                "Read an advisory domain-knowledge skill listed in the system "
                "prompt. Returns conventions, definitions, and established "
                "facts. It is knowledge only: it never establishes readiness, "
                "approval, terminal state, or an accuracy claim, and never "
                "substitutes for a typed host receipt."
            ),
            {"skill_id": _public_identifier()},
            ("skill_id",),
        ),
        _tool(
            "plan_command_workflow",
            (
                "Build a typed command DAG after binding scientific identity "
                "to every initial geometry. Every node needs at least one "
                "expected output. Future producer outputs remain unresolved. "
                "A null scientific_workflow_plan means the binding must be "
                "repaired and this tool called again."
            ),
            {
                "workflow_id": _string(),
                "task_spec_id": _string(),
                "nodes": {
                    "type": "array",
                    "items": _workflow_node_schema(),
                },
            },
            ("workflow_id", "task_spec_id", "nodes"),
        ),
        _tool(
            "synthesize_command",
            "Compile scientific intent to canonical argv through live Click.",
            {
                "node_id": _string(),
                "program": program,
                "jobtype": _string(),
                "project_artifact_id": _string(),
                "input_artifact_id": _string(),
                "scientific_identity_sha256": digest,
                "charge": {"type": "integer"},
                "multiplicity": {"type": "integer", "minimum": 1},
                "capability_receipt_sha256": digest,
                "engine_binding_sha256": digest,
                "project_validation_receipt_sha256": digest,
            },
            (
                "node_id",
                "program",
                "jobtype",
                "project_artifact_id",
                "input_artifact_id",
                "scientific_identity_sha256",
                "charge",
                "multiplicity",
                "capability_receipt_sha256",
                "engine_binding_sha256",
            ),
        ),
        _tool(
            "repair_command",
            "Apply one constrained counterexample patch without changing bindings.",
            {
                "invocation_sha256": digest,
                "counterexample_id": _string(),
            },
            ("invocation_sha256", "counterexample_id"),
        ),
        _tool(
            "preview_command",
            "Execute exact compiled argv in isolated fake/test mode and hash outputs.",
            {"invocation_sha256": digest},
            ("invocation_sha256",),
        ),
        _tool(
            "preflight_program_node",
            "Cross-check all receipt ancestry and effective scientific settings.",
            {
                "node_id": _string(),
                "capability_receipt_sha256": digest,
                "program_binding_sha256": digest,
                "engine_binding_sha256": digest,
                "geometry_artifact_sha256": digest,
                "scientific_identity_sha256": digest,
                "charge": {"type": "integer"},
                "multiplicity": {"type": "integer", "minimum": 1},
                "project_validation_receipt_sha256": digest,
                "invocation_sha256": digest,
                "command_inspection_receipt_sha256": digest,
                "safe_preview_receipt_sha256": digest,
            },
            (
                "node_id",
                "capability_receipt_sha256",
                "program_binding_sha256",
                "engine_binding_sha256",
                "geometry_artifact_sha256",
                "scientific_identity_sha256",
                "charge",
                "multiplicity",
                "invocation_sha256",
                "command_inspection_receipt_sha256",
            ),
        ),
        _tool(
            "record_scientific_decision",
            "Record public chemical rationale, alternatives, and uncertainty.",
            {
                "decision_id": _string(),
                "task_spec_sha256": digest,
                "assumptions": {"type": "array", "items": _string()},
                "method_rationale": _string(),
                "alternatives": {"type": "array", "items": _string()},
                "uncertainties": {"type": "array", "items": _string()},
                "diagnostics": {"type": "array", "items": _string()},
                "stage_order": {"type": "array", "items": _string()},
                "evidence_refs": {"type": "array", "items": _string()},
                "postprocessing_receipt_sha256s": {
                    "type": "array",
                    "items": digest,
                    "maxItems": 64,
                    "description": (
                        "Exact receipt_sha256 values returned by quantity "
                        "extraction, thermochemistry, or quantity-expression "
                        "tools. The host validates and canonicalizes them; do "
                        "not embed receipt IDs in free-form evidence strings."
                    ),
                },
            },
            (
                "decision_id",
                "task_spec_sha256",
                "assumptions",
                "method_rationale",
                "alternatives",
                "uncertainties",
                "diagnostics",
                "stage_order",
                "evidence_refs",
            ),
        ),
        _tool(
            "inspect_calculation_artifact",
            "Run a deterministic verifier over a host-bound result and run receipt.",
            {
                "program": program,
                "artifact_id": _string(),
                "project_artifact_id": _string(),
                "settings_id": _string(),
                "run_receipt_id": _string(),
            },
            (
                "program",
                "artifact_id",
                "project_artifact_id",
                "settings_id",
                "run_receipt_id",
            ),
        ),
        _tool(
            "extract_result_quantities",
            (
                "Parse selected numerical or scientific fields from a trusted "
                "host-bound result artifact. The model supplies semantic selectors, "
                "never a file path."
            ),
            {
                "program": structured_result_program,
                "artifact_id": _string(),
                "selectors": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 32,
                    "items": {
                        "type": "object",
                        "properties": {
                            "quantity_id": _public_identifier(),
                            "selector": {
                                "type": "string",
                                "enum": [
                                    "energy",
                                    "energies",
                                    "positions",
                                    "symbols",
                                    "vibrational_frequencies",
                                    "homo",
                                    "lumo",
                                    "gap",
                                    "charge",
                                    "multiplicity",
                                    "method",
                                    "basis",
                                ],
                            },
                        },
                        "required": ["quantity_id", "selector"],
                        "additionalProperties": False,
                    },
                },
            },
            ("program", "artifact_id", "selectors"),
        ),
        _tool(
            "derive_thermochemistry",
            (
                "Derive RRHO thermochemistry from a trusted Hessian result at "
                "an explicit temperature and pressure using ChemSmart's common engine."
            ),
            {
                "program": structured_result_program,
                "artifact_id": _string(),
                "temperature_k": {"type": "number", "exclusiveMinimum": 0},
                "pressure_atm": {"type": "number", "exclusiveMinimum": 0},
            },
            ("program", "artifact_id", "temperature_k", "pressure_atm"),
        ),
        _tool(
            "evaluate_quantity_expression",
            (
                "Evaluate a bounded dimension-aware expression DAG over prior "
                "extraction or thermochemistry receipts, or over typed literal "
                "nodes when inputs is empty; Python and formula strings are not "
                "accepted. Expression inputs use local input_id aliases. "
                "For converted outputs, source_value/source_unit are the requested "
                "display pair; value/unit remain the canonical arithmetic pair."
            ),
            {
                "expression_id": _public_identifier(),
                "inputs": {
                    "type": "array",
                    "minItems": 0,
                    "maxItems": 64,
                    "items": {
                        "type": "object",
                        "properties": {
                            "input_id": _public_identifier(),
                            "semantic_role": _public_identifier(),
                            "receipt_sha256": digest,
                            "quantity_id": _public_identifier(),
                        },
                        "required": [
                            "input_id",
                            "receipt_sha256",
                            "quantity_id",
                        ],
                        "additionalProperties": False,
                    },
                },
                "nodes": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 128,
                    "items": _quantity_expression_node_schema(),
                },
                "output_node_ids": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 128,
                    "items": _public_identifier(),
                },
            },
            ("expression_id", "inputs", "nodes", "output_node_ids"),
        ),
        _tool(
            "record_analysis_claims",
            (
                "Bind reportable numerical claims to exact typed receipt "
                "quantities. Supply identifiers and display units only; the "
                "host copies and converts the values."
            ),
            {
                "task_spec_sha256": digest,
                "claims": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 64,
                    "items": {
                        "type": "object",
                        "properties": {
                            "claim_id": _public_identifier(),
                            "receipt_sha256": digest,
                            "quantity_id": _public_identifier(),
                            "display_unit": _string(),
                        },
                        "required": [
                            "claim_id",
                            "receipt_sha256",
                            "quantity_id",
                            "display_unit",
                        ],
                        "additionalProperties": False,
                    },
                },
            },
            ("task_spec_sha256", "claims"),
        ),
    )
    if not skills_enabled():
        # Skills off restores the exact pre-skill tool-schema digest, so a
        # skills-off arm reproduces a recorded baseline byte for byte.
        tools = tuple(
            item
            for item in tools
            if item["function"]["name"] != "consult_domain_skill"
        )
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="command_compiled_preview",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def build_approved_execution_tool_surface(
    registry: ProgramCapabilityRegistryV1 | None = None,
) -> AgentToolSurfaceV1:
    """Command-compiled surface plus one host-resolved execution action.

    The model never supplies argv, paths, resources, approval material, or
    dependency artifacts.  It can only request execution of a previously
    compiled and approved node by its stable identifier.
    """

    full = build_command_compiled_tool_surface(registry)
    tools = full.tool_definitions + (
        _tool(
            "execute_approved_program_node",
            "Execute one host-compiled node only when its workflow approval and dependencies are green.",
            {"node_id": _string()},
            ("node_id",),
        ),
    )
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="command_compiled_approved_execution",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def build_single_agent_baseline_tool_surface(
    registry: ProgramCapabilityRegistryV1 | None = None,
) -> AgentToolSurfaceV1:
    """H0 surface with host-prebound capability/program/environment state."""

    full = build_command_compiled_tool_surface(registry)
    excluded = {
        "inspect_program_capability",
        "inspect_program_environment",
        "assess_program_candidate",
    }
    tools = tuple(
        item
        for item in full.tool_definitions
        if item["function"]["name"] not in excluded
    )
    return AgentToolSurfaceV1(
        schema_version="chemsmart.agent-tool-surface.v1",
        profile="single_agent_typed_baseline",
        tool_definitions=tools,
        tool_schema_sha256=canonical_sha256(tools),
    )


def _string() -> dict:
    return {"type": "string"}


def _public_identifier() -> dict:
    return {
        "type": "string",
        "pattern": "^[a-z][a-z0-9_.-]*$",
        "description": (
            "Lower-case public identifier; use dots or dashes instead of "
            "spaces, parentheses, or placeholder syntax."
        ),
    }


def _workflow_node_schema() -> dict:
    return {
        "type": "object",
        "properties": {
            "node_id": _string(),
            "program": _string(),
            "jobtype": _string(),
            "node_kind": {
                "type": "string",
                "enum": list(WORKFLOW_NODE_KINDS),
                "description": (
                    "Omit or use 'program_call' to invoke a program. Use "
                    "'aggregate' for a stage that combines finished results "
                    "into the requested number: set program to "
                    f"'{AGGREGATE_NODE_PROGRAM}' and jobtype to "
                    f"'{AGGREGATE_NODE_STAGE}'. An observable that needs such "
                    "a stage has no other producer, so declare it rather than "
                    "dropping it; the arithmetic itself is supplied later "
                    "as a quantity expression."
                ),
            },
            "project_role": _string(),
            "dependencies": {"type": "array", "items": _string()},
            "inputs": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "binding_id": _string(),
                        "artifact_id": _string(),
                        "artifact_class": _string(),
                        "producer_node_id": _string(),
                        "producer_output_id": _string(),
                    },
                    "required": [
                        "binding_id",
                        "artifact_class",
                        "producer_node_id",
                        "producer_output_id",
                    ],
                    "additionalProperties": False,
                },
            },
            "expected_outputs": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "output_id": _string(),
                        "artifact_class": _string(),
                    },
                    "required": ["output_id", "artifact_class"],
                    "additionalProperties": False,
                },
            },
            "unresolved_fields": {
                "type": "array",
                "items": _public_identifier(),
            },
        },
        "required": [
            "node_id",
            "program",
            "jobtype",
            "project_role",
            "dependencies",
            "inputs",
            "expected_outputs",
            "unresolved_fields",
        ],
        "additionalProperties": False,
    }


def _quantity_expression_node_schema() -> dict:
    return {
        "type": "object",
        "properties": {
            "node_id": _public_identifier(),
            "operation": {
                "type": "string",
                "enum": [
                    "ref",
                    "literal",
                    "add",
                    "subtract",
                    "multiply",
                    "divide",
                    "scale",
                    "abs",
                    "sqrt",
                    "power",
                    "exp",
                    "log",
                    "sum",
                    "mean",
                    "min",
                    "max",
                    "distance",
                    "angle",
                    "convert",
                    "linear_fit_slope",
                    "linear_fit_intercept",
                ],
            },
            "input_ids": {
                "type": "array",
                "items": _public_identifier(),
                "description": (
                    "Inputs for arithmetic, bounded transforms, reductions, "
                    "linear fits, distance, angle, and conversion. For ref, "
                    "prefer reference instead."
                ),
            },
            "reference": {
                "type": "string",
                "description": (
                    "For ref, identify one expression input alias or earlier node; "
                    "omit input_ids. Other operations use input_ids and omit reference."
                ),
            },
            "indices": {
                "type": "array",
                "items": {"type": "integer", "minimum": 0},
                "description": (
                    "For ref, select nested zero-based indices. Create one indexed "
                    "ref node per coordinate vector before distance or angle."
                ),
            },
            "literal_value": {
                "oneOf": [
                    {"type": "number"},
                    {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 1,
                        "maxItems": 64,
                    },
                ],
                "description": (
                    "Finite scalar/vector for literal, or finite scalar exponent "
                    "for power."
                ),
            },
            "literal_unit": _string(),
            "scale_factor": {"type": "number"},
            "target_unit": _string(),
        },
        "required": ["node_id", "operation"],
        "additionalProperties": False,
    }


def _tool(
    name: str,
    description: str,
    properties: dict,
    required: tuple[str, ...],
) -> dict:
    return {
        "type": "function",
        "function": {
            "name": name,
            "description": description,
            "parameters": {
                "type": "object",
                "properties": properties,
                "required": list(required),
                "additionalProperties": False,
            },
        },
    }


__all__ = [
    "AgentToolSurfaceV1",
    "build_approved_execution_tool_surface",
    "build_command_compiled_tool_surface",
    "build_single_agent_baseline_tool_surface",
]
