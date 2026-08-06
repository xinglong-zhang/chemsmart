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
from chemsmart.analysis.quantity_expressions import OPERATION_DESCRIPTIONS
from chemsmart.analysis.result_quantities import SUPPORTED_SELECTORS
from chemsmart.analysis.result_readers import registered_reader_programs


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
    result_programs = tuple(sorted({"pyscf", *registered_reader_programs()}))
    structured_result_program = {
        "type": "string",
        "enum": list(result_programs),
    }
    result_selector = {
        "type": "string",
        "enum": sorted(SUPPORTED_SELECTORS),
        "description": (
            "Program-neutral semantic selector. Support is resolved by the "
            "registered parser for the bound program artifact; a selector that "
            "the program or result does not provide remains explicitly blocked."
        ),
    }
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
            "plan_scientific_workflow",
            (
                "Plan one connected paper-level tool chain containing both "
                "program calculations and deterministic analysis stages. "
                "Future analysis inputs name producer node/output pairs; they "
                "do not require artifact or receipt hashes before execution. "
                "Keep unsupported requested analyses as blocked_unsupported nodes."
            ),
            {
                "plan_id": _public_identifier(),
                "workflow_id": _public_identifier(),
                "task_spec_id": _string(),
                "calculation_nodes": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 64,
                    "items": _scientific_workflow_node_schema(),
                },
                "analysis_nodes": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 128,
                    "items": _analysis_intent_node_schema(),
                },
                "required_output_ids": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 64,
                    "items": _public_identifier(),
                },
            },
            (
                "plan_id",
                "workflow_id",
                "task_spec_id",
                "calculation_nodes",
                "analysis_nodes",
                "required_output_ids",
            ),
        ),
        _tool(
            "inspect_workflow_frontier",
            (
                "Inspect the current calculation-and-analysis frontier for a "
                "previously planned scientific workflow."
            ),
            {"workflow_id": _public_identifier()},
            ("workflow_id",),
        ),
        _tool(
            "prepare_program_node",
            (
                "Prepare and safe-preview one actionable calculation node "
                "from a scientific workflow. The host resolves its program, "
                "project, input, electronic state, capability, and engine "
                "bindings from the typed workflow; do not copy receipt hashes."
            ),
            {
                "workflow_id": _public_identifier(),
                "node_id": _public_identifier(),
            },
            ("workflow_id", "node_id"),
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
            (
                "Apply one constrained counterexample patch without changing "
                "bindings. PRECONDITION: a counterexample must already exist, "
                "and only the host creates one -- when a compiled command "
                "fails inspection, safe preview, or program validation. There "
                "is nothing to repair until something has failed, so calling "
                "this after a successful preview cannot succeed. Fix a "
                "rejected argument by calling the tool again with it "
                "corrected, not by repairing it."
            ),
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
                            "selector": result_selector,
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
    tools = _describe_tool_definitions(tools)
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
    tools = _describe_tool_definitions(tools)
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


#: What each recurring argument name means, applied wherever that name appears.
#:
#: Measured before this existed: 7 of 93 required arguments across the whole
#: model-visible surface carried a description, so a model had to infer 86
#: mandatory fields from their spelling.  Two live sessions failed on
#: ``record_scientific_decision`` for exactly that reason -- one omitted
#: ``alternatives``, another wrote a sentence into ``stage_order``.
#:
#: Keying by argument name rather than by tool is deliberate: the same field
#: means the same thing everywhere it appears, and one entry then describes it
#: in every tool that takes it.
ARGUMENT_DESCRIPTIONS: dict[str, str] = {
    "basis_mode": (
        "How the basis is specified: a single set, or split by element class."
    ),
    "claims": (
        "The reported values, each bound to the extraction or expression "
        "receipt that produced it and carrying its display unit."
    ),
    "command_inspection_receipt_sha256": (
        "Digest of the inspection receipt for the compiled command."
    ),
    "constraint_kinds": "The geometric constraints this request needs.",
    "counterexample_id": (
        "ID of a counterexample the host produced when a compiled command "
        "failed inspection, safe preview or program validation. One exists "
        "only after such a failure."
    ),
    "decision_id": (
        "Stable identifier for this decision, lower case. Reuse it when you "
        "revise the same decision so the record supersedes rather than "
        "duplicates."
    ),
    "invocation_sha256": (
        "Digest of the canonical command invocation this acts on."
    ),
    "method_family": (
        "The broad method class, for example dft, hartree_fock or post_hf."
    ),
    "method_name": "The specific functional or method literal.",
    "render_receipt_sha256": (
        "Digest of the render receipt for the project document being promoted."
    ),
    "requires_double_hybrid": (
        "Whether the protocol needs a double hybrid, which several programs "
        "and engines do not support."
    ),
    "requires_post_hf": (
        "Whether the protocol needs a correlated wavefunction method beyond "
        "Hartree-Fock."
    ),
    "sections": (
        "The project YAML body, keyed by the program's own section names. A "
        "phase-keyed program uses gas/solv; a stage-keyed program uses its "
        "job-type names."
    ),
    "source_claim_sha256s": (
        "Digests of the claims this assessment is derived from."
    ),
    "alternatives": (
        "The other scientifically defensible options you considered and did "
        "not take, each with the reason. Required: a decision with no "
        "alternatives is a preference, not a decision."
    ),
    "analysis_nodes": (
        "The extraction, validation, mathematics and reporting stages that "
        "turn finished results into the requested values."
    ),
    "artifact_id": "The host-bound ID of an artifact already recorded.",
    "assumptions": (
        "What you are taking as given and did not verify, stated so a reader "
        "can check them independently."
    ),
    "calculation_nodes": (
        "The program invocations of the workflow, each naming its program, "
        "job type and project role."
    ),
    "capability_receipt_sha256": (
        "Digest of the capability receipt returned by "
        "inspect_program_capability for this program and engine."
    ),
    "charge": "Total molecular charge as an integer.",
    "diagnostics": (
        "The checks you will use to tell whether this decision was right, "
        "stated before the results exist."
    ),
    "engine": "Execution engine: 'cpu' or 'gpu'.",
    "engine_binding_sha256": (
        "Digest of the engine binding from inspect_program_environment."
    ),
    "evidence_refs": (
        "Digests of the host receipts this decision rests on. Every claim in "
        "the decision must be traceable to one."
    ),
    "geometry_artifact_sha256": (
        "Digest of the coordinate bytes this node consumes."
    ),
    "input_artifact_id": "The host-bound ID of the input geometry artifact.",
    "inputs": (
        "The measured quantities the expression consumes, each bound to an "
        "extraction receipt."
    ),
    "job_families": "The job types this request covers.",
    "jobtype": (
        "The target program's ChemSmart CLI job form, not a program-neutral "
        "label."
    ),
    "method_rationale": (
        "Why this method and these settings answer the question, in the terms "
        "the protocol being reproduced uses."
    ),
    "multiplicity": "Spin multiplicity 2S+1 as an integer, not PySCF's spin.",
    "node_id": "Stable identifier of this node within the workflow.",
    "nodes": "The workflow's nodes, in the order you intend them to run.",
    "output_node_ids": (
        "Which expression nodes are the reported outputs; the rest are "
        "intermediates."
    ),
    "pressure_atm": "Standard-state pressure in atmospheres.",
    "program": (
        "The executable program name as ChemSmart registers it, lower case."
    ),
    "program_binding_sha256": (
        "Digest of the program binding from inspect_program_environment."
    ),
    "project_artifact_id": (
        "The host-bound ID of the promoted project YAML this node uses."
    ),
    "required_output_ids": (
        "The observables the task asked for. A workflow that cannot produce "
        "one of these is incomplete, whatever else it computes."
    ),
    "requested_engine": "The engine the task implies, before selection.",
    "requested_program": "The program the task implies, before selection.",
    "request_id": "Stable identifier for this request.",
    "run_receipt_id": "The host-bound ID of the execution receipt.",
    "scientific_identity_sha256": (
        "Digest of the approved molecular identity binding this node uses."
    ),
    "selected_engine": "The engine you chose, which may differ from requested.",
    "selected_program": (
        "The program you chose, which may differ from requested."
    ),
    "selectors": (
        "Which registered quantities to read from the result file, by name."
    ),
    "settings_id": "The host-bound ID of the validated settings object.",
    "stage_order": (
        "The stages in the order they run, each a lower-case identifier and "
        "nothing else. Do not put dependency prose here -- the host states "
        "dependencies itself in the workflow frontier."
    ),
    "task_spec_id": "Identifier of the task specification being planned.",
    "task_spec_sha256": "Digest of the task specification this binds to.",
    "temperature_k": "Temperature in kelvin.",
    "uncertainties": (
        "What could still make this wrong, and what would resolve it."
    ),
    "workflow_id": "Stable identifier for this workflow.",
}


def _describe(name: str, schema: dict) -> dict:
    """Give ``name`` its meaning, keeping any format rule already stated.

    The two are different things and a caller needs both: the shared entry says
    what the argument is for, while an existing description usually says what
    shape it must take.  Composing them in a fixed order keeps one argument
    name reading identically everywhere it appears.
    """

    meaning = ARGUMENT_DESCRIPTIONS.get(name)
    if meaning is None:
        return schema
    existing = str(schema.get("description") or "").strip()
    if not existing:
        return {**schema, "description": meaning}
    if meaning in existing:
        return schema
    return {**schema, "description": f"{meaning} {existing}"}


#: Arguments the runtime passes through its public-identifier validator.
#:
#: Measured: 36 of 37 identifier-guarded arguments were exposed as
#: unconstrained strings, so the rule existed only after submission.  Two live
#: sessions wrote whole sentences into identifier fields -- a dependency phrase
#: into ``stage_order``, a method description into ``project_role`` -- and only
#: learned the rule from the rejection.  Declaring the pattern states it before
#: the call instead of after it.
IDENTIFIER_ARGUMENTS = frozenset(
    {
        "artifact_class",
        "artifact_id",
        "counterexample_id",
        "jobtype",
        "node_id",
        "output_id",
        "program",
        "project_role",
        "workflow_id",
    }
)

#: Fields that are an identifier *or* deliberately empty.  An input bound to an
#: initial artifact has no producer, and says so with "".  Constraining these to
#: the plain identifier pattern would forbid the commonest edge in any workflow.
OPTIONAL_IDENTIFIER_ARGUMENTS = frozenset(
    {"producer_node_id", "producer_output_id"}
)

_PUBLIC_IDENTIFIER_PATTERN = "^[a-z][a-z0-9_.-]*$"
_OPTIONAL_IDENTIFIER_PATTERN = "^$|^[a-z][a-z0-9_.-]*$"


def _constrain(name: str, schema: dict) -> dict:
    """Declare the identifier rule the runtime will enforce anyway."""

    if not isinstance(schema, dict):
        return schema
    if name in OPTIONAL_IDENTIFIER_ARGUMENTS:
        pattern = _OPTIONAL_IDENTIFIER_PATTERN
    elif name in IDENTIFIER_ARGUMENTS:
        pattern = _PUBLIC_IDENTIFIER_PATTERN
    else:
        return schema
    if schema.get("type") != "string" or schema.get("enum") or schema.get("pattern"):
        return schema
    return {**schema, "pattern": pattern}


def _describe_tool_definitions(definitions: tuple[dict, ...]) -> tuple[dict, ...]:
    """Describe every argument the surface exposes, by argument name."""

    described = []
    for item in definitions:
        function = dict(item["function"])
        parameters = dict(function.get("parameters") or {})
        properties = parameters.get("properties")
        if isinstance(properties, dict):
            parameters["properties"] = {
                name: _describe(name, _walk_constrain(name, schema))
                if isinstance(schema, dict)
                else schema
                for name, schema in properties.items()
            }
            function["parameters"] = parameters
        described.append({**item, "function": function})
    return tuple(described)


def _walk_constrain(name: str, schema: dict) -> dict:
    """Apply the identifier rule at every depth, including inside arrays."""

    schema = _constrain(name, schema)
    if not isinstance(schema, dict):
        return schema
    updated = dict(schema)
    properties = updated.get("properties")
    if isinstance(properties, dict):
        updated["properties"] = {
            key: _describe(key, _walk_constrain(key, value))
            if isinstance(value, dict)
            else value
            for key, value in properties.items()
        }
    items = updated.get("items")
    if isinstance(items, dict):
        updated["items"] = _walk_constrain(name, items)
    return updated


def _string() -> dict:
    return {"type": "string"}


def _nullable_positive_number() -> dict:
    """A positive number, or an explicit null where the concept does not apply.

    A field the node contract types as ``float | None`` must accept null on the
    wire.  Otherwise omitting the key succeeds and saying null fails, which
    makes the explicit statement the rejected one.
    """

    return {
        "type": ["number", "null"],
        "exclusiveMinimum": 0,
        "description": (
            "Positive value, or null when this stage has no thermodynamic "
            "state."
        ),
    }


def _public_identifier() -> dict:
    return {
        "type": "string",
        "pattern": "^[a-z][a-z0-9_.-]*$",
        "description": (
            "Lower-case public identifier; use dots, dashes, or underscores "
            "instead of spaces, parentheses, hashes, or placeholder syntax."
        ),
    }


def _workflow_node_schema() -> dict:
    return {
        "type": "object",
        "properties": {
            "node_id": _string(),
            "program": _string(),
            "jobtype": {
                "type": "string",
                "description": (
                    "Use the target program's ChemSmart CLI job form, not a "
                    "program-neutral conceptual label. In particular, ORCA "
                    "harmonic frequencies are requested by freq: true in an "
                    "opt project and remain one opt node; ORCA has no hess CLI "
                    "jobtype. PySCF uses a separate hess node."
                ),
            },
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
                        "binding_id": _public_identifier(),
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
            "produces_observables": {
                "type": "array",
                "items": _public_identifier(),
            },
            "support_state": {
                "type": "string",
                "enum": ["planned", "blocked_unsupported"],
            },
            "blocked_reason": _string(),
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


def _scientific_workflow_node_schema() -> dict:
    """Calculation node schema with explicit scientific output semantics."""

    schema = _workflow_node_schema()
    schema["required"] = list(schema["required"]) + [
        "produces_observables",
        "support_state",
        "blocked_reason",
    ]
    return schema


def _analysis_intent_node_schema() -> dict:
    """Planning-only analysis node; artifacts are bound after producers run."""

    return {
        "type": "object",
        "properties": {
            "node_id": _public_identifier(),
            "analysis_kind": {
                "type": "string",
                "enum": [
                    "claim_rendering",
                    "quantity_expression",
                    "result_extraction",
                    "scientific_validation",
                    "thermochemistry",
                    "unsupported_external",
                ],
                "description": (
                    "Only result_extraction carries selectors; only "
                    "quantity_expression carries expression_nodes. Put a "
                    "numerical check in a quantity_expression producer and "
                    "feed its output to a scientific_validation node."
                ),
            },
            "dependencies": {
                "type": "array",
                "items": _public_identifier(),
            },
            "inputs": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "input_id": _public_identifier(),
                        "source_kind": {
                            "type": "string",
                            "enum": ["analysis_output", "program_output"],
                        },
                        "producer_node_id": _public_identifier(),
                        "producer_output_id": _public_identifier(),
                    },
                    "required": [
                        "input_id",
                        "source_kind",
                        "producer_node_id",
                        "producer_output_id",
                    ],
                    "additionalProperties": False,
                },
            },
            "selectors": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "quantity_id": _public_identifier(),
                        "selector": _public_identifier(),
                    },
                    "required": ["quantity_id", "selector"],
                    "additionalProperties": False,
                },
            },
            "outputs": {
                "type": "array",
                "minItems": 1,
                "items": {
                    "type": "object",
                    "properties": {
                        "output_id": _public_identifier(),
                        "quantity_kind": _public_identifier(),
                        "unit": _string(),
                    },
                    "required": ["output_id", "quantity_kind", "unit"],
                    "additionalProperties": False,
                },
            },
            "expression_nodes": {
                "type": "array",
                "items": _quantity_expression_node_schema(),
            },
            "expression_output_node_ids": {
                "type": "array",
                "items": _public_identifier(),
            },
            # Most analysis kinds have no thermodynamic state, and the node
            # contract accepts None for them.  Refusing an explicit null while
            # accepting an omitted key made the honest way to say "this stage
            # has no temperature" the one the schema rejected -- a live session
            # sent null on a result-extraction node and was refused for it.
            "temperature_k": _nullable_positive_number(),
            "pressure_atm": _nullable_positive_number(),
            "support_state": {
                "type": "string",
                "enum": ["planned", "blocked_unsupported"],
            },
            "blocked_reason": _string(),
        },
        "required": [
            "node_id",
            "analysis_kind",
            "dependencies",
            "inputs",
            "selectors",
            "outputs",
            "expression_nodes",
            "expression_output_node_ids",
            "support_state",
            "blocked_reason",
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
                "enum": sorted(OPERATION_DESCRIPTIONS),
                "description": (
                    "Pick the operation that owns the step. Where a named "
                    "operation exists for a scientific convention, use it "
                    "rather than rebuilding the convention from arithmetic "
                    "primitives: the named one carries the convention, its "
                    "validity conditions, and its provenance. "
                    + " | ".join(
                        f"{name}: {text}"
                        for name, text in sorted(
                            OPERATION_DESCRIPTIONS.items()
                        )
                    )
                ),
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
            "cardinal_numbers": {
                "type": "array",
                "items": {"type": "integer", "minimum": 2},
                "minItems": 2,
                "maxItems": 2,
                "description": (
                    "Increasing lower/higher basis cardinal numbers for a "
                    "two-point CBS operation."
                ),
            },
            "extrapolation_exponent": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Method/protocol-derived positive exponent for a two-point "
                    "SCF exponential or correlation inverse-power CBS limit. "
                    "This is a number you supply, not one the host measured, "
                    "so the receipt records it as model-authored and it is "
                    "auditable as such. Supply it only when the protocol you "
                    "are reproducing states it; when the protocol just says "
                    "the energy was extrapolated exponentially and you have "
                    "three successive cardinal numbers, prefer "
                    "exponential_cbs_limit, which fits the decay from the "
                    "data and introduces no constant of your own."
                ),
            },
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
