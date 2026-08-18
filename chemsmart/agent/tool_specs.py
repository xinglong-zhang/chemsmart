"""Approved model-visible command-compiled tool surface."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.skills import skills_enabled
from chemsmart.agent.capabilities import (
    ProgramCapabilityRegistryV1,
    load_program_capabilities,
)
from chemsmart.analysis.quantity_expressions import OPERATION_DESCRIPTIONS
from chemsmart.analysis.result_quantities import SUPPORTED_SELECTORS
from chemsmart.analysis.result_readers import (
    registered_reader_programs,
    registered_reader_selectors,
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


#: ``repair_command`` is implemented internally but is not advertised because
#: the public Runtime does not create its prerequisite counterexample record.
_REPAIR_COMMAND_TOOL_SOURCE = """_tool(
            "repair_command",
            (
                "Apply one constrained counterexample patch without changing "
                "bindings. PRECONDITION: a counterexample must already"""


def build_command_compiled_tool_surface(
    registry: ProgramCapabilityRegistryV1 | None = None,
) -> AgentToolSurfaceV1:
    registry = registry or load_program_capabilities()
    programs = [item.program for item in registry.programs]
    program = {"type": "string", "enum": programs}
    result_programs = tuple(sorted({"pyscf", *registered_reader_programs()}))
    reader_selector_inventory = "; ".join(
        f"{name}: {', '.join(selectors)}"
        for name, selectors in registered_reader_selectors().items()
    )
    structured_result_program = {
        "type": "string",
        "enum": list(result_programs),
        "description": (
            "Select the parser matching the registered artifact. Current "
            "program-wide reader selector union (not a promise for every job "
            "type): "
            f"{reader_selector_inventory}. PySCF uses its structured HDF5 "
            "result registry. Query inspect_program_capability for job-scoped "
            "parser support where declared; the selected method/settings must "
            "still emit the quantity."
        ),
    }
    thermochemistry_program = {
        "type": "string",
        # A geometry-only XYZ artifact can provide coordinates and an
        # embedded electronic energy, but never a Hessian thermochemistry
        # result.  Keep that truthful distinction in the model surface.
        "enum": [item for item in result_programs if item != "xyz"],
    }
    result_selector = {
        "type": "string",
        "enum": sorted(SUPPORTED_SELECTORS),
        "description": (
            "Program-neutral semantic selector. Support is resolved by the "
            "registered parser for the bound program artifact; a selector that "
            "the program or result does not provide remains explicitly blocked. "
            "The connectivity selector returns binary geometry-perceived "
            "adjacency in source atom order from covalent radii; it is not an "
            "electronic bond-order assignment."
        ),
    }
    digest = {"type": "string", "pattern": "^[0-9a-f]{64}$"}
    tools = (
        _tool(
            "inspect_program_capability",
            (
                "Join core capability, active support overlay, and live Click "
                "schema. For an exact job type, the receipt also returns "
                "parser-established job_result_selector_coverage when declared; "
                "coverage_semantics=parser_supported_when_emitted means the "
                "selected method/settings must still emit those quantities. "
                "This is not evidence that a particular run completed."
            ),
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
                "write it. Sections must be an object mapping section names "
                "to setting-name/value objects, never a list. PySCF and xTB "
                "use job sections. Gaussian and ORCA retain gas/solv phase "
                "sections plus optional job overrides: SP reads solv when "
                "present, otherwise gas, and an explicit sp override wins; "
                "the section name solv does not enable solvation by itself. "
                "For ORCA, ab_initio is the method field for HF-family and "
                "correlated wave-function methods; reference only selects the "
                "SCF determinant and does not replace ab_initio or functional. "
                "CPU count and memory belong to the ChemSmart run/server "
                "layer, not project additional_route_parameters; use that "
                "escape hatch only for source-required scientific keywords "
                "that refine an otherwise supported typed method. It is not "
                "a way to introduce, duplicate, or override an electronic-"
                "structure method or basis that is absent from the program "
                "capability. Keep such a source-exact stage explicitly "
                "unsupported and label any typed, supported surrogate."
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
            (
                "Validate a bound project through the checked-out settings "
                "loader and report which YAML sections actually feed the "
                "requested job. Loader-valid with no explicit applied settings "
                "requires a project-section repair before planning."
            ),
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
            ("workflow_id", "nodes"),
        ),
        _tool(
            "plan_scientific_workflow",
            (
                "Plan one connected scientific tool chain containing any "
                "required program calculations and deterministic analysis "
                "stages. For an analysis-only task over registered results, "
                "calculation_nodes may be empty; do not invent a documentary "
                "or blocked calculation placeholder. "
                "Future analysis inputs name producer node/output pairs; they "
                "do not require artifact or receipt hashes before execution. "
                "A result-extraction or thermochemistry root may instead "
                "consume one existing host-registered result by artifact_id. "
                "Keep unsupported requested analyses as blocked_unsupported nodes."
            ),
            {
                "plan_id": _public_identifier(),
                "workflow_id": _public_identifier(),
                "task_spec_id": _string(),
                "calculation_nodes": {
                    "type": "array",
                    "minItems": 0,
                    "maxItems": 64,
                    "items": _scientific_workflow_node_schema(),
                    "description": (
                        "Program calculations required by the task. Use an empty "
                        "array when every scientific root is an existing "
                        "host-registered result."
                    ),
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
                "analysis_nodes",
                "required_output_ids",
            ),
        ),
        _tool(
            "amend_scientific_workflow",
            (
                "Repair how part of the latest scientific workflow is "
                "expressed, without resubmitting the whole DAG. Use it after a "
                "rejection that names a field: an identifier's case, a missing "
                "or wrong unit, a declared quantity kind the operation does not "
                "derive, a selector the result does not resolve, or a project "
                "role whose corrected project has since been promoted and "
                "validated. The host preserves every node, binding, dependency "
                "and receipt you do not name. It refuses anything that changes "
                "the science rather than its expression: molecular identity, "
                "state, program, job type, an analysis kind, which producer an "
                "input reads from, thermochemical conditions, and validation "
                "thresholds are all a new workflow and need their own review."
            ),
            {
                "workflow_id": _public_identifier(),
                "project_replacements": {
                    "type": "array",
                    "description": (
                        "Project-only repairs. Each item binds one existing "
                        "calculation node to a newly promoted project role; "
                        "all scientific inputs, outputs, and dependencies are "
                        "preserved by the host."
                    ),
                    "maxItems": 64,
                    "items": {
                        "type": "object",
                        "properties": {
                            "node_id": _public_identifier(),
                            "project_role": _public_identifier(),
                        },
                        "required": ["node_id", "project_role"],
                        "additionalProperties": False,
                    },
                },
                "analysis_repairs": {
                    "type": "array",
                    "description": (
                        "Expression repairs to named analysis nodes. Each "
                        "entry addresses elements that already exist on that "
                        "node and replaces only the fields named."
                    ),
                    "maxItems": 128,
                    "items": {
                        "type": "object",
                        "properties": {
                            "node_id": _string(),
                            "outputs": {
                                "type": "array",
                                "maxItems": 64,
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        # An address into what the node
                                        # already declares, not a name being
                                        # authored, so the naming rule does
                                        # not apply and is not repeated here.
                                        "output_id": _string(),
                                        "new_output_id": _public_identifier(),
                                        "quantity_kind": _string(),
                                        "unit": _unit_string(
                                            "Corrected unit for this output."
                                        ),
                                    },
                                    "required": ["output_id"],
                                    "additionalProperties": False,
                                },
                            },
                            "selectors": {
                                "type": "array",
                                "maxItems": 64,
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "quantity_id": _string(),
                                        "selector": _string(),
                                    },
                                    "required": ["quantity_id", "selector"],
                                    "additionalProperties": False,
                                },
                            },
                            "inputs": {
                                "type": "array",
                                "maxItems": 64,
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "input_id": _string(),
                                        "producer_output_id": _string(),
                                    },
                                    "required": [
                                        "input_id",
                                        "producer_output_id",
                                    ],
                                    "additionalProperties": False,
                                },
                            },
                        },
                        "required": ["node_id"],
                        "additionalProperties": False,
                    },
                },
            },
            ("workflow_id",),
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
                "bindings from the typed workflow; do not copy receipt hashes. "
                "For a multi-file CLI job, the workflow uses binding_id "
                "'filename' for the primary molecular geometry and the exact "
                "live ChemSmart job-option name for each additional registered "
                "artifact (for example 'ending_xyzfile')."
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
                        "extraction, thermochemistry, quantity-expression, "
                        "scientific-validation, or analysis-claim "
                        "tools. The host validates and canonicalizes them; do "
                        "not embed receipt IDs in free-form evidence strings."
                    ),
                },
            },
            (
                "decision_id",
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
                "Derive harmonic RRHO and, when requested, Grimme/Truhlar "
                "quasi-harmonic thermochemistry from a trusted frequency result "
                "using ChemSmart's common engine. A supplied concentration "
                "defines the translational standard state instead of pressure. "
                "Grimme or Truhlar requires entropy_cutoff_cm1; an enthalpy "
                "cutoff independently enables Head-Gordon qRRHO enthalpy. The "
                "receipt distinguishes thermal_enthalpy_correction = H(T) - "
                "E_electronic, which includes ZPE, from "
                "enthalpy_increment_above_zero_point = H(T) - E_electronic - "
                "ZPE. For Grimme/Truhlar or Head-Gordon treatment it also "
                "provides quasi_harmonic_thermal_gibbs_correction = G_qh(T) "
                "- E_electronic; use that quantity, not the harmonic "
                "thermal_gibbs_correction, when composing a high-level "
                "electronic energy with low-level qRRHO thermochemistry. Use "
                "enthalpy_increment_above_zero_point when adding a finite-"
                "temperature increment to an already ZPE-corrected 0 K "
                "quantity."
            ),
            {
                "program": thermochemistry_program,
                "artifact_id": _string(),
                "temperature_k": {"type": "number", "exclusiveMinimum": 0},
                "pressure_atm": {"type": "number", "exclusiveMinimum": 0},
                "concentration_mol_l": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "description": (
                        "Optional solution standard-state concentration in mol/L; "
                        "when supplied, pressure remains recorded but is not used "
                        "for the translational partition function."
                    ),
                },
                "entropy_method": {
                    "type": "string",
                    "enum": ["rrho", "grimme", "truhlar"],
                    "description": (
                        "Entropy treatment; omitted means harmonic RRHO. Grimme "
                        "and Truhlar require entropy_cutoff_cm1."
                    ),
                },
                "entropy_cutoff_cm1": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "description": (
                        "Positive low-frequency entropy cutoff in cm^-1, required "
                        "for Grimme or Truhlar and invalid for harmonic RRHO."
                    ),
                },
                "enthalpy_cutoff_cm1": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "description": (
                        "Optional positive Head-Gordon qRRHO enthalpy cutoff in "
                        "cm^-1; omission retains harmonic enthalpy."
                    ),
                },
                "alpha": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "Damping exponent for Grimme entropy and Head-Gordon "
                        "enthalpy corrections; omitted means 4."
                    ),
                },
                "use_weighted_mass": {
                    "type": "boolean",
                    "description": (
                        "Use natural-abundance weighted isotope masses; omitted "
                        "means the backward-compatible most-abundant masses."
                    ),
                },
                "frequency_scale_factor": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "description": (
                        "Positive multiplicative scale applied to every "
                        "vibrational frequency before harmonic or "
                        "quasi-harmonic thermochemistry; omitted means 1.0."
                    ),
                },
            },
            ("program", "artifact_id", "temperature_k", "pressure_atm"),
        ),
        _tool(
            "evaluate_quantity_expression",
            (
                "Evaluate a bounded dimension-aware expression DAG over prior "
                "extraction, thermochemistry, quantity-expression, or "
                "scientific-validation receipts, "
                "or over typed literal nodes when inputs is empty; Python and "
                "formula strings are not accepted. Expression inputs use local "
                "input_id aliases, so derived outputs can feed later expressions. "
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
                            "semantic_role": _semantic_role_identifier(),
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
            "evaluate_scientific_validation",
            (
                "Evaluate the rules already sealed into one planned "
                "scientific-validation node against exact upstream typed "
                "quantities. Supply only receipt bindings; predicates and "
                "thresholds cannot be restated or weakened at execution."
            ),
            {
                "workflow_id": _public_identifier(),
                "node_id": _public_identifier(),
                "inputs": {
                    "type": "array",
                    "minItems": 1,
                    "maxItems": 64,
                    "items": {
                        "type": "object",
                        "properties": {
                            "input_id": _public_identifier(),
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
            },
            ("workflow_id", "node_id", "inputs"),
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
                            "display_unit": _unit_string(
                                "Unit to display the bound quantity in."
                            ),
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
            ("claims",),
        ),
    )
    if not skills_enabled():
        # Keep the model-visible surface aligned with the host feature state.
        tools = tuple(
            item
            for item in tools
            if item["function"]["name"] != "consult_domain_skill"
        )
    # Advertise only what this runtime can actually deliver.  The handlers and
    # contracts stay, so restoring one of these is a producer away.
    tools = tuple(
        item for item in tools if not _requires_an_unbound_registry(item)
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
    # ``inspect_calculation_artifact`` belongs to the legacy externally
    # seeded verifier surface: it requires separate settings-object and run
    # receipt IDs. Runtime V2 execution instead creates a typed program result
    # validation receipt and registers the resulting artifact directly; it
    # never binds those legacy IDs. Advertising the verifier here made live
    # models guess impossible identifiers after an otherwise valid run.
    tools = tuple(
        item
        for item in full.tool_definitions
        if item["function"]["name"] != "inspect_calculation_artifact"
    ) + (
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


#: What each recurring argument name means, applied wherever that name appears.
#: Keying by argument name keeps the same field definition consistent across
#: every tool that accepts it.
ARGUMENT_DESCRIPTIONS: dict[str, str] = {
    "basis_mode": (
        "How the basis is specified: a single set, or split by element class."
    ),
    "claims": (
        "The reported values, each bound to the typed analysis "
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
        "The typed quantities this operation consumes, each bound to an "
        "exact upstream receipt."
    ),
    "job_families": "The job types this request covers.",
    "jobtype": (
        "The target program's ChemSmart CLI job form, not a program-neutral "
        "label. Use sp when the supplied geometry must remain fixed, opt only "
        "when a minimum geometry search is intended, and ts only when a "
        "transition-state search is intended. Project frequency or VPT2 "
        "settings request properties and do not change this geometry operation."
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
    "task_spec_id": (
        "Identifier of the task specification being planned. Omit it when "
        "the host has exactly one active task; multi-task hosts require the "
        "exact identifier."
    ),
    "task_spec_sha256": (
        "Digest of the task specification this binds to. Omit it when the "
        "host has exactly one active task; multi-task hosts require the exact "
        "digest."
    ),
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


#: Arguments the runtime passes through its public-identifier validator.  The
#: schema publishes the constraint so callers can construct valid requests.
IDENTIFIER_ARGUMENTS = frozenset(
    {
        "artifact_class",
        "counterexample_id",
        "jobtype",
        "node_id",
        "output_id",
        "program",
        "project_role",
        "workflow_id",
    }
)

#: Fields that are an identifier *or* deliberately empty.  The two halves of a
#: workflow edge are mutually exclusive: an input bound to an initial artifact
#: names artifact_id and leaves the producer fields empty, while an input fed by
#: an upstream node names the producer fields and leaves artifact_id empty.
#: Constraining either half to the plain identifier pattern forbids the other.
#:
#: Membership here is not a judgement call.  It is what the runtime validators
#: actually accept, checked by probing them rather than by reading the code --
#: artifact_id was misclassed on the first attempt and a live session hit it.
OPTIONAL_IDENTIFIER_ARGUMENTS = frozenset(
    {"artifact_id", "producer_node_id", "producer_output_id"}
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
    if (
        schema.get("type") != "string"
        or schema.get("enum")
        or schema.get("pattern")
    ):
        return schema
    return {**schema, "pattern": pattern}


#: What binds each kind of host-owned object, phrased to complete the sentence
#: "no <label> is bound yet; one is bound <producer>".  A caller told only
#: that a registry is empty cannot act; a caller told what fills it can.
REGISTRY_PRODUCERS: dict[str, str] = {
    "canonical invocation": "by compiling a prepared program node",
    "capability receipt": "by inspecting a program capability",
    "command context": "by preparing a program node",
    "command inspection receipt": "by compiling a program node",
    "counterexample": (
        "by the host when a compiled command fails inspection, safe preview, "
        "or program validation -- do not reference one before a failure has "
        "produced it"
    ),
    "engine binding": "by inspecting the program environment",
    "functional equivalence receipt": "by validating a project document",
    "program binding": "by inspecting the program environment",
    "program validator receipt": "by safe-previewing a compiled command",
    "project render receipt": "by rendering a project YAML document",
    "project validation receipt": "by validating a rendered project",
    "run receipt": "by executing an approved node",
    "safe preview receipt": "by safe-previewing a compiled command",
    "scientific claim evidence": "by extracting quantities from a result",
    "scientific identity": "by binding an approved molecular identity",
    "settings object": "by validating a project document",
    "trusted artifact": "by recording a workspace file as an artifact",
}


#: Which host registry each late-bound argument indexes.  A tool taking one of
#: these cannot succeed until something else has run, and saying so only in the
#: rejection means the model learns it by failing.  Observed across six live
#: sessions: repair_command called with no counterexample bound (four times)
#: and assess_program_candidate called with no claim evidence bound.
LATE_BOUND_ARGUMENTS: dict[str, str] = {
    "command_inspection_receipt_sha256": "command inspection receipt",
    "counterexample_id": "counterexample",
    "invocation_sha256": "canonical invocation",
    "render_receipt_sha256": "project render receipt",
    "run_receipt_id": "run receipt",
    "settings_id": "settings object",
    "source_claim_sha256s": "scientific claim evidence",
}


#: Registries the externally-seeded V1 surface filled but Runtime V2 never
#: binds: nothing writes them during a session, and no live entry point seeds
#: them at construction.  A tool requiring one cannot succeed here at all, so
#: advertising it can only invite a guess at an identifier that will never
#: exist -- the same defect ``repair_command`` had, and the same rule the
#: capability registry applies when it declares an unsupported job type instead
#: of offering it as runnable.
#:
#: This is deliberately a property of the *registry* rather than a list of tool
#: names, so a tool becomes reachable again by giving its registry a producer,
#: not by editing a second list that can drift from the first.
UNBOUND_RUNTIME_V2_REGISTRIES: frozenset[str] = frozenset(
    {
        "run receipt",
        "scientific claim evidence",
        "settings object",
    }
)


def _requires_an_unbound_registry(definition: dict) -> bool:
    """Whether this tool can never succeed on the active runtime."""

    parameters = definition["function"].get("parameters") or {}
    for name in parameters.get("required", ()):
        label = LATE_BOUND_ARGUMENTS.get(name)
        if label in UNBOUND_RUNTIME_V2_REGISTRIES:
            return True
    return False


def _precondition_sentence(properties) -> str:
    """State what must already exist before this tool can succeed."""

    parts = []
    for name in sorted(properties):
        label = LATE_BOUND_ARGUMENTS.get(name)
        if not label:
            continue
        producer = REGISTRY_PRODUCERS.get(label)
        if producer:
            parts.append(f"{name} names a {label}, which is bound {producer}")
    if not parts:
        return ""
    return (
        " PRECONDITION: "
        + "; ".join(parts)
        + ". Calling this before that has happened cannot succeed."
    )


def _describe_tool_definitions(
    definitions: tuple[dict, ...],
) -> tuple[dict, ...]:
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
        sentence = _precondition_sentence(parameters.get("properties") or {})
        if sentence and "PRECONDITION" not in function.get("description", ""):
            function["description"] = function["description"] + sentence
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


def _unit_string(lead: str) -> dict:
    """A unit field that states the convention where the unit is written.

    The typed vocabulary takes units, not quantity names and not rescalings,
    and a dimensionless quantity is spelled ``1``.  Chemists write "percent",
    "count" and "mole fraction" by habit, so say so here rather than only in
    the refusal: this text is read on every call, while a refusal is only read
    once the call has already been spent.
    """

    return {
        "type": "string",
        "description": (
            f"{lead} Give a unit, not a quantity name and not a rescaling. "
            "A dimensionless quantity -- a count, population, mole fraction, "
            "branching ratio, equilibrium constant or oscillator strength -- "
            "takes '1'; 'percent' is a rescaling of a dimensionless value, "
            "not a unit, so report the fraction and describe it as a "
            "percentage in prose. Energies accept hartree, eV, kJ/mol or "
            "kcal/mol; frequencies cm^-1; temperatures K."
        ),
    }


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
            "instead of spaces, parentheses, hashes, or placeholder syntax. "
            "It must begin with a letter, so a name taken from a compound "
            "whose locants come first needs a leading word: "
            "'dfe-12-rotamers', not '12-difluoroethane'. Chemical notation is "
            "mixed case and this field is not: unit symbols and quantity "
            "names must be folded down, so write 'gap-adiab-ev' not "
            "'gap-adiab-eV', 'delta-e' not 'dE', and 'ddg-compose' not "
            "'compose-ddG'. Fold the case; do not drop the letters."
        ),
    }


def _semantic_role_identifier() -> dict:
    """An optional readable label for one occurrence of a repeated quantity.

    An expression identifies its inputs by the quantity they came from, and
    the same named quantity routinely arrives twice -- the two conformers of a
    population, the reactant and product of a difference, the two states of a
    gap.  The host resolves that repetition itself, from each input's own id,
    which the expression contract already requires to be unique.

    This field therefore buys readability, not correctness.  Demanding it
    instead cost five cycles of refusals: the requirement was stated in terms
    of receipt internals the model cannot see, and each attempt to explain it
    better halved the failure without clearing it.  Describe what supplying a
    role is *for*, and let the host derive what it can derive.
    """

    identifier = _public_identifier()
    return {
        **identifier,
        "description": (
            "Optional. Which occurrence this input is, when the same source "
            "quantity is drawn more than once. You do not have to supply it: "
            "the host falls back to this input's own id, which is unique "
            "within the expression, so an omitted role is never ambiguous. "
            "Supply one only to label an occurrence more readably than its "
            "input id does, as 'gauche-gibbs' rather than 'in7'. If you do "
            "supply roles, keep them distinct -- two inputs sharing one role "
            "would collapse onto a single slot and make the evidence "
            "reference ambiguous, so that is refused. Never use a receipt "
            "hash as a role. " + identifier["description"]
        ),
    }


def _internal_coordinates_schema() -> dict:
    """Which internal coordinates this node scans or holds fixed.

    A scanned dihedral or a frozen bond is a fact about *this molecule in this
    calculation*, the same class of fact as charge and multiplicity, so it
    belongs on the node rather than frozen into a reusable project.  The
    specification here is physical and program-neutral: the host renders it
    into each program's own idiom, which genuinely differ -- ORCA takes an
    absolute range, Gaussian an increment.
    """

    atoms = {
        "type": "array",
        "items": {"type": "integer", "minimum": 1},
        "minItems": 2,
        "maxItems": 4,
        "description": (
            "The atoms defining this coordinate, numbered from 1 in the "
            "order of the bound geometry: two for a bond, three for an "
            "angle, four for a dihedral."
        ),
    }
    kind = {
        "type": "string",
        "enum": ["bond", "angle", "dihedral"],
        "description": "Which internal coordinate these atoms define.",
    }
    return {
        "type": "object",
        "description": (
            "Internal coordinates this node scans or constrains. Required by "
            "a scan or a constrained optimisation and meaningless without "
            "one; the geometry itself, the method and the basis come from the "
            "bound artifact and the project, not from here."
        ),
        "properties": {
            "scan": {
                "type": "object",
                "description": (
                    "The one coordinate driven across a range. State the "
                    "range physically -- angstrom for a bond, degrees for an "
                    "angle or dihedral -- and the host renders it as the "
                    "target program expects."
                ),
                "properties": {
                    "kind": kind,
                    "atoms": atoms,
                    "start": {
                        "type": "number",
                        "description": "First value of the driven coordinate.",
                    },
                    "stop": {
                        "type": "number",
                        "description": "Last value of the driven coordinate.",
                    },
                    "points": {
                        "type": "integer",
                        "minimum": 2,
                        "description": (
                            "How many values are computed, endpoints "
                            "included."
                        ),
                    },
                },
                "required": ["kind", "atoms", "start", "stop", "points"],
                "additionalProperties": False,
            },
            "constrained": {
                "type": "array",
                "maxItems": 32,
                "description": (
                    "Coordinates held fixed while everything else relaxes. "
                    "These are the constraint of a constrained optimisation, "
                    "and may also accompany a scan."
                ),
                "items": {
                    "type": "object",
                    "properties": {"kind": kind, "atoms": atoms},
                    "required": ["kind", "atoms"],
                    "additionalProperties": False,
                },
            },
        },
        "additionalProperties": False,
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
                    "For ORCA, freq and vpt2 remain project properties and "
                    "ORCA exposes no hess CLI jobtype. PySCF uses a separate "
                    "hess node."
                ),
            },
            "node_kind": {
                "type": "string",
                # Legacy Runtime V2 events may still replay an aggregate
                # command node, but new model proposals use the scientific
                # analysis DAG as the single post-processing authority.
                "enum": ["program_call"],
                "description": (
                    "Omit or use 'program_call' to invoke a program. Declare "
                    "post-processing with the scientific toolchain's typed "
                    "analysis nodes, not a second command-level aggregate "
                    "plane."
                ),
            },
            "project_role": _string(),
            "charge": {
                "type": "integer",
                "description": (
                    "Optional explicit charge for this node. Supply it only "
                    "together with multiplicity. On an optimized-geometry "
                    "data-edge consumer this deliberately reuses the exact "
                    "producer geometry on another electronic-state surface."
                ),
            },
            "multiplicity": {
                "type": "integer",
                "minimum": 1,
                "description": (
                    "Optional explicit spin multiplicity for this node. "
                    "Supply it only together with charge; omission inherits "
                    "the task-bound state of the molecular input."
                ),
            },
            "internal_coordinates": _internal_coordinates_schema(),
            "dependencies": {"type": "array", "items": _string()},
            "inputs": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "binding_id": {
                            **_public_identifier(),
                            "description": (
                                "Semantic input role. For a single input any "
                                "unique public role is valid. For a multi-file "
                                "program call use 'filename' for the primary "
                                "geometry and the exact live ChemSmart job-option "
                                "parameter for additional artifacts, such as "
                                "'ending_xyzfile' for an ORCA NEB product. "
                                "An ORCA IRC that reads the final transition-"
                                "state Hessian has exactly two producer inputs: "
                                "'filename'/geometry_xyz and "
                                "'hess_filename'/orca_hessian."
                            ),
                        },
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
                "description": (
                    "Typed program artifacts made available to downstream "
                    "nodes. Future typed result extraction or thermochemistry "
                    "must consume a reader artifact: orca_output, "
                    "gaussian_output, xtb_output, or pyscf_hdf5 for the "
                    "matching program. Declare it as a separate expected "
                    "output. Native handoff artifacts have different roles: "
                    "for example, an ORCA TS with frequencies may also "
                    "declare final_hessian/orca_hessian for a downstream IRC, "
                    "and ChemSmart selects the Hessian bound to the validated "
                    "final TS rather than guessing a filename."
                ),
            },
            "unresolved_fields": {
                "type": "array",
                "items": _public_identifier(),
            },
            "produces_observables": {
                "type": "array",
                "items": _public_identifier(),
                "description": (
                    "Scientific quantities produced by this calculation. Match "
                    "the loader-effective settings returned by "
                    "validate_project_yaml. In particular, when that tool says "
                    "the project already requests frequencies, put "
                    "vibrational_frequencies on this node instead of scheduling "
                    "a duplicate Hessian at the same geometry and method."
                ),
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
            "artifact_id": {
                "type": "string",
                "description": (
                    "For a result_extraction or thermochemistry root, name "
                    "one existing registered result instead of a future "
                    "program output. Leave inputs empty. Do not supply a "
                    "path, hash, or program-native text."
                ),
            },
            "inputs": {
                "type": "array",
                "description": (
                    "Typed producer edges. A planned result_extraction or "
                    "thermochemistry node must name the producer's registered "
                    "result-reader output, not a native geometry or Hessian "
                    "handoff output."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "input_id": _public_identifier(),
                        "source_kind": {
                            "type": "string",
                            "enum": ["analysis_output", "program_output"],
                            "description": (
                                "Describe the immediate producer: use "
                                "program_output when producer_node_id names a "
                                "calculation node, and analysis_output when it "
                                "names an analysis node. A planned extraction "
                                "or thermochemistry input must be the "
                                "calculation's typed result-reader output."
                            ),
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
                        "unit": _unit_string(
                            "Physical unit this output is declared in."
                        ),
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
            "validation_rules": {
                "type": "array",
                "description": (
                    "Program-neutral validation predicates over named typed "
                    "inputs. Use minimum_greater_equal with threshold 0 cm-1 "
                    "for a no-imaginary-frequency requirement; do not hide "
                    "criteria only in prose. A scientific_validation node "
                    "declares exactly one dimensionless verdict output."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "rule_id": _public_identifier(),
                        "predicate": {
                            "type": "string",
                            "enum": [
                                "all_equal",
                                "all_finite",
                                "count_equals",
                                "maximum_absolute_less_equal",
                                "minimum_greater_equal",
                                "symmetric_within",
                            ],
                        },
                        "input_ids": {
                            "type": "array",
                            "minItems": 1,
                            "items": _public_identifier(),
                        },
                        "threshold": {"type": "number"},
                        "expected_count": {"type": "integer", "minimum": 0},
                        "unit": _string(),
                    },
                    "required": ["rule_id", "predicate", "input_ids"],
                    "additionalProperties": False,
                },
            },
            # Most analysis kinds have no thermodynamic state, and the node
            # contract accepts None for them. Refusing an explicit null while
            # accepting an omitted key makes an honest stateless analysis fail.
            "temperature_k": _nullable_positive_number(),
            "pressure_atm": _nullable_positive_number(),
            "concentration_mol_l": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Thermochemistry-only solution standard state in mol/L. "
                    "Omit for the pressure-defined ideal-gas standard state."
                ),
            },
            "entropy_method": {
                "type": "string",
                "enum": ["rrho", "grimme", "truhlar"],
                "description": (
                    "Thermochemistry-only entropy treatment; omitted means "
                    "harmonic RRHO."
                ),
            },
            "entropy_cutoff_cm1": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Thermochemistry-only entropy cutoff required for "
                    "Grimme or Truhlar treatment."
                ),
            },
            "enthalpy_cutoff_cm1": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Thermochemistry-only Head-Gordon qRRHO enthalpy cutoff."
                ),
            },
            "alpha": {
                "type": "integer",
                "minimum": 1,
                "description": (
                    "Thermochemistry-only damping exponent; omitted means 4."
                ),
            },
            "use_weighted_mass": {
                "type": "boolean",
                "description": (
                    "Thermochemistry-only isotope-mass convention; omitted "
                    "uses most-abundant isotopes."
                ),
            },
            "frequency_scale_factor": {
                "type": "number",
                "exclusiveMinimum": 0,
                "description": (
                    "Thermochemistry-only positive multiplicative frequency "
                    "scale; omitted means 1.0."
                ),
            },
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
                    "SCF exponential, SCF inverse-power, or correlation "
                    "inverse-power CBS limit. "
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
]
