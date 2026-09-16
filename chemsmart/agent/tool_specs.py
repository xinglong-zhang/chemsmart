"""Approved model-visible command-compiled tool surface."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.capabilities import (
    ProgramCapabilityRegistryV1,
    load_program_capabilities,
)
from chemsmart.agent.execution import EDITABLE_COORDINATE_OPERATIONS
from chemsmart.agent.rules import render_rules
from chemsmart.agent.scientific_toolchain import (
    ANALYSIS_VALIDATION_PREDICATES,
)
from chemsmart.analysis.literature_constants import LITERATURE_CONSTANTS
from chemsmart.analysis.quantity_expressions import OPERATION_DESCRIPTIONS
from chemsmart.analysis.result_quantities import (
    SUPPORTED_SELECTORS,
    derivable_thermochemistry_quantities,
)
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


def _legacy_tool_definitions(
    registry: ProgramCapabilityRegistryV1 | None = None,
    *,
    operations: tuple[str, ...] | None = None,
) -> tuple[dict, ...]:
    """Every tool the host implements, one definition each, before the
    model-facing merge. The provider-free executor drives nodes through
    these names; the planning surface the model reads is built from them
    by :func:`_merge_planning_tools`."""

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
            f"{reader_selector_inventory}. Query inspect_program for "
            "job-scoped parser support where declared; the selected "
            "method/settings must still emit the quantity."
        ),
    }
    structured_result_program_brief = {
        "type": "string",
        "enum": list(result_programs),
        "description": (
            "Select the parser matching the registered artifact. Which "
            "selectors each program serves, per job type, is listed on "
            "inspect_run.program (with artifact_id) and on the capability "
            "receipt's coverage; the selected method/settings must still "
            "emit the quantity."
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
            "electronic bond-order assignment. The "
            "vibrational_mode_atom_participation selector is host-derived, "
            "not printed by any program: it is each atom's share of a mode's "
            "squared displacement, computed from the displacement vectors "
            "the program did print and renormalised so a row sums to one, "
            "which is what makes it comparable across programs. It says how "
            "much of a mode an atom carries, never what the motion is; and "
            "inside a degenerate set the individual vectors are an arbitrary "
            "basis, so read vibrational_mode_degeneracy_group before "
            "assigning motion to a single mode."
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
                "layer, not project additional_route_parameters. That "
                "escape hatch carries two things and nothing else: "
                "source-required scientific keywords that refine an "
                "otherwise supported typed method, and print directives "
                "that change no method and only make the program report "
                "more of what it already computed -- ORCA's Hirshfeld "
                "population analysis being the case that matters, because "
                "the condensed-Fukui literature asks for that scheme and "
                "ORCA prints it only when told to. What goes here is a "
                "route-line token, not a block setting: for ORCA these are "
                "the words that follow '!', so the Hirshfeld analysis is "
                "requested by the single token Hirshfeld. A block-style "
                "directive such as Print[P_Hirshfeld] 1 belongs inside a "
                "%output or %scf block, and ORCA rejects the entire input "
                "before any calculation starts when one appears on the "
                "route line. Every token you put here is displayed to the "
                "reviewer verbatim. It is not "
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
            "establish_project",
            (
                "Render, promote, and validate one project in a single call, "
                "returning all three receipts. Prefer this whenever a node "
                "needs a project: it is the same three steps in the same "
                "order, and doing them separately costs three turns per node "
                "for no additional evidence. Section rules are exactly those "
                "of project_yaml(action=render)."
            ),
            {
                "program": program,
                "sections": {
                    "type": "object",
                    "additionalProperties": {"type": "object"},
                },
                "artifact_id": _string(),
                "capability_receipt_sha256": digest,
            },
            (
                "program",
                "sections",
                "artifact_id",
                "capability_receipt_sha256",
            ),
        ),
        _tool(
            "fetch_pubchem_geometry",
            (
                "Bring in a molecule this workspace does not hold, by name, "
                "CID, or SMILES: you name the identifier, the host fetches "
                "the record and owns the bytes. Use it for a molecule the "
                "question needs and nobody supplied -- a reference computed "
                "at your own level so its systematic cancels, a calibration "
                "standard. What arrives is a depositor's conformer carrying "
                "that depositor's symmetry, not a relaxed structure, and it "
                "binds no electronic state: bind charge and multiplicity "
                "explicitly, and the consuming stage is a new workflow."
            ),
            {
                "artifact_id": _string(),
                "identifier": _string(),
            },
            ("artifact_id", "identifier"),
        ),
        _tool(
            "bind_scientific_identity",
            (
                "Bind an explicit charge and multiplicity to the exact host "
                "geometry artifact. input_artifact_id must identify a "
                "geometry_xyz artifact, never a project YAML or result. "
                "File names, XYZ comments, element lists, project settings, "
                "and preview artifacts do not establish molecular identity: "
                "measure the structure before binding what it is."
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
            "bind_scan_point_geometry",
            (
                "Carry one point of a completed relaxed scan forward as a "
                "geometry input. You choose the point -- read the surface "
                "first through scan_point_indices, scan_coordinate_values "
                "and scan_energies, then name the 1-based index of the "
                "structure you want; the host records which result and which "
                "point it came from, at what coordinate and energy. Using "
                "the returned geometry is a changed molecular input, so the "
                "stage that consumes it is a new workflow needing its own "
                "review. When the point you want is the surface's "
                "minimum-energy sample, you do not need this tool or a "
                "second workflow: declare the consumer's geometry input as a "
                "producer edge from the scan node, and one approval covers "
                "scan and consumer. artifact_id must identify a completed "
                "scan result already registered in this workspace."
            ),
            {
                "artifact_id": _string(),
                "point_index": {"type": "integer", "minimum": 1},
                "program": program,
            },
            ("artifact_id", "point_index"),
        ),
        _tool(
            "compose_molecular_arrangement",
            (
                "Place two identity-bound geometry artifacts into one "
                "arrangement at an explicit atomic contact. You choose the "
                "fragments, the 1-based contact atoms, and the contact "
                "distance in angstrom; the host owns the placement "
                "mathematics (the named pair held at the distance, every "
                "other interfragment pair kept outside covalent radii plus "
                "buffer, remaining freedom maximising separation) and "
                "writes the composed geometry into the workspace with full "
                "parent lineage. Composition never infers an electronic "
                "state: bind the arrangement's charge and multiplicity "
                "explicitly afterwards, and the stage that consumes it is "
                "a new workflow needing its own review. Both fragments must "
                "already carry a scientific identity."
            ),
            {
                "composed_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the composed "
                        "geometry artifact."
                    ),
                },
                "fragment_a_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact whose "
                        "coordinates are kept fixed; its atoms come first "
                        "in the composed file."
                    ),
                },
                "fragment_b_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact the host "
                        "places against fragment A."
                    ),
                },
                "fragment_a_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": ("1-based contact atom within fragment A."),
                },
                "fragment_b_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": ("1-based contact atom within fragment B."),
                },
                "distance_angstrom": {
                    "type": "number",
                    "minimum": 0.5,
                    "maximum": 10.0,
                    "description": (
                        "Contact distance between the two named atoms, in "
                        "angstrom."
                    ),
                },
                "fragment_b_atoms": {
                    "type": "array",
                    "items": {"type": "integer", "minimum": 1},
                    "minItems": 2,
                    "description": (
                        "1-based atoms of B all held at "
                        "distance_angstrom from fragment_a_atom, and "
                        "exempt from the clash floor for that pair "
                        "only. For a contact spread over several atoms "
                        "at once. Include fragment_b_atom. The set must "
                        "lie on a circle and be coplanar; the receipt "
                        "reports the distance achieved per atom."
                    ),
                },
                "fragment_a_atom_2": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "Optional second contact atom within fragment A. "
                        "Single-contact placement points the unpinned ends "
                        "apart (remaining freedom maximises separation), so "
                        "a doubly-contacted motif -- a cyclic "
                        "hydrogen-bonded dimer -- is asked for with a "
                        "second contact, never approximated by scanning a "
                        "contact closed. Give all three _2 fields together; "
                        "the host solves both distances simultaneously."
                    ),
                },
                "fragment_b_atom_2": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "Optional second contact atom within fragment B."
                    ),
                },
                "distance_angstrom_2": {
                    "type": "number",
                    "minimum": 0.5,
                    "maximum": 10.0,
                    "description": (
                        "Optional second contact distance, in angstrom."
                    ),
                },
            },
            (
                "composed_artifact_id",
                "fragment_a_artifact_id",
                "fragment_b_artifact_id",
                "fragment_a_atom",
                "fragment_b_atom",
                "distance_angstrom",
            ),
        ),
        _tool(
            "derive_molecular_species",
            (
                "Make a new geometry from one identity-bound parent by "
                "keeping an ordered subset of its atoms -- the operation "
                "underneath homolysis, deprotonation and pulling a fragment "
                "out of a larger structure. Give exactly one of "
                "removed_atoms (natural when one atom leaves a large "
                "molecule) or kept_atoms (natural when extracting a "
                "fragment, and it fixes the new atom order); the host records "
                "both either way. Coordinates are copied unchanged, so the "
                "result is a starting structure and not a relaxed one -- "
                "optimise it in the consuming stage. Derivation never infers "
                "an electronic state: removing a hydrogen gives a radical or "
                "an anion depending on where its electron went, so bind the "
                "new charge and multiplicity explicitly afterwards, and the "
                "stage that consumes this geometry is a new workflow needing "
                "its own review. The parent must already carry a scientific "
                "identity."
            ),
            {
                "derived_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the derived "
                        "geometry artifact."
                    ),
                },
                "parent_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact to derive from."
                    ),
                },
                "removed_atoms": {
                    "type": "array",
                    "items": {"type": "integer", "minimum": 1},
                    "minItems": 1,
                    "description": (
                        "1-based parent atoms to leave out. Mutually "
                        "exclusive with kept_atoms."
                    ),
                },
                "kept_atoms": {
                    "type": "array",
                    "items": {"type": "integer", "minimum": 1},
                    "minItems": 1,
                    "description": (
                        "1-based parent atoms to keep, in the order they "
                        "should appear in the new geometry. Mutually "
                        "exclusive with removed_atoms."
                    ),
                },
            },
            ("derived_artifact_id", "parent_artifact_id"),
        ),
        _tool(
            "break_symmetry",
            (
                "Perturb every atom of an identity-bound geometry by a "
                "seed you name, by at most the amplitude you name, so an "
                "exactly symmetric starting structure can relax away from "
                "the saddle its symmetry pins it to. The host draws the "
                "displacement deterministically from the seed, removes the "
                "net translation, records the largest step it actually "
                "took beside the one asked, and states the point-group "
                "estimate before and after. Atom count, order and formula "
                "are unchanged, so "
                "parent atom i is perturbed atom i. The result is a "
                "STARTING structure with no electronic state bound; bind "
                "charge and multiplicity afterwards, and let the "
                "optimisation that consumes it decide whether the step "
                "escaped the saddle -- an amplitude is never refused on "
                "merit. Refusals are structural: a non-integer seed, an "
                "amplitude outside (0, 0.5] Å, a parent with no bound "
                "identity. An amplitude of 0.02-0.05 Å is a perturbation; "
                "0.2 Å and above is a different starting guess."
            ),
            {
                "perturbed_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the perturbed "
                        "geometry artifact."
                    ),
                },
                "input_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact to perturb."
                    ),
                },
                "seed": {
                    "type": "integer",
                    "minimum": 0,
                    "description": (
                        "Integer seed of the displacement; the same seed "
                        "and amplitude on the same parent give the same "
                        "bytes, and a second seed is a second perturbation."
                    ),
                },
                "amplitude_angstrom": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "maximum": 0.5,
                    "description": (
                        "Largest displacement any atom may take, in "
                        "angstrom; the receipt records the largest it did."
                    ),
                },
            },
            (
                "perturbed_artifact_id",
                "input_artifact_id",
                "seed",
                "amplitude_angstrom",
            ),
        ),
        _tool(
            "edit_molecular_geometry",
            (
                "Set one internal coordinate of an identity-bound geometry "
                "to a value you choose. The host owns the arithmetic: it "
                "moves the side of the coordinate you name as one rigid "
                "piece, leaves every other atom exactly where it was, "
                "measures the coordinate before and after, and records which "
                "atoms moved. Use it to reach a structure a plain "
                "optimisation cannot -- another conformer, the far side of a "
                "rotational barrier, a deliberately stretched bond. Atom "
                "count, atom order and formula are unchanged, so parent atom "
                "i is edited atom i. Which side moves is your scientific "
                "choice and must be named. The result is a STARTING "
                "structure: the value you requested is what you asked for, "
                "and only an optimisation and its validation verdict say "
                "what the coordinate really is, so plan the consuming "
                "workflow to find out. Editing does not change or infer "
                "electronic state -- bind charge and multiplicity explicitly "
                "afterwards. Refusals are structural (an axis that is not a "
                "bond, a ring a rigid motion would tear, collinear atoms); "
                "for a coordinate inside a ring, use a constrained "
                "optimisation or a relaxed scan instead. The two atoms "
                "must already be bonded to each other in the perceived "
                "connectivity, because the axis is what tells the host "
                "which side to carry: setting a distance between two "
                "separate fragments -- a nucleophile and its substrate, "
                "the approach in a transition-state guess -- is not an "
                "edit at all, and compose_molecular_arrangement is the "
                "operation that places unbound pieces at a chosen "
                "contact and distance."
            ),
            {
                "edited_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the edited "
                        "geometry artifact."
                    ),
                },
                "input_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact to edit."
                    ),
                },
                "operation": {
                    "type": "string",
                    "enum": sorted(EDITABLE_COORDINATE_OPERATIONS),
                    "description": (
                        "Which internal coordinate to set: set_bond_length "
                        "takes 2 atoms, set_angle 3, set_dihedral 4 -- the "
                        "same coordinates a scan drives and a modred holds."
                    ),
                },
                "atoms": {
                    "type": "array",
                    "items": {"type": "integer", "minimum": 1},
                    "minItems": 2,
                    "maxItems": 4,
                    "description": (
                        "1-based atoms defining the coordinate, in bonded "
                        "order: i-j for a bond, i-j-k for an angle with the "
                        "vertex in the middle, i-j-k-l for a torsion about "
                        "the j-k bond."
                    ),
                },
                "moving_side_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "Which of the coordinate's own atoms sits on the "
                        "side that moves. Both choices reach the same "
                        "measured value and give different molecules, so "
                        "there is no default: name the side whose motion is "
                        "the chemistry you mean."
                    ),
                },
                "target_value": {
                    "type": "number",
                    "description": (
                        "The value to set, in the coordinate's own unit: "
                        "angstrom for set_bond_length, degrees for set_angle "
                        "(strictly between 0 and 180) and set_dihedral "
                        "(-180 to 180, signed by the IUPAC convention)."
                    ),
                },
            },
            (
                "edited_artifact_id",
                "input_artifact_id",
                "operation",
                "atoms",
                "moving_side_atom",
                "target_value",
            ),
        ),
        _tool(
            "bind_reached_geometry",
            (
                "Carry the structure a run actually reached into a new "
                "geometry input. This is the ordinary answer to an "
                "optimisation that ran out of iterations or wall time: a "
                "fresh start from the original coordinates repeats the "
                "same path, while the reached structure is many steps "
                "further down the surface. The host reads the geometry "
                "through the same selector plane every quantity uses, "
                "writes the bytes itself, and records which result it "
                "came from, what ending this workspace recorded for that "
                "run, and whether the program terminated normally. It "
                "binds no electronic state -- what the structure is "
                "depends on the question you ask next -- so bind charge "
                "and multiplicity explicitly afterwards, and the stage "
                "that optimises it is a new workflow needing its own "
                "review. Nothing here upgrades the source: a node that "
                "failed still satisfies no producer edge, and this "
                "geometry is a starting structure, graded by the "
                "optimisation that consumes it."
            ),
            {
                "artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "The result whose reached geometry you want, "
                        "already registered in this workspace."
                    ),
                },
                "reached_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "The id to give the new geometry artifact."
                    ),
                },
                "program": {
                    **_string(),
                    "description": (
                        "The program that wrote the result; its geometry "
                        "is read with that program's reader."
                    ),
                },
            },
            ("artifact_id", "reached_artifact_id", "program"),
        ),
        _tool(
            "characterise_stationary_point",
            (
                "State what kind of stationary point a completed result "
                "actually is, and have the host check that statement "
                "against the frequencies the program itself printed. Reach "
                "for this when a search converged onto something other "
                "than what it promised and you judge the structure it "
                "found to be worth reporting: a minimum search that lands "
                "on a first-order saddle has found a transition state, and "
                "its energy and thermochemistry are readable like any "
                "other result's. The node keeps its own verdict -- the "
                "promise it was launched under was still not met -- and "
                "this receipt stands beside that, so a number you deliver "
                "from those bytes can say which structure it belongs to. "
                "The host refuses an order the printed modes do not "
                "support, counting a mode as imaginary below -20 cm^-1, "
                "the same convention the validator uses. What the "
                "structure means is yours to say."
            ),
            {
                "result_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "The completed, frequency-bearing result to "
                        "characterise."
                    ),
                },
                "program": {
                    **_string(),
                    "description": (
                        "The program that wrote the result; its "
                        "frequencies are read with that program's reader."
                    ),
                },
                "order_claimed": {
                    "type": "integer",
                    "minimum": 0,
                    "description": (
                        "The number of imaginary modes you say this "
                        "structure has: 0 for a minimum, 1 for a "
                        "transition state, more for a higher-order saddle."
                    ),
                },
                "node_id": {
                    **_string(),
                    "description": (
                        "Optional: the workflow node whose run produced "
                        "this result, when you know it."
                    ),
                },
                "anomaly_receipt_sha256": {
                    **_string(),
                    "description": (
                        "Optional: the anomaly observation this answers, "
                        "by its receipt digest."
                    ),
                },
            },
            ("result_artifact_id", "program", "order_claimed"),
        ),
        _tool(
            "displace_along_vibrational_mode",
            (
                "Step a completed result's geometry along one of its own "
                "printed normal modes, producing a new starting structure. "
                "This is the move to reach for when an optimisation "
                "converges onto a saddle rather than a minimum, or when a "
                "transition-state search returns the wrong number of "
                "imaginary modes: read which mode is wrong, step along it, "
                "and relax again. The displacement vectors are the ones the "
                "program itself printed and the host owns the arithmetic; "
                "you own which mode and how far. Mode 1 is the lowest "
                "printed mode, which is the imaginary one on a saddle. Atom "
                "count and order are preserved, so parent atom i is "
                "displaced atom i. Displacing never infers an electronic "
                "state, so bind charge and multiplicity explicitly "
                "afterwards, and the optimisation that consumes the "
                "structure is a new workflow needing its own review -- "
                "whether the step escaped the saddle is decided there, by "
                "physics, not here."
            ),
            {
                "displaced_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Identifier for the displaced geometry this "
                        "produces."
                    ),
                },
                "result_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "The completed, frequency-bearing result whose "
                        "printed modes are used. A result with no "
                        "frequencies carries no modes to step along."
                    ),
                },
                "program": {
                    **_string(),
                    "description": (
                        "The program that wrote the result; its modes are "
                        "read with that program's reader."
                    ),
                },
                "mode_index": {
                    "type": "integer",
                    "description": (
                        "Which printed mode to step along, 1-based. Mode 1 "
                        "is the lowest, which is the imaginary mode on a "
                        "saddle."
                    ),
                },
                "amplitude_angstrom": {
                    "type": "number",
                    "description": (
                        "Largest atomic displacement, in angstrom, after "
                        "the mode is scaled; the sign chooses the "
                        "direction, so +0.3 and -0.3 are the two sides of "
                        "a saddle along its imaginary mode. A step is a "
                        "starting guess: the consuming optimisation grades "
                        "it, and nothing here refuses a value on "
                        "scientific merit."
                    ),
                },
            },
            (
                "displaced_artifact_id",
                "result_artifact_id",
                "program",
                "mode_index",
                "amplitude_angstrom",
            ),
        ),
        _tool(
            "append_molecular_atom",
            (
                "Add one atom to an identity-bound parent, placed by the "
                "three internal coordinates that say where it sits: a bond "
                "length to the atom it attaches to, an angle to a second "
                "atom, and a dihedral against a third. This is derivation's "
                "mirror -- protonation, hydrogenation, capping a radical, "
                "deuteration -- and the host owns the placement arithmetic "
                "and the bytes. Parent atom indices are unchanged and the "
                "appended atom is last. Appending never infers an electronic "
                "state: adding a hydrogen gives a cation or a radical "
                "depending on whether it brought an electron, so bind charge "
                "and multiplicity explicitly afterwards. The result is a "
                "starting structure and the consuming stage is a new "
                "workflow needing its own review."
            ),
            {
                "appended_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the appended "
                        "geometry artifact."
                    ),
                },
                "input_artifact_id": {
                    **_string(),
                    "description": (
                        "Identity-bound geometry_xyz artifact to add to."
                    ),
                },
                "element": {
                    **_string(),
                    "description": (
                        "Element symbol of the atom to add, e.g. 'H'."
                    ),
                },
                "anchor_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "1-based parent atom the new atom bonds to; the "
                        "bond length is measured to this atom."
                    ),
                },
                "angle_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "1-based parent atom the angle new-anchor-angle is "
                        "measured against."
                    ),
                },
                "dihedral_atom": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "1-based parent atom the torsion "
                        "new-anchor-angle-dihedral is measured against; it "
                        "must not be collinear with the other two."
                    ),
                },
                "bond_length_angstrom": {
                    "type": "number",
                    "minimum": 0.5,
                    "maximum": 5.0,
                    "description": (
                        "Distance from the anchor atom, in angstrom."
                    ),
                },
                "angle_degrees": {
                    "type": "number",
                    "exclusiveMinimum": 0.0,
                    "exclusiveMaximum": 180.0,
                    "description": ("Angle at the anchor atom, in degrees."),
                },
                "dihedral_degrees": {
                    "type": "number",
                    "description": (
                        "Torsion about the anchor-angle bond, in degrees, "
                        "signed by the IUPAC convention."
                    ),
                },
            },
            (
                "appended_artifact_id",
                "input_artifact_id",
                "element",
                "anchor_atom",
                "angle_atom",
                "dihedral_atom",
                "bond_length_angstrom",
                "angle_degrees",
                "dihedral_degrees",
            ),
        ),
        _tool(
            "inspect_database_records",
            (
                "Enumerate the records of a workspace chemsmart .db "
                "artifact: record ids, formulas, structure counts, and the "
                "record's own stored fields (charge, multiplicity, energy, "
                "optimized flag). Stored fields are observations from the "
                "database's provenance, never identity bindings -- a record "
                "may store no electronic state at all. An optional query "
                "filters records with the database query language (FIELD "
                "OPERATOR VALUE joined by AND/OR, ~ for contains). To use a "
                "record's geometry, admit it with "
                "extract_database_record_geometry and then bind charge and "
                "multiplicity explicitly with bind_scientific_identity."
            ),
            {
                "database_artifact_id": {
                    **_string(),
                    "description": (
                        "Workspace chemsmart_db artifact to enumerate."
                    ),
                },
                "query": {
                    **_string(),
                    "description": (
                        "Optional record filter, e.g. "
                        "\"program = 'xtb' AND normal_termination = 1\"; "
                        "an invalid query is refused naming the supported "
                        "fields and operators."
                    ),
                },
                "limit": {
                    "type": "integer",
                    "minimum": 1,
                    "maximum": 200,
                    "description": (
                        "Maximum records to return (default 50); the total "
                        "count is always reported."
                    ),
                },
            },
            ("database_artifact_id",),
        ),
        _tool(
            "extract_database_record_geometry",
            (
                "Copy one database record's stored coordinates into a new "
                "host-owned workspace geometry artifact with full lineage "
                "(database digest, record id, structure selection). "
                "Coordinates are copied unchanged, and execution never "
                "reads the .db again -- the extracted artifact is an "
                "ordinary geometry input. It carries NO electronic state: "
                "the record's stored charge and multiplicity are returned "
                "as observations, and binding is yours -- call "
                "bind_scientific_identity explicitly afterwards, stating "
                "your own charge and multiplicity for the calculation you "
                "intend. Give exactly one of record_index or record_id; "
                "structure_index is required when the record stores more "
                "than one structure."
            ),
            {
                "extracted_artifact_id": {
                    **_public_identifier(),
                    "description": (
                        "Workspace-unique identifier for the extracted "
                        "geometry artifact."
                    ),
                },
                "database_artifact_id": {
                    **_string(),
                    "description": (
                        "Workspace chemsmart_db artifact to extract from."
                    ),
                },
                "record_index": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "1-based record index. Mutually exclusive with "
                        "record_id."
                    ),
                },
                "record_id": {
                    **_string(),
                    "description": (
                        "Record id or unique prefix. Mutually exclusive "
                        "with record_index."
                    ),
                },
                "structure_index": {
                    "type": "integer",
                    "minimum": 1,
                    "description": (
                        "1-based structure within the record, in the "
                        "record's own order; required when the record "
                        "stores several structures."
                    ),
                },
            },
            ("extracted_artifact_id", "database_artifact_id"),
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
            "open_guide",
            (
                "Open one guide named in the system prompt -- a family unit "
                "of tools, operations and guidance (structure, scan, "
                "constants, cbs, ensemble, spectroscopy, database, "
                "crossprogram, recovery, saddle) -- or one advisory "
                "domain-knowledge skill. Returns the text and, for a guide, "
                "the tools and operations that join the surface on the next "
                "turn. Guidance only: it never establishes readiness, "
                "approval, terminal state, or an accuracy claim, and never "
                "substitutes for a typed host receipt."
            ),
            {"guide_id": _public_identifier()},
            ("guide_id",),
        ),
        _tool(
            "declare_requested_observable",
            (
                "Restate the task's requested observables as your first "
                "typed act, before planning: an identifier, the unit the "
                "answer will be reported in, and one sentence of meaning "
                "each. The host verifies the unit parses; the completion "
                "gate later requires, for every declared observable, a "
                "delivered claim carrying its id -- as the claim's "
                "claim_id, or as the quantity id of the receipt it stands "
                "on -- in the declared dimension; values are never "
                "checked -- and an "
                "undelivered declared observable is named in the "
                "completion receipt like a plan output the chain could "
                "not fulfil. A declaration cannot be rebound to a "
                "different unit; later calls may add observables. You "
                "may also record the sign you expect and what it rests "
                "on; the host prints it beside the delivered number so "
                "a reader sees both, and grades neither. To print them "
                "beside each other the host must know which claim "
                "answers which observable, and it joins them on the "
                "claim's ``claim_id`` or on the receipt quantity id the "
                "claim stands on: give the claim that answers an "
                "observable that observable's own id as its claim_id "
                "(in a planned claim node, as the input_id that binds "
                "it). A dimension is never used to guess: two energies "
                "or three potentials -- the ordinary shape of a "
                "comparison -- cannot be told apart by it, and the row "
                "reports the expectation with no number beside it "
                "rather than guess."
            ),
            {
                "observables": {
                    "type": "array",
                    "minItems": 1,
                    "description": (
                        "The observables the task asks for, one entry " "each."
                    ),
                    "items": {
                        "type": "object",
                        "properties": {
                            "observable_id": {
                                **_public_identifier(),
                                "description": (
                                    "Stable identifier for this "
                                    "observable; it cannot be rebound "
                                    "to a different unit later."
                                ),
                            },
                            "unit": {
                                **_string(),
                                "description": (
                                    "Unit the answer will be reported "
                                    "in, from the typed unit vocabulary "
                                    "(e.g. 'kcal/mol', 'eV', "
                                    "'angstrom', '1' for a count)."
                                ),
                            },
                            "meaning": {
                                **_string(),
                                "description": (
                                    "One sentence saying what this "
                                    "observable is, in scientific "
                                    "terms."
                                ),
                            },
                            "expected_sign": {
                                "type": "string",
                                "enum": ["positive", "negative", ""],
                                "description": (
                                    "Optional. The sign you expect the "
                                    "delivered value to carry, recorded "
                                    "before the evidence exists. The "
                                    "host restates it beside the "
                                    "delivered number and nothing else: "
                                    "a prediction is displayed, never "
                                    "scored, and a diverging one settles "
                                    "nothing -- being wrong about the "
                                    "chemistry is a result, not a "
                                    "defect."
                                ),
                            },
                            "expected_low": {
                                "type": "number",
                                "description": (
                                    "Optional, with expected_high: the "
                                    "range you expect the delivered "
                                    "magnitude to fall in, in this "
                                    "observable's own unit. Direction "
                                    "is the easy half of a prediction "
                                    "and scale is the half that gets "
                                    "away, so a range is worth "
                                    "recording even when the sign feels "
                                    "certain. Displayed, never scored."
                                ),
                            },
                            "expected_high": {
                                "type": "number",
                                "description": (
                                    "The upper end of that range, in "
                                    "the same unit; both ends are given "
                                    "together or neither is. A point "
                                    "expectation -- exactly one imaginary "
                                    "mode, a yes-or-no verdict -- is "
                                    "written with both ends equal."
                                ),
                            },
                            "expectation_basis": {
                                **_string(),
                                "description": (
                                    "Required when a sign or a range is "
                                    "given: what the expectation rests "
                                    "on (a mechanism, a trend, a "
                                    "measured fact). An expectation "
                                    "without a reason is a coin flip "
                                    "and is refused."
                                ),
                            },
                            "supersedes_observable_id": {
                                **_public_identifier(),
                                "description": (
                                    "Optional: an earlier declaration of "
                                    "yours whose unit was wrong, retired by "
                                    "this one. Both stay on the record and "
                                    "the retired id stops being owed."
                                ),
                            },
                            "role": {
                                "type": "string",
                                "enum": ["requested", "diagnostic"],
                                "description": (
                                    "Optional, default requested. A "
                                    "diagnostic is your own prediction "
                                    "about the route -- which stationary "
                                    "point a search reaches, which spin "
                                    "state lies lower -- scored like any "
                                    "expectation and never owed: "
                                    "undelivered is no limitation, "
                                    "diverged is an observation the word "
                                    "carries."
                                ),
                            },
                            "failure_update_rule": {
                                **_string(),
                                "description": (
                                    "Required with role diagnostic: what "
                                    "you change about the route if this "
                                    "is falsified. A prediction whose "
                                    "failure changes nothing is not a "
                                    "diagnostic."
                                ),
                            },
                            "required_tolerance": {
                                "type": "number",
                                "minimum": 0,
                                "description": (
                                    "The error the task says it can "
                                    "tolerate, in this unit -- 'good to "
                                    "about 0.2 V' is 0.2. Restate it here "
                                    "when the task states one; the host "
                                    "compares it with the uncertainty your "
                                    "claim states, and never with the "
                                    "value. It is frozen at this "
                                    "declaration."
                                ),
                            },
                            "tolerance_basis": {
                                **_string(),
                                "description": (
                                    "Required with required_tolerance: the "
                                    "task's own words that fix it, or that "
                                    "the task fixes none and this is your "
                                    "reading of what the answer must "
                                    "support."
                                ),
                            },
                            "tolerance_origin": {
                                "type": "string",
                                "enum": ["task", "session"],
                                "description": (
                                    "State it with required_tolerance. "
                                    "'task' means you restated a precision "
                                    "the task states; 'session' means the "
                                    "task fixes none for this quantity and "
                                    "the tolerance is your own reading of "
                                    "what the decision needs. Both are "
                                    "legitimate -- declaring the margin a "
                                    "decision actually turns on is a route "
                                    "the host offers -- and neither is "
                                    "refused, and omitting it is recorded "
                                    "as 'unstated' rather than refused. It "
                                    "is asked because a "
                                    "reader cannot otherwise tell a "
                                    "requirement the requester set from "
                                    "one you set, and a session tolerance "
                                    "chosen after the number is known is "
                                    "a different thing from a restated "
                                    "one."
                                ),
                            },
                            "method_resolution": {
                                "type": "number",
                                "description": (
                                    "Optional. A deadband on the delivered "
                                    "magnitude: a value smaller than this "
                                    "prints indeterminate instead of "
                                    "agreeing or diverging. It is not your "
                                    "answer's uncertainty -- that goes on "
                                    "the claim."
                                ),
                            },
                        },
                        "required": ["observable_id", "unit", "meaning"],
                        "additionalProperties": False,
                    },
                },
            },
            ("observables",),
        ),
        _tool(
            "plan_scientific_workflow",
            (
                "Plan one connected scientific tool chain containing any "
                "required program calculations and deterministic analysis "
                "stages. For an analysis-only task over registered results, "
                "calculation_nodes may be empty; do not invent a documentary "
                "or blocked calculation placeholder. For a calculation-only "
                "task, analysis_nodes and required_output_ids may be empty. "
                "Bind scientific identity to every initial geometry first; "
                "every calculation node needs at least one expected output, "
                "and future producer outputs remain unresolved. A null "
                "scientific_workflow_plan means the binding must be repaired "
                "and this tool called again. "
                "Future analysis inputs name producer node/output pairs; they "
                "do not require artifact or receipt hashes before execution. "
                "A result-extraction or thermochemistry root may instead "
                "consume one existing host-registered result by artifact_id. "
                "Keep unsupported requested analyses as blocked_unsupported nodes."
            ),
            {
                "plan_id": _public_identifier(),
                "workflow_id": _public_identifier(spelling_rule=True),
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
                    "minItems": 0,
                    "maxItems": 128,
                    "items": _analysis_intent_node_schema(
                        operations=operations
                    ),
                },
                "required_output_ids": {
                    "type": "array",
                    "minItems": 0,
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
                "support_repairs": {
                    "type": "array",
                    "description": (
                        "Declare named calculation stages non-executable "
                        "scientific intent, with the reason -- for example a "
                        "functional the program's validator refused. One "
                        "direction only: a declared stage stays displayed "
                        "with the workflow, is excluded from approval, and "
                        "is never launched, so this can only narrow what "
                        "runs. Reversing it is a new workflow."
                    ),
                    "maxItems": 64,
                    "items": {
                        "type": "object",
                        "properties": {
                            "node_id": _public_identifier(),
                            "blocked_reason": _string(),
                        },
                        "required": ["node_id", "blocked_reason"],
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
                                        "uncertainty_producer_output_id": {
                                            **_string(),
                                            "description": (
                                                "On a claim: the analysis "
                                                "output that computes this "
                                                "number's uncertainty. You "
                                                "cannot know an estimate "
                                                "before the results exist, "
                                                "but you can plan its "
                                                "estimator; the host "
                                                "evaluates it after "
                                                "execution and the claim "
                                                "carries it as measured, "
                                                "cited to that output. "
                                                "Without one the delivery "
                                                "states no uncertainty and "
                                                "the goal reopens to assess "
                                                "it."
                                            ),
                                        },
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
                "previously planned scientific workflow. The response also "
                "names, as artifacts_without_structural_read, every "
                "registered geometry or result artifact this session has not "
                "yet read with a structural selector (positions or "
                "connectivity) -- an id on that list is a measurement not "
                "yet made, not a verdict."
            ),
            {"workflow_id": _public_identifier()},
            ("workflow_id",),
        ),
        _tool(
            "select_execution_wave",
            (
                "Name the calculations to run together as one wave, in "
                "the order you want them."
            ),
            {
                "workflow_id": _public_identifier(),
                "node_ids": {
                    "type": "array",
                    "minItems": 1,
                    "items": _public_identifier(),
                },
            },
            ("workflow_id", "node_ids"),
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
                # Optional, and only meaningful for a family that drives or
                # holds a coordinate. The executor rebuilds an approved node
                # through this tool, so the coordinate has to be expressible
                # here or it cannot survive into the launched command.
                "internal_coordinates": _internal_coordinates_schema(),
                # Optional: the digest of a host-recorded anomaly this node
                # investigates. Charged to the envelope's excursion line, never
                # to the engine-call budget; it may feed no untagged node, so the
                # asked observable stays owed to the budget.
                "excursion": digest,
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
            (
                "Record public chemical rationale, alternatives, and "
                "uncertainty. A doubt that cites its receipt binds to the "
                "claim standing on that receipt: add the evidence reference "
                "doubt:{receipt_sha256} and the completion gate will not "
                "certify past it (passed becomes partial naming the doubted "
                "quantity, and the goal returns to the human). A doubt kept "
                "in prose alone binds to nothing."
            ),
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
                "unreachable_observable_ids": {
                    "type": "array",
                    "maxItems": 32,
                    "description": (
                        "The typed refusal: each declared observable this "
                        "goal cannot deliver from admissible evidence, with "
                        "the producer it would need and the receipts that "
                        "show the gap. The host verifies what it can -- a "
                        "selector no program in the envelope declares for "
                        "the job type, or a blocked_unsupported analysis "
                        "node in this session's plan whose output_id is the "
                        "observable -- and only a verified refusal settles "
                        "the goal unreachable_from_evidence; an unverified "
                        "one returns to the human naming it. Never a "
                        "shortcut past computing what can be computed."
                    ),
                    "items": {
                        "type": "object",
                        "properties": {
                            "observable_id": _public_identifier(),
                            "statement": _string(),
                            "selector": _public_identifier(),
                            "jobtype": _public_identifier(),
                            "blocked_node_id": _public_identifier(),
                            "receipt_sha256s": {
                                "type": "array",
                                "items": digest,
                                "minItems": 1,
                                "maxItems": 16,
                            },
                        },
                        "required": [
                            "observable_id",
                            "statement",
                            "receipt_sha256s",
                        ],
                        "additionalProperties": False,
                    },
                },
                "menu_route_dispositions": {
                    "type": "array",
                    "description": (
                        "What you did with each route the wake's "
                        "repair_menu offered -- taken, rejected or deferred "
                        "-- with the mechanism and its receipts. The host "
                        "verifies the route was offered and the receipts "
                        "are its own, never the choice; the next wake shows "
                        "them."
                    ),
                    "items": {
                        "type": "object",
                        "properties": {
                            "route": _public_identifier(),
                            "disposition": {
                                "type": "string",
                                "enum": ["taken", "rejected", "deferred"],
                            },
                            "reason": _string(),
                            "receipt_sha256s": {
                                "type": "array",
                                "items": digest,
                            },
                        },
                        "required": ["route", "disposition", "reason"],
                        "additionalProperties": False,
                    },
                },
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
            "inspect_result_selectors",
            (
                "List the selectors one completed result actually resolves, "
                "by probing its parser. A declaration describes a job type; "
                "what a particular result carries is narrower, because the "
                "method and settings decide what the program printed -- a "
                "spin-restricted run prints no <S^2>, a single point has no "
                "Hessian. Use it when planning analysis over results that "
                "already exist, so the shape of an artifact is learned once "
                "rather than one refused selector at a time. It cannot help "
                "before a calculation has run: there is no artifact to probe "
                "until the engine has produced one."
            ),
            {
                "program": structured_result_program,
                "artifact_id": _string(),
            },
            ("program", "artifact_id"),
        ),
        _tool(
            "inspect_run_outcome",
            (
                "Read how a recorded run in this workspace ended, as one "
                "typed vocabulary: per-node terminal states (validated, "
                "failed_nonconverged_scan_step, timeout_terminated, "
                "interrupted_mid_engine, ...) with the program-native "
                "findings, the engine's own redacted words, convergence "
                "and scan reached-versus-planned facts, wall seconds, and "
                "the evidence digests a revision cites. Call it with no "
                "run reference to list the runs this workspace records; "
                "name one to read its full outcome. Read the outcome "
                "before planning a revision of failed work: a plan that "
                "answers a failure it never read is answering a guess."
            ),
            {
                "run": _string(),
            },
            (),
        ),
        _tool(
            "extract_result_quantities",
            (
                "Parse selected numerical or scientific fields from a trusted "
                "host-bound result artifact. The model supplies semantic selectors, "
                "never a file path."
            ),
            {
                "program": structured_result_program_brief,
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
                "quasi-harmonic thermochemistry from a trusted frequency "
                "result. A supplied concentration "
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
                "quantity. The receipt also carries near_zero_mode_count: "
                "modes within 20 cm-1 of zero are treated as numerical "
                "noise, not evidence of a saddle point -- the documented "
                "convention and the default acceptance criterion when the "
                "task states none; an explicit criterion in the task text "
                "always overrides it."
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
                        "Use natural-abundance weighted isotope masses; "
                        "omitted means most-abundant masses."
                    ),
                },
                "frequency_scale_factor": {
                    "type": "number",
                    "exclusiveMinimum": 0,
                    "description": (
                        "Positive scale applied to every vibrational "
                        "frequency before thermochemistry; omitted is 1.0."
                    ),
                },
                "reaction_coordinate_mode": {
                    "type": "integer",
                    "minimum": 0,
                    "description": (
                        "1-based reaction-coordinate mode of a "
                        "characterised saddle; omit for a minimum."
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
                    "items": _quantity_expression_node_schema(
                        operations=operations
                    ),
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
                "host copies and converts the values. Where a claim "
                "answers an observable you declared, set its "
                "``claim_id`` to that observable's id -- the "
                "expectation you recorded is printed beside the "
                "delivered number only when those two meet, and "
                "``quantity_id`` names the receipt quantity, not the "
                "observable."
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
                            "uncertainty": {
                                "type": "number",
                                "minimum": 0,
                                "description": (
                                    "The error you attribute to this "
                                    "number, in its display unit: every "
                                    "term you would state as a limitation, "
                                    "together. State it whenever the "
                                    "observable declared a "
                                    "required_tolerance -- silence is not "
                                    "sufficiency, and the host reads the "
                                    "two together."
                                ),
                            },
                            "uncertainty_basis": {
                                "type": "string",
                                "enum": [
                                    "measured",
                                    "inferred",
                                    "asserted",
                                ],
                                "description": (
                                    "Required with uncertainty. measured "
                                    "or inferred name evidence and need "
                                    "uncertainty_reference; asserted is "
                                    "your judgement and takes none. An "
                                    "assertion is never penalised and is "
                                    "often the honest number, but it does "
                                    "not on its own discharge a tolerance "
                                    "the task set -- it routes you to "
                                    "measure it, to show the decision does "
                                    "not turn on it, or to say it cannot "
                                    "be established here."
                                ),
                            },
                            "uncertainty_reference": {
                                **_string(),
                                "description": (
                                    "Required with basis measured or "
                                    "inferred. For measured, "
                                    "'<receipt_sha256>:<quantity_id>': the "
                                    "host reads that quantity and checks "
                                    "your number against it, as it copies "
                                    "the value you claim. For inferred, the "
                                    "receipt or registered constant your "
                                    "judgement rests on. Only a magnitude "
                                    "the host can check discharges a "
                                    "tolerance. The host follows an "
                                    "expression to its roots and reports "
                                    "what it finds -- a coefficient of "
                                    "yours in the chain, a spread over one "
                                    "receipt, a magnitude of zero -- "
                                    "beside the assessment, and refuses "
                                    "none of them: whether such a number "
                                    "is the right uncertainty for this "
                                    "claim is yours to argue and the "
                                    "reader's to judge."
                                ),
                            },
                            "approximates": {
                                "type": "object",
                                "additionalProperties": False,
                                "properties": {
                                    "observable_id": _string(),
                                    "relationship": _string(),
                                    "basis": _string(),
                                },
                                "description": (
                                    "State this when your number "
                                    "approximates a declared observable "
                                    "rather than being it -- a 0 K gap "
                                    "for a Gibbs one. Name the id and "
                                    "the relationship, in your words. "
                                    "Never refused, still delivered; the "
                                    "expectation row stops reading as an "
                                    "exact match."
                                ),
                            },
                            "uncertainty_combination": {
                                "type": "object",
                                "additionalProperties": False,
                                "properties": {
                                    "rule": _string(),
                                    "input_meaning": _string(),
                                    "coverage": _string(),
                                    "dependence": _string(),
                                },
                                "description": (
                                    "How you combined your components: "
                                    "the rule, what its inputs mean, the "
                                    "coverage the total claims, what you "
                                    "assume about dependence. Copied and "
                                    "never graded, but it rides the "
                                    "verdict: the same components can "
                                    "fall inside or outside one tolerance "
                                    "on this choice alone."
                                ),
                            },
                            "uncertainty_components": {
                                "type": "array",
                                "maxItems": 16,
                                "description": (
                                    "Optional: the terms your uncertainty "
                                    "is made of, in your own words. The "
                                    "host reads their provenance, never "
                                    "their meaning, and never combines "
                                    "them: how they add is your judgement "
                                    "and the total is the one you state. "
                                    "Name a term you cannot quantify "
                                    "rather than omitting it -- it holds "
                                    "the requirement open, which is what "
                                    "such a term means."
                                ),
                                "items": {
                                    "type": "object",
                                    "properties": {
                                        "meaning": {
                                            **_string(),
                                            "description": (
                                                "What this term is, in one "
                                                "sentence of your own."
                                            ),
                                        },
                                        "magnitude": {
                                            "type": "number",
                                            "minimum": 0,
                                            "description": (
                                                "Its size in the display "
                                                "unit; omit when "
                                                "unquantified."
                                            ),
                                        },
                                        "basis": {
                                            "type": "string",
                                            "enum": [
                                                "measured",
                                                "inferred",
                                                "asserted",
                                                "unquantified",
                                            ],
                                        },
                                        "reference": {
                                            **_string(),
                                            "description": (
                                                "Required with measured or "
                                                "inferred: the receipt, or "
                                                "'<receipt>:<quantity_id>' "
                                                "to name the number itself. "
                                                "The host resolves what you "
                                                "name and never grades the "
                                                "magnitude: how the terms "
                                                "add, and at what coverage "
                                                "each is stated, is yours."
                                            ),
                                        },
                                    },
                                    "required": ["meaning", "basis"],
                                    "additionalProperties": False,
                                },
                            },
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
    return tools


#: Legacy tool names each merged planning tool replaces, in surface order.
MERGED_PLANNING_TOOLS: dict[str, tuple[str, ...]] = {
    "inspect_program": (
        "inspect_program_capability",
        "inspect_program_environment",
        "assess_program_candidate",
    ),
    "project_yaml": (
        "render_project_yaml",
        "promote_project_yaml",
        "establish_project",
        "read_project_yaml",
        "validate_project_yaml",
    ),
    "compile_command": (
        "prepare_program_node",
        "synthesize_command",
        "preview_command",
        "preflight_program_node",
    ),
    "inspect_run": ("inspect_run_outcome", "inspect_result_selectors"),
}

PROJECT_YAML_ACTIONS = ("establish", "render", "promote", "read", "validate")

#: What each project_yaml action needs, beyond ``action`` itself.
PROJECT_YAML_ACTION_ARGUMENTS: dict[str, tuple[str, ...]] = {
    "establish": (
        "program",
        "sections",
        "artifact_id",
        "capability_receipt_sha256",
    ),
    "render": ("program", "sections"),
    "promote": ("render_receipt_sha256", "artifact_id"),
    "read": ("program", "project_artifact_id"),
    "validate": ("project_artifact_id", "capability_receipt_sha256"),
}


def _merge_planning_tools(tools: tuple[dict, ...]) -> tuple[dict, ...]:
    """The model-facing surface: four gate-feeding groups become one tool
    each. Eleven of the withdrawn tools existed mainly so the model could
    carry a receipt digest from one call to the next; the merged tool does
    the join and returns every receipt it produced, so a red preview costs
    one turn instead of four.
    """

    by_name = {item["function"]["name"]: item for item in tools}

    def parameters(name: str) -> dict:
        return dict(by_name[name]["function"]["parameters"])

    def properties(name: str) -> dict:
        return dict(parameters(name).get("properties") or {})

    capability = by_name["inspect_program_capability"]["function"]
    inspect_program = _tool(
        "inspect_program",
        (
            "Inspect one program, job type and engine in a single call: "
            "the capability query (what the registry, the live CLI and the "
            "conformance overlay say, with the coverage cell naming which "
            "typed axes and validity rules the host reaches) and the "
            "environment query (which installed program and engine bind "
            "here). Returns the capability receipt, the environment receipt "
            "and the program and engine bindings; the capability receipt "
            "digest is what project_yaml(validate/establish) takes. "
            + str(capability.get("description") or "")
        ),
        properties("inspect_program_capability"),
        tuple(parameters("inspect_program_capability").get("required", ())),
    )
    project_properties: dict = {
        "action": {
            "type": "string",
            "enum": list(PROJECT_YAML_ACTIONS),
            "description": (
                "establish = render + promote + validate in one turn, the "
                "ordinary route for a new node; render/promote/validate are "
                "the same three steps one at a time, for repairing one; "
                "read returns an existing project artifact's document. "
                + "; ".join(
                    f"{action} needs {', '.join(names)}"
                    for action, names in PROJECT_YAML_ACTION_ARGUMENTS.items()
                )
                + "."
            ),
        }
    }
    for legacy in MERGED_PLANNING_TOOLS["project_yaml"]:
        for key, schema in properties(legacy).items():
            project_properties.setdefault(key, schema)
    project_yaml = _tool(
        "project_yaml",
        (
            "The one channel for program settings: a project YAML document "
            "per node, rendered from typed sections, promoted into the "
            "workspace as an artifact, and validated against the capability "
            "receipt for its program and job type. "
            + str(
                by_name["establish_project"]["function"].get("description")
                or ""
            )
        ),
        project_properties,
        ("action",),
    )
    prepare = by_name["prepare_program_node"]["function"]
    compile_command = _tool(
        "compile_command",
        (
            "Compile one planned calculation node: the host joins the "
            "recorded capability, environment, project, identity and "
            "artifact records, compiles argv through the live Click schema, "
            "runs the safe preview, and preflights the node, returning every "
            "receipt in one reply. A red finding names the field that needs "
            "repair; a node whose producer has not run yet returns "
            "waiting_for_artifact with the producer named. "
            + str(prepare.get("description") or "")
        ),
        properties("prepare_program_node"),
        tuple(parameters("prepare_program_node").get("required", ())),
    )
    selectors = by_name["inspect_result_selectors"]["function"]
    outcome = by_name["inspect_run_outcome"]["function"]
    inspect_run_properties = {
        **properties("inspect_run_outcome"),
        **properties("inspect_result_selectors"),
    }
    inspect_run = _tool(
        "inspect_run",
        (
            "Read what a run left behind. Without arguments it lists the "
            "recorded runs; with run it returns that run's typed outcome, "
            "how each node ended and why; with program and artifact_id it "
            "returns which selectors one finished result actually resolves. "
            + str(outcome.get("description") or "")
            + " "
            + str(selectors.get("description") or "")
        ),
        inspect_run_properties,
        (),
    )
    merged = {
        "inspect_program": inspect_program,
        "project_yaml": project_yaml,
        "compile_command": compile_command,
        "inspect_run": inspect_run,
    }
    withdrawn = {
        legacy for names in MERGED_PLANNING_TOOLS.values() for legacy in names
    }
    first_member = {
        names[0]: merged_name
        for merged_name, names in MERGED_PLANNING_TOOLS.items()
    }
    result: list[dict] = []
    for item in tools:
        name = item["function"]["name"]
        if name in first_member:
            result.append(merged[first_member[name]])
        elif name not in withdrawn:
            result.append(item)
    return tuple(result)


def stem_operations(guides: tuple[str, ...] = ()) -> tuple[str, ...]:
    """The operation vocabulary the surface exposes: every operation that
    belongs to no leaf, plus the operations of the open guides."""

    from chemsmart.agent.guides import LEAF_OPERATIONS

    active = set(guides)
    return tuple(
        sorted(
            name
            for name in OPERATION_DESCRIPTIONS
            if LEAF_OPERATIONS.get(name) is None
            or LEAF_OPERATIONS[name] in active
        )
    )


def build_command_compiled_tool_surface(
    registry: ProgramCapabilityRegistryV1 | None = None,
    *,
    guides: tuple[str, ...] = (),
) -> AgentToolSurfaceV1:
    """The planning surface the model reads: the stem, plus the tools and
    operations of every open guide."""

    from chemsmart.agent.guides import LEAF_TOOLS

    active = set(guides)
    tools = _merge_planning_tools(
        _legacy_tool_definitions(registry, operations=stem_operations(guides))
    )
    tools = tuple(
        item
        for item in tools
        if LEAF_TOOLS.get(item["function"]["name"]) is None
        or LEAF_TOOLS[item["function"]["name"]] in active
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

    # The executor drives nodes through the legacy names one step at a
    # time (executor.PROGRAM_NODE_SEQUENCE); it never reads this surface as
    # a prompt, so the model-facing merge does not apply here.
    # ``inspect_calculation_artifact`` belongs to the legacy externally
    # seeded verifier surface: it requires separate settings-object and run
    # receipt IDs. Runtime V2 execution instead creates a typed program result
    # validation receipt and registers the resulting artifact directly; it
    # never binds those legacy IDs. Advertising the verifier here made live
    # models guess impossible identifiers after an otherwise valid run.
    tools = tuple(
        item
        for item in _legacy_tool_definitions(registry)
        if item["function"]["name"] != "inspect_calculation_artifact"
        and not _requires_an_unbound_registry(item)
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
    "blocked_reason": (
        "Why this stage cannot run here, stated for the human review -- for "
        "example the program validator's refusal of its functional."
    ),
    "point_index": (
        "1-based position of the scan point whose converged geometry to "
        "carry forward, counted in step order along the driven coordinate. "
        "Read scan_point_indices with scan_coordinate_values and "
        "scan_energies first; the choice is yours and is recorded as such."
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
        "inspect_program for this program and engine."
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
    "canonical invocation": "by compile_command",
    "capability receipt": "by inspect_program",
    "command context": "by compile_command",
    "command inspection receipt": "by compile_command",
    "engine binding": "by inspect_program",
    "functional equivalence receipt": "by validating a project document",
    "program binding": "by inspect_program",
    "program validator receipt": "by compile_command",
    "project render receipt": "by project_yaml(action=render)",
    "project validation receipt": "by project_yaml(action=validate or establish)",
    "run receipt": "by executing an approved node",
    "safe preview receipt": "by compile_command",
    "scientific claim evidence": "by extracting quantities from a result",
    "scientific identity": "by binding an approved molecular identity",
    "settings object": "by validating a project document",
    "trusted artifact": "by recording a workspace file as an artifact",
}


#: Which host registry each late-bound argument indexes.  A tool taking one of
#: these cannot succeed until something else has run, and saying so only in the
#: rejection means the model learns it by failing.  Observed live:
#: assess_program_candidate called with no claim evidence bound.
LATE_BOUND_ARGUMENTS: dict[str, str] = {
    "command_inspection_receipt_sha256": "command inspection receipt",
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
                name: (
                    _describe(name, _walk_constrain(name, schema))
                    if isinstance(schema, dict)
                    else schema
                )
                for name, schema in properties.items()
            }
            function["parameters"] = parameters
        sentence = _precondition_sentence(parameters.get("properties") or {})
        if sentence and "PRECONDITION" not in function.get("description", ""):
            function["description"] = function["description"] + sentence
        placed = render_rules(f"tool:{function['name']}")
        if placed and placed not in function.get("description", ""):
            function["description"] = function["description"] + " " + placed
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
            key: (
                _describe(key, _walk_constrain(key, value))
                if isinstance(value, dict)
                else value
            )
            for key, value in properties.items()
        }
    items = updated.get("items")
    if isinstance(items, dict):
        updated["items"] = _walk_constrain(name, items)
    return updated


def _describe_string(description: str) -> dict:
    """A free-form string that states what it must join to.

    The plan-time joins are enforced and were unstated, which is how a
    finished nine-species profile lost every claim: the rule existed,
    the refusal named it, and the field where the id gets chosen said
    nothing. Point of use is where a sentence changes behaviour.
    """

    return {"type": "string", "description": description}


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


#: The spelling rule every public identifier follows. Stated in full on
#: one field (plan_scientific_workflow.workflow_id); every other field
#: points there. It appeared verbatim 42 times, 12 KB per turn.
_IDENTIFIER_SPELLING_RULE = (
    "Lower-case public identifier; use dots, dashes, or underscores "
    "instead of spaces, parentheses, hashes, or placeholder syntax. "
    "It must begin with a letter, so a name taken from a compound "
    "whose locants come first needs a leading word: "
    "'dfe-12-rotamers', not '12-difluoroethane'. Chemical notation is "
    "mixed case and this field is not: unit symbols and quantity "
    "names must be folded down, so write 'gap-adiab-ev' not "
    "'gap-adiab-eV', 'delta-e' not 'dE', and 'ddg-compose' not "
    "'compose-ddG'. Fold the case; do not drop the letters."
)


#: The quantity kinds a thermochemistry node can declare as outputs.
#: Two live sessions declared a Gibbs correction as kind "energy", were
#: refused when planned with the list in the message, and repeated it
#: in a later session; the list belongs on the field.
_THERMOCHEMISTRY_KINDS = tuple(
    sorted(derivable_thermochemistry_quantities("rrho"))
)


def _public_identifier(
    joins: str | None = None, *, spelling_rule: bool = False
) -> dict:
    """A public identifier, optionally stating what it must join to.

    The spelling rule is shared by every identifier on this surface.
    What an id must *match* is not: an extraction output names a
    selector on its own node, a producer_output_id names an output on
    another node, a validation rule names an input on its own node.
    Those joins are enforced when the plan is checked, and a session
    that learns of one from a refusal has already spent the engines.

    A live nine-species reaction profile validated every node and then
    lost every claim to exactly this, so the join belongs on the field
    where the id is chosen rather than only in the error after it.
    """

    description = (
        _IDENTIFIER_SPELLING_RULE
        if spelling_rule
        else (
            # The pattern below states the shape machine-readably on
            # every one of these fields, require_identifier's refusal
            # quotes the full rule when one is malformed, and the rule
            # itself lives on plan_scientific_workflow.workflow_id. A
            # prose restatement on all 46 of them, and then even a
            # pointer to it, were paying 7.6 kB and 2.5 kB of a surface
            # at its ceiling for what the schema already enforces.
            "Lower-case id, case folded."
        )
    )
    if joins:
        description = f"{description} {joins}"
    return {
        "type": "string",
        "pattern": "^[a-z][a-z0-9_.-]*$",
        "description": description,
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
            # Optional: the digest of a host-recorded anomaly this node
            # investigates. Charged to the envelope's excursion line, never
            # to the engine-call budget; it may feed no untagged node, so the
            # asked observable stays owed to the budget.
            "excursion": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
            "dependencies": {
                "type": "array",
                "description": (
                    "Node ids this node runs after. Two rules bind "
                    "here and both are refused when planned: list the "
                    "nodes in topological order, so every dependency "
                    "appears earlier in this list than the node naming "
                    "it; and every producer_node_id on this node's "
                    "inputs must ALSO appear here, because a data edge "
                    "does not imply the ordering edge -- a producer "
                    "that is not a direct dependency is refused."
                ),
                "items": _string(),
            },
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
                                "'hess_filename'/orca_hessian. An ORCA TS "
                                "search that starts from a producer's Hessian "
                                "(any frequency-bearing ORCA stage; any "
                                "imaginary-mode count) likewise has "
                                "'filename'/geometry_xyz plus "
                                "'inhess_filename'/orca_hessian, and its "
                                "project ts section must set inhess: true. A "
                                "'filename'/geometry_xyz input fed by an ORCA "
                                "scan node carries the scan's minimum-energy "
                                "sampled point under one approval -- declare "
                                "the edge and the consumer defers until the "
                                "surface exists; carrying any other point is "
                                "the bind_scan_point_geometry route with its "
                                "own new workflow."
                            ),
                        },
                        "artifact_id": _string(),
                        "artifact_class": _string(),
                        "producer_node_id": _describe_string(
                            "The node_id whose output this binding "
                            "consumes. It must be a node earlier in "
                            "this plan and must also be listed in this "
                            "node's dependencies; naming a producer "
                            "without the matching dependency is "
                            "refused when planned."
                        ),
                        "producer_output_id": _describe_string(
                            "Which of that producer's expected_outputs "
                            "this binding reads, by its output_id. It "
                            "must be an output_id the named producer "
                            "itself declares -- the edge resolves "
                            "against that node's own list, so a name "
                            "invented for the quantity is refused."
                        ),
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
                        "output_id": _describe_string(
                            "Name for one artifact this node produces. A "
                            "consumer cites it as producer_output_id, and "
                            "that edge resolves against this list, so a "
                            "downstream node naming anything else is "
                            "refused when planned."
                        ),
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
                    "project_yaml(action=validate). In particular, when that tool says "
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


def _analysis_intent_node_schema(
    *, operations: tuple[str, ...] | None = None
) -> dict:
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
                "description": (
                    "Ordering-only edges, named by node_id. Every id must "
                    "be a node in this same plan -- a calculation node or "
                    "another analysis node -- and an id that matches "
                    "nothing is refused when planned. Data edges belong in "
                    "inputs; use this only to sequence work that carries "
                    "no value between the nodes."
                ),
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
                        "producer_node_id": _public_identifier(
                            "The node_id producing this input. It must "
                            "name a node in this same plan, and "
                            "source_kind must match what that node is: "
                            "program_output for a calculation node, "
                            "analysis_output for an analysis node. A "
                            "producer nothing declares is refused when "
                            "planned."
                        ),
                        "producer_output_id": _public_identifier(
                            "Which of that producer's declared outputs "
                            "this input reads, named by its output_id. It "
                            "must be an output the named producer itself "
                            "declares -- naming the quantity you want "
                            "rather than the output that carries it is "
                            "refused when planned, because the edge is "
                            "resolved against the producer's own list."
                        ),
                        "uncertainty_producer_node_id": _string(),
                        "uncertainty_producer_output_id": _describe_string(
                            "On a claim_rendering input: the analysis "
                            "output that computes this number's "
                            "uncertainty, with its node in "
                            "uncertainty_producer_node_id. You cannot "
                            "know an estimate before the results exist, "
                            "but you can plan its estimator -- a spread "
                            "over methods, a conformer range -- and the "
                            "host evaluates it in the same run and "
                            "records the claim as measured, cited to "
                            "that output. Leave it out and the delivery "
                            "states no uncertainty, which reopens the "
                            "goal to assess it."
                        ),
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
                "description": (
                    "What this node produces, one entry per quantity a "
                    "later node or claim can name. At least one is "
                    "required, because a node that declares nothing "
                    "produces nothing another node can cite and the "
                    "plan is refused. An extraction node's outputs are "
                    "the quantities its selectors read; a "
                    "thermochemistry node's are the terms it derives; "
                    "an expression node's are the values it computes. "
                    "Downstream inputs cite these by output_id."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "output_id": _public_identifier(
                            "Public name for this output. When the node "
                            "carries more than one selector, each "
                            "output_id must equal one of this node's own "
                            "selector quantity_ids: with several "
                            "selectors the host cannot tell which "
                            "extracted quantity an unmatched name meant, "
                            "so it is refused when planned and the "
                            "message names the selectors available. "
                            "Downstream inputs cite this id as "
                            "producer_output_id."
                        ),
                        "quantity_kind": _public_identifier(
                            "For a thermochemistry node, one of the kinds "
                            "the derivation produces -- "
                            + ", ".join(_THERMOCHEMISTRY_KINDS)
                            + " -- so a Gibbs correction is "
                            "thermal_gibbs_correction and a free energy is "
                            "gibbs_free_energy, never a bare 'energy'. "
                            "Extraction and expression outputs name their "
                            "own kinds."
                        ),
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
                "description": (
                    "The typed expression DAG. Give the nodes in any order: "
                    "the host orders them so each follows what it reads, "
                    "and refuses when planned only a name that no analysis "
                    "input or expression node provides, or nodes that read "
                    "each other in a cycle, naming the node and the name."
                ),
                "items": _quantity_expression_node_schema(
                    compact=True, operations=operations
                ),
            },
            "expression_output_node_ids": {
                "type": "array",
                "description": (
                    "Which of this node's own expression nodes are the "
                    "reported outputs; the rest are intermediates. Every "
                    "id here must be the node_id of a node listed above -- "
                    "the receipt keys each produced quantity by the node "
                    "that computed it, so an output naming anything else "
                    "cannot be recorded and the whole expression is "
                    "refused. Name the observable on the node that "
                    "produces it rather than inventing a separate "
                    "reporting name; a live run computed a complete "
                    "nine-species reaction profile and lost every claim "
                    "to this, because its outputs were named after the "
                    "quantities it wanted rather than after the nodes it "
                    "had written."
                ),
                "items": _public_identifier(),
            },
            "validation_rules": {
                "type": "array",
                "description": (
                    "Program-neutral validation predicates over named typed "
                    "inputs; do not hide criteria only in prose. A "
                    "scientific_validation node declares exactly one "
                    "dimensionless verdict output. For imaginary-mode "
                    "criteria there are two different questions -- say which "
                    "you asked. Strict: minimum_greater_equal with threshold "
                    "0 cm-1 fails on ANY negative mode; use it when the task "
                    "demands zero imaginary modes without tolerance. "
                    "Convention: ChemSmart's thermochemistry treats a mode "
                    "within 20 cm-1 of zero as numerical noise "
                    "(near_zero_mode_count in the derive_thermochemistry "
                    "receipt counts them), so minimum_greater_equal with "
                    "threshold -20 cm-1 asks the question the "
                    "thermochemistry answers. Default to the 20 cm-1 "
                    "convention when the task states no criterion; an "
                    "explicit criterion in the task text always overrides "
                    "the default."
                ),
                "items": {
                    "type": "object",
                    "properties": {
                        "rule_id": _public_identifier(),
                        "predicate": {
                            "type": "string",
                            # Derived, not hand-listed. The eight names were
                            # typed here as literals while the constant they
                            # mirror lived in the toolchain module, with
                            # nothing keeping the two lists together -- so a
                            # ninth predicate would have been invisible to the
                            # model until somebody noticed. Every other
                            # model-facing vocabulary here is derived from its
                            # source of truth; this one now is too.
                            "enum": sorted(ANALYSIS_VALIDATION_PREDICATES),
                        },
                        "input_ids": {
                            "type": "array",
                            "minItems": 1,
                            "description": (
                                "Which of this validation node's own "
                                "declared inputs the predicate reads, by "
                                "input_id. Every id must be an input_id "
                                "listed on this same node -- not a "
                                "quantity name, not an upstream output_id "
                                "-- and a rule naming anything else is "
                                "refused when planned, naming the rule "
                                "and the unknown ids. Bind the value as "
                                "an input first, then cite that input_id "
                                "here."
                            ),
                            "items": _public_identifier(),
                        },
                        "threshold": {
                            "type": "number",
                            "description": (
                                "Only for maximum_absolute_less_equal, "
                                "minimum_greater_equal and symmetric_within, "
                                "each of which also requires unit."
                            ),
                        },
                        # No schema-level minimum: count_equals enforces
                        # non-negativity in its own contract, while
                        # integer_equals exists precisely for negative
                        # state labels (an anion's charge of -1) -- a bound
                        # here made that case unreachable at the tool
                        # surface.
                        "expected_count": {
                            "type": "integer",
                            "description": (
                                "Only for count_equals (non-negative) and "
                                "integer_equals (any integer, e.g. an "
                                "anion's charge of -1)."
                            ),
                        },
                        "unit": {
                            **_string(),
                            "description": (
                                "The unit of threshold; required with it "
                                "and accepted by nothing else. all_equal, "
                                "all_equal_text and all_finite take "
                                "input_ids alone."
                            ),
                        },
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


def _quantity_expression_node_schema(
    *, compact: bool = False, operations: tuple[str, ...] | None = None
) -> dict:
    """The expression node, in full on evaluate_quantity_expression and
    compact on the planner: the same fields and enums, with the operation
    semantics and the constant purposes stated once. The full form was
    serialised twice, 31 KB of the surface every turn."""

    operation_description = (
        "Operation name; its meaning, arity, and conventions are stated "
        "on evaluate_quantity_expression.nodes.operation, and the named "
        "convention operations there are preferred over rebuilding a "
        "convention from arithmetic."
        if compact
        else (
            "Pick the operation that owns the step. Where a named "
            "operation exists for a scientific convention, use it "
            "rather than rebuilding the convention from arithmetic "
            "primitives: the named one carries the convention, its "
            "validity conditions, and its provenance. "
            + " | ".join(
                f"{name}: {text}"
                for name, text in sorted(OPERATION_DESCRIPTIONS.items())
                if operations is None or name in set(operations)
            )
        )
    )
    exposed = (
        sorted(OPERATION_DESCRIPTIONS)
        if operations is None
        else sorted(
            name for name in OPERATION_DESCRIPTIONS if name in set(operations)
        )
    )
    return {
        "type": "object",
        "properties": {
            "node_id": _public_identifier(),
            "operation": {
                "type": "string",
                "enum": exposed,
                "description": operation_description,
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
            "constant_name": {
                "type": "string",
                "enum": sorted(LITERATURE_CONSTANTS),
                "description": (
                    "For constant, the registered literature-constant name; "
                    "units, convention families and purposes are listed on "
                    "evaluate_quantity_expression.nodes.constant_name."
                    if compact
                    else "For constant, the registered literature-constant name. "
                    "The host owns the value, unit, and standard-state "
                    "convention; an unregistered name is refused when "
                    "planned. Other operations omit this. Registered names, "
                    "each with its unit, the convention family it may be "
                    "combined within, and what it is for: "
                    + " || ".join(
                        f"{name} [{entry.unit}, {entry.convention_family}]"
                        + (f" -- {entry.purpose}" if entry.purpose else "")
                        for name, entry in sorted(LITERATURE_CONSTANTS.items())
                    )
                    + ". Read the purpose before composing several entries by "
                    "hand: a family says which scale an entry sits on and "
                    "says nothing about its standard state, and two entries "
                    "on one scale at different standard states still need "
                    "the term that bridges them. Where a finished composed "
                    "value is registered, prefer it. Two entries sharing a "
                    "family combine freely. Two "
                    "families in one chain is not refused, but it is "
                    "displayed to the reviewer, because constants determined "
                    "on different scales are not interchangeable even when "
                    "each is correct on its own -- pick the electrode "
                    "potential belonging to the solvation scale you used. "
                    "An entry marked independent combines with any family. "
                    "The values themselves are deliberately not listed here: "
                    "select by meaning, and let the host supply the number."
                ),
            },
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
                    "Positive exponent for a two-point CBS limit. You "
                    "supply it, so the receipt records it as "
                    "model-authored; the cbs guide says when to."
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
    "MERGED_PLANNING_TOOLS",
    "stem_operations",
    "PROJECT_YAML_ACTIONS",
    "PROJECT_YAML_ACTION_ARGUMENTS",
    "AgentToolSurfaceV1",
    "build_approved_execution_tool_surface",
    "build_command_compiled_tool_surface",
]
