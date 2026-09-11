"""The human-approval display: everything the single /approve covers.

Charter step 3: the terminal shows the complete plan -- molecular identity,
electronic state, effective settings, CLI operations, dependencies,
environment, resources -- and, because the one approval also covers them,
the typed analysis chain carried with the workflow and the full lineage of
any composed arrangement. Every panel is a projection of the same typed
review object; the identical canonical bytes live on in the run evidence,
and no digest is ever displayed to the human.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Mapping

from rich.console import Group
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1

from chemsmart.agent.voice import human_identity_evidence, human_state

from .presentation import _canonical_tool_results, human_cli_operation


def _walk_render_records(value: Any) -> Iterable[Mapping[str, Any]]:
    if isinstance(value, Mapping):
        if "rendered_yaml" in value and "rendered_sha256" in value:
            yield value
        for item in value.values():
            yield from _walk_render_records(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _walk_render_records(item)


def resolve_project_yaml_texts(
    review: "WorkflowExecutionReviewV1",
    *,
    public_transcript: Iterable[Mapping[str, Any]] = (),
    workspace: Path | None = None,
) -> dict[str, str]:
    """Readable project YAML per reviewed node, from session evidence.

    Promotion enforces that a node's project artifact digest equals the
    digest of the rendered YAML bytes, so the join is exact: first over the
    session transcript's render/promotion records, then over the promoted
    YAML files kept under the planning runs, and a node that resolves
    nowhere keeps its canonical-settings fallback.
    """

    wanted = {
        digest: item.node_id
        for item in review.node_reviews
        if (digest := str(getattr(item, "project_artifact_sha256", "")))
    }
    texts: dict[str, str] = {}
    for _tool, record in _canonical_tool_results(public_transcript):
        for found in _walk_render_records(record):
            digest = str(found.get("rendered_sha256") or "")
            node_id = wanted.get(digest)
            rendered = found.get("rendered_yaml")
            if node_id and node_id not in texts and isinstance(rendered, str):
                texts[node_id] = rendered
    missing = set(wanted.values()) - set(texts)
    if workspace is not None and missing:
        pattern = ".chemsmart-agent/runs/*/projects/*.yaml"
        for path in sorted(Path(workspace).glob(pattern)):
            try:
                data = path.read_bytes()
            except OSError:
                continue
            node_id = wanted.get(hashlib.sha256(data).hexdigest())
            if node_id and node_id not in texts:
                texts[node_id] = data.decode("utf-8", errors="replace")
    return texts


def _overview_table(review: "WorkflowExecutionReviewV1") -> Table:
    review_by_id = {item.node_id: item for item in review.node_reviews}
    deferred = set(review.non_executable_node_ids)
    overview = Table(title="ChemSmart workflow awaiting human approval")
    overview.add_column("Node", style="bold cyan")
    overview.add_column("Program / engine")
    overview.add_column("Stage")
    overview.add_column("Molecule")
    overview.add_column("Charge / multiplicity")
    overview.add_column("Execution")
    overview.add_column("Reason")
    for planned in review.scientific_plan.nodes:
        item = review_by_id.get(planned.node_id)
        if item is not None:
            identity = item.molecular_identity
            approved_names = identity.get("approved_names") or ()
            name = str(approved_names[0]) if approved_names else ""
            formula = str(identity.get("formula") or "unknown formula")
            atom_order = identity.get("atom_order") or ()
            atom_summary = "-".join(str(symbol) for symbol in atom_order)
            molecule = name or formula
            if atom_summary:
                molecule += f" · {atom_summary}"
            program_engine = f"{item.program} / {item.engine}"
            stage = item.stage
            molecular_state = (
                f"{identity.get('charge')} / "
                f"{identity.get('multiplicity')}"
            )
            execution_state = (
                "Excursion" if planned.excursion else "Executable"
            )
            reason = (
                f"investigates anomaly {planned.excursion[:8]}; charged "
                "to the excursion line"
                if planned.excursion
                else "reviewed below"
            )
        else:
            molecule = planned.project_role
            program_engine = f"{planned.program} / {planned.engine}"
            stage = planned.stage
            if planned.charge is not None and planned.multiplicity is not None:
                molecular_state = f"{planned.charge} / {planned.multiplicity}"
            else:
                source_ids = tuple(
                    edge.source_node_id
                    for edge in review.scientific_plan.edges
                    if edge.target_node_id == planned.node_id
                )
                source_review = next(
                    (
                        review_by_id[source_id]
                        for source_id in source_ids
                        if source_id in review_by_id
                    ),
                    None,
                )
                if source_review is not None:
                    source_identity = source_review.molecular_identity
                    molecular_state = (
                        f"{source_identity.get('charge')} / "
                        f"{source_identity.get('multiplicity')} "
                        f"(from {source_review.node_id})"
                    )
                else:
                    molecular_state = "not separately reviewed"
            execution_state = (
                "Deferred" if planned.node_id in deferred else "Not executable"
            )
            reason = planned.blocked_reason or "no reviewed command"
        overview.add_row(
            planned.node_id,
            program_engine,
            stage,
            molecule,
            molecular_state,
            execution_state,
            reason,
        )
    return overview


def _environment_table(review: "WorkflowExecutionReviewV1") -> Table:
    environments = Table(title="Observed execution environments")
    environments.add_column("Node", style="bold cyan")
    environments.add_column("Status")
    environments.add_column("Target")
    environments.add_column("Version")
    environments.add_column("Observed by")
    for item in review.node_reviews:
        summary = item.environment_summary
        dependencies = {
            str(name): str(version)
            for name, version in summary.get("dependency_versions", ())
        }
        version = str(summary.get("observed_version") or "")
        if not version:
            version = dependencies.get(item.program, "")
        if not version and dependencies:
            version = ", ".join(
                f"{name} {value}"
                for name, value in sorted(dependencies.items())
            )
        environments.add_row(
            item.node_id,
            str(summary.get("status") or "unknown"),
            str(summary.get("target_kind") or "unknown"),
            version or "not reported",
            str(summary.get("observation_method") or "host probe"),
        )
    return environments


def _bounds_panel(review: "WorkflowExecutionReviewV1") -> Panel:
    resources = review.execution_resources
    excursions = review.execution_envelope.get("max_excursion_calls")
    lines = [
        f"cores: {resources.cores}",
        f"memory: {resources.memory_gb:g} GB",
        f"GPUs: {resources.gpu_count}",
        f"node timeout: {resources.node_timeout_seconds} s",
        f"engine-call limit: "
        f"{review.execution_envelope.get('max_engine_calls')}",
    ]
    if excursions:
        lines.append(
            f"excursion-call limit: {excursions} (investigation of a "
            "recorded anomaly; never the deliverable)"
        )
    return Panel("\n".join(lines), title="Execution bounds")


def _edges_table(review: "WorkflowExecutionReviewV1") -> Table:
    edge_table = Table(title="Scientific data flow")
    edge_table.add_column("Producer")
    edge_table.add_column("Artifact")
    edge_table.add_column("Consumer")
    for edge in review.scientific_plan.edges:
        edge_table.add_row(
            edge.source_node_id,
            edge.artifact_class or edge.edge_kind,
            edge.target_node_id,
        )
    return edge_table


def _composition_panels(review: "WorkflowExecutionReviewV1") -> list[Panel]:
    panels: list[Panel] = []
    for item in review.node_reviews:
        record = item.molecular_identity.get("composition")
        if not isinstance(record, Mapping):
            continue
        contact = record.get("placement", {}).get("contact", {})
        status = item.molecular_identity.get("identity_evidence_status")
        lines = [
            f"parent A: {record.get('fragment_a_artifact_id')}",
            f"parent B: {record.get('fragment_b_artifact_id')}",
            f"contact: atom {contact.get('fragment_a_atom')} of A to atom "
            f"{contact.get('fragment_b_atom')} of B (1-based)",
            f"requested distance: {contact.get('distance_angstrom')} Å · "
            f"achieved "
            f"{record.get('achieved_contact_distance_angstrom')} Å · "
            f"closest interfragment contact "
            f"{record.get('min_interfragment_distance_angstrom')} Å",
            f"composed: {record.get('formula')} "
            f"({record.get('atom_count')} atoms); "
            f"{record.get('atom_order_note')}",
            f"evidence: {human_identity_evidence(str(status))}",
        ]
        panels.append(
            Panel(
                "\n".join(lines),
                title=(
                    f"{item.node_id} · composed arrangement lineage "
                    "(covered by this approval)"
                ),
            )
        )
    return panels


def _route_directive_panel(review: "WorkflowExecutionReviewV1") -> Any:
    """Name every project route-hatch token, or return ``None``.

    The hatch carries two kinds of thing.  A source-required scientific
    keyword changes the calculation, and a print directive changes only what
    the program reports -- ORCA's ``Hirshfeld`` being the case that matters,
    because the condensed-Fukui literature asks for that scheme and ORCA
    prints it only when told to.  Both are already inside the effective
    project settings the reviewer can read per node, but a single token in a
    settings dump is exactly what a reader skims past, and this is the one
    channel where free text reaches a program's input.  So it is named on
    its own, next to the node it belongs to.

    Which kind a token is stays the reviewer's reading.  The host does not
    classify a keyword it did not define.
    """

    import json

    rows = []
    for item in review.node_reviews:
        try:
            settings = json.loads(item.project_settings_text)
        except (TypeError, ValueError):
            continue
        if not isinstance(settings, Mapping):
            continue
        tokens = settings.get("additional_route_parameters")
        if not tokens:
            continue
        if isinstance(tokens, str):
            rendered = tokens
        else:
            rendered = " ".join(str(token) for token in tokens)
        rows.append((item.node_id, rendered.strip()))
    if not rows:
        return None
    table = Table(
        title=(
            "Project route directives (free text reaching the program input)"
        )
    )
    table.add_column("Node", style="bold cyan")
    table.add_column("Tokens")
    for node_id, rendered in rows:
        table.add_row(node_id, rendered)
    return table


def _levels_of_theory(
    review: "WorkflowExecutionReviewV1",
) -> dict[str, str]:
    """Map each calculation node to a one-line level of theory.

    A typed value carries its unit and its dimension and not the method that
    produced it, and the arithmetic checks only those two.  So a GFN2 energy
    minus a hybrid-DFT energy is accepted, and so is a Mulliken charge from
    one program minus a Mulliken charge from another at a different basis.

    Refusing that would be wrong in both directions.  A high-level single
    point on a low-level geometry is the most ordinary multi-program
    protocol there is, and composite methods mix levels on purpose; two
    energies from one program at different basis sets are just as
    unsubtractable and no program check would catch them.  What was missing
    is not a rule but a sightline: the level was already known per node and
    simply absent from the surface where the arithmetic is approved.
    """

    import json

    levels: dict[str, str] = {}
    for item in review.node_reviews:
        parts = [str(item.program)]
        try:
            settings = json.loads(item.project_settings_text)
        except (TypeError, ValueError):
            settings = {}
        if isinstance(settings, Mapping):
            method = " ".join(
                str(settings[key]).strip()
                for key in ("functional", "ab_initio", "method")
                if settings.get(key)
            ).strip()
            basis = str(settings.get("basis") or "").strip()
            if method and basis:
                parts.append(f"{method}/{basis}")
            elif method or basis:
                parts.append(method or basis)
            solvent = str(settings.get("solvent_id") or "").strip()
            model = str(settings.get("solvent_model") or "").strip()
            if solvent and model:
                parts.append(f"{model}({solvent})")
            elif solvent or model:
                parts.append(solvent or model)
        levels[item.node_id] = " · ".join(part for part in parts if part)
    return levels


def _levels_carried_by_analysis_node(
    plan: Any, levels: Mapping[str, str]
) -> dict[str, set[str]]:
    """Map each analysis node to every level of theory its values came from.

    Resolved through the analysis DAG rather than one hop, because the
    mixture that matters most is usually a hop away: an expression reading
    two other expressions is exactly where one program's number meets
    another's, and a column that goes blank there would invite a reader to
    take the blank for "nothing mixed here".
    """

    carried: dict[str, set[str]] = {}
    nodes = tuple(getattr(plan, "analysis_nodes", ()) or ())
    for _pass in range(len(nodes) + 1):
        changed = False
        for node in nodes:
            seen = set(carried.get(node.node_id, ()))
            for item in node.inputs:
                producer_id = getattr(item, "producer_node_id", None)
                if not producer_id:
                    continue
                direct = levels.get(producer_id)
                if direct:
                    seen.add(direct)
                else:
                    seen |= carried.get(producer_id, set())
            if seen != carried.get(node.node_id, set()):
                carried[node.node_id] = seen
                changed = True
        if not changed:
            break
    return carried


def _analysis_chain_renderable(review: "WorkflowExecutionReviewV1") -> Any:
    plan = review.scientific_toolchain_plan
    if plan is None:
        return Panel(
            "No typed analysis chain is planned with this workflow. "
            "Approval covers the calculation nodes only; analysis of the "
            "completed results is a later explicit request.",
            title="Analysis chain",
        )
    table = Table(
        title=(
            "Typed analysis chain (runs provider-free after every approved "
            "calculation node validates; covered by this approval)"
        )
    )
    table.add_column("Node", style="bold cyan")
    table.add_column("Kind")
    table.add_column("Inputs")
    table.add_column("Level of theory")
    table.add_column("Outputs")
    table.add_column("Conditions")
    table.add_column("State")
    carried = _levels_carried_by_analysis_node(plan, _levels_of_theory(review))
    for node in plan.analysis_nodes:
        inputs = ", ".join(
            (
                f"{item.producer_node_id}.{item.producer_output_id}"
                if hasattr(item, "producer_node_id")
                # A registered-result input names an artifact, not a
                # producer node; rendering it used to crash the review.
                else f"artifact:{item.artifact_id}"
            )
            for item in node.inputs
        )
        # Which levels of theory feed this node's arithmetic.  Displayed,
        # never refused: mixing them can be exactly right or exactly wrong,
        # and that is the reviewer's call to make with it in front of them.
        producers = sorted(carried.get(node.node_id, ()))
        outputs = ", ".join(
            f"{item.output_id} ({item.unit})" for item in node.outputs
        )
        conditions = []
        if node.temperature_k is not None:
            conditions.append(f"{node.temperature_k:g} K")
        concentration = getattr(node, "concentration_mol_l", None)
        if concentration is not None:
            # The standard state the receipt will bind: never display a
            # solution free energy under a gas-phase pressure label.
            conditions.append(f"{concentration:g} mol/L")
        elif node.pressure_atm is not None:
            conditions.append(f"{node.pressure_atm:g} atm")
        entropy_method = getattr(node, "entropy_method", "rrho")
        if node.temperature_k is not None and entropy_method != "rrho":
            conditions.append(entropy_method)
        scale = getattr(node, "frequency_scale_factor", 1.0)
        if node.temperature_k is not None and scale != 1.0:
            conditions.append(f"scale {scale:g}")
        state = human_state(node.support_state)
        if node.blocked_reason:
            state += f": {node.blocked_reason}"
        table.add_row(
            node.node_id,
            node.analysis_kind,
            inputs or "—",
            "\n".join(producers) or "—",
            outputs or "—",
            ", ".join(conditions) or "—",
            state,
        )
    return table


def _literature_constants_renderable(
    review: "WorkflowExecutionReviewV1",
) -> Any:
    """Show every host-owned constant the planned chain selects, or None.

    The reviewer approves a cycle whose answer moves with these numbers, so
    the value, unit, and standard-state convention appear at the decision
    surface -- resolved from the registry by name, never restated by the
    session that planned the chain.

    The convention family is here for a narrower reason.  Constants that
    look independent are often matched pairs -- an absolute electrode
    potential means one thing beside the proton solvation free energy
    determined on the same scale and another beside a different one -- and
    the literature circulates the halves separately.  A crossed pair fails
    silently: both values are correct, the units are right, nothing
    diverges, and the answer is wrong.  So when a chain draws on more than
    one family the table says so, and the reviewer decides; mixing can be
    deliberate, and it is sometimes the only way to reproduce a published
    cycle.
    """

    plan = review.scientific_toolchain_plan
    if plan is None:
        return None
    names: list[str] = []
    for node in plan.analysis_nodes:
        if node.analysis_kind != "quantity_expression":
            continue
        for item in node.expression_nodes:
            if str(item.get("operation", "")) != "constant":
                continue
            name = str(item.get("constant_name", ""))
            if name and name not in names:
                names.append(name)
    if not names:
        return None
    from chemsmart.analysis.literature_constants import (
        UnknownLiteratureConstantError,
        literature_constant,
    )

    table = Table(
        title=(
            "Literature constants (host-owned values the chain selects "
            "by name)"
        )
    )
    table.add_column("Constant", style="bold cyan")
    table.add_column("Value", justify="right")
    table.add_column("Unit")
    table.add_column("Family")
    table.add_column("Convention")
    families: list[str] = []
    for name in names:
        try:
            entry = literature_constant(name)
        except UnknownLiteratureConstantError:
            table.add_row(name, "—", "—", "—", "no longer registered")
            continue
        if (
            entry.convention_family != "independent"
            and entry.convention_family not in families
        ):
            families.append(entry.convention_family)
        table.add_row(
            entry.name,
            f"{entry.value:g}",
            entry.unit,
            entry.convention_family,
            entry.convention,
        )
    if len(families) > 1:
        return Group(
            table,
            Panel(
                Text(
                    "This chain draws on "
                    f"{len(families)} convention families: "
                    + ", ".join(sorted(families))
                    + ". Constants determined on different scales are not "
                    "interchangeable even when each is correct on its own, "
                    "and a crossed pair shifts the answer without failing. "
                    "Confirm the combination is the one you intend.",
                    style="bold yellow",
                ),
                title="Mixed constant conventions",
                border_style="yellow",
            ),
        )
    return table


def _short_record_id(record_id: Any) -> str:
    """Truncate a content-derived database record id for human reading.

    The database family's own tables truncate record ids to twelve
    characters and every database CLI resolves unique prefixes, so the
    decision surface follows the same convention: the first live batch
    review printed six 64-character ids a human cannot read.  The full
    id stays in the review packet's typed lineage record.
    """

    text = str(record_id or "")
    return text[:12] + "…" if len(text) > 12 else text


def _batch_records_renderable(
    review: "WorkflowExecutionReviewV1",
) -> Any | None:
    """One row per record when the plan carries several disconnected ones.

    The batch state is derived from the per-record rows by reading them,
    never stored beside them, and no row is ever elided: each row carries
    that molecule's load-bearing facts (identity, explicit state, origin),
    so approving the batch means having seen every record.
    """

    from chemsmart.agent.workflows import workflow_record_components

    components = workflow_record_components(
        review.scientific_plan.nodes, review.scientific_plan.edges
    )
    if len(components) <= 1:
        return None
    review_by_id = {item.node_id: item for item in review.node_reviews}
    table = Table(
        title=(
            f"Batch records ({len(components)}) -- each row is one "
            "molecule's sub-workflow under this single approval"
        )
    )
    table.add_column("Record", style="bold cyan")
    table.add_column("Molecule")
    table.add_column("Charge / multiplicity")
    table.add_column("Stages")
    table.add_column("Origin")
    for index, group in enumerate(components, start=1):
        reviewed = [
            review_by_id[node_id]
            for node_id in group
            if node_id in review_by_id
        ]
        first = reviewed[0] if reviewed else None
        if first is not None:
            identity = first.molecular_identity
            approved_names = identity.get("approved_names") or ()
            molecule = (
                str(approved_names[0])
                if approved_names
                else str(identity.get("formula") or "unknown formula")
            )
            state = (
                f"{identity.get('charge')} / {identity.get('multiplicity')}"
            )
            extraction = identity.get("database_extraction")
            derivation = identity.get("derivation")
            if isinstance(extraction, Mapping):
                origin = (
                    f"db {extraction.get('database_filename')} record "
                    f"{_short_record_id(extraction.get('record_id'))}"
                )
                stored_charge = extraction.get("stored_charge")
                stored_multiplicity = extraction.get("stored_multiplicity")
                stored = []
                if stored_charge is not None:
                    stored.append(f"charge {stored_charge}")
                if stored_multiplicity is not None:
                    stored.append(f"mult {stored_multiplicity}")
                if stored:
                    origin += " (stored: " + ", ".join(stored) + ")"
                if (
                    stored_charge is not None
                    and identity.get("charge") is not None
                    and int(stored_charge) != int(identity["charge"])
                ) or (
                    stored_multiplicity is not None
                    and identity.get("multiplicity") is not None
                    and int(stored_multiplicity)
                    != int(identity["multiplicity"])
                ):
                    origin += (
                        " [bold red]!= bound state -- deliberate?[/bold red]"
                    )
            elif isinstance(derivation, Mapping):
                origin = (
                    f"derived from {derivation.get('parent_artifact_id')}, "
                    f"removed atoms "
                    f"{list(derivation.get('removed_atoms') or ())}"
                )
            else:
                coordinate = identity.get("coordinate_identity") or {}
                origin = str(coordinate.get("kind") or "workspace geometry")
        else:
            molecule = "not separately reviewed"
            state = ""
            origin = ""
        stage_by_id = {
            node.node_id: node.stage for node in review.scientific_plan.nodes
        }
        stages = " -> ".join(
            f"{node_id} ({stage_by_id.get(node_id, '?')})" for node_id in group
        )
        table.add_row(str(index), molecule, state, stages, origin)
    return table


def _derivation_panels(review: "WorkflowExecutionReviewV1") -> list[Panel]:
    """A derived species' kept/removed atoms, on the decision surface.

    A wrong-atom removal was once caught only from the artifact bytes;
    the lineage the host already records belongs on the page the human
    approves from.
    """

    panels: list[Panel] = []
    for item in review.node_reviews:
        record = item.molecular_identity.get("derivation")
        if not isinstance(record, Mapping):
            continue
        lines = [
            f"parent: {record.get('parent_artifact_id')} "
            f"({record.get('parent_formula')}, "
            f"{record.get('parent_atom_count')} atoms)",
            f"removed atoms (1-based): "
            f"{list(record.get('removed_atoms') or ())}",
            f"kept atoms (1-based): {list(record.get('kept_atoms') or ())}",
            f"derived: {record.get('formula')} "
            f"({record.get('atom_count')} atoms); "
            f"{record.get('fragment_count')} connected piece(s); "
            f"{record.get('atom_order_note')}",
            "electronic state was deliberately unbound at derivation; "
            f"this node binds charge "
            f"{item.molecular_identity.get('charge')}, multiplicity "
            f"{item.molecular_identity.get('multiplicity')} explicitly",
        ]
        panels.append(
            Panel(
                "\n".join(lines),
                title=(
                    f"{item.node_id} · derived species lineage "
                    "(covered by this approval)"
                ),
            )
        )
    return panels


def _pubchem_geometry_panels(
    review: "WorkflowExecutionReviewV1",
) -> list[Panel]:
    """What arrived when a molecule entered by public identifier.

    `pubchem_geometries` had no reader anywhere in the tree: the host
    minted the receipt, returned the formula, atom count and fragment
    count to the session, and none of it reached the page a human
    approves from -- while a derived species carried its full panel. A
    live session named two numeric CIDs from prior knowledge and got
    two unrelated molecules, one of 102 atoms in 27 pieces, each
    registered under a confident artifact id of the model's choosing.
    Neither reached a plan, so nothing was approved on them; that was
    the session reading the formula, not the review showing it.

    These are measurements and never refusals. A salt, an ion pair, a
    solvate and a metallocene as deposited are all legitimately more
    than one piece, and `compose_molecular_arrangement` exists to
    consume fragments -- so the count is stated and the scientist
    judges it.
    """

    panels: list[Panel] = []
    for item in review.node_reviews:
        record = item.molecular_identity.get("pubchem_geometry")
        if not isinstance(record, Mapping):
            continue
        pieces = record.get("fragment_count")
        lines = [
            f"requested identifier: {record.get('identifier')!r} "
            f"({record.get('identifier_kind')})",
            f"what arrived: {record.get('formula')} "
            f"({record.get('atom_count')} atoms); "
            f"{pieces} connected piece(s)",
            "a depositor's conformer, not a relaxed structure, and it "
            "carries that depositor's symmetry",
            "electronic state was deliberately unbound at fetch; this "
            f"node binds charge {item.molecular_identity.get('charge')}, "
            f"multiplicity {item.molecular_identity.get('multiplicity')} "
            "explicitly",
        ]
        if isinstance(pieces, int) and pieces > 1:
            lines.append(
                f"the record converted to {pieces} disconnected pieces -- "
                "read as an observation, not a verdict: an ion pair or a "
                "solvate is legitimately more than one"
            )
        panels.append(
            Panel(
                "\n".join(lines),
                title=(
                    f"{item.node_id} · molecule entered by identifier "
                    "(covered by this approval)"
                ),
            )
        )
    return panels


def _database_extraction_panels(
    review: "WorkflowExecutionReviewV1",
) -> list[Panel]:
    """A record-sourced geometry's lineage and stored-field observations."""

    panels: list[Panel] = []
    for item in review.node_reviews:
        record = item.molecular_identity.get("database_extraction")
        if not isinstance(record, Mapping):
            continue
        stored = (
            f"stored charge {record.get('stored_charge')}, "
            f"multiplicity {record.get('stored_multiplicity')}, "
            f"energy {record.get('stored_energy')}, "
            f"optimized {record.get('stored_is_optimized')}"
        )
        lines = [
            f"database: {record.get('database_filename')}",
            f"record: {_short_record_id(record.get('record_id'))} "
            f"(structure {record.get('structure_index')} of "
            f"{record.get('structure_count')})",
            f"extracted: {record.get('formula')} "
            f"({record.get('atom_count')} atoms); coordinates unchanged",
            stored + " -- observations from the record's own provenance, "
            "never bindings",
            f"this node binds charge "
            f"{item.molecular_identity.get('charge')}, multiplicity "
            f"{item.molecular_identity.get('multiplicity')} explicitly",
        ]
        panels.append(
            Panel(
                "\n".join(lines),
                title=(
                    f"{item.node_id} · database record lineage "
                    "(covered by this approval)"
                ),
            )
        )
    return panels


def _decision_panel(review: "WorkflowExecutionReviewV1") -> Panel:
    executable_text = ", ".join(item.node_id for item in review.node_reviews)
    deferred = set(review.non_executable_node_ids)
    decision = (
        "Review the molecule/state, settings, CLI operations, DAG, "
        "observed execution environments and resource bounds above. "
        "Enter /approve once to execute the reviewed nodes: "
        f"{executable_text}."
    )
    if review.scientific_toolchain_plan is not None:
        decision += (
            " The displayed analysis chain executes provider-free in the "
            "same run."
        )
    if deferred:
        decision += (
            " Deferred stages remain unapproved and will not launch: "
            + ", ".join(sorted(deferred))
            + "."
        )
    decision += (
        " A changed scientific request must be replanned. The full review "
        "record is kept in the run evidence."
    )
    return Panel(decision, title="Human decision", border_style="yellow")


def _advisory_knowledge_panel(
    review: "WorkflowExecutionReviewV1",
) -> Any | None:
    """Provenance of consulted skills; never a gate, never authority."""

    records = getattr(review, "consulted_domain_knowledge", ())
    if not records:
        return None
    lines = [
        f"{item.get('skill_id', '?')} {item.get('skill_version', '')}".strip()
        + f" ({item.get('origin', 'builtin')})"
        for item in records
    ]
    return Panel(
        Text("\n".join(lines)),
        title=(
            "Advisory domain knowledge consulted "
            "(carries no readiness or accuracy authority)"
        ),
    )


def _declared_observable_panel(
    review: "WorkflowExecutionReviewV1",
) -> Any | None:
    """What this workflow is for, and what the session expected.

    These are the session's own commitments, and unlike the consultation
    records they are consumed: the completion gate requires a delivered
    claim of matching dimension for every entry and states an unmatched
    one as a limitation.  Showing them beside the plan puts a recorded
    expectation in front of the reviewer *before* the grant, which is the
    only moment at which a wrong premise is still cheap.  An expectation
    is never a criterion -- only the validation verdicts are -- so it is
    displayed and never refused.
    """

    records = getattr(review, "requested_observable_declarations", ())
    if not records:
        return None
    lines = []
    for item in records:
        line = (
            f"{item.get('observable_id', '?')} "
            f"[{item.get('unit', '?')}] -- {item.get('meaning', '')}"
        )
        # The precision the task asked for, restated by the session and
        # frozen at this decision. It reached no reviewer at all, so a
        # model-proposed tolerance became a user-authorised one with no
        # human ever seeing it -- and an *absent* tolerance, on a task
        # that states a number in plain words, was the cheapest escape
        # from the contract and the least visible. Both are rendered
        # here, at the one moment a wrong premise is still cheap.
        tolerance = item.get("required_tolerance")
        if str(item.get("role") or "requested") == "requested":
            if tolerance is None:
                line += (
                    "\n    required precision: none declared -- if the "
                    "task states one, it is not on this contract"
                )
            else:
                line += (
                    f"\n    required precision: +/-{tolerance} "
                    f"{item.get('unit', '')}".rstrip()
                )
                basis = str(item.get("tolerance_basis") or "").strip()
                origin = str(item.get("tolerance_origin") or "").strip()
                # "from the task" was printed over every tolerance,
                # including one the session formulated for a decision of
                # its own -- a host word that can be false about whose
                # requirement the reader is approving.
                whose = {
                    "task": "from the task",
                    "session": "the session's own reading (the task "
                    "fixes none here)",
                }.get(origin, "source not stated")
                if basis:
                    line += f"\n    {whose}: {basis}"
                elif origin:
                    line += f"\n    {whose}"
        low, high = item.get("expected_low"), item.get("expected_high")
        sign = item.get("expected_sign", "")
        if sign or low is not None:
            expected = []
            if sign:
                expected.append(str(sign))
            if low is not None:
                expected.append(f"{low} to {high}")
            line += "\n    expected: " + ", ".join(expected)
            basis = item.get("expectation_basis", "")
            if basis:
                line += f"\n    basis: {basis}"
        lines.append(line)
    return Panel(
        Text("\n".join(lines)),
        title=(
            "Requested observables declared (expectations are displayed, "
            "never scored; a required precision is the contract)"
        ),
    )


def _shown_value(value: Any) -> Any:
    """Display a measured coordinate without a signed zero.

    A torsion set to zero is reached as something like -1e-17, and a record
    written before the receipts normalised it still carries -0.0. Reading
    "achieved -0.0 degree" on the page you approve from looks like a defect,
    so the display collapses the two zeros as well.
    """

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return float(value) + 0.0
    return value


def _geometry_lineage_panels(
    review: "WorkflowExecutionReviewV1",
) -> list[Panel]:
    """Every hop of a built geometry's edit chain, on the decision surface.

    A built structure may be several operations deep -- the cis rotamer with
    both methyls staggered is three edits -- and the hop that decides what
    the molecule IS can sit at the root of the chain.  The first chained
    build displayed only the final hop, which hid exactly the edit a
    reviewer would need to catch, so every hop renders in the order it was
    performed.  The atoms carry their element symbols because that is what
    makes a wrong-atom edit catchable by reading.
    """

    panels: list[Panel] = []
    for item in review.node_reviews:
        chain = item.molecular_identity.get("geometry_lineage")
        if not isinstance(chain, (list, tuple)) or not chain:
            continue
        total = len(chain)
        for position, record in enumerate(chain, start=1):
            if not isinstance(record, Mapping):
                continue
            step = f"step {position} of {total}" if total > 1 else "one step"
            if record.get("kind") == "atom_append":
                anchors = list(record.get("anchor_atoms") or ())
                symbols = list(record.get("anchor_symbols") or ())
                named = ", ".join(
                    f"{symbol}{index}"
                    for symbol, index in zip(symbols, anchors)
                )
                contacts = list(record.get("close_contact_pairs") or ())
                lines = [
                    f"parent: {record.get('parent_artifact_id')} "
                    f"({record.get('parent_formula')}, "
                    f"{record.get('parent_atom_count')} atoms)",
                    f"appended {record.get('element')} as atom "
                    f"{record.get('appended_atom_index')}, against {named}",
                    f"bond {_shown_value(record.get('bond_length_angstrom'))} "
                    f"angstrom, angle "
                    f"{_shown_value(record.get('angle_degrees'))} deg, "
                    f"dihedral {_shown_value(record.get('dihedral_degrees'))} "
                    f"deg (achieved "
                    f"{_shown_value(record.get('achieved_bond_length_angstrom'))}, "
                    f"{_shown_value(record.get('achieved_angle_degrees'))}, "
                    f"{_shown_value(record.get('achieved_dihedral_degrees'))})",
                    f"{record.get('parent_formula')} -> "
                    f"{record.get('formula')}; "
                    f"{record.get('fragment_count')} connected piece(s); "
                    f"{record.get('atom_order_note')}",
                    "closest contact "
                    f"{_shown_value(record.get('min_interatomic_distance_angstrom'))}"
                    f" angstrom; {len(contacts)} pair(s) inside covalent "
                    "contact -- observations, not verdicts",
                ]
                title_kind = "appended atom"
            elif record.get("kind") == "mode_displacement":
                moved = list(record.get("moved_atoms") or ())
                leading = list(record.get("leading_atoms") or ())
                contacts = list(record.get("close_contact_pairs") or ())
                character = (
                    "imaginary" if record.get("mode_is_imaginary") else "real"
                )
                lines = [
                    f"from result {record.get('result_artifact_id')} "
                    f"({record.get('program')}; {record.get('formula')}, "
                    f"{record.get('atom_count')} atoms)",
                    f"stepped along mode {record.get('mode_index')} "
                    f"({_shown_value(record.get('mode_frequency_cm_1'))} "
                    f"cm^-1, {character}) by "
                    f"{_shown_value(record.get('amplitude_angstrom'))} "
                    "angstrom (largest atom displacement achieved "
                    f"{_shown_value(record.get('achieved_max_displacement_angstrom'))})",
                    f"{len(moved)} atom(s) moved; leading atoms "
                    f"{', '.join(str(index) for index in leading)}; "
                    f"{record.get('atom_order_note')}",
                    "closest contact "
                    f"{_shown_value(record.get('min_interatomic_distance_angstrom'))}"
                    f" angstrom; {len(contacts)} pair(s) inside covalent "
                    "contact; connectivity "
                    + (
                        "changed"
                        if record.get("connectivity_changed")
                        else "unchanged"
                    )
                    + " -- observations, not verdicts",
                ]
                title_kind = "mode displacement"
            elif record.get("kind") == "symmetry_break":
                contacts = list(record.get("close_contact_pairs") or ())
                lines = [
                    f"parent: {record.get('parent_artifact_id')} "
                    f"({record.get('formula')}, "
                    f"{record.get('atom_count')} atoms)",
                    f"every atom perturbed by seed {record.get('seed')} "
                    "within "
                    f"{_shown_value(record.get('amplitude_angstrom'))} "
                    "angstrom (largest step achieved "
                    f"{_shown_value(record.get('max_displacement_angstrom'))}"
                    ", rms "
                    f"{_shown_value(record.get('rms_displacement_angstrom'))}"
                    ")",
                    f"point group estimate {record.get('point_group_before')}"
                    f" -> {record.get('point_group_after')}; "
                    f"{record.get('atom_order_note')}",
                    "closest contact "
                    f"{_shown_value(record.get('min_interatomic_distance_angstrom'))}"
                    f" angstrom; {len(contacts)} pair(s) inside covalent "
                    "contact; connectivity "
                    + (
                        "changed"
                        if record.get("connectivity_changed")
                        else "unchanged"
                    )
                    + " -- observations, not verdicts",
                ]
                title_kind = "symmetry break"
            else:
                atoms = list(record.get("coordinate_atoms") or ())
                symbols = list(record.get("coordinate_symbols") or ())
                named = " - ".join(
                    f"{symbol}{index}" for symbol, index in zip(symbols, atoms)
                )
                unit = record.get("value_unit")
                moved = list(record.get("moved_atoms") or ())
                contacts = list(record.get("close_contact_pairs") or ())
                lines = [
                    f"parent: {record.get('parent_artifact_id')}",
                    f"{record.get('operation')} on {named}",
                    f"before {_shown_value(record.get('value_before'))} "
                    f"{unit} -> requested "
                    f"{_shown_value(record.get('value_requested'))} {unit} "
                    f"-> achieved "
                    f"{_shown_value(record.get('value_achieved'))} {unit}",
                    f"moved the side of atom "
                    f"{record.get('moving_side_atom')}: {moved} "
                    f"({len(moved)} atoms); every other atom unchanged",
                    f"{record.get('formula')} "
                    f"({record.get('atom_count')} atoms); "
                    f"{record.get('atom_order_note')}",
                    "closest contact "
                    f"{_shown_value(record.get('min_interatomic_distance_angstrom'))}"
                    f" angstrom; {len(contacts)} pair(s) inside covalent "
                    "contact; connectivity "
                    + (
                        "changed"
                        if record.get("connectivity_changed")
                        else "unchanged"
                    )
                    + " -- observations, not verdicts",
                ]
                title_kind = "edited geometry"
            if position == total:
                lines.append(str(record.get("starting_structure_role")))
                lines.append(
                    "electronic state was deliberately unbound while "
                    "building; this node binds charge "
                    f"{item.molecular_identity.get('charge')}, multiplicity "
                    f"{item.molecular_identity.get('multiplicity')} "
                    "explicitly"
                )
            panels.append(
                Panel(
                    "\n".join(lines),
                    title=(
                        f"{item.node_id} · {title_kind} lineage, {step} "
                        "(covered by this approval)"
                    ),
                )
            )
    return panels


def render_review_blocks(
    review: "WorkflowExecutionReviewV1",
    *,
    project_yaml: Mapping[str, str] | None = None,
) -> tuple[Any, ...]:
    """Every renderable of the approval display, in reading order."""

    yaml_texts = project_yaml or {}
    blocks: list[Any] = [_overview_table(review), _environment_table(review)]
    records_table = _batch_records_renderable(review)
    if records_table is not None:
        blocks.append(records_table)
    blocks.append(_bounds_panel(review))
    if review.scientific_plan.edges:
        blocks.append(_edges_table(review))
    blocks.extend(_composition_panels(review))
    blocks.extend(_derivation_panels(review))
    blocks.extend(_pubchem_geometry_panels(review))
    blocks.extend(_database_extraction_panels(review))
    blocks.extend(_geometry_lineage_panels(review))
    blocks.append(_analysis_chain_renderable(review))
    constants_table = _literature_constants_renderable(review)
    if constants_table is not None:
        blocks.append(constants_table)
    declared_panel = _declared_observable_panel(review)
    if declared_panel is not None:
        blocks.append(declared_panel)
    knowledge_panel = _advisory_knowledge_panel(review)
    if knowledge_panel is not None:
        blocks.append(knowledge_panel)
    route_directives = _route_directive_panel(review)
    if route_directives is not None:
        blocks.append(route_directives)
    for item in review.node_reviews:
        yaml_text = yaml_texts.get(item.node_id)
        if yaml_text:
            settings = Syntax(
                yaml_text.rstrip() + "\n",
                "yaml",
                word_wrap=True,
                background_color="default",
            )
            settings_title = f"{item.node_id} · project settings (YAML)"
        else:
            settings = Syntax(
                item.project_settings_text,
                "json",
                word_wrap=True,
                background_color="default",
            )
            settings_title = f"{item.node_id} · effective project settings"
        blocks.append(Panel(settings, title=settings_title))
        blocks.append(
            Panel(
                Text(human_cli_operation(item.real_execution_argv)),
                title=f"{item.node_id} · ChemSmart CLI operation",
            )
        )
        stated = tuple(
            line
            for entry in getattr(review, "node_observations", ())
            if entry.get("node_id") == item.node_id
            for line in entry.get("observations", ())
        )
        if stated:
            blocks.append(
                Panel(
                    Text("\n".join(f"- {line}" for line in stated)),
                    title=f"{item.node_id} · host observations",
                )
            )
    blocks.append(_decision_panel(review))
    return tuple(blocks)


__all__ = ["render_review_blocks", "resolve_project_yaml_texts"]
