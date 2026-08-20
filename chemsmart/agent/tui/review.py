"""The human-approval display: everything the single /approve covers.

Charter step 3: the terminal shows the complete plan -- molecular identity,
electronic state, effective settings, CLI operations, dependencies,
environment, resources -- and, because the one approval also covers them,
the typed analysis chain carried with the workflow and the full lineage of
any composed arrangement. Every panel is a projection of the same typed
review object; ``render_raw_review`` shows the identical canonical bytes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Mapping

from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text

from chemsmart.agent._contracts import canonical_data, canonical_json

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1

from .presentation import human_cli_operation


def _short(digest: str, length: int = 12) -> str:
    return digest[:length] + "…" if digest else ""


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
            execution_state = "Executable"
            reason = "reviewed below"
        else:
            molecule = planned.project_role
            program_engine = f"{planned.program} / {planned.engine}"
            stage = planned.stage
            if (
                planned.charge is not None
                and planned.multiplicity is not None
            ):
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
    return Panel(
        f"cores: {resources.cores}\n"
        f"memory: {resources.memory_gb:g} GB\n"
        f"GPUs: {resources.gpu_count}\n"
        f"node timeout: {resources.node_timeout_seconds} s\n"
        f"engine-call limit: "
        f"{review.execution_envelope.get('max_engine_calls')}",
        title="Execution bounds",
    )


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
            f"parent A: {record.get('fragment_a_artifact_id')} "
            f"(geometry {_short(str(record.get('fragment_a_sha256') or ''))}, "
            f"identity {_short(str(record.get('fragment_a_identity_sha256') or ''))})",
            f"parent B: {record.get('fragment_b_artifact_id')} "
            f"(geometry {_short(str(record.get('fragment_b_sha256') or ''))}, "
            f"identity {_short(str(record.get('fragment_b_identity_sha256') or ''))})",
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
            f"evidence status: {status}",
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
    table.add_column("Outputs")
    table.add_column("Conditions")
    table.add_column("State")
    for node in plan.analysis_nodes:
        inputs = ", ".join(
            f"{item.producer_node_id}.{item.producer_output_id}"
            for item in node.inputs
        )
        outputs = ", ".join(
            f"{item.output_id} ({item.unit})" for item in node.outputs
        )
        conditions = []
        if node.temperature_k is not None:
            conditions.append(f"{node.temperature_k:g} K")
        if node.pressure_atm is not None:
            conditions.append(f"{node.pressure_atm:g} atm")
        state = node.support_state
        if node.blocked_reason:
            state += f": {node.blocked_reason}"
        table.add_row(
            node.node_id,
            node.analysis_kind,
            inputs or "—",
            outputs or "—",
            ", ".join(conditions) or "—",
            state,
        )
    return table


def _digest_line(review: "WorkflowExecutionReviewV1") -> Text:
    return Text(
        f"review sha256 {_short(review.review_sha256, 16)} · plan sha256 "
        f"{_short(review.scientific_plan.plan_sha256, 16)} · /raw shows the "
        "canonical bytes this display projects",
        style="dim",
    )


def _decision_panel(review: "WorkflowExecutionReviewV1") -> Panel:
    executable_text = ", ".join(
        item.node_id for item in review.node_reviews
    )
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
    decision += " A changed scientific request must be replanned."
    return Panel(decision, title="Human decision", border_style="yellow")


def render_review_blocks(
    review: "WorkflowExecutionReviewV1",
) -> tuple[Any, ...]:
    """Every renderable of the approval display, in reading order."""

    blocks: list[Any] = [_overview_table(review), _environment_table(review)]
    blocks.append(_bounds_panel(review))
    if review.scientific_plan.edges:
        blocks.append(_edges_table(review))
    blocks.extend(_composition_panels(review))
    blocks.append(_analysis_chain_renderable(review))
    for item in review.node_reviews:
        blocks.append(
            Panel(
                Syntax(
                    item.project_settings_text,
                    "json",
                    word_wrap=True,
                    background_color="default",
                ),
                title=f"{item.node_id} · effective project settings",
            )
        )
        blocks.append(
            Panel(
                Text(human_cli_operation(item.real_execution_argv)),
                title=f"{item.node_id} · ChemSmart CLI operation",
            )
        )
    blocks.append(_digest_line(review))
    blocks.append(_decision_panel(review))
    return tuple(blocks)


def render_raw_review(review: "WorkflowExecutionReviewV1") -> Panel:
    """The canonical bytes behind the humanized panels, unchanged."""

    return Panel(
        Syntax(
            canonical_json({"workflow_execution_review": canonical_data(review)}),
            "json",
            word_wrap=True,
            background_color="default",
        ),
        title=f"Canonical review · sha256 {review.review_sha256}",
    )


__all__ = ["render_raw_review", "render_review_blocks"]
