"""Host-authored mermaid projection of an approved workflow.

Written into the run directory at approve time so the DAG the human saw in
the terminal also exists as a file artifact next to the run's evidence. It
is a projection of the reviewed plan bytes -- deterministic, no styling
beyond executable/deferred/analysis classes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from chemsmart.agent.execution import WorkflowExecutionReviewV1


def _identifier(node_id: str) -> str:
    return "".join(
        char if char.isalnum() or char == "_" else "_" for char in node_id
    )


def render_workflow_mermaid(review: "WorkflowExecutionReviewV1") -> str:
    lines = ["flowchart TD"]
    deferred = set(review.non_executable_node_ids)
    reviewed = {item.node_id for item in review.node_reviews}
    for node in review.scientific_plan.nodes:
        label = f"{node.node_id}<br/>{node.program} {node.stage}"
        shape = f'{_identifier(node.node_id)}["{label}"]'
        lines.append(f"    {shape}")
        if node.node_id in deferred or (node.node_id not in reviewed):
            lines.append(f"    class {_identifier(node.node_id)} deferred")
    for edge in review.scientific_plan.edges:
        arrow = f"-->|{edge.artifact_class}|" if edge.artifact_class else "-->"
        lines.append(
            f"    {_identifier(edge.source_node_id)} {arrow} "
            f"{_identifier(edge.target_node_id)}"
        )
    plan = review.scientific_toolchain_plan
    if plan is not None:
        calculation_ids = {
            node.node_id for node in review.scientific_plan.nodes
        }
        for node in plan.analysis_nodes:
            label = f"{node.node_id}<br/>{node.analysis_kind}"
            lines.append(f'    {_identifier(node.node_id)}(["{label}"])')
            lines.append(f"    class {_identifier(node.node_id)} analysis")
            producers = {
                item.producer_node_id
                for item in node.inputs
                if item.producer_node_id
            } | set(node.dependencies)
            for producer in sorted(producers):
                if producer in calculation_ids or any(
                    other.node_id == producer for other in plan.analysis_nodes
                ):
                    lines.append(
                        f"    {_identifier(producer)} -.-> "
                        f"{_identifier(node.node_id)}"
                    )
    lines.append("    classDef deferred stroke-dasharray: 5 5,opacity:0.6;")
    lines.append("    classDef analysis stroke-width:1px,opacity:0.85;")
    return "\n".join(lines) + "\n"


__all__ = ["render_workflow_mermaid"]
