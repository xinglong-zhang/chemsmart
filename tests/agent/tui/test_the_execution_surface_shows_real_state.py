"""The waiting UI shows what really runs, and the DAG shows every status.

States come from host events and typed results only. The mermaid artifact
is a deterministic projection of the reviewed plan, written beside the
run's evidence.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.mermaid import render_workflow_mermaid  # noqa: E402
from chemsmart.agent.tui.panels import dag_rows, state_glyph  # noqa: E402
from chemsmart.agent.tui.runs import list_runs  # noqa: E402


def _flatten(renderable) -> str:
    console = Console(record=True, width=200)
    console.print(renderable)
    return " ".join(console.export_text().split())


def test_dag_rows_show_every_status_with_minimal_glyphs():
    nodes = {
        "opt": {
            "kind": "calc",
            "label": "xtb opt · H2O",
            "state": "validated",
        },
        "sp": {"kind": "calc", "label": "xtb sp · H2O", "state": "running"},
        "irc": {
            "kind": "calc",
            "label": "orca irc",
            "state": "deferred",
            "detail": "IRC execution is not release-qualified",
        },
        "extract": {
            "kind": "analysis",
            "label": "result_extraction",
            "state": "queued",
        },
        "verdict": {
            "kind": "analysis",
            "label": "scientific_validation",
            "state": "failed",
        },
    }
    text = _flatten(dag_rows(nodes))

    assert "[✓] opt" in text
    assert "[•] sp" in text
    assert "⊘ irc" in text
    assert "↳ IRC execution is not release-qualified" in text
    assert "[ ] extract" in text
    assert "[✗] verdict" in text
    assert state_glyph("engine_complete")[0] == "[»]"


def test_the_mermaid_artifact_projects_the_reviewed_plan():
    calc = SimpleNamespace(
        node_id="opt-water", stage="opt", program="xtb"
    )
    consumer = SimpleNamespace(node_id="sp-water", stage="sp", program="xtb")
    review = SimpleNamespace(
        scientific_plan=SimpleNamespace(
            nodes=(calc, consumer),
            edges=(
                SimpleNamespace(
                    source_node_id="opt-water",
                    target_node_id="sp-water",
                    artifact_class="geometry_xyz",
                ),
            ),
        ),
        node_reviews=(
            SimpleNamespace(node_id="opt-water"),
            SimpleNamespace(node_id="sp-water"),
        ),
        non_executable_node_ids=(),
        scientific_toolchain_plan=SimpleNamespace(
            analysis_nodes=(
                SimpleNamespace(
                    node_id="extract-e",
                    analysis_kind="result_extraction",
                    dependencies=("sp-water",),
                    inputs=(
                        SimpleNamespace(producer_node_id="sp-water"),
                    ),
                ),
            )
        ),
    )

    mermaid = render_workflow_mermaid(review)

    assert mermaid.startswith("flowchart TD")
    assert 'opt_water["opt-water<br/>xtb opt"]' in mermaid
    assert "opt_water -->|geometry_xyz| sp_water" in mermaid
    assert 'extract_e(["extract-e<br/>result_extraction"])' in mermaid
    assert "sp_water -.-> extract_e" in mermaid


def test_runs_are_listed_from_their_durable_evidence(tmp_path: Path):
    import json

    execution = tmp_path / ".chemsmart-agent" / "executions" / "tui-abc"
    execution.mkdir(parents=True)
    (execution / "events.jsonl").write_text(
        json.dumps(
            {
                "kind": "runtime_terminated",
                "payload": {"terminal_state": "complete"},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    report_dir = execution / "analysis"
    report_dir.mkdir()
    (report_dir / "completed-analysis-report.md").write_text(
        "# report\n", encoding="utf-8"
    )
    replay = tmp_path / ".chemsmart-agent" / "replays" / "replay-1" / "run"
    replay.mkdir(parents=True)
    (replay / "events.jsonl").write_text("", encoding="utf-8")

    summaries = list_runs(tmp_path)

    by_name = {summary.name: summary for summary in summaries}
    assert by_name["tui-abc"].terminal_state == "complete"
    assert by_name["tui-abc"].report_path is not None
    assert by_name["replay-1"].kind == "replay"
    assert by_name["replay-1"].terminal_state == "in progress"
