"""A chemist reads project settings as YAML, not canonical JSON.

Promotion enforces that a node's project artifact digest equals the digest
of its rendered YAML bytes, so the display can resolve the readable text
exactly: from the session transcript first, from the promoted files under
the planning runs second, and only a node that resolves nowhere keeps the
canonical-settings fallback.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.review import (  # noqa: E402
    render_review_blocks,
    resolve_project_yaml_texts,
)

_YAML = "project:\n  functional: gfn2\n  jobtype: sp\n"
_SHA = hashlib.sha256(_YAML.encode()).hexdigest()


def _review():
    node = SimpleNamespace(
        node_id="sp-water",
        program="xtb",
        engine="cpu",
        stage="sp",
        project_artifact_sha256=_SHA,
        molecular_identity={
            "approved_names": ("water",),
            "formula": "H2O",
            "atom_order": ("O", "H", "H"),
            "charge": 0,
            "multiplicity": 1,
        },
        environment_summary={"status": "available"},
        project_settings_text='{"canonical": true}',
        real_execution_argv=("chemsmart", "run", "xtb", "sp"),
    )
    return SimpleNamespace(
        execution_resources=SimpleNamespace(
            cores=4, memory_gb=8.0, gpu_count=0, node_timeout_seconds=600
        ),
        execution_envelope={"max_engine_calls": 1},
        node_reviews=(node,),
        non_executable_node_ids=(),
        scientific_plan=SimpleNamespace(
            nodes=(node,), edges=(), plan_sha256="b" * 64
        ),
        scientific_toolchain_plan=None,
        review_sha256="a" * 64,
    )


def _transcript_with_render():
    record = {
        "schema_version": "chemsmart.tool-result.v1",
        "tool": "render_project_yaml",
        "status": "ok",
        "result": {"rendered_yaml": _YAML, "rendered_sha256": _SHA},
    }
    return (
        {
            "role": "assistant",
            "content": "",
            "tool_calls": [
                {
                    "id": "c1",
                    "type": "function",
                    "function": {
                        "name": "render_project_yaml",
                        "arguments": "{}",
                    },
                }
            ],
        },
        {"role": "tool", "tool_call_id": "c1", "content": json.dumps(record)},
    )


def test_the_transcript_route_resolves_the_yaml():
    texts = resolve_project_yaml_texts(
        _review(), public_transcript=_transcript_with_render()
    )
    assert texts == {"sp-water": _YAML}


def test_the_promoted_file_is_the_fallback(tmp_path: Path):
    projects = tmp_path / ".chemsmart-agent" / "runs" / "live-x" / "projects"
    projects.mkdir(parents=True)
    (projects / "sp-water.yaml").write_text(_YAML, encoding="utf-8")

    texts = resolve_project_yaml_texts(_review(), workspace=tmp_path)

    assert texts == {"sp-water": _YAML}


def test_an_unresolvable_node_keeps_the_canonical_fallback():
    review = _review()
    texts = resolve_project_yaml_texts(review)
    assert texts == {}

    console = Console(record=True, width=200)
    for block in render_review_blocks(review, project_yaml=texts):
        console.print(block)
    text = console.export_text()
    assert "effective project settings" in text


def test_the_yaml_panel_is_titled_and_rendered():
    review = _review()
    texts = resolve_project_yaml_texts(
        review, public_transcript=_transcript_with_render()
    )
    console = Console(record=True, width=200)
    for block in render_review_blocks(review, project_yaml=texts):
        console.print(block)
    text = console.export_text()
    assert "project settings (YAML)" in text
    assert "functional: gfn2" in text
    assert "effective project settings" not in text
