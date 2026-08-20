"""A terminal-level smoke observation for the production Textual shell."""

from __future__ import annotations

import asyncio
from pathlib import Path
from types import SimpleNamespace

import pytest
from rich.console import Console

pytest.importorskip("textual")

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
)


def test_tui_launches_and_exposes_explicit_approval_chain(tmp_path: Path):
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    envelope = tmp_path / "execution-envelope.json"
    envelope.write_text("{}\n", encoding="utf-8")
    app = ChemSmartAgentApp(
        AgentTuiController(
            AgentSessionConfigV1(
                workspace=tmp_path,
                secret_file=secret,
                execution_envelope_file=envelope,
                review_file=tmp_path / "review-output.json",
            )
        ),
        plain=True,
    )

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            composer = app.query_one("#composer")
            composer.value = "/help"
            await pilot.press("enter")
            await pilot.pause()

            transcript = app.query_one("#transcript")
            assert len(transcript.lines) > 10
            assert app.controller.phase.value == "ready"

    asyncio.run(observe())


def test_approval_surface_displays_deferred_stage_and_reason(tmp_path: Path):
    secret = tmp_path / "unused.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    envelope = tmp_path / "unused.json"
    envelope.write_text("{}\n", encoding="utf-8")
    app = ChemSmartAgentApp(
        AgentTuiController(
            AgentSessionConfigV1(
                workspace=tmp_path,
                secret_file=secret,
                execution_envelope_file=envelope,
                review_file=tmp_path / "review.json",
            )
        ),
        plain=True,
    )
    rendered = []
    app._write = rendered.append  # type: ignore[method-assign]
    executable = SimpleNamespace(
        node_id="ts",
        program="orca",
        engine="cpu",
        stage="ts",
        molecular_identity={
            "approved_names": ("reaction transition structure",),
            "formula": "C2H2",
            "atom_order": ("C", "C", "H", "H"),
            "charge": 0,
            "multiplicity": 1,
        },
        environment_summary={
            "status": "available",
            "target_kind": "executable",
            "observed_version": "6.1.1",
            "observation_method": "host probe",
        },
        project_settings_text="{}",
        real_execution_argv=("chemsmart", "run", "orca", "ts"),
    )
    deferred = SimpleNamespace(
        node_id="irc",
        program="orca",
        engine="cpu",
        stage="irc",
        project_role="reaction-path",
        charge=0,
        multiplicity=1,
        blocked_reason="IRC execution is not release-qualified",
    )
    review = SimpleNamespace(
        execution_resources=SimpleNamespace(
            cores=4,
            memory_gb=8.0,
            gpu_count=0,
            node_timeout_seconds=600,
        ),
        execution_envelope={"max_engine_calls": 1},
        node_reviews=(executable,),
        non_executable_node_ids=("irc",),
        scientific_plan=SimpleNamespace(nodes=(executable, deferred), edges=()),
    )

    app._present_request(review)
    console = Console(record=True, width=220)
    for item in rendered:
        console.print(item)
    transcript = console.export_text()

    flattened = " ".join(transcript.replace("\u2502", " ").split())
    assert "irc" in flattened
    assert "Deferred" in flattened
    assert "IRC execution is not release-qualified" in flattened
    assert (
        "Deferred stages remain unapproved and will not launch: irc"
        in flattened
    )
