"""The terminal shell: explicit approval chain, honest guards, live chrome."""

from __future__ import annotations

import asyncio
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
)


def _app(tmp_path: Path) -> ChemSmartAgentApp:
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    envelope = tmp_path / "execution-envelope.json"
    envelope.write_text("{}\n", encoding="utf-8")
    return ChemSmartAgentApp(
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


def _transcript_text(app: ChemSmartAgentApp) -> str:
    transcript = app.query_one("#transcript")
    console = Console(record=True, width=200)
    for renderable in transcript.recorder.entries:
        console.print(renderable)
    return " ".join(
        console.export_text().replace("│", " ").split()
    )


async def _submit(app, pilot, text: str) -> None:
    composer = app.query_one("#composer")
    composer.value = text
    composer.focus()
    await pilot.press("enter")
    await pilot.pause()


def test_tui_launches_and_exposes_explicit_approval_chain(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            await _submit(app, pilot, "/help")
            text = _transcript_text(app)
            assert "/approve" in text
            assert "provider-free ChemSmart executor" in text
            assert app.controller.phase.value == "ready"

    asyncio.run(observe())


def test_a_mistyped_command_is_answered_with_the_nearest_real_one(
    tmp_path: Path,
):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            await _submit(app, pilot, "/aprove")
            text = _transcript_text(app)
            assert "Unknown command: /aprove" in text
            assert "Did you mean /approve?" in text

    asyncio.run(observe())


def test_approve_refuses_before_any_banner_when_nothing_is_reviewed(
    tmp_path: Path,
):
    """The guard speaks first; an approval banner over nothing is a lie."""

    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            await _submit(app, pilot, "/approve")
            text = _transcript_text(app)
            assert "Approve and run" not in text
            assert "Approval stopped" in text
            assert "review the displayed ChemSmart workflow" in text

    asyncio.run(observe())


def test_revise_declines_and_hands_the_task_back_for_editing(tmp_path: Path):
    app = _app(tmp_path)

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            app.controller.task = "optimize the water dimer with xtb"
            await _submit(app, pilot, "/revise")
            composer = app.query_one("#composer")
            assert composer.value == "optimize the water dimer with xtb"
            text = _transcript_text(app)
            assert "revision requested" in text

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

    flattened = " ".join(transcript.replace("│", " ").split())
    assert "irc" in flattened
    assert "Deferred" in flattened
    assert "IRC execution is not release-qualified" in flattened
    assert (
        "Deferred stages remain unapproved and will not launch: irc"
        in flattened
    )
