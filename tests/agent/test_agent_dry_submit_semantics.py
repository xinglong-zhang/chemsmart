from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok

REQUEST = "single-point energy on water.xyz at B3LYP/6-31G(d) Gaussian"


def _single_point_plan(filepath: str, include_submit: bool) -> dict:
    steps = [
        {"tool": "build_molecule", "args": {"filepath": filepath}},
        {
            "tool": "build_gaussian_settings",
            "args": {"functional": "b3lyp", "basis": "6-31g(d)"},
        },
        {
            "tool": "build_job",
            "args": {
                "kind": "gaussian.sp",
                "molecule": "$step1",
                "settings": "$step2",
                "label": "sp_case",
            },
        },
        {"tool": "dry_run_input", "args": {"job": "$step3"}},
        {
            "tool": "validate_runtime",
            "args": {"job": "$step3", "server": None},
        },
    ]
    if include_submit:
        steps.append({"tool": "submit_hpc", "args": {"job": "$step3"}})
    return {"steps": steps, "rationale": "Prepare a SP workflow."}


def _patch_runtime_partial(monkeypatch) -> None:
    import chemsmart.agent.tools as agent_tools

    monkeypatch.setattr(
        agent_tools,
        "validate_runtime",
        lambda job, server=None: {
            "ok": "partial",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [
                "server.queue required",
                "server.account required",
                "server.scratch_dir required",
                "server.modules_or_executable_path required",
            ],
        },
    )


def _dry_submit_provider(filepath: str) -> FakeProvider:
    return FakeProvider(
        [
            _single_point_plan(filepath, include_submit=True),
            {
                "verdict": "warn",
                "confidence": 0.55,
                "issues": [
                    "geometry handoff missing",
                    "server.queue required",
                ],
                "rationale": "Flagged geometry handoff and server setup.",
            },
        ]
    )


def test_dry_submit_single_point_exits_clean(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    _patch_runtime_partial(monkeypatch)

    def fake_submit_hpc(job, server=None, transport=None, execute=False):
        raise ValueError(
            "submit_hpc requires server when no configured servers are "
            "available"
        )

    monkeypatch.setattr(agent_tools, "submit_hpc", fake_submit_hpc)

    session = AgentSession(
        provider=_dry_submit_provider(single_molecule_xyz_file),
        registry=ToolRegistry.default(),
        session_root=tmp_path / "session",
    )
    result = session.run(REQUEST, dry_submit=True)

    assert result["blocked"] is False
    assert result["critic_verdict"].verdict == "ok"
    assert result["critic_verdict"].issues == []
    assert result["preview_submit"]["skipped"] is True

    provider = _dry_submit_provider(single_molecule_xyz_file)
    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", lambda: provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", lambda: provider)
    monkeypatch.setattr(
        "chemsmart.agent.core._default_session_root",
        lambda: str(tmp_path / "cli-sessions"),
    )

    cli_result = CliRunner().invoke(
        agent,
        ["run", "--dry-submit", REQUEST],
        catch_exceptions=False,
    )

    assert cli_result.exit_code == 0, cli_result.output
    assert "Error: critic gating blocked execution" not in cli_result.output


def test_execute_still_blocks_when_server_missing(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    _patch_runtime_partial(monkeypatch)

    session = AgentSession(
        provider=FakeProvider(
            [
                _single_point_plan(
                    single_molecule_xyz_file, include_submit=False
                ),
                critic_ok(),
            ]
        ),
        registry=ToolRegistry.default(),
        session_root=tmp_path / "execute",
    )
    result = session.run(REQUEST, dry_submit=False)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "warn"
    assert "server.queue required" in result["critic_verdict"].issues
