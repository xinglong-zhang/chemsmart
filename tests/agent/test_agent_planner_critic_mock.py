from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from chemsmart.agent.core import AgentSession, Plan
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def _opt_freq_plan(filepath: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Load the structure.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Use a simple DFT level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.opt",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "opt_case",
                },
                "rationale": "Build the optimization job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Inspect the optimization input.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.freq",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "freq_case",
                },
                "rationale": "Build the frequency job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step5"},
                "rationale": "Inspect the frequency input.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step5", "server": None},
                "rationale": "Check runtime prerequisites.",
            },
        ],
        "rationale": "Prepare optimization and frequency dry runs.",
        "estimated_cost": "low",
    }


def _gaussian_irc_plan(filepath: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Load the structure.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Use a simple DFT level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.irc",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "irc_case",
                },
                "rationale": "Build the IRC job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Inspect the IRC input.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Check runtime prerequisites.",
            },
        ],
        "rationale": "Prepare a Gaussian IRC workflow.",
        "estimated_cost": "low",
    }


@pytest.mark.parametrize(
    ("user_request", "label"),
    [
        ("optimize examples/h2o.xyz", "plan_h2o"),
        ("optimize then frequency examples/co2.xyz", "plan_co2"),
        ("optimize methane with Gaussian", "plan_ch4"),
    ],
)
def test_planner_emits_valid_plan_for_canonical_requests(
    single_molecule_xyz_file,
    tmp_path: Path,
    user_request: str,
    label: str,
):
    provider = FakeProvider([planner_plan(single_molecule_xyz_file, label)])
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )

    plan = session._planner_call(user_request)

    assert isinstance(plan, Plan)
    assert plan.steps[0].tool == "build_molecule"
    assert plan.steps[-1].tool == "validate_runtime"


def test_critic_blocks_malformed_input(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "malformed_case", "run_local"
            ),
            critic_ok(),
        ]
    )

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_dry_run_input(job):
        return {
            "inputfile": str(Path(job.folder) / "malformed.com"),
            "content": "%chk=bad.chk\nmissing route line\n",
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "dry_run_input", fake_dry_run_input)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("optimize malformed input", dry_submit=True)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "warn"
    assert result["critic_verdict"].confidence == pytest.approx(0.95)
    assert (
        "Gaussian route line missing or malformed"
        in result["critic_verdict"].issues
    )


def test_critic_refuses_partial_without_allow_remote_unknown(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [
            planner_plan(
                single_molecule_xyz_file, "partial_case", "run_local"
            ),
            critic_ok(),
        ]
    )
    run_local_calls = {"count": 0}

    def fake_run_local(job):
        run_local_calls["count"] += 1
        return {
            "ok": True,
            "returncode": 0,
            "stdout_path": str(Path(job.folder) / "run.stdout"),
            "stderr_path": str(Path(job.folder) / "run.stderr"),
            "output_summary": {},
        }

    monkeypatch.setattr(agent_tools, "run_local", fake_run_local)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("optimize with partial runtime", dry_submit=True)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "warn"
    assert run_local_calls["count"] == 0
    assert "server.queue required" in result["critic_verdict"].issues


def test_critic_receives_all_dry_run_inputs(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [_opt_freq_plan(single_molecule_xyz_file), critic_ok()]
    )

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("optimize then frequency", dry_submit=True)

    assert result["blocked"] is False
    assert len(result["dry_run_results"]) == 2

    critic_payload = json.loads(provider.calls[1]["messages"][1]["content"])
    assert "dry_run_inputs" in critic_payload
    assert len(critic_payload["dry_run_inputs"]) == 2
    assert critic_payload["dry_run_inputs"][0]["inputfile"].endswith(
        "opt_case.com"
    )
    assert critic_payload["dry_run_inputs"][1]["inputfile"].endswith(
        "freq_case.com"
    )


def test_critic_hard_rejects_irc_without_irc_keyword(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [_gaussian_irc_plan(single_molecule_xyz_file), critic_ok()]
    )

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    def fake_dry_run_input(job):
        return {
            "inputfile": str(Path(job.folder) / "irc_case.com"),
            "content": "#p b3lyp/6-31g* opt\n\nirc test\n\n0 1\nO 0.0 0.0 0.0\n",
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)
    monkeypatch.setattr(agent_tools, "dry_run_input", fake_dry_run_input)

    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path,
    )
    result = session.run("run IRC without irc keyword", dry_submit=True)

    assert result["blocked"] is True
    assert result["critic_verdict"].verdict == "reject"
    assert (
        "Gaussian IRC input missing irc= keyword"
        in result["critic_verdict"].issues
    )


def test_identical_queries_produce_identical_plan_hashes(
    single_molecule_xyz_file,
    tmp_path: Path,
):
    query = "optimize examples/h2o.xyz"
    explicit_defaults_plan = {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": single_molecule_xyz_file, "index": "-1"},
                "rationale": "Load the neutral closed-shell structure.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {
                    "functional": "B3LYP",
                    "basis": "6-31G*",
                    "charge": 0,
                    "multiplicity": 1,
                    "freq": False,
                    "numfreq": False,
                },
                "rationale": "B3LYP/6-31G* is a cost-effective starting point for a small closed-shell neutral molecule.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.opt",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "hash_case",
                },
                "rationale": "A geometry optimization is appropriate because the request is to relax the structure before any higher-level follow-up.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "The dry run verifies the chosen route line before any execution.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Runtime validation confirms the local prerequisites for this small Gaussian optimization.",
            },
        ],
        "rationale": "Prepare a Gaussian optimization dry run for a small neutral molecule.",
        "estimated_cost": "low",
    }
    omitted_defaults_plan = {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": single_molecule_xyz_file},
                "rationale": "Load the neutral closed-shell structure.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "B3LYP", "basis": "6-31G*"},
                "rationale": "B3LYP/6-31G* is a cost-effective starting point for a small closed-shell neutral molecule.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.opt",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": "hash_case",
                },
                "rationale": "A geometry optimization is appropriate because the request is to relax the structure before any higher-level follow-up.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "The dry run verifies the chosen route line before any execution.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3"},
                "rationale": "Runtime validation confirms the local prerequisites for this small Gaussian optimization.",
            },
        ],
        "rationale": "Prepare a Gaussian optimization dry run for a small neutral molecule.",
        "estimated_cost": "low",
    }

    plan_hashes = []
    for response in (explicit_defaults_plan, omitted_defaults_plan):
        provider = FakeProvider([response])
        session = AgentSession(
            provider=provider,
            registry=ToolRegistry.default(),
            session_root=tmp_path,
        )
        plan = session._planner_call(query)
        payload = json.dumps(plan.model_dump(), sort_keys=True)
        plan_hashes.append(hashlib.sha256(payload.encode("utf-8")).hexdigest())

    assert plan_hashes[0] == plan_hashes[1]
