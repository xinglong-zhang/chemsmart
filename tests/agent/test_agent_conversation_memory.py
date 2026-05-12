from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.services.conversation_memory import (
    ConversationMemory,
    EntityMemory,
)

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def _single_point_plan(filepath: str, label: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Reload the same molecule for the follow-up.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Reuse the same lightweight DFT level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.sp",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": label,
                },
                "rationale": "Prepare the single-point follow-up job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Preview the single-point input file.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Check local prerequisites before execution.",
            },
        ],
        "rationale": "Follow the earlier setup with a single-point job.",
        "estimated_cost": "low",
    }


def _request_entry(request: str) -> dict:
    return {"kind": "request", "payload": {"request": request}}


def _tool_use_entries(
    *,
    step: int,
    tool: str,
    args: dict,
    result: dict,
    status: str = "ok",
) -> list[dict]:
    return [
        {
            "kind": "tool_use_request",
            "payload": {
                "step": step,
                "tool": tool,
                "args": args,
            },
        },
        {
            "kind": "tool_use_result",
            "payload": {
                "step": step,
                "tool": tool,
                "status": status,
                "payload": result,
            },
        },
    ]


def _tool_result_entries(
    *,
    step_index: int,
    tool: str,
    args: dict,
    result: dict,
) -> list[dict]:
    return [
        {
            "kind": "tool_call",
            "payload": {
                "step_index": step_index,
                "tool": tool,
                "args": args,
            },
        },
        {
            "kind": "tool_result",
            "payload": {
                "step_index": step_index,
                "tool": tool,
                "payload": result,
            },
        },
    ]


def test_resumed_session_rebuilds_memory_from_decision_log(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    first_request = f"build water from {single_molecule_xyz_file}"
    second_request = "now run a single-point on it"
    session_root = tmp_path / "sessions"

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)

    first_session = AgentSession(
        provider=FakeProvider(
            [planner_plan(single_molecule_xyz_file, "memory_opt"), critic_ok()]
        ),
        registry=ToolRegistry.default(),
        session_root=session_root,
    )
    first_result = first_session.run(first_request, dry_submit=True)

    resumed_provider = FakeProvider(
        [
            _single_point_plan(single_molecule_xyz_file, "memory_sp"),
            critic_ok(),
        ]
    )
    resumed_session = AgentSession.load(
        first_result["session_id"],
        provider=resumed_provider,
        registry=ToolRegistry.default(),
        session_root=session_root,
    )

    assert [
        turn.request for turn in resumed_session.conversation_history.turns
    ] == [first_request]

    second_result = resumed_session.run(second_request, dry_submit=True)
    planner_payload = json.loads(
        resumed_provider.calls[0]["messages"][1]["content"]
    )
    history = planner_payload["conversation_history"]["recent_turns"]

    assert second_result["session_id"] == first_result["session_id"]
    assert history[0]["request"] == first_request
    assert any(
        "dry_run_input wrote" in item
        for item in history[0]["reusable_results"]
    )


def test_conversation_memory_tracks_ssh_probe_scheduler_entities():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("detect the scheduler"),
            *_tool_use_entries(
                step=1,
                tool="ssh_probe",
                args={
                    "server": "cluster-a",
                    "probe_name": "scheduler.detect_kind",
                },
                result={"scheduler": "slurm"},
            ),
        ]
    )

    assert memory.entities.last_server == "cluster-a"
    assert memory.entities.last_scheduler == "slurm"


def test_conversation_memory_tracks_scheduler_query_entities():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("check job state"),
            *_tool_use_entries(
                step=1,
                tool="scheduler_query",
                args={"server": "cluster-a", "job_id": "12345"},
                result={
                    "scheduler": "slurm",
                    "job_id": "12345",
                    "state": "RUNNING",
                    "queue": "debug",
                },
            ),
        ]
    )

    assert memory.entities.last_server == "cluster-a"
    assert memory.entities.last_scheduler == "slurm"
    assert memory.entities.last_job_id == "12345"


def test_conversation_memory_tracks_log_tail_entities_and_error():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("tail the remote log"),
            *_tool_use_entries(
                step=1,
                tool="log_tail",
                args={"server": "cluster-a", "path": "/scratch/job.out"},
                result={
                    "lines_returned": 12,
                    "errors": [
                        {
                            "kind": "oom_killed",
                            "line": "Detected 1 oom-kill event(s)",
                        }
                    ],
                },
            ),
        ]
    )

    assert memory.entities.last_server == "cluster-a"
    assert memory.entities.last_log_path == "/scratch/job.out"
    assert memory.entities.last_probe_error == "Detected 1 oom-kill event(s)"


def test_conversation_memory_read_log_path_updates_last_log_path():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("read the log"),
            *_tool_use_entries(
                step=1,
                tool="read",
                args={"path": "/tmp/foo.log"},
                result={
                    "path": "/tmp/foo.log",
                    "start_line": 10,
                    "end_line": 16,
                    "total_lines": 200,
                },
            ),
        ]
    )

    assert memory.entities.last_log_path == "/tmp/foo.log"


def test_conversation_memory_read_non_log_path_does_not_update_last_log_path():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("read the notes"),
            *_tool_use_entries(
                step=1,
                tool="read",
                args={"path": "/tmp/foo.txt"},
                result={
                    "path": "/tmp/foo.txt",
                    "start_line": 1,
                    "end_line": 4,
                    "total_lines": 20,
                },
            ),
        ]
    )

    assert memory.entities.last_log_path is None


def test_conversation_memory_preserves_prior_server_when_next_turn_omits_it():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("check a queued job"),
            *_tool_use_entries(
                step=1,
                tool="scheduler_query",
                args={"server": "cluster-a", "job_id": "12345"},
                result={"scheduler": "slurm", "job_id": "12345"},
            ),
            _request_entry("read the latest log"),
            *_tool_use_entries(
                step=1,
                tool="log_tail",
                args={"path": "/scratch/job.out"},
                result={"lines_returned": 8, "errors": []},
            ),
        ]
    )

    assert memory.entities.last_server == "cluster-a"
    assert memory.entities.last_log_path == "/scratch/job.out"


def test_prompt_context_includes_entities_only_when_populated():
    empty_context = ConversationMemory().prompt_context()
    populated_context = ConversationMemory(
        entities=EntityMemory(last_server="cluster-a", last_scheduler="slurm")
    ).prompt_context(token_budget=1)

    assert "entities" not in empty_context
    assert populated_context["entities"] == {
        "last_server": "cluster-a",
        "last_scheduler": "slurm",
    }


def test_conversation_memory_reusable_results_include_new_mva_summaries():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("inspect remote state"),
            *_tool_use_entries(
                step=1,
                tool="ssh_probe",
                args={
                    "server": "cluster-a",
                    "probe_name": "scheduler.detect_kind",
                },
                result={"scheduler": "slurm"},
            ),
            *_tool_use_entries(
                step=2,
                tool="scheduler_query",
                args={"server": "cluster-a", "job_id": "12345"},
                result={
                    "scheduler": "slurm",
                    "job_id": "12345",
                    "state": "RUNNING",
                    "queue": "debug",
                },
            ),
            *_tool_use_entries(
                step=3,
                tool="log_tail",
                args={"server": "cluster-a", "path": "/scratch/job.err"},
                result={
                    "path": "/scratch/job.err",
                    "lines_returned": 12,
                    "errors": [
                        {
                            "kind": "oom_killed",
                            "line": "Detected 1 oom-kill event(s)",
                        }
                    ],
                },
            ),
            *_tool_use_entries(
                step=4,
                tool="read",
                args={"path": "/scratch/job.err"},
                result={
                    "path": "/scratch/job.err",
                    "start_line": 20,
                    "end_line": 40,
                    "total_lines": 400,
                },
            ),
        ]
    )

    assert memory.turns[0].reusable_results[-4:] == [
        "ssh_probe scheduler.detect_kind on cluster-a → slurm",
        "scheduler_query on cluster-a: job 12345 state=RUNNING queue=debug",
        "log_tail /scratch/job.err (12L, 1 errors: oom_killed)",
        "read /scratch/job.err L20-40/400",
    ]


def test_conversation_memory_updates_entities_for_tool_result_branch():
    memory = ConversationMemory.from_entries(
        [
            _request_entry("check scheduler through legacy logging"),
            *_tool_result_entries(
                step_index=1,
                tool="scheduler_query",
                args={"server": "cluster-b", "job_id": "7788"},
                result={"scheduler": "pbs", "job_id": "7788"},
            ),
        ]
    )

    assert memory.entities.last_server == "cluster-b"
    assert memory.entities.last_scheduler == "pbs"
    assert memory.entities.last_job_id == "7788"
