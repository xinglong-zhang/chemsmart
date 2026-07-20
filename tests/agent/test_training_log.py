from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from chemsmart.agent.harness.terminal_state import (
    assertion,
    build_terminal_state,
)
from chemsmart.agent.training_log import (
    TrainingEpisodeWriter,
    TrainingLogConfig,
    load_training_log_config,
)


def _load_export_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "training"
        / "export_sft.py"
    )
    spec = importlib.util.spec_from_file_location("export_sft", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _read_jsonl(path: Path):
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def test_default_training_log_is_disabled_under_pytest():
    assert load_training_log_config().enabled is False


def test_training_episode_writer_masks_paths_and_stores_prompt(tmp_path):
    cwd = tmp_path / "work"
    cwd.mkdir()
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )

    written = writer.write_episode(
        session_id="sess-1",
        turn=1,
        provider_name="openai",
        model="gpt-test",
        cwd=str(cwd),
        messages=[
            {"role": "system", "content": "system prompt"},
            {
                "role": "user",
                "content": f"optimize {cwd}/examples/h2o.xyz",
            },
        ],
        tool_records=[
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result": {
                    "status": "ready",
                    "command": (
                        f"chemsmart run gaussian -f {cwd}/examples/h2o.xyz "
                        "-c 0 -m 1 opt"
                    ),
                    "schema_variant": "run/gaussian[opt,sp]",
                    "semantic": {
                        "verdict": "ok",
                        "generated_inputs": [
                            {
                                "path": "/tmp/chemsmart-gate/h2o_opt.com",
                                "route": "# opt b3lyp def2svp",
                            }
                        ],
                    },
                },
            }
        ],
        approvals_count=1,
        denials_count=0,
        available_tools=["synthesize_command", "repair_command"],
        final_answer=f"Prepared {cwd}/examples/h2o.xyz",
    )

    assert written is not None
    record = _read_jsonl(written)[0]
    assert record["system_prompt_sha"]
    assert (tmp_path / "training" / "prompts").is_dir()
    assert str(cwd) not in json.dumps(record)
    assert "./examples/h2o.xyz" in json.dumps(record)
    assert record["synthesis"]["schema_variant"] == "run/gaussian[opt,sp]"
    assert record["outcome"]["gate"] == "ok"
    assert record["v"] == 2
    assert record["available_tools"] == [
        "repair_command",
        "synthesize_command",
    ]
    assert record["invoked_tools"] == ["synthesize_command"]
    assert record["tools"] == ["synthesize_command"]
    assert record["tool_events"][0]["semantic_verdict"] == "ok"
    assert record["final_answer"] == "Prepared ./examples/h2o.xyz"


def test_training_episode_writer_drops_hidden_provider_message_fields(
    tmp_path,
):
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-public-messages",
        turn=1,
        provider_name="deepseek",
        model="deepseek-test",
        messages=[
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "prepare a job",
                "reasoning_content": "hidden provider trace",
                "annotations": ["transport-only"],
                "audio": {"id": "transport-only"},
            },
        ],
        tool_records=[],
        final_answer="Ready.",
    )

    assert written is not None
    record = _read_jsonl(written)[0]
    assert record["messages"] == [{"role": "user", "content": "prepare a job"}]
    assert "reasoning_content" not in json.dumps(record)
    assert "transport-only" not in json.dumps(record)


def test_training_episode_writer_masks_terminal_state_and_exports_receipt(
    tmp_path,
):
    export_sft = _load_export_module()
    cwd = tmp_path / "workspace"
    cwd.mkdir()
    terminal = build_terminal_state(
        action="submit_job",
        command="chemsmart sub -s mock-pbs gaussian ...",
        returncode=0,
        server="mock-pbs",
        scheduler="PBS",
        assertions=[
            assertion("server.yaml_exists", expected=True, observed=True),
            assertion("submit.returncode", expected=0, observed=0),
        ],
        artifacts=[{"kind": "server_yaml", "path": str(cwd / "server.yaml")}],
    )
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-terminal",
        turn=1,
        provider_name="deepseek",
        model="deepseek-v4-pro",
        cwd=str(cwd),
        messages=[
            {"role": "system", "content": "public system"},
            {"role": "user", "content": "Submit to mock-pbs."},
            {"role": "assistant", "content": "Submitted."},
        ],
        tool_records=[
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result": {
                    "status": "ready",
                    "command": "chemsmart sub -s mock-pbs gaussian ...",
                    "semantic": {
                        "verdict": "ok",
                        "generated_inputs": [
                            {"path": "/tmp/mock.com", "route": "# opt"}
                        ],
                    },
                },
            },
            {
                "tool": "execute_chemsmart_command",
                "status": "ok",
                "result": {"ok": True, "returncode": 0},
            },
        ],
        final_answer="Submitted.",
        terminal_state=terminal,
    )

    assert written is not None
    record = _read_jsonl(written)[0]
    assert record["terminal_state"]["all_passed"] is True
    assert str(cwd) not in json.dumps(record)

    counts = export_sft.export_sft(
        episode_paths=[written],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["tool_loop_written"] == 1
    assert counts["terminal_state_written"] == 1
    receipt = _read_jsonl(tmp_path / "out" / "terminal_state_assertions.jsonl")
    assert receipt[0]["terminal_state"]["status"] == "passed"


def test_export_routes_submission_without_terminal_state_to_review(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-submit-without-state",
        "turn": 1,
        "provider": {"name": "deepseek", "model": "deepseek-v4-pro"},
        "messages": [
            {"role": "user", "content": "Submit to the cluster."},
            {"role": "assistant", "content": "Submitted."},
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart sub -s cluster gaussian -f h2o.xyz -c 0 -m 1 opt",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/h2o.com", "route": "# opt"}
            ],
        },
        "outcome": {"gate": "ok", "execute_rc": 0},
        "final_answer": "Submitted.",
    }
    source = tmp_path / "episodes.jsonl"
    source.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[source],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["tool_loop_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["skipped"]["tool_loop:missing_terminal_state"] == 1


def test_export_requires_one_system_prompt_and_strips_hidden_fields(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-public-export",
        "turn": 1,
        "provider": {"name": "deepseek", "model": "deepseek-test"},
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "prepare h2o",
                "reasoning_content": "hidden trace",
                "annotations": ["transport-only"],
            },
            {
                "role": "assistant",
                "content": "Done.",
                "reasoning_content": "hidden trace",
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart run gaussian -f h2o.xyz -c 0 -m 1 opt",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/h2o.com", "route": "# opt"}
            ],
        },
        "outcome": {"gate": "ok", "execute_rc": None},
        "tool_events": [{"tool": "synthesize_command", "status": "ok"}],
        "final_answer": "Done.",
    }
    source = tmp_path / "episodes.jsonl"
    source.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[source],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["tool_loop_written"] == 1
    exported = _read_jsonl(tmp_path / "out" / "tool_loop_sft.jsonl")[0]
    assert sum(m["role"] == "system" for m in exported["messages"]) == 1
    assert all(
        "reasoning_content" not in message and "annotations" not in message
        for message in exported["messages"]
    )

    episode["session_id"] = "sess-missing-system"
    episode["messages"] = [
        {"role": "user", "content": "prepare h2o"},
        {"role": "assistant", "content": "Done."},
    ]
    source.write_text(json.dumps(episode) + "\n", encoding="utf-8")
    counts = export_sft.export_sft(
        episode_paths=[source],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out-missing-system",
    )
    assert counts["tool_loop_written"] == 0
    manifest = json.loads(
        (tmp_path / "out-missing-system" / "manifest.json").read_text()
    )
    assert manifest["skipped"]["tool_loop:missing_system_prompt"] == 1


def test_export_sft_splits_tool_loop_and_verified_compact_spec(tmp_path):
    export_sft = _load_export_module()
    training_dir = tmp_path / "training"
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=training_dir, mask_paths=False)
    )

    written = writer.write_episode(
        session_id="sess-2",
        turn=1,
        provider_name="openai",
        model="gpt-test",
        cwd=str(tmp_path),
        messages=[
            {"role": "system", "content": "tool loop prompt"},
            {
                "role": "user",
                "content": "Prepare Gaussian opt for examples/h2o.xyz.",
            },
        ],
        tool_records=[
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result": {
                    "status": "ready",
                    "command": (
                        "chemsmart run gaussian -f examples/h2o.xyz "
                        "-c 0 -m 1 opt"
                    ),
                    "schema_variant": "run/gaussian[opt,sp]",
                    "raw_response": json.dumps(
                        {
                            "intent": "workflow",
                            "jobs": [
                                {
                                    "id": 1,
                                    "kind": "gaussian.opt",
                                    "file": "examples/h2o.xyz",
                                    "charge": 0,
                                    "mult": 1,
                                }
                            ],
                        }
                    ),
                    "semantic": {
                        "verdict": "ok",
                        "generated_inputs": [
                            {
                                "path": "/tmp/chemsmart-gate/h2o_opt.com",
                                "route": "# opt b3lyp def2svp",
                            }
                        ],
                    },
                },
            }
        ],
        approvals_count=0,
        denials_count=0,
        available_tools=["synthesize_command"],
        final_answer="Prepared a Gaussian optimization command.",
    )
    assert written is not None

    counts = export_sft.export_sft(
        episode_paths=[written],
        training_dir=training_dir,
        out_dir=tmp_path / "out",
    )

    assert counts["episodes_seen"] == 1
    assert counts["tool_loop_written"] == 1
    assert counts["compact_spec_written"] == 1
    assert counts["compact_spec_skipped"] == 0
    assert counts["command_answer_written"] == 1
    tool_records = _read_jsonl(tmp_path / "out" / "tool_loop_sft.jsonl")
    compact_records = _read_jsonl(tmp_path / "out" / "compact_spec_sft.jsonl")
    command_records = _read_jsonl(
        tmp_path / "out" / "command_answer_sft.jsonl"
    )
    assert tool_records[0]["messages"][0]["role"] == "system"
    assert compact_records[0]["messages"][-1]["content"].startswith(
        '{"intent":"workflow"'
    )
    assert compact_records[0]["meta"]["commands"] == [
        "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 opt"
    ]
    assert command_records[0]["messages"][-1]["content"] == (
        "Prepared a Gaussian optimization command."
    )
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["written"]["tool_loop"] == 1
    assert manifest["diversity"]["distinct_query_skeletons"] == 1
    assert manifest["output_dir"] == "../out"
    assert manifest["files"]["tool_loop"] == "../out/tool_loop_sft.jsonl"


def test_export_sft_recovers_compact_spec_from_semantic_command(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-command",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "system_prompt_sha": "",
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "Prepare Gaussian opt for examples/h2o.xyz.",
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -f examples/h2o.xyz -c 0 -m 1 opt"
            ),
            "schema_variant": "run/gaussian[opt,sp]",
            "semantic_verdict": "ok",
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["compact_spec_written"] == 1
    compact_records = _read_jsonl(tmp_path / "out" / "compact_spec_sft.jsonl")
    spec = json.loads(compact_records[0]["messages"][-1]["content"])
    assert spec["jobs"][0]["kind"] == "gaussian.opt"
    assert spec["jobs"][0]["file"] == "examples/h2o.xyz"
    assert compact_records[0]["meta"]["source"] == (
        "deterministic_command_reverse"
    )


def test_export_sft_project_yaml_family(tmp_path):
    export_sft = _load_export_module()
    yaml_text = "gas:\n  functional: b3lyp\n  basis: def2svp\n  freq: true\n"
    episode = {
        "v": 2,
        "session_id": "sess-yaml",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "Build a project YAML for B3LYP/def2-SVP.",
            },
        ],
        "tool_events": [
            {
                "tool": "render_project_yaml",
                "status": "ok",
                "result_summary": {"yaml_text": yaml_text},
            },
            {
                "tool": "validate_project_yaml",
                "status": "ok",
                "result_summary": {"verdict": "ok"},
            },
        ],
        "outcome": {"gate": "none"},
        "final_answer": yaml_text,
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["project_yaml_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "project_yaml_sft.jsonl")
    assert records[0]["meta"]["yaml_text"] == yaml_text
    assert records[0]["meta"]["validation_verdict"] == "ok"


def test_export_sft_command_answer_uses_episode_workspace_label(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-command-answer",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "Prepare Gaussian opt for examples/h2o.xyz.",
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p test -f examples/h2o.xyz "
                "-c 0 -m 1 opt"
            ),
            "schema_variant": "run/gaussian[opt,sp]",
            "semantic_verdict": "ok",
        },
        "outcome": {"gate": "ok"},
        "cwd_masked": ".",
        "workspace": {
            "yaml_loaded": True,
            "project": "test",
            "program": "gaussian",
            "path": "./.chemsmart/gaussian/test.yaml",
        },
        "tool_events": [
            {
                "tool": "read_project_yaml",
                "status": "ok",
                "result_summary": {
                    "program": "gaussian",
                    "project_name": "test",
                    "yaml_text": (
                        "gas:\n  functional: b3lyp\n"
                        "  basis: def2svp\n  freq: true\n"
                        "solv:\n  functional: b3lyp\n"
                        "  basis: def2svp\n  freq: false\n"
                    ),
                },
            }
        ],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["command_answer_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "command_answer_sft.jsonl")
    deterministic = records[0]["meta"]["deterministic_answer"]
    assert "- workspace: `.`" in deterministic
    assert (
        "- resolved method: functional `b3lyp`, basis `def2svp`"
        in deterministic
    )
    assert "project settings could not be resolved" not in deterministic
    assert str(Path.cwd()) not in deterministic


def test_export_sft_masks_runtime_temp_paths_at_write(tmp_path):
    export_sft = _load_export_module()
    temp_path = (
        "/private/var/folders/fl/example/T/"
        "chemsmart-gate-abcd1234/h2o_opt_fake.com"
    )
    episode = {
        "v": 2,
        "session_id": "sess-temp-path",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "Prepare Gaussian opt for examples/h2o.xyz.",
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -f examples/h2o.xyz " "-c 0 -m 1 opt"
            ),
            "schema_variant": "run/gaussian[opt,sp]",
            "semantic_verdict": "ok",
        },
        "outcome": {"gate": "ok"},
        "tool_events": [
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result_summary": {
                    "command": (
                        "chemsmart run gaussian -f examples/h2o.xyz "
                        "-c 0 -m 1 opt"
                    ),
                    "semantic": {
                        "verdict": "ok",
                        "generated_inputs": [{"path": temp_path}],
                    },
                },
            }
        ],
        "final_answer": f"Generated input file: {temp_path}",
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["tool_loop_written"] == 1
    assert counts["command_answer_written"] == 1
    exported = (tmp_path / "out" / "tool_loop_sft.jsonl").read_text() + (
        tmp_path / "out" / "command_answer_sft.jsonl"
    ).read_text()
    assert "/private/var/folders" not in exported
    assert "<runtime-temp-path>" in exported


def test_export_sft_read_project_yaml_is_not_project_yaml_family(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-read-yaml",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {
                "role": "user",
                "content": "Prepare Gaussian opt with active project.",
            }
        ],
        "tool_events": [
            {
                "tool": "read_project_yaml",
                "status": "ok",
                "result_summary": {
                    "yaml_text": (
                        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
                    ),
                    "verdict": "ok",
                },
            }
        ],
        "outcome": {"gate": "none"},
        "final_answer": "The active YAML is valid.",
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["project_yaml_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["skipped"]["project_yaml:not_project_yaml_workflow"] == 1


def test_export_sft_repair_pair_family(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-repair",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [{"role": "user", "content": "Repair this command."}],
        "tool_events": [
            {
                "tool": "synthesize_command",
                "status": "ok",
                "result_summary": {
                    "command": "chemsmart run gaussian opt",
                    "semantic": {
                        "verdict": "reject",
                        "failed_rule_ids": [
                            "cmd.semantic.generated_input_missing"
                        ],
                    },
                },
            },
            {
                "tool": "repair_command",
                "status": "ok",
                "result_summary": {
                    "command": (
                        "chemsmart run gaussian -f examples/h2o.xyz "
                        "-c 0 -m 1 opt"
                    ),
                    "semantic": {"verdict": "ok", "failed_rule_ids": []},
                },
            },
        ],
        "outcome": {"gate": "reject"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["repair_pairs_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "repair_pairs.jsonl")
    assert records[0]["rejected"]["command"] == "chemsmart run gaussian opt"
    assert records[0]["chosen"]["semantic"]["verdict"] == "ok"


def test_export_sft_excludes_rejected_tool_loop_by_default(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 1,
        "session_id": "sess-3",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "system_prompt_sha": "",
        "messages": [{"role": "user", "content": "bad command"}],
        "synthesis": {
            "schema_variant": "run/gaussian[opt,sp]",
            "raw_response": "{}",
        },
        "outcome": {"gate": "reject"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["episodes_seen"] == 1
    assert counts["tool_loop_written"] == 0
    assert counts["compact_spec_written"] == 0


def test_export_sft_dedups_paused_then_resumed_turn(tmp_path):
    # An ask_user pause snapshots the turn; the resumed completion re-writes
    # the same (session_id, turn). Export must keep only the completion.
    export_sft = _load_export_module()

    def episode(paused: bool, messages: list[dict]) -> dict:
        return {
            "v": 2,
            "session_id": "sess-pause",
            "turn": 3,
            "provider": {"name": "openai", "model": "gpt-test"},
            "system_prompt_sha": "",
            "messages": [
                {"role": "system", "content": "public system"},
                *messages,
            ],
            "synthesis": None,
            "outcome": {"gate": "none"},
            "tool_events": [],
            "paused": paused,
            "final_answer": "" if paused else "Done.",
        }

    partial = episode(
        True,
        [{"role": "user", "content": "optimize my molecule"}],
    )
    complete = episode(
        False,
        [
            {"role": "user", "content": "optimize my molecule"},
            {"role": "assistant", "content": "Which file?"},
            {"role": "user", "content": "h2o.xyz"},
            {"role": "assistant", "content": "Done."},
        ],
    )
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        json.dumps(partial) + "\n" + json.dumps(complete) + "\n",
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["episodes_seen"] == 1  # deduped before stats
    records = _read_jsonl(tmp_path / "out" / "tool_loop_sft.jsonl")
    assert len(records) == 1
    assert len(records[0]["messages"]) == 5  # system + completed conversation


def test_export_sft_reconstructs_correction_chain_and_repair_pair(tmp_path):
    export_sft = _load_export_module()
    rejected_command = "chemsmart run gaussian -p demo -f bridge.xyz modred"
    repaired_command = (
        "chemsmart run gaussian -p demo -f bridge.xyz -c 0 -m 1 "
        "modred -c '[[4,8]]' -j opt"
    )
    first = {
        "v": 2,
        "session_id": "sess-correction",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher"},
        "messages": [
            {"role": "system", "content": "public system"},
            {"role": "user", "content": "Keep bond 4-8 fixed."},
            {"role": "assistant", "content": "I need charge and spin."},
        ],
        "synthesis": {
            "status": "needs_clarification",
            "command": rejected_command,
            "semantic_verdict": "reject",
            "failed_rule_ids": ["cmd.semantic.missing_charge"],
        },
        "outcome": {"gate": "reject"},
        "tool_events": [],
        "final_answer": "I need charge and spin.",
        "paused": False,
    }
    second = {
        "v": 2,
        "session_id": "sess-correction",
        "turn": 2,
        "provider": {"name": "openai", "model": "teacher"},
        "messages": [
            {"role": "system", "content": "public system"},
            {
                "role": "user",
                "content": "Use charge 0 and multiplicity 1; make it modred.",
            },
            {
                "role": "assistant",
                "content": "Prepared the corrected command.",
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": repaired_command,
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/bridge.com", "route": "# opt=modredundant"}
            ],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
        "final_answer": "Prepared the corrected command.",
        "paused": False,
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        json.dumps(first) + "\n" + json.dumps(second) + "\n",
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["episodes_seen"] == 2
    assert counts["tool_loop_written"] == 1
    assert counts["repair_pairs_written"] == 1
    tool_loop = _read_jsonl(tmp_path / "out" / "tool_loop_sft.jsonl")[0]
    users = [
        message["content"]
        for message in tool_loop["messages"]
        if message["role"] == "user"
    ]
    assert users == [
        "Keep bond 4-8 fixed.",
        "Use charge 0 and multiplicity 1; make it modred.",
    ]
    assert tool_loop["meta"]["source_episode_count"] == 2
    assert tool_loop["meta"]["user_turn_count"] == 2
    assert tool_loop["meta"]["recovered_reject"] is True
    repair = _read_jsonl(tmp_path / "out" / "repair_pairs.jsonl")[0]
    assert repair["rejected"]["command"] == rejected_command
    assert repair["chosen"]["command"] == repaired_command
    assert repair["meta"]["source"] == "cross_turn"


def test_export_sft_routes_terminal_nosynth_chain_to_review(tmp_path):
    export_sft = _load_export_module()
    episodes = [
        {
            "v": 2,
            "session_id": "sess-nosynth",
            "turn": turn,
            "provider": {"name": "openai", "model": "teacher"},
            "messages": [{"role": "user", "content": content}],
            "synthesis": None,
            "outcome": {"gate": "none"},
            "tool_events": [],
            "final_answer": "",
            "paused": False,
        }
        for turn, content in (
            (1, "Optimize then derive RESP charges."),
            (2, "Prepare only the RESP command."),
        )
    ]
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        "".join(json.dumps(episode) + "\n" for episode in episodes),
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["tool_loop_written"] == 0
    review = _read_jsonl(tmp_path / "out" / "tool_loop_review.jsonl")
    assert len(review) == 1
    assert review[0]["meta"]["skip_reason"] == "terminal_no_response"
    assert review[0]["meta"]["user_turn_count"] == 2


def test_export_sft_routes_terminal_tool_error_to_review(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-tool-error",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher"},
        "messages": [
            {"role": "user", "content": "Prepare an optimization."},
            {
                "role": "assistant",
                "content": "I will synthesize the command.",
                "tool_calls": [{"id": "call-1", "type": "function"}],
            },
            {
                "role": "tool",
                "tool_call_id": "call-1",
                "content": '{"ok": false, "error": "provider failed"}',
            },
        ],
        "synthesis": None,
        "outcome": {"gate": "none"},
        "tool_events": [
            {
                "tool": "synthesize_command",
                "status": "error",
                "result_summary": {"ok": False, "error": "provider failed"},
            }
        ],
        "final_answer": "I will synthesize the command.",
        "paused": False,
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["tool_loop_written"] == 0
    review = _read_jsonl(tmp_path / "out" / "tool_loop_review.jsonl")
    assert review[0]["meta"]["skip_reason"] == "terminal_tool_error"


def test_writer_records_paused_flag(tmp_path):
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-p",
        turn=1,
        provider_name="openai",
        model="gpt-test",
        messages=[{"role": "user", "content": "hello"}],
        tool_records=[],
        paused=True,
    )
    assert written is not None
    record = _read_jsonl(written)[0]
    assert record["paused"] is True


def test_run_loop_wires_paused_and_completed_episodes(tmp_path):
    # End-to-end wiring: an ask_user pause writes a paused episode, the
    # resumed completion writes the same (session_id, turn) unpaused.
    from chemsmart.agent.core import AgentSession
    from chemsmart.agent.loop import ASK_USER_TOOL_NAME

    from ._agent_session_helpers import FakeProvider
    from ._loop_helpers import (
        ScriptedRegistry,
        openai_final_response,
        openai_tool_call_response,
        tool_call,
    )

    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call(
                        "call_ask",
                        ASK_USER_TOOL_NAME,
                        {"question": "which file?", "options": []},
                    )
                )
            },
            {"__raw_response__": openai_final_response("Done with h2o.xyz.")},
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ScriptedRegistry({"recommend_method": {"method": "b3lyp"}}),
        session_root=tmp_path / "sessions",
    )
    session._training_writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )

    session.run_loop("Optimize my molecule.")
    session.run_loop("h2o.xyz")

    episode_files = sorted(
        (tmp_path / "training" / "episodes").glob("*.jsonl")
    )
    assert episode_files
    records = [
        json.loads(line)
        for path in episode_files
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert len(records) == 2
    assert records[0]["paused"] is True
    assert records[1]["paused"] is False
    assert records[0]["session_id"] == records[1]["session_id"]
    assert records[0]["turn"] == records[1]["turn"]
    assert records[1]["final_answer"] == "Done with h2o.xyz."
    assert records[1]["system_prompt_sha"]
    # The completed episode carries the full resumed conversation.
    blob = json.dumps(records[1]["messages"])
    assert "which file?" in blob and "h2o.xyz" in blob


def test_export_sft_filters_canonical_violations(tmp_path):
    # freq is runtime-owned: commands smuggling it via route params must not
    # reach positive SFT families.
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-canon",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "user", "content": "optimize and get frequencies"}
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt -r 'freq'"
            ),
            "semantic_verdict": "ok",
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["command_answer_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["skipped"]["command_answer:canonical_freq_in_route"] == 1


def test_export_sft_filters_short_form_runtime_owned_method_flags(tmp_path):
    # -x=--functional, -b=--basis are real CLI flags (they parse), so the
    # long-form-only guard let them poison positives. The parser-based guard
    # must reject short forms from every positive family.
    export_sft = _load_export_module()
    for command, reason in (
        (
            "chemsmart run gaussian -p demo -f a.xyz -c 0 -m 1 -x b3lyp opt",
            "runtime_owned_method_flag",
        ),
        (
            "chemsmart run gaussian -p demo -f a.xyz -c 0 -m 1 -b def2svp opt",
            "runtime_owned_basis_flag",
        ),
        (
            "chemsmart run orca -p demo -f a.xyz -c 0 -m 1 sp -si dmso",
            "runtime_owned_solvent_flag",
        ),
    ):
        episode = {
            "v": 2,
            "session_id": f"sess-leak-{reason}",
            "turn": 1,
            "provider": {"name": "openai", "model": "gpt-test"},
            "messages": [
                {"role": "user", "content": "prepare a job for a.xyz"}
            ],
            "synthesis": {
                "status": "ready",
                "command": command,
                "semantic_verdict": "ok",
                "reasoning": "public trace.",
                "reasoning_provenance": "public_decision_trace",
                "generated_input_evidence": [
                    {"path": "/tmp/a.com", "route": "# route"}
                ],
            },
            "outcome": {"gate": "ok"},
            "tool_events": [],
        }
        source = tmp_path / f"ep_{reason}.jsonl"
        source.write_text(json.dumps(episode) + "\n", encoding="utf-8")
        out = tmp_path / f"out_{reason}"
        counts = export_sft.export_sft(
            episode_paths=[source],
            training_dir=tmp_path / "training",
            out_dir=out,
        )
        # No positive family admits a runtime-owned leak.
        assert counts["command_answer_written"] == 0
        assert counts["reasoning_synthesis_written"] == 0
        assert counts["tool_loop_written"] == 0
        manifest = json.loads((out / "manifest.json").read_text())
        assert manifest["skipped"][f"command_answer:canonical_{reason}"] == 1
        # A clean command of the same shape still exports.
        clean = dict(episode)
        clean["session_id"] = f"sess-clean-{reason}"
        clean["synthesis"] = dict(episode["synthesis"])
        clean["synthesis"][
            "command"
        ] = "chemsmart run gaussian -p demo -f a.xyz -c 0 -m 1 opt"
        source.write_text(json.dumps(clean) + "\n", encoding="utf-8")
        out2 = tmp_path / f"out2_{reason}"
        counts2 = export_sft.export_sft(
            episode_paths=[source],
            training_dir=tmp_path / "training",
            out_dir=out2,
        )
        assert counts2["command_answer_written"] == 1


def test_export_sft_scan_numeric_x_flag_is_not_a_method_leak(tmp_path):
    # scan's `-x <number>` is dist_start (a coordinate), NOT --functional. The
    # parser-based guard must not false-positive on it.
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-scan-x",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [{"role": "user", "content": "scan the bond in a.xyz"}],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p demo -f a.xyz -c 0 -m 1 "
                "scan -c '[[1,2]]' -x 1.4 -n 10 -s 0.05"
            ),
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/a.com", "route": "# opt=modredundant"}
            ],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    source = tmp_path / "scan.jsonl"
    source.write_text(json.dumps(episode) + "\n", encoding="utf-8")
    counts = export_sft.export_sft(
        episode_paths=[source],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["command_answer_written"] == 1  # not falsely rejected


def test_export_sft_filters_long_route_canonical_violations(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-canon-long",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {
                "role": "user",
                "content": (
                    "Confirm a Gaussian transition state with one "
                    "imaginary frequency."
                ),
            }
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -f a.xyz -c 0 -m 1 "
                "--additional-route-parameters freq sp"
            ),
            "semantic_verdict": "ok",
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["command_answer_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["skipped"]["command_answer:canonical_freq_in_route"] == 1


def test_export_sft_filters_empty_qmmm_layer_override(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-empty-qmmm-layer",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher"},
        "messages": [
            {
                "role": "user",
                "content": "Prepare a Gaussian QMMM optimization.",
            }
        ],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt qmmm "
                '-ha 1-8 -la 9-16 -mx "" -mb ""'
            ),
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/a.com", "route": "oniom"}
            ],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["command_answer_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert (
        manifest["skipped"][
            "command_answer:canonical_empty_qmmm_layer_override"
        ]
        == 1
    )


def test_export_sft_repair_pair_from_repair_tool_only(tmp_path):
    # "fix this command" flows call repair_command directly; the result's
    # original_command supplies the rejected side of the contrast pair.
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-repair2",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "user", "content": "fix this broken command please"}
        ],
        "tool_events": [
            {
                "tool": "repair_command",
                "status": "ok",
                "result_summary": {
                    "command": (
                        "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt"
                    ),
                    "original_command": "chemsmart run gaussian opt -f a.xyz",
                    "repaired": True,
                    "semantic": {"verdict": "ok", "failed_rule_ids": []},
                },
            }
        ],
        "outcome": {"gate": "none"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["repair_pairs_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "repair_pairs.jsonl")
    assert records[0]["rejected"]["command"] == (
        "chemsmart run gaussian opt -f a.xyz"
    )
    assert records[0]["chosen"]["command"].endswith("-m 1 opt")


def test_export_sft_resolves_prompts_from_runs_subdirs(tmp_path):
    # Parallel harness episodes store prompt hashes under runs/<model>/.
    export_sft = _load_export_module()
    training_dir = tmp_path / "training"
    prompt_dir = training_dir / "runs" / "model-a" / "prompts"
    prompt_dir.mkdir(parents=True)
    (prompt_dir / "abc123def456.txt").write_text(
        "canonical system prompt", encoding="utf-8"
    )
    episode = {
        "v": 2,
        "session_id": "sess-prompt",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "system_prompt_sha": "abc123def456",
        "messages": [{"role": "user", "content": "optimize a.xyz (0/1)"}],
        "synthesis": None,
        "outcome": {"gate": "none"},
        "tool_events": [],
        "final_answer": "Ready to help.",
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=training_dir,
        out_dir=tmp_path / "out",
    )
    records = _read_jsonl(tmp_path / "out" / "tool_loop_sft.jsonl")
    assert records[0]["messages"][0] == {
        "role": "system",
        "content": "canonical system prompt",
    }


def test_synthesis_summary_preserves_public_reasoning_only():
    from chemsmart.agent.training_log import _synthesis_summary

    records = [
        {
            "tool": "synthesize_command",
            "status": "ok",
            "result": {
                "status": "ready",
                "command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt",
                "reasoning": "opt of a.xyz, gaussian, neutral singlet.",
                "reasoning_provenance": "public_decision_trace",
                "raw_response": "<think>private provider reasoning</think>",
                "semantic": {"verdict": "ok"},
            },
        }
    ]
    summary = _synthesis_summary(records)
    assert summary["reasoning"] == "opt of a.xyz, gaussian, neutral singlet."
    assert summary["reasoning_provenance"] == "public_decision_trace"
    assert "raw_response" not in summary


def test_synthesis_summary_keeps_only_compact_spec_raw_response():
    from chemsmart.agent.training_log import _synthesis_summary

    spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.qmmm",
                "file": "enzyme.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {"parent_job": "opt", "high_level_atoms": [1, 2]},
            }
        ],
    }
    summary = _synthesis_summary(
        [
            {
                "tool": "synthesize_command",
                "result": {
                    "status": "ready",
                    "command": "chemsmart run gaussian ...",
                    "raw_response": json.dumps(spec),
                },
            }
        ]
    )

    assert summary["raw_response"] == spec


def test_export_reasoning_synthesis_family(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-reason",
        "turn": 1,
        "provider": {"name": "qwen", "model": "qwen3-235b-thinking"},
        "messages": [
            {"role": "user", "content": "optimize h2o.xyz with gaussian, 0/1"},
            {
                "role": "assistant",
                "content": "",
                "reasoning_content": (
                    "Geometry optimization requested for h2o.xyz. Gaussian "
                    "engine, neutral singlet, opt subcommand."
                ),
                "tool_calls": [
                    {
                        "id": "c1",
                        "type": "function",
                        "function": {
                            "name": "synthesize_command",
                            "arguments": '{"request": "optimize h2o.xyz"}',
                        },
                    }
                ],
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart run gaussian -p demo -f h2o.xyz -c 0 -m 1 opt",
            "explanation": "Prepare a Gaussian optimization.",
            "reasoning": (
                "Geometry optimization requested for h2o.xyz. Gaussian "
                "engine, neutral singlet, opt subcommand."
            ),
            "reasoning_provenance": "public_decision_trace",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/chemsmart-gate/h2o_opt.com", "route": "# opt"}
            ],
            "missing_info": [],
            "alternatives": [],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["reasoning_synthesis_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "reasoning_synthesis_sft.jsonl")
    assistant = records[0]["messages"][-1]["content"]
    assert assistant.startswith("<think>")
    assert "</think>" in assistant
    # The trace was marked public at synthesis time, not recovered from a
    # provider-private assistant payload.
    assert "neutral singlet" in assistant
    # command JSON follows the think block
    json_part = assistant.split("</think>", 1)[1].strip()
    payload = json.loads(json_part)
    assert payload["command"].endswith("opt")
    assert payload["status"] == "ready"
    assert records[0]["messages"][0]["role"] == "system"
    assert (
        records[0]["meta"]["reasoning_provenance"] == "public_decision_trace"
    )


def test_export_reasoning_synthesis_routes_message_reasoning_to_review(
    tmp_path,
):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-review",
        "turn": 1,
        "provider": {"name": "qwen", "model": "qwen3-235b-thinking"},
        "messages": [
            {"role": "user", "content": "optimize h2o.xyz with gaussian, 0/1"},
            {
                "role": "assistant",
                "content": "",
                "reasoning_content": (
                    "This is post-tool presentation reasoning, not the "
                    "synthesis subcall reasoning."
                ),
            },
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart run gaussian -p demo -f h2o.xyz -c 0 -m 1 opt",
            "explanation": "Prepare a Gaussian optimization.",
            "reasoning": "",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/chemsmart-gate/h2o_opt.com", "route": "# opt"}
            ],
            "missing_info": [],
            "alternatives": [],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["reasoning_synthesis_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert (
        manifest["skipped"][
            "reasoning_synthesis:untrusted_reasoning:assistant_message"
        ]
        == 1
    )
    review = _read_jsonl(tmp_path / "out" / "reasoning_synthesis_review.jsonl")
    assert review[0]["meta"]["reasoning_source"] == "assistant_message"
    assert review[0]["meta"]["skip_reason"].startswith("untrusted_reasoning")
    assert "reasoning" not in review[0]


def test_export_reasoning_synthesis_skips_without_reasoning(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-noreason",
        "turn": 1,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": [
            {"role": "user", "content": "sp of a.xyz gaussian 0/1"},
            {"role": "assistant", "content": "", "tool_calls": []},
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart run gaussian -p demo -f a.xyz -c 0 -m 1 sp",
            "reasoning": "",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/chemsmart-gate/a_sp.com", "route": "# sp"}
            ],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")
    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    assert counts["reasoning_synthesis_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert manifest["skipped"]["reasoning_synthesis:no_reasoning_trace"] == 1


def test_export_skips_maxstep_request_mapped_to_maxcycle(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-maxstep-mismatch",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher"},
        "messages": [
            {
                "role": "user",
                "content": (
                    "gussian TS opt for proton_style.xyz, max step 8 plz, "
                    "neutral singlett. use loaded project."
                ),
            }
        ],
        "invoked_tools": ["synthesize_command"],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p demo -f proton_style.xyz "
                "-c 0 -m 1 -r maxcycle=8 ts"
            ),
            "reasoning": "The user asks for a TS command.",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {
                    "path": "/tmp/chemsmart-gate/proton_style.com",
                    "route": "# opt=(ts,calcfc,noeigentest)",
                }
            ],
        },
        "outcome": {"gate": "ok"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["reasoning_synthesis_written"] == 0
    assert counts["compact_spec_written"] == 0
    assert counts["command_answer_written"] == 0
    manifest = json.loads((tmp_path / "out" / "manifest.json").read_text())
    assert (
        manifest["skipped"][
            "reasoning_synthesis:semantic_maxstep_mapped_to_maxcycle"
        ]
        == 1
    )


def test_export_wrong_route_contrast_family(tmp_path):
    export_sft = _load_export_module()
    wrong = {
        "v": 2,
        "session_id": "sess-wrong-route",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher-a"},
        "messages": [
            {
                "role": "user",
                "content": (
                    "Using the loaded Gaussian project, create the command "
                    "for Wiberg bond indices on cationic_cmd.xyz, charge +1 "
                    "multiplicity 1."
                ),
            }
        ],
        "invoked_tools": [
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
        ],
        "tool_events": [
            {
                "tool": "extract_project_protocol",
                "status": "ok",
                "result_summary": {"ok": True},
            }
        ],
        "outcome": {"gate": "none"},
    }
    corrected = {
        "v": 2,
        "session_id": "sess-correct-route",
        "turn": 1,
        "provider": {"name": "openai", "model": "teacher-b"},
        "messages": [
            {
                "role": "user",
                "content": (
                    "Using the loaded Gaussian project, create the command "
                    "for Wiberg bond indices on cationic_cmd.xyz, charge +1 "
                    "multiplicity 1."
                ),
            }
        ],
        "invoked_tools": ["synthesize_command"],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p demo -f cationic_cmd.xyz "
                "-c 1 -m 1 wbi"
            ),
            "reasoning": "This is a command request for Gaussian WBI.",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {
                    "path": "/tmp/chemsmart-gate/wbi.com",
                    "route": "# pop=nboread",
                }
            ],
        },
        "outcome": {"gate": "ok"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        json.dumps(wrong) + "\n" + json.dumps(corrected) + "\n",
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["wrong_route_contrast_written"] == 1
    records = _read_jsonl(tmp_path / "out" / "wrong_route_contrast.jsonl")
    assert records[0]["rejected"]["route"] == "project_yaml"
    assert records[0]["chosen"]["route"] == "synthesize_command"
    assert records[0]["chosen"]["command"].endswith("wbi")
    assert records[0]["meta"]["contrast_key"] == "gaussian:wbi"


def test_export_wrong_route_contrast_accepts_oniom_yaml_then_cli_request(
    tmp_path,
):
    export_sft = _load_export_module()
    prompt = (
        "Create a project YAML for a two-layer Gaussian ONIOM optimization "
        "and then prepare the actual CLI command for enzyme.xyz."
    )
    wrong = {
        "v": 2,
        "session_id": "sess-oniom-wrong-route",
        "turn": 1,
        "messages": [{"role": "user", "content": prompt}],
        "invoked_tools": ["extract_project_protocol", "render_project_yaml"],
        "outcome": {"gate": "none"},
    }
    corrected = {
        "v": 2,
        "session_id": "sess-oniom-correct-route",
        "turn": 1,
        "messages": [{"role": "user", "content": prompt}],
        "invoked_tools": ["synthesize_command"],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run gaussian -p demo -f enzyme.xyz "
                "-c 0 -m 1 opt qmmm -ct 0 -mt 1 -ch 0 -mh 1 -ha 1-8 -la 9-16"
            ),
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/chemsmart-gate/enzyme.com", "route": "oniom"}
            ],
        },
        "outcome": {"gate": "ok"},
    }
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        json.dumps(wrong) + "\n" + json.dumps(corrected) + "\n",
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["wrong_route_contrast_written"] == 1
    record = _read_jsonl(tmp_path / "out" / "wrong_route_contrast.jsonl")[0]
    assert record["meta"]["contrast_key"] == "gaussian:qmmm"


def test_write_episode_records_sanitized_dataset_provenance(tmp_path):
    fake_secret = "sk-" + "abcdef0123456789abcdefghij"
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-prov",
        turn=1,
        provider_name="deepseek",
        model="deepseek-reasoner",
        messages=[{"role": "user", "content": "prepare a job"}],
        tool_records=[],
        final_answer="Ready.",
        dataset_provenance={
            "scenario_id": "orca_modred_multiturn_03",
            "scenario_family": "orca.modred.multiturn",
            "batch_id": "rk-batch-8",
            "teacher_id": "deepseek-reasoner",
            "fixture_id": "abc123def456",
            "schema_version": 8,
            "collector_version": "agentic_workflow_accum@1",
            "unknown_field": "dropped",
            "leaked_key": fake_secret,
        },
    )
    assert written is not None
    record = _read_jsonl(written)[0]
    provenance = record["dataset_provenance"]
    assert provenance["scenario_family"] == "orca.modred.multiturn"
    assert provenance["schema_version"] == "8"  # coerced to string
    assert "unknown_field" not in provenance  # whitelist only
    assert "leaked_key" not in provenance  # not a whitelisted key
    assert fake_secret[:21] not in json.dumps(record)  # secret masked


def test_write_episode_reads_dataset_provenance_from_env(
    tmp_path, monkeypatch
):
    from chemsmart.agent.training_log import DATASET_PROVENANCE_ENV

    monkeypatch.setenv(
        DATASET_PROVENANCE_ENV,
        json.dumps({"scenario_family": "gaussian.crest.multiturn"}),
    )
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-env-prov",
        turn=1,
        provider_name="deepseek",
        model="deepseek-reasoner",
        messages=[{"role": "user", "content": "prepare a job"}],
        tool_records=[],
        final_answer="Ready.",
    )
    assert written is not None
    record = _read_jsonl(written)[0]
    assert (
        record["dataset_provenance"]["scenario_family"]
        == "gaussian.crest.multiturn"
    )


def test_write_episode_omits_dataset_provenance_when_absent(tmp_path):
    # Legacy behaviour: no provenance argument and no env var leaves the field
    # off entirely so old readers and old rows stay identical.
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=tmp_path / "training")
    )
    written = writer.write_episode(
        session_id="sess-noprov",
        turn=1,
        provider_name="openai",
        model="gpt-test",
        messages=[{"role": "user", "content": "prepare a job"}],
        tool_records=[],
        final_answer="Ready.",
    )
    assert written is not None
    record = _read_jsonl(written)[0]
    assert "dataset_provenance" not in record


def test_export_meta_propagates_dataset_provenance(tmp_path):
    export_sft = _load_export_module()
    episode = {
        "v": 2,
        "session_id": "sess-prov-export",
        "turn": 1,
        "provider": {"name": "deepseek", "model": "deepseek-reasoner"},
        "dataset_provenance": {"scenario_family": "orca.qrc.multiturn"},
        "messages": [
            {
                "role": "user",
                "content": "prepare an ORCA QRC job for a.xyz 0/1",
            }
        ],
        "synthesis": {
            "status": "ready",
            "command": "chemsmart run orca -p demo -f a.xyz -c 0 -m 1 qrc",
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/a.inp", "route": "! qrc"}
            ],
        },
        "outcome": {"gate": "ok"},
        "tool_events": [],
    }
    source = tmp_path / "episodes.jsonl"
    source.write_text(json.dumps(episode) + "\n", encoding="utf-8")

    export_sft.export_sft(
        episode_paths=[source],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )
    command_records = _read_jsonl(
        tmp_path / "out" / "command_answer_sft.jsonl"
    )
    assert command_records[0]["meta"]["dataset_provenance"] == {
        "scenario_family": "orca.qrc.multiturn"
    }


def test_export_wrong_route_prefers_same_session_correction(tmp_path):
    export_sft = _load_export_module()
    prompt = (
        "Create an ORCA project YAML for QM/MM and then prepare the actual "
        "CLI command for route_gold.xyz, high-level atoms 1-12."
    )
    wrong = {
        "v": 2,
        "session_id": "same-session-route",
        "turn": 1,
        "messages": [{"role": "user", "content": prompt}],
        "invoked_tools": ["extract_project_protocol", "render_project_yaml"],
        "outcome": {"gate": "none"},
    }
    same_session = {
        "v": 2,
        "session_id": "same-session-route",
        "turn": 2,
        "messages": [{"role": "user", "content": "Correction: use CLI."}],
        "invoked_tools": ["synthesize_command"],
        "synthesis": {
            "status": "ready",
            "command": (
                "chemsmart run orca -p demo -f route_gold.xyz -c 0 -m 1 "
                "opt qmmm -j QMMM -ha 1-12 -ct 0 -mt 1 -lm AMBER=HardFirst"
            ),
            "semantic_verdict": "ok",
            "generated_input_evidence": [
                {"path": "/tmp/route.inp", "route": "! Opt QMMM"}
            ],
        },
        "outcome": {"gate": "ok"},
    }
    other = json.loads(json.dumps(same_session))
    other["session_id"] = "other-session-route"
    other["synthesis"]["command"] = (
        "chemsmart run orca -p demo -f other.xyz -c 0 -m 1 "
        "opt qmmm -j QMMM -ha 1-5 -ct 0 -mt 1 -lm UFF"
    )
    episode_path = tmp_path / "episodes.jsonl"
    episode_path.write_text(
        "\n".join(json.dumps(row) for row in (wrong, other, same_session))
        + "\n",
        encoding="utf-8",
    )

    counts = export_sft.export_sft(
        episode_paths=[episode_path],
        training_dir=tmp_path / "training",
        out_dir=tmp_path / "out",
    )

    assert counts["wrong_route_contrast_written"] == 1
    record = _read_jsonl(tmp_path / "out" / "wrong_route_contrast.jsonl")[0]
    assert record["chosen"]["command"].endswith("-lm AMBER=HardFirst")
    assert record["meta"]["wrong_session_id"] == "same-session-route"
