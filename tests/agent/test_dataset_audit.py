from __future__ import annotations

import importlib.util
import json
from pathlib import Path


def _load_audit_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "training"
        / "audit_dataset.py"
    )
    spec = importlib.util.spec_from_file_location("audit_dataset", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _episode(
    *,
    session="sess",
    turn=1,
    user="optimize h2o.xyz with gaussian charge 0 mult 1",
    command="chemsmart run gaussian -f h2o.xyz -c 0 -m 1 opt",
    gate="ok",
    execute_rc=None,
    denied=False,
    paused=False,
    final_answer="Done.",
    tool_events=None,
    extra_messages=None,
):
    events = (
        tool_events
        if tool_events is not None
        else [
            {
                "tool": "synthesize_command",
                "status": "ok",
                "args": {"request": user},
                "normalized_args": {"request": user},
                "result_summary": {
                    "command": command,
                    "semantic": {
                        "verdict": gate,
                        "generated_inputs": [
                            {"path": "/tmp/job.inp", "route": "# test"}
                        ],
                    },
                },
            }
        ]
    )
    messages = [{"role": "user", "content": user}]
    messages.extend(extra_messages or [])
    return {
        "v": 2,
        "session_id": session,
        "turn": turn,
        "provider": {"name": "openai", "model": "gpt-test"},
        "messages": messages,
        "synthesis": {
            "status": "ready",
            "command": command,
            "semantic_verdict": gate,
        },
        "outcome": {
            "gate": gate,
            "execute_rc": execute_rc,
            "denied": denied,
        },
        "tool_events": events,
        "invoked_tools": [e["tool"] for e in events],
        "paused": paused,
        "final_answer": final_answer,
    }


def _write(tmp_path, episodes):
    path = tmp_path / "episodes.jsonl"
    path.write_text(
        "".join(json.dumps(e) + "\n" for e in episodes), encoding="utf-8"
    )
    return path


def test_pipeline_counts_by_trajectory_and_kind(tmp_path):
    audit = _load_audit_module()
    episodes = [
        _episode(session="a", turn=1),
        _episode(
            session="a",
            turn=2,
            user="single point of x.xyz with gaussian charge 0 mult 1",
            command="chemsmart run gaussian -f x.xyz -c 0 -m 1 sp",
        ),
        _episode(
            session="b",
            turn=1,
            user="ts search of g.xyz with orca charge 0 mult 1",
            command="chemsmart run orca -f g.xyz -c 0 -m 1 ts",
        ),
    ]
    report = audit.audit([_write(tmp_path, episodes)])
    assert report["episodes_after_dedup"] == 3
    assert report["pipelines"]["by_trajectory"]["synthesize_command"] == 3
    assert report["pipelines"]["by_job_kind"] == {
        "gaussian.opt": 1,
        "gaussian.sp": 1,
        "orca.ts": 1,
    }


def test_quality_flags_for_bad_episodes(tmp_path):
    audit = _load_audit_module()
    episodes = [
        _episode(session="ok", turn=1),
        _episode(session="rej", turn=1, gate="reject"),
        _episode(session="rc", turn=1, execute_rc=1),
        _episode(
            session="cmd",
            turn=1,
            command="rm -rf /",
        ),
        _episode(session="deny", turn=1, denied=True),
        {
            "v": 2,
            "session_id": "nosynth",
            "turn": 1,
            "provider": {"name": "openai", "model": "gpt-test"},
            "messages": [{"role": "user", "content": "do the job"}],
            "synthesis": None,
            "outcome": {"gate": "none"},
            "tool_events": [],
            "invoked_tools": [],
            "paused": False,
            "final_answer": "",
        },
    ]
    report = audit.audit([_write(tmp_path, episodes)])
    counts = report["quality"]["flag_counts"]
    assert counts["gate_reject_unrepaired"] == 1
    assert counts["execute_failed"] == 1
    assert counts["malformed_command"] == 1
    assert counts["user_denied"] == 1
    assert counts["terminal_nosynth"] == 1
    assert report["quality"]["clean_episodes"] == 1


def test_gate_reject_with_positive_repair_is_not_flagged(tmp_path):
    audit = _load_audit_module()
    events = [
        {
            "tool": "synthesize_command",
            "status": "ok",
            "result_summary": {
                "command": "chemsmart run gaussian opt",
                "semantic": {"verdict": "reject"},
            },
        },
        {
            "tool": "repair_command",
            "status": "ok",
            "result_summary": {
                "command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt",
                "semantic": {"verdict": "ok"},
            },
        },
    ]
    episode = _episode(gate="reject", tool_events=events)
    report = audit.audit([_write(tmp_path, [episode])])
    assert "gate_reject_unrepaired" not in report["quality"]["flag_counts"]


def test_cross_turn_reject_is_retained_when_later_turn_repairs_it(tmp_path):
    audit = _load_audit_module()
    rejected = _episode(
        session="repair-chain",
        turn=1,
        command="chemsmart run gaussian -f bridge.xyz modred",
        gate="reject",
    )
    repaired = _episode(
        session="repair-chain",
        turn=2,
        user="Use charge 0 and multiplicity 1.",
        command=(
            "chemsmart run gaussian -f bridge.xyz -c 0 -m 1 "
            "modred -c '[[4,8]]' -j opt"
        ),
        gate="ok",
    )

    report = audit.audit([_write(tmp_path, [rejected, repaired])])

    assert report["quality"]["recovered_cross_turn_rejects"] == 1
    assert "gate_reject_unrepaired" not in report["quality"]["flag_counts"]
    assert report["quality"]["clean_episodes"] == 2


def test_over_frequent_skeletons_flagged_and_capped(tmp_path):
    audit = _load_audit_module()
    episodes = [
        _episode(
            session=f"s{i}",
            turn=1,
            user=f"optimize mol{i}.xyz with gaussian charge 0 mult 1",
        )
        for i in range(8)
    ] + [
        _episode(
            session="uniq",
            turn=1,
            user="wiberg bond index of complex.xyz with gaussian",
            command="chemsmart run gaussian -f complex.xyz -c 0 -m 1 wbi",
        )
    ]
    report = audit.audit(
        [_write(tmp_path, episodes)],
        max_per_skeleton=3,
        max_skeleton_share=0.5,
    )
    over = report["diversity"]["over_frequent_skeletons"]
    assert len(over) == 1
    assert list(over.values()) == [8]
    # 8 same-skeleton episodes capped to 3: 5 capped out.
    assert report["diversity"]["capped_out_episodes"] == 5
    assert report["clean_after_cap"] == 4  # 3 kept + 1 unique


def test_stability_flags(tmp_path):
    audit = _load_audit_module()
    naked_execute = _episode(
        session="exe",
        turn=1,
        tool_events=[
            {
                "tool": "execute_chemsmart_command",
                "status": "ok",
                "result_summary": {"returncode": 0},
            }
        ],
    )
    thrash_call = {
        "tool": "read_project_yaml",
        "status": "ok",
        "args": {"project": "demo"},
        "normalized_args": {"project": "demo"},
        "result_summary": {"ok": True},
    }
    thrash = _episode(
        session="loop",
        turn=1,
        tool_events=[dict(thrash_call) for _ in range(3)],
    )
    multi = _episode(
        session="multi",
        turn=1,
        extra_messages=[
            {"role": "assistant", "content": "Which basis?"},
            {"role": "user", "content": "def2svp"},
        ],
    )
    report = audit.audit([_write(tmp_path, [naked_execute, thrash, multi])])
    counts = report["stability"]["flag_counts"]
    assert counts["execute_without_synthesize"] == 1
    assert counts["repeated_identical_tool_call"] == 1
    assert report["stability"]["multi_turn_episodes"] == 1
    assert report["stability"]["unresolved_ask_user_pauses"] == 0


def test_session_level_multi_turn_metrics_reconstruct_separate_turns(tmp_path):
    audit = _load_audit_module()
    episodes = [
        _episode(session="chain", turn=1, user="optimize bridge.xyz"),
        _episode(
            session="chain",
            turn=2,
            user="Correction: constrain bond 4-8 instead.",
            command=(
                "chemsmart run gaussian -f bridge.xyz -c 0 -m 1 "
                "modred -c '[[4,8]]' -j opt"
            ),
        ),
        _episode(session="single", turn=1, user="single point of a.xyz"),
    ]

    report = audit.audit([_write(tmp_path, episodes)])

    stability = report["stability"]
    assert stability["session_count"] == 2
    assert stability["multi_turn_sessions"] == 1
    assert stability["multi_turn_session_share"] == 0.5
    assert stability["episodes_in_multi_turn_sessions"] == 2
    assert stability["record_level_multi_turn_episodes"] == 0
    assert report["quality"]["trainable_session_chains"] == 2


def test_job_kind_audit_canonicalizes_gaussian_td_to_tddft(tmp_path):
    audit = _load_audit_module()
    episode = _episode(
        session="td",
        user="Gaussian TD-DFT for dye.xyz, ten singlet states",
        command=(
            "chemsmart run gaussian -f dye.xyz -c 0 -m 1 "
            "td -n 10 -s singlets"
        ),
    )

    report = audit.audit([_write(tmp_path, [episode])])

    assert report["pipelines"]["by_job_kind"] == {"gaussian.tddft": 1}
    coverage = report["pipelines"]["canonical_kind_coverage"]
    assert "gaussian.tddft" in coverage["covered"]
    assert "gaussian.td" not in coverage["noncanonical"]


def test_command_coverage_counts_nested_qmmm_and_excludes_freq_cli(tmp_path):
    audit = _load_audit_module()
    episode = _episode(
        session="qmmm",
        user="Run a Gaussian QM/MM optimization for enzyme.xyz.",
        command=(
            "chemsmart run gaussian -f enzyme.xyz -c 0 -m 1 "
            "opt qmmm --high-level-atoms 1-8 --low-level-atoms 9-16"
        ),
    )

    report = audit.audit([_write(tmp_path, [episode])])

    assert report["pipelines"]["by_job_kind"] == {"gaussian.qmmm": 1}
    coverage = report["pipelines"]["canonical_kind_coverage"]
    assert "gaussian.qmmm" in coverage["covered"]
    assert "gaussian.freq" not in coverage["missing"]
    assert coverage["project_yaml_only"] == ["gaussian.freq", "orca.freq"]


def test_paused_superseded_episode_is_deduped_not_counted(tmp_path):
    audit = _load_audit_module()
    paused = _episode(session="p", turn=1, paused=True, final_answer="")
    resumed = _episode(
        session="p",
        turn=1,
        extra_messages=[
            {"role": "assistant", "content": "which file?"},
            {"role": "user", "content": "h2o.xyz"},
        ],
    )
    report = audit.audit([_write(tmp_path, [paused, resumed])])
    assert report["episodes_raw"] == 2
    assert report["episodes_after_dedup"] == 1
    assert report["stability"]["unresolved_ask_user_pauses"] == 0


def test_health_failures_and_strict_exit(tmp_path):
    audit = _load_audit_module()
    # 10 episodes, all the same skeleton -> ratio 0.1 < 0.35.
    episodes = [
        _episode(session=f"s{i}", turn=1, user="optimize 1.xyz")
        for i in range(10)
    ]
    path = _write(tmp_path, episodes)
    report = audit.audit([path])
    failures = audit.health_failures(report)
    assert any("skeleton_ratio" in failure for failure in failures)

    rc = audit.main([str(path), "--strict", "--out", str(tmp_path / "r.json")])
    assert rc == 1


def test_main_include_runs_audits_model_specific_stores(tmp_path):
    audit = _load_audit_module()
    training_dir = tmp_path / "agent-training"
    run_episode_dir = training_dir / "runs" / "teacher-a" / "episodes"
    run_episode_dir.mkdir(parents=True)
    episode_path = run_episode_dir / "202607.jsonl"
    episode_path.write_text(
        json.dumps(_episode(session="run-store")) + "\n", encoding="utf-8"
    )
    out = tmp_path / "audit.json"

    rc = audit.main(
        [
            "--training-dir",
            str(training_dir),
            "--include-runs",
            "--out",
            str(out),
        ]
    )

    assert rc == 0
    report = json.loads(out.read_text(encoding="utf-8"))
    assert report["episodes_after_dedup"] == 1
    assert report["pipelines"]["by_provider"] == {"openai:gpt-test": 1}


def test_write_clean_applies_quality_filter_and_cap(tmp_path):
    audit = _load_audit_module()
    episodes = [
        _episode(session="good", turn=1),
        _episode(session="bad", turn=1, execute_rc=2),
    ]
    path = _write(tmp_path, episodes)
    clean = tmp_path / "clean.jsonl"
    rc = audit.main([str(path), "--write-clean", str(clean)])
    assert rc == 0
    lines = [
        json.loads(line)
        for line in clean.read_text(encoding="utf-8").splitlines()
    ]
    assert len(lines) == 1
    assert lines[0]["session_id"] == "good"


def test_trajectory_fingerprint_deterministic_and_order_preserving():
    audit = _load_audit_module()
    ep_ab = {"tool_events": [{"tool": "a"}, {"tool": "b"}]}
    ep_ab2 = {"tool_events": [{"tool": "a"}, {"tool": "b"}]}
    ep_ba = {"tool_events": [{"tool": "b"}, {"tool": "a"}]}
    ep_empty = {"tool_events": []}
    ep_one = {"tool_events": [{"tool": "synthesize_command"}]}
    fp = audit.trajectory_fingerprint
    # deterministic + stable string hash
    assert fp(ep_ab) == fp(ep_ab2)
    assert isinstance(fp(ep_ab), str)
    # order preserved (a->b != b->a)
    assert fp(ep_ab) != fp(ep_ba)
    # empty and single-tool are DISTINCT (no collapse)
    assert fp(ep_empty) != fp(ep_one)
    assert fp(ep_empty) != fp(ep_ab)


def test_schema_efficiency_none_for_no_tools_and_error_ratio():
    audit = _load_audit_module()
    assert audit.schema_efficiency_score({"tool_events": []}) is None
    clean = {
        "tool_events": [
            {"tool": "x", "status": "ok"},
            {"tool": "y", "status": "ok"},
        ]
    }
    assert audit.schema_efficiency_score(clean) == 1.0
    mixed = {
        "tool_events": [
            {"tool": "x", "status": "ok"},
            {"tool": "y", "status": "error"},
        ]
    }
    assert audit.schema_efficiency_score(mixed) == 0.5


def test_sota_metrics_no_empty_bigram_inflation(tmp_path):
    audit = _load_audit_module()
    # 3 single-tool synthesize episodes + 1 no-tool: previously all 4 would
    # collide into the empty-bigram bucket and report 3 duplicates. Now
    # single-step is not counted as multi-step flooding => 0.
    episodes = [_episode(session=f"s{i}", turn=1) for i in range(3)] + [
        {
            "v": 2,
            "session_id": "chat",
            "turn": 1,
            "provider": {"name": "openai", "model": "m"},
            "messages": [{"role": "user", "content": "hi"}],
            "tool_events": [],
            "outcome": {"gate": "none"},
            "final_answer": "hello",
        }
    ]
    report = audit.audit([_write(tmp_path, episodes)])
    d = report["diversity"]
    assert d["semantic_trajectory_duplicates"] == 0
    # no-tool episode excluded from efficiency average
    assert report["stability"]["schema_efficiency_sampled_episodes"] == 3
    assert report["stability"]["average_schema_efficiency_score"] == 1.0
