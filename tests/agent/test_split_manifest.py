"""Deterministic tests for the leakage-safe train/eval split builder.

Each test targets one observed defect in the previous first-turn-only,
force-eval split construction. No network, no provider, no real training store.
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts.training.build_split_manifest import build_split_manifest


def _write_export(tmp_path: Path, families: dict[str, list[dict]]) -> Path:
    export = tmp_path / "export"
    export.mkdir(parents=True)
    for family, rows in families.items():
        (export / f"{family}.jsonl").write_text(
            "".join(json.dumps(row) + "\n" for row in rows),
            encoding="utf-8",
        )
    return export


def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _tool_loop(
    session_id: str,
    *turns: str,
    recovered: bool = False,
    terminal: bool = False,
    provenance: dict | None = None,
) -> dict:
    messages: list[dict] = []
    for index, turn in enumerate(turns):
        messages.append({"role": "user", "content": turn})
        messages.append({"role": "assistant", "content": f"ok {index}"})
    meta: dict = {
        "family": "tool_loop",
        "session_id": session_id,
        "recovered_reject": recovered,
        "terminal_state": {"status": "passed"} if terminal else None,
    }
    if provenance is not None:
        meta["dataset_provenance"] = provenance
    return {"messages": messages, "meta": meta}


def _repair_pair(session_id: str, prompt: str) -> dict:
    return {
        "prompt": prompt,
        "rejected": {"command": "chemsmart run gaussian opt"},
        "chosen": {"command": "chemsmart run gaussian -f a.xyz -c 0 -m 1 opt"},
        "meta": {"family": "repair_pair", "session_id": session_id},
    }


def _wrong_route(session_id: str, prompt: str) -> dict:
    return {
        "prompt": prompt,
        "rejected": {"route": "project_yaml"},
        "chosen": {
            "route": "synthesize_command",
            "command": "chemsmart run gaussian -f b.xyz -c 0 -m 1 wbi",
        },
        "meta": {
            "family": "wrong_route_contrast",
            "session_id": session_id,
            "wrong_session_id": session_id,
        },
    }


def _terminal_assertion(session_id: str) -> dict:
    return {
        "terminal_state": {
            "status": "passed",
            "command": "chemsmart sub -s mock gaussian opt",
        },
        "meta": {"family": "terminal_state", "session_id": session_id},
    }


def _splits(export: Path) -> Path:
    return export / "splits"


def test_zero_all_turn_skeleton_leakage(tmp_path):
    # Two distinct sessions whose FIRST turns differ but whose LATER turn is an
    # identical query shape. First-turn-only splitting let the shared later
    # skeleton straddle train/eval; all-turn union must co-locate them.
    export = _write_export(
        tmp_path,
        {
            "tool_loop_sft": [
                _tool_loop("s1", "optimize alpha.xyz", "freeze bond 3 and 5"),
                _tool_loop("s2", "single point on beta.xyz", "freeze bond 3 and 5"),
            ]
        },
    )
    manifest = build_split_manifest(export, eval_fraction=50)
    assert manifest["leakage"]["all_turn_skeleton_leak_count"] == 0
    assert manifest["production_gate"]["no_all_turn_skeleton_cross_split"]
    rows_train = _read(_splits(export) / "train" / "tool_loop_sft.jsonl")
    rows_eval = _read(_splits(export) / "eval" / "tool_loop_sft.jsonl")
    # Both rows share a component and therefore a split.
    assert (len(rows_train), len(rows_eval)) in {(2, 0), (0, 2)}


def test_zero_scenario_family_leakage(tmp_path):
    prov = {"scenario_family": "gaussian.modred.multiturn", "scenario_id": "x"}
    export = _write_export(
        tmp_path,
        {
            "tool_loop_sft": [
                _tool_loop("s1", "totally distinct request one", provenance=prov),
                _tool_loop("s2", "an unrelated phrasing two", provenance=prov),
            ]
        },
    )
    manifest = build_split_manifest(export, eval_fraction=50)
    assert manifest["leakage"]["scenario_family_leak_count"] == 0
    assert manifest["production_gate"]["no_scenario_family_cross_split"]
    assert manifest["provenance_coverage"]["scenario_family_ratio"] == 1.0
    rows_train = _read(_splits(export) / "train" / "tool_loop_sft.jsonl")
    rows_eval = _read(_splits(export) / "eval" / "tool_loop_sft.jsonl")
    assert (len(rows_train), len(rows_eval)) in {(2, 0), (0, 2)}


def test_repaired_positive_chain_can_land_in_train(tmp_path):
    # A recovered-reject positive chain plus its repair pair share a session.
    # The old builder force-routed that whole component to eval; now it follows
    # the deterministic split and can be train.
    export = _write_export(
        tmp_path,
        {
            "tool_loop_sft": [
                _tool_loop("s1", "keep bond 4-8 fixed then opt", recovered=True)
            ],
            "repair_pairs": [_repair_pair("s1", "keep bond 4-8 fixed then opt")],
        },
    )
    manifest = build_split_manifest(export, eval_fraction=0)  # all -> train
    assert manifest["recovered_rejects_by_split"]["train"] == 1
    assert manifest["recovered_rejects_by_split"]["eval"] == 0
    assert _read(_splits(export) / "train" / "tool_loop_sft.jsonl")
    # Preference data stays out of the positive stream but under train.
    assert _read(
        _splits(export) / "train" / "preference" / "repair_pairs.jsonl"
    )
    assert not (_splits(export) / "train" / "repair_pairs.jsonl").exists()


def test_terminal_success_chain_can_land_in_train(tmp_path):
    export = _write_export(
        tmp_path,
        {
            "tool_loop_sft": [
                _tool_loop("s1", "submit an opt to the cluster", terminal=True)
            ],
            "terminal_state_assertions": [_terminal_assertion("s1")],
        },
    )
    manifest = build_split_manifest(export, eval_fraction=0)  # all -> train
    terminal = manifest["terminal_state_by_split"]
    assert terminal["tool_loop_terminal_success_chains"]["train"] == 1
    assert terminal["terminal_state_assertions"]["train"] == 1
    assert _read(
        _splits(export) / "train" / "evidence" / "terminal_state_assertions.jsonl"
    )


def test_contrast_data_is_not_force_routed_to_eval(tmp_path):
    # Contrast/evidence no longer force their component to eval, so an
    # unrelated positive component is not dragged along either.
    export = _write_export(
        tmp_path,
        {
            "tool_loop_sft": [
                _tool_loop("pos", "an unrelated positive optimization request")
            ],
            "wrong_route_contrast": [
                _wrong_route("wr", "a separate wiberg bond index request")
            ],
        },
    )
    manifest = build_split_manifest(export, eval_fraction=0)  # all -> train
    assert _read(_splits(export) / "train" / "tool_loop_sft.jsonl")
    assert _read(
        _splits(export) / "train" / "preference" / "wrong_route_contrast.jsonl"
    )
    # Nothing forced onto the eval side.
    assert manifest["counts"]["eval"] == {}


def test_legacy_rows_without_provenance_split_deterministically(tmp_path):
    families = {
        "tool_loop_sft": [
            _tool_loop(f"s{i}", f"distinct legacy request number {i}")
            for i in range(8)
        ]
    }
    export_a = _write_export(tmp_path / "a", families)
    export_b = _write_export(tmp_path / "b", families)
    manifest_a = build_split_manifest(export_a, eval_fraction=40)
    manifest_b = build_split_manifest(export_b, eval_fraction=40)
    assert manifest_a["counts"] == manifest_b["counts"]
    assert manifest_a["provenance_coverage"]["dataset_provenance_ratio"] == 0.0
    # Re-running the same export is byte-stable per split file.
    train_a = (_splits(export_a) / "train" / "tool_loop_sft.jsonl").read_text()
    train_b = (_splits(export_b) / "train" / "tool_loop_sft.jsonl").read_text()
    assert train_a == train_b
    manifest_a2 = build_split_manifest(export_a, eval_fraction=40)
    assert manifest_a2["counts"] == manifest_a["counts"]
    train_a2 = (_splits(export_a) / "train" / "tool_loop_sft.jsonl").read_text()
    assert train_a2 == train_a


def test_positive_export_records_are_preserved_verbatim(tmp_path):
    row = _tool_loop("s1", "prepare a gaussian optimization for water")
    export = _write_export(tmp_path, {"tool_loop_sft": [row]})
    build_split_manifest(export, eval_fraction=0)
    written = _read(_splits(export) / "train" / "tool_loop_sft.jsonl")
    assert written == [row]  # content round-trips unchanged
    # Diagnostics required by the manifest contract are present.
    manifest = build_split_manifest(export, eval_fraction=0)
    assert "per_family_eval_shares" in manifest
    assert "command_template_overlap" in manifest
    assert "provenance_coverage" in manifest
