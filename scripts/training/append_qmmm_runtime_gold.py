"""Append runtime-verified QMMM gold episodes to the local training ledger.

The records are deliberately labelled as ``chemsmart-runtime`` rather than as
frontier-model output.  They provide exact CLI/harness supervision for cases
where a teacher response was semantically incomplete.  Every positive command
is checked by the production semantic gate before it is appended.
"""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path
from typing import Any

from chemsmart.agent.harness.command_semantics import (
    evaluate_command_semantics,
)
from chemsmart.agent.training_log import (
    TrainingEpisodeWriter,
    TrainingLogConfig,
)

REPO = Path(__file__).resolve().parents[2]


def _semantic_payload(result: Any) -> dict[str, Any]:
    payload = result.to_dict()
    return {
        "verdict": payload.get("verdict"),
        "failed_rule_ids": payload.get("failed_rule_ids") or [],
        "generated_inputs": [
            {
                "path": "<runtime-temp-path>",
                "route": item.get("route"),
            }
            for item in payload.get("generated_inputs") or []
            if isinstance(item, dict)
        ],
        "issues": [
            {
                "rule_id": issue.get("rule_id"),
                "severity": issue.get("severity"),
                "message": issue.get("message"),
            }
            for issue in payload.get("issues") or []
            if isinstance(issue, dict)
        ],
    }


def _synthesis_record(
    command: str,
    semantic: dict[str, Any],
    *,
    status: str,
    explanation: str,
    compact_spec: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "tool": "synthesize_command",
        "status": "ok",
        "args": {"request": "runtime-verified QMMM fixture"},
        "normalized_args": {"request": "runtime-verified QMMM fixture"},
        "result": {
            "status": status,
            "action": "synthesize_command",
            "command": command,
            "schema_variant": "runtime-gold/qmmm",
            "explanation": explanation,
            "reasoning": (
                "Public decision trace: select nested qmmm parent; preserve "
                "1-based regions and explicit layer state; verify generated input."
            ),
            "reasoning_provenance": "public_decision_trace",
            "semantic": semantic,
            **(
                {
                    "raw_response": json.dumps(
                        compact_spec,
                        ensure_ascii=False,
                        separators=(",", ":"),
                    )
                }
                if compact_spec is not None
                else {}
            ),
        },
    }


def _workspace(program: str) -> dict[str, Any]:
    return {
        "project": "demo",
        "program": program,
        "path": f"./.chemsmart/{program}/demo.yaml",
        "candidates": [f"./.chemsmart/{program}/demo.yaml"],
        "yaml_loaded": True,
        "message": f"Loaded workspace project YAML: ./.chemsmart/{program}/demo.yaml",
    }


def _assistant(command: str, semantic: dict[str, Any]) -> str:
    route = ", ".join(
        str(item.get("route") or "").strip()
        for item in semantic.get("generated_inputs") or []
        if item.get("route")
    )
    return (
        "Validated ChemSmart CLI command:\n\n"
        f"{command}\n\n"
        f"Semantic verdict: {semantic.get('verdict')}.\n"
        f"Generated route evidence: {route or '<none>'}."
    )


def _append(
    writer: TrainingEpisodeWriter,
    *,
    session_id: str,
    turn: int,
    program: str,
    user: str,
    assistant: str,
    tool_records: list[dict[str, Any]],
) -> None:
    episode_file = writer.config.dir / "episodes" / "202607.jsonl"
    if episode_file.exists():
        for line in episode_file.read_text(encoding="utf-8").splitlines():
            if line.strip():
                try:
                    existing = json.loads(line)
                    if (
                        existing.get("session_id") == session_id
                        and existing.get("turn") == turn
                    ):
                        return
                except json.JSONDecodeError:
                    continue
    writer.write_episode(
        session_id=session_id,
        turn=turn,
        provider_name="chemsmart-runtime",
        model="semantic-gate-qmmm-v1",
        messages=[
            {"role": "user", "content": user},
            {"role": "assistant", "content": assistant},
        ],
        tool_records=tool_records,
        approvals_count=1,
        denials_count=0,
        cwd=str(REPO),
        available_tools=[
            "extract_project_protocol",
            "render_project_yaml",
            "validate_project_yaml",
            "repair_command",
            "synthesize_command",
        ],
        final_answer=assistant,
        workspace=_workspace(program),
    )


def _checked_command(
    command: str, program: str, request: str
) -> dict[str, Any]:
    from scripts.training import reasoning_accum

    with tempfile.TemporaryDirectory(prefix="qmmm-runtime-gold-") as raw:
        work = Path(raw)
        reasoning_accum._prepare_workspace(work, request, program)
        previous = Path.cwd()
        try:
            os.chdir(work)
            result = evaluate_command_semantics(command)
        finally:
            os.chdir(previous)
    semantic = _semantic_payload(result)
    if semantic["verdict"] not in {"ok", "warn"}:
        raise RuntimeError(
            f"runtime-gold command failed: {command}; "
            f"rules={semantic['failed_rule_ids']}"
        )
    if not semantic["generated_inputs"]:
        raise RuntimeError(
            f"runtime-gold command had no generated input: {command}"
        )
    return semantic


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="qmmm-runtime-gold-20260711")
    args = parser.parse_args()
    training_dir = REPO / "var" / "agent-training" / "runs" / args.run_id
    writer = TrainingEpisodeWriter(
        TrainingLogConfig(enabled=True, dir=training_dir)
    )

    gaussian_command = (
        "chemsmart run gaussian -p demo -f enzyme_gaussian_gold.xyz "
        "-c 0 -m 1 opt qmmm -ha 1-8 -ma 9-12 -la 13-16 "
        "-ct 0 -mt 1 -ch 0 -mh 1 -ci 0 -mi 1"
    )
    gaussian_request = (
        "three-layer Gaussian ONIOM optimization of enzyme_gaussian_gold.xyz "
        "with high atoms 1-8, medium atoms 9-12, low atoms 13-16, all "
        "neutral singlets"
    )
    gaussian_spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "gaussian.qmmm",
                "file": "enzyme_gaussian_gold.xyz",
                "charge": 0,
                "mult": 1,
                "settings": {
                    "parent_job": "opt",
                    "high_level_atoms": list(range(1, 9)),
                    "medium_level_atoms": list(range(9, 13)),
                    "low_level_atoms": list(range(13, 17)),
                    "charge_total": 0,
                    "mult_total": 1,
                    "charge_intermediate": 0,
                    "mult_intermediate": 1,
                    "charge_high": 0,
                    "mult_high": 1,
                },
            }
        ],
    }
    orca_command = (
        "chemsmart run orca -p demo -f enzyme_orca_gold.xyz -c -1 -m 1 "
        "opt qmmm -j QMMM -ha 1-12 -ct -1 -mt 1 -lm AMBER=HardFirst"
    )
    orca_request = (
        "ORCA opt qmmm of enzyme_orca_gold.xyz with high-level atoms 1-12, "
        "total charge -1, multiplicity 1, and AMBER=HardFirst"
    )
    orca_spec = {
        "intent": "workflow",
        "jobs": [
            {
                "id": 1,
                "kind": "orca.qmmm",
                "file": "enzyme_orca_gold.xyz",
                "charge": -1,
                "mult": 1,
                "settings": {
                    "parent_job": "opt",
                    "high_level_atoms": list(range(1, 13)),
                    "jobtype": "QMMM",
                    "low_level_method": "AMBER=HardFirst",
                    "charge_total": -1,
                    "mult_total": 1,
                },
            }
        ],
    }

    gaussian_semantic = _checked_command(
        gaussian_command, "gaussian", gaussian_request
    )
    orca_semantic = _checked_command(orca_command, "orca", orca_request)

    _append(
        writer,
        session_id="runtime-gold-qmmm-positive-gaussian-001",
        turn=1,
        program="gaussian",
        user=gaussian_request,
        assistant=_assistant(gaussian_command, gaussian_semantic),
        tool_records=[
            _synthesis_record(
                gaussian_command,
                gaussian_semantic,
                status="ready",
                explanation="Three-layer atom regions and layer state are explicit.",
                compact_spec=gaussian_spec,
            )
        ],
    )
    _append(
        writer,
        session_id="runtime-gold-qmmm-positive-orca-001",
        turn=1,
        program="orca",
        user=orca_request,
        assistant=_assistant(orca_command, orca_semantic),
        tool_records=[
            _synthesis_record(
                orca_command,
                orca_semantic,
                status="ready",
                explanation="ORCA QMMM parent, total state, QM region, and force field are explicit.",
                compact_spec=orca_spec,
            )
        ],
    )

    bad_gaussian = (
        "chemsmart run gaussian -p demo -f enzyme_gaussian_repair_gold.xyz "
        "-c 0 -m 1 opt qmmm -ha 1-8 -ma 9-12 -la 13-16"
    )
    good_gaussian = gaussian_command.replace(
        "enzyme_gaussian_gold.xyz", "enzyme_gaussian_repair_gold.xyz"
    )
    good_gaussian_spec = json.loads(
        json.dumps(gaussian_spec).replace(
            "enzyme_gaussian_gold.xyz", "enzyme_gaussian_repair_gold.xyz"
        )
    )
    bad_semantic = _semantic_payload(
        _evaluate_in_workspace(
            bad_gaussian,
            "gaussian",
            "three-layer Gaussian ONIOM repair fixture for "
            "enzyme_gaussian_repair_gold.xyz, atoms 1-16",
        )
    )
    if bad_semantic["verdict"] != "reject":
        raise RuntimeError(
            "repair fixture did not produce the expected reject"
        )
    good_semantic = _checked_command(
        good_gaussian,
        "gaussian",
        "three-layer Gaussian ONIOM repair fixture for "
        "enzyme_gaussian_repair_gold.xyz, atoms 1-16",
    )
    _append(
        writer,
        session_id="runtime-gold-qmmm-repair-gaussian-001",
        turn=1,
        program="gaussian",
        user=(
            "Attempt this Gaussian QMMM command and report the semantic gate; "
            f"do not repair it yet: {bad_gaussian}"
        ),
        assistant="The draft is rejected because the nested QMMM layer state is incomplete.",
        tool_records=[
            _synthesis_record(
                bad_gaussian,
                bad_semantic,
                status="needs_clarification",
                explanation="Nested QMMM charge/multiplicity fields are missing.",
            )
        ],
    )
    _append(
        writer,
        session_id="runtime-gold-qmmm-repair-gaussian-001",
        turn=2,
        program="gaussian",
        user=(
            "Repair the rejected command without changing its atom regions. "
            "Add the explicit total, high, and medium layer state and validate "
            "the generated input."
        ),
        assistant=_assistant(good_gaussian, good_semantic),
        tool_records=[
            _synthesis_record(
                good_gaussian,
                good_semantic,
                status="ready",
                explanation="The missing QMMM layer state is restored; route is unchanged.",
                compact_spec=good_gaussian_spec,
            )
        ],
    )

    route_command = orca_command.replace(
        "enzyme_orca_gold.xyz", "enzyme_orca_route_gold.xyz"
    )
    route_semantic = _checked_command(
        route_command,
        "orca",
        "ORCA QMMM YAML route correction fixture for "
        "enzyme_orca_route_gold.xyz, atoms 1-12",
    )
    route_spec = json.loads(
        json.dumps(orca_spec).replace(
            "enzyme_orca_gold.xyz", "enzyme_orca_route_gold.xyz"
        )
    )
    session_id = "runtime-gold-qmmm-route-orca-001"
    route_user = (
        "Create an ORCA project YAML for a QM/MM enzyme calculation and then "
        "prepare the actual job command for enzyme_orca_route_gold.xyz with "
        "high-level atoms 1-12, total charge -1, multiplicity 1."
    )
    _append(
        writer,
        session_id=session_id,
        turn=1,
        program="orca",
        user=route_user,
        assistant="I routed the first request to project-YAML authoring; no CLI command is claimed yet.",
        tool_records=[
            {
                "tool": "extract_project_protocol",
                "status": "ok",
                "args": {"text": route_user},
                "result": {
                    "status": "ready",
                    "facts": {"program": "orca", "kind": "qmmm"},
                },
            },
            {
                "tool": "render_project_yaml",
                "status": "ok",
                "args": {"project_name": "route-demo"},
                "result": {
                    "status": "ready",
                    "project_name": "route-demo",
                    "yaml_text": "gas: {}",
                },
            },
            {
                "tool": "validate_project_yaml",
                "status": "ok",
                "args": {"project_name": "route-demo"},
                "result": {"validation_verdict": "ok", "issues": []},
            },
        ],
    )
    _append(
        writer,
        session_id=session_id,
        turn=2,
        program="orca",
        user=(
            "Correction: the workspace demo project is already loaded. Do not "
            "author YAML; synthesize the actual ORCA opt qmmm CLI command and "
            "validate its generated input."
        ),
        assistant=_assistant(route_command, route_semantic),
        tool_records=[
            _synthesis_record(
                route_command,
                route_semantic,
                status="ready",
                explanation="The request is routed to CLI synthesis after the YAML clarification.",
                compact_spec=route_spec,
            )
        ],
    )

    variant_cases = [
        {
            "session_id": "runtime-gold-qmmm-variant-gaussian-ts-001",
            "program": "gaussian",
            "filename": "enzyme_gaussian_ts.xyz",
            "command": (
                "chemsmart run gaussian -p demo -f enzyme_gaussian_ts.xyz "
                "-c 0 -m 1 ts qmmm -ha 1-8 -ma 9-12 -la 13-16 "
                "-ct 0 -mt 1 -ch 0 -mh 1 -ci 0 -mi 1"
            ),
            "request": (
                "Gaussian three-layer ONIOM transition-state optimization of "
                "enzyme_gaussian_ts.xyz with high atoms 1-8, medium atoms 9-12, "
                "low atoms 13-16, all neutral singlets"
            ),
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.qmmm",
                        "file": "enzyme_gaussian_ts.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "parent_job": "ts",
                            "high_level_atoms": list(range(1, 9)),
                            "medium_level_atoms": list(range(9, 13)),
                            "low_level_atoms": list(range(13, 17)),
                            "charge_total": 0,
                            "mult_total": 1,
                            "charge_intermediate": 0,
                            "mult_intermediate": 1,
                            "charge_high": 0,
                            "mult_high": 1,
                        },
                    }
                ],
            },
        },
        {
            "session_id": "runtime-gold-qmmm-variant-gaussian-sp-001",
            "program": "gaussian",
            "filename": "enzyme_gaussian_sp.xyz",
            "command": (
                "chemsmart run gaussian -p demo -f enzyme_gaussian_sp.xyz "
                "-c 0 -m 1 sp qmmm -ha 1-8 -ma 9-12 -la 13-16 "
                "-ct 0 -mt 1 -ch 0 -mh 1 -ci 0 -mi 1"
            ),
            "request": (
                "Gaussian three-layer ONIOM single point for "
                "enzyme_gaussian_sp.xyz with high atoms 1-8, medium atoms "
                "9-12, low atoms 13-16, neutral singlet"
            ),
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "gaussian.qmmm",
                        "file": "enzyme_gaussian_sp.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "parent_job": "sp",
                            "high_level_atoms": list(range(1, 9)),
                            "medium_level_atoms": list(range(9, 13)),
                            "low_level_atoms": list(range(13, 17)),
                            "charge_total": 0,
                            "mult_total": 1,
                            "charge_intermediate": 0,
                            "mult_intermediate": 1,
                            "charge_high": 0,
                            "mult_high": 1,
                        },
                    }
                ],
            },
        },
        {
            "session_id": "runtime-gold-qmmm-variant-orca-ts-001",
            "program": "orca",
            "filename": "enzyme_orca_ts.xyz",
            "command": (
                "chemsmart run orca -p demo -f enzyme_orca_ts.xyz -c 0 -m 1 "
                "ts qmmm -j QMMM -ha 1-12 -ct 0 -mt 1 -lm AMBER=HardFirst"
            ),
            "request": (
                "ORCA QMMM transition-state optimization of enzyme_orca_ts.xyz "
                "with high atoms 1-12, total charge 0, multiplicity 1, and "
                "AMBER=HardFirst"
            ),
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "orca.qmmm",
                        "file": "enzyme_orca_ts.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "parent_job": "ts",
                            "high_level_atoms": list(range(1, 13)),
                            "jobtype": "QMMM",
                            "low_level_method": "AMBER=HardFirst",
                            "charge_total": 0,
                            "mult_total": 1,
                        },
                    }
                ],
            },
        },
        {
            "session_id": "runtime-gold-qmmm-variant-orca-sp-001",
            "program": "orca",
            "filename": "enzyme_orca_sp.xyz",
            "command": (
                "chemsmart run orca -p demo -f enzyme_orca_sp.xyz -c 0 -m 1 "
                "sp qmmm -j QMMM -ha 1-12 -ct 0 -mt 1 -lm AMBER=HardFirst"
            ),
            "request": (
                "ORCA QMMM single point for enzyme_orca_sp.xyz with high atoms "
                "1-12, total charge 0, multiplicity 1, and AMBER=HardFirst"
            ),
            "spec": {
                "intent": "workflow",
                "jobs": [
                    {
                        "id": 1,
                        "kind": "orca.qmmm",
                        "file": "enzyme_orca_sp.xyz",
                        "charge": 0,
                        "mult": 1,
                        "settings": {
                            "parent_job": "sp",
                            "high_level_atoms": list(range(1, 13)),
                            "jobtype": "QMMM",
                            "low_level_method": "AMBER=HardFirst",
                            "charge_total": 0,
                            "mult_total": 1,
                        },
                    }
                ],
            },
        },
    ]
    for case in variant_cases:
        semantic = _checked_command(
            case["command"], case["program"], case["request"]
        )
        _append(
            writer,
            session_id=case["session_id"],
            turn=1,
            program=case["program"],
            user=case["request"],
            assistant=_assistant(case["command"], semantic),
            tool_records=[
                _synthesis_record(
                    case["command"],
                    semantic,
                    status="ready",
                    explanation=(
                        "Distinct QMMM parent and layer topology is preserved "
                        "through the production semantic gate."
                    ),
                    compact_spec=case["spec"],
                )
            ],
        )

    print(
        json.dumps(
            {
                "training_dir": str(training_dir),
                "episodes": 10,
                "validated_positive_commands": 7,
                "repair_pairs_expected": 1,
                "wrong_route_contrasts_expected": 1,
                "commands": [gaussian_command, orca_command, route_command],
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0


def _evaluate_in_workspace(command: str, program: str, request: str) -> Any:
    from scripts.training import reasoning_accum

    with tempfile.TemporaryDirectory(prefix="qmmm-runtime-gold-") as raw:
        work = Path(raw)
        reasoning_accum._prepare_workspace(work, request, program)
        previous = Path.cwd()
        try:
            os.chdir(work)
            return evaluate_command_semantics(command)
        finally:
            os.chdir(previous)


if __name__ == "__main__":
    raise SystemExit(main())
