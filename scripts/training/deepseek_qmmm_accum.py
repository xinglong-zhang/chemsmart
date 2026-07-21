"""Collect a small, grounded DeepSeek QMMM supplement.

This is intentionally a separate collector rather than a modification of the
older reasoning/agentic collectors.  It keeps the batch bounded at four
scenarios and asks the teacher for observable CLI decisions only.  The
production adapter and semantic gate remain the source of truth; provider
responses are never promoted without generated-input evidence.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import signal
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

from dotenv import dotenv_values

REPO = Path(__file__).resolve().parents[2]
BASE_URL = "https://api.deepseek.com/v1"

SCENARIOS = (
    {
        "name": "gaussian_qmmm_three_layer_positive",
        "program": "gaussian",
        "turns": (
            "The Gaussian project YAML demo is already loaded in this workspace. "
            "Prepare only a real chemsmart CLI command for a three-layer ONIOM "
            "optimization of enzyme_gaussian_supplement.xyz. The QM/high-level "
            "region is atoms 1-8, the medium-level region is atoms 9-12, the "
            "MM/low-level region is atoms 13-16, and the total system is charge "
            "0, multiplicity 1. The high and medium layers are neutral singlets. "
            "Do not create or modify YAML and do "
            "not use a direct internal tool path. Include the nested qmmm "
            "parent and return the command after semantic input validation.",
        ),
    },
    {
        "name": "orca_qmmm_two_layer_positive",
        "program": "orca",
        "turns": (
            "Use the already loaded ORCA project YAML demo. Prepare only a "
            "chemsmart CLI command for a QMMM geometry optimization of "
            "enzyme_orca_supplement.xyz. The total charge is -1 and total "
            "multiplicity is 1; the explicit high-level QM region is atoms "
            "1-12. The project supplies the high-level method and MM force "
            "field. Keep this as ORCA opt qmmm, not plain opt, and validate the "
            "generated input before returning the command.",
        ),
    },
    {
        "name": "gaussian_qmmm_missing_layer_state_repair",
        "program": "gaussian",
        "turns": (
            "Repair this failed Gaussian QM/MM command without changing the "
            "molecule or atom regions. The intended job is a three-layer ONIOM "
            "optimization for enzyme_gaussian_repair_supplement.xyz with high "
            "atoms 1-8, medium atoms 9-12, low atoms 13-16, neutral singlet "
            "total/high/medium system. The "
            "failed draft is:\n"
            "chemsmart run gaussian -p demo -f enzyme_gaussian_repair_supplement.xyz "
            "-c 0 -m 1 opt qmmm -ha 1-8 -ma 9-12 -la 13-16\n"
            "The nested qmmm layer state is missing. Use repair_command or the "
            "equivalent grounded repair path, and return the repaired command "
            "only after generated-input semantic validation.",
        ),
    },
    {
        "name": "orca_qmmm_yaml_route_contrast",
        "program": "orca",
        "turns": (
            "Create an ORCA project YAML for a QM/MM enzyme calculation and then "
            "prepare the actual job for enzyme_orca_route_supplement.xyz, with "
            "high-level atoms 1-10, total charge 0, multiplicity 1. This request "
            "is deliberately ambiguous: first show the project-YAML route and "
            "do not silently claim a calculation command.",
            "Correction: the demo project YAML is already loaded and must be "
            "used. Do not author or write YAML now. Route this as the actual "
            "ORCA opt qmmm command, include the explicit QMMM job type and "
            "high-level region, and return the grounded CLI command after input "
            "validation.",
        ),
    },
)


def _last_synthesis(output: dict[str, Any]) -> dict[str, Any] | None:
    for outcome in reversed(output.get("tool_outcomes") or []):
        if getattr(outcome, "name", "") == "synthesize_command":
            raw = getattr(outcome, "raw_result", None)
            if isinstance(raw, dict):
                return raw
    return None


def _semantic_fields(
    synthesis: dict[str, Any],
) -> tuple[str, list[dict[str, Any]]]:
    """Read normalized or raw synthesis payloads without losing gate evidence."""

    semantic = synthesis.get("semantic")
    if not isinstance(semantic, dict):
        semantic = {}
    verdict = str(
        synthesis.get("semantic_verdict") or semantic.get("verdict") or ""
    )
    evidence = synthesis.get("generated_input_evidence")
    if not isinstance(evidence, list):
        evidence = semantic.get("generated_inputs")
    return verdict, [item for item in evidence or [] if isinstance(item, dict)]


def _kind_present(command: str) -> bool:
    tokens = command.lower().replace("\\\n", " ").split()
    return "qmmm" in tokens or "oniom" in command.lower()


class _ScenarioTimeout(RuntimeError):
    pass


class _ScenarioDeadline:
    def __init__(self, seconds: float):
        self.seconds = seconds
        self.previous = None

    def __enter__(self):
        if self.seconds <= 0:
            return self
        self.previous = signal.getsignal(signal.SIGALRM)

        def raise_timeout(signum, frame):
            del signum, frame
            raise _ScenarioTimeout(f"scenario exceeded {self.seconds:g}s")

        signal.signal(signal.SIGALRM, raise_timeout)
        signal.setitimer(signal.ITIMER_REAL, self.seconds)
        return self

    def __exit__(self, exc_type, exc, tb):
        if self.seconds > 0:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, self.previous)
        return False


def _prepare_qmmm_workspace(
    reasoning_accum: Any,
    work: Path,
    request_text: str,
    program: str,
) -> None:
    """Prepare a project whose layers match the requested QMMM topology."""

    reasoning_accum._prepare_workspace(work, request_text, program)
    # The generic fixture intentionally carries a medium layer.  The prompts
    # above therefore use a complete three-layer topology rather than asking
    # the teacher to emit empty ``-mx/-mb`` overrides.


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="deepseek-v4-pro")
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument(
        "--batch-id",
        default="prod-readiness-deepseek-qmmm-supplement-20260711",
    )
    parser.add_argument("--timeout-s", type=float, default=150.0)
    args = parser.parse_args()
    if not 1 <= args.limit <= 4:
        raise SystemExit("--limit must be between 1 and 4")

    sys.path.insert(0, str(REPO))
    os.chdir(REPO)
    # Provider HTTP logs are not training evidence and make a bounded batch
    # difficult to audit in a terminal.  The API result is recorded through
    # the normal AgentSession ledger instead.
    logging.disable(logging.CRITICAL)

    import reasoning_accum

    import chemsmart.agent.providers as providers_mod
    import chemsmart.agent.tools_command as tools_command
    from chemsmart.agent.core import AgentSession
    from chemsmart.agent.permissions import (
        ApprovalDecision,
        PermissionMode,
        PermissionPolicy,
    )
    from chemsmart.agent.providers import OpenAIProvider
    from chemsmart.agent.registry import ToolRegistry

    key = dotenv_values(REPO / "api.env").get("DEEPSEEK-api-key", "")
    if not key:
        raise SystemExit("api.env is missing DEEPSEEK-api-key")

    class DeepSeekProvider(OpenAIProvider):
        def chat(
            self,
            messages: list[Any],
            tools: list[Any] | None = None,
            timeout_s: float = 30,
        ) -> dict[str, Any]:
            kwargs: dict[str, Any] = {
                "model": self.default_model,
                "messages": messages,
                "timeout": max(float(timeout_s), 90.0),
            }
            if tools:
                kwargs["tools"] = tools
            return self._client.chat.completions.create(**kwargs).model_dump()

    provider = DeepSeekProvider(key, model=args.model, base_url=BASE_URL)
    providers_mod.get_provider = lambda *a, **k: provider
    run_dir = REPO / "var" / "agent-training" / "runs" / args.model
    os.environ["CHEMSMART_AGENT_TRAINING_DIR"] = str(run_dir)
    harness_id = os.environ.get(
        "CHEMSMART_AGENT_HARNESS_ID", "uncommitted-current"
    )
    plan_dir = run_dir / "collection_plans"
    plan_dir.mkdir(parents=True, exist_ok=True)
    plan_path = plan_dir / f"{args.batch_id}.jsonl"
    with plan_path.open("w", encoding="utf-8") as handle:
        for index, scenario in enumerate(SCENARIOS[: args.limit], start=1):
            handle.write(
                json.dumps(
                    {
                        "batch_id": args.batch_id,
                        "case_index": index,
                        "model": args.model,
                        "scenario": scenario["name"],
                        "program": scenario["program"],
                        "turns": list(scenario["turns"]),
                        "workspace_yaml_policy": "one demo YAML per task workspace",
                        "quality_policy": (
                            "semantic ok|warn plus generated-input evidence; "
                            "otherwise review only"
                        ),
                    },
                    ensure_ascii=False,
                    sort_keys=True,
                )
                + "\n"
            )

    result_path = run_dir / "qmmm_supplement_results.jsonl"
    results: list[dict[str, Any]] = []
    for index, scenario in enumerate(SCENARIOS[: args.limit], start=1):
        work = Path(tempfile.mkdtemp(prefix=f"qmmm-supplement-{index:02d}-"))
        request_text = " ".join(scenario["turns"])
        _prepare_qmmm_workspace(
            reasoning_accum, work, request_text, scenario["program"]
        )
        original_cwd = Path.cwd()
        os.chdir(work)
        tools_command.reset_command_tools_state()
        started = time.time()
        row: dict[str, Any] = {
            "batch_id": args.batch_id,
            "case_index": index,
            "model": args.model,
            "harness_id": harness_id,
            "scenario": scenario["name"],
            "program": scenario["program"],
            "turn_count": len(scenario["turns"]),
        }
        try:
            session = AgentSession(
                provider=provider,
                registry=ToolRegistry.default(
                    groups=["synthesis", "project_yaml"]
                ),
                session_root=work / "sessions" / f"case-{index:02d}",
                stage_prompt="unified_agent.md",
            )
            policy = PermissionPolicy(
                mode=PermissionMode.DRIVING, prompt_risky=True
            )
            output: dict[str, Any] = {}
            turn_rows: list[dict[str, Any]] = []
            with _ScenarioDeadline(args.timeout_s):
                for turn_index, turn in enumerate(scenario["turns"], start=1):
                    output = session.run_loop(
                        turn,
                        policy=policy,
                        approver=lambda req: ApprovalDecision.ALLOW_ONCE,
                    )
                    synthesis = _last_synthesis(output) or {}
                    turn_verdict, turn_evidence = _semantic_fields(synthesis)
                    turn_rows.append(
                        {
                            "turn": turn_index,
                            "status": synthesis.get("status"),
                            "semantic_verdict": turn_verdict,
                            "command": synthesis.get("command", ""),
                            "generated_input_evidence": turn_evidence,
                            "ask_user": bool(output.get("ask_user_question")),
                        }
                    )
            synthesis = _last_synthesis(output) or {}
            command = str(synthesis.get("command") or "").strip()
            verdict, generated_input_evidence = _semantic_fields(synthesis)
            row.update(
                {
                    "session_id": output.get("session_id"),
                    "turns": turn_rows,
                    "command": command,
                    "semantic_verdict": verdict,
                    "generated_input_evidence": generated_input_evidence,
                    "invoked_tools": output.get("invoked_tools") or [],
                    "ask_user": bool(output.get("ask_user_question")),
                    "latency_s": round(time.time() - started, 1),
                }
            )
            if (
                command.startswith("chemsmart ")
                and _kind_present(command)
                and verdict in {"ok", "warn"}
                and row["generated_input_evidence"]
            ):
                row["grade"] = "PASS_QMMM"
            elif row["ask_user"]:
                row["grade"] = "ASK"
            else:
                row["grade"] = "WRONG"
        except Exception as exc:  # pragma: no cover - live API path
            row.update(
                {
                    "grade": "ERR",
                    "error": f"{type(exc).__name__}: {exc}"[:300],
                    "latency_s": round(time.time() - started, 1),
                }
            )
        finally:
            os.chdir(original_cwd)
        results.append(row)
        print(json.dumps(row, ensure_ascii=False, sort_keys=True), flush=True)

    with result_path.open("a", encoding="utf-8") as handle:
        for row in results:
            handle.write(
                json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n"
            )
    print(
        json.dumps(
            {
                "batch_id": args.batch_id,
                "plan_path": str(plan_path),
                "result_path": str(result_path),
                "count": len(results),
                "grades": {
                    grade: sum(row.get("grade") == grade for row in results)
                    for grade in sorted({row.get("grade") for row in results})
                },
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
