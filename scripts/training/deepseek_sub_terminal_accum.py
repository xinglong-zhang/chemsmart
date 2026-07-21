"""Collect DeepSeek-authored ``chemsmart sub`` terminal-state trajectories.

The workspace contains only a mock PBS server YAML, fake Gaussian/ORCA
executables, and a marker-writing qsub replacement.  The real agent loop still
has to synthesize a valid ChemSmart CLI command, pass the semantic gate, invoke
the command, and produce a terminal receipt.  No real scheduler or chemistry
program is contacted.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from dotenv import dotenv_values

REPO = Path(__file__).resolve().parents[2]
BASE_URL = "https://api.deepseek.com/v1"

SCENARIOS = (
    (
        "submit_wrong_server_then_repair",
        (
            "Prepare and execute a Gaussian optimization of examples/h2o.xyz "
            "using the workspace project mock and server missing-pbs. Do not use "
            "build_job or submit_hpc directly; the user-facing artifact must "
            "be a chemsmart sub command. Charge and multiplicity are intentionally "
            "not given yet, so request them rather than guessing."
        ),
        "Correction: missing-pbs was a typo; H2O is neutral singlet (charge 0, "
        "multiplicity 1). Use workspace project mock and configured mock-pbs, "
        "keep Gaussian opt, then execute the corrected chemsmart sub command "
        "with test=false so the mock scheduler produces a terminal receipt.",
    ),
    (
        "submit_gaussian_opt_actual",
        (
            "Use workspace project mock and the configured mock-pbs server to execute a Gaussian geometry "
            "optimization for examples/h2o.xyz, neutral singlet. Generate and "
            "execute the chemsmart sub command, not a direct internal tool path."
            " Use execute_chemsmart_command with test=false so the mock scheduler "
            "is actually invoked."
        ),
        None,
    ),
    (
        "submit_orca_sp_actual",
        (
            "Use workspace project mock and the configured mock-pbs server to execute an ORCA single-point "
            "calculation for examples/h2o.xyz, neutral singlet. The final "
            "user-facing artifact must be a valid chemsmart sub command and it "
            "must be executed against the mock scheduler with test=false."
        ),
        "The ORCA command is valid. Continue the workflow by calling the execution tool now; do not stop after synthesis, and preserve the mock-pbs server and neutral singlet settings.",
    ),
    (
        "submit_missing_server_then_repair",
        (
            "I need a real ChemSmart CLI submission for a Gaussian optimization "
            "of examples/h2o.xyz, but the server name is not certain. Do not "
            "invent a server or use a direct tool path; ask for the missing "
            "server information."
        ),
        "The workspace project is mock and server is mock-pbs. Use them, keep neutral singlet Gaussian "
        "opt, and execute the valid chemsmart sub command.",
    ),
)


def _latest_terminal_state(
    episode_file: Path, session_id: str
) -> dict[str, Any] | None:
    latest = None
    if not episode_file.exists():
        return None
    for line in episode_file.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("session_id") == session_id and row.get("terminal_state"):
            latest = row["terminal_state"]
    return latest


def _last_command(output: dict[str, Any]) -> str:
    for outcome in reversed(output.get("tool_outcomes") or []):
        raw = getattr(outcome, "raw_result", None)
        if isinstance(raw, dict) and raw.get("command"):
            return str(raw["command"])
    return ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default="deepseek-v4-pro")
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument(
        "--batch-id", default="prod-readiness-deepseek-sub-terminal-20260711"
    )
    args = parser.parse_args()

    sys.path.insert(0, str(REPO))
    os.chdir(REPO)
    logging.disable(logging.INFO)
    from scripts.training.collect_sub_terminal_fixtures import _write_fixture

    run_dir = (REPO / "var" / "agent-training" / "runs" / args.model).resolve()
    batch_dir = run_dir / "sub_terminal_batches" / args.batch_id
    os.environ["CHEMSMART_AGENT_TRAINING_DIR"] = str(run_dir)
    harness_id = os.environ.get(
        "CHEMSMART_AGENT_HARNESS_ID", "uncommitted-current"
    )

    key = dotenv_values(REPO / "api.env").get("DEEPSEEK-api-key", "")
    if not key:
        raise SystemExit("api.env is missing DEEPSEEK-api-key")

    import chemsmart.agent.providers as providers_mod
    import chemsmart.agent.tools_command as tools_command
    from chemsmart.agent.core import AgentSession
    from chemsmart.agent.harness.terminal_state import (
        terminal_state_is_positive,
    )
    from chemsmart.agent.permissions import (
        ApprovalDecision,
        PermissionMode,
        PermissionPolicy,
    )
    from chemsmart.agent.providers import OpenAIProvider
    from chemsmart.agent.registry import ToolRegistry

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
    episode_file = (
        run_dir
        / "episodes"
        / f"{datetime.now(timezone.utc).strftime('%Y%m')}.jsonl"
    )
    results: list[dict[str, Any]] = []
    for index, (name, first, second) in enumerate(SCENARIOS[: args.limit], 1):
        # Each case gets an isolated workspace.  The resolver intentionally
        # refuses to choose between Gaussian and ORCA YAMLs, so sharing one
        # fixture would turn a valid case into an artificial clarification.
        program = "orca" if name.startswith("submit_orca") else "gaussian"
        paths = _write_fixture(
            batch_dir / "fixtures" / f"case-{index:02d}",
            program=program,
        )
        os.environ["HOME"] = str(paths["workspace"].parent / "home")
        tools_command.reset_command_tools_state()
        session = AgentSession(
            provider=provider,
            registry=ToolRegistry.default(
                groups=["synthesis", "project_yaml", "execution"]
            ),
            session_root=batch_dir / "sessions" / f"case-{index:02d}",
            stage_prompt="unified_agent.md",
        )
        policy = PermissionPolicy(
            mode=PermissionMode.DRIVING, prompt_risky=True
        )
        started = time.time()
        row: dict[str, Any] = {
            "batch_id": args.batch_id,
            "harness_id": harness_id,
            "scenario": name,
            "model": args.model,
        }
        original_cwd = Path.cwd()
        try:
            # The fixture owns the workspace YAML and relative structure path.
            # Running from the repository root silently bypasses both and
            # causes the teacher to enter a project-creation loop.
            os.chdir(paths["workspace"])
            out = session.run_loop(
                first,
                policy=policy,
                approver=lambda req: ApprovalDecision.ALLOW_ONCE,
            )
            if second:
                out = session.run_loop(
                    second,
                    policy=policy,
                    approver=lambda req: ApprovalDecision.ALLOW_ONCE,
                )
            state = _latest_terminal_state(episode_file, out["session_id"])
            command = _last_command(out)
            row.update(
                {
                    "session_id": out["session_id"],
                    "turn_count": 2 if second else 1,
                    "command": command,
                    "terminal_state": state,
                    "completed_steps": out.get("completed_steps"),
                    "ask_user": bool(out.get("ask_user_question")),
                    "latency_s": round(time.time() - started, 1),
                }
            )
            if (
                state
                and terminal_state_is_positive(state)
                and command.startswith("chemsmart sub")
            ):
                row["grade"] = "PASS_SUB_TERMINAL"
            elif out.get("ask_user_question"):
                row["grade"] = "ASK"
            else:
                row["grade"] = "WRONG"
        except Exception as exc:
            row.update(
                {"grade": "ERR", "error": f"{type(exc).__name__}: {exc}"[:240]}
            )
        finally:
            os.chdir(original_cwd)
        results.append(row)
        print(json.dumps(row, ensure_ascii=False, sort_keys=True), flush=True)

    result_path = run_dir / "sub_terminal_results.jsonl"
    with result_path.open("a", encoding="utf-8") as handle:
        for row in results:
            handle.write(
                json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n"
            )
    print(
        json.dumps(
            {"result_path": str(result_path), "count": len(results)},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
