"""Evaluate repeated agent workflow trials from observable result ledgers.

The evaluator deliberately consumes the small JSONL result files produced by
the DeepSeek sub-terminal and QMMM collectors.  A trial is positive only when
the collector's deterministic terminal or generated-input gate says so.  A
provider grade without evidence is never promoted to success.

``pass_at_1`` is the observed one-trial success rate.  ``pass_power_k`` is the
independence approximation ``pass_at_1 ** k`` used by tau-bench-style
reporting.  ``all_trials_pass`` is the stricter empirical result for every
trial in a group.  Groups with fewer than ``k`` trials are explicitly marked
as insufficient rather than being treated as reliable.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from chemsmart.agent.harness.terminal_state import terminal_state_is_positive

JsonDict = dict[str, Any]


def _read_rows(paths: Iterable[Path]) -> list[JsonDict]:
    rows: list[JsonDict] = []
    for path in paths:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            value = json.loads(line)
            if isinstance(value, dict):
                rows.append(value)
    return rows


def _is_positive(
    row: JsonDict,
    source: str,
    *,
    require_intent: bool = False,
) -> bool:
    if source == "sub_terminal":
        terminal_positive = (
            row.get("grade") == "PASS_SUB_TERMINAL"
            and terminal_state_is_positive(row.get("terminal_state"))
        )
        if not terminal_positive:
            return False
        if not require_intent:
            return True
        assertions = row.get("intent_assertions")
        return (
            isinstance(assertions, list)
            and bool(assertions)
            and all(item.get("status") == "pass" for item in assertions)
        )
    if source == "qmmm":
        evidence = row.get("generated_input_evidence")
        return (
            row.get("grade") == "PASS_QMMM"
            and row.get("semantic_verdict") in {"ok", "warn"}
            and isinstance(evidence, list)
            and bool(evidence)
            and str(row.get("command") or "").startswith("chemsmart ")
        )
    raise ValueError(f"unsupported source: {source}")


def _configuration_id(row: JsonDict) -> str:
    """Use an explicit harness id when present; otherwise isolate each batch."""

    return str(
        row.get("harness_id")
        or row.get("configuration_id")
        or row.get("batch_id")
        or "unknown"
    )


def evaluate(
    rows: Iterable[JsonDict],
    *,
    source: str,
    k: int = 3,
    harness_id: str | None = None,
    require_intent: bool = False,
) -> JsonDict:
    if k < 1:
        raise ValueError("k must be positive")
    groups: dict[tuple[str, str], list[JsonDict]] = defaultdict(list)
    excluded = 0
    for row in rows:
        if harness_id is not None and _configuration_id(row) != harness_id:
            excluded += 1
            continue
        scenario = str(row.get("scenario") or "unknown")
        groups[(scenario, _configuration_id(row))].append(row)

    reports: list[JsonDict] = []
    for (scenario, config), trials in sorted(groups.items()):
        positives = sum(
            _is_positive(row, source, require_intent=require_intent)
            for row in trials
        )
        trial_count = len(trials)
        pass_at_1 = positives / trial_count if trial_count else 0.0
        enough_trials = trial_count >= k
        reports.append(
            {
                "source": source,
                "scenario": scenario,
                "harness_id": config,
                "trials": trial_count,
                "positive_trials": positives,
                "failed_trials": trial_count - positives,
                "pass_at_1": round(pass_at_1, 6),
                "pass_power_k": round(pass_at_1**k, 6)
                if enough_trials
                else None,
                "all_trials_pass": positives == trial_count
                if enough_trials
                else None,
                "status": "measured" if enough_trials else "insufficient_trials",
                "batch_ids": sorted(
                    {str(row.get("batch_id") or "unknown") for row in trials}
                ),
                "session_ids": sorted(
                    {str(row.get("session_id") or "unknown") for row in trials}
                ),
            }
        )

    measured = [row for row in reports if row["status"] == "measured"]
    return {
        "schema_version": 1,
        "source": source,
        "k": k,
        "require_intent": require_intent,
        "groups": reports,
        "summary": {
            "groups": len(reports),
            "measured_groups": len(measured),
            "insufficient_groups": len(reports) - len(measured),
            "excluded_rows": excluded,
            "all_measured_groups_pass": bool(measured)
            and all(row["all_trials_pass"] for row in measured),
            "mean_pass_at_1": round(
                sum(row["pass_at_1"] for row in measured) / len(measured), 6
            )
            if measured
            else None,
            "mean_pass_power_k": round(
                sum(row["pass_power_k"] for row in measured) / len(measured), 6
            )
            if measured
            else None,
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sub-results", action="append", default=[])
    parser.add_argument("--qmmm-results", action="append", default=[])
    parser.add_argument("--out", required=True)
    parser.add_argument("--k", type=int, default=3)
    parser.add_argument("--harness-id")
    parser.add_argument(
        "--require-intent",
        action="store_true",
        help="For sub trials, require request-level intent assertions.",
    )
    args = parser.parse_args(argv)
    if not args.sub_results and not args.qmmm_results:
        parser.error("provide --sub-results and/or --qmmm-results")

    output: JsonDict = {"schema_version": 1, "k": args.k, "evaluations": []}
    if args.sub_results:
        output["evaluations"].append(
            evaluate(
                _read_rows(Path(path) for path in args.sub_results),
                source="sub_terminal",
                k=args.k,
                harness_id=args.harness_id,
                require_intent=args.require_intent,
            )
        )
    if args.qmmm_results:
        output["evaluations"].append(
            evaluate(
                _read_rows(Path(path) for path in args.qmmm_results),
                source="qmmm",
                k=args.k,
                harness_id=args.harness_id,
                require_intent=args.require_intent,
            )
        )
    Path(args.out).write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(output["evaluations"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
