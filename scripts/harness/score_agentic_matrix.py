#!/usr/bin/env python3
"""Score repeated provider results against the frozen high-risk matrix."""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from chemsmart.agent.harness.evaluation import (
    classify_agent_result,
    load_case_matrix,
    reliability_metrics,
)


def score(matrix_path: Path, result_paths: list[Path]) -> dict[str, Any]:
    cases = {case.case_id: case for case in load_case_matrix(matrix_path)}
    rows: list[dict[str, Any]] = []
    outcomes: Counter[str] = Counter()
    by_provider: dict[str, list[dict[str, Any]]] = defaultdict(list)
    by_family: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for path in result_paths:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            raw = json.loads(line)
            case_id = str(raw.get("case_id") or "")
            case = cases.get(case_id)
            if case is None:
                raise ValueError(
                    f"result references unknown case_id: {case_id}"
                )
            outcome = classify_agent_result(
                status=str(raw.get("status") or ""),
                expected_outcome=case.expected_outcome,
                semantic_rule_ids=raw.get("semantic_rule_ids") or (),
                intent_rule_ids=raw.get("intent_rule_ids") or (),
                repaired=bool(raw.get("repaired")),
            )
            passed = outcome.value == case.expected_outcome or (
                case.expected_outcome == "direct_pass"
                and outcome.value in {"direct_pass", "repair_pass"}
            )
            row = {
                "case_id": case_id,
                "family": case.family,
                "provider": str(raw.get("provider") or "unknown"),
                "model": str(raw.get("model") or "unknown"),
                "trial": int(raw.get("trial") or 1),
                "outcome": outcome.value,
                "expected_outcome": case.expected_outcome,
                "passed": passed,
                "tokens": raw.get("tokens"),
                "latency_s": raw.get("latency_s"),
            }
            rows.append(row)
            outcomes[outcome.value] += 1
            provider_key = f"{row['provider']}:{row['model']}"
            by_provider[provider_key].append(row)
            by_family[case.family].append(row)
    return {
        "schema_version": 1,
        "matrix_case_count": len(cases),
        "result_count": len(rows),
        "outcomes": dict(outcomes.most_common()),
        "overall": reliability_metrics(rows, k=3),
        "by_provider": {
            key: reliability_metrics(value, k=3)
            for key, value in sorted(by_provider.items())
        },
        "by_family": {
            key: reliability_metrics(value, k=3)
            for key, value in sorted(by_family.items())
        },
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--matrix",
        type=Path,
        default=Path("tests/agent/harness/fixtures/high_risk_matrix.json"),
    )
    parser.add_argument("results", nargs="+", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = score(args.matrix, args.results)
    rendered = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
