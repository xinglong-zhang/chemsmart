#!/usr/bin/env python3
"""Replay exported review chains into the deterministic failure taxonomy."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path
from typing import Any

from chemsmart.agent.harness.failure_taxonomy import classify_runtime_failure


def replay(paths: list[Path]) -> dict[str, Any]:
    outcomes: Counter[str] = Counter()
    rules: Counter[str] = Counter()
    providers: Counter[str] = Counter()
    variants: Counter[str] = Counter()
    rows_seen = 0
    for path in paths:
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            rows_seen += 1
            meta = row.get("meta") if isinstance(row.get("meta"), dict) else {}
            provider = (
                meta.get("provider")
                if isinstance(meta.get("provider"), dict)
                else {}
            )
            provider_name = str(
                provider.get("model") or provider.get("name") or "unknown"
            )
            providers[provider_name] += 1
            variants[str(meta.get("schema_variant") or "unknown")] += 1
            row_rules, outcome = _classify_row(row)
            outcomes[outcome] += 1
            rules.update(row_rules or [f"outcome.{outcome}"])
    return {
        "schema_version": 1,
        "rows_seen": rows_seen,
        "outcomes": dict(outcomes.most_common()),
        "failed_rule_ids": dict(rules.most_common()),
        "providers": dict(providers.most_common()),
        "schema_variants": dict(variants.most_common()),
        "source_files": [str(path) for path in paths],
        "raw_ledger_modified": False,
    }


def _classify_row(row: dict[str, Any]) -> tuple[set[str], str]:
    rules: set[str] = set()
    statuses: list[str] = []
    tool_names = _tool_names(row)
    for message in row.get("messages") or []:
        if not isinstance(message, dict) or message.get("role") != "tool":
            continue
        try:
            payload = json.loads(str(message.get("content") or ""))
        except json.JSONDecodeError:
            rules.add("tool.result.non_json")
            continue
        if not isinstance(payload, dict):
            continue
        statuses.append(
            str(payload.get("status") or payload.get("verdict") or "")
        )
        semantic = (
            payload.get("semantic")
            if isinstance(payload.get("semantic"), dict)
            else {}
        )
        for issue in semantic.get("issues") or []:
            if not isinstance(issue, dict):
                continue
            rule_id = str(issue.get("rule_id") or issue.get("id") or "")
            if rule_id == "cmd.semantic.safe_execution_failed":
                evidence = (
                    issue.get("evidence")
                    if isinstance(issue.get("evidence"), dict)
                    else {}
                )
                failure = classify_runtime_failure(
                    stdout=str(evidence.get("stdout_tail") or ""),
                    stderr=str(evidence.get("stderr_tail") or ""),
                    returncode=evidence.get("returncode"),
                )
                rule_id = failure.rule_id
            if rule_id:
                rules.add(rule_id)
        for issue in payload.get("issues") or []:
            if isinstance(issue, dict):
                rule_id = str(issue.get("rule_id") or issue.get("id") or "")
                if rule_id:
                    rules.add(rule_id)
    skip = str((row.get("meta") or {}).get("skip_reason") or "")
    if any(rule.startswith("input.") for rule in rules):
        return rules, "generated_input_invariant_failure"
    if any(rule.startswith("intent.") for rule in rules):
        return rules, "intent_drift"
    if any("project" in rule or "yaml" in rule for rule in rules):
        return rules, "yaml_state_failure"
    if any("terminal" in rule or "server" in rule for rule in rules):
        return rules, "terminal_environment_failure"
    if skip == "unresolved_clarification" or "needs_clarification" in statuses:
        return rules, "unresolved_clarification"
    if skip == "unresolved_pause":
        return rules, "unresolved_pause"
    if skip in {"terminal_tool_error", "terminal_execute_failed"}:
        return rules, "cli_runtime_failure"
    if skip == "missing_generated_input_evidence":
        return rules, "generated_input_evidence_missing"
    if not tool_names:
        return rules, "no_tool_trajectory"
    return rules, skip or "unclassified_review"


def _tool_names(row: dict[str, Any]) -> list[str]:
    names: list[str] = []
    for message in row.get("messages") or []:
        if not isinstance(message, dict):
            continue
        for call in message.get("tool_calls") or []:
            if not isinstance(call, dict):
                continue
            function = (
                call.get("function")
                if isinstance(call.get("function"), dict)
                else {}
            )
            name = str(function.get("name") or "")
            if name:
                names.append(name)
    return names


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = replay(args.inputs)
    rendered = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
