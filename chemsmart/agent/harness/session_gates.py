"""Compose deterministic runtime and generated-input gates for a Plan."""

from __future__ import annotations

import re
from typing import Any, Protocol

from chemsmart.agent.harness.models import HarnessResult
from chemsmart.agent.harness.runner import evaluate_harness
from chemsmart.agent.models import CriticVerdict, Plan, Step
from chemsmart.agent.services.result_codec import REFERENCE_RE

_GAUSSIAN_ROUTE_RE = re.compile(r"^\s*#\s*\S+", re.MULTILINE)
_ORCA_ROUTE_RE = re.compile(r"^\s*!\s*\S+", re.MULTILINE)


class GateHost(Protocol):
    decision_log: Any | None
    _last_harness_result: HarnessResult | None


def apply_deterministic_gates(
    host: GateHost,
    *,
    plan: Plan,
    verdict: CriticVerdict,
    runtime_result: dict[str, Any] | None,
    dry_run_results: list[dict[str, Any]],
    preview_submit: dict[str, Any] | None,
    dry_submit: bool,
) -> CriticVerdict:
    issues = filter_critic_issues(
        plan=plan, issues=verdict.issues, dry_submit=dry_submit
    )
    rationale = [verdict.rationale] if verdict.rationale else []
    final = verdict.verdict
    final = _apply_runtime_result(
        final, issues, rationale, runtime_result, dry_submit
    )
    final = _apply_route_checks(
        final, issues, rationale, plan, dry_run_results
    )
    final = _apply_software_harness(
        host, final, issues, rationale, plan, dry_run_results
    )
    final = _apply_duplicate_check(final, issues, rationale, preview_submit)
    if final == "warn" and not issues:
        final = "ok"
    return CriticVerdict(
        verdict=final,
        confidence=verdict.confidence,
        issues=dedupe_strings(issues),
        rationale="; ".join(part for part in rationale if part),
    )


def _apply_runtime_result(
    verdict: str,
    issues: list[str],
    rationale: list[str],
    runtime_result: dict[str, Any] | None,
    dry_submit: bool,
) -> str:
    if runtime_result is None:
        return verdict
    runtime_ok = runtime_result.get("ok")
    if runtime_ok == "fail":
        issues.extend(runtime_result.get("local_issues", []))
        rationale.append("validate_runtime returned fail")
        return "reject"
    if runtime_ok == "partial" and not dry_submit:
        issues.extend(runtime_result.get("remote_unknown", []))
        rationale.append("validate_runtime returned partial")
        return "warn" if verdict == "ok" else verdict
    return verdict


def _apply_route_checks(
    verdict: str,
    issues: list[str],
    rationale: list[str],
    plan: Plan,
    dry_run_results: list[dict[str, Any]],
) -> str:
    malformed = [
        issue
        for result in dry_run_results
        if (issue := malformed_input_issue(result)) is not None
    ]
    if malformed:
        issues.extend(malformed)
        rationale.append("dry-run input route line failed basic validation")
        if verdict != "reject":
            verdict = "warn"
    irc_issues = missing_irc_keyword_issues(plan, dry_run_results)
    if irc_issues:
        issues.extend(irc_issues)
        rationale.append("IRC input route line missing required keyword")
        verdict = "reject"
    return verdict


def _apply_software_harness(
    host: GateHost,
    verdict: str,
    issues: list[str],
    rationale: list[str],
    plan: Plan,
    dry_run_results: list[dict[str, Any]],
) -> str:
    result = evaluate_harness(plan=plan, dry_run_results=dry_run_results)
    host._last_harness_result = result
    if host.decision_log is not None:
        host.decision_log.write(
            "harness_result",
            result.to_dict(),
            rationale=f"runtime harness verdict: {result.verdict}",
        )
    if result.verdict == "reject":
        issues.extend(harness_issue_messages(result))
        rationale.append("software invariant harness rejected generated input")
        return "reject"
    if result.verdict == "warn":
        issues.extend(harness_issue_messages(result))
        rationale.append(
            "software invariant harness warned on generated input"
        )
        return "warn" if verdict == "ok" else verdict
    return verdict


def _apply_duplicate_check(
    verdict: str,
    issues: list[str],
    rationale: list[str],
    preview_submit: dict[str, Any] | None,
) -> str:
    if preview_submit is None:
        return verdict
    duplicate = preview_submit.get("duplicate_check", {})
    if not duplicate.get("duplicate"):
        return verdict
    issues.append(duplicate.get("message") or "duplicate submission detected")
    rationale.append("submit_hpc duplicate check rejected the plan")
    return "reject"


def block_reason(
    verdict: CriticVerdict,
    dry_submit: bool,
    allow_remote_unknown: bool,
    allow_critic_override: bool,
) -> str | None:
    if verdict.verdict == "reject":
        return "critic_reject"
    if verdict.verdict == "ok":
        return None
    has_remote = any(is_remote_unknown_issue(i) for i in verdict.issues)
    has_other = any(not is_remote_unknown_issue(i) for i in verdict.issues)
    if not verdict.issues:
        has_other = True
    if has_remote and not dry_submit and not allow_remote_unknown:
        return "critic_warn_remote_unknown"
    if has_other and not allow_critic_override:
        return "critic_warn_no_override"
    return None


def filter_critic_issues(
    *, plan: Plan, issues: list[str], dry_submit: bool
) -> list[str]:
    needs_handoff = plan_requires_geometry_handoff(plan)
    filtered: list[str] = []
    for issue in issues:
        if dry_submit and is_remote_unknown_issue(issue):
            continue
        if "geometry handoff missing" in issue.lower() and not needs_handoff:
            continue
        filtered.append(issue)
    return filtered


def plan_requires_geometry_handoff(plan: Plan) -> bool:
    for index, step in enumerate(plan.steps):
        if step.tool != "build_job":
            continue
        kind = step.args.get("kind")
        if not (isinstance(kind, str) and kind.endswith(".sp")):
            continue
        source = source_step_for_ref(plan, step.args.get("molecule"))
        if source is None or source.tool != "build_molecule":
            continue
        if any(
            _is_geometry_source(candidate) for candidate in plan.steps[:index]
        ):
            return True
    return False


def _is_geometry_source(step: Step) -> bool:
    kind = step.args.get("kind")
    return (
        step.tool == "build_job"
        and isinstance(kind, str)
        and kind.endswith((".opt", ".ts", ".irc"))
    )


def source_step_for_ref(plan: Plan, value: Any) -> Step | None:
    if not isinstance(value, str):
        return None
    match = REFERENCE_RE.match(value)
    if match is None or match.group("path"):
        return None
    index = int(match.group("index")) - 1
    return plan.steps[index] if 0 <= index < len(plan.steps) else None


def is_remote_unknown_issue(issue: str) -> bool:
    return (
        issue.startswith("server.")
        or issue.endswith("on HPC")
        or issue == "ssh login reachable"
    )


def malformed_input_issue(result: dict[str, Any] | None) -> str | None:
    if not result:
        return None
    content = result.get("content")
    inputfile = str(result.get("inputfile", "")).lower()
    if not isinstance(content, str):
        return None
    if inputfile.endswith(".inp"):
        return (
            None
            if _ORCA_ROUTE_RE.search(content)
            else ("ORCA route line missing or malformed")
        )
    return (
        None
        if _GAUSSIAN_ROUTE_RE.search(content)
        else ("Gaussian route line missing or malformed")
    )


def missing_irc_keyword_issues(
    plan: Plan, dry_run_results: list[dict[str, Any]]
) -> list[str]:
    issues: list[str] = []
    steps = [step for step in plan.steps if step.tool == "dry_run_input"]
    for step, result in zip(steps, dry_run_results):
        kind = dry_run_job_kind(plan, step)
        route = route_text(result)
        if kind == "gaussian.irc" and (
            route is None or re.search(r"\birc\s*=", route, re.I) is None
        ):
            issues.append("Gaussian IRC input missing irc= keyword")
        elif kind == "orca.irc" and (
            route is None or re.search(r"\birc\b", route, re.I) is None
        ):
            issues.append("ORCA IRC input missing IRC keyword")
    return issues


def dry_run_job_kind(plan: Plan, dry_run_step: Step) -> str | None:
    job_ref = dry_run_step.args.get("job")
    if not isinstance(job_ref, str):
        return None
    match = REFERENCE_RE.match(job_ref)
    if match is None:
        return None
    index = int(match.group("index")) - 1
    if index < 0 or index >= len(plan.steps):
        return None
    job_step = plan.steps[index]
    kind = job_step.args.get("kind")
    if job_step.tool != "build_job" or not isinstance(kind, str):
        return None
    return kind.strip().lower()


def route_text(result: dict[str, Any] | None) -> str | None:
    if not result or not isinstance(result.get("content"), str):
        return None
    content = result["content"]
    if str(result.get("inputfile", "")).lower().endswith(".inp"):
        matches = re.findall(r"^\s*!\s*(.+)$", content, re.MULTILINE)
    else:
        matches = re.findall(r"^\s*#.*$", content, re.MULTILINE)
    return " ".join(match.strip() for match in matches) if matches else None


def dedupe_strings(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def harness_issue_messages(result: HarnessResult) -> list[str]:
    return [f"{issue.rule_id}: {issue.message}" for issue in result.issues]


__all__ = [
    "apply_deterministic_gates",
    "block_reason",
    "dedupe_strings",
    "dry_run_job_kind",
    "filter_critic_issues",
    "harness_issue_messages",
    "is_remote_unknown_issue",
    "malformed_input_issue",
    "missing_irc_keyword_issues",
    "plan_requires_geometry_handoff",
    "route_text",
    "source_step_for_ref",
]
