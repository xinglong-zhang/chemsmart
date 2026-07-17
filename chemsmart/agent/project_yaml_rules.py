"""Deterministic static and protocol-alignment rules for project YAML."""

from __future__ import annotations

from typing import Any

from chemsmart.agent.harness.basis_sets import resolve_basis_name
from chemsmart.agent.project_protocol import render_method_block
from chemsmart.agent.project_yaml_values import string_list, string_or_none


def static_project_yaml_issues(
    parsed: dict[str, Any],
    program: str,
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    if "gas" not in parsed and "solv" not in parsed:
        issues.append(
            issue(
                "yaml.project_phase_missing",
                "reject",
                "project YAML must define at least one of gas or solv.",
            )
        )
    if "gas" in parsed and "solv" not in parsed:
        issues.append(
            issue(
                "yaml.solv_block_required_for_loader",
                "reject",
                (
                    "current chemsmart loader requires solv when gas is "
                    "present so sp settings can be built."
                ),
            )
        )
    allowed_top = {"gas", "solv", "td", "qmmm"}
    for key in sorted(set(parsed) - allowed_top):
        issues.append(
            issue(
                "yaml.unknown_top_level_key",
                "warn",
                (
                    f"unknown top-level key {key!r} will not be used by the "
                    "project settings loader."
                ),
            )
        )
    issues.extend(_phase_issues(parsed, program))
    if program == "gaussian":
        issues.extend(_gaussian_static_issues(parsed))
    return issues


def protocol_alignment_issues(
    parsed: Any,
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    if not isinstance(parsed, dict):
        return []
    method = protocol.get("method")
    if not isinstance(method, dict):
        return []
    block = parsed.get("gas") if isinstance(parsed.get("gas"), dict) else {}
    issues = _method_alignment_issues(block, method)
    td_method = protocol.get("td")
    if isinstance(td_method, dict):
        issues.extend(_td_alignment_issues(parsed, td_method))
    return issues


def issue(rule_id: str, severity: str, message: str) -> dict[str, Any]:
    return {"rule_id": rule_id, "severity": severity, "message": message}


def _phase_issues(
    parsed: dict[str, Any],
    program: str,
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    for phase in ("gas", "solv", "td", "qmmm"):
        block = parsed.get(phase)
        if block is None:
            continue
        if not isinstance(block, dict):
            issues.append(
                issue(
                    "yaml.phase_not_mapping",
                    "reject",
                    f"{phase} must be a mapping or null.",
                )
            )
        else:
            issues.extend(_basis_catalog_issues(block, phase, program))
    return issues


def _basis_catalog_issues(
    block: dict[str, Any],
    phase: str,
    program: str,
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    for key in (
        "basis",
        "heavy_elements_basis",
        "light_elements_basis",
        "aux_basis",
    ):
        value = string_or_none(block.get(key))
        if value is None or value.lower() in {"gen", "genecp"}:
            continue
        result = resolve_basis_name(value, program=program)
        if result.verdict == "ok":
            continue
        issues.append(
            issue(
                "yaml.basis.program_unsupported"
                if result.canonical_name
                else "yaml.basis.unrecognized",
                "reject",
                f"{phase}.{key}={value!r}: {result.message}",
            )
        )
    return issues


def _gaussian_static_issues(
    parsed: dict[str, Any],
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    for phase in ("gas", "solv", "td"):
        block = parsed.get(phase)
        if not isinstance(block, dict):
            continue
        issues.extend(_gaussian_method_issues(block, phase))
    return issues


def _gaussian_method_issues(
    block: dict[str, Any],
    phase: str,
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    basis = string_or_none(block.get("basis"))
    if string_or_none(block.get("functional")) is None:
        issues.append(
            issue(
                "yaml.method_missing_functional",
                "reject",
                f"{phase} must define functional for a generated project YAML.",
            )
        )
    if basis is None:
        issues.append(
            issue(
                "yaml.method_missing_basis",
                "reject",
                f"{phase} must define basis for a generated project YAML.",
            )
        )
    heavy_basis = string_or_none(block.get("heavy_elements_basis"))
    light_basis = string_or_none(block.get("light_elements_basis"))
    heavy_elements = string_list(block.get("heavy_elements"))
    if (heavy_basis or light_basis or heavy_elements) and basis not in {
        "gen",
        "genecp",
    }:
        issues.append(
            issue(
                "yaml.gaussian.mixed_basis_without_gen",
                "reject",
                (
                    f"{phase} defines mixed-basis fields but basis is not "
                    "gen/genecp."
                ),
            )
        )
    if basis in {"gen", "genecp"} and not (heavy_basis or light_basis):
        issues.append(
            issue(
                "yaml.gaussian.gen_without_basis_sections",
                "reject",
                (
                    f"{phase} uses {basis} but does not define heavy/light "
                    "basis sections."
                ),
            )
        )
    return issues


def _method_alignment_issues(
    block: dict[str, Any],
    method: dict[str, Any],
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    expected_route = string_or_none(method.get("functional_route"))
    if expected_route and block.get("functional") != expected_route:
        issues.append(
            issue(
                "critic.functional_mismatch",
                "reject",
                (
                    f"gas.functional should be {expected_route!r} for the "
                    "reported method."
                ),
            )
        )
    expected_basis = string_or_none(method.get("basis"))
    if expected_basis and block.get("basis") != expected_basis:
        issues.append(
            issue(
                "critic.basis_mismatch",
                "reject",
                (
                    f"gas.basis should be {expected_basis!r} for the reported "
                    "method."
                ),
            )
        )
    expected_heavy = string_list(method.get("heavy_elements"))
    if (
        expected_heavy
        and string_list(block.get("heavy_elements")) != expected_heavy
    ):
        issues.append(
            issue(
                "critic.heavy_elements_mismatch",
                "reject",
                f"gas.heavy_elements should be {expected_heavy!r}.",
            )
        )
    issues.extend(_basis_alignment_issues(block, method))
    if method.get("freq") is True and block.get("freq") is not True:
        issues.append(
            issue(
                "critic.freq_missing",
                "warn",
                (
                    "reported minima/TS confirmation requires harmonic "
                    "frequency analysis; set gas.freq: true."
                ),
            )
        )
    return issues


def _basis_alignment_issues(
    block: dict[str, Any],
    method: dict[str, Any],
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    for key in ("heavy_elements_basis", "light_elements_basis"):
        expected = string_or_none(method.get(key))
        if expected and block.get(key) != expected:
            issues.append(
                issue(
                    f"critic.{key}_mismatch",
                    "reject",
                    f"gas.{key} should be {expected!r}.",
                )
            )
    return issues


def _td_alignment_issues(
    parsed: dict[str, Any],
    td_method: dict[str, Any],
) -> list[dict[str, Any]]:
    td_block = parsed.get("td") if isinstance(parsed.get("td"), dict) else {}
    expected_td = render_method_block(td_method, "gaussian")
    if not td_block:
        return [
            issue(
                "critic.td_block_missing",
                "reject",
                "reported TD-DFT method requires a top-level td block.",
            )
        ]
    issues: list[dict[str, Any]] = []
    for key in ("functional", "basis"):
        if td_block.get(key) != expected_td.get(key):
            issues.append(
                issue(
                    f"critic.td_{key}_mismatch",
                    "reject",
                    (
                        f"td.{key} should be "
                        f"{expected_td.get(key)!r} for the reported TD-DFT method."
                    ),
                )
            )
    return issues


__all__ = ["issue", "protocol_alignment_issues", "static_project_yaml_issues"]
