"""Project YAML builder, validator, and critic tools for the agent."""

from __future__ import annotations

import json
import tempfile
from copy import deepcopy
from difflib import unified_diff
from pathlib import Path
from typing import Any, Literal

import yaml

from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_workspace_project,
)
from chemsmart.agent.project_protocol import (
    extract_project_protocol,
    render_project_document,
)
from chemsmart.agent.project_yaml_rules import issue as _issue
from chemsmart.agent.project_yaml_rules import (
    protocol_alignment_issues as _protocol_alignment_issues,
)
from chemsmart.agent.project_yaml_rules import (
    static_project_yaml_issues as _static_project_yaml_issues,
)
from chemsmart.agent.project_yaml_values import (
    normalize_program as _normalize_program,
)
from chemsmart.agent.project_yaml_values import (
    normalize_project_name as _normalize_project_name,
)
from chemsmart.agent.project_yaml_values import (
    normalize_yaml_text as _normalize_yaml_text,
)
from chemsmart.settings.workspace_project import (
    resolve_workspace_project,
    workspace_project_path,
)

ProjectProgram = Literal["gaussian", "orca"]


def render_project_yaml(
    protocol: dict[str, Any],
    project_name: str | None = None,
    program: ProjectProgram = "gaussian",
) -> dict[str, Any]:
    """Render a chemsmart project YAML candidate from extracted protocol facts."""

    candidate = render_project_document(protocol, project_name, program)
    normalized_program = str(candidate["program"])
    name = str(candidate["project_name"])
    yaml_text = str(candidate["yaml_text"])
    validation = validate_project_yaml(
        yaml_text,
        program=normalized_program,
        project_name=name,
    )
    return {
        "ok": validation["verdict"] in {"ok", "warn"},
        "project_name": name,
        "program": normalized_program,
        "yaml_text": yaml_text,
        "validation": validation,
        "unsupported_yaml_features": candidate["unsupported_yaml_features"],
    }


def _coerce_yaml_text(value: Any) -> str:
    """Accept either a YAML string or a render_project_yaml result dict.

    The tool-loop model frequently chains ``render_project_yaml`` straight into
    ``validate_project_yaml`` and passes the whole render result object. Unwrap
    the ``yaml_text`` field so the harness accepts what the model naturally
    produces instead of raising an opaque type error.
    """

    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        for key in ("yaml_text", "yaml"):
            inner = value.get(key)
            if isinstance(inner, str) and inner.strip():
                return inner
    raise ValueError(
        "yaml_text must be a YAML string or a render_project_yaml result "
        "containing a 'yaml_text' field"
    )


_VALIDATION_CACHE: dict[tuple[str, str, str], dict[str, Any]] = {}


def validate_project_yaml(
    yaml_text: str | dict[str, Any],
    program: ProjectProgram = "gaussian",
    project_name: str = "candidate",
) -> dict[str, Any]:
    """Validate project YAML by loading it through chemsmart project settings.

    Results are memoized per (yaml_text, program, project_name) so the build-mode
    tool loop can re-request validation of an unchanged candidate without
    repeating the runtime loader work (dedup guard against re-validation loops).
    """

    yaml_text = _coerce_yaml_text(yaml_text)
    normalized_program = _normalize_program(program)
    name = _normalize_project_name(project_name)
    cache_key = (yaml_text, normalized_program, name)
    cached = _VALIDATION_CACHE.get(cache_key)
    if cached is not None:
        repeat = deepcopy(cached)
        # Signal the tool loop that this identical candidate was already
        # validated, so it stops re-validating and reports the result.
        repeat["revalidation_skipped"] = True
        return repeat
    result = _validate_project_yaml_uncached(
        yaml_text, normalized_program, name
    )
    if len(_VALIDATION_CACHE) >= 256:
        _VALIDATION_CACHE.clear()
    _VALIDATION_CACHE[cache_key] = deepcopy(result)
    return deepcopy(result)


def _validate_project_yaml_uncached(
    yaml_text: str,
    normalized_program: str,
    name: str,
) -> dict[str, Any]:
    issues: list[dict[str, Any]] = []
    try:
        parsed = yaml.safe_load(yaml_text) or {}
    except yaml.YAMLError as exc:
        return _validation_result(
            ok=False,
            verdict="reject",
            program=normalized_program,
            project_name=name,
            issues=[
                _issue(
                    "yaml.parse",
                    "reject",
                    f"YAML could not be parsed: {exc}",
                )
            ],
        )
    if not isinstance(parsed, dict):
        issues.append(
            _issue("yaml.root", "reject", "project YAML must be a mapping")
        )
    else:
        issues.extend(_static_project_yaml_issues(parsed, normalized_program))

    if any(issue["severity"] == "reject" for issue in issues):
        return _validation_result(
            ok=False,
            verdict="reject",
            program=normalized_program,
            project_name=name,
            parsed=parsed,
            issues=issues,
        )

    try:
        summary = _load_project_yaml_via_runtime(
            yaml_text=yaml_text,
            program=normalized_program,
            project_name=name,
        )
    except Exception as exc:
        issues.append(
            _issue(
                "yaml.runtime_load",
                "reject",
                f"chemsmart project settings loader rejected YAML: {exc}",
            )
        )
        return _validation_result(
            ok=False,
            verdict="reject",
            program=normalized_program,
            project_name=name,
            parsed=parsed,
            issues=issues,
        )

    verdict = "warn" if issues else "ok"
    return _validation_result(
        ok=verdict == "ok",
        verdict=verdict,
        program=normalized_program,
        project_name=name,
        parsed=parsed,
        runtime_summary=summary,
        issues=issues,
    )


def critic_project_yaml(
    yaml_text: str | dict[str, Any],
    protocol: dict[str, Any] | None = None,
    program: ProjectProgram = "gaussian",
    project_name: str = "candidate",
) -> dict[str, Any]:
    """Critique whether YAML matches a literature protocol and chemsmart use."""

    yaml_text = _coerce_yaml_text(yaml_text)
    validation = validate_project_yaml(
        yaml_text=yaml_text,
        program=program,
        project_name=project_name,
    )
    issues = list(validation.get("issues") or [])
    try:
        parsed = yaml.safe_load(yaml_text) or {}
    except yaml.YAMLError:
        parsed = {}
    if isinstance(protocol, dict):
        issues.extend(_protocol_alignment_issues(parsed, protocol))
        for feature in protocol.get("unsupported_yaml_features") or []:
            issues.append(
                _issue(
                    "protocol.unsupported_yaml_feature",
                    "warn",
                    f"{feature} is part of the reported protocol but is not a chemsmart project-YAML field.",
                )
            )

    severities = {issue.get("severity") for issue in issues}
    verdict = (
        "reject" if "reject" in severities else "warn" if issues else "ok"
    )
    return {
        "ok": verdict == "ok",
        "verdict": verdict,
        "project_name": project_name,
        "program": _normalize_program(program),
        "issues": issues,
        "validation": validation,
        "summary": _critic_summary(verdict, issues),
    }


def write_project_yaml(
    project_name: str,
    yaml_text: str | dict[str, Any],
    program: ProjectProgram = "gaussian",
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write a validated project YAML into the current workspace."""

    yaml_text = _coerce_yaml_text(yaml_text)
    normalized_program = _normalize_program(program)
    name = _normalize_project_name(project_name)
    validation = validate_project_yaml(
        yaml_text=yaml_text,
        program=normalized_program,
        project_name=name,
    )
    if validation["verdict"] == "reject":
        return {
            "ok": False,
            "project_name": name,
            "program": normalized_program,
            "written_path": None,
            "validation": validation,
            "error": "project YAML failed validation",
        }
    target = workspace_project_path(name, normalized_program)
    if target.exists() and not overwrite:
        return {
            "ok": False,
            "project_name": name,
            "program": normalized_program,
            "written_path": str(target),
            "validation": validation,
            "error": "target exists; pass overwrite=True to replace it",
        }
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(_normalize_yaml_text(yaml_text), encoding="utf-8")
    state_delta = select_workspace_project(name, normalized_program)
    return {
        "ok": True,
        "project_name": name,
        "program": normalized_program,
        "written_path": str(target),
        "overwrite": overwrite,
        "validation": validation,
        "state_delta": {"project": state_delta},
    }


def read_project_yaml(
    project_name: str = "",
    program: str = "",
) -> dict[str, Any]:
    """Read the active workspace project YAML and summarize runtime settings."""

    resolved = _resolve_project_yaml_target(project_name, program)
    if resolved["path"] is None:
        return {
            "ok": False,
            "project_name": resolved["project_name"],
            "program": resolved["program"],
            "path": None,
            "candidates": resolved["candidates"],
            "message": resolved["message"],
        }
    path = Path(str(resolved["path"]))
    yaml_text = path.read_text(encoding="utf-8")
    validation = validate_project_yaml(
        yaml_text,
        program=resolved["program"],
        project_name=resolved["project_name"],
    )
    state_delta = (
        select_workspace_project(resolved["project_name"], resolved["program"])
        if validation.get("verdict") != "reject"
        else {"selected": False, "rule_id": "workflow.project.invalid"}
    )
    return {
        "ok": True,
        "project_name": resolved["project_name"],
        "program": resolved["program"],
        "path": str(path),
        "yaml_text": yaml_text,
        "parsed": yaml.safe_load(yaml_text) or {},
        "validation": validation,
        "runtime_summary": validation.get("runtime_summary"),
        "state_delta": {"project": state_delta},
    }


def update_project_yaml(
    updates: dict[str, Any] | str | None = None,
    project_name: str = "",
    program: str = "",
    unset: list[str] | None = None,
) -> dict[str, Any]:
    """Patch an existing workspace project YAML after validation."""

    updates, decode_error = _decode_project_updates(
        updates, project_name, program
    )
    if decode_error is not None:
        return decode_error

    resolved = _resolve_project_yaml_target(project_name, program)
    if resolved["path"] is None:
        return {
            "ok": False,
            "project_name": resolved["project_name"],
            "program": resolved["program"],
            "written_path": None,
            "error": resolved["message"]
            or "workspace project YAML is not loaded",
        }
    path = Path(str(resolved["path"]))
    before_text, document, load_error = _load_project_document(path, resolved)
    if load_error is not None:
        return load_error
    patched = _patch_project_document(document, updates, unset)

    after_text = yaml.safe_dump(patched, sort_keys=False)
    validation = validate_project_yaml(
        after_text,
        program=resolved["program"],
        project_name=resolved["project_name"],
    )
    diff = "".join(
        unified_diff(
            before_text.splitlines(keepends=True),
            after_text.splitlines(keepends=True),
            fromfile=str(path),
            tofile=str(path),
        )
    )
    if validation["verdict"] == "reject":
        return {
            "ok": False,
            "project_name": resolved["project_name"],
            "program": resolved["program"],
            "written_path": str(path),
            "validation": validation,
            "diff": diff,
            "error": "project YAML update failed validation",
        }

    path.write_text(_normalize_yaml_text(after_text), encoding="utf-8")
    state_delta = select_workspace_project(
        resolved["project_name"], resolved["program"]
    )
    return {
        "ok": True,
        "project_name": resolved["project_name"],
        "program": resolved["program"],
        "written_path": str(path),
        "validation": validation,
        "diff": diff,
        "state_delta": {"project": state_delta},
    }


def _decode_project_updates(
    updates: dict[str, Any] | str | None,
    project_name: str,
    program: str,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    if not isinstance(updates, str):
        return updates, None
    try:
        decoded = json.loads(updates)
    except json.JSONDecodeError as exc:
        return None, _project_update_error(
            project_name,
            program,
            rule_id="yaml.update.stringified_json_invalid",
            error=f"stringified updates could not be decoded: {exc}",
        )
    if not isinstance(decoded, dict):
        return None, _project_update_error(
            project_name,
            program,
            rule_id="yaml.update.mapping_required",
            error="decoded updates must be a mapping",
        )
    return decoded, None


def _project_update_error(
    project_name: str,
    program: str,
    *,
    rule_id: str,
    error: str,
) -> dict[str, Any]:
    return {
        "ok": False,
        "project_name": _normalize_project_name(project_name),
        "program": str(program or ""),
        "written_path": None,
        "rule_id": rule_id,
        "repair_hint": "Pass updates as a JSON object/dictionary.",
        "error": error,
    }


def _load_project_document(
    path: Path, resolved: dict[str, Any]
) -> tuple[str, dict[str, Any], dict[str, Any] | None]:
    before_text = path.read_text(encoding="utf-8")
    try:
        document = yaml.safe_load(before_text) or {}
    except yaml.YAMLError as exc:
        return (
            before_text,
            {},
            {
                "ok": False,
                "project_name": resolved["project_name"],
                "program": resolved["program"],
                "written_path": str(path),
                "error": f"existing YAML could not be parsed: {exc}",
            },
        )
    if not isinstance(document, dict):
        return (
            before_text,
            {},
            {
                "ok": False,
                "project_name": resolved["project_name"],
                "program": resolved["program"],
                "written_path": str(path),
                "error": "existing YAML root is not a mapping",
            },
        )
    return before_text, document, None


def _patch_project_document(
    document: dict[str, Any],
    updates: dict[str, Any] | None,
    unset: list[str] | None,
) -> dict[str, Any]:
    patched = deepcopy(document)
    for dotted, value in (updates or {}).items():
        _set_dotted_path(patched, str(dotted), value)
    for dotted in unset or []:
        _unset_dotted_path(patched, str(dotted))
    return patched


def _resolve_project_yaml_target(
    project_name: str = "",
    program: str = "",
) -> dict[str, Any]:
    requested_program = str(program or "").strip()
    requested_project = str(project_name or "").strip()
    if requested_program or requested_project:
        normalized_program = _normalize_program(
            requested_program or "gaussian"
        )
        normalized_project = _normalize_project_name(
            requested_project or "project"
        )
        path = workspace_project_path(normalized_project, normalized_program)
        if not path.exists():
            status = resolve_workspace_project()
            return {
                "project_name": normalized_project,
                "program": normalized_program,
                "path": None,
                "candidates": [str(item) for item in status.candidates],
                "message": (
                    f"Workspace project YAML not found: {path}. "
                    f"{status.message}"
                ),
            }
        return {
            "project_name": normalized_project,
            "program": normalized_program,
            "path": str(path),
            "candidates": [],
            "message": f"Loaded workspace project YAML: {path}",
        }

    selected = current_workflow_state().project
    if selected is not None:
        selected_path = Path(selected.path)
        if selected_path.is_file():
            return {
                "project_name": selected.name,
                "program": selected.program,
                "path": str(selected_path),
                "candidates": [],
                "message": f"Loaded selected workspace project YAML: {selected_path}",
            }

    status = resolve_workspace_project()
    return {
        "project_name": status.project,
        "program": status.program,
        "path": str(status.path) if status.path is not None else None,
        "candidates": [str(item) for item in status.candidates],
        "message": status.message,
    }


def _set_dotted_path(
    document: dict[str, Any],
    dotted: str,
    value: Any,
) -> None:
    parts = [part for part in dotted.split(".") if part]
    if not parts:
        raise ValueError("update path cannot be empty")
    cursor: dict[str, Any] = document
    for part in parts[:-1]:
        existing = cursor.get(part)
        if not isinstance(existing, dict):
            existing = {}
            cursor[part] = existing
        cursor = existing
    cursor[parts[-1]] = value


def _unset_dotted_path(document: dict[str, Any], dotted: str) -> None:
    parts = [part for part in dotted.split(".") if part]
    if not parts:
        return
    cursor: Any = document
    for part in parts[:-1]:
        if not isinstance(cursor, dict):
            return
        cursor = cursor.get(part)
    if isinstance(cursor, dict):
        cursor.pop(parts[-1], None)


def _load_project_yaml_via_runtime(
    *,
    yaml_text: str,
    program: str,
    project_name: str,
) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="chemsmart-project-yaml-") as tmp:
        path = Path(tmp) / f"{project_name}.yaml"
        path.write_text(_normalize_yaml_text(yaml_text), encoding="utf-8")
        if program == "gaussian":
            from chemsmart.settings.gaussian import (
                YamlGaussianProjectSettings,
            )

            settings = YamlGaussianProjectSettings.from_yaml(str(path))
        else:
            from chemsmart.settings.orca import YamlORCAProjectSettings

            settings = YamlORCAProjectSettings.from_yaml(str(path))
        return _settings_summary(settings)


def _settings_summary(settings: Any) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for name in ("opt", "ts", "sp", "td", "qmmm"):
        method = getattr(settings, f"{name}_settings", None)
        if not callable(method):
            continue
        item = method()
        if item is None:
            continue
        summary[name] = {
            "functional": getattr(item, "functional", None),
            "basis": getattr(item, "basis", None),
            "freq": getattr(item, "freq", None),
            "solvent_model": getattr(item, "solvent_model", None),
            "solvent_id": getattr(item, "solvent_id", None),
            "heavy_elements": getattr(item, "heavy_elements", None),
            "heavy_elements_basis": getattr(
                item, "heavy_elements_basis", None
            ),
            "light_elements_basis": getattr(
                item, "light_elements_basis", None
            ),
        }
    return summary


def _validation_result(
    *,
    ok: bool,
    verdict: str,
    program: str,
    project_name: str,
    issues: list[dict[str, Any]],
    parsed: Any | None = None,
    runtime_summary: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "ok": ok,
        "verdict": verdict,
        "program": program,
        "project_name": project_name,
        "issues": issues,
        "parsed": parsed,
        "runtime_summary": runtime_summary or {},
    }


def _critic_summary(verdict: str, issues: list[dict[str, Any]]) -> str:
    if verdict == "ok":
        return "Project YAML matches the extracted protocol and loads in chemsmart."
    failed = ", ".join(str(issue["rule_id"]) for issue in issues[:5])
    return f"Project YAML critic verdict is {verdict}; issues: {failed}."


__all__ = [
    "critic_project_yaml",
    "extract_project_protocol",
    "read_project_yaml",
    "render_project_yaml",
    "update_project_yaml",
    "validate_project_yaml",
    "write_project_yaml",
]
