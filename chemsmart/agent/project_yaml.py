"""Project YAML builder, validator, and critic tools for the agent."""

from __future__ import annotations

import json
import re
import tempfile
from copy import deepcopy
from difflib import unified_diff
from pathlib import Path
from typing import Any, Literal

import yaml

from chemsmart.agent.harness.basis_sets import resolve_basis_name
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    select_workspace_project,
)
from chemsmart.settings.workspace_project import (
    resolve_workspace_project,
    workspace_project_path,
)

ProjectProgram = Literal["gaussian", "orca"]

_KNOWN_PROGRAMS = {"gaussian", "orca"}
_D3BJ_ALIASES = (
    "d3bj",
    "d3-bj",
    "d3 bj",
    "becke-johnson",
    "becke johnson",
)
_CREST_MARKERS = ("crest", "mtd", "metadynamics", "gfn2", "xtb")
_DEF2_BASIS_PATTERN = (
    r"def2[- ]?(?:svp|svpd|tzvp|tzvpd|tzvpp|tzvppd|qzvp|qzvpd)"
)
_SOLVENT_ALIASES = {
    "acetonitrile": "acetonitrile",
    "chloroform": "chloroform",
    "cyclohexane": "cyclohexane",
    "dichloroethane": "dichloroethane",
    "dichloromethane": "dichloromethane",
    "dmso": "dmso",
    "ethanol": "ethanol",
    "methanol": "methanol",
    "thf": "thf",
    "tetrahydrofuran": "thf",
    "toluene": "toluene",
    "water": "water",
}


def extract_project_protocol(
    text: str,
    project_name: str = "co2",
    program: ProjectProgram = "gaussian",
) -> dict[str, Any]:
    """Extract project-YAML-relevant method facts from a literature protocol."""

    normalized_program = _normalize_program(program)
    name = _normalize_project_name(project_name)
    source = text or ""
    lowered = source.lower()
    functional = _extract_functional(lowered)
    dispersion = _extract_dispersion(lowered)
    heavy_basis = _extract_heavy_basis(source)
    light_basis = _extract_light_basis(source)
    basis = _canonical_basis_for_yaml(
        heavy_basis, light_basis, _extract_first_basis(source)
    )
    solvent = _extract_solvent(lowered)
    heavy_elements = sorted(heavy_basis)
    protocol_notes = _extract_protocol_notes(source)
    freq = _mentions_frequency_confirmation(lowered)
    result = {
        "ok": True,
        "project_name": name,
        "program": normalized_program,
        "source_excerpt": source[:1200],
        "method": {
            "functional": functional,
            "dispersion": dispersion,
            "functional_route": _functional_route(functional, dispersion),
            "basis": basis,
            "heavy_elements": heavy_elements,
            "heavy_elements_basis": (
                _first_or_none(heavy_basis.values()) if heavy_basis else None
            ),
            "light_elements_basis": light_basis,
            "solvent_model": solvent.get("solvent_model"),
            "solvent_id": solvent.get("solvent_id"),
            "freq": freq,
        },
        "protocol_notes": protocol_notes,
        "unsupported_yaml_features": _unsupported_protocol_features(lowered),
    }
    td_method = _extract_td_method(source)
    if td_method is not None:
        result["td"] = td_method
    return result


def render_project_yaml(
    protocol: dict[str, Any],
    project_name: str | None = None,
    program: ProjectProgram = "gaussian",
) -> dict[str, Any]:
    """Render a chemsmart project YAML candidate from extracted protocol facts."""

    normalized_program = _normalize_program(
        str(protocol.get("program") or program)
    )
    name = _normalize_project_name(
        project_name or str(protocol.get("project_name") or "project")
    )
    method = _method_from_protocol(protocol)

    block = _render_method_block(method, normalized_program)
    solv_freq = bool(method.get("solv_freq", False))

    # The current chemsmart project loader uses `gas` for opt/ts/scan/etc. and
    # `solv` for sp jobs. A gas-phase project still needs a solv-shaped method
    # block with no solvent fields so all job settings can be constructed.
    gas_block = deepcopy(block)
    solv_block = deepcopy(block)
    solv_block["freq"] = solv_freq
    solvent_model = _string_or_none(method.get("solvent_model"))
    solvent_id = _string_or_none(method.get("solvent_id"))
    if solvent_model and solvent_id:
        solv_block["solvent_model"] = solvent_model
        solv_block["solvent_id"] = solvent_id

    document = {"gas": gas_block, "solv": solv_block}
    td_method = protocol.get("td")
    if normalized_program == "gaussian" and isinstance(td_method, dict):
        # Gaussian's TD-DFT subcommand reads a distinct top-level `td` block.
        # Preserve an explicitly extracted excited-state method instead of
        # forcing the tool loop to patch YAML by trial and error.
        document["td"] = _render_method_block(td_method, normalized_program)
    yaml_text = yaml.safe_dump(document, sort_keys=False)
    return {
        "ok": True,
        "project_name": name,
        "program": normalized_program,
        "yaml_text": yaml_text,
        "unsupported_yaml_features": protocol.get(
            "unsupported_yaml_features", []
        ),
    }


def _render_method_block(
    method: dict[str, Any],
    program: str,
) -> dict[str, Any]:
    """Normalize one gas/solv/TD method section for project YAML output."""

    functional = _string_or_none(method.get("functional_route"))
    if functional is None:
        normalized_functional, normalized_dispersion = (
            _normalize_functional_and_dispersion(
                _string_or_none(method.get("functional")),
                _string_or_none(method.get("dispersion")),
            )
        )
        functional = _functional_route(
            normalized_functional,
            normalized_dispersion,
        )
    else:
        functional = functional.lower()
    basis = (
        _normalize_basis_if_known(_string_or_none(method.get("basis")))
        or "def2svp"
    )
    block: dict[str, Any] = {
        "functional": functional,
        "basis": basis,
        "freq": bool(method.get("freq", True)),
    }

    heavy_elements = _string_list(method.get("heavy_elements"))
    heavy_basis = _normalize_basis_if_known(
        _string_or_none(method.get("heavy_elements_basis"))
    )
    light_basis = _normalize_basis_if_known(
        _string_or_none(method.get("light_elements_basis"))
    )
    if heavy_elements and heavy_basis and basis not in {"gen", "genecp"}:
        light_basis = light_basis or basis
        basis = "gen"
        block["basis"] = basis
    if basis in {"gen", "genecp"}:
        if heavy_elements:
            block["heavy_elements"] = heavy_elements
        if heavy_basis:
            block["heavy_elements_basis"] = heavy_basis
        if light_basis:
            block["light_elements_basis"] = light_basis

    if program == "orca":
        dispersion = _string_or_none(method.get("dispersion"))
        if dispersion == "d3bj":
            block["dispersion"] = "D3BJ"
        elif dispersion == "d3":
            block["dispersion"] = "D3"
        if basis == "gen":
            block["basis"] = light_basis or "def2-svp"
    return block


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
        select_workspace_project(
            resolved["project_name"], resolved["program"]
        )
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

    if isinstance(updates, str):
        try:
            decoded = json.loads(updates)
        except json.JSONDecodeError as exc:
            return {
                "ok": False,
                "project_name": _normalize_project_name(project_name),
                "program": str(program or ""),
                "written_path": None,
                "rule_id": "yaml.update.stringified_json_invalid",
                "repair_hint": "Pass updates as a JSON object/dictionary.",
                "error": f"stringified updates could not be decoded: {exc}",
            }
        if not isinstance(decoded, dict):
            return {
                "ok": False,
                "project_name": _normalize_project_name(project_name),
                "program": str(program or ""),
                "written_path": None,
                "rule_id": "yaml.update.mapping_required",
                "repair_hint": "Pass updates as a JSON object/dictionary.",
                "error": "decoded updates must be a mapping",
            }
        updates = decoded

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
    before_text = path.read_text(encoding="utf-8")
    try:
        document = yaml.safe_load(before_text) or {}
    except yaml.YAMLError as exc:
        return {
            "ok": False,
            "project_name": resolved["project_name"],
            "program": resolved["program"],
            "written_path": str(path),
            "error": f"existing YAML could not be parsed: {exc}",
        }
    if not isinstance(document, dict):
        return {
            "ok": False,
            "project_name": resolved["project_name"],
            "program": resolved["program"],
            "written_path": str(path),
            "error": "existing YAML root is not a mapping",
        }

    patched = deepcopy(document)
    for dotted, value in (updates or {}).items():
        _set_dotted_path(patched, str(dotted), value)
    for dotted in unset or []:
        _unset_dotted_path(patched, str(dotted))

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


def _extract_functional(lowered: str) -> str | None:
    matches: list[tuple[int, str]] = []
    for pattern, canonical in (
        (r"cam[- ]?b3lyp", "camb3lyp"),
        (r"m06[- ]?2x|m062x", "m062x"),
        (r"\bpbe0\b", "pbe0"),
        (r"\bb3lyp\b", "b3lyp"),
        (r"\bpbe\b", "pbe"),
    ):
        match = re.search(pattern, lowered)
        if match:
            matches.append((match.start(), canonical))
    return min(matches, default=(0, None))[1]


def _extract_td_method(source: str) -> dict[str, Any] | None:
    """Extract a separately reported Gaussian TD-DFT method when present."""

    marker = re.search(
        r"(?i)\b(?:td[- ]?dft|tddft|time[- ]dependent|td)\b",
        source,
    )
    if marker is None:
        return None
    excerpt = source[marker.start() : marker.start() + 320]
    functional = _extract_functional(excerpt.lower())
    basis = _extract_first_basis(excerpt)
    if functional is None and basis is None:
        return None
    normalized_functional, dispersion = _normalize_functional_and_dispersion(
        functional,
        _extract_dispersion(excerpt.lower()),
    )
    return {
        "functional": normalized_functional,
        "dispersion": dispersion,
        "functional_route": _functional_route(
            normalized_functional,
            dispersion,
        ),
        "basis": basis,
        "freq": True,
    }


def _extract_dispersion(lowered: str) -> str | None:
    if any(alias in lowered for alias in _D3BJ_ALIASES):
        return "d3bj"
    if "d3" in lowered:
        return "d3"
    return None


def _normalize_functional_and_dispersion(
    functional: str | None,
    dispersion: str | None,
) -> tuple[str | None, str | None]:
    if functional is None:
        return None, dispersion
    lowered = functional.lower().replace("_", "-")
    inferred_dispersion = dispersion
    if any(alias in lowered for alias in _D3BJ_ALIASES):
        inferred_dispersion = "d3bj"
    elif re.search(
        r"(?:\bgd3\b|-d3\b|\bd3\b|empiricaldispersion=gd3\b)", lowered
    ):
        # Plain Grimme D3 (zero-damping), e.g. "b3lyp-d3" or
        # "b3lyp empiricaldispersion=gd3" — distinct from D3BJ and must not
        # be silently dropped from the rendered route.
        inferred_dispersion = "d3"
    extracted = _extract_functional(lowered)
    return (
        extracted
        or lowered.replace("-d3bj", "").replace("-d3", "").strip("- "),
        inferred_dispersion,
    )


def _extract_heavy_basis(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for basis, element_blob in re.findall(
        rf"(?i)({_DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?\s+[^.;,]{{0,120}}?\bfor\s+([A-Z][a-z]?(?:\s*(?:,|and)\s*[A-Z][a-z]?)*)\s*(?:atom|atoms|element|elements)?",
        text,
    ):
        basis_norm = _normalize_basis_name(basis)
        for symbol in re.findall(r"\b[A-Z][a-z]?\b", element_blob):
            if symbol not in {"DFT", "PES", "MTD", "GC", "BS"}:
                result[symbol] = basis_norm
    return result


def _extract_light_basis(text: str) -> str | None:
    patterns = (
        rf"(?i)({_DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?(?:\s+basis\s+set)?\s+for\s+all\s+other\s+atoms",
        rf"(?i)({_DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?(?:\s+basis\s+set)?\s+for\s+light\s+atoms",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return _normalize_basis_name(match.group(1))
    # A bare single basis is the whole-molecule basis, not a mixed-basis
    # "light elements" section. Returning the first basis here would spuriously
    # mark a single-basis protocol as mixed and break rendering/critic. The
    # general basis is resolved separately via _extract_first_basis.
    return None


def _extract_first_basis(text: str) -> str | None:
    match = re.search(rf"(?i)\b({_DEF2_BASIS_PATTERN})\b", text)
    if match:
        return _normalize_basis_name(match.group(1))
    return None


def _extract_solvent(lowered: str) -> dict[str, str | None]:
    model = None
    if "smd" in lowered:
        model = "smd"
    elif "cpcm" in lowered:
        model = "cpcm"
    elif "pcm" in lowered:
        model = "pcm"
    if model is None:
        return {"solvent_model": None, "solvent_id": None}

    solvent_id = None
    direct = re.search(
        r"(?i)\b(?:smd|cpcm|pcm)\s*\(\s*([a-z][a-z0-9 -]+?)\s*\)",
        lowered,
    )
    if direct:
        solvent_id = _normalize_solvent_id(direct.group(1))
    if solvent_id is None:
        for alias, canonical in _SOLVENT_ALIASES.items():
            if re.search(rf"\b{re.escape(alias)}\b", lowered):
                solvent_id = canonical
                break
    return {"solvent_model": model, "solvent_id": solvent_id}


def _method_from_protocol(protocol: Any) -> dict[str, Any]:
    if not isinstance(protocol, dict):
        return {}
    method = protocol.get("method")
    if isinstance(method, dict):
        return method
    method_keys = {
        "basis",
        "dispersion",
        "freq",
        "functional",
        "functional_route",
        "heavy_elements",
        "heavy_elements_basis",
        "light_elements_basis",
        "solvent_model",
        "solvent_id",
    }
    if method_keys.intersection(protocol):
        return protocol
    return {}


def _canonical_basis_for_yaml(
    heavy_basis: dict[str, str],
    light_basis: str | None,
    fallback_basis: str | None = None,
) -> str | None:
    if heavy_basis and light_basis:
        heavy_values = set(heavy_basis.values())
        if len(heavy_values) == 1 and next(iter(heavy_values)) != light_basis:
            return "gen"
    return (
        light_basis or _first_or_none(heavy_basis.values()) or fallback_basis
    )


def _functional_route(
    functional: str | None, dispersion: str | None
) -> str | None:
    if functional is None:
        return None
    if dispersion == "d3bj":
        return f"{functional} empiricaldispersion=gd3bj"
    if dispersion == "d3":
        return f"{functional} empiricaldispersion=gd3"
    return functional


def _extract_protocol_notes(text: str) -> list[str]:
    lowered = text.lower()
    notes = []
    if "crest" in lowered:
        notes.append("CREST conformational sampling was reported.")
    if "gfn2" in lowered or "xtb" in lowered:
        notes.append("GFN2-xTB semiempirical sampling level was reported.")
    if "gaussian16" in lowered or "gaussian 16" in lowered:
        notes.append("Gaussian 16 was reported for DFT refinements.")
    if "frequency" in lowered:
        notes.append(
            "Frequency analysis was reported for minima/TS confirmation."
        )
    return notes


def _unsupported_protocol_features(lowered: str) -> list[str]:
    features = []
    if any(marker in lowered for marker in _CREST_MARKERS):
        features.append("CREST/GFN2-xTB conformer sampling workflow")
    if "ten of the lowest" in lowered or "lowest energy" in lowered:
        features.append(
            "selection of lowest-energy conformers for downstream DFT"
        )
    return features


def _mentions_frequency_confirmation(lowered: str) -> bool:
    return "frequency" in lowered or "imaginary" in lowered


def _static_project_yaml_issues(
    parsed: dict[str, Any],
    program: str,
) -> list[dict[str, Any]]:
    issues: list[dict[str, Any]] = []
    if "gas" not in parsed and "solv" not in parsed:
        issues.append(
            _issue(
                "yaml.project_phase_missing",
                "reject",
                "project YAML must define at least one of gas or solv.",
            )
        )
    if "gas" in parsed and "solv" not in parsed:
        issues.append(
            _issue(
                "yaml.solv_block_required_for_loader",
                "reject",
                "current chemsmart loader requires solv when gas is present so sp settings can be built.",
            )
        )
    allowed_top = {"gas", "solv", "td", "qmmm"}
    unknown = sorted(set(parsed) - allowed_top)
    for key in unknown:
        issues.append(
            _issue(
                "yaml.unknown_top_level_key",
                "warn",
                f"unknown top-level key {key!r} will not be used by the project settings loader.",
            )
        )
    for phase in ("gas", "solv", "td", "qmmm"):
        block = parsed.get(phase)
        if block is None:
            continue
        if not isinstance(block, dict):
            issues.append(
                _issue(
                    "yaml.phase_not_mapping",
                    "reject",
                    f"{phase} must be a mapping or null.",
                )
            )
        else:
            issues.extend(_basis_catalog_issues(block, phase, program))
    if program == "gaussian":
        issues.extend(_gaussian_static_issues(parsed))
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
        value = _string_or_none(block.get(key))
        if value is None or value.lower() in {"gen", "genecp"}:
            continue
        result = resolve_basis_name(value, program=program)
        if result.verdict == "ok":
            continue
        issues.append(
            _issue(
                "yaml.basis.program_unsupported"
                if result.canonical_name
                else "yaml.basis.unrecognized",
                "reject",
                f"{phase}.{key}={value!r}: {result.message}",
            )
        )
    return issues


def _gaussian_static_issues(parsed: dict[str, Any]) -> list[dict[str, Any]]:
    issues = []
    for phase in ("gas", "solv", "td"):
        block = parsed.get(phase)
        if not isinstance(block, dict):
            continue
        basis = _string_or_none(block.get("basis"))
        functional = _string_or_none(block.get("functional"))
        if functional is None:
            issues.append(
                _issue(
                    "yaml.method_missing_functional",
                    "reject",
                    f"{phase} must define functional for a generated project YAML.",
                )
            )
        if basis is None:
            issues.append(
                _issue(
                    "yaml.method_missing_basis",
                    "reject",
                    f"{phase} must define basis for a generated project YAML.",
                )
            )
        heavy_basis = _string_or_none(block.get("heavy_elements_basis"))
        light_basis = _string_or_none(block.get("light_elements_basis"))
        heavy_elements = _string_list(block.get("heavy_elements"))
        if (heavy_basis or light_basis or heavy_elements) and basis not in {
            "gen",
            "genecp",
        }:
            issues.append(
                _issue(
                    "yaml.gaussian.mixed_basis_without_gen",
                    "reject",
                    f"{phase} defines mixed-basis fields but basis is not gen/genecp.",
                )
            )
        if basis in {"gen", "genecp"} and not (heavy_basis or light_basis):
            issues.append(
                _issue(
                    "yaml.gaussian.gen_without_basis_sections",
                    "reject",
                    f"{phase} uses {basis} but does not define heavy/light basis sections.",
                )
            )
    return issues


def _protocol_alignment_issues(
    parsed: Any,
    protocol: dict[str, Any],
) -> list[dict[str, Any]]:
    if not isinstance(parsed, dict):
        return []
    method = protocol.get("method")
    if not isinstance(method, dict):
        return []
    block = parsed.get("gas") if isinstance(parsed.get("gas"), dict) else {}
    issues = []
    expected_route = _string_or_none(method.get("functional_route"))
    if expected_route and block.get("functional") != expected_route:
        issues.append(
            _issue(
                "critic.functional_mismatch",
                "reject",
                f"gas.functional should be {expected_route!r} for the reported method.",
            )
        )
    expected_basis = _string_or_none(method.get("basis"))
    if expected_basis and block.get("basis") != expected_basis:
        issues.append(
            _issue(
                "critic.basis_mismatch",
                "reject",
                f"gas.basis should be {expected_basis!r} for the reported method.",
            )
        )
    expected_heavy = _string_list(method.get("heavy_elements"))
    if (
        expected_heavy
        and _string_list(block.get("heavy_elements")) != expected_heavy
    ):
        issues.append(
            _issue(
                "critic.heavy_elements_mismatch",
                "reject",
                f"gas.heavy_elements should be {expected_heavy!r}.",
            )
        )
    for key in ("heavy_elements_basis", "light_elements_basis"):
        expected = _string_or_none(method.get(key))
        if expected and block.get(key) != expected:
            issues.append(
                _issue(
                    f"critic.{key}_mismatch",
                    "reject",
                    f"gas.{key} should be {expected!r}.",
                )
            )
    if method.get("freq") is True and block.get("freq") is not True:
        issues.append(
            _issue(
                "critic.freq_missing",
                "warn",
                "reported minima/TS confirmation requires harmonic frequency analysis; set gas.freq: true.",
            )
        )
    td_method = protocol.get("td")
    if isinstance(td_method, dict):
        td_block = parsed.get("td") if isinstance(parsed.get("td"), dict) else {}
        expected_td = _render_method_block(td_method, "gaussian")
        if not td_block:
            issues.append(
                _issue(
                    "critic.td_block_missing",
                    "reject",
                    "reported TD-DFT method requires a top-level td block.",
                )
            )
        else:
            for key in ("functional", "basis"):
                if td_block.get(key) != expected_td.get(key):
                    issues.append(
                        _issue(
                            f"critic.td_{key}_mismatch",
                            "reject",
                            (
                                f"td.{key} should be "
                                f"{expected_td.get(key)!r} for the reported TD-DFT method."
                            ),
                        )
                    )
    return issues


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


def _issue(rule_id: str, severity: str, message: str) -> dict[str, Any]:
    return {"rule_id": rule_id, "severity": severity, "message": message}


def _normalize_program(program: str) -> str:
    normalized = (program or "").strip().lower()
    if normalized not in _KNOWN_PROGRAMS:
        raise ValueError(f"unsupported project YAML program: {program!r}")
    return normalized


def _normalize_project_name(project_name: str) -> str:
    stem = Path(project_name or "project").stem.lower()
    stem = re.sub(r"[^a-z0-9_.-]+", "_", stem).strip("._-")
    return stem or "project"


def _normalize_basis_name(value: str) -> str:
    return value.lower().replace("-", "").replace(" ", "")


def _normalize_basis_if_known(value: str | None) -> str | None:
    if value is None:
        return None
    if re.fullmatch(rf"(?i){_DEF2_BASIS_PATTERN}", value.strip()):
        return _normalize_basis_name(value)
    return value.strip()


def _normalize_solvent_id(value: str) -> str:
    key = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    return _SOLVENT_ALIASES.get(key, key.replace(" ", ""))


def _normalize_yaml_text(yaml_text: str) -> str:
    text = yaml_text.rstrip() + "\n"
    return text


def _string_or_none(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _string_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        value = value.replace(",", " ").split()
    if not isinstance(value, list):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


def _first_or_none(values: Any) -> Any:
    for value in values:
        return value
    return None
