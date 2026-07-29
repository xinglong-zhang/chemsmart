"""Extract literature protocols and render project-YAML method documents."""

from __future__ import annotations

import re
from copy import deepcopy
from typing import Any, Literal

import yaml

from chemsmart.agent.project_yaml_values import (
    DEF2_BASIS_PATTERN,
    SOLVENT_ALIASES,
    first_or_none,
    normalize_basis_if_known,
    normalize_basis_name,
    normalize_program,
    normalize_project_name,
    normalize_solvent_id,
    string_list,
    string_or_none,
)
from chemsmart.utils.periodictable import PeriodicTable

_D3BJ_ALIASES = (
    "d3bj",
    "d3-bj",
    "d3 bj",
    "becke-johnson",
    "becke johnson",
)
_CREST_MARKERS = ("crest", "mtd", "metadynamics", "gfn2", "xtb")
ProjectProgram = Literal["gaussian", "orca"]
_ELEMENT_SYMBOLS = frozenset(PeriodicTable.PERIODIC_TABLE)


def extract_project_protocol(
    text: str,
    project_name: str = "co2",
    program: ProjectProgram = "gaussian",
) -> dict[str, Any]:
    """Extract project-YAML-relevant facts from a literature protocol."""

    normalized_program = normalize_program(program)
    name = normalize_project_name(project_name)
    source = text or ""
    lowered = source.lower()
    functional = _extract_functional(lowered)
    dispersion = _extract_dispersion(lowered)
    explicit_heavy_basis = _extract_named_basis(source, "heavy_elements_basis")
    explicit_light_basis = _extract_named_basis(source, "light_elements_basis")
    explicit_heavy_elements = _extract_named_elements(source)
    heavy_basis = _extract_heavy_basis(source)
    if explicit_heavy_basis:
        heavy_elements = (
            explicit_heavy_elements
            or _extract_contextual_heavy_elements(source)
        )
        heavy_basis.update(
            {element: explicit_heavy_basis for element in heavy_elements}
        )
    light_basis = explicit_light_basis or _extract_light_basis(source)
    explicit_basis = _extract_named_basis(source, "basis")
    fallback_basis = explicit_basis
    if fallback_basis is None and not (
        explicit_heavy_basis or explicit_light_basis
    ):
        fallback_basis = _extract_first_basis(source)
    basis = _canonical_basis_for_yaml(
        heavy_basis,
        light_basis,
        fallback_basis,
    )
    solvent = _extract_solvent(lowered)
    mixed_basis_intent = bool(
        explicit_heavy_basis
        or explicit_light_basis
        or explicit_heavy_elements
        or (heavy_basis and light_basis)
    )
    missing_fields = _missing_extracted_fields(
        functional=functional,
        basis=basis,
        mixed_basis_intent=mixed_basis_intent,
        heavy_elements=sorted(heavy_basis),
        heavy_basis=first_or_none(heavy_basis.values()),
        light_basis=light_basis,
        fallback_basis=fallback_basis,
    )
    result = {
        "ok": True,
        "complete": not missing_fields,
        "missing_fields": missing_fields,
        "requires_clarification": bool(missing_fields),
        "project_name": name,
        "program": normalized_program,
        "source_excerpt": source[:1200],
        "method": {
            "functional": functional,
            "dispersion": dispersion,
            "functional_route": functional_route(functional, dispersion),
            "basis": basis,
            "heavy_elements": sorted(heavy_basis),
            "heavy_elements_basis": (
                first_or_none(heavy_basis.values()) if heavy_basis else None
            ),
            "light_elements_basis": light_basis,
            "solvent_model": solvent.get("solvent_model"),
            "solvent_id": solvent.get("solvent_id"),
            "freq": _mentions_frequency_confirmation(lowered),
        },
        "protocol_notes": _extract_protocol_notes(source),
        "unsupported_yaml_features": _unsupported_protocol_features(lowered),
    }
    td_method = _extract_td_method(source)
    if td_method is not None:
        result["td"] = td_method
    return result


def render_project_document(
    protocol: dict[str, Any],
    project_name: str | None = None,
    program: ProjectProgram = "gaussian",
) -> dict[str, Any]:
    """Render an unvalidated project-YAML document and its metadata."""

    normalized_program = normalize_program(
        str(protocol.get("program") or program)
    )
    name = normalize_project_name(
        project_name or str(protocol.get("project_name") or "project")
    )
    method = method_from_protocol(protocol)
    block = render_method_block(method, normalized_program)
    gas_block = deepcopy(block)
    solv_block = deepcopy(block)
    solv_block["freq"] = bool(method.get("solv_freq", False))
    solvent_model = string_or_none(method.get("solvent_model"))
    solvent_id = string_or_none(method.get("solvent_id"))
    if solvent_model and solvent_id:
        solv_block["solvent_model"] = solvent_model
        solv_block["solvent_id"] = solvent_id

    document = {"gas": gas_block, "solv": solv_block}
    td_method = protocol.get("td")
    if normalized_program == "gaussian" and isinstance(td_method, dict):
        document["td"] = render_method_block(td_method, normalized_program)
    return {
        "project_name": name,
        "program": normalized_program,
        "yaml_text": yaml.safe_dump(document, sort_keys=False),
        "unsupported_yaml_features": protocol.get(
            "unsupported_yaml_features",
            [],
        ),
    }


def render_method_block(
    method: dict[str, Any],
    program: str,
) -> dict[str, Any]:
    """Normalize one gas, solv, or TD method section."""

    functional = _render_functional(method)
    basis = (
        normalize_basis_if_known(string_or_none(method.get("basis")))
        or "def2svp"
    )
    block: dict[str, Any] = {
        "functional": functional,
        "basis": basis,
        "freq": bool(method.get("freq", True)),
    }
    basis, light_basis = _apply_mixed_basis(block, method)
    if program == "orca":
        _apply_orca_method(block, method, basis, light_basis)
    return block


def _render_functional(method: dict[str, Any]) -> str | None:
    functional = string_or_none(method.get("functional_route"))
    if functional is not None:
        return functional.lower()
    normalized_functional, normalized_dispersion = (
        normalize_functional_and_dispersion(
            string_or_none(method.get("functional")),
            string_or_none(method.get("dispersion")),
        )
    )
    return functional_route(normalized_functional, normalized_dispersion)


def _apply_mixed_basis(
    block: dict[str, Any],
    method: dict[str, Any],
) -> tuple[str, str | None]:
    basis = str(block["basis"])

    heavy_elements = string_list(method.get("heavy_elements"))
    heavy_basis = normalize_basis_if_known(
        string_or_none(method.get("heavy_elements_basis"))
    )
    light_basis = normalize_basis_if_known(
        string_or_none(method.get("light_elements_basis"))
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
    return basis, light_basis


def _apply_orca_method(
    block: dict[str, Any],
    method: dict[str, Any],
    basis: str,
    light_basis: str | None,
) -> None:
    dispersion = string_or_none(method.get("dispersion"))
    if dispersion == "d3bj":
        block["dispersion"] = "D3BJ"
    elif dispersion == "d3":
        block["dispersion"] = "D3"
    if basis == "gen":
        block["basis"] = light_basis or "def2-svp"


def method_from_protocol(protocol: Any) -> dict[str, Any]:
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
    return protocol if method_keys.intersection(protocol) else {}


def functional_route(
    functional: str | None,
    dispersion: str | None,
) -> str | None:
    if functional is None:
        return None
    if dispersion == "d3bj":
        return f"{functional} empiricaldispersion=gd3bj"
    if dispersion == "d3":
        return f"{functional} empiricaldispersion=gd3"
    return functional


def normalize_functional_and_dispersion(
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
        r"(?:\bgd3\b|-d3\b|\bd3\b|empiricaldispersion=gd3\b)",
        lowered,
    ):
        inferred_dispersion = "d3"
    extracted = _extract_functional(lowered)
    return (
        extracted
        or lowered.replace("-d3bj", "").replace("-d3", "").strip("- "),
        inferred_dispersion,
    )


def _extract_functional(lowered: str) -> str | None:
    matches: list[tuple[int, str]] = []
    for pattern, canonical in (
        (r"\bmn15\b", "mn15"),
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
    normalized_functional, dispersion = normalize_functional_and_dispersion(
        functional,
        _extract_dispersion(excerpt.lower()),
    )
    return {
        "functional": normalized_functional,
        "dispersion": dispersion,
        "functional_route": functional_route(
            normalized_functional,
            dispersion,
        ),
        "basis": basis,
        "freq": True,
    }


def _extract_dispersion(lowered: str) -> str | None:
    if any(alias in lowered for alias in _D3BJ_ALIASES):
        return "d3bj"
    return "d3" if "d3" in lowered else None


def _extract_heavy_basis(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    pattern = (
        rf"(?i)({DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?\s+"
        rf"(?:(?!{DEF2_BASIS_PATTERN})[^.;,]){{0,120}}?\bfor\s+"
        rf"([A-Z][a-z]?\b"
        rf"(?:\s*(?:,|and)\s*[A-Z][a-z]?\b)*)\s*"
        rf"(?:atom|atoms|element|elements)?"
    )
    for basis, element_blob in re.findall(pattern, text):
        basis_norm = normalize_basis_name(basis)
        for raw_symbol in re.findall(r"(?i)\b[a-z]{1,2}\b", element_blob):
            symbol = _canonical_element_symbol(raw_symbol)
            if symbol is not None:
                result[symbol] = basis_norm
    reverse_pattern = (
        rf"(?i)\bheavy\s+(?:atom|atoms|element|elements)"
        rf"(?:\s+such\s+as)?\s+"
        rf"([a-z]{{1,2}}\b"
        rf"(?:\s*(?:,|and)\s*[a-z]{{1,2}}\b)*)"
        rf"[^.;]{{0,80}}?\b(?:use|uses|using|with)\b[^.;]{{0,32}}?"
        rf"({DEF2_BASIS_PATTERN})\b"
    )
    for element_blob, basis in re.findall(reverse_pattern, text):
        basis_norm = normalize_basis_name(basis)
        for raw_symbol in re.findall(r"(?i)\b[a-z]{1,2}\b", element_blob):
            symbol = _canonical_element_symbol(raw_symbol)
            if symbol is not None:
                result[symbol] = basis_norm
    return result


def _extract_light_basis(text: str) -> str | None:
    patterns = (
        rf"(?i)({DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?"
        rf"(?:\s+basis\s+set)?\s+for\s+all\s+other\s+atoms",
        rf"(?i)({DEF2_BASIS_PATTERN})(?:\s*\[[^\]]+\])?"
        rf"(?:\s+basis\s+set)?\s+for\s+light\s+atoms",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            return normalize_basis_name(match.group(1))
    return None


def _extract_first_basis(text: str) -> str | None:
    match = re.search(rf"(?i)\b({DEF2_BASIS_PATTERN})\b", text)
    return normalize_basis_name(match.group(1)) if match else None


def _extract_named_basis(text: str, field: str) -> str | None:
    aliases = {field, field.replace("_", " ")}
    for alias in aliases:
        match = re.search(
            rf"(?i)(?<![a-z0-9_]){re.escape(alias)}\s*[:=]\s*"
            rf"({DEF2_BASIS_PATTERN})\b",
            text,
        )
        if match:
            return normalize_basis_name(match.group(1))
    return None


def _extract_named_elements(text: str) -> list[str]:
    match = re.search(
        r"(?i)(?<![a-z0-9_])heavy(?:_elements|\s+elements)" r"\s*[:=]\s*",
        text,
    )
    if match is None:
        return []
    remainder = text[match.end() :]
    boundary = re.search(
        r"(?i)(?:[;\n]|,\s*(?:heavy_elements_basis|"
        r"light_elements_basis|functional|basis)\s*[:=])",
        remainder,
    )
    value = remainder[: boundary.start()] if boundary else remainder
    return _canonical_elements(value)


def _extract_contextual_heavy_elements(text: str) -> list[str]:
    patterns = (
        r"(?i)\bheavy\s+(?:atom|atoms|element|elements)"
        r"(?:\s+such\s+as)?\s+"
        r"([a-z]{1,2}\b(?:\s*(?:,|and)\s*[a-z]{1,2}\b)*)",
        r"(?i)\bfor\s+"
        r"([a-z]{1,2}\b(?:\s*(?:,|and)\s*[a-z]{1,2}\b)*)"
        r"\s*(?:atom|atoms|element|elements)?(?:\s|[.;,]|$)",
    )
    for pattern in patterns:
        match = re.search(pattern, text)
        if match:
            elements = _canonical_elements(match.group(1))
            if elements:
                return elements
    return []


def _canonical_elements(value: str) -> list[str]:
    result: list[str] = []
    for raw_symbol in re.findall(r"(?i)\b[a-z]{1,2}\b", value):
        symbol = _canonical_element_symbol(raw_symbol)
        if symbol is not None and symbol not in result:
            result.append(symbol)
    return result


def _canonical_element_symbol(value: str) -> str | None:
    symbol = value[:1].upper() + value[1:].lower()
    return symbol if symbol in _ELEMENT_SYMBOLS else None


def _extract_solvent(lowered: str) -> dict[str, str | None]:
    model = (
        "smd"
        if "smd" in lowered
        else (
            "cpcm"
            if "cpcm" in lowered
            else "pcm" if "pcm" in lowered else None
        )
    )
    if model is None:
        return {"solvent_model": None, "solvent_id": None}
    solvent_id = None
    direct = re.search(
        r"(?i)\b(?:smd|cpcm|pcm)\s*\(\s*([a-z][a-z0-9 -]+?)\s*\)",
        lowered,
    )
    if direct:
        solvent_id = normalize_solvent_id(direct.group(1))
    if solvent_id is None:
        solvent_id = next(
            (
                canonical
                for alias, canonical in SOLVENT_ALIASES.items()
                if re.search(rf"\b{re.escape(alias)}\b", lowered)
            ),
            None,
        )
    return {"solvent_model": model, "solvent_id": solvent_id}


def _canonical_basis_for_yaml(
    heavy_basis: dict[str, str],
    light_basis: str | None,
    fallback_basis: str | None = None,
) -> str | None:
    if heavy_basis and light_basis:
        heavy_values = set(heavy_basis.values())
        if len(heavy_values) == 1 and next(iter(heavy_values)) != light_basis:
            return "gen"
    return light_basis or first_or_none(heavy_basis.values()) or fallback_basis


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
    return bool(
        "frequency" in lowered
        or "imaginary" in lowered
        or re.search(r"\btransition[- ]?state\b", lowered)
        or re.search(r"\bts\b", lowered)
    )


def _missing_extracted_fields(
    *,
    functional: str | None,
    basis: str | None,
    mixed_basis_intent: bool,
    heavy_elements: list[str],
    heavy_basis: str | None,
    light_basis: str | None,
    fallback_basis: str | None,
) -> list[str]:
    missing: list[str] = []
    if functional is None:
        missing.append("functional")
    if mixed_basis_intent:
        if not heavy_elements:
            missing.append("heavy_elements")
        if heavy_basis is None:
            missing.append("heavy_elements_basis")
        if light_basis is None and fallback_basis is None:
            missing.append("light_elements_basis")
    elif basis is None:
        missing.append("basis")
    return missing


__all__ = [
    "extract_project_protocol",
    "functional_route",
    "method_from_protocol",
    "normalize_functional_and_dispersion",
    "render_method_block",
    "render_project_document",
]
