"""SPEC-level invariants shared by runtime eval and dataset gates."""

from __future__ import annotations

import re
from typing import Any

from chemsmart.agent import v8_kind_index as KI
from chemsmart.agent.harness.invariants.gaussian_ts import _ALLOWED_EXTRA
from chemsmart.agent.harness.models import InvariantIssue

CANON_KINDS = set(KI.KIND_SETTINGS)

RUNTIME_OWNED = {
    "project",
    "functional",
    "basis",
    "semiempirical",
    "ab_initio",
    "aux_basis",
    "extrapolation_basis",
    "dispersion",
    "defgrid",
    "scf_tol",
    "scf_algorithm",
    "scf_maxiter",
    "scf_convergence",
    "solvent_model",
    "solvent_id",
    "custom_solvent",
    "heavy_elements_basis",
    "light_elements_basis",
    "gen_genecp_file",
}
_RUNTIME_OWNED_PREFIX = ("mdci_",)

TS_RUNTIME_ROUTE = {"ts", "calcfc", "noeigentest"}

KIND_CANON = {
    "gaussian.opt+freq": "gaussian.opt",
    "gaussian.opt.freq": "gaussian.opt",
    "gaussian.opt_freq": "gaussian.opt",
    "orca.opt+freq": "orca.opt",
    "orca.opt.freq": "orca.opt",
    "orca.opt_freq": "orca.opt",
    "gaussian.qm/mm": "gaussian.qmmm",
    "orca.qm/mm": "orca.qmmm",
}

REQUIRED_SLOT = {
    "gaussian.scan": "scan_definition",
    "orca.scan": "scan_definition",
    "gaussian.modred": "modred",
    "orca.modred": "modred",
    "gaussian.dias": "fragment_indices",
    "gaussian.qmmm": "high_level_atoms",
    "orca.qmmm": "high_level_atoms",
}

DB_SELECTOR_KEYS = {
    "record_index",
    "record_id",
    "structure_index",
    "structure_id",
    "molecule_id",
}
DB_RECORD_SELECTOR_KEYS = {"record_index", "record_id"}
DB_TOP_LEVEL_SELECTOR_KEYS = {"record_index", "record_id", "structure_id"}

_CUE = re.compile(
    r"\b(atom|atoms|fragment|fragments|bond|dihedral|angle|torsion|distance|"
    r"length|freeze|frozen|constrain|constraint|between|indices|index|"
    r"region|high[\s-]*level|qm[\s-]*region|reaction\s+coordinate|"
    r"coordinate)\b",
    re.IGNORECASE,
)
_ATOM_PAIR = re.compile(r"\b\d+\s*[-,、]\s*\d+\b")
_COORD_SPEC = re.compile(r"\b[BAD]\s+\d+(?:\s+\d+)+\b")


def query_specifies_atoms(query: str) -> bool:
    text = (query or "").replace("_", " ")
    return bool(_CUE.search(text) or _ATOM_PAIR.search(text) or _COORD_SPEC.search(text))


def check_spec(spec: dict[str, Any], query: str) -> list[InvariantIssue]:
    if not isinstance(spec, dict):
        return [_issue("spec.parse", "spec is not an object", {})]
    intent = spec.get("intent")
    if intent != "workflow":
        if intent in {"advisory", "decline", "chitchat"} and not spec.get("message"):
            return [_issue("spec.message.missing", f"intent={intent} requires a message", {})]
        return []

    jobs = spec.get("jobs") or []
    if not jobs:
        return [_issue("spec.jobs.empty", "workflow intent with no jobs", {})]

    issues: list[InvariantIssue] = []
    for job in jobs:
        issues.extend(check_job(job, query))
    return issues


def check_job(job: dict[str, Any], query: str) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    kind = job.get("kind")
    evidence = {"kind": kind}

    if kind not in CANON_KINDS:
        canon = KIND_CANON.get(str(kind).replace(" ", ""))
        issues.append(
            _issue(
                "spec.kind.canonical",
                f"non-canonical kind {kind!r}" + (f" -> use {canon!r}" if canon else ""),
                {**evidence, "suggested": canon},
            )
        )
        return issues

    if job.get("project") not in (None, ""):
        issues.append(
            _issue(
                "spec.project.runtime_owned",
                "project selection is runtime-owned and must not be model-emitted",
                {**evidence, "project": job.get("project")},
            )
        )

    issues.extend(_db_selector_issues(job, str(kind), evidence))

    settings = job.get("settings") or {}
    leaked = sorted(key for key in settings if _is_runtime_owned(key))
    if leaked:
        issues.append(
            _issue(
                "spec.runtime_owned.leak",
                f"runtime-owned keys must not be model-emitted: {leaked}",
                {**evidence, "leaked": leaked},
            )
        )

    bad = sorted(
        key
        for key in settings
        if key not in KI.KIND_SETTINGS.get(kind, {}) and not _is_runtime_owned(key)
    )
    if bad:
        issues.append(
            _issue(
                "spec.settings.allowed",
                f"settings not allowed for {kind}: {bad}",
                {**evidence, "disallowed": bad},
            )
        )

    if isinstance(kind, str) and kind.endswith(".ts"):
        if settings.get("freq") is not None:
            issues.append(
                _issue(
                    "spec.ts.runtime_freq",
                    "freq is runtime-owned for *.ts; do not emit it",
                    evidence,
                )
            )
        extras = settings.get("additional_opt_options_in_route") or []
        if isinstance(extras, str):
            extras = re.split(r"[,\s]+", extras)
        tokens = [str(token).strip().lower() for token in extras if str(token).strip()]
        canon_in_extras = sorted(set(tokens) & TS_RUNTIME_ROUTE)
        if canon_in_extras:
            issues.append(
                _issue(
                    "spec.ts.runtime_route",
                    f"runtime-owned TS route token(s) emitted as extras: {canon_in_extras}",
                    {**evidence, "tokens": canon_in_extras},
                )
            )
        junk = sorted(
            token
            for token in tokens
            if token not in TS_RUNTIME_ROUTE and not _ALLOWED_EXTRA.match(token)
        )
        if junk:
            issues.append(
                _issue(
                    "spec.ts.bad_extra",
                    f"non-allowlisted TS opt extra(s): {junk}",
                    {**evidence, "tokens": junk},
                )
            )

    slot = REQUIRED_SLOT.get(kind)
    if slot:
        present = slot in settings and settings.get(slot) not in (None, [], "")
        if not query_specifies_atoms(query):
            issues.append(
                _issue(
                    "spec.decline_contract",
                    f"{kind} needs {slot} but the query specifies no atoms/fragments",
                    {**evidence, "slot": slot},
                )
            )
        elif not present:
            issues.append(
                _issue(
                    "spec.required_present",
                    f"{kind} query specifies atoms but {slot} is missing",
                    {**evidence, "slot": slot},
                )
            )

    label = job.get("label")
    if label and "label" not in (query or "").lower():
        issues.append(
            _issue(
                "spec.label.runtime_owned",
                "label is runtime-owned; omit unless the query asks for one",
                {**evidence, "label": label},
                severity="warn",
            )
        )
    return issues


def _db_selector_issues(
    job: dict[str, Any],
    kind: str,
    evidence: dict[str, Any],
) -> list[InvariantIssue]:
    issues: list[InvariantIssue] = []
    source = job.get("file")
    is_db = isinstance(source, str) and source.lower().endswith(".db")
    present = {
        key: job.get(key)
        for key in DB_SELECTOR_KEYS
        if job.get(key) not in (None, "")
    }
    if not is_db and not present:
        return issues

    if present.get("molecule_id") is not None:
        issues.append(
            _issue(
                "spec.db.molecule_id_job",
                "molecule_id is not supported for Gaussian/ORCA job submission",
                {**evidence, "selectors": sorted(present)},
            )
        )

    if not is_db:
        issues.append(
            _issue(
                "spec.db.selector_without_db",
                "database selectors require a .db source file",
                {**evidence, "file": source, "selectors": sorted(present)},
            )
        )
        return issues

    top_level = [key for key in DB_TOP_LEVEL_SELECTOR_KEYS if key in present]
    if len(top_level) != 1:
        issues.append(
            _issue(
                "spec.db.selector_cardinality",
                (
                    ".db job input must select exactly one of record_index, "
                    "record_id, or structure_id"
                ),
                {**evidence, "selectors": sorted(present)},
            )
        )

    if "structure_index" in present and not (
        DB_RECORD_SELECTOR_KEYS & set(present)
    ):
        issues.append(
            _issue(
                "spec.db.structure_index_requires_record",
                "structure_index requires record_index or record_id",
                {**evidence, "selectors": sorted(present)},
            )
        )
    return issues


def _is_runtime_owned(key: str) -> bool:
    return key in RUNTIME_OWNED or key.startswith(_RUNTIME_OWNED_PREFIX)


def _issue(
    rule_id: str,
    message: str,
    evidence: dict[str, Any],
    severity: str = "reject",
) -> InvariantIssue:
    return InvariantIssue(
        rule_id=rule_id,
        severity=severity,
        message=message,
        evidence=evidence,
    )


__all__ = [
    "CANON_KINDS",
    "DB_SELECTOR_KEYS",
    "KIND_CANON",
    "REQUIRED_SLOT",
    "RUNTIME_OWNED",
    "TS_RUNTIME_ROUTE",
    "check_job",
    "check_spec",
    "query_specifies_atoms",
]
