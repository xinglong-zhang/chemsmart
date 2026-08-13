"""Domain-knowledge skills for the ChemSmart agent.

A skill is a markdown document with YAML frontmatter, mirroring the convention
used by the repository's Claude-facing skills.  It carries **knowledge**:
conventions, definitions, and established facts.  It carries no accuracy target,
no error budget, and no readiness authority.

Two authority tiers are deliberately separated:

* **Tier A -- advisory.**  The ``SKILL.md`` body.  Resolved
  ``~/.chemsmart/skills/<skill_id>/SKILL.md`` first, then the packaged builtin,
  matching how :meth:`chemsmart.settings.pyscf.PySCFProjectSettings.from_project`
  resolves project YAML.  A user overlay may change what the model reads; it has
  no effect on any host decision.
* **Tier B -- deterministic.**  :mod:`chemsmart.agent.skills.rules` convention
  rules, declared in version-controlled Python only.  An overlay file can never
  introduce or alter one, so a user-editable document cannot move a scientific
  verdict.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

import yaml

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    require_identifier,
)
from chemsmart.agent.skills.rules import (
    CONVENTION_SCOPES,
    DeterministicConventionRuleV1,
    build_convention_rule,
    rules_for_scope,
)

BUILTIN_SKILL_ROOT = Path(__file__).resolve().parent
DEFAULT_OVERLAY_ROOT = Path("~/.chemsmart/skills")

_FRONTMATTER_DELIMITER = "---"


@dataclass(frozen=True)
class SkillDocumentV1:
    """One resolved advisory skill document."""

    schema_version: str
    skill_id: str
    skill_version: str
    description: str
    body: str
    origin: str
    body_sha256: str
    document_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.domain-skill-document.v1":
            raise ContractError("unsupported skill document schema")
        require_identifier(self.skill_id, "skill_id")
        if self.origin not in ("builtin", "overlay"):
            raise ContractError("skill origin must be builtin or overlay")
        if not str(self.description).strip():
            raise ContractError("skill description must not be empty")
        if not str(self.body).strip():
            raise ContractError("skill body must not be empty")
        if self.body_sha256 != canonical_sha256({"body": self.body}):
            raise ContractError("skill body digest mismatch")
        if self.document_sha256 != canonical_sha256(self._body()):
            raise ContractError("skill document digest mismatch")

    def _body(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "skill_id": self.skill_id,
            "skill_version": self.skill_version,
            "description": self.description,
            "body_sha256": self.body_sha256,
            "origin": self.origin,
        }

    def index_entry(self) -> str:
        """Return the one-line prompt index entry for progressive disclosure."""

        return f"{self.skill_id}: {self.description}"


def _split_frontmatter(text: str) -> tuple[dict[str, object], str]:
    lines = text.splitlines()
    if not lines or lines[0].strip() != _FRONTMATTER_DELIMITER:
        raise ContractError("skill document must open with YAML frontmatter")
    for index in range(1, len(lines)):
        if lines[index].strip() == _FRONTMATTER_DELIMITER:
            header = "\n".join(lines[1:index])
            body = "\n".join(lines[index + 1 :]).strip()
            parsed = yaml.safe_load(header) or {}
            if not isinstance(parsed, dict):
                raise ContractError("skill frontmatter must be a mapping")
            return parsed, body
    raise ContractError("skill frontmatter is not terminated")


def _load_document(path: Path, origin: str) -> SkillDocumentV1:
    header, body = _split_frontmatter(path.read_text(encoding="utf-8"))
    body_sha256 = canonical_sha256({"body": body})
    fields = {
        "schema_version": "chemsmart.domain-skill-document.v1",
        "skill_id": require_identifier(
            str(header.get("name", "")), "skill_id"
        ),
        "skill_version": str(header.get("version", "0.0.0")),
        "description": str(header.get("description", "")).strip(),
        "body_sha256": body_sha256,
        "origin": origin,
    }
    return SkillDocumentV1(
        **fields,
        body=body,
        document_sha256=canonical_sha256(fields),
    )


def skills_enabled() -> bool:
    """Return whether domain-knowledge skills are surfaced to the agent.

    Defaults to enabled. ``CHEMSMART_AGENT_SKILLS=0`` removes both the prompt
    index and the ``consult_domain_skill`` tool for historical clients that
    expect the smaller surface.
    """

    return os.environ.get("CHEMSMART_AGENT_SKILLS", "1").strip() not in (
        "0",
        "false",
        "no",
    )


def _overlay_root(overlay_root: str | Path | None = None) -> Path:
    if overlay_root is not None:
        return Path(overlay_root).expanduser()
    env_root = os.environ.get("CHEMSMART_SKILL_ROOT")
    if env_root:
        return Path(env_root).expanduser()
    return DEFAULT_OVERLAY_ROOT.expanduser()


def _candidate_paths(
    skill_id: str, overlay_root: str | Path | None
) -> tuple[tuple[Path, str], ...]:
    normalized = require_identifier(skill_id, "skill_id")
    return (
        (_overlay_root(overlay_root) / normalized / "SKILL.md", "overlay"),
        (BUILTIN_SKILL_ROOT / normalized / "SKILL.md", "builtin"),
    )


def resolve_skill(
    skill_id: str, *, overlay_root: str | Path | None = None
) -> SkillDocumentV1 | None:
    """Return the overlay document if present, else the builtin, else ``None``."""

    for path, origin in _candidate_paths(skill_id, overlay_root):
        if path.is_file():
            return _load_document(path, origin)
    return None


def available_skill_ids(
    *, overlay_root: str | Path | None = None
) -> tuple[str, ...]:
    """Return every skill id resolvable from the overlay or the builtins."""

    found: set[str] = set()
    for root in (_overlay_root(overlay_root), BUILTIN_SKILL_ROOT):
        if not root.is_dir():
            continue
        for child in root.iterdir():
            if (child / "SKILL.md").is_file():
                found.add(child.name)
    return tuple(sorted(found))


def resolve_skills(
    skill_ids: tuple[str, ...], *, overlay_root: str | Path | None = None
) -> tuple[SkillDocumentV1, ...]:
    """Resolve several skills, silently skipping ids with no document."""

    resolved = []
    for skill_id in skill_ids:
        document = resolve_skill(skill_id, overlay_root=overlay_root)
        if document is not None:
            resolved.append(document)
    return tuple(resolved)


__all__ = [
    "BUILTIN_SKILL_ROOT",
    "CONVENTION_SCOPES",
    "DEFAULT_OVERLAY_ROOT",
    "DeterministicConventionRuleV1",
    "SkillDocumentV1",
    "available_skill_ids",
    "build_convention_rule",
    "resolve_skill",
    "resolve_skills",
    "rules_for_scope",
    "skills_enabled",
]
