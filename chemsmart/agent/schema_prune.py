"""Request-aware pruning of the ChemSmart CLI schema for synthesis prompts.

The full CLI schema serializes to ~545k chars (~99% of the synthesis system
prompt) and blocks small-context providers entirely. Synthesis only ever
emits ``chemsmart run|sub <program> <jobkind> …``, so for one request the
model needs the entry group(s), one engine (or both when ambiguous), and the
job subcommands the request could plausibly mean. Everything else is dead
weight in the prompt.

:func:`prune_schema_for_request` builds a filtered copy of the schema along
those axes. It never mutates the input; unchanged subtrees are shared by
reference (schema nodes are read-only downstream). The node shape
``{name, description, options, subcommands}`` is preserved at every level so
schema walkers such as ``_validate_tokens_against_schema`` keep working.
"""

from __future__ import annotations

import re
from typing import Any

JsonDict = dict[str, Any]

# Top-level groups synthesis can target. `agent`, `config`, and `update` are
# never legal synthesis outputs and are always dropped.
_ENTRY_KEYS = ("run", "sub")

_SUBMIT_PATTERN = re.compile(
    r"\b(submit|sbatch|qsub|slurm|pbs|hpc|cluster|queue|walltime)\b|제출|클러스터",
    re.IGNORECASE,
)
_LOCAL_PATTERN = re.compile(
    r"\b(locally|local run|on my (machine|laptop|mac)|without submission)\b"
    r"|로컬",
    re.IGNORECASE,
)

_PROGRAM_PATTERNS = {
    "gaussian": re.compile(r"\bgaussian\b|\bg16\b|\bg09\b", re.IGNORECASE),
    "orca": re.compile(r"\borca\b", re.IGNORECASE),
    # Requests name xTB by the program, by the Hamiltonian they want
    # (GFN2-xTB, GFN-FF), by the method family, or by the pre-optimization
    # role it plays before an expensive DFT job. Missing a cue here would
    # prune the xtb subcommand out of the schema, leaving the model unable
    # to emit it at all.
    "xtb": re.compile(
        r"\bxtb\b|\bgfn-?ff\b|\bgfn-?[012]\b|semi-?empirical"
        r"|tight[- ]binding|pre-?opt\w*|사전\s?최적화|반경험",
        re.IGNORECASE,
    ),
}

# Non-engine programs under run/sub, kept only when the request names them.
_AUX_PROGRAM_PATTERNS = {
    # (?!\.) keeps filenames like mol.xyz from dragging in the mol program.
    "mol": re.compile(r"\bmol\b(?!\.)|\bpymol\b|visuali[sz]", re.IGNORECASE),
    "grouper": re.compile(r"\bgrouper\b|\bgroup(ing)?\b", re.IGNORECASE),
    "database": re.compile(r"\bdatabase\b", re.IGNORECASE),
    "nciplot": re.compile(r"\bnciplot\b|\.wfn\b|\.cube\b", re.IGNORECASE),
    "thermochemistry": re.compile(
        r"\bthermochem\w*\b|boltzmann", re.IGNORECASE
    ),
    "iterate": re.compile(r"\biterate\b|\biteration\b", re.IGNORECASE),
}

# Job-kind cues. A hit keeps that subcommand (when the program has it) in
# addition to the {opt, sp} safety floor. Keys cover the union of Gaussian
# and ORCA jobkinds; per-program intersection happens during pruning.
_KIND_PATTERNS = {
    "opt": re.compile(
        r"\bopt\b|optimi[sz]\w*|geometry optim|최적화", re.IGNORECASE
    ),
    "sp": re.compile(r"\bsp\b|single[- ]?point|\benergy\b", re.IGNORECASE),
    "ts": re.compile(
        r"\bts\b|transition[- ]?state|saddle point|전이상태", re.IGNORECASE
    ),
    "irc": re.compile(r"\birc\b|intrinsic reaction", re.IGNORECASE),
    # xTB spells its frequency leaf "hess"; users ask for either word.
    "hess": re.compile(
        r"\bhess\w*\b|\bfreq\w*\b|vibrational|harmonic", re.IGNORECASE
    ),
    "scan": re.compile(
        r"\bscan\b|\bpes\b|potential energy surface", re.IGNORECASE
    ),
    "modred": re.compile(
        r"\bmodred\w*\b|freez\w+|frozen|constrain\w*", re.IGNORECASE
    ),
    "td": re.compile(
        r"\btd\b|\btd-?dft\b|excited state|uv[- ]?vis|absorption spectr",
        re.IGNORECASE,
    ),
    "wbi": re.compile(
        r"\bwbi\b|wiberg|\bnbo\b|bond (index|order)", re.IGNORECASE
    ),
    "nci": re.compile(r"\bnci\b|non[- ]?covalent", re.IGNORECASE),
    "resp": re.compile(
        r"\bresp\b|esp charge|electrostatic potential", re.IGNORECASE
    ),
    "crest": re.compile(r"\bcrest\b|conformer", re.IGNORECASE),
    "dias": re.compile(
        r"\bdias\b|distortion[/ -]?interaction|activation strain",
        re.IGNORECASE,
    ),
    "traj": re.compile(r"\btraj\w*\b", re.IGNORECASE),
    # The codebase spells QRC both "quick" (job classes) and "quasi"
    # (model_command_parser); users also describe it by what it does.
    "qrc": re.compile(
        r"\bqrc\b|qu(?:ick|asi)[- ]?reaction coordinate|imaginary mode",
        re.IGNORECASE,
    ),
    "neb": re.compile(r"\bneb\b|nudged elastic", re.IGNORECASE),
    "com": re.compile(r"\bcom file\b|\.com\b", re.IGNORECASE),
    "inp": re.compile(r"\binp\b|\.inp\b", re.IGNORECASE),
    "link": re.compile(r"\blink\s?job\b|\blink\b", re.IGNORECASE),
    "userjob": re.compile(r"\buser\s?job\b", re.IGNORECASE),
}

# Always-kept jobkind floor: cheap, and covers the most common intent when
# keyword matching misses.
_KIND_FLOOR = frozenset({"opt", "sp"})

# Confusable sibling kinds the disambiguator or the model may legitimately swap
# to. When one is cued we keep the other(s) too, so the pruned schema never
# hides the corrected target: a "scan" request the user actually means as a
# frozen ``modred`` (or vice versa), or a ``dias`` request the disambiguator
# rewrites to ``wbi``. Per-program intersection in ``_prune_program`` drops any
# sibling the selected program does not have, so this is always safe.
_CONFUSABLE_SIBLINGS = {
    "scan": frozenset({"modred"}),
    "modred": frozenset({"scan"}),
    "dias": frozenset({"wbi"}),
}

_QMMM_PATTERN = re.compile(r"\bqm/?mm\b|\boniom\b", re.IGNORECASE)


def prune_schema_for_request(
    schema: JsonDict,
    request: str,
    *,
    workspace_program: str | None = None,
) -> JsonDict:
    """Return a request-relevant filtered copy of the CLI ``schema``.

    Falls back to returning ``schema`` unchanged (same object) when the tree
    does not look like the full ChemSmart CLI schema — e.g. the minimal
    ``{"subcommands": {}}`` stubs used in tests — so callers can detect a
    no-op prune by identity.
    """

    subcommands = schema.get("subcommands")
    if not isinstance(subcommands, dict):
        return schema
    entry_names = [key for key in _ENTRY_KEYS if key in subcommands]
    if not entry_names:
        return schema

    text = str(request or "")
    kept_entries = _select_entries(entry_names, text)
    programs = _select_programs(text, workspace_program)
    kinds = _select_kinds(text)
    keep_qmmm = bool(_QMMM_PATTERN.search(text))
    aux_programs = {
        name
        for name, pattern in _AUX_PROGRAM_PATTERNS.items()
        if pattern.search(text)
    }

    pruned = {
        key: value for key, value in schema.items() if key != "subcommands"
    }
    pruned["subcommands"] = {
        entry: _prune_entry(
            subcommands[entry], programs, aux_programs, kinds, keep_qmmm
        )
        for entry in kept_entries
    }
    return pruned


def schema_variant_id(schema: JsonDict) -> str:
    """Short stable label of what a (pruned) schema contains.

    Used for decision-log / training-episode metadata and A/B reporting,
    e.g. ``"run/gaussian[opt,sp]"`` or ``"full"`` for an unpruned tree.
    """

    subcommands = schema.get("subcommands")
    if not isinstance(subcommands, dict) or not subcommands:
        return "empty"
    names = sorted(subcommands)
    if any(name not in _ENTRY_KEYS for name in names):
        return "full"
    parts = []
    for entry in names:
        programs = (subcommands[entry] or {}).get("subcommands") or {}
        program_bits = []
        for program in sorted(programs):
            kinds = (programs[program] or {}).get("subcommands") or {}
            program_bits.append(f"{program}[{','.join(sorted(kinds))}]")
        parts.append(f"{entry}/" + "+".join(program_bits))
    return " ".join(parts)


def _select_entries(entry_names: list[str], text: str) -> list[str]:
    wants_submit = bool(_SUBMIT_PATTERN.search(text))
    wants_local = bool(_LOCAL_PATTERN.search(text))
    if wants_submit and not wants_local and "sub" in entry_names:
        return ["sub"]
    if not wants_submit and "run" in entry_names:
        return ["run"]
    return list(entry_names)


def _select_programs(text: str, workspace_program: str | None) -> list[str]:
    named = [
        name
        for name, pattern in _PROGRAM_PATTERNS.items()
        if pattern.search(text)
    ]
    if named:
        return named
    workspace = (workspace_program or "").strip().lower()
    if workspace in _PROGRAM_PATTERNS:
        return [workspace]
    return list(_PROGRAM_PATTERNS)


def _select_kinds(text: str) -> set[str]:
    matched = {
        kind
        for kind, pattern in _KIND_PATTERNS.items()
        if pattern.search(text)
    }
    if not matched:
        # No cue means the model should ask or choose only the conservative
        # common floor.  Exposing every kind recreates the full-schema prompt
        # and encourages unsupported guesses.
        return set(_KIND_FLOOR)
    expanded = set(matched)
    for kind in matched:
        expanded |= _CONFUSABLE_SIBLINGS.get(kind, frozenset())
    return expanded | set(_KIND_FLOOR)


def _prune_entry(
    entry_node: JsonDict,
    programs: list[str],
    aux_programs: set[str],
    kinds: set[str],
    keep_qmmm: bool,
) -> JsonDict:
    out = {
        key: value for key, value in entry_node.items() if key != "subcommands"
    }
    available = entry_node.get("subcommands") or {}
    kept: dict[str, JsonDict] = {}
    for program in programs:
        if program in available:
            kept[program] = _prune_program(
                available[program], kinds, keep_qmmm
            )
    for aux in sorted(aux_programs):
        if aux in available and aux not in kept:
            kept[aux] = available[aux]
    out["subcommands"] = kept
    return out


def _prune_program(
    program_node: JsonDict, kinds: set[str], keep_qmmm: bool
) -> JsonDict:
    out = {
        key: value
        for key, value in program_node.items()
        if key != "subcommands"
    }
    available = program_node.get("subcommands") or {}
    if kinds:
        selected = {
            name: node for name, node in available.items() if name in kinds
        }
        if not selected:
            selected = dict(available)
    else:
        selected = dict(available)
    out["subcommands"] = {
        name: _prune_kind(node, keep_qmmm) for name, node in selected.items()
    }
    return out


def _prune_kind(kind_node: JsonDict, keep_qmmm: bool) -> JsonDict:
    children = kind_node.get("subcommands") or {}
    if keep_qmmm or "qmmm" not in children:
        return kind_node
    out = {
        key: value for key, value in kind_node.items() if key != "subcommands"
    }
    out["subcommands"] = {
        name: node for name, node in children.items() if name != "qmmm"
    }
    return out
