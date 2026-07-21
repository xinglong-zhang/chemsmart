"""Tool-grounded, model-composed answers for command questions.

The deterministic command parser and the runtime semantic gate are treated as
*tools*: they emit ground-truth facts about a synthesized ``chemsmart`` command.
The frontier provider is then given freedom to compose the actual user-facing
answer in the user's own language and register, constrained by those facts.

A grounding harness verifies that the composed answer did not contradict the
ground-truth facts; on any contradiction, provider error, or malformed output we
fall back to the deterministic explanation. This keeps agent behaviour flexible
(natural, multilingual, question-appropriate) while remaining faithful to what
the command actually does.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any

from chemsmart.agent.model_command_parser import (
    ParsedModelCommand,
    format_parsed_model_command,
)
from chemsmart.agent.provider_adapter import extract_response_text
from chemsmart.settings.workspace_project import PROJECT_PROGRAMS

JsonDict = dict[str, Any]

# Ground-truth fields the model must not contradict. If the model echoes any of
# these in ``facts_used`` with a value that disagrees with the parsed command,
# the composed answer is rejected as hallucinated.
_HARD_FACT_KEYS = (
    "action",
    "program",
    "job",
    "project",
    "filename",
    "charge",
    "multiplicity",
)

# Fields worth surfacing to the model as ground truth (compact prompt).
_FACT_KEYS = (
    "action",
    "program",
    "job",
    "server",
    "dry_run",
    "project",
    "filename",
    "label",
    "charge",
    "multiplicity",
    "functional",
    "ab_initio",
    "basis",
    "aux_basis",
    "extrapolation_basis",
    "defgrid",
    "scf_tol",
    "scf_algorithm",
    "solvent_model",
    "solvent_id",
    "route_parameters",
    "opt_options",
    "record_index",
    "record_id",
    "structure_index",
    "structure_id",
    "molecule_id",
)

_ANSWER_SYSTEM_PROMPT = """You are the ChemSmart command assistant. You answer a user's question about ONE already-synthesized `chemsmart` command.

You are given GROUND_TRUTH: a JSON object of deterministically parsed facts about the command, plus an optional runtime semantic verdict. GROUND_TRUTH is the ONLY source of truth. You must NOT invent programs, job types, files, projects, methods, flags, charges, or multiplicities that are not in GROUND_TRUTH. If a fact is null/absent, say it is runtime-derived or unspecified rather than guessing.

Write for the user, not for a machine. Respond in the SAME language the user wrote their request in, and honour any explicit language request (e.g. a request that says "in Chinese"/"중국어"/"in Korean" must be answered in that language). Be concise and natural — a short paragraph, not a field dump. Match the question:
- explain_command: say what the command does and why, in plain language.
- critique_command: judge whether it is valid/safe/correctly grounded, citing the semantic verdict.
- repair_command: state whether repair is needed and, if so, what must change (only using legal facts/flags).

Return ONLY one JSON object:
{
  "answer": "<user-facing prose in the requested language>",
  "facts_used": {"program": "...", "job": "...", "action": "...", "project": "...", "filename": "...", "charge": "...", "multiplicity": "..."},
  "reasoning": ["short observable reasoning step", ...],
  "caveats": ["short caveat if any", ...]
}

`facts_used` MUST echo the corresponding GROUND_TRUTH values verbatim (use the same string, or null when GROUND_TRUTH is null). `reasoning` is a public, user-auditable rationale (NOT hidden chain-of-thought): 1-4 short observable steps. Do not include any text outside the JSON object."""

_MISSING_INFO_SYSTEM_PROMPT = """You are the ChemSmart intake reasoner. A user asked to prepare/run a computational chemistry job, but the request may be missing information required to synthesize a concrete `chemsmart` command.

Reason about what is PRESENT and what is MISSING, then produce the minimal set of clarifying questions needed to proceed. Only ask for information that is genuinely required and not already implied by the request. Prefer 1-3 sharp questions over a long checklist. Respond to the user in the SAME language they wrote in (honour explicit language requests).

Return ONLY one JSON object:
{
  "present": ["fact the request already provides", ...],
  "missing": ["required fact that is absent", ...],
  "questions": ["clarifying question in the user's language", ...],
  "reasoning": ["short observable reasoning step", ...]
}

`reasoning` is a public, user-auditable rationale (NOT hidden chain-of-thought). Do not include any text outside the JSON object."""


@dataclass(frozen=True)
class ComposedAnswer:
    """Result of composing a grounded answer to a command question."""

    answer: str
    reasoning: tuple[str, ...] = ()
    caveats: tuple[str, ...] = ()
    grounded: bool = True
    fallback_used: bool = False
    fallback_reason: str = ""
    facts_used: JsonDict = field(default_factory=dict)


@dataclass(frozen=True)
class MissingInfoReasoning:
    """Chain-of-thought result for judging insufficient job information."""

    questions: tuple[str, ...] = ()
    reasoning: tuple[str, ...] = ()
    present: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()


def build_command_facts(
    parsed: ParsedModelCommand,
    *,
    semantic_summary: JsonDict | None = None,
) -> JsonDict:
    """Return a compact GROUND_TRUTH dict for one parsed command."""

    data = parsed.to_dict()
    facts: JsonDict = {"command": parsed.command}
    for key in _FACT_KEYS:
        value = data.get(key)
        if value not in (None, "", {}, []):
            facts[key] = value
    if parsed.parse_error:
        facts["parse_error"] = parsed.parse_error
    structural = data.get("structural_options") or {}
    if structural:
        facts["job_specific_options"] = structural
    if parsed.warnings:
        facts["parser_warnings"] = list(parsed.warnings)
    if semantic_summary:
        facts["semantic_verdict"] = semantic_summary.get("verdict")
        issues = semantic_summary.get("issues") or []
        if issues:
            facts["semantic_issues"] = [
                {
                    "rule_id": issue.get("rule_id"),
                    "message": issue.get("message"),
                }
                for issue in issues
                if isinstance(issue, dict)
            ]
        missing = semantic_summary.get("missing_info") or []
        if missing:
            facts["semantic_missing_info"] = list(missing)
    return facts


def compose_command_answer(
    provider: Any,
    *,
    request: str,
    action: str,
    parsed: ParsedModelCommand,
    semantic_summary: JsonDict | None = None,
    timeout_s: float = 30.0,
) -> ComposedAnswer:
    """Compose a grounded, language-appropriate answer to a command question.

    Falls back to the deterministic explanation on provider error, malformed
    output, or a grounding violation.
    """

    deterministic = format_parsed_model_command(parsed)
    facts = build_command_facts(parsed, semantic_summary=semantic_summary)
    messages = [
        {"role": "system", "content": _ANSWER_SYSTEM_PROMPT},
        {
            "role": "user",
            "content": json.dumps(
                {
                    "action": action,
                    "request": request,
                    "ground_truth": facts,
                },
                ensure_ascii=False,
            ),
        },
    ]
    try:
        raw = provider.chat(messages, timeout_s=timeout_s)
    except TypeError:
        # Providers whose chat() does not accept timeout_s.
        try:
            raw = provider.chat(messages)
        except Exception as exc:  # pragma: no cover - defensive
            return _fallback_answer(deterministic, f"provider error: {exc}")
    except Exception as exc:
        return _fallback_answer(deterministic, f"provider error: {exc}")

    try:
        parsed_response = _load_json_object(extract_response_text(raw))
    except (ValueError, TypeError, json.JSONDecodeError) as exc:
        return _fallback_answer(deterministic, f"malformed response: {exc}")

    answer = str(parsed_response.get("answer") or "").strip()
    if not answer:
        return _fallback_answer(deterministic, "empty answer")

    facts_used = parsed_response.get("facts_used")
    if not isinstance(facts_used, dict):
        facts_used = {}
    contradiction = _grounding_contradiction(facts_used, facts, answer)
    if contradiction:
        return _fallback_answer(deterministic, contradiction)

    return ComposedAnswer(
        answer=answer,
        reasoning=_string_tuple(parsed_response.get("reasoning")),
        caveats=_string_tuple(parsed_response.get("caveats")),
        grounded=True,
        fallback_used=False,
        facts_used={k: v for k, v in facts_used.items()},
    )


def reason_missing_info(
    provider: Any,
    *,
    request: str,
    hints: tuple[str, ...] = (),
    timeout_s: float = 30.0,
) -> MissingInfoReasoning | None:
    """Run the missing-information chain-of-thought; None on failure."""

    messages = [
        {"role": "system", "content": _MISSING_INFO_SYSTEM_PROMPT},
        {
            "role": "user",
            "content": json.dumps(
                {"request": request, "known_hints": list(hints)},
                ensure_ascii=False,
            ),
        },
    ]
    try:
        raw = provider.chat(messages, timeout_s=timeout_s)
    except TypeError:
        try:
            raw = provider.chat(messages)
        except Exception:  # pragma: no cover - defensive
            return None
    except Exception:
        return None

    try:
        parsed_response = _load_json_object(extract_response_text(raw))
    except (ValueError, TypeError, json.JSONDecodeError):
        return None

    questions = _string_tuple(parsed_response.get("questions"))
    reasoning = _string_tuple(parsed_response.get("reasoning"))
    present = _string_tuple(parsed_response.get("present"))
    missing = _string_tuple(parsed_response.get("missing"))
    if not questions and not missing:
        return None
    return MissingInfoReasoning(
        questions=questions,
        reasoning=reasoning,
        present=present,
        missing=missing,
    )


def _fallback_answer(deterministic: str, reason: str) -> ComposedAnswer:
    return ComposedAnswer(
        answer=deterministic,
        reasoning=(),
        caveats=(),
        grounded=False,
        fallback_used=True,
        fallback_reason=reason,
    )


def _grounding_contradiction(
    facts_used: JsonDict, facts: JsonDict, answer: str
) -> str:
    """Return a non-empty reason string when the answer contradicts facts."""

    # Primary check: any echoed hard fact must match ground truth.
    for key in _HARD_FACT_KEYS:
        if key not in facts_used:
            continue
        echoed = facts_used.get(key)
        truth = facts.get(key)
        if not _facts_equal(echoed, truth):
            return (
                f"grounding violation: facts_used.{key}={echoed!r} "
                f"contradicts ground truth {truth!r}"
            )

    # Secondary check: proper nouns survive translation, so a wrong program
    # name appearing in prose (while the correct one is absent) is a red flag.
    program = str(facts.get("program") or "").lower()
    other_programs = set(PROJECT_PROGRAMS) - {program}
    lowered_answer = answer.lower()
    if program and any(name in lowered_answer for name in other_programs):
        if program not in lowered_answer:
            return (
                "grounding violation: answer references a different program "
                f"than ground truth {program!r}"
            )
    return ""


def _facts_equal(echoed: Any, truth: Any) -> bool:
    if echoed is None and truth is None:
        return True
    if echoed is None or truth is None:
        return False
    return str(echoed).strip().lower() == str(truth).strip().lower()


def _string_tuple(value: Any) -> tuple[str, ...]:
    if not isinstance(value, list):
        return ()
    return tuple(str(item).strip() for item in value if str(item).strip())


def _load_json_object(text: str) -> JsonDict:
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    parsed = json.loads(text)
    if not isinstance(parsed, dict):
        raise ValueError("expected a JSON object")
    return parsed
