"""The one vocabulary authority for every human-facing string.

The terminal speaks to a chemist. Machine state tokens, protocol words, and
Python class names stay in the durable records; the display translates them
once, here, so no surface invents its own dialect. Unknown values pass
through verbatim -- host-authored sentences are already written for humans.

This module imports no UI framework so the CLI may use it too.
"""

from __future__ import annotations

#: Machine state -> what a chemist reads. Identity entries are listed so a
#: grep for a state token finds its display form here first.
HUMAN_WORDS: dict[str, str] = {
    # Session terminal states.
    "waiting_for_approval": "ready for your review",
    "planned": "planned",
    "complete": "complete",
    "failed": "failed",
    "blocked": "blocked",
    "cancelled": "cancelled",
    # TUI phases.
    "ready": "ready",
    "planning": "planning",
    "preview-ready": "preview complete",
    "request-reviewed": "awaiting your approval",
    "executing": "executing",
    "error": "error",
    # Workflow node run states.
    "pending": "queued",
    "queued": "queued",
    "deferred": "deferred",
    "running": "running",
    "engine_complete": "validating result",
    "validated": "validated",
    "ambiguous": "ambiguous",
    # Analysis node states.
    "executed": "executed",
    "skipped": "skipped",
    "blocked_unsupported": "not executable in this release",
    "blocked_upstream": "waiting on a blocked stage",
    "waiting_for_artifact": "waiting for its input",
    "actionable": "ready",
    # Analysis intent support states.
    "resolvable": "planned",
    "unresolved_future": "planned",
}

#: Host refusal reasons that read as protocol, mapped to plain speech.
TURN_BLOCKED_REASONS: dict[str, str] = {
    "provider token budget exceeded": (
        "the model's token budget for this session is used up"
    ),
    "public provider tool envelope failed atomic decoding": (
        "the model's reply could not be decoded safely"
    ),
}

#: Composition evidence statuses, spoken.
IDENTITY_EVIDENCE_WORDS: dict[str, str] = {
    "composed-from-approved-parents": (
        "built from two approved parent structures"
    ),
    "composed-task-bound": (
        "built from two task-bound parent structures"
    ),
}


def human_state(value: str) -> str:
    return HUMAN_WORDS.get(str(value), str(value))


def human_turn_blocked(reason: str) -> str:
    return TURN_BLOCKED_REASONS.get(str(reason), str(reason))


def human_identity_evidence(status: str) -> str:
    return IDENTITY_EVIDENCE_WORDS.get(str(status), str(status))


__all__ = [
    "HUMAN_WORDS",
    "IDENTITY_EVIDENCE_WORDS",
    "TURN_BLOCKED_REASONS",
    "human_identity_evidence",
    "human_state",
    "human_turn_blocked",
]
