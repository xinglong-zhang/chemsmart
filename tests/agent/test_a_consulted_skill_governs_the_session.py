"""Fetching knowledge and being governed by it are not the same event.

Anthropic's account of the pattern is that an invoked skill's body "becomes
loaded context, not inert data". ChemSmart returned it stamped
`advisory_only`, `readiness_authority: False`, `accuracy_authority: False` and
nothing else. Those flags exist for a real reason and must not move: a skill
may never establish readiness, approval, terminal state or an accuracy claim.

But they were carrying a second meaning they were never meant to carry. "This
establishes no scientific status" and "you need not follow this" are different
statements, and a payload that says only the first is read as both. A session
consults a skill in order to be governed by it; the payload should say so.

So this check pins both axes at once, because the failure mode is fixing one by
weakening the other: adoption language must appear, and every status flag must
still read exactly as before.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.skills import resolve_skill

_SKILL = "method-adequacy"


@pytest.fixture()
def host(tmp_path):
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="s1"
    )
    return CommandCompiledToolHostV1(event_store=store)


def _consult(host, skill_id=_SKILL):
    """Through `dispatch`, which is what a turn actually calls.

    The handler's return value is wrapped in a status envelope before it
    reaches the model, so a check that called the private method directly
    would assert against a shape no session ever sees.
    """

    envelope = host.dispatch(
        turn_id="t1",
        tool_name="consult_domain_skill",
        arguments={"skill_id": skill_id},
    )
    assert envelope["tool"] == "consult_domain_skill"
    assert envelope["status"] == "ok"
    return envelope["result"]


def test_the_status_flags_are_exactly_what_they_were(host):
    """The governance guarantee. If this loosens, the fix went wrong."""

    payload = _consult(host)

    assert payload["advisory_only"] is True
    assert payload["readiness_authority"] is False
    assert payload["accuracy_authority"] is False


def test_the_payload_says_the_session_now_works_under_it(host):
    payload = _consult(host)

    assert payload["applies_to"] == "the rest of this session"
    assert "work under it" in payload["guidance"]


def test_adoption_never_outranks_a_receipt(host):
    """Guidance the session follows must still lose to typed evidence."""

    guidance = _consult(host)["guidance"]

    assert "receipt" in guidance
    for status in ("readiness", "approval", "terminal state", "accuracy"):
        assert status in guidance


def test_the_body_is_delivered_verbatim_and_digest_bound(host):
    """Replay must reconstruct the exact text the model was given."""

    document = resolve_skill(_SKILL)
    payload = _consult(host)

    assert payload["body"] == document.body
    assert payload["body_sha256"] == document.body_sha256
    assert payload["document_sha256"] == document.document_sha256


def test_consultation_is_recorded_against_the_session(host):
    document = resolve_skill(_SKILL)

    _consult(host)

    assert document.document_sha256 in host.consulted_skills


def test_an_unknown_skill_is_refused_rather_than_improvised(host):
    with pytest.raises(ContractError, match="unknown domain skill"):
        _consult(host, skill_id="method-adequacy-v2")
