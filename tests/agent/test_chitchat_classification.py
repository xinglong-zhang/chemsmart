from __future__ import annotations

import pytest

from chemsmart.agent.core import (
    Plan,
    _is_chitchat_request,
    _resolve_plan_intent,
)


@pytest.mark.parametrize(
    "user_request",
    [
        "hi, what is your name?",
        "what is your name",
        "who are you",
        "what are you",
        "introduce yourself",
        "tell me about yourself",
        "which model are you",
        "hi",
        "thanks",
    ],
)
def test_is_chitchat_request_accepts_identity_and_smalltalk(
    user_request: str,
):
    assert _is_chitchat_request(user_request) is True


@pytest.mark.parametrize(
    "user_request",
    [
        "optimize h2o.xyz",
        "transition state for carbene",
    ],
)
def test_is_chitchat_request_rejects_chemistry_requests(
    user_request: str,
):
    assert _is_chitchat_request(user_request) is False


@pytest.mark.parametrize(
    "user_request",
    [
        "hi, what is your name?",
        "what is your name",
        "who are you",
        "what are you",
        "introduce yourself",
        "tell me about yourself",
        "which model are you",
        "hi",
        "thanks",
    ],
)
def test_resolve_plan_intent_marks_zero_step_smalltalk_as_chitchat(
    user_request: str,
):
    plan = Plan(
        steps=[],
        rationale="I'm the chemsmart agent.",
        estimated_cost="none",
    )

    assert _resolve_plan_intent(user_request, plan) == "chitchat"
