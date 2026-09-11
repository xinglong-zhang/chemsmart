"""Three timeouts recorded nothing about what they were asked to carry.

SUFFICIENCY-1's cycle 1 died on three `inter_event_timeout` failures at
510, 463 and 598 seconds. Each recorded `input_tokens: 0` -- correctly,
because the provider never billed them -- and so said nothing at all
about the request behind them. The cause was a context that had grown
past what the silence deadline tolerated, and it was invisible until
the last *successful* attempt beside them was read by hand.

A failure that says nothing about its own request cannot calibrate the
deadline that killed it. The size travels from the same public
projection the request digest is taken from, so it carries no content
and no secret -- only how large the thing was.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.loop import ToolLoopRunner
from chemsmart.agent.request_context import build_provider_attempt_receipt

pytestmark = pytest.mark.capability("rule:*")

_DIGESTS = {
    "request_context_sha256": "a" * 64,
    "provider_budget_sha256": "b" * 64,
    "request_sha256": "c" * 64,
}


def _attempt(**extra):
    return build_provider_attempt_receipt(
        attempt_id="turn-1.provider.1",
        provider="alibaba-token-plan",
        endpoint_origin="https://example.invalid/v1",
        **_DIGESTS,
        **extra,
    )


def test_a_timeout_carries_its_request_size():
    attempt = _attempt(
        status="timeout",
        latency_ms=510_000,
        nonsecret_error_class="inter_event_timeout",
        request_bytes=880_123,
    )
    assert attempt.input_tokens == 0
    assert attempt.request_bytes == 880_123
    # And it is inside the digest, so the number is evidence.
    assert attempt.receipt_sha256


def test_the_size_is_a_counter_like_any_other():
    with pytest.raises(Exception):
        _attempt(status="timeout", request_bytes=-1)


def test_every_attempt_of_a_live_loop_carries_its_size(tmp_path):
    """The loop records the size on the success too, not only the failure.

    The diagnosis that earned this field read the last successful
    attempt beside three timeouts -- 219,005 input tokens at 374 s --
    and it was exactly the attempt that recorded nothing. The test this
    replaces asserted two source strings were present in
    ``inspect.getsource(loop)`` and never ran the loop at all.
    """

    from chemsmart.agent.runtime.alibaba import (
        Qwen38MaxConfigV1,
        Qwen38MaxToolSession,
    )
    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.runtime.events import EventKind
    from tests.agent.provider_fakes import _DispatchSpyHost, _run_contracts
    from tests.agent.test_a_flaky_provider_response_is_asked_again import (
        _FINAL,
        _MALFORMED,
    )

    responses = iter((_MALFORMED, _FINAL))
    config = Qwen38MaxConfigV1()
    session = Qwen38MaxToolSession(
        transport=lambda _payload: next(responses),
        messages=[{"role": "user", "content": "Plan the workflow."}],
        config=config,
    )
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="protocol-session"
    )
    host = _DispatchSpyHost()
    envelope, request_context, network = _run_contracts(host, config)
    ToolLoopRunner(host=host, event_store=store).run(
        session=session,
        envelope=envelope,
        request_context=request_context,
        provider_budget=network,
    )

    attempts = [
        event.payload
        for event in store.read_events()
        if event.kind == EventKind.API_ATTEMPT_OBSERVED.value
    ]
    assert [a["status"] for a in attempts] == ["protocol_failed", "succeeded"]
    assert all(a["request_bytes"] > 0 for a in attempts)
