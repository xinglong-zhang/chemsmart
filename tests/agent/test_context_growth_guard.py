"""A conversation that outgrows the context must be refused, not paid for.

Measured across fourteen live sessions: input tokens rose on *every single
turn* -- 51 of 51, 78 of 78 -- from about 16k at the start to 146k and 249k at
the end.  Nothing is ever pruned from the session history, so total cost scales
with the square of the turn count, and a long enough workflow reaches the
provider's context limit by construction rather than by accident.

A budget check already existed, but it ran on the provider *receipt*: the
request had been sent and charged before anyone could notice it was too large.
These tests pin the pre-request guard, which spends nothing and says why.

The underlying fix -- selecting which preceding work to carry forward instead
of carrying all of it -- is not here.  chemsmart.agent.dependency_context
already implements that policy with a budget, and the live loop does not call
it; see the archived observation record.
"""

from chemsmart.agent.loop import estimate_request_input_tokens


def _request(messages, tools=()):
    return {"messages": list(messages), "tools": list(tools)}


def test_the_estimate_grows_with_the_conversation():
    short = _request([{"role": "user", "content": "hello"}])
    long = _request(
        [{"role": "user", "content": "hello"}]
        + [{"role": "tool", "content": "x" * 10_000} for _ in range(20)]
    )
    assert estimate_request_input_tokens(long) > 50 * (
        estimate_request_input_tokens(short) + 1
    )


def test_the_estimate_counts_the_tool_schema_too():
    """The schema is resent every turn and is not small."""

    bare = _request([{"role": "user", "content": "hi"}])
    with_tools = _request(
        [{"role": "user", "content": "hi"}],
        tools=[{"function": {"x": "y" * 5000}}],
    )
    assert (
        estimate_request_input_tokens(with_tools)
        > estimate_request_input_tokens(bare) + 1000
    )


def test_the_estimate_is_conservative_rather_than_optimistic():
    """Overestimating refuses a working request; underestimating pays for a
    doomed one.  The ratio must not be so high that ordinary text is refused.
    """

    text = "The optimized geometry has three equivalent N-H bonds. " * 200
    estimate = estimate_request_input_tokens(
        _request([{"role": "assistant", "content": text}])
    )
    # English prose is roughly four characters per token; the estimate must sit
    # near that, not wildly above it.
    assert len(text) / 5.0 < estimate < len(text) / 2.0


def test_a_non_string_field_is_still_counted():
    """Tool-call payloads arrive as structures, and they are not free."""

    structured = _request(
        [{"role": "assistant", "tool_calls": [{"id": "a" * 4000}]}]
    )
    assert estimate_request_input_tokens(structured) > 500


def test_a_malformed_message_does_not_crash_the_guard():
    """A guard that raises is worse than no guard."""

    assert estimate_request_input_tokens(_request(["not a mapping"])) >= 0
    assert estimate_request_input_tokens({}) == 0


def test_the_growth_the_guard_exists_for_is_reproduced_here():
    """A twenty-turn conversation that never prunes, priced turn by turn."""

    history = [{"role": "system", "content": "s" * 4000}]
    sizes = []
    for _ in range(20):
        history.append({"role": "assistant", "content": "a" * 2000})
        history.append({"role": "tool", "content": "t" * 8000})
        sizes.append(estimate_request_input_tokens(_request(history)))
    assert all(b > a for a, b in zip(sizes, sizes[1:])), "must grow every turn"
    # Total spend over the session is superlinear in the number of turns.
    assert sum(sizes) > 8 * sizes[-1]
