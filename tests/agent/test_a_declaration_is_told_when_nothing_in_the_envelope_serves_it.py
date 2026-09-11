"""A declaration is told, at once, when no envelope program serves it.

NOVEL-3 ino3 (2026-09-05): a session declared the cation's spin
distribution as an observable, ran four cycles, and learned only at the
end that no route to it existed inside the envelope. The readers' own
selector vocabulary says what a meaning names; the envelope says which
programs may run; the declaration reply says when the two do not meet,
and names the two routes. Warn and route, never refuse.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest


def _host(tmp_path, programs):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    host = _host(tmp_path, approved_requested_observable_declarations=[])
    host.bounded_execution_envelope = SimpleNamespace(
        allowed_program_engines=tuple(
            (program, ("cpu",)) for program in programs
        ),
        resources=SimpleNamespace(resource_sha256="a" * 64),
    )
    return host


def _declare(host, meaning, observable_id="spin-ni"):
    return host.dispatch(
        turn_id="t1",
        tool_name="declare_requested_observable",
        arguments={
            "observables": [
                {
                    "observable_id": observable_id,
                    "unit": "1",
                    "meaning": meaning,
                }
            ]
        },
    )["result"]


@pytest.mark.capability("tool:declare_requested_observable")
def test_a_kind_no_envelope_program_declares_is_named_with_its_routes(
    tmp_path,
):
    host = _host(tmp_path, ("xtb",))
    reply = _declare(
        host, "Mulliken spin population on the nickel atom of the cation"
    )
    assert reply["declared"][0]["observable_id"] == "spin-ni"
    warning = reply["reach_warnings"]["spin-ni"]
    assert "mulliken_atomic_spin_populations" in warning
    assert "no envelope program (xtb) declares" in warning
    assert "orca" in warning
    assert "unreachable_observable_ids" in warning
    assert "Nothing is refused here." in warning


def test_a_kind_the_envelope_serves_or_none_it_can_judge_says_nothing(
    tmp_path,
):
    served = _declare(
        _host(tmp_path / "orca", ("orca",)),
        "Mulliken spin population on the nickel atom of the cation",
    )
    assert "reach_warnings" not in served
    unjudged = _declare(
        _host(tmp_path / "xtb", ("xtb",)),
        "difference of two Gibbs energies at 298 K",
        observable_id="ddg",
    )
    assert "reach_warnings" not in unjudged
