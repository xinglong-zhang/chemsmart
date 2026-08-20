"""Some waits a session can end, and some it cannot. It should be told which.

A consumer whose geometry comes from an `opt` or `ts` producer is deferrable:
that stage ends at one stationary structure, so the consumer can sit inside the
same approval and take the geometry when it exists.

A consumer whose geometry comes from a relaxed `scan` cannot be. A scan ends at
a surface, not a structure, and which point to carry forward is a scientific
judgement the surface has to inform. The host must not settle that on the
scientist's behalf by silently picking, say, the lowest point.

What the host did instead was worse than either: it told the session to
"materialize the declared workflow inputs", which is impossible here, and left
the node blocking approval for ever with no reason given. Every observation of
the cycle-038 paper task died this way, because the paper's own protocol is scan
then reoptimise:

    cal-scan-psi [scan] --geometry_xyz--> cal-conf-psimin [opt]

so `blocks_approval = not previewed and not deferred and not non_executable`
stayed true no matter what the session did.
"""

from __future__ import annotations

from chemsmart.agent.execution import DEFERRABLE_GEOMETRY_PRODUCER_STAGES
from chemsmart.agent.tool_runtime import _undeferrable_producer_finding


def _waiting(stage, deferrable, node_id="producer"):
    return {
        "binding_id": "filename",
        "producer_node_id": node_id,
        "producer_output_id": "out",
        "producer_stage": stage,
        "deferrable_within_one_approval": deferrable,
    }


def test_a_scan_is_not_a_deferrable_geometry_producer():
    """The premise. A scan ends at a surface, so it is deliberately absent."""

    assert "scan" not in DEFERRABLE_GEOMETRY_PRODUCER_STAGES
    assert {"opt", "ts"} <= DEFERRABLE_GEOMETRY_PRODUCER_STAGES


def test_waiting_on_an_optimisation_keeps_the_ordinary_advice():
    finding = _undeferrable_producer_finding([_waiting("opt", True)])

    assert finding["next_action"] == "materialize the declared workflow inputs"
    assert "finding" not in finding


def test_waiting_on_a_scan_names_the_producer_and_the_reason():
    finding = _undeferrable_producer_finding(
        [_waiting("scan", False, node_id="cal-scan-psi")]
    )

    assert "cal-scan-psi" in finding["finding"]
    assert "scan" in finding["finding"]
    assert "scientific choice" in finding["finding"]


def test_the_advice_offers_the_two_routes_that_actually_exist():
    """Retain it as declared intent, or plan it once a structure is chosen."""

    action = _undeferrable_producer_finding([_waiting("scan", False)])[
        "next_action"
    ]

    assert "non-executable intent" in action
    assert "chosen" in action


def test_one_undeferrable_producer_is_enough_to_change_the_advice():
    """A mixed wait must not be reported as though it were ordinary."""

    finding = _undeferrable_producer_finding(
        [_waiting("opt", True, "a"), _waiting("scan", False, "b")]
    )

    assert "finding" in finding
    assert "b (scan)" in finding["finding"]
    assert "a (" not in finding["finding"]
