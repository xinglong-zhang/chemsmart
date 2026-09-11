"""A declared observable is delivered by id and by dimension, once.

OPEN-1 ino3 (2026-09-07) declared six spin observables in ``e`` while
ChemSmart represents spin populations as dimensionless counts. The
completion gate recorded six declared-observable limitations naming the
mismatch; the settlement read the workspace record by id alone, called
the same six delivered in an earlier cycle, and settled ``achieved``.
Two auditors reproduced it independently on the live stream. The
predicate now lives in one module both readers call, and a mistaken
unit is repairable by declaring the right one and naming what it
replaces.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent.delivery import (
    observable_is_delivered,
    superseded_observable_ids,
)
from chemsmart.agent.driver import _AnalysisDelivery
from chemsmart.agent.tool_runtime import ContractError

from .test_a_guide_opens_when_something_asks import _host

_SPIN = {
    "observable_id": "spin-ni",
    "unit": "e",
    "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
    "meaning": "spin population on nickel",
}


def _claim(unit, dimension):
    return SimpleNamespace(
        claim_id="spin-ni",
        quantity_id="spin-ni",
        dimension=dimension,
        display_value=0.787,
        display_unit=unit,
    )


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_a_charge_unit_is_not_answered_by_a_dimensionless_count():
    row = {"claim_id": "spin-ni", "unit": "1", "dimension": (0,) * 6}
    assert not observable_is_delivered(_SPIN, row)
    assert observable_is_delivered(_SPIN, {"claim_id": "spin-ni", "unit": "e"})


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_the_settlement_projection_keeps_the_gate_verdict():
    delivery = _AnalysisDelivery(
        anomaly_output_ids=(),
        unanswered_verdicts=(),
        completion_status="passed",
        limitation_output_ids=("declared_observable:spin-ni",),
        claims=1,
        decisions=1,
        receipt_sha256s=(),
        declared_observables=(_SPIN,),
        goal_delivered={
            "spin-ni": {"cycle": 1, "unit": "1", "dimension": (0,) * 6}
        },
    )
    assert delivery.undelivered_declared_ids == ("spin-ni",)
    assert delivery.delivered_in_earlier_cycles == ()

    answered = _AnalysisDelivery(
        anomaly_output_ids=(),
        unanswered_verdicts=(),
        completion_status="passed",
        limitation_output_ids=("declared_observable:spin-ni",),
        claims=1,
        decisions=1,
        receipt_sha256s=(),
        declared_observables=(_SPIN,),
        goal_delivered={
            "spin-ni": {
                "cycle": 1,
                "unit": "e",
                "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
            }
        },
    )
    assert answered.undelivered_declared_ids == ()
    assert answered.delivered_in_earlier_cycles == (
        "spin-ni (delivered in cycle 1)",
    )


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_a_mistaken_unit_is_repaired_by_declaring_its_replacement(tmp_path):
    host = _host(tmp_path)

    def declare(item):
        return host.dispatch(
            turn_id=f"t-{item['observable_id']}",
            tool_name="declare_requested_observable",
            arguments={"observables": [item]},
        )["result"]

    declare(
        {
            "observable_id": "spin-ni",
            "unit": "e",
            "meaning": "spin population on nickel",
        }
    )
    with pytest.raises(ContractError, match="names no observable"):
        declare(
            {
                "observable_id": "spinpop-ni",
                "unit": "1",
                "meaning": "spin population on nickel",
                "supersedes_observable_id": "never-declared",
            }
        )
    corrected = declare(
        {
            "observable_id": "spinpop-ni",
            "unit": "1",
            "meaning": "spin population on nickel, dimensionless",
            "supersedes_observable_id": "spin-ni",
        }
    )
    (record,) = corrected["declared"]
    assert record["supersedes_observable_id"] == "spin-ni"
    assert superseded_observable_ids(
        tuple(host.requested_observable_declarations.values())
    ) == frozenset({"spin-ni"})

    task = next(iter(host.task_spec_sha256s))
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256=task, claims=(_claim("1", (0,) * 6),)
    )
    misses, limitations = host._declared_observable_completion(
        task_spec_sha256=task
    )
    # The retired id is not owed; the corrected one is, and is open
    # until a claim carries its id.
    assert "spin-ni" not in " ".join(limitations)
    assert limitations == ("declared_observable:spinpop-ni",)


def test_the_rewake_reads_delivery_the_way_the_settlement_does(tmp_path):
    """A third reader had grown alone.

    The repair above gave the completion gate and the settlement one
    function. ``_rewake`` -- the organ that decides whether a goal gets
    another cycle at all -- still built its own id-only union from the
    same stream. So a claim in the wrong dimension looked delivered
    there and undelivered here: the cycle that could have repaired it
    was never opened, and the goal settled naming it. That is the
    defect this module exists for, one organ over.
    """

    import json

    from chemsmart.agent.driver import _analysis_delivery

    declared = {
        "observable_id": "spin-ni",
        "unit": "e",
        "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
        "meaning": "spin population on nickel",
    }
    stream = tmp_path / "events.jsonl"
    stream.write_text(
        json.dumps(
            {
                "kind": "requested_observable_declared",
                "payload": {"observables": [declared]},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    # The record row carries the id and the wrong dimension.
    delivery = _analysis_delivery(
        stream,
        goal_delivered_ids={
            "spin-ni": {"unit": "1", "dimension": (0,) * 9, "value": 0.81}
        },
        declared_observables=(declared,),
    )
    assert not delivery._delivered_by_the_goal("spin-ni")

    # And in the dimension it was declared in, it is delivered.
    delivered = _analysis_delivery(
        stream,
        goal_delivered_ids={
            "spin-ni": {
                "unit": "e",
                "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
                "value": 0.81,
            }
        },
        declared_observables=(declared,),
    )
    assert delivered._delivered_by_the_goal("spin-ni")


def test_a_current_cycle_claim_is_judged_in_its_declared_dimension(tmp_path):
    """ac127b41 repaired half of the re-wake's join.

    Record rows from earlier cycles went through the shared predicate;
    claims of the current stream were then added on their id alone. So a
    claim of this cycle in the wrong dimension counted as delivered
    here and undelivered at settlement -- and the cycle that could have
    repaired it was never opened, which is the same defect this module
    exists for, in the same organ, one source over.

    A claim written before dimensions travelled carries only its display
    unit; the predicate resolves that through the unit table, so this
    must not refuse it.
    """

    import json

    from chemsmart.agent.driver import _analysis_delivery

    declared = {
        "observable_id": "spin-ni",
        "unit": "e",
        "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
        "meaning": "spin population on nickel",
    }

    def delivery_for(claim):
        stream = tmp_path / f"events-{claim['display_unit']}.jsonl"
        stream.write_text(
            json.dumps(
                {
                    "kind": "requested_observable_declared",
                    "payload": {"observables": [declared]},
                }
            )
            + "\n"
            + json.dumps(
                {
                    "kind": "analysis_claims_recorded",
                    "payload": {
                        "receipt_sha256": "3" * 64,
                        "record": {"claims": [claim]},
                    },
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return _analysis_delivery(stream, declared_observables=(declared,))

    # Right dimension, carried explicitly: delivered.
    right = delivery_for(
        {
            "claim_id": "spin-ni",
            "quantity_id": "spin-ni",
            "display_value": 0.81,
            "display_unit": "e",
            "dimension": (0, 0, 0, 0, 0, 0, 0, 0, 1),
        }
    )
    assert right.answers_declaration("spin-ni")

    # Wrong dimension: not delivered, so the wake can open.
    wrong = delivery_for(
        {
            "claim_id": "spin-ni",
            "quantity_id": "spin-ni",
            "display_value": 0.81,
            "display_unit": "1",
            "dimension": (0,) * 9,
        }
    )
    assert not wrong.answers_declaration("spin-ni")
