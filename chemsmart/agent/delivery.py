"""One judgement of whether a declared observable was delivered.

The completion gate and the settlement are two readers of the same
goal, and they disagreed. OPEN-1 ino3 (2026-09-07) declared six spin
observables in ``e`` and claimed them as dimensionless counts: the
completion gate recorded six declared-observable limitations naming the
dimension mismatch, and the settlement, reading the workspace record by
id alone, called the same six delivered in an earlier cycle and settled
``achieved``. Two auditors reproduced it independently.

The predicate lives here so both call it. It compares what the gate has
always compared -- the id and the dimension, zero-padded, because
dimension vectors from different eras differ only by trailing bases --
and nothing else: a value is never checked, because a declaration is a
contract about meaning, not about the answer.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

#: Dimension vectors are compared zero-padded to this width, so a
#: six-base row written by an earlier era matches a nine-base one.
_DIMENSION_WIDTH = 9


def padded_dimension(dimension: Sequence[Any]) -> tuple[int, ...]:
    values = tuple(int(value) for value in dimension)
    return values + (0,) * (_DIMENSION_WIDTH - len(values))


def _dimension_of(record: Mapping[str, Any]) -> tuple[int, ...] | None:
    """The dimension a declaration or a claim row carries.

    A row written before this repair carries only its display unit, so
    the unit is resolved through the same table the analysis layer
    uses. A unit the table does not know yields None, and an unknown
    dimension never satisfies a declaration.
    """

    raw = record.get("dimension")
    if isinstance(raw, (list, tuple)) and raw:
        try:
            return padded_dimension(raw)
        except (TypeError, ValueError):
            return None
    unit = str(record.get("unit") or record.get("display_unit") or "").strip()
    if not unit:
        return None
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionError,
        unit_dimension,
    )

    try:
        return padded_dimension(unit_dimension(unit))
    except (QuantityExpressionError, KeyError, ValueError):
        return None


def observable_is_delivered(
    declaration: Mapping[str, Any], claim_row: Mapping[str, Any] | None
) -> bool:
    """True when this claim row answers this declaration.

    The row must exist and carry the declaration's dimension. The id is
    the caller's join: both callers look the row up by the declared id.
    """

    if not claim_row:
        return False
    declared = _dimension_of(declaration)
    delivered = _dimension_of(claim_row)
    if declared is None or delivered is None:
        return False
    return declared == delivered


def superseded_observable_ids(
    declarations: Sequence[Mapping[str, Any]],
) -> frozenset[str]:
    """Ids a later declaration explicitly retired.

    A mistaken unit must be repairable without re-running physics: the
    session declares the corrected observable and names the one it
    replaces, and the retired id stops being owed. The correction is
    recorded, so a reader sees both.
    """

    retired = {
        str(item.get("supersedes_observable_id") or "")
        for item in declarations
        if item.get("supersedes_observable_id")
    }
    return frozenset(item for item in retired if item)


def declarations_by_id(
    declarations: Sequence[Mapping[str, Any]],
) -> dict[str, Mapping[str, Any]]:
    return {
        str(item.get("observable_id") or ""): item
        for item in declarations
        if item.get("observable_id")
    }


#: How a requirement stands once a number exists. Only ``met`` is
#: resolved; ``attested``, ``short`` and ``unstated`` are open, and
#: an open requirement with budget in hand is what re-opens a goal.
SUFFICIENCY_STATES = ("met", "attested", "short", "unstated")


def judge_sufficiency(
    declaration: Mapping[str, Any] | None,
    claim_row: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    """Whether a delivered number answers the precision it was asked for.

    Arithmetic on numbers the model wrote and the host copied; nothing
    here judges chemistry, and nothing here checks a value against an
    expectation. Three things can be true of a requirement:

    ``met``
        the stated uncertainty is within the tolerance the task asked
        for, **and** the uncertainty rests on evidence the host
        resolved with no component left unquantified. Only then is the
        requirement discharged.
    ``attested``
        the stated uncertainty is within the tolerance and rests on the
        session's own word, or leaves a component unquantified. The
        number stands and the requirement stays open: SUFFICIENCY-1
        closed three requirements on ``uncertainty: 0.2`` against
        ``required_tolerance: 0.2``, basis ``asserted``, while its own
        decision named a second term it had not quantified. A
        self-report cannot discharge a contract that judges the
        self-report. Nothing here punishes an assertion -- an asserted
        0.2 V was the honest number, and the session declined a
        measured 0.013 V spread to give it -- it simply does not close
        the obligation on its own.
    ``short``
        the stated uncertainty exceeds the tolerance. The requirement is
        open, and the number stays delivered.

    There is deliberately no state for "the decision this number
    supports is unaffected". A ``decision_boundary`` field once carried
    that, and it was the wrong shape twice over: a model-typed literal
    where the plane already offers a host-evaluated quantity, and a
    slot that anticipated an idea rather than absorbing one.
    SUFFICIENCY-1's session showed the better route without being asked
    -- it composed the margin to ferrocenium as a quantity and declared
    it as its own observable with its own tolerance, which travels
    through the ordinary machinery with real provenance. A decision
    question is a quantity with a tolerance, and needs nothing else.

    A claim that states no uncertainty at all is ``unstated`` and open
    too: silence is not sufficiency. A declaration with no
    ``required_tolerance`` asks nothing and returns ``None``.
    """

    if not declaration:
        return None
    tolerance = declaration.get("required_tolerance")
    if tolerance is None:
        return None
    row: Mapping[str, Any] = claim_row or {}
    uncertainty = row.get("uncertainty")
    verdict: dict[str, Any] = {
        "observable_id": str(declaration.get("observable_id") or ""),
        "unit": str(declaration.get("unit") or ""),
        "required_tolerance": float(tolerance),
        # Whose requirement this is. A tolerance the session formulated
        # for a decision of its own is legitimate and is not the
        # precision the requester asked for, and every reader that
        # decides anything read the number without that fact
        # (SUFFICIENCY-5, 2026-09-10).
        "tolerance_origin": str(
            declaration.get("tolerance_origin") or "unstated"
        ),
        "uncertainty": (None if uncertainty is None else float(uncertainty)),
        "uncertainty_basis": str(row.get("uncertainty_basis") or ""),
    }
    if uncertainty is None:
        verdict["state"] = "unstated"
        return verdict
    verdict["meets_tolerance"] = float(uncertainty) <= float(tolerance)
    # What the host established about the uncertainty itself, written
    # by the claim handler: whether its magnitude resolved through a
    # reference the host owns, and whether the session left any
    # component of its own budget unquantified.
    evidence_backed = bool(row.get("uncertainty_evidence_backed"))
    unquantified = tuple(row.get("unquantified_components") or ())
    verdict["uncertainty_evidence_backed"] = evidence_backed
    # And what the host observed about the cited magnitude without
    # ruling on it. Copied onto the verdict explicitly because this
    # function builds a fresh dict and returns early on `met` -- the
    # one state where the observation matters most, since a `met`
    # standing on a zero spread or a single receipt is exactly what the
    # two removed guards used to refuse (owner ruling, 2026-09-10).
    observations = tuple(
        str(item) for item in row.get("uncertainty_observations") or ()
    )
    if observations:
        verdict["uncertainty_observations"] = list(observations)
    # Which rule turned these components into this total, what its
    # inputs mean, what coverage it claims and what it assumes about
    # dependence. Copied and never graded: identical components gave
    # 0.249 by root-sum-square and 0.427 by linear sum against the same
    # tolerance, and the verdict named neither. Absence is reported as
    # absence rather than assumed independent.
    combination = row.get("uncertainty_combination")
    verdict["uncertainty_combination"] = (
        dict(combination) if isinstance(combination, Mapping) else None
    )
    if unquantified:
        verdict["unquantified_components"] = list(unquantified)
    if verdict["meets_tolerance"]:
        verdict["state"] = (
            "met" if evidence_backed and not unquantified else "attested"
        )
        if verdict["state"] == "met":
            return verdict
    if verdict.get("state") == "attested":
        return verdict
    verdict["state"] = "short"
    return verdict


#: An assessment in one of these states has not answered its
#: requirement. ``unstated`` is open because silence is not
#: sufficiency; ``attested`` is open because a self-report cannot
#: discharge a contract that judges the self-report.
_OPEN_STATES = frozenset({"short", "unstated", "attested"})


def current_assessments(
    verdicts: Sequence[Mapping[str, Any]],
) -> dict[str, Mapping[str, Any]]:
    """The latest assessment of each requirement, in first-seen order.

    History is evidence; it is not state. The first version of this
    module kept every row and called an id open if *any* row was open,
    so a requirement that was assessed ``unstated``, re-claimed with a
    proper uncertainty, and assessed ``met`` stayed open for ever --
    which is the one route the re-wake offers, unable to close the wake
    that offered it. A repair whose own condition cannot be discharged
    is the second defect of that repair.
    """

    latest: dict[str, Mapping[str, Any]] = {}
    for verdict in verdicts:
        observable_id = str(verdict.get("observable_id") or "")
        if not observable_id:
            continue
        latest[observable_id] = verdict
    return latest


def unresolved_requirement_ids(
    verdicts: Sequence[Mapping[str, Any]],
) -> tuple[str, ...]:
    """Requirements whose current assessment has not answered them.

    A state this vocabulary does not know is open, not resolved. Two
    readers used to disagree in opposite directions -- this one closed
    an unknown state silently while the menu refused one loudly -- and
    the rows are read back out of a durable workspace record that a
    previous commit on this branch wrote ``separated`` into. The
    fail-open half is the dangerous half: it discharges a contract
    nothing assessed.
    """

    return tuple(
        observable_id
        for observable_id, verdict in current_assessments(verdicts).items()
        if str(verdict.get("state") or "") != "met"
    )


def refused_requirement_ids(
    verdicts: Sequence[Mapping[str, Any]],
    verified_unreachable_ids: Sequence[str],
) -> tuple[str, ...]:
    """Open requirements the session refused and the host verified.

    Route three of the sufficiency menu, and the charter calls its
    ending a deliverable. The subtraction that closes these ids was
    written before anything read the difference, so a verified
    precision refusal removed the requirement from every projection and
    the goal settled ``achieved`` -- the refusal bought the best word
    in the vocabulary. They are named here so a settlement can say
    ``unreachable_from_evidence`` over them instead.
    """

    refused = set(verified_unreachable_ids)
    return tuple(
        observable_id
        for observable_id in unresolved_requirement_ids(verdicts)
        if observable_id in refused
    )
