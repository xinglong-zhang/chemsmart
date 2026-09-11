"""A self-report cannot discharge a contract that judges the self-report.

SUFFICIENCY-1 closed three requirements on `uncertainty: 0.2` against
`required_tolerance: 0.2`, basis `asserted`, while its own recorded
decision named a second term it had not quantified. The arithmetic was
right and the contract was empty.

Nothing here penalises an assertion. That session had a *measured*
spread available -- 0.013 V across three vibrational treatments -- and
declined it to assert an honest 0.2 V, which was the better science.
An assertion simply does not close an evidence obligation on its own:
it routes to measuring it, to showing the decision does not turn on it,
or to saying it cannot be established here.

The ontology stays the session's. The host reads provenance and never
meaning, never combines components, and imposes no vocabulary of error
kinds.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.delivery import judge_sufficiency
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_KJ_PER_MOL_IN_HARTREE = 2625.4996394798254


pytestmark = pytest.mark.capability("tool:record_analysis_claims")

_TASK = "a" * 64
_RECEIPT = "b" * 64


def _host(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir(exist_ok=True)
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=(_TASK,),
        approved_workspace=workspace,
    )
    host.quantity_extractions[_RECEIPT] = SimpleNamespace(
        quantities=(
            SimpleNamespace(
                quantity_id="dg",
                value=-10.0,
                unit="kJ/mol",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="9" * 64,
            ),
            # The spread a `measured` uncertainty cites: the host reads
            # this number where the session says it got it.
            SimpleNamespace(
                quantity_id="dg_spread",
                # Canonical, as a real receipt carries it: the host
                # converts to the claim's display unit exactly as it
                # does for the value being claimed.
                value=1.0 / _KJ_PER_MOL_IN_HARTREE,
                unit="hartree",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="7" * 64,
            ),
            # A second quantity no declaration puts a tolerance on, so
            # the quiet case is quiet for the right reason.
            SimpleNamespace(
                quantity_id="dg_other",
                value=-11.0,
                unit="kJ/mol",
                dimension=(1, 0, 0, 0, 0, 0),
                data_kind="scalar",
                value_sha256="8" * 64,
            ),
        )
    )
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "dg",
                    "unit": "kJ/mol",
                    "meaning": "solvation free energy",
                    "required_tolerance": 2.0,
                    "tolerance_basis": "the author's design threshold",
                }
            ]
        },
    )
    return host


def _claim(host, **extra):
    return host._record_analysis_claims(
        "t1",
        {
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    **extra,
                }
            ],
        },
    )


def _rows(host, tmp_path):
    import json

    return [
        row
        for line in (tmp_path / "events.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
        for event in (json.loads(line),)
        if event["kind"] == "analysis_claims_recorded"
        for row in event["payload"].get("sufficiency") or ()
    ]


def test_a_word_that_names_evidence_must_name_the_evidence(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="needs uncertainty_reference"):
        _claim(host, uncertainty=1.0, uncertainty_basis="measured")
    with pytest.raises(ContractError, match="takes no uncertainty_reference"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="asserted",
            uncertainty_reference=f"{_RECEIPT}:dg_spread",
        )
    with pytest.raises(ContractError, match="does not resolve"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="inferred",
            uncertainty_reference="c" * 64,
        )


def test_a_resolved_reference_discharges_and_an_assertion_does_not(tmp_path):
    host = _host(tmp_path)
    _claim(host, uncertainty=1.0, uncertainty_basis="asserted")
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="measured",
        uncertainty_reference=f"{_RECEIPT}:dg_spread",
    )
    states = [row["state"] for row in _rows(host, tmp_path)]
    assert states == ["attested", "met"]


def test_a_named_unquantified_term_holds_the_requirement_open(tmp_path):
    """Naming what you cannot quantify must not be the losing move.

    It holds the requirement open, which is what such a term means --
    and the bounded refusal is what makes that honest rather than a
    trap.
    """

    host = _host(tmp_path)
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="measured",
        uncertainty_reference=f"{_RECEIPT}:dg_spread",
        uncertainty_components=[
            {
                "meaning": "vibrational treatment spread",
                "magnitude": 1.0,
                "basis": "measured",
                "reference": _RECEIPT,
            },
            {
                "meaning": "basin coverage of the cation; no alternative "
                "conformer was searched",
                "basis": "unquantified",
            },
        ],
    )
    (row,) = _rows(host, tmp_path)
    assert row["state"] == "attested"
    assert row["unquantified_components"] == [
        "basin coverage of the cation; no alternative conformer was searched"
    ]


def test_a_component_states_a_magnitude_or_says_it_cannot(tmp_path):
    host = _host(tmp_path)
    with pytest.raises(ContractError, match="needs a magnitude"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="asserted",
            uncertainty_components=[
                {"meaning": "functional systematic", "basis": "asserted"}
            ],
        )
    with pytest.raises(ContractError, match="carries no magnitude"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="asserted",
            uncertainty_components=[
                {
                    "meaning": "functional systematic",
                    "magnitude": 0.4,
                    "basis": "unquantified",
                }
            ],
        )


def test_provenance_follows_the_chain_to_its_roots(tmp_path):
    """An expression receipt is host arithmetic; its leaves need not be.

    A value the session supplied as a literal stays the session's
    however many receipts sit above it, so a budget composed over
    literals cannot launder itself into evidence.
    """

    host = _host(tmp_path)
    # Evaluated through the tool. The version of this test that
    # hand-populated `quantity_expression_requests` passed over a
    # production path that restored no request at all, which is the
    # fixture-supplies-what-production-lacks pattern it was written to
    # close.
    reply = host.dispatch(
        turn_id="t1",
        tool_name="evaluate_quantity_expression",
        arguments={
            "expression_id": "budget",
            "inputs": [],
            "nodes": [
                {
                    "node_id": "n1",
                    "operation": "literal",
                    "literal_value": 0.2,
                    "literal_unit": "kJ/mol",
                }
            ],
            "output_node_ids": ["n1"],
        },
    )
    literal_expression = reply["result"]["receipt_sha256"]
    # The chain is still walked to its roots and the authorship is
    # still read off the receipt. What changed on 2026-09-10 is what
    # the host does with the answer: it reports it and refuses nothing,
    # because the refusal read spelling rather than value and taught
    # receipt partitioning rather than derivation.
    backed, detail = host._uncertainty_evidence(literal_expression)
    assert backed is True
    assert host._expression_is_model_authored(literal_expression) is True

    # An extraction receipt is the host's own reading and does resolve.
    assert host._uncertainty_evidence(_RECEIPT)[0] is True

    # And a host rebuilt over the same durable stream -- every
    # continuation is one -- reaches the same word. The receipt carries
    # node values and no operations, so a restored receipt left the
    # guard nothing to walk and the literal became host-owned evidence.
    resumed = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=(_TASK,),
        approved_workspace=tmp_path / "workspace",
    )
    assert literal_expression in resumed.quantity_expression_receipts
    assert resumed._uncertainty_evidence(literal_expression)[0] is True
    # The authorship itself is what must survive a restart: the receipt
    # carries node values and no operations, so a restored receipt once
    # left the walk nothing to read and the literal became host-owned
    # without observation. It is the report that rests on this now.
    assert resumed._expression_is_model_authored(literal_expression) is True


def test_the_declaration_and_the_claim_agree_on_nothing_but_provenance():
    """The host reads provenance and never meaning: two budgets with the
    same numbers and different sentences are the same to it."""

    declaration = {
        "observable_id": "dg",
        "unit": "kJ/mol",
        "required_tolerance": 2.0,
    }
    for meaning in ("continuum solvation", "a term I will not name here"):
        row = {
            "uncertainty": 1.0,
            "uncertainty_basis": "measured",
            "uncertainty_evidence_backed": True,
            "unquantified_components": [],
            "meaning": meaning,
        }
        assert judge_sufficiency(declaration, row)["state"] == "met"


def test_the_host_tells_the_session_what_it_made_of_its_own_number(tmp_path):
    """The verdict reached the stream and never the reply.

    A session could not learn in the turn it claimed that its number
    missed the tolerance it had itself declared: it had to spend a whole
    cycle to be told through the wake, while two of the three routes
    that answer cost no engine call. An independent audit named this the
    cheapest available improvement in the mechanism, and it is, because
    the reply is where this model is actually taught.
    """

    host = _host(tmp_path)
    reply = host.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    "uncertainty": 6.0,
                    "uncertainty_basis": "asserted",
                }
            ],
        },
    )
    (row,) = reply["observations"]
    assert row["observable_id"] == "dg"
    assert row["state"] == "short"
    assert row["required_tolerance"] == 2.0

    # A claim on an observable that asked for no precision says nothing.
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "free",
                    "unit": "kJ/mol",
                    "meaning": "an observable with no stated tolerance",
                }
            ]
        },
    )
    quiet = host.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "free",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg_other",
                    "display_unit": "kJ/mol",
                }
            ],
        },
    )
    assert "observations" not in quiet

    # But a claim that *delivers* a tolerance-bearing observable is
    # assessed however it is labelled. The gate credits a claim under
    # both its claim_id and its host-minted quantity_id, and this join
    # took only the first: claiming the quantity `dg`, which carries a
    # tolerance, under the tolerance-free label `free` delivered both
    # ids and was assessed against neither. The completion certified,
    # no sufficiency row existed, and the requirement left the question
    # rather than standing open. An earlier version of this test
    # asserted that silence.
    aliased = host.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "free",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                }
            ],
        },
    )
    (aliased_row,) = aliased["observations"]
    assert aliased_row["observable_id"] == "dg"
    assert aliased_row["state"] == "unstated"


def test_every_open_state_has_a_route_and_an_unknown_one_is_refused():
    """A state without a route is a state the host cannot ask about.

    The menu fell through to the "compute it" route for anything it did
    not recognise -- which is the route the function exists to stop
    offering blindly, since a session that has named no term cannot
    narrow one. `attested` shipped in this round with no route at all
    and re-opened goals with wording that was false of it.
    """

    from chemsmart.agent._contracts import ContractError as DriverError
    from chemsmart.agent.delivery import SUFFICIENCY_STATES
    from chemsmart.agent.driver import (
        SUFFICIENCY_ATTESTED_ROUTE,
        SUFFICIENCY_SHORT_ROUTE,
        SUFFICIENCY_UNSTATED_ROUTE,
        sufficiency_menu,
    )

    assert sufficiency_menu(["unstated"]) == SUFFICIENCY_UNSTATED_ROUTE
    assert sufficiency_menu(["attested"]) == SUFFICIENCY_ATTESTED_ROUTE
    assert sufficiency_menu(["short"]) == SUFFICIENCY_SHORT_ROUTE

    # Every open state the vocabulary admits is answerable.
    from chemsmart.agent.delivery import _OPEN_STATES

    for state in _OPEN_STATES:
        assert sufficiency_menu([state])
    assert set(_OPEN_STATES) <= set(SUFFICIENCY_STATES)

    with pytest.raises(DriverError, match="no route answers"):
        sufficiency_menu(["separated"])


def test_an_accuracy_this_envelope_cannot_reach_is_a_deliverable(tmp_path):
    """The third way an observable is unreachable.

    The verifier knew two kinds of absence: a selector no envelope
    program declares, and a blocked_unsupported node. Neither fits a
    precision insufficiency, where the producer exists and the number
    was computed and no method on the shelf establishes it to what was
    asked. So a session facing that would have had to invent an absence
    to say a true thing -- while the menu offered refusal as one of its
    three routes.

    This is the pair to the rule that an asserted uncertainty does not
    discharge. Without it, naming a term you cannot quantify holds the
    requirement open for ever and the cheapest response is to stop
    naming such terms, which is a trap rather than a design.
    """

    host = _host(tmp_path)
    _claim(host, uncertainty=6.0, uncertainty_basis="asserted")

    (row,) = host._verify_unreachable_observables(
        [
            {
                "observable_id": "dg",
                "statement": "no functional or basis available in this "
                "envelope establishes a solvation free energy to 2 kJ/mol; "
                "the continuum systematic alone exceeds it",
                "receipt_sha256s": [_RECEIPT],
            }
        ]
    )
    assert row["verified"] is True
    assert "stands open on this goal's own record" in row["basis"]
    # The host verifies the open requirement and claims nothing about
    # what no calculation could ever do.
    assert "never that no calculation could reach it" in row["basis"]


def test_a_precision_refusal_needs_the_requirement_to_be_open(tmp_path):
    """Refusing must not be a shortcut past claiming.

    A session that never states what its number is worth cannot refuse
    the precision of a number nobody has assessed.
    """

    host = _host(tmp_path)
    (row,) = host._verify_unreachable_observables(
        [
            {
                "observable_id": "dg",
                "statement": "cannot be established here",
                "receipt_sha256s": [_RECEIPT],
            }
        ]
    )
    assert row["verified"] is False
    assert "no open assessment stands against it" in row["basis"]

    # And once it is met, there is nothing left to refuse.
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="measured",
        uncertainty_reference=f"{_RECEIPT}:dg_spread",
    )
    (met_row,) = host._verify_unreachable_observables(
        [
            {
                "observable_id": "dg",
                "statement": "cannot be established here",
                "receipt_sha256s": [_RECEIPT],
            }
        ]
    )
    assert met_row["verified"] is False


def test_the_evidence_a_word_names_is_kept_on_the_claim(tmp_path):
    """The reference and the budget ride inside the claim's digest.

    They were resolved, consumed into one boolean and dropped, so the
    only durable trace of `measured` was `uncertainty_evidence_backed:
    true` -- and nothing afterwards could ask which receipt backed the
    number that discharged the contract. Inside the digest for the same
    reason the uncertainty is: re-claiming with a corrected budget must
    mint a second receipt rather than collide on the idempotency key,
    which is the route the attested wake offers.
    """

    host = _host(tmp_path)
    reference = f"{next(iter(host.quantity_extractions))}:dg_spread"
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="measured",
        uncertainty_reference=reference,
        uncertainty_components=[
            {
                "meaning": "vibrational treatment spread",
                "magnitude": 0.6,
                "basis": "measured",
                "reference": reference,
            }
        ],
    )
    (record,) = tuple(host.analysis_claim_records.values())
    (claim,) = record.claims
    assert claim.uncertainty_reference == reference
    assert claim.uncertainty_components[0]["meaning"] == (
        "vibrational treatment spread"
    )
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="measured",
        uncertainty_reference=reference,
        uncertainty_components=[
            {
                "meaning": "vibrational treatment spread",
                "magnitude": 0.6,
                "basis": "measured",
                "reference": reference,
            },
            {
                "meaning": "cation basin coverage, searched and null",
                "magnitude": 0.0,
                "basis": "measured",
                "reference": reference,
            },
        ],
    )
    assert len(host.analysis_claim_records) == 2


def test_a_components_reference_resolves_as_the_claims_does(tmp_path):
    """The schema promises the host reads a component's provenance.

    It read only that the string was non-empty, so a component could
    name evidence that does not exist under a word that claims it does.
    """

    host = _host(tmp_path)
    with pytest.raises(ContractError) as excinfo:
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="asserted",
            uncertainty_components=[
                {
                    "meaning": "vibrational treatment spread",
                    "magnitude": 0.4,
                    "basis": "measured",
                    "reference": "not-a-receipt",
                }
            ],
        )
    assert "does not resolve" in str(excinfo.value)
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="asserted",
        uncertainty_components=[
            {
                "meaning": "vibrational treatment spread",
                "magnitude": 0.4,
                "basis": "measured",
                "reference": _RECEIPT,
            }
        ],
    )


def test_a_supersession_carries_the_tolerance_it_replaces(tmp_path):
    """A supersession corrects the observable, never the requirement.

    Requiring the replacement to state *a* tolerance preserved its
    presence and not the obligation: a 2 kJ/mol requirement was
    replaced by 200 kJ/mol under the same meaning, the same unit and
    the same quoted source, and the host retired the original. The
    model authored its finish line one declaration later. The
    comparison is made at the precision the session wrote, so 0.05 eV
    restated as 1.15 kcal/mol is the same requirement and 200 for 2 is
    a different one.
    """

    host = _host(tmp_path)
    with pytest.raises(ContractError, match="how good the answer has to be"):
        host._declare_requested_observable(
            "t1",
            {
                "observables": [
                    {
                        "observable_id": "dg_loose",
                        "unit": "kJ/mol",
                        "meaning": "solvation free energy",
                        "required_tolerance": 200.0,
                        "tolerance_basis": "the author's design threshold",
                        "supersedes_observable_id": "dg",
                    }
                ]
            },
        )
    # A genuine unit correction still walks: 2 kJ/mol is 0.478 kcal/mol.
    host._declare_requested_observable(
        "t1",
        {
            "observables": [
                {
                    "observable_id": "dg_kcal",
                    "unit": "kcal/mol",
                    "meaning": "solvation free energy, corrected unit",
                    "required_tolerance": 0.478,
                    "tolerance_basis": "2 kJ/mol restated in kcal/mol",
                    "supersedes_observable_id": "dg",
                }
            ]
        },
    )
    record = host.requested_observable_declarations["dg_kcal"]
    assert record["superseded_required_tolerance"] == 2.0
    assert record["superseded_tolerance_unit"] == "kJ/mol"


def test_only_a_host_checkable_magnitude_discharges(tmp_path):
    """`met` had meant a small number beside a resolvable citation.

    An uncertainty of 0.01 kJ/mol citing a receipt whose quantity was an
    unrelated standard-state correction reached it. Binding a number to
    its numerical source is provenance and the host owns that; judging
    whether it estimates the relevant scientific error is chemistry and
    stays the session's. ChemSmart already draws that line for the
    claimed value, which the model never types, and the uncertainty is
    now held to the same rule (owner ruling, 2026-09-09).
    """

    host = _host(tmp_path)

    # A bare receipt no longer discharges: `measured` names the quantity.
    with pytest.raises(ContractError, match="names the quantity too"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="measured",
            uncertainty_reference=_RECEIPT,
        )

    # Understating what the cited quantity reads is refused, because
    # that is the escape the check exists to close.
    with pytest.raises(ContractError, match="means the host reads"):
        _claim(
            host,
            uncertainty=0.1,
            uncertainty_basis="measured",
            uncertainty_reference=f"{_RECEIPT}:dg_spread",
        )

    # Cite the quantity and omit the number: the host copies it in,
    # exactly as it copies the value being claimed.
    reply = host.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    "uncertainty_basis": "measured",
                    "uncertainty_reference": f"{_RECEIPT}:dg_spread",
                }
            ],
        },
    )
    (row,) = reply["observations"]
    assert row["state"] == "met"
    assert round(float(row["uncertainty"]), 6) == 1.0

    # An inference stands on a receipt and is still the session's
    # judgement, so it attests rather than discharges.
    (tmp_path / "two").mkdir(exist_ok=True)
    host_two = _host(tmp_path / "two")
    inferred = host_two.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    "uncertainty": 1.0,
                    "uncertainty_basis": "inferred",
                    "uncertainty_reference": _RECEIPT,
                }
            ],
        },
    )
    (inferred_row,) = inferred["observations"]
    assert inferred_row["state"] == "attested"


def test_the_host_reports_what_it_saw_and_rules_on_none_of_it(tmp_path):
    """Three refusals became one report (owner ruling, 2026-09-10).

    A zero magnitude, a spread over a single receipt, and a chain
    naming a coefficient the session supplied were each refused. Each
    check was computable and none was scientifically defensible: a
    spread of exactly zero is a real observation, a variance over many
    samples inside one receipt compares plenty, and an equivalent
    coefficient does not become worth less because operators spelled
    it -- the vocabulary is rational-complete over any non-zero
    quantity, so the authored-constant rule read spelling and not
    value. None of the three stopped the substitution it was aimed at.

    So the host reports what it observed beside the assessment and
    refuses none of it. The requirement can now be discharged by such a
    number, and a reader of that `met` is told what it rests on -- on
    the assessment, which is what the session, the workspace record and
    the settlement all read.
    """

    host = _host(tmp_path)
    zero = SimpleNamespace(
        quantity_id="fake_spread",
        value=0.0,
        unit="hartree",
        dimension=(1, 0, 0, 0, 0, 0),
        data_kind="scalar",
        value_sha256="6" * 64,
    )
    host.quantity_expression_receipts["f" * 64] = SimpleNamespace(
        outputs=(zero,),
        output_dependencies=(
            SimpleNamespace(
                output_id="fake_spread",
                source_receipt_sha256s=(_RECEIPT, "e" * 64),
                model_authored_constants=(),
            ),
        ),
    )
    _claim(
        host,
        uncertainty_basis="measured",
        uncertainty_reference=f"{'f' * 64}:fake_spread",
    )
    row = _rows(host, tmp_path)[-1]
    assert row["state"] == "met"
    observed = row["uncertainty_observations"]
    assert "uncertainty.zero_magnitude" in observed
    # And the report says what the evidence is rather than staying
    # silent: an empty list must never read as "assumption-free".
    assert any(o.startswith("uncertainty.source_receipts=") for o in observed)

    # A spread over one receipt, and a chain carrying a coefficient of
    # the session's own, are both admitted and both reported. The
    # authored-constant fact used to refuse the citation outright.
    lonely = SimpleNamespace(
        quantity_id="lonely",
        value=1.0 / _KJ_PER_MOL_IN_HARTREE,
        unit="hartree",
        dimension=(1, 0, 0, 0, 0, 0),
        data_kind="scalar",
        value_sha256="5" * 64,
    )
    host.quantity_expression_receipts["a" * 63 + "c"] = SimpleNamespace(
        outputs=(lonely,),
        output_dependencies=(
            SimpleNamespace(
                output_id="lonely",
                source_receipt_sha256s=(_RECEIPT,),
                # The real dataclass shape: a bare string here made
                # the test unfaithful to production, which is the
                # fixture-supplies-what-production-lacks pattern.
                model_authored_constants=(
                    SimpleNamespace(
                        node_id="lonely",
                        role="scale_factor",
                        value="0.5",
                    ),
                ),
            ),
        ),
    )
    _claim(
        host,
        uncertainty_basis="measured",
        uncertainty_reference=f"{'a' * 63 + 'c'}:lonely",
    )
    row = _rows(host, tmp_path)[-1]
    assert row["state"] == "met"
    observed = row["uncertainty_observations"]
    # The coefficient is named with its role and its value, because a
    # half that defines a half-range and a typed solvation budget are
    # both the session's and are not the same statement.
    assert any(
        o.startswith("uncertainty.model_authored_constants[")
        and "scale_factor=" in o
        for o in observed
    ), observed
    assert "uncertainty.source_receipts=1" in observed

    # An honest spread over two receipts carries no observation at all,
    # so the report distinguishes rather than decorating everything.
    clean = SimpleNamespace(
        quantity_id="clean",
        value=1.0 / _KJ_PER_MOL_IN_HARTREE,
        unit="hartree",
        dimension=(1, 0, 0, 0, 0, 0),
        data_kind="scalar",
        value_sha256="4" * 64,
    )
    host.quantity_expression_receipts["b" * 63 + "d"] = SimpleNamespace(
        outputs=(clean,),
        output_dependencies=(
            SimpleNamespace(
                output_id="clean",
                source_receipt_sha256s=(_RECEIPT, "e" * 64),
                model_authored_constants=(),
            ),
        ),
    )
    _claim(
        host,
        uncertainty_basis="measured",
        uncertainty_reference=f"{'b' * 63 + 'd'}:clean",
    )
    row = _rows(host, tmp_path)[-1]
    assert row["state"] == "met"
    observed = row["uncertainty_observations"]
    # An honest two-receipt spread names no coefficient of the
    # session's, and still says how many receipts it compared.
    assert not any("model_authored_constants" in o for o in observed)
    assert "uncertainty.source_receipts=2" in observed


def test_a_measured_claim_records_the_quantity_it_claims(tmp_path):
    """The `measured` branch rebound this loop's `quantity`.

    Every claim on the one path that can discharge a tolerance recorded
    the *uncertainty's* quantity_id, canonical value and value digest,
    so the claim said two different numbers at once and its digest no
    longer hashed the delivered one. It also reintroduced the alias
    bypass on that path alone: the declaration join read the
    uncertainty's id, so a tolerance-bearing observable delivered under
    an alternate label produced no assessment at all.
    """

    host = _host(tmp_path)
    host.dispatch(
        turn_id="t1",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": _TASK,
            "claims": [
                {
                    "claim_id": "dg",
                    "receipt_sha256": _RECEIPT,
                    "quantity_id": "dg",
                    "display_unit": "kJ/mol",
                    "uncertainty_basis": "measured",
                    "uncertainty_reference": f"{_RECEIPT}:dg_spread",
                }
            ],
        },
    )
    (record,) = tuple(host.analysis_claim_records.values())
    (claim,) = record.claims
    assert claim.quantity_id == "dg"
    assert claim.canonical_value == -10.0
    assert claim.quantity_value_sha256 == "9" * 64


def test_a_citation_is_read_at_the_value_it_names(tmp_path):
    """One receipt's other outputs are other numbers.

    The authored-constant check read every output row of the cited
    receipt, so a clean, fully derived term was refused because an
    unrelated output of the same expression carried a literal. Live: a
    session put the task's own three oxidant potentials in the same
    expression as its uncertainty budget, and its functional term --
    derived end to end -- was rejected, twice. What it learned was to
    partition receipts, not to derive anything. Forty lines below, the
    two-source guard already scoped by output id: two readers of one
    citation, one of them wrong.
    """

    host = _host(tmp_path)
    clean = SimpleNamespace(
        quantity_id="clean",
        value=1.0 / _KJ_PER_MOL_IN_HARTREE,
        unit="hartree",
        dimension=(1, 0, 0, 0, 0, 0),
        data_kind="scalar",
        value_sha256="4" * 64,
    )
    authored = SimpleNamespace(
        quantity_id="authored",
        value=0.27,
        unit="hartree",
        dimension=(1, 0, 0, 0, 0, 0),
        data_kind="scalar",
        value_sha256="5" * 64,
    )
    host.quantity_expression_receipts["d" * 64] = SimpleNamespace(
        outputs=(clean, authored),
        output_dependencies=(
            SimpleNamespace(
                output_id="clean",
                source_receipt_sha256s=(_RECEIPT, "e" * 64),
                model_authored_constants=(),
            ),
            SimpleNamespace(
                output_id="authored",
                source_receipt_sha256s=(_RECEIPT, "e" * 64),
                model_authored_constants=(
                    SimpleNamespace(role="literal_value"),
                ),
            ),
        ),
    )
    # The scoping outlived the refusal it was built for: since
    # 2026-09-10 an authored coefficient is reported and not refused,
    # so what the citation's own rows decide is whether the *report*
    # names this output or its neighbour. Reading the whole receipt
    # would now mark every clean term in a shared expression.
    clean_observed = host._uncertainty_observations("d" * 64, "clean")
    assert not any(
        "model_authored_constants" in o for o in clean_observed
    ), clean_observed
    authored_observed = host._uncertainty_observations("d" * 64, "authored")
    assert any(
        o.startswith("uncertainty.model_authored_constants[")
        for o in authored_observed
    ), authored_observed
    # Naming no output at all stays conservative: the host cannot see
    # which output was read, so it reports the receipt's authorship.
    assert host._expression_is_model_authored("d" * 64) is True
    # And resolution itself no longer turns on any of this.
    assert host._uncertainty_evidence("d" * 64, "clean")[0] is True
    assert host._uncertainty_evidence("d" * 64, "authored")[0] is True


def test_a_components_reference_names_a_quantity_that_exists(tmp_path):
    """The schema said the host resolves a component's reference as the
    claim's own. It resolved the receipt half and never asked whether
    the quantity existed, so a component could name evidence that does
    not exist under a word claiming it does. The magnitude stays
    ungraded: how the terms add is the science.
    """

    host = _host(tmp_path)
    with pytest.raises(ContractError, match="carries no quantity"):
        _claim(
            host,
            uncertainty=1.0,
            uncertainty_basis="asserted",
            uncertainty_components=[
                {
                    "meaning": "a term that names nothing",
                    "magnitude": 0.4,
                    "basis": "measured",
                    "reference": f"{_RECEIPT}:no-such-quantity",
                }
            ],
        )
    # And a component whose magnitude differs from the quantity it names
    # is admitted: a term stated at 2 sigma from a 1 sigma receipt is
    # ordinary, and grading it would be the host judging chemistry.
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="asserted",
        uncertainty_components=[
            {
                "meaning": "the same term at a different coverage",
                "magnitude": 2.0,
                "basis": "measured",
                "reference": f"{_RECEIPT}:dg_spread",
            }
        ],
    )


@pytest.mark.capability("tool:record_analysis_claims")
def test_a_component_says_what_the_quantity_it_cites_reads(tmp_path):
    """A component's magnitude is not graded, and is not hidden either.

    The claim's own ``measured`` magnitude is checked against its
    citation. A component's deliberately is not: a term restated at two
    sigma from a one sigma receipt is ordinary science, and the host
    cannot tell that apart from a wrong citation without judging which
    arithmetic was meant -- which is the session's. ``ino3-r17`` cycle 3
    (2026-09-10) stated a 0.214 V component whose cited quantity holds
    -2.103 V, an endpoint rather than any spread, under a rule reading
    "the measured SVP-to-TZVP spread ... rounded down to 0.15 V". No
    refusal belongs there. What belongs there is the number, on the row,
    so the human who reads the claim can see the two do not match.
    """

    host = _host(tmp_path)
    _claim(
        host,
        uncertainty=1.0,
        uncertainty_basis="asserted",
        uncertainty_components=[
            {
                "meaning": "the same term at a different coverage",
                "magnitude": 2.0,
                "basis": "measured",
                "reference": f"{_RECEIPT}:dg_spread",
            }
        ],
    )
    (record,) = tuple(host.analysis_claim_records.values())
    (claim,) = record.claims
    observed = tuple(claim.uncertainty_components[0].get("observations") or ())
    reads = [o for o in observed if "cited_quantity_reads" in o]
    assert reads, (
        "the component's cited quantity is not reported, so a reader "
        f"cannot see it differs from the stated magnitude: {observed}"
    )
    assert "dg_spread" in reads[0], reads
    # And it is an observation, never a verdict: the claim stands.
    assert claim.uncertainty_components[0]["magnitude"] == 2.0
