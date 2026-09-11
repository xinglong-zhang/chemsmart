"""Planning contracts for connected calculation and analysis tool chains.

This module models both future data flow and analysis of results that ChemSmart
has already registered. Future values bind producer node/output names. Existing
results enter the DAG as registered artifact inputs, without a fictitious
calculation producer or a model-authored path.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Sequence

from chemsmart.agent._contracts import (
    DIMENSIONLESS_UNIT_HINT,
    ContractError,
    canonical_sha256,
    require_identifier,
)
from chemsmart.agent.postprocessing import typed_result_artifact_kind
from chemsmart.agent.workflows import CommandNodeIntentV1
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionError,
    normalize_numeric_value,
)
from chemsmart.analysis.result_quantities import (
    DIMENSIONLESS,
    canonical_thermochemistry_quantity,
    derivable_thermochemistry_quantities,
)

ANALYSIS_INTENT_KINDS = (
    "claim_rendering",
    "quantity_expression",
    "result_extraction",
    "scientific_validation",
    "thermochemistry",
    "unsupported_external",
)
ANALYSIS_INTENT_SUPPORT_STATES = ("blocked_unsupported", "planned")
# Registered analysis operations for these kinds read a typed program result
# artifact; a derived quantity from another analysis node cannot stand in.
RESULT_ARTIFACT_ANALYSIS_KINDS = ("result_extraction", "thermochemistry")
ANALYSIS_VALIDATION_PREDICATES = (
    "all_equal",
    "all_equal_text",
    "all_finite",
    "count_equals",
    "integer_equals",
    "maximum_absolute_less_equal",
    "minimum_greater_equal",
    "symmetric_within",
)


#: The closed refusal-cause vocabulary. Every entry was earned by an
#: observed loss class; a refusal outside the set raises without a
#: cause rather than guessing one, and the set grows only when a new
#: class has taught us its shape.
TOOLCHAIN_REFUSAL_CAUSES = frozenset(
    {
        "expression_shape",
        "unregistered_constant",
        "expression_read_order",
        "extraction_unit_dimension",
        "expression_unit_dimension",
        "required_output_has_no_producer",
        "claim_output_unrendered",
    }
)


#: The invariant each typed cause protects, stated once. A refusal that
#: types its cause is a failure report (gate, invariant, diagnosis,
#: route, cost) and rides the tool rejection as ``failure_report``; in
#: NOVEL-3 the reports rendered 0 of 31 times because the causes typed
#: here never reached that field (2026-09-05).
TOOLCHAIN_REFUSAL_INVARIANTS: Mapping[str, str] = {
    "expression_shape": (
        "an expression node takes the input count and shape its "
        "operation's arity table admits."
    ),
    "unregistered_constant": (
        "a literature constant is selected by a registered name the host "
        "resolves; a value the model types is a literal, never a constant."
    ),
    "expression_read_order": (
        "expression nodes evaluate in the order given; a node reads only "
        "what an earlier node or an analysis input provides."
    ),
    "extraction_unit_dimension": (
        "an extracted selector carries the dimension its reader declares; "
        "a unit of another dimension names a different quantity."
    ),
    "expression_unit_dimension": (
        "arithmetic across the analysis DAG is dimension-checked; an "
        "operation combines only what its dimension rules admit."
    ),
    "required_output_has_no_producer": (
        "every required output id is produced by some node of the plan, or "
        "declared blocked by a node naming the producer it would need."
    ),
    "claim_output_unrendered": (
        "a claim node renders only what one of its inputs carries; a claim "
        "output no input carries would be a number the host never computed."
    ),
}


class ScientificToolchainContractError(ContractError):
    """Raised when a proposed paper-level tool chain is inconsistent.

    A refusal may carry a typed cause from the closed vocabulary and a
    next-legal-route directive naming what would be admissible instead.
    Both surface as structured fields on the tool rejection: a refusal
    message is an affordance, and the steering belongs in a field the
    session can act on, not only in prose. A route without a cause is
    refused -- the class is what makes the steering auditable.
    """

    def __init__(
        self,
        message: str,
        *,
        cause: str = "",
        next_legal_route: str = "",
    ) -> None:
        if cause and cause not in TOOLCHAIN_REFUSAL_CAUSES:
            raise ContractError(
                f"unknown toolchain refusal cause {cause!r}; the closed "
                f"set is {sorted(TOOLCHAIN_REFUSAL_CAUSES)}"
            )
        if next_legal_route and not cause:
            raise ContractError("a routed refusal names its cause")
        super().__init__(message)
        self.cause = cause
        self.next_legal_route = next_legal_route
        if cause:
            self.failure_report = {
                "gate": cause,
                "invariant": TOOLCHAIN_REFUSAL_INVARIANTS[cause],
                "diagnosis": message,
                "route": next_legal_route
                or "no route is named; the diagnosis says what the plan "
                "must change",
                "cost": "no engine call",
            }
        else:
            # Sixty-six admission refusals type no cause; a session still
            # meets them as failure reports, under one general gate, so
            # the stream counts them and the reply carries the shape
            # (REACH-1: 21 of 26 refusals were bare).
            self.failure_report = {
                "gate": "plan.admission",
                "invariant": (
                    "a scientific workflow is admitted only when every "
                    "node, edge, output and unit resolves against the "
                    "host's typed vocabulary."
                ),
                "diagnosis": message,
                "route": "change the plan as the diagnosis says.",
                "cost": "no engine call",
            }


def _identifier(value: str, field: str) -> str:
    try:
        return require_identifier(value, field)
    except ContractError as exc:
        raise ScientificToolchainContractError(str(exc)) from exc


def _expression_references(item: Mapping[str, Any]) -> tuple[str, ...]:
    """Every name one expression node reads: its inputs and its reference."""

    referenced = [
        str(value) for value in (item.get("input_ids") or ()) if str(value)
    ]
    reference = str(item.get("reference", "") or "")
    if reference:
        referenced.append(reference)
    return tuple(referenced)


def _order_expression_nodes(
    nodes: Sequence[Mapping[str, Any]], *, available: set[str]
) -> tuple[Mapping[str, Any], ...]:
    """Order expression nodes so each follows everything it reads.

    The given order is kept wherever it already satisfies that, so a
    plan written in evaluation order is returned unchanged. A node whose
    reads can never be satisfied is refused with the name it reads; a set
    of nodes that read each other is refused as a cycle.
    """

    provided = set(available)
    pending = list(nodes)
    ordered: list[Mapping[str, Any]] = []
    while pending:
        progressed = False
        for item in list(pending):
            if all(name in provided for name in _expression_references(item)):
                ordered.append(item)
                provided.add(str(item.get("node_id", "?")))
                pending.remove(item)
                progressed = True
        if progressed:
            continue
        pending_ids = {str(item.get("node_id", "?")) for item in pending}
        unprovided = sorted(
            {
                name
                for item in pending
                for name in _expression_references(item)
                if name not in provided and name not in pending_ids
            }
        )
        if unprovided:
            readers = sorted(
                str(item.get("node_id", "?"))
                for item in pending
                if any(
                    name in unprovided for name in _expression_references(item)
                )
            )
            raise ScientificToolchainContractError(
                f"expression node(s) {readers} read {unprovided}, which no "
                "expression node or analysis input provides; "
                f"available: {sorted(provided)}",
                cause="expression_read_order",
                next_legal_route=(
                    "add the missing name as an analysis input or as an "
                    "expression node that computes it"
                ),
            )
        raise ScientificToolchainContractError(
            f"expression nodes {sorted(pending_ids)} read each other in a "
            "cycle, so no evaluation order exists",
            cause="expression_read_order",
            next_legal_route=(
                "break the cycle: each value must be computed from analysis "
                "inputs and nodes that do not read it back"
            ),
        )
    return tuple(ordered)


def _sorted_unique(values: Sequence[str], field: str) -> tuple[str, ...]:
    normalized = tuple(sorted(set(values)))
    if tuple(values) != normalized:
        raise ScientificToolchainContractError(
            f"{field} must be sorted and unique"
        )
    return normalized


@dataclass(frozen=True)
class AnalysisInputIntentV1:
    """A future program output or analysis quantity required by an analysis."""

    input_id: str
    source_kind: str
    producer_node_id: str
    producer_output_id: str
    #: For a claim, the analysis output that computes this number's
    #: uncertainty. An estimate cannot be known before the results
    #: exist; its *estimator* can, so a plan may name the expression
    #: that will produce the spread and the host evaluates it after
    #: execution like any other. The claim it renders is `measured` and
    #: cites that output by id, which is host-checkable by
    #: construction -- the number is the host's own. Without it an
    #: executed delivery arrives `unstated` and costs a further cycle
    #: to assess (owner ruling, 2026-09-09).
    uncertainty_producer_node_id: str = ""
    uncertainty_producer_output_id: str = ""

    def __post_init__(self) -> None:
        _identifier(self.input_id, "analysis input_id")
        if self.source_kind not in {"analysis_output", "program_output"}:
            raise ScientificToolchainContractError(
                "analysis source_kind must be analysis_output or program_output"
            )
        _identifier(self.producer_node_id, "analysis producer_node_id")
        _identifier(self.producer_output_id, "analysis producer_output_id")
        named = bool(self.uncertainty_producer_node_id) + bool(
            self.uncertainty_producer_output_id
        )
        if named == 1:
            raise ScientificToolchainContractError(
                "a planned uncertainty names both its producer node and "
                "its output, or neither"
            )
        if named:
            _identifier(
                self.uncertainty_producer_node_id,
                "uncertainty_producer_node_id",
            )
            _identifier(
                self.uncertainty_producer_output_id,
                "uncertainty_producer_output_id",
            )


@dataclass(frozen=True)
class RegisteredResultInputIntentV1:
    """A host-registered result consumed directly by an analysis root."""

    input_id: str
    artifact_id: str
    source_kind: str = "registered_result"

    def __post_init__(self) -> None:
        _identifier(self.input_id, "analysis input_id")
        _identifier(self.artifact_id, "registered result artifact_id")
        if self.source_kind != "registered_result":
            raise ScientificToolchainContractError(
                "registered result source_kind must be registered_result"
            )


@dataclass(frozen=True)
class AnalysisSelectorIntentV1:
    """Program-neutral quantity requested from a future result artifact."""

    quantity_id: str
    selector: str

    def __post_init__(self) -> None:
        _identifier(self.quantity_id, "analysis selector quantity_id")
        _identifier(self.selector, "analysis selector")


@dataclass(frozen=True)
class AnalysisOutputIntentV1:
    """A named quantity promised by a planned analysis node."""

    output_id: str
    quantity_kind: str
    unit: str

    def __post_init__(self) -> None:
        _identifier(self.output_id, "analysis output_id")
        _identifier(self.quantity_kind, "analysis quantity_kind")


@dataclass(frozen=True)
class AnalysisValidationRuleIntentV1:
    """One program-neutral numerical or categorical validation predicate."""

    rule_id: str
    predicate: str
    input_ids: tuple[str, ...]
    threshold: float | None = None
    expected_count: int | None = None
    unit: str = ""

    def __post_init__(self) -> None:
        _identifier(self.rule_id, "analysis validation rule_id")
        if self.predicate not in ANALYSIS_VALIDATION_PREDICATES:
            raise ScientificToolchainContractError(
                "unsupported analysis validation predicate"
            )
        object.__setattr__(self, "input_ids", tuple(self.input_ids))
        _sorted_unique(self.input_ids, "analysis validation input IDs")
        if not self.input_ids:
            raise ScientificToolchainContractError(
                "analysis validation rule requires a typed input"
            )
        if self.predicate == "count_equals":
            if self.expected_count is None or self.expected_count < 0:
                raise ScientificToolchainContractError(
                    "count_equals requires a non-negative expected_count"
                )
            if self.threshold is not None or self.unit:
                raise ScientificToolchainContractError(
                    "count_equals does not accept threshold or unit"
                )
        elif self.predicate == "integer_equals":
            # A count is non-negative by definition; an integer state label
            # is not -- an anion's charge of -1 was uncertifiable because
            # the only integer predicate refused negatives.
            if self.expected_count is None:
                raise ScientificToolchainContractError(
                    "integer_equals requires expected_count as the "
                    "expected integer"
                )
            if self.threshold is not None or self.unit:
                raise ScientificToolchainContractError(
                    "integer_equals does not accept threshold or unit"
                )
        elif self.predicate in {
            "maximum_absolute_less_equal",
            "minimum_greater_equal",
            "symmetric_within",
        }:
            if self.threshold is None or not str(self.unit).strip():
                raise ScientificToolchainContractError(
                    f"{self.predicate} requires threshold and unit"
                )
            if self.expected_count is not None:
                raise ScientificToolchainContractError(
                    f"{self.predicate} does not accept expected_count"
                )
        elif (
            self.threshold is not None
            or self.expected_count is not None
            or self.unit
        ):
            raise ScientificToolchainContractError(
                f"{self.predicate} accepts only input_ids"
            )


@dataclass(frozen=True)
class AnalysisNodeIntentV1:
    """One extraction, validation, mathematics, or reporting intent."""

    node_id: str
    analysis_kind: str
    dependencies: tuple[str, ...]
    inputs: tuple[AnalysisInputIntentV1 | RegisteredResultInputIntentV1, ...]
    selectors: tuple[AnalysisSelectorIntentV1, ...]
    outputs: tuple[AnalysisOutputIntentV1, ...]
    expression_nodes: tuple[dict[str, object], ...]
    expression_output_node_ids: tuple[str, ...]
    temperature_k: float | None
    pressure_atm: float | None
    support_state: str
    blocked_reason: str
    validation_rules: tuple[AnalysisValidationRuleIntentV1, ...] = ()
    concentration_mol_l: float | None = None
    entropy_method: str = "rrho"
    entropy_cutoff_cm1: float | None = None
    enthalpy_cutoff_cm1: float | None = None
    alpha: int = 4
    use_weighted_mass: bool = False
    frequency_scale_factor: float = 1.0
    #: Which printed mode the session names as the reaction coordinate,
    #: 1-based. Absent from this schema entirely, so the selection could
    #: not travel an approved DAG: it was reachable only from a direct
    #: analysis call, and the executor had nothing to forward.
    reaction_coordinate_mode: int = 0

    def __post_init__(self) -> None:
        _identifier(self.node_id, "analysis node_id")
        if int(self.reaction_coordinate_mode or 0) < 0:
            raise ContractError(
                "reaction_coordinate_mode is a 1-based mode index"
            )
        if self.analysis_kind not in ANALYSIS_INTENT_KINDS:
            raise ScientificToolchainContractError(
                "unsupported analysis intent kind"
            )
        if self.support_state not in ANALYSIS_INTENT_SUPPORT_STATES:
            raise ScientificToolchainContractError(
                "unsupported analysis intent support_state"
            )
        object.__setattr__(self, "dependencies", tuple(self.dependencies))
        object.__setattr__(self, "inputs", tuple(self.inputs))
        object.__setattr__(self, "selectors", tuple(self.selectors))
        object.__setattr__(self, "outputs", tuple(self.outputs))
        object.__setattr__(
            self, "validation_rules", tuple(self.validation_rules)
        )
        object.__setattr__(
            self,
            "expression_nodes",
            tuple(dict(item) for item in self.expression_nodes),
        )
        object.__setattr__(
            self,
            "expression_output_node_ids",
            tuple(self.expression_output_node_ids),
        )
        _sorted_unique(self.dependencies, "analysis dependencies")
        _sorted_unique(
            tuple(item.input_id for item in self.inputs),
            "analysis input IDs",
        )
        _sorted_unique(
            tuple(item.quantity_id for item in self.selectors),
            "analysis selector IDs",
        )
        _sorted_unique(
            tuple(item.output_id for item in self.outputs),
            "analysis output IDs",
        )
        _sorted_unique(
            tuple(item.rule_id for item in self.validation_rules),
            "analysis validation rule IDs",
        )
        if not self.outputs:
            raise ScientificToolchainContractError(
                "analysis intent must preserve at least one output"
            )
        if (
            self.analysis_kind == "claim_rendering"
            and self.support_state == "planned"
            and not self.inputs
        ):
            # The executor renders claims from this node's inputs and from
            # nothing else, so a claim node that binds none can never render
            # one -- determinable from the node alone, with no execution.
            # Found live: such a node previewed green, both engines ran, and
            # it then failed with "claims has 0 items". Dependencies alone do
            # not feed it; each claimed value must be bound as an input
            # naming its producer node and output.
            raise ScientificToolchainContractError(
                f"claim rendering node {self.node_id!r} binds no inputs, so "
                "it can never render a claim; bind each claimed value as an "
                "input naming its producer node and output"
            )
        if (
            self.analysis_kind == "claim_rendering"
            and self.support_state == "planned"
        ):
            # The executor names each claim after the input that binds
            # it, and a claim node's outputs are the claims it renders,
            # so an output id no input carries is a claim nothing will
            # ever render under that name. A live plan declared the six
            # declared observable ids as this node's outputs and bound
            # them under six short input labels; the plan-time gate
            # accepted the outputs, the executor rendered the labels,
            # and the completion gate found none of the six (NOVEL-2
            # po2, 2026-09-04). Provable from the node alone.
            input_ids = {item.input_id for item in self.inputs}
            output_ids = {item.output_id for item in self.outputs}
            unrendered = sorted(output_ids.difference(input_ids))
            if unrendered:
                raise ScientificToolchainContractError(
                    f"claim rendering node {self.node_id!r} declares "
                    f"output(s) {unrendered} that no input of the node "
                    "carries; a claim is rendered under its input_id, "
                    "so name each input by the id the claim must carry "
                    "-- for a declared observable, its observable_id "
                    "verbatim -- and declare that same id as the output",
                    cause="claim_output_unrendered",
                    next_legal_route=(
                        "rename the claim node's input_id to the id the "
                        "claim must carry and keep the output_id equal to "
                        "it; no engine call"
                    ),
                )
        if self.analysis_kind == "result_extraction":
            # The executor already refuses an extraction output that names no
            # selector, and it is right to: with more than one selector the
            # mapping from outputs to quantities is genuinely ambiguous. But
            # it refused *after* the workflow was approved and the engine had
            # finished. Observed four times in one campaign, and the cost is
            # the whole point -- one session ran an optimisation, a minimum
            # check, RRHO thermochemistry at 298.15 K and a kcal/mol
            # conversion, and then lost every number to an output it had
            # called `homo-e` instead of `homo`. The run exited 0 and reported
            # nothing.
            #
            # Both facts the rule needs are here at planning time: the
            # declared outputs, and the selectors on this same node. So the
            # rule belongs here, where the session can still fix it, and the
            # message names the selectors it could have used.
            selector_ids = sorted(
                {selector.quantity_id for selector in self.selectors}
            )
            if len(selector_ids) > 1:
                unmatched = [
                    output.output_id
                    for output in self.outputs
                    if output.output_id not in set(selector_ids)
                ]
                if unmatched:
                    raise ScientificToolchainContractError(
                        f"extraction output(s) {unmatched} name no selector "
                        f"quantity on node {self.node_id!r}; with "
                        f"{len(selector_ids)} selectors the mapping is "
                        "ambiguous, so name each output after the selector "
                        f"it carries -- available here: {selector_ids}"
                    )
        else:
            for output in self.outputs:
                if not str(output.unit).strip():
                    raise ScientificToolchainContractError(
                        f"analysis output {output.output_id!r} "
                        f"({output.quantity_kind}) declares no unit; "
                        + DIMENSIONLESS_UNIT_HINT
                    )
                if str(output.unit).strip().lower() == "count":
                    raise ScientificToolchainContractError(
                        f"analysis output {output.output_id!r} "
                        f"({output.quantity_kind}) declares unsupported unit "
                        f"{output.unit!r}; counts use physical unit '1'"
                    )
        if self.support_state == "blocked_unsupported":
            if not str(self.blocked_reason).strip():
                raise ScientificToolchainContractError(
                    "blocked analysis intent requires a localized reason"
                )
        elif str(self.blocked_reason).strip():
            raise ScientificToolchainContractError(
                "blocked_reason applies only to blocked_unsupported"
            )
        if self.analysis_kind == "result_extraction":
            if len(self.inputs) != 1 or not self.selectors:
                raise ScientificToolchainContractError(
                    "result extraction needs one result input and selectors"
                )
        elif self.analysis_kind == "thermochemistry":
            if self.support_state == "planned" and len(self.inputs) != 1:
                raise ScientificToolchainContractError(
                    "planned thermochemistry needs exactly one result input"
                )
            if self.selectors:
                raise ScientificToolchainContractError(
                    "selectors apply only to result extraction"
                )
            # The RRHO engine writes a fixed vocabulary of quantity IDs.  A
            # node that names an output outside it is asking for something its
            # own receipt will never carry, and nothing downstream can bind:
            # the node stays unmatched and starves every dependent expression,
            # validation and claim, so the failure would surface far away from
            # its cause.  Refuse it here and name what is available, the way
            # the extraction plane answers an unknown selector.
            if self.support_state == "planned":
                available = derivable_thermochemistry_quantities(
                    self.entropy_method
                )
                for output in self.outputs:
                    if (
                        canonical_thermochemistry_quantity(
                            output.quantity_kind
                        )
                        in available
                    ):
                        continue
                    raise ScientificToolchainContractError(
                        f"analysis output {output.output_id!r} declares "
                        f"quantity_kind {output.quantity_kind!r}, which "
                        f"thermochemistry does not derive; it derives "
                        f"{list(available)}"
                    )
        elif any(
            isinstance(item, RegisteredResultInputIntentV1)
            for item in self.inputs
        ):
            raise ScientificToolchainContractError(
                "registered result inputs apply only to result extraction "
                "or thermochemistry"
            )
        elif self.selectors:
            raise ScientificToolchainContractError(
                "selectors apply only to result extraction"
            )
        if (
            self.analysis_kind != "scientific_validation"
            and self.validation_rules
        ):
            raise ScientificToolchainContractError(
                "validation_rules apply only to scientific_validation"
            )
        if (
            self.analysis_kind == "scientific_validation"
            and self.support_state == "planned"
        ):
            if not self.validation_rules:
                raise ScientificToolchainContractError(
                    "scientific validation requires declared rules"
                )
            if len(self.outputs) != 1:
                raise ScientificToolchainContractError(
                    "scientific validation requires one verdict output"
                )
            try:
                _value, _unit, dimension = normalize_numeric_value(
                    0.0, self.outputs[0].unit
                )
            except (QuantityExpressionError, ValueError) as exc:
                raise ScientificToolchainContractError(
                    "scientific validation verdict unit is invalid"
                ) from exc
            if tuple(dimension) != DIMENSIONLESS:
                raise ScientificToolchainContractError(
                    "scientific validation verdict must be dimensionless"
                )
        input_ids = {item.input_id for item in self.inputs}
        for rule in self.validation_rules:
            unknown = set(rule.input_ids).difference(input_ids)
            if unknown:
                raise ScientificToolchainContractError(
                    f"validation rule {rule.rule_id!r} references unknown inputs "
                    f"{sorted(unknown)}"
                )
        if self.analysis_kind == "thermochemistry":
            if (
                self.temperature_k is None
                or not math.isfinite(float(self.temperature_k))
                or self.temperature_k <= 0
            ):
                raise ScientificToolchainContractError(
                    "thermochemistry requires a positive temperature"
                )
            if (
                self.pressure_atm is None
                or not math.isfinite(float(self.pressure_atm))
                or self.pressure_atm <= 0
            ):
                raise ScientificToolchainContractError(
                    "thermochemistry requires a positive pressure"
                )
            method = str(self.entropy_method).strip().lower()
            object.__setattr__(self, "entropy_method", method)
            if method not in {"rrho", "grimme", "truhlar"}:
                raise ScientificToolchainContractError(
                    "thermochemistry entropy_method must be rrho, grimme, "
                    "or truhlar"
                )
            for field in (
                "concentration_mol_l",
                "entropy_cutoff_cm1",
                "enthalpy_cutoff_cm1",
            ):
                value = getattr(self, field)
                if value is not None and (
                    not math.isfinite(float(value)) or float(value) <= 0.0
                ):
                    raise ScientificToolchainContractError(
                        f"thermochemistry {field} must be finite and positive"
                    )
            if method == "rrho" and self.entropy_cutoff_cm1 is not None:
                raise ScientificToolchainContractError(
                    "thermochemistry entropy_cutoff_cm1 requires grimme or "
                    "truhlar entropy"
                )
            if method in {"grimme", "truhlar"} and (
                self.entropy_cutoff_cm1 is None
            ):
                raise ScientificToolchainContractError(
                    f"thermochemistry {method} entropy requires "
                    "entropy_cutoff_cm1"
                )
            if (
                isinstance(self.alpha, bool)
                or not isinstance(self.alpha, int)
                or self.alpha <= 0
            ):
                raise ScientificToolchainContractError(
                    "thermochemistry alpha must be a positive integer"
                )
            if not isinstance(self.use_weighted_mass, bool):
                raise ScientificToolchainContractError(
                    "thermochemistry use_weighted_mass must be boolean"
                )
            if (
                not math.isfinite(float(self.frequency_scale_factor))
                or float(self.frequency_scale_factor) <= 0.0
            ):
                raise ScientificToolchainContractError(
                    "thermochemistry frequency_scale_factor must be finite "
                    "and positive"
                )
        elif self.temperature_k is not None or self.pressure_atm is not None:
            raise ScientificToolchainContractError(
                "temperature and pressure apply only to thermochemistry"
            )
        elif (
            self.concentration_mol_l is not None
            or self.entropy_method != "rrho"
            or self.entropy_cutoff_cm1 is not None
            or self.enthalpy_cutoff_cm1 is not None
            or self.alpha != 4
            or self.use_weighted_mass is not False
            or self.frequency_scale_factor != 1.0
        ):
            raise ScientificToolchainContractError(
                "thermochemistry controls apply only to thermochemistry"
            )
        if self.analysis_kind == "quantity_expression":
            if (
                not self.expression_nodes
                or not self.expression_output_node_ids
            ):
                raise ScientificToolchainContractError(
                    "quantity expression requires a typed expression DAG"
                )
            # Constructing the typed node here makes the evaluator's own
            # shape contract -- known operation, admissible input count,
            # per-operation fields -- a planning refusal. The gap was found
            # live: a six-species pKa workflow wrote an isodesmic exchange
            # as one four-input `subtract`, passed planning and approval,
            # ran six solvated opt+freq jobs, and lost every number to
            # "subtract requires two inputs" -- pure arithmetic on the plan
            # text, checked only after the engines had finished.
            from chemsmart.analysis.quantity_expressions import (
                expression_node_from_plan,
            )
            from chemsmart.analysis.result_quantities import (
                QuantityContractError,
            )

            for item in self.expression_nodes:
                try:
                    expression_node_from_plan(item)
                except (
                    QuantityContractError,
                    QuantityExpressionError,
                ) as exc:
                    raise ScientificToolchainContractError(
                        f"expression node "
                        f"{str(item.get('node_id', '?'))!r}: {exc}",
                        cause="expression_shape",
                        next_legal_route=(
                            "restate the arithmetic as nodes each taking "
                            "its operation's own input count -- e.g. a "
                            "many-term sum as chained two-input add nodes "
                            "-- using only operations from the typed "
                            "vocabulary"
                        ),
                    ) from exc
            # A constant node names a host-owned registry entry; resolving
            # the name here makes an unregistered constant a plan-time
            # refusal that spells out the available set, instead of an
            # execution failure after the engines have already run.
            from chemsmart.analysis.literature_constants import (
                UnknownLiteratureConstantError,
                literature_constant,
            )

            for item in self.expression_nodes:
                if str(item.get("operation", "")) != "constant":
                    continue
                try:
                    literature_constant(str(item.get("constant_name", "")))
                except UnknownLiteratureConstantError as exc:
                    raise ScientificToolchainContractError(
                        f"expression node "
                        f"{str(item.get('node_id', '?'))!r}: {exc.args[0]}",
                        cause="unregistered_constant",
                        next_legal_route=(
                            "select a registered constant from the set "
                            "the message lists, or compute the value from "
                            "delivered quantities; a literal stays "
                            "recorded as model-authored"
                        ),
                    ) from exc
            # An expression node's inputs are resolved in list order while
            # the run evaluates, so a node placed before the node it reads
            # used to raise "references an unavailable prior value" -- and
            # only then, once every engine had finished. The order is a
            # structural fact about the DAG, fully determinable now, and
            # the host can settle it itself: the nodes are put into an
            # order in which every node follows what it reads, and the
            # refusal is kept only for what no order can fix -- a name
            # that nothing provides, or nodes that read each other.
            object.__setattr__(
                self,
                "expression_nodes",
                _order_expression_nodes(
                    self.expression_nodes,
                    available={item.input_id for item in self.inputs},
                ),
            )
            declared = sorted({item.output_id for item in self.outputs})
            exported = sorted(set(self.expression_output_node_ids))
            if declared != exported:
                raise ScientificToolchainContractError(
                    f"expression outputs {declared} on node "
                    f"{self.node_id!r} do not name the exported expression "
                    f"nodes {exported}; the execution receipt keys every "
                    "produced quantity by its expression node id, so a "
                    "claim citing a differently named output cannot "
                    "resolve after the engines have run -- name each "
                    "output after the expression node it exports"
                )
        elif self.expression_nodes or self.expression_output_node_ids:
            raise ScientificToolchainContractError(
                "expression fields apply only to quantity_expression"
            )


@dataclass(frozen=True)
class ScientificToolchainPlanV1:
    """Connected calculation and postprocessing plan proposed in one tool call."""

    schema_version: str
    plan_id: str
    workflow_id: str
    command_workflow_draft_sha256: str
    calculation_node_ids: tuple[str, ...]
    calculation_observables: tuple[tuple[str, tuple[str, ...]], ...]
    analysis_nodes: tuple[AnalysisNodeIntentV1, ...]
    required_output_ids: tuple[str, ...]
    node_order: tuple[str, ...]
    status: str
    plan_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.scientific-toolchain-plan.v1":
            raise ScientificToolchainContractError(
                "unsupported scientific toolchain schema"
            )
        _identifier(self.plan_id, "toolchain plan_id")
        _identifier(self.workflow_id, "toolchain workflow_id")
        if len(self.command_workflow_draft_sha256) != 64:
            raise ScientificToolchainContractError(
                "toolchain requires a command workflow draft digest"
            )
        _sorted_unique(self.calculation_node_ids, "calculation node IDs")
        _sorted_unique(self.required_output_ids, "required output IDs")
        if self.status != "planned":
            raise ScientificToolchainContractError(
                "scientific toolchain status must be planned"
            )
        body = scientific_toolchain_plan_body(self)
        if self.plan_sha256 != canonical_sha256(body):
            raise ScientificToolchainContractError(
                "scientific toolchain plan digest mismatch"
            )


def scientific_toolchain_plan_body(
    plan: ScientificToolchainPlanV1,
) -> dict[str, object]:
    return {
        "schema_version": plan.schema_version,
        "plan_id": plan.plan_id,
        "workflow_id": plan.workflow_id,
        "command_workflow_draft_sha256": plan.command_workflow_draft_sha256,
        "calculation_node_ids": plan.calculation_node_ids,
        "calculation_observables": plan.calculation_observables,
        "analysis_nodes": plan.analysis_nodes,
        "required_output_ids": plan.required_output_ids,
        "node_order": plan.node_order,
        "status": plan.status,
    }


def build_scientific_toolchain_plan(
    *,
    plan_id: str,
    workflow_id: str,
    command_workflow_draft_sha256: str,
    calculation_nodes: Iterable[CommandNodeIntentV1],
    calculation_observables: Mapping[str, Iterable[str]],
    analysis_nodes: Iterable[AnalysisNodeIntentV1],
    required_output_ids: Iterable[str],
) -> ScientificToolchainPlanV1:
    """Validate and topologically order a connected paper-level plan."""

    calculations = tuple(calculation_nodes)
    analyses = tuple(analysis_nodes)
    calculation_by_id = {node.node_id: node for node in calculations}
    analysis_by_id = {node.node_id: node for node in analyses}
    if len(calculation_by_id) != len(calculations):
        raise ScientificToolchainContractError(
            "calculation node IDs must be unique"
        )
    if len(analysis_by_id) != len(analyses):
        raise ScientificToolchainContractError(
            "analysis node IDs must be unique"
        )
    overlap = set(calculation_by_id).intersection(analysis_by_id)
    if overlap:
        raise ScientificToolchainContractError(
            "calculation and analysis node IDs must be distinct"
        )
    all_ids = set(calculation_by_id) | set(analysis_by_id)
    calculation_output_classes: dict[str, dict[str, str]] = {
        node.node_id: {
            item.output_id: item.artifact_class
            for item in node.expected_outputs
        }
        for node in calculations
    }
    produced: dict[str, set[str]] = {
        node_id: set(outputs)
        for node_id, outputs in calculation_output_classes.items()
    }
    produced.update(
        {
            node.node_id: {item.output_id for item in node.outputs}
            for node in analyses
        }
    )
    dependencies: dict[str, set[str]] = {node_id: set() for node_id in all_ids}
    for node in calculations:
        for dependency in node.dependencies:
            if dependency not in calculation_by_id:
                raise ScientificToolchainContractError(
                    "calculation dependency must name a calculation node"
                )
            dependencies[node.node_id].add(dependency)
    for node in analyses:
        # The plan-time half of the declaration gate.  A selector declared
        # for a jobtype is a semantic claim; an extraction that asks a
        # calculation's result for an undeclared selector would either be
        # refused after the engines have run (the late-check class that lost
        # completed work before) or, worse, deliver a value whose meaning
        # for that jobtype was never audited -- an IRC log's starting
        # structure labeled as path endpoints.  Both facts are on hand when
        # the plan is built, so the refusal happens here, naming what IS
        # declared so the session can repair in the same turn.
        if node.analysis_kind != "result_extraction":
            continue
        producers = {
            item.producer_node_id
            for item in node.inputs
            if not isinstance(item, RegisteredResultInputIntentV1)
            and item.producer_node_id in calculation_by_id
        }
        for producer_id in sorted(producers):
            calculation = calculation_by_id[producer_id]
            from chemsmart.analysis.result_readers import reader_for

            reader = reader_for(calculation.program)
            if reader is None or not reader.jobtype_selectors:
                continue
            declared = reader.selectors_for_jobtype(calculation.jobtype)
            if declared is None:
                raise ScientificToolchainContractError(
                    f"extraction node {node.node_id!r} reads "
                    f"{calculation.program!r} jobtype "
                    f"{calculation.jobtype!r}, which declares no selector "
                    "coverage; declared jobtypes: "
                    f"{sorted(j for j, _ in reader.jobtype_selectors)}"
                )
            undeclared = sorted(
                {
                    item.selector
                    for item in node.selectors
                    if item.selector not in declared
                }
            )
            if undeclared:
                from chemsmart.analysis.result_quantities import (
                    thermochemistry_route_hint,
                )

                raise ScientificToolchainContractError(
                    f"extraction node {node.node_id!r} requests selector(s) "
                    f"{undeclared} that are not declared for "
                    f"{calculation.program!r} jobtype "
                    f"{calculation.jobtype!r}; a declaration is a semantic "
                    "claim about what the value means for this job type. "
                    f"Declared here: {sorted(declared)}."
                    + thermochemistry_route_hint(undeclared)
                )
    # The dimensional half of the plan-time gate. Observed live (C5,
    # twice): a relaxed distance scan's coordinate values extract
    # dimensionless -- a scan coordinate may be a length or an angle,
    # so the selector cannot fix a physical dimension -- while the plan
    # declared the output in angstrom; the length-declared chain
    # produced a dimensionless quantity, and the claim render died
    # after every engine had finished, on two independent sessions with
    # correctly ordered expressions. Both halves of the check are
    # derivable from the plan alone: an extraction output's declared
    # unit against the selector's own fixed dimension, and an
    # expression output's declared unit against the chain's derivable
    # result dimension. Unknown rules and unknown dimensions skip --
    # this gate never guesses.
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionError,
        canonical_unit_for_dimension,
        operation_result_dimension,
        unit_dimension,
    )
    from chemsmart.analysis.result_readers import selector_dimension

    for node in analyses:
        if node.analysis_kind != "result_extraction":
            continue
        selector_by_quantity = {
            item.quantity_id: item.selector for item in node.selectors
        }
        for output in node.outputs:
            selector = selector_by_quantity.get(output.output_id)
            if not selector:
                continue
            fixed = selector_dimension(selector)
            if fixed is None:
                continue
            try:
                declared_dimension = unit_dimension(output.unit)
            except QuantityExpressionError:
                continue
            if canonical_unit_for_dimension(
                declared_dimension
            ) != canonical_unit_for_dimension(fixed):
                raise ScientificToolchainContractError(
                    f"extraction output {output.output_id!r} on node "
                    f"{node.node_id!r} declares unit {output.unit!r}, but "
                    f"selector {selector!r} carries dimension "
                    f"{canonical_unit_for_dimension(fixed)!r}. A scan's "
                    "coordinate values, for example, are positional numbers "
                    "in the scan's own coordinate unit; a physical distance "
                    "or angle is measured from a geometry with the "
                    "distance/angle/dihedral operations instead.",
                    cause="extraction_unit_dimension",
                    next_legal_route=(
                        "declare the output in the selector's own "
                        "dimension, or measure the physical coordinate "
                        "from the delivered geometry with the "
                        "distance/angle/dihedral operations"
                    ),
                )

    analysis_output_units = {
        (node.node_id, output.output_id): output.unit
        for node in analyses
        for output in node.outputs
    }
    for node in analyses:
        if node.analysis_kind != "quantity_expression":
            continue
        dimensions: dict[str, tuple | None] = {}
        for item in node.inputs:
            unit = analysis_output_units.get(
                (
                    getattr(item, "producer_node_id", ""),
                    getattr(item, "producer_output_id", ""),
                )
            )
            try:
                dimensions[item.input_id] = (
                    unit_dimension(unit) if unit else None
                )
            except QuantityExpressionError:
                dimensions[item.input_id] = None
        for expression in node.expression_nodes:
            expression_id = str(expression.get("node_id", ""))
            referenced = [
                str(value)
                for value in (expression.get("input_ids") or ())
                if str(value)
            ]
            reference = str(expression.get("reference", "") or "")
            if reference:
                referenced.append(reference)
            try:
                dimensions[expression_id] = operation_result_dimension(
                    str(expression.get("operation", "")),
                    tuple(dimensions.get(name) for name in referenced),
                    literal_unit=str(expression.get("unit", "") or ""),
                    target_unit=str(expression.get("target_unit", "") or ""),
                    constant_name=str(
                        expression.get("constant_name", "") or ""
                    ),
                )
            except QuantityExpressionError as exc:
                raise ScientificToolchainContractError(
                    f"expression node {expression_id!r} on "
                    f"{node.node_id!r}: {exc}",
                    cause="expression_unit_dimension",
                    next_legal_route=(
                        "align the operation's input dimensions with its "
                        "dimension rule, or route the value through an "
                        "explicit conversion or domain-conversion "
                        "operation"
                    ),
                ) from exc
        for output in node.outputs:
            computed = dimensions.get(output.output_id)
            if computed is None:
                continue
            try:
                declared_dimension = unit_dimension(output.unit)
            except QuantityExpressionError:
                continue
            if canonical_unit_for_dimension(
                declared_dimension
            ) != canonical_unit_for_dimension(computed):
                raise ScientificToolchainContractError(
                    f"expression output {output.output_id!r} on node "
                    f"{node.node_id!r} declares unit {output.unit!r} but "
                    "the expression chain produces a quantity of dimension "
                    f"{canonical_unit_for_dimension(computed)!r}. The "
                    "declared unit is checkable from the plan alone; the "
                    "alternative is a claim render failing after every "
                    "engine has finished.",
                    cause="expression_unit_dimension",
                    next_legal_route=(
                        "declare the output in a unit of the computed "
                        "dimension the message names, or change the chain "
                        "so it deliberately produces the declared "
                        "dimension through an explicit conversion"
                    ),
                )
    normalized_analyses: list[AnalysisNodeIntentV1] = []
    for node in analyses:
        effective_dependencies = set(node.dependencies)
        requires_result_artifact = (
            node.support_state == "planned"
            and node.analysis_kind in RESULT_ARTIFACT_ANALYSIS_KINDS
        )
        binds_result_artifact = False
        for item in node.inputs:
            if isinstance(item, RegisteredResultInputIntentV1):
                binds_result_artifact = True
                continue
            producer = item.producer_node_id
            if producer not in all_ids:
                raise ScientificToolchainContractError(
                    "analysis input references an unknown producer"
                )
            expected_kind = (
                "program_output"
                if producer in calculation_by_id
                else "analysis_output"
            )
            if item.source_kind != expected_kind:
                producer_class = (
                    "calculation"
                    if producer in calculation_by_id
                    else "analysis"
                )
                raise ScientificToolchainContractError(
                    f"analysis input {item.input_id!r} declares source_kind "
                    f"{item.source_kind!r}, but producer {producer!r} is a "
                    f"{producer_class} node and requires {expected_kind!r}"
                )
            if item.producer_output_id not in produced[producer]:
                raise ScientificToolchainContractError(
                    "analysis input references an unknown producer output"
                )
            if producer in calculation_by_id and requires_result_artifact:
                calculation = calculation_by_id[producer]
                try:
                    required_artifact_class = typed_result_artifact_kind(
                        calculation.program
                    )
                except ContractError as exc:
                    raise ScientificToolchainContractError(str(exc)) from exc
                declared_artifact_class = calculation_output_classes[producer][
                    item.producer_output_id
                ]
                if declared_artifact_class != required_artifact_class:
                    raise ScientificToolchainContractError(
                        f"analysis node {node.node_id!r} "
                        f"({node.analysis_kind}) consumes "
                        f"{producer}.{item.producer_output_id}, declared as "
                        f"{declared_artifact_class!r}; typed "
                        f"{calculation.program} analysis requires a "
                        f"{required_artifact_class!r} calculation output. "
                        "Declare that result output separately; native "
                        "geometry or Hessian outputs remain available for "
                        "program-to-program handoffs."
                    )
                binds_result_artifact = True
            effective_dependencies.add(producer)
        if requires_result_artifact and not binds_result_artifact:
            raise ScientificToolchainContractError(
                f"analysis node {node.node_id!r} ({node.analysis_kind}) "
                "binds no complete typed program result. Use a registered "
                "result root or consume the immediate calculation's "
                "program_output; an analysis_output quantity cannot replace "
                "the result artifact."
            )
        unknown = effective_dependencies.difference(all_ids)
        if unknown:
            raise ScientificToolchainContractError(
                "analysis dependency references an unknown node"
            )
        effective_dependencies.discard(node.node_id)
        normalized = AnalysisNodeIntentV1(
            node_id=node.node_id,
            analysis_kind=node.analysis_kind,
            dependencies=tuple(sorted(effective_dependencies)),
            inputs=node.inputs,
            selectors=node.selectors,
            outputs=node.outputs,
            expression_nodes=node.expression_nodes,
            expression_output_node_ids=node.expression_output_node_ids,
            temperature_k=node.temperature_k,
            pressure_atm=node.pressure_atm,
            support_state=node.support_state,
            blocked_reason=node.blocked_reason,
            validation_rules=node.validation_rules,
            concentration_mol_l=node.concentration_mol_l,
            entropy_method=node.entropy_method,
            entropy_cutoff_cm1=node.entropy_cutoff_cm1,
            enthalpy_cutoff_cm1=node.enthalpy_cutoff_cm1,
            alpha=node.alpha,
            use_weighted_mass=node.use_weighted_mass,
            frequency_scale_factor=node.frequency_scale_factor,
        )
        normalized_analyses.append(normalized)
        dependencies[node.node_id].update(effective_dependencies)
    children: dict[str, set[str]] = {node_id: set() for node_id in all_ids}
    for node_id, parents in dependencies.items():
        for parent in parents:
            children[parent].add(node_id)
    ready = sorted(
        node_id for node_id, parents in dependencies.items() if not parents
    )
    ordered: list[str] = []
    while ready:
        node_id = ready.pop(0)
        ordered.append(node_id)
        for child in sorted(children[node_id]):
            dependencies[child].discard(node_id)
            if not dependencies[child] and child not in ordered + ready:
                ready.append(child)
                ready.sort()
    if len(ordered) != len(all_ids):
        raise ScientificToolchainContractError(
            "scientific toolchain contains a dependency cycle"
        )
    all_output_ids = {
        output_id
        for output_ids in produced.values()
        for output_id in output_ids
    }
    required = tuple(sorted(set(required_output_ids)))
    missing = sorted(set(required).difference(all_output_ids))
    if missing:
        # A refusal names the route. A live session declared the ground
        # state's multiplicity as a required output, met this refusal
        # bare, and built the label itself from min, typed all_equal
        # verdicts and a one-hot sum -- correctly, and at the cost of a
        # plan (NOVEL-2 ino1, 2026-09-04); coordinate_at_minimum was in
        # the vocabulary and nothing said so.
        raise ScientificToolchainContractError(
            f"required output(s) have no producer: {missing}; every "
            "required output id is an output_id of some node in this "
            "plan",
            cause="required_output_has_no_producer",
            next_legal_route=(
                "declare each missing id as an output of the node that "
                "computes it -- for a declared observable, the claim "
                "node's input_id and output_id, verbatim; for a label "
                "such as which state lies lowest, an expression node "
                "using coordinate_at_minimum over (energies, labels), or "
                "min with all_equal verdicts; for an observable no "
                "producer can deliver, a blocked_unsupported analysis node "
                "whose output_id is that id, with blocked_reason naming the "
                "producer it would need, so the refusal is typed; no "
                "engine call"
            ),
        )
    observable_rows = tuple(
        sorted(
            (
                node_id,
                tuple(sorted(set(calculation_observables.get(node_id, ())))),
            )
            for node_id in calculation_by_id
        )
    )
    analysis_by_id = {node.node_id: node for node in normalized_analyses}
    body = {
        "schema_version": "chemsmart.scientific-toolchain-plan.v1",
        "plan_id": plan_id,
        "workflow_id": workflow_id,
        "command_workflow_draft_sha256": command_workflow_draft_sha256,
        "calculation_node_ids": tuple(sorted(calculation_by_id)),
        "calculation_observables": observable_rows,
        "analysis_nodes": tuple(
            analysis_by_id[node_id]
            for node_id in ordered
            if node_id in analysis_by_id
        ),
        "required_output_ids": required,
        "node_order": tuple(ordered),
        "status": "planned",
    }
    return ScientificToolchainPlanV1(
        **body, plan_sha256=canonical_sha256(body)
    )


#: The tool a session calls to execute each analysis kind. The frontier
#: names it on every actionable analysis node, so "actionable" never
#: reads as "awaiting approval"; under a goal wake the host executes an
#: analysis-only plan itself and these nodes arrive completed.
ANALYSIS_NEXT_TOOL: Mapping[str, str] = {
    "result_extraction": "extract_result_quantities",
    "thermochemistry": "derive_thermochemistry",
    "quantity_expression": "evaluate_quantity_expression",
    "scientific_validation": "evaluate_scientific_validation",
    "claim_rendering": "record_analysis_claims",
}


def project_scientific_toolchain_frontier(
    plan: ScientificToolchainPlanV1,
    *,
    actionable_calculation_node_ids: Iterable[str] = (),
    unresolved_calculation_node_ids: Iterable[str] = (),
    completed_calculation_node_ids: Iterable[str] = (),
    completed_analysis_node_ids: Iterable[str] = (),
    non_executable_calculation_node_ids: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Return model-useful next actions without pretending future data exist."""

    actionable = set(actionable_calculation_node_ids)
    unresolved = set(unresolved_calculation_node_ids)
    completed = set(completed_calculation_node_ids)
    completed_analysis = set(completed_analysis_node_ids)
    # A live session read "actionable / compile_and_preview" on a stage that
    # approval readiness was simultaneously excluding as non-executable, and
    # had to flag the contradiction itself. The readiness side already owns
    # the truth (declared blocked_unsupported plus its data-edge cascade);
    # the frontier now states the same truth instead of a second opinion.
    non_executable = dict(non_executable_calculation_node_ids or {})
    states: dict[str, str] = {}
    nodes: list[dict[str, object]] = []
    analysis_by_id = {node.node_id: node for node in plan.analysis_nodes}
    # Which nodes consume each producer.  The model submitted the edges; the
    # host can therefore state the branching instead of leaving the model to
    # recall its own DAG, which is where a stage carrying the answer is lost.
    dependents: dict[str, set[str]] = {}
    for node in plan.analysis_nodes:
        for edge in node.inputs:
            if isinstance(edge, RegisteredResultInputIntentV1):
                continue
            dependents.setdefault(edge.producer_node_id, set()).add(
                node.node_id
            )
    for node_id in plan.node_order:
        waiting_on: tuple[str, ...] = ()
        unsatisfied: tuple[dict[str, str], ...] = ()
        if node_id in plan.calculation_node_ids:
            if node_id in completed:
                state = "completed"
            elif node_id in non_executable:
                state = "non_executable"
            elif node_id in actionable:
                state = "actionable"
            else:
                state = "waiting"
            if node_id in unresolved and state not in {
                "completed",
                "non_executable",
            }:
                state = "waiting"
            if state == "completed":
                reason = "validated execution result already exists"
            elif state == "non_executable":
                reason = non_executable[node_id]
            elif state == "actionable":
                reason = "compile_and_preview"
            else:
                reason = "resolve_inputs_or_scientific_settings"
        else:
            node = analysis_by_id[node_id]
            blocked_parents = tuple(
                sorted(
                    parent
                    for parent in node.dependencies
                    if states.get(parent)
                    in {
                        "blocked_unsupported",
                        "blocked_upstream",
                        "non_executable",
                    }
                )
            )
            if node.support_state == "blocked_unsupported":
                state = "blocked_unsupported"
                reason = node.blocked_reason
            elif node_id in completed_analysis:
                state = "completed"
                reason = "typed analysis evidence already exists"
            elif blocked_parents:
                state = "blocked_upstream"
                waiting_on = blocked_parents
                reason = (
                    "blocked because "
                    + ", ".join(blocked_parents)
                    + (" is" if len(blocked_parents) == 1 else " are")
                    + " blocked"
                )
            elif node.dependencies and not all(
                states.get(parent) == "completed"
                for parent in node.dependencies
            ):
                state = "waiting_for_artifact"
                future_inputs = tuple(
                    edge
                    for edge in node.inputs
                    if isinstance(edge, AnalysisInputIntentV1)
                )
                # Name the exact producer outputs.  "await producer outputs"
                # is the same sentence for every waiting node in every
                # workflow, so it carries none of the dependency structure the
                # host already holds.
                unsatisfied = tuple(
                    {
                        "input_id": edge.input_id,
                        "producer_node_id": edge.producer_node_id,
                        "producer_output_id": edge.producer_output_id,
                        "source_kind": edge.source_kind,
                    }
                    for edge in sorted(
                        future_inputs,
                        key=lambda item: (
                            item.producer_node_id,
                            item.producer_output_id,
                            item.input_id,
                        ),
                    )
                )
                waiting_on = tuple(sorted(set(node.dependencies)))
                named = ", ".join(
                    f"{edge['producer_node_id']}.{edge['producer_output_id']}"
                    for edge in unsatisfied
                ) or ", ".join(waiting_on)
                reason = f"await {named}"
            else:
                state = "actionable"
                # Name who acts. A woken session read "execute the
                # registered analysis operation" on fourteen thermochemistry
                # nodes, took it for a stage awaiting approval, and stopped
                # with every number in hand (NOVEL-2 po2, 2026-09-04). The
                # host knows the tool each kind is executed with.
                reason = (
                    "this session calls "
                    + ANALYSIS_NEXT_TOOL.get(
                        node.analysis_kind, "the registered analysis tool"
                    )
                    + " for this node"
                )
        states[node_id] = state
        entry: dict[str, object] = {
            "node_id": node_id,
            "state": state,
            "next_action": reason,
            "feeds": tuple(sorted(dependents.get(node_id, ()))),
        }
        if waiting_on:
            entry["waiting_on"] = waiting_on
        if unsatisfied:
            entry["unsatisfied_inputs"] = unsatisfied
        if (
            node_id not in plan.calculation_node_ids
            and state == "actionable"
            and analysis_by_id[node_id].analysis_kind in ANALYSIS_NEXT_TOOL
        ):
            entry["next_tool"] = ANALYSIS_NEXT_TOOL[
                analysis_by_id[node_id].analysis_kind
            ]
        nodes.append(entry)
    return {
        "schema_version": "chemsmart.scientific-toolchain-frontier.v1",
        "workflow_id": plan.workflow_id,
        "nodes": tuple(nodes),
        "actionable_node_ids": tuple(
            item["node_id"] for item in nodes if item["state"] == "actionable"
        ),
        "completed_node_ids": tuple(
            item["node_id"] for item in nodes if item["state"] == "completed"
        ),
        "waiting_node_ids": tuple(
            item["node_id"]
            for item in nodes
            if item["state"] in {"waiting", "waiting_for_artifact"}
        ),
        "blocked_node_ids": tuple(
            item["node_id"]
            for item in nodes
            if item["state"] in {"blocked_unsupported", "blocked_upstream"}
        ),
    }


__all__ = [
    "ANALYSIS_INTENT_KINDS",
    "ANALYSIS_INTENT_SUPPORT_STATES",
    "AnalysisInputIntentV1",
    "AnalysisNodeIntentV1",
    "AnalysisOutputIntentV1",
    "RegisteredResultInputIntentV1",
    "AnalysisSelectorIntentV1",
    "AnalysisValidationRuleIntentV1",
    "ANALYSIS_VALIDATION_PREDICATES",
    "ScientificToolchainContractError",
    "ScientificToolchainPlanV1",
    "build_scientific_toolchain_plan",
    "project_scientific_toolchain_frontier",
    "scientific_toolchain_plan_body",
]
