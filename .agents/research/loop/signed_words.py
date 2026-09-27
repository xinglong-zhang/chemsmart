"""The classes of word the CHEMSMART host signs, each with its one signer.

A word the host signs -- a settlement, a verified refusal, a certification,
a stationarity verdict, a finding's standing, a category's answer -- is
computed by one function and read by others. This table names that
function per class, the record the word lands in, and how a census can
check it: RE-SIGNED (a replay or a pure function recomputes it on the
imported tree), READER (an independent reader checks it against the records
it cites; the word was signed inside a session with host state no replay
rebuilds), or CENSUS (counted and classified, not verified).

    PYTHONPATH=<tree> python .agents/research/loop/signed_words.py [--json]

prints the table and imports every named signer from the tree on
``PYTHONPATH``, so a renamed or deleted signer turns the census red rather
than leaving a sentence about a function that no longer exists. Written for
R11 episode `truth`; the enumeration was read from the code at 9185770e.
"""

from __future__ import annotations

import importlib
import json
import sys

#: (id, word, signers "module:qualname", record, method)
CLASSES: tuple[tuple[str, str, tuple[str, ...], str, str], ...] = (
    (
        "W1",
        "settlement: achieved | achieved_with_observations | "
        "unreachable_from_evidence | exhausted | returned_to_human, with reasons",
        (
            "chemsmart.agent.driver:GoalDriver._settle",
            "chemsmart.agent.driver:_delivery_settlement",
            "chemsmart.agent.driver:GoalDriver._settle_after_reading",
            "chemsmart.agent.driver:_typed_error_settlement",
            "chemsmart.agent.driver:_achieved_word",
        ),
        "ledger goal_settled {state, reasons, evidence}",
        "RE-SIGNED (settle replay, final step) + READER",
    ),
    (
        "W2",
        "held words: recovery_opened, reading_opened, rewake_opened",
        ("chemsmart.agent.driver:GoalDriver._settle",),
        "ledger rows",
        "RE-SIGNED where it is the goal's final settle-type row",
    ),
    (
        "W3",
        "qualified capability rows",
        ("chemsmart.agent.driver:_record_goal_qualification",),
        "ledger qualified",
        "RE-SIGNED (settle replay) + READER",
    ),
    (
        "W4",
        "executor analysis word: completed | partial | ''",
        (
            "chemsmart.agent.executor:"
            "ApprovedWorkflowExecutor._run_analysis_phase",
        ),
        "run_recorded.analysis_status; execution-result.json",
        "RE-SIGNED (walk) for chains that ran no node; READER otherwise",
    ),
    (
        "W5",
        "certification: completion passed | partial, limitation, anomaly, "
        "falsified-expectation and failed-criterion ids",
        (
            "chemsmart.agent.tool_runtime:"
            "CommandCompiledToolHostV1._record_toolchain_completion",
        ),
        "analysis_completion_evaluated",
        "READER",
    ),
    (
        "W6",
        "expectation verdict: agreed | diverged | indeterminate | "
        "not_comparable | agreed_as_approximation",
        (
            "chemsmart.agent.tool_runtime:"
            "CommandCompiledToolHostV1._declared_observable_predictions",
        ),
        "completion declared_observable_predictions",
        "READER (arithmetic recomputed from the row)",
    ),
    (
        "W7",
        "verified refusal: verified + basis",
        (
            "chemsmart.agent.tool_runtime:"
            "CommandCompiledToolHostV1._verify_unreachable_observables",
            "chemsmart.agent.tool_runtime:refusal_read_against_results",
            "chemsmart.agent.driver:"
            "GoalDriver._refusals_the_results_now_answer",
        ),
        "scientific_decision_recorded.unreachable_observables",
        "READER; the final cycle's re-read is inside the W1 replay",
    ),
    (
        "W8",
        "finding standing: answers | on_the_request | unrequested; "
        "relation truth",
        (
            "chemsmart.agent.tool_runtime:"
            "CommandCompiledToolHostV1._verify_findings",
            "chemsmart.agent.analysis_claims:evaluate_finding_relation",
        ),
        "scientific_decision_recorded.findings",
        "READER",
    ),
    (
        "W9",
        "category answer: the word the host read",
        (
            "chemsmart.agent.tool_runtime:"
            "CommandCompiledToolHostV1._categorical_answer_row",
        ),
        "finding answer[]; completion declared_categorical_answers",
        "READER",
    ),
    (
        "W10",
        "stationary-point order and stationarity",
        (
            "chemsmart.agent.execution:build_stationary_point_characterisation",
            "chemsmart.analysis.result_quantities:structure_stationarity",
        ),
        "stationary_point_characterised",
        "RE-SIGNED (pure function over the archived result file)",
    ),
    (
        "W11",
        "a free energy stands on a stationary structure of a named surface",
        (
            "chemsmart.analysis.result_quantities:structure_stationarity",
            "chemsmart.analysis.result_quantities:free_energy_surface",
        ),
        "thermochemistry_derived",
        "RE-SIGNED (stationarity of the archived source result) + READER",
    ),
    (
        "W12",
        "node terminal state and result validity",
        ("chemsmart.agent.terminal_states:derive_run_outcome",),
        "program_result_verified; workflow_node_state_changed",
        "READER-lite (stationary-point rule from the verification record)",
    ),
    (
        "W13",
        "anomaly standing: unreplicated | replicated | refuted",
        ("chemsmart.agent.execution:anomaly_standing",),
        "anomaly_observed; ledger anomalies_observed",
        "READER (as observations the first word must carry)",
    ),
    (
        "W14",
        "sufficiency: met | attested | short | unstated",
        ("chemsmart.agent.delivery:judge_sufficiency",),
        "claim requirement assessment rows",
        "READER where the rows carry the operands",
    ),
    (
        "W15",
        "a failed acceptance criterion answered or not",
        ("chemsmart.agent.goal:failed_criteria",),
        "derived into W1 and W5",
        "READER (inside W1 and W5)",
    ),
    (
        "W16",
        "revision admitted | returned",
        ("chemsmart.agent.goal:admit_revision",),
        "ledger revision_admitted / revision_returned",
        "not checked (authority, not truth of evidence)",
    ),
    (
        "W17",
        "refusal: gate, invariant, diagnosis, route",
        ("chemsmart.agent._contracts:RoutedContractError",),
        "tool_failed",
        "CENSUS (grouped by tool, error class and gate)",
    ),
)


def resolve(target: str) -> bool:
    module_name, _, qualname = target.partition(":")
    try:
        obj = importlib.import_module(module_name)
    except Exception:  # noqa: BLE001 - reported, never hidden
        return False
    for part in qualname.split("."):
        obj = getattr(obj, part, None)
        if obj is None:
            return False
    return True


def main() -> int:
    rows = []
    missing = []
    for cid, word, signers, record, method in CLASSES:
        found = {signer: resolve(signer) for signer in signers}
        missing.extend(signer for signer, ok in found.items() if not ok)
        rows.append(
            {
                "id": cid,
                "word": word,
                "signers": found,
                "record": record,
                "method": method,
            }
        )
    if "--json" in sys.argv:
        print(json.dumps(rows, indent=1))
    else:
        import chemsmart

        print("chemsmart imported from", chemsmart.__file__)
        for row in rows:
            marks = ", ".join(
                f"{name.split(':', 1)[1]}{'' if ok else ' (MISSING)'}"
                for name, ok in row["signers"].items()
            )
            print(f"{row['id']:4} {row['method'][:26]:26} {marks}")
    if missing:
        print("missing signers:", ", ".join(missing), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
