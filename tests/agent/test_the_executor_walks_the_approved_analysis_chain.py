"""The approved analysis chain executes provider-free after the engines.

Through four audited benchmark executions the subject's planned analysis
chains never ran: `agent run` stopped at engine-complete, and the
reviewing scientist extracted every number outside the product. The
executor now walks the digest-bound chain the same /approve covered,
dispatching the session's own analysis tools with host-computed arguments
-- extraction, thermochemistry, expressions, validation verdicts, claims
-- then records the toolchain completion receipt and renders the report,
with zero provider calls by construction.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.executor import ApprovedWorkflowExecutor
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.scientific_toolchain import (
    AnalysisInputIntentV1,
    AnalysisNodeIntentV1,
    AnalysisOutputIntentV1,
    AnalysisSelectorIntentV1,
    AnalysisValidationRuleIntentV1,
    build_scientific_toolchain_plan,
)
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.agent.workflows import (
    ArtifactInputIntentV1,
    ArtifactOutputIntentV1,
    CommandNodeIntentV1,
)

_RESULT = Path("tests/data/ORCATests/outputs/CO2.out")


def _calculation():
    return CommandNodeIntentV1(
        node_id="sp",
        program="orca",
        jobtype="sp",
        project_role="r",
        dependencies=(),
        inputs=(
            ArtifactInputIntentV1(
                binding_id="geometry",
                artifact_class="geometry_xyz",
                artifact_id="start",
                producer_node_id="",
                producer_output_id="",
            ),
        ),
        expected_outputs=(
            ArtifactOutputIntentV1(
                output_id="sp-out", artifact_class="orca_output"
            ),
        ),
        unresolved_fields=(),
    )


def _analysis_node(node_id, kind, **overrides):
    base = dict(
        node_id=node_id,
        analysis_kind=kind,
        dependencies=(),
        inputs=(),
        selectors=(),
        outputs=(),
        expression_nodes=(),
        expression_output_node_ids=(),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
        validation_rules=(),
    )
    base.update(overrides)
    return AnalysisNodeIntentV1(**base)


def _chain(selector="energy", validation_threshold=None):
    extraction = _analysis_node(
        "extract-sp",
        "result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector=selector),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    expression = _analysis_node(
        "to-kcal",
        "quantity_expression",
        dependencies=("extract-sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="e-in",
                source_kind="analysis_output",
                producer_node_id="extract-sp",
                producer_output_id="e",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e-kcal",
                quantity_kind="energy",
                unit="kcal/mol",
            ),
        ),
        expression_nodes=(
            {
                "node_id": "e-kcal",
                "operation": "convert",
                "input_ids": ("e-in",),
                "target_unit": "kcal/mol",
            },
        ),
        expression_output_node_ids=("e-kcal",),
    )
    rules = (
        AnalysisValidationRuleIntentV1(
            rule_id="finite",
            predicate="all_finite",
            input_ids=("val-in",),
        ),
    )
    if validation_threshold is not None:
        rules = (
            AnalysisValidationRuleIntentV1(
                rule_id="above-threshold",
                predicate="minimum_greater_equal",
                input_ids=("val-in",),
                threshold=validation_threshold,
                unit="kcal/mol",
            ),
        )
    validation = _analysis_node(
        "check",
        "scientific_validation",
        dependencies=("to-kcal",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="val-in",
                source_kind="analysis_output",
                producer_node_id="to-kcal",
                producer_output_id="e-kcal",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="verdict", quantity_kind="count", unit="1"
            ),
        ),
        validation_rules=rules,
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("to-kcal",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="to-kcal",
                producer_output_id="e-kcal",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="kcal/mol",
            ),
        ),
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, expression, validation, claims),
        required_output_ids=("final-energy",),
    )


def _executor(tmp_path, toolchain):
    resolved = _RESULT.resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="result.sp.1",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="exec"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
        approved_scientific_toolchain_plan=toolchain,
    )
    host.artifacts[artifact.artifact_id] = artifact
    host.execution_receipts["sp"] = SimpleNamespace(validated=True)
    run_directory = tmp_path / "run"
    run_directory.mkdir()
    return ApprovedWorkflowExecutor(
        host=host,
        plan=SimpleNamespace(
            workflow_id="w",
            plan_sha256="b" * 64,
            nodes=(SimpleNamespace(node_id="sp", program="orca"),),
        ),
        approval=SimpleNamespace(node_bindings=()),
        frozen_approval=SimpleNamespace(approval_sha256="c" * 64),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256="a" * 64,
        run_directory=run_directory,
        execution_bundle=SimpleNamespace(non_executable_node_ids=()),
        approval_workspace=tmp_path / "workspace",
        claim_workspace_bundle=False,
    )


def test_the_chain_executes_completes_and_renders(tmp_path):
    toolchain = _chain()
    executor = _executor(tmp_path, toolchain)

    nodes, status, completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    states = {node.node_id: node.state for node in nodes}
    assert states == {
        "extract-sp": "executed",
        "to-kcal": "executed",
        "check": "executed",
        "claims": "executed",
    }
    assert status == "completed"
    assert len(completions) == 1
    report = Path(report_path).read_text(encoding="utf-8")
    assert "Host-rendered numerical claims" in report
    assert "final-energy" in report
    assert "Scientific decision: not recorded" in report

    kinds = [event.kind for event in executor.host.event_store.read_events()]
    assert kinds.count(EventKind.WORKFLOW_ANALYSIS_NODE_SETTLED.value) == 4
    assert EventKind.WORKFLOW_ANALYSIS_REPORT_RENDERED.value in kinds


def test_a_failed_verdict_is_still_a_completed_analysis(tmp_path):
    toolchain = _chain(validation_threshold=1.0e9)
    executor = _executor(tmp_path, toolchain)

    nodes, status, completions, _report = executor._run_analysis_phase(
        toolchain
    )

    check = next(node for node in nodes if node.node_id == "check")
    assert check.state == "executed"
    assert "=0" in check.reason.replace(" ", "")
    assert status == "completed"
    assert completions


def test_a_broken_analysis_input_skips_its_dependents(tmp_path):
    # solvation_electrostatic_energy is declared for orca sp (so the plan
    # gate admits it), carries the energy dimension the chain's convert
    # needs (so the dimensional gate admits it), and is absent at runtime
    # on this gas-phase result -- absence is meaning, and the walk must
    # starve dependents rather than crash.
    toolchain = _chain(selector="solvation_electrostatic_energy")
    executor = _executor(tmp_path, toolchain)

    nodes, status, completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    states = {node.node_id: node.state for node in nodes}
    # The node ran and delivered nothing it was asked for, so it settles
    # executed with the output named absent -- the walk did run it -- while
    # the consumer that named that output still does not run.  Which of the
    # two the state records is the round's whole point: what the walk did,
    # against what the evidence contains.
    assert states["extract-sp"] == "executed"
    assert states["to-kcal"] == "skipped"
    assert states["check"] == "skipped"
    assert status == "partial"
    extraction = next(node for node in nodes if node.node_id == "extract-sp")
    assert extraction.absent_output_ids == ("e",)
    starved = next(node for node in nodes if node.node_id == "to-kcal")
    assert "extract-sp.e" in starved.reason
    # The partial envelope delivers the failure to the reader instead of
    # leaving an empty run directory.
    assert report_path.endswith("partial-analysis-report.md")
    assert len(completions) == 1


def test_blocked_analysis_intent_stays_visible_not_executed(tmp_path):
    extraction = _analysis_node(
        "extract-sp",
        "result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    blocked = _analysis_node(
        "unsupported",
        "unsupported_external",
        support_state="blocked_unsupported",
        blocked_reason="this release has no parser for the companion output",
        outputs=(
            AnalysisOutputIntentV1(
                output_id="companion",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
    )
    toolchain = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, blocked),
        required_output_ids=("e",),
    )
    executor = _executor(tmp_path, toolchain)

    nodes, status, completions, _report = executor._run_analysis_phase(
        toolchain
    )

    states = {node.node_id: node for node in nodes}
    assert states["unsupported"].state == "blocked_unsupported"
    assert (
        states["unsupported"].reason
        == "this release has no parser for the companion output"
    )
    assert states["extract-sp"].state == "executed"
    assert status == "completed"
    assert completions


def _chain_with_selector_and_no_dependents(selector):
    """One extraction plus a dependent claim -- the l3b shape, minimized."""

    extraction = _analysis_node(
        "extract-spin",
        "result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="s2", selector=selector),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="s2", quantity_kind="count", unit="1"
            ),
        ),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("extract-spin",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="s2-final",
                source_kind="analysis_output",
                producer_node_id="extract-spin",
                producer_output_id="s2",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="s2-final", quantity_kind="count", unit="1"
            ),
        ),
    )
    return build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, claims),
        required_output_ids=("s2-final",),
    )


def test_an_absent_quantity_settles_the_node_instead_of_crashing(tmp_path):
    """Six benchmark runs died here: <S^2> asked of a closed-shell result.

    The extraction layer raised its own typed error with an exemplary
    message -- naming the cause and every quantity the result does resolve --
    and the walk, catching only ContractError, let it kill the executor, so
    nobody ever read that message and every other receipt was lost with it.
    """

    toolchain = _chain_with_selector_and_no_dependents("spin_square")
    executor = _executor(tmp_path, toolchain)

    nodes, status, receipts, report_path = executor._run_analysis_phase(
        toolchain
    )

    settled = {record.node_id: record for record in nodes}
    # The refusal is no longer fatal to the node, and the claim that mattered
    # here is unchanged: a consumer that names the refused output does not
    # run, nothing crashes, and the exemplary message still reaches a reader.
    assert settled["extract-spin"].state == "executed"
    assert settled["extract-spin"].absent_output_ids == ("s2",)
    assert settled["claims"].state == "skipped"
    assert "extract-spin.s2" in settled["claims"].reason
    assert status == "partial"
    # The envelope delivers the exemplary refusal to the reader: findings
    # name what did not run and the footer names the recovery act, while no
    # claims heading appears because nothing reached claim standing.
    assert report_path.endswith("partial-analysis-report.md")
    report = Path(report_path).read_text()
    assert "Partial analysis: 1 finding(s) over 2 analysis nodes" in report
    # The exemplary refusal still reaches the reader -- now under its own
    # heading, as a quantity the result does not carry, rather than as a
    # node that did not execute.
    from chemsmart.agent.report_format import ABSENCES_HEADING

    assert ABSENCES_HEADING in report
    assert "no <S^2>" in report
    assert "Recovery:" in report
    from chemsmart.agent.report_format import CLAIMS_HEADING

    assert CLAIMS_HEADING not in report


def test_several_claim_records_complete_green(tmp_path):
    """Two claim-rendering nodes are a whole delivery, not a broken one.

    This shape used to crash, was repaired into a durable refusal, and the
    refusal itself was still wrong: a live chain delivering two reduction
    potentials and their difference recorded three claim records with every
    node executed and every verdict passed, and the binding rule downgraded
    it to "partial" over a finding about arity rather than chemistry -- in
    a report that rendered all three records regardless. A chain that
    breaks is named by its findings; the number of claim stages a healthy
    plan chose says nothing about whether it broke.
    """

    from chemsmart.agent.scientific_toolchain import (
        build_scientific_toolchain_plan as _build,
    )

    base = _chain()
    second_claims = _analysis_node(
        "claims-again",
        "claim_rendering",
        dependencies=("to-kcal",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy-again",
                source_kind="analysis_output",
                producer_node_id="to-kcal",
                producer_output_id="e-kcal",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy-again",
                quantity_kind="energy",
                unit="kcal/mol",
            ),
        ),
    )
    toolchain = _build(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=tuple(base.analysis_nodes) + (second_claims,),
        required_output_ids=("final-energy", "final-energy-again"),
    )
    executor = _executor(tmp_path, toolchain)

    nodes, status, receipts, report_path = executor._run_analysis_phase(
        toolchain
    )

    assert all(
        record.state in {"executed", "blocked_unsupported"} for record in nodes
    )
    assert status == "completed"
    refusals = [
        event
        for event in executor.host.event_store.read_events()
        if event.kind == EventKind.WORKFLOW_ANALYSIS_COMPLETION_REFUSED.value
    ]
    assert refusals == []
    # Both validated claim records render, each under its own label, in a
    # report that no longer calls a whole delivery partial.
    assert report_path.endswith("completed-analysis-report.md")
    report = Path(report_path).read_text()
    assert "at most one claim record" not in report
    assert "Partial analysis" not in report
    assert report.count("Claim record:") == 2
    assert len(receipts) == 1


def test_a_successful_report_shows_its_validation_verdicts(tmp_path):
    """Verdicts used to survive only as prose on settled events; a reader
    never saw them even when every rule passed. The report carries them."""

    toolchain = _chain()
    executor = _executor(tmp_path, toolchain)

    nodes, status, receipts, report_path = executor._run_analysis_phase(
        toolchain
    )

    assert status == "completed"
    report = Path(report_path).read_text()
    from chemsmart.agent.report_format import VERDICTS_HEADING

    assert VERDICTS_HEADING in report
    assert "all_finite" in report
    assert "| passed |" in report


def test_a_chain_using_a_constant_renders_its_provenance_table(tmp_path):
    # The expression applies the 1 atm -> 1 mol/L standard-state
    # correction as a registry constant: the value is host-resolved, and
    # the report and its terminal projection must show name, value, unit,
    # and convention at the surface a reader actually sees.
    extraction = _analysis_node(
        "extract-sp",
        "result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    expression = _analysis_node(
        "apply-standard-state",
        "quantity_expression",
        dependencies=("extract-sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="e-in",
                source_kind="analysis_output",
                producer_node_id="extract-sp",
                producer_output_id="e",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e-corrected",
                quantity_kind="energy",
                unit="kcal/mol",
            ),
        ),
        expression_nodes=(
            {
                "node_id": "ss-correction",
                "operation": "constant",
                "constant_name": ("standard_state_correction_1atm_to_1M_298K"),
            },
            {
                "node_id": "shifted",
                "operation": "add",
                "input_ids": ("e-in", "ss-correction"),
            },
            {
                "node_id": "e-corrected",
                "operation": "convert",
                "input_ids": ("shifted",),
                "target_unit": "kcal/mol",
            },
        ),
        expression_output_node_ids=("e-corrected",),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("apply-standard-state",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="apply-standard-state",
                producer_output_id="e-corrected",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="kcal/mol",
            ),
        ),
    )
    toolchain = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, expression, claims),
        required_output_ids=("final-energy",),
    )
    executor = _executor(tmp_path, toolchain)

    nodes, status, _completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    assert status == "completed"
    assert all(node.state == "executed" for node in nodes)
    report = Path(report_path).read_text(encoding="utf-8")
    assert "## Literature constants" in report
    assert "standard_state_correction_1atm_to_1M_298K" in report
    assert "1.894" in report
    assert "RT ln(24.4654)" in report

    from rich.console import Console

    from chemsmart.agent.tui.report import render_report_for_humans

    console = Console(record=True, width=200)
    console.print(render_report_for_humans(report))
    projected = console.export_text()
    assert "Literature constants" in projected
    assert "standard_state_correction_1atm_to_1M_298K" in projected


def test_a_kernel_refusal_is_a_finding_not_a_crash(tmp_path, monkeypatch):
    # Observed live: an acetate optimization converged onto its
    # methyl-torsion saddle, the thermochemistry kernel refused the
    # meaningless pure-RRHO correction with a bare ValueError, and the
    # refusal escaped the delivery envelope -- the whole run crashed and
    # every surviving receipt was lost. The kernel's refusal must settle
    # the node and reach the partial report as a finding.
    extraction = _analysis_node(
        "extract-sp",
        "result_extraction",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    thermo = _analysis_node(
        "thermo",
        "thermochemistry",
        dependencies=("sp",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="raw",
                source_kind="program_output",
                producer_node_id="sp",
                producer_output_id="sp-out",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="gcorr",
                quantity_kind="thermal_gibbs_correction",
                unit="hartree",
            ),
        ),
        temperature_k=298.15,
        pressure_atm=1.0,
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("thermo",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-gcorr",
                source_kind="analysis_output",
                producer_node_id="thermo",
                producer_output_id="gcorr",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-gcorr",
                quantity_kind="thermal_gibbs_correction",
                unit="kcal/mol",
            ),
        ),
    )
    toolchain = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, thermo, claims),
        required_output_ids=("final-gcorr",),
    )
    executor = _executor(tmp_path, toolchain)

    import chemsmart.agent.tool_runtime as tool_runtime_module

    def _kernel_refusal(**_kwargs):
        raise ValueError(
            "!! ERROR: Detected imaginary frequencies in geometry "
            "optimization. A valid optimized geometry should not contain "
            "imaginary frequencies."
        )

    monkeypatch.setattr(
        tool_runtime_module,
        "derive_trusted_thermochemistry",
        _kernel_refusal,
    )

    nodes, status, completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    states = {node.node_id: node.state for node in nodes}
    assert states == {
        "extract-sp": "executed",
        "thermo": "failed",
        "claims": "skipped",
    }
    thermo_record = next(node for node in nodes if node.node_id == "thermo")
    assert "imaginary frequencies" in thermo_record.reason
    assert status == "partial"
    assert completions
    report = Path(report_path).read_text(encoding="utf-8")
    assert "Partial analysis" in report
    assert "imaginary frequencies" in report
    # The extraction's receipt survives as evidence at its rung.
    assert "parsed" in report


def test_no_analysis_node_may_end_the_goal(tmp_path, monkeypatch):
    """Whatever a kernel raises, the node fails and the run continues.

    Widening the caught tuple one exception class at a time was the
    wrong repair, and the third strike proved it: a kernel reached into
    a result's absent thermochemistry, ``float + None`` raised
    TypeError, and the escape killed a goal outright -- three hours
    twenty minutes of engine time, five converged spin states, a ledger
    three entries long with no run_recorded and no settlement.

    The invariant is not "these errors settle". It is that no analysis
    node may end the goal, so this drives the exception classes that are
    *not* typed refusals and asserts the same outcome for each.
    """

    import chemsmart.agent.tool_runtime as tool_runtime_module

    for error in (
        TypeError("unsupported operand type(s) for +: 'float' and 'NoneType'"),
        KeyError("total_internal_energy"),
        AttributeError("'NoneType' object has no attribute 'real'"),
        ZeroDivisionError("float division by zero"),
    ):
        extraction = _analysis_node(
            "extract-sp",
            "result_extraction",
            dependencies=("sp",),
            inputs=(
                AnalysisInputIntentV1(
                    input_id="raw",
                    source_kind="program_output",
                    producer_node_id="sp",
                    producer_output_id="sp-out",
                ),
            ),
            selectors=(
                AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
            ),
            outputs=(
                AnalysisOutputIntentV1(
                    output_id="e", quantity_kind="energy", unit="hartree"
                ),
            ),
        )
        thermo = _analysis_node(
            "thermo",
            "thermochemistry",
            dependencies=("sp",),
            inputs=(
                AnalysisInputIntentV1(
                    input_id="raw",
                    source_kind="program_output",
                    producer_node_id="sp",
                    producer_output_id="sp-out",
                ),
            ),
            outputs=(
                AnalysisOutputIntentV1(
                    output_id="gcorr",
                    quantity_kind="thermal_gibbs_correction",
                    unit="hartree",
                ),
            ),
            temperature_k=298.15,
            pressure_atm=1.0,
        )
        toolchain = build_scientific_toolchain_plan(
            plan_id="p",
            workflow_id="w",
            command_workflow_draft_sha256="9" * 64,
            calculation_nodes=(_calculation(),),
            calculation_observables={"sp": ("sp-out",)},
            analysis_nodes=(extraction, thermo),
            required_output_ids=("gcorr",),
        )
        executor = _executor(
            tmp_path / f"run-{type(error).__name__}", toolchain
        )

        def _raise(**_kwargs):
            raise error

        monkeypatch.setattr(
            tool_runtime_module, "derive_trusted_thermochemistry", _raise
        )
        nodes, status, _completions, _report = executor._run_analysis_phase(
            toolchain
        )
        states = {node.node_id: node.state for node in nodes}
        assert states["thermo"] == "failed", type(error).__name__
        # The sibling still ran: a failure settles its own node only.
        assert states["extract-sp"] == "executed", type(error).__name__
        assert status == "partial", type(error).__name__
        record = next(n for n in nodes if n.node_id == "thermo")
        # Recorded as a host defect, not dressed up as a scientific
        # finding -- the reason carries the exception's own type.
        assert "host defect" in record.reason, type(error).__name__
        assert type(error).__name__ in record.reason


def test_a_result_without_thermochemistry_is_refused_not_computed():
    """A run that never reached its frequency step has no thermochemistry.

    The kernel used to reach straight into the absent values, and
    ``electronic_energy + total_internal_energy`` is a TypeError when
    the second is None. Driven here on real archived output rather than
    a fake: the refusal must name what is absent.
    """

    import hashlib

    from chemsmart.analysis.result_quantities import (
        QuantityExtractionError,
        ThermochemistryRequestV1,
        derive_result_thermochemistry,
    )

    # An archived optimisation with no frequency step: the same shape a
    # run that hit its iteration cap leaves behind.
    path = Path("tests/data/ORCATests/outputs/ethanol_fixed_bond.out")
    request = ThermochemistryRequestV1(
        schema_version="chemsmart.thermochemistry-request.v1",
        artifact_id="orca-result-nonconverged",
        artifact_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        program="orca",
        temperature_k=298.15,
        pressure_atm=1.0,
        concentration_mol_l=None,
        entropy_method="rrho",
        entropy_cutoff_cm1=None,
        enthalpy_cutoff_cm1=None,
        alpha=4,
        use_weighted_mass=True,
        frequency_scale_factor=1.0,
    )
    with pytest.raises(QuantityExtractionError, match="no thermochemistry"):
        derive_result_thermochemistry(request=request, artifact_path=str(path))


def test_a_registered_result_input_resolves_from_the_workspace(tmp_path):
    # Observed live: the first composed chain to reuse a registered acid
    # result crashed with a bare KeyError -- the fresh executor host held
    # only this run's artifacts, and the approved registered input was
    # never resolved. The id is content-derived, so the approved
    # workspace rebuilds exactly the artifact the planning session
    # registered.
    from chemsmart.agent.scientific_toolchain import (
        RegisteredResultInputIntentV1,
    )

    workspace_nodes = tmp_path / "workspace" / "nodes" / "prior-sp"
    workspace_nodes.mkdir(parents=True)
    staged = workspace_nodes / "CO2.out"
    staged.write_bytes(_RESULT.read_bytes())
    registered_id = "orca-result-8ae1cdc683f8eb7d"

    extraction = _analysis_node(
        "extract-registered",
        "result_extraction",
        dependencies=(),
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="raw", artifact_id=registered_id
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("extract-registered",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="extract-registered",
                producer_output_id="e",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
    )
    toolchain = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, claims),
        required_output_ids=("final-energy",),
    )
    executor = _executor(tmp_path, toolchain)

    nodes, status, _completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    assert status == "completed"
    assert all(node.state == "executed" for node in nodes)
    report = Path(report_path).read_text(encoding="utf-8")
    assert "final-energy" in report


def test_a_missing_registered_result_is_a_finding_not_a_crash(tmp_path):
    from chemsmart.agent.scientific_toolchain import (
        RegisteredResultInputIntentV1,
    )

    extraction = _analysis_node(
        "extract-registered",
        "result_extraction",
        dependencies=(),
        inputs=(
            RegisteredResultInputIntentV1(
                input_id="raw",
                artifact_id="orca-result-0000000000000000",
            ),
        ),
        selectors=(
            AnalysisSelectorIntentV1(quantity_id="e", selector="energy"),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="e", quantity_kind="energy", unit="hartree"
            ),
        ),
    )
    claims = _analysis_node(
        "claims",
        "claim_rendering",
        dependencies=("extract-registered",),
        inputs=(
            AnalysisInputIntentV1(
                input_id="final-energy",
                source_kind="analysis_output",
                producer_node_id="extract-registered",
                producer_output_id="e",
            ),
        ),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="final-energy",
                quantity_kind="energy",
                unit="hartree",
            ),
        ),
    )
    toolchain = build_scientific_toolchain_plan(
        plan_id="p",
        workflow_id="w",
        command_workflow_draft_sha256="9" * 64,
        calculation_nodes=(_calculation(),),
        calculation_observables={"sp": ("sp-out",)},
        analysis_nodes=(extraction, claims),
        required_output_ids=("final-energy",),
    )
    executor = _executor(tmp_path, toolchain)

    nodes, status, completions, report_path = executor._run_analysis_phase(
        toolchain
    )

    states = {node.node_id: node.state for node in nodes}
    assert states == {"extract-registered": "failed", "claims": "skipped"}
    failed = next(
        node for node in nodes if node.node_id == "extract-registered"
    )
    assert "orca-result-0000000000000000" in failed.reason
    assert status == "partial"
    assert completions
    report = Path(report_path).read_text(encoding="utf-8")
    assert "not present in the approved workspace" in report
