"""A converged search on a stationary point of the wrong order is a
typed failure whatever program found it, and the capability cell says
whether the host can judge it.

Before this, only an ORCA transition-state search was checked, and no
program's minimum was: an opt+freq that landed on a saddle validated.
The rule is program-neutral because the physics is -- the declared
jobtype promises a count of imaginary modes and the program's own
frequencies deliver one -- and the 20 cm-1 convention is the one the
thermochemistry already uses for numerical noise.
"""

from __future__ import annotations

import math

from chemsmart.agent.capabilities import (
    COVERAGE_AXES,
    CapabilityQueryV1,
    coverage_for,
    load_program_capabilities,
    query_capability,
)
from chemsmart.agent.terminal_states import (
    STATIONARY_POINT_ORDER_FINDING,
    _classify_failure,
    consequential_imaginary_mode_count,
    expected_imaginary_mode_count,
    stationary_point_order_finding,
)
from chemsmart.analysis.result_readers import reader_for


def test_the_declared_jobtype_promises_a_count():
    assert expected_imaginary_mode_count("ts") == 1
    assert expected_imaginary_mode_count("opt") == 0
    assert expected_imaginary_mode_count("hess") == 0
    assert expected_imaginary_mode_count("sp") is None
    assert expected_imaginary_mode_count("scan") is None


def test_noise_below_the_convention_is_not_an_imaginary_mode():
    assert consequential_imaginary_mode_count((-15.0, 40.0, 1200.0)) == 0
    assert consequential_imaginary_mode_count((-45.2, 40.0)) == 1
    assert consequential_imaginary_mode_count((-300.0, -25.0, 10.0)) == 2
    assert consequential_imaginary_mode_count(()) is None
    assert consequential_imaginary_mode_count(None) is None
    assert consequential_imaginary_mode_count((math.nan, 10.0)) is None


def test_the_finding_fires_only_on_a_mismatch_with_a_claim():
    assert stationary_point_order_finding("ts", 0) == (
        STATIONARY_POINT_ORDER_FINDING
    )
    assert stationary_point_order_finding("ts", 2) == (
        STATIONARY_POINT_ORDER_FINDING
    )
    assert stationary_point_order_finding("ts", 1) == ""
    assert stationary_point_order_finding("opt", 1) == (
        STATIONARY_POINT_ORDER_FINDING
    )
    assert stationary_point_order_finding("opt", 0) == ""
    # No frequencies printed: an absent fact, never a verdict.
    assert stationary_point_order_finding("opt", None) == ""
    # No claim made by the jobtype: nothing to judge.
    assert stationary_point_order_finding("sp", 3) == ""


def test_the_finding_classifies_as_a_wrong_stationary_point():
    state = _classify_failure(
        jobtype="opt",
        findings=(STATIONARY_POINT_ORDER_FINDING,),
        native_class="",
        converged=True,
        reached=None,
        planned=None,
    )
    assert state == "failed_wrong_stationary_point"
    assert (
        _classify_failure(
            jobtype="ts",
            findings=("orca.result.ts_imaginary_mode_count",),
            native_class="",
            converged=True,
            reached=None,
            planned=None,
        )
        == "failed_wrong_stationary_point"
    )


def _cell(program, jobtype):
    reader = reader_for(program)
    selectors = reader.selectors_for_jobtype(jobtype)
    assert selectors is not None, (program, jobtype)
    axes, rules = coverage_for(program, jobtype, selectors)
    return dict(axes), rules


def test_the_coverage_cell_says_what_the_host_can_judge():
    orca_ts_axes, orca_ts_rules = _cell("orca", "ts")
    assert set(orca_ts_axes) == set(COVERAGE_AXES)
    assert orca_ts_axes["identity"] == "validated"
    assert orca_ts_axes["thermochemistry"] == "readable"
    assert "stationary_point_order" in orca_ts_rules
    assert "spin_square_observed" in orca_ts_rules

    _gaussian_opt_axes, gaussian_opt_rules = _cell("gaussian", "opt")
    assert "stationary_point_order" in gaussian_opt_rules

    xtb_hess_axes, xtb_hess_rules = _cell("xtb", "hess")
    assert "stationary_point_order" in xtb_hess_rules
    assert xtb_hess_axes["spin"] == "unsupported"

    # This asserted ``not in`` while PySCF's verification branch fed the
    # rule nothing; the cell was honest and nothing read it. The rule is
    # now fed through every program's reader, so the cell derives from
    # the declaration alone and PySCF's Hessian is judged like xTB's.
    pyscf_hess_axes, pyscf_hess_rules = _cell("pyscf", "hess")
    assert "stationary_point_order" in pyscf_hess_rules
    assert "spin_square_observed" in pyscf_hess_rules
    assert pyscf_hess_axes["identity"] == "validated"

    _orca_sp_axes, orca_sp_rules = _cell("orca", "sp")
    assert "convergence" not in orca_sp_rules
    assert "stationary_point_order" not in orca_sp_rules


def test_the_capability_receipt_carries_the_cell():
    registry = load_program_capabilities()
    receipt = query_capability(
        CapabilityQueryV1(program="orca", jobtype="ts", engine="cpu"),
        registry=registry,
    )
    coverage = receipt.job_result_selector_coverage
    assert coverage is not None
    assert "stationary_point_order" in coverage.validity_rules
    assert dict(coverage.axes)["geometry"] == "readable"


def test_the_xtb_sidecar_is_not_a_second_log_and_the_log_gives_the_count():
    """Observed live (W1, methanol xTB opt+hess): g98.out registered as a
    second xtb_output and the analysis walker refused to choose; and the
    result receipt carries no frequencies, so the count came back None
    and the stationary-point rule never reached xTB. The log is the
    source, read by the typed layer's own parser."""

    from pathlib import Path

    from chemsmart.agent.tool_runtime import (
        _output_artifact_kind,
        _xtb_log_frequencies,
    )

    assert _output_artifact_kind("xtb", Path("run/g98.out")) == (
        "program_output"
    )
    assert _output_artifact_kind("xtb", Path("run/meoh_hess.out")) == (
        "xtb_output"
    )
    assert _output_artifact_kind("orca", Path("run/a.hess")) == "orca_hessian"
    assert _output_artifact_kind("orca", Path("run/a.out")) == "orca_output"

    log = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "XTBTests"
        / ("outputs/co2_ohess/co2_ohess.out")
    )
    frequencies = _xtb_log_frequencies((log,))
    assert frequencies, "the archived hess log prints frequencies"
    assert consequential_imaginary_mode_count(frequencies) == 0
