"""A validated saddle with a soft imaginary mode is an observation.

REACH-1 po3 cycle 4 (2026-09-06): a transition-state search relaxed to
a van der Waals complex -- azide 2.9 to 3.3 A from an unreacted alkyne
-- with one imaginary mode at -22.8 cm-1. The stationary-point rule
counts modes below the 20 cm-1 noise convention, so it certified the
structure and the host's word was `validated`; the session had named
the forming-bond distances as its own check, in prose. A rule at a
threshold certifies noise on the far side of the threshold. The number
is now recorded as an anomaly observation with standing, never a
verdict: the same archived SN2 saddle is no anomaly at -407.58 cm-1
and one at -22.80.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_SN2 = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "ORCATests"
    / "outputs"
    / "sn2_ts.out"
)
_HARD = "  -407.58 cm**-1  ***imaginary mode***"
_SOFT = "   -22.80 cm**-1  ***imaginary mode***"


def _artifact(path: Path) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id="result.sn2",
        kind="orca_output",
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


def _evaluate(path: Path):
    return CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="orca",
        jobtype="ts",
        charge=-1,
        multiplicity=1,
        output_artifacts=(_artifact(path),),
        exit_status=0,
    )


def _signals(evaluation) -> dict[str, dict]:
    return {item["signal_id"]: item for item in evaluation.anomalies}


@pytest.mark.capability("predicate:stationary_point.imaginary_mode_lt_50")
def test_a_bond_forming_saddle_is_no_soft_mode_anomaly():
    signals = _signals(_evaluate(_SN2))
    assert "stationary_point.imaginary_mode_lt_50" not in signals
    assert "stationary_point.unexpected_order" not in signals


@pytest.mark.capability("predicate:stationary_point.imaginary_mode_lt_50")
def test_a_saddle_inside_the_band_is_observed_with_its_number(tmp_path):
    text = _SN2.read_text()
    assert text.count(_HARD) == 1
    soft = tmp_path / "sn2_soft.out"
    soft.write_text(text.replace(_HARD, _SOFT))
    evaluation = _evaluate(soft)
    assert "result.stationary_point_order" not in evaluation.findings
    signals = _signals(evaluation)
    anomaly = signals["stationary_point.imaginary_mode_lt_50"]
    assert anomaly["imaginary_mode_cm1"] == -22.8
    assert anomaly["noise_convention_cm1"] == 20.0
    assert anomaly["soft_band_cm1"] == 50.0
