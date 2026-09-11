"""A thermochemistry derivation names the modes below 50 cm-1 and the RRHO
entropy they carry, beside its receipt.

Two enantiomeric rotamers differed by 0.06 kJ/mol electronically and
0.40 kJ/mol in Gibbs energy from 8-16 cm-1 modes under harmonic RRHO, on
a 2 kJ/mol design decision (NOVEL-1/2 po2, 2026-09-04). Invisible from
inside one run unless the host says it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.analysis.result_quantities import low_frequency_mode_entropy

pytestmark = pytest.mark.capability("tool:derive_thermochemistry")


def test_the_entropy_of_low_modes_is_the_harmonic_oscillator_term():
    summary = low_frequency_mode_entropy(
        (12.0, 40.0, 600.0), temperature_k=298.15
    )
    assert summary["low_modes_cm1"] == (12.0, 40.0)
    # A 12 cm-1 oscillator at 298 K carries about 33 J/mol/K.
    assert 30.0 < summary["entropy_j_per_mol_k"] < 60.0
    assert summary["entropy_term_kj_per_mol"] > 9.0
    assert (
        low_frequency_mode_entropy((600.0,), temperature_k=298.15)[
            "low_modes_cm1"
        ]
        == ()
    )


def _derive(tmp_path, fixture):
    resolved = Path(fixture).resolve()
    artifact = TrustedArtifactRefV1(
        artifact_id="result",
        kind="orca_output",
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="thermo"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path,
    )
    host.artifacts[artifact.artifact_id] = artifact
    reply = host.dispatch(
        turn_id="t1",
        tool_name="derive_thermochemistry",
        arguments={
            "artifact_id": "result",
            "program": "orca",
            "temperature_k": 298.15,
            "pressure_atm": 1.0,
        },
    )
    return host, reply


def test_a_result_with_low_modes_is_observed(tmp_path):
    host, reply = _derive(
        tmp_path / "low", "tests/data/ORCATests/outputs/L2_ts1_opt_pka_A.out"
    )
    assert reply["status"] == "ok"
    (observation,) = reply["observations"]
    assert observation["kind"] == "low_frequency_modes_under_rrho"
    assert observation["low_modes_cm1"][0] == pytest.approx(31.54, abs=0.01)
    assert "entropy_method" in observation["meaning"]
    derived = [
        event
        for event in host.event_store.read_events()
        if event.kind == "thermochemistry_derived"
    ]
    assert derived[-1].payload["observations"][0]["kind"] == (
        "low_frequency_modes_under_rrho"
    )


def test_a_result_with_no_low_modes_is_silent(tmp_path):
    _host, reply = _derive(
        tmp_path / "high", "tests/data/ORCATests/outputs/CO2.out"
    )
    assert reply["status"] == "ok"
    assert "observations" not in reply
