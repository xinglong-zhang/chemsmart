"""A capability cell's thermochemistry axis says what the host's derivation serves.

A session learns whether it can have a free energy from a stage's results
before it plans the stage: ``inspect_program_capability`` returns the
stage's coverage cell, and its thermochemistry axis is the host's word on
the question. Since R10 Q27 the derivation serves the free energy of a
held surface -- ``derive_thermochemistry`` with ``projected_coordinates``
on a modred result that carries its Hessian -- yet the cell was derived
from the reader's frequency and Gibbs selectors alone and said
``unsupported`` for orca/modred and gaussian/modred. A live session cited
exactly those capability receipts when it gave up the free energy of a
point on a torsional profile (R10 Q21 g1-hooh, CUHK 2153623), which the
host now derives at 0.346 kcal/mol (R11 truth, item 3; truth-3, item 2).

Asked here through the public tool surface, for every stage with an
archived result in this repository: the cell says ``readable`` exactly
where the derivation, asked through its own tool, returns a free energy.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.analysis.result_readers import reader_for

pytestmark = [
    pytest.mark.capability("tool:inspect_program"),
    pytest.mark.capability("tool:derive_thermochemistry"),
]

_DATA = Path(__file__).resolve().parents[1] / "data"
#: One archived result per stage, and the coordinates it held (one-based).
STAGES = {
    # R10 Q21 g1-hooh node mod90 (CUHK 2153623): ORCA ``! Opt Freq`` with
    # H-O-O-H held at 90 deg.
    ("orca", "modred"): (
        _DATA
        / "ORCATests/constrained_dihedral"
        / "h2o2_b3lypg_d3bj_def2svp_hooh90_freq.out",
        ((3, 1, 2, 4),),
    ),
    # R10 q7 oracle O2 (CUHK 2150076): Gaussian ``opt=modredundant freq``
    # holding the same torsion at 90 deg.
    ("gaussian", "modred"): (
        _DATA
        / "GaussianTests/constrained_dihedral/h2o2_b3lyp_def2svp_hooh90.log",
        ((3, 1, 2, 4),),
    ),
    # R10 Q21 g2-hooh: the relaxed equilibrium, a stationary point.
    ("orca", "opt"): (_DATA / "ORCATests/hooh_torsion/h2o2_opt_opt.out", ()),
    # A relaxed scan, which prints no Hessian: no free energy of it.
    ("orca", "scan"): (
        _DATA / "ORCATests/scan_completed/h2o2_o_scan.out",
        ((3, 1, 2, 4),),
    ),
}


def _host(tmp_path):
    artifacts = {}
    for (program, jobtype), (path, _held) in STAGES.items():
        resolved = path.resolve()
        artifact_id = f"{program}-{jobtype}-result"
        artifacts[artifact_id] = TrustedArtifactRefV1(
            artifact_id=artifact_id,
            kind=reader_for(program).artifact_kind,
            sha256=hashlib.sha256(resolved.read_bytes()).hexdigest(),
            size_bytes=resolved.stat().st_size,
            path=str(resolved),
            cli_value=str(resolved),
        )
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s"
        ),
        artifacts=artifacts,
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def _cell(host, program, jobtype):
    """The cell the session reads: ``inspect_program``'s capability."""

    reply = host.dispatch(
        turn_id="turn-1",
        tool_name="inspect_program",
        arguments={"program": program, "jobtype": jobtype, "engine": "cpu"},
    )
    assert reply["status"] == "ok", reply
    coverage = reply["result"]["capability"]["job_result_selector_coverage"]
    return dict(coverage["axes"])["thermochemistry"]


def _derives(host, program, jobtype, held):
    """Whether the derivation, asked with the coordinates the result held,
    returns a free energy (a refusal is raised to the loop that runs it)."""

    arguments = {
        "program": program,
        "artifact_id": f"{program}-{jobtype}-result",
        "temperature_k": 298.15,
        "pressure_atm": 1.0,
    }
    if held:
        arguments["projected_coordinates"] = [list(item) for item in held]
    try:
        reply = host.dispatch(
            turn_id="turn-2",
            tool_name="derive_thermochemistry",
            arguments=arguments,
        )
    except ContractError:
        return False
    return reply["status"] == "ok" and any(
        item["quantity_id"] == "gibbs_free_energy"
        for item in reply["result"]["quantities"]
    )


@pytest.mark.parametrize("stage", sorted(STAGES), ids="/".join)
def test_the_cell_says_readable_where_the_derivation_serves(tmp_path, stage):
    program, jobtype = stage
    _path, held = STAGES[stage]
    host = _host(tmp_path)
    served = _derives(host, program, jobtype, held)
    assert _cell(host, program, jobtype) == (
        "readable" if served else "unsupported"
    ), (stage, served)
