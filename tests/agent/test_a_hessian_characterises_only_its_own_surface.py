"""A Hessian characterises the geometry it consumed, on its own surface.

The host credits an optimisation's stationary point to the Hessian that
consumed its reached geometry, so a two-node minimum gets the word a
one-node opt+freq gets. It credited that Hessian whatever surface it was
computed on. A ground-state Hessian at the geometry an excited-root
optimisation reached, or a cheap Hessian at a correlated minimum, is a
real number about a real structure that says nothing about whether that
structure is a stationary point of the surface the optimisation walked
on -- and round 2 delivered exactly that shape of number as "the S1
geometry" while the better structure sat unclaimed.

The records here are the host's own: each fixture is evaluated through
the organ the goal driver calls, so the surface in every observation
block was written by the reader and not by this test. Three cases, three
answers -- agreement, disagreement, and a comparison that cannot be made
because an older contract recorded no surface -- and the third is not
silently the first.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.driver import _analysis_delivery
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _output_artifact_kind,
)
from chemsmart.io.pyscf.output import read_pyscf_h5

FIXTURES = (
    Path(__file__).resolve().parents[1] / "data" / "PySCFTests" / "outputs"
)

pytestmark = pytest.mark.capability("program_jobtype:pyscf:cpu:hess")

#: producer directory, consumer directory, and what the host should say
#: about whether the consumer characterises the producer's structure.
_PAIRS = (
    ("water_dft_opt_v6", "water_dft_hess_on_dft_geometry", True),
    ("water_mp2_opt_v6", "water_dft_hess_on_mp2_geometry", False),
    # Written under contract v5, which records no surface at all.
    ("water_opt", "water_hess", None),
)


def _artifacts(case: str) -> tuple[str, tuple[TrustedArtifactRefV1, ...]]:
    directory = FIXTURES / case
    label = sorted(directory.glob("*.h5"))[0].stem
    found = []
    for path in sorted(directory.glob(f"{label}*")):
        if path.suffix == ".json" and "reference" in path.name:
            continue
        kind = (
            "pyscf_hdf5"
            if path.suffix == ".h5"
            else _output_artifact_kind("pyscf", path)
        )
        found.append(
            TrustedArtifactRefV1(
                artifact_id=f"result.{path.name}",
                kind=kind,
                sha256=file_sha256(path),
                size_bytes=path.stat().st_size,
                path=str(path),
                cli_value=str(path),
            )
        )
    return label, tuple(found)


def _verified_record(case: str, node_id: str) -> dict:
    """One result as the host's own evaluator describes it."""

    label, outputs = _artifacts(case)
    spec, _provenance, _status, _results = read_pyscf_h5(
        str(FIXTURES / case / f"{label}.h5")
    )
    evaluation = CommandCompiledToolHostV1._evaluate_execution_outputs(
        program="pyscf",
        jobtype=str(spec["jobtype"]),
        charge=int(spec["charge"]),
        multiplicity=int(spec["multiplicity"]),
        expected_settings={
            "jobtype": spec["jobtype"],
            **{
                field: spec[field]
                for field in ("ab_initio", "excited_state_root", "frozen_core")
                if spec.get(field) is not None
            },
        },
        output_artifacts=outputs,
        exit_status=0,
    )
    return {
        "node_id": node_id,
        "state": "valid",
        "jobtype": str(spec["jobtype"]),
        "observations": dict(evaluation.observations),
        "output_artifacts": [
            {"sha256": artifact.sha256} for artifact in outputs
        ],
    }


def _stream(tmp_path, producer_case, consumer_case, name):
    rows = [
        {
            "kind": "program_result_verified",
            "payload": {"record": _verified_record(producer_case, "producer")},
        },
        {
            "kind": "program_result_verified",
            "payload": {"record": _verified_record(consumer_case, "consumer")},
        },
        {
            "kind": "optimized_geometry_handed_off",
            "payload": {
                "status": "validated_handoff",
                "producer_node_id": "producer",
                "consumer_node_id": "consumer",
            },
        },
    ]
    path = tmp_path / f"{name}.jsonl"
    path.write_text(
        "\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8"
    )
    return path


@pytest.mark.parametrize("producer,consumer,agree", _PAIRS)
def test_only_a_same_surface_hessian_characterises_the_geometry(
    producer, consumer, agree, tmp_path
):
    stream = _stream(tmp_path, producer, consumer, "run")
    delivery = _analysis_delivery(stream)
    mismatched = delivery.surface_mismatched_characterisations
    if agree is False:
        assert mismatched == (("producer", "consumer"),), (
            f"{consumer} was computed on another surface from {producer} "
            "and the host credited it anyway"
        )
    else:
        # Agreement, and a comparison that cannot be made, both credit
        # the Hessian -- a program that records no surface keeps the
        # behaviour it had -- and neither is recorded as a mismatch.
        assert mismatched == ()


def test_the_surfaces_the_fixtures_carry_are_the_ones_compared():
    """The discrimination is in the recorded surfaces, not in this test."""

    from chemsmart.analysis.result_readers import reader_for, surfaces_agree
    from chemsmart.io.pyscf.output import PySCFOutput

    reader = reader_for("pyscf")
    surfaces = {}
    for case in {name for pair in _PAIRS for name in pair[:2]}:
        label, _outputs = _artifacts(case)
        surfaces[case] = reader.surface_for_output(
            PySCFOutput(str(FIXTURES / case / f"{label}.h5"))
        )
    for producer, consumer, agree in _PAIRS:
        assert (
            surfaces_agree(surfaces[producer], surfaces[consumer]) is agree
        ), (producer, consumer)
