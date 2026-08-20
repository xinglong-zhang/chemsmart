"""A bound identity told the session nothing about the molecule it bound.

In the cycle-042 observations, one of four sessions discovered that a
bound geometry artifact is measurable through the typed layer and verified
its own conformer labels against the coordinates; the other three assumed
their atom numbering and disclosed the assumption for the reviewer to
check. The capability existed -- the xyz reader and the
distance/angle/dihedral operations -- but nothing at the binding site said
so. The binding result now carries the geometry facts and names the
measurement route, while the binding record itself stays digest-frozen.
"""

from __future__ import annotations

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.commands import build_scientific_identity_binding
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_WATER = "3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.467\nH 0.0 -0.757 -0.467\n"


def test_the_binding_result_names_the_molecule_and_the_route(tmp_path):
    xyz = tmp_path / "water.xyz"
    xyz.write_text(_WATER)
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    artifact = TrustedArtifactRefV1(
        artifact_id="geometry-water",
        kind="geometry_xyz",
        sha256=file_sha256(xyz),
        size_bytes=xyz.stat().st_size,
        path=str(xyz),
        cli_value=str(xyz),
    )
    host.artifacts[artifact.artifact_id] = artifact

    envelope = host.dispatch(
        turn_id="t1",
        tool_name="bind_scientific_identity",
        arguments={
            "input_artifact_id": "geometry-water",
            "charge": 0,
            "multiplicity": 1,
        },
    )
    result = envelope["result"]

    assert result["geometry"]["atom_count"] == 3
    assert result["geometry"]["formula"] == "H2O"
    assert "extract_result_quantities" in result["measurement_route"]
    assert "dihedral" in result["measurement_route"]

    # The digest-frozen binding is unchanged by the wrapper: the recorded
    # binding and the surfaced digest both equal the direct build.
    expected = build_scientific_identity_binding(
        task_spec_sha256="a" * 64,
        geometry_artifact=artifact,
        charge=0,
        multiplicity=1,
    )
    assert result["binding_sha256"] == expected.binding_sha256
    assert result["binding"]["binding_sha256"] == expected.binding_sha256
    assert (
        host.scientific_identities[expected.binding_sha256].binding_sha256
        == expected.binding_sha256
    )
