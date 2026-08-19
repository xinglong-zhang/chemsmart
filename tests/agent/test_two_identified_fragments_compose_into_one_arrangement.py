"""Nothing could join two approved monomers; now the host can, with lineage.

Four cycle-039 sessions ended with the paper's observable uncomputable
because no affordance builds a two-molecule arrangement -- one session
probed four invented option names looking for it. The host now owns the
placement mathematics (contact pair at the requested distance, everything
else clash-free and maximally separated, reusing the iterate module's
deterministic SLSQP machinery); the model owns the scientific choices and
must bind the electronic state explicitly afterwards.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.runtime.events import EventKind
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

_WATER = "3\nwater\nO 0.0 0.0 0.117\nH 0.0 0.757 -0.467\nH 0.0 -0.757 -0.467\n"
_AMMONIA = (
    "4\nammonia\nN 0.0 0.0 0.0\nH 0.937 0.0 -0.381\n"
    "H -0.468 0.811 -0.381\nH -0.468 -0.811 -0.381\n"
)


def _host_with_fragments(tmp_path, bind_identities=True):
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    (tmp_path / "workspace").mkdir(exist_ok=True)
    for name, text in (("water", _WATER), ("ammonia", _AMMONIA)):
        path = tmp_path / f"{name}.xyz"
        path.write_text(text)
        artifact = TrustedArtifactRefV1(
            artifact_id=f"geometry-{name}",
            kind="geometry_xyz",
            sha256=file_sha256(path),
            size_bytes=path.stat().st_size,
            path=str(path),
            cli_value=str(path),
        )
        host.artifacts[artifact.artifact_id] = artifact
        if bind_identities:
            host.dispatch(
                turn_id="t0",
                tool_name="bind_scientific_identity",
                arguments={
                    "input_artifact_id": artifact.artifact_id,
                    "charge": 0,
                    "multiplicity": 1,
                },
            )
    return host


def test_the_arrangement_lands_in_the_workspace_with_lineage(tmp_path):
    host = _host_with_fragments(tmp_path)

    result = host.dispatch(
        turn_id="t1",
        tool_name="compose_molecular_arrangement",
        arguments={
            "composed_artifact_id": "water-ammonia-complex",
            "fragment_a_artifact_id": "geometry-water",
            "fragment_b_artifact_id": "geometry-ammonia",
            "fragment_a_atom": 1,
            "fragment_b_atom": 1,
            "distance_angstrom": 2.8,
        },
    )["result"]

    composition = result["composition"]
    artifact = result["artifact"]
    composed_path = Path(
        str(tmp_path / "workspace" / "artifacts" / "water-ammonia-complex.xyz")
    )
    assert composed_path.exists()
    assert composition["atom_count"] == 7
    assert composition["formula"] == "H5NO"
    assert abs(
        composition["achieved_contact_distance_angstrom"] - 2.8
    ) < 0.05
    assert composition["min_interfragment_distance_angstrom"] > 0.5
    assert composition["fragment_a_artifact_id"] == "geometry-water"
    assert composition["fragment_b_artifact_id"] == "geometry-ammonia"
    assert composition["atom_order_note"].startswith("fragment A")
    assert "electronic state deliberately unbound" in composed_path.read_text()
    assert "bind charge and multiplicity explicitly" in result["next_action"]
    assert artifact["kind"] == "geometry_xyz"

    kinds = [event.kind for event in host.event_store.read_events()]
    assert EventKind.MOLECULAR_ARRANGEMENT_COMPOSED.value in kinds

    # The composed artifact itself carries NO identity until the model
    # binds one -- composition never invents an electronic state.
    assert not any(
        binding.geometry_artifact_sha256 == artifact["sha256"]
        for binding in host.scientific_identities.values()
    )


def test_an_unidentified_fragment_cannot_compose(tmp_path):
    host = _host_with_fragments(tmp_path, bind_identities=False)

    with pytest.raises(ContractError) as refusal:
        host.dispatch(
            turn_id="t1",
            tool_name="compose_molecular_arrangement",
            arguments={
                "composed_artifact_id": "water-ammonia-complex",
                "fragment_a_artifact_id": "geometry-water",
                "fragment_b_artifact_id": "geometry-ammonia",
                "fragment_a_atom": 1,
                "fragment_b_atom": 1,
                "distance_angstrom": 2.8,
            },
        )

    message = str(refusal.value)
    assert "carries no scientific identity" in message
    assert "bind_scientific_identity" in message
