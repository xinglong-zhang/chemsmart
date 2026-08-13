from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256, file_sha256
from chemsmart.agent.identity import (
    build_approved_molecular_identity,
    validate_identity_for_geometry,
)
from chemsmart.agent.live_session import (
    _scan_xyz_artifacts,
    _task_spec_sha256,
    _validated_identity_records,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _identity(xyz_sha256: str):
    return build_approved_molecular_identity(
        identity_id="nist-h2o-experimental-cartesian",
        approved_names=("H2O", "water"),
        geometry_sha256=xyz_sha256,
        coordinate_units="angstrom",
        atom_order=("O", "H", "H"),
        source_locator="NIST CCCBDB geometric data, Cartesians",
        source_record_sha256=canonical_sha256(
            {"source": "NIST CCCBDB", "identity": "H2O"}
        ),
    )


def test_approved_identity_binds_names_but_not_electronic_state(tmp_path):
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\nuntrusted comment\nO 0 0 0\nH 0 0.7 0.5\nH 0 -0.7 0.5\n",
        encoding="utf-8",
    )
    digest = file_sha256(xyz)
    identity = _identity(digest)
    record = identity.public_record()
    assert record["approved_names"] == ["H2O", "water"]
    assert record["electronic_state_status"] == (
        "not_established_by_identity_record"
    )
    assert record["evidence_ref"] == f"molecular_identity:{identity.identity_sha256}"


def test_approved_identity_rejects_geometry_or_atom_order_substitution():
    identity = _identity("a" * 64)
    with pytest.raises(ContractError, match="another geometry"):
        validate_identity_for_geometry(
            identity, geometry_sha256="b" * 64, atom_order=("O", "H", "H")
        )
    with pytest.raises(ContractError, match="atom order"):
        validate_identity_for_geometry(
            identity, geometry_sha256="a" * 64, atom_order=("H", "O", "H")
        )


def test_multiple_state_geometries_keep_distinct_approved_identities(tmp_path):
    high_spin = tmp_path / "high-spin.xyz"
    low_spin = tmp_path / "low-spin.xyz"
    high_spin.write_text(
        "3\nhigh-spin source frame\nFe 0 0 0\nO 0 0 2.1\nO 0 0 -2.1\n",
        encoding="utf-8",
    )
    low_spin.write_text(
        "3\nlow-spin source frame\nFe 0 0 0\nO 0 0 1.9\nO 0 0 -1.9\n",
        encoding="utf-8",
    )
    observations = _scan_xyz_artifacts(tmp_path)
    by_name = {Path(item.artifact.path).name: item for item in observations}
    identities = tuple(
        build_approved_molecular_identity(
            identity_id=identity_id,
            approved_names=(approved_name,),
            geometry_sha256=by_name[filename].artifact.sha256,
            coordinate_units="angstrom",
            atom_order=("Fe", "O", "O"),
            source_locator=f"supporting information {approved_name}",
            source_record_sha256=canonical_sha256(
                {"source": "supporting information", "state": approved_name}
            ),
        )
        for identity_id, approved_name, filename in (
            ("fe-aquo-high-spin", "high-spin state", "high-spin.xyz"),
            ("fe-aquo-low-spin", "low-spin state", "low-spin.xyz"),
        )
    )

    records = _validated_identity_records(observations, identities)
    task_sha256 = _task_spec_sha256("compare spin states", observations, identities)

    assert tuple(record["identity_id"] for record in records) == (
        "fe-aquo-high-spin",
        "fe-aquo-low-spin",
    )
    assert len({record["geometry_sha256"] for record in records}) == 2
    assert len(task_sha256) == 64


def test_host_rejects_identity_reference_absent_from_approved_registry(tmp_path):
    task_sha256 = canonical_sha256("identity-task")
    identity = _identity("a" * 64)
    arguments = {
        "decision_id": "identity-decision",
        "task_spec_sha256": task_sha256,
        "assumptions": ["The approved system is water."],
        "method_rationale": "No method inference is made.",
        "alternatives": [],
        "uncertainties": ["Electronic state remains separate."],
        "diagnostics": [],
        "stage_order": ["sp"],
        "evidence_refs": [identity.evidence_ref],
    }
    unbound = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "unbound.jsonl", session_id="unbound"
        ),
        task_spec_sha256s=(task_sha256,),
    )
    with pytest.raises(ContractError, match="unapproved molecular identity"):
        unbound.dispatch(
            turn_id="turn-1",
            tool_name="record_scientific_decision",
            arguments=arguments,
        )

    bound = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "bound.jsonl", session_id="bound"
        ),
        task_spec_sha256s=(task_sha256,),
        approved_molecular_identities={identity.identity_sha256: identity},
    )
    result = bound.dispatch(
        turn_id="turn-1",
        tool_name="record_scientific_decision",
        arguments=arguments,
    )
    assert result["status"] == "ok"


def test_host_rejects_ungrounded_or_unknown_functional_convention(tmp_path):
    task_sha256 = canonical_sha256("functional-task")
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "functional.jsonl", session_id="functional"
        ),
        task_spec_sha256s=(task_sha256,),
    )
    arguments = {
        "decision_id": "functional-decision",
        "task_spec_sha256": task_sha256,
        "assumptions": ["The applied functional uses the VWN3 convention."],
        "method_rationale": "Use the selected DFT method.",
        "alternatives": [],
        "uncertainties": [],
        "diagnostics": [],
        "stage_order": ["sp"],
        "evidence_refs": [],
    }
    with pytest.raises(ContractError, match="require a host resolution"):
        host.dispatch(
            turn_id="turn-1",
            tool_name="record_scientific_decision",
            arguments=arguments,
        )

    arguments["evidence_refs"] = ["functional_resolution:" + "f" * 64]
    with pytest.raises(ContractError, match="unknown functional resolution"):
        host.dispatch(
            turn_id="turn-2",
            tool_name="record_scientific_decision",
            arguments=arguments,
        )
