from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, canonical_sha256, file_sha256
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    approved_xyz_source,
    approved_xyz_source_from_ledger,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import qwen_pyscf_cases_v1
from chemsmart.agent.experiments.qwen_pyscf_grading import (
    grade_qwen_pyscf_episode,
)
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
    source = approved_xyz_source(
        xyz,
        expected_sha256=digest,
        molecular_identity=identity,
    )

    record = source.public_record()["molecular_identity"]
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


def test_nist_coordinate_ledger_builds_path_free_approved_identity(tmp_path):
    ledger = (
        Path(__file__).parents[1]
        / "data"
        / "agent"
        / "pyscf_24h"
        / "coordinate-sources.json"
    ).resolve()
    source = approved_xyz_source_from_ledger(
        ledger, artifact_id="nist-h2o-experimental-cartesian"
    )

    public = source.public_record()
    assert public["molecular_identity"]["approved_names"] == ["H2O", "water"]
    assert "path" not in json.dumps(public).lower()

    copied = tmp_path / "approved.xyz"
    copied.write_bytes(source.path.read_bytes())
    observations = _scan_xyz_artifacts(tmp_path)
    records = _validated_identity_records(
        observations, source.molecular_identity
    )
    with_identity = _task_spec_sha256(
        "preview task", observations, source.molecular_identity
    )
    without_identity = _task_spec_sha256("preview task", observations)
    assert records == (public["molecular_identity"],)
    assert with_identity != without_identity


def _decision_message(*, evidence_refs: tuple[str, ...]) -> dict:
    result = {
        "assumptions": (
            "The approved molecular system is water; its electronic state remains missing.",
        ),
        "method_rationale": "No method is selected while the state is unresolved.",
        "alternatives": (),
        "uncertainties": ("Charge and multiplicity are missing.",),
        "diagnostics": (),
        "evidence_refs": evidence_refs,
    }
    return {
        "role": "tool",
        "content": json.dumps(
            {
                "schema_version": "chemsmart.tool-result.v1",
                "tool": "record_scientific_decision",
                "status": "ok",
                "result": result,
            }
        ),
    }


def test_grader_requires_identity_record_reference_when_public_name_is_used():
    case = next(
        item for item in qwen_pyscf_cases_v1() if item.case_id == "QP-TR-001"
    )
    identity = _identity(case.source_sha256s[0])
    artifact_records = (identity.public_record(),)

    ungrounded = grade_qwen_pyscf_episode(
        case=case,
        live_result={
            "terminal_state": "blocked",
            "public_transcript": (_decision_message(evidence_refs=()),),
            "artifact_records": artifact_records,
            "successful_tool_calls": 1,
            "failed_tool_calls": 0,
        },
        legacy_transcript_fallback=True,
    )
    grounded = grade_qwen_pyscf_episode(
        case=case,
        live_result={
            "terminal_state": "blocked",
            "public_transcript": (
                _decision_message(evidence_refs=(identity.evidence_ref,)),
            ),
            "artifact_records": artifact_records,
            "successful_tool_calls": 1,
            "failed_tool_calls": 0,
        },
        legacy_transcript_fallback=True,
    )

    assert ungrounded.verdict == "fail"
    assert grounded.verdict == "pass"


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
