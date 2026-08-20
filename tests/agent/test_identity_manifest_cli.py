from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner
import pytest
import yaml

from chemsmart.agent._contracts import (
    ContractError,
    canonical_sha256,
    file_sha256,
)
from chemsmart.agent.identity import load_approved_molecular_input_manifest
from chemsmart.agent.live_session import (
    _approved_input_state_bindings,
    _scan_xyz_artifacts,
    _task_spec_sha256,
)
from chemsmart.cli.agent import agent


class _Result:
    terminal_state = "planned"
    successful_tool_calls = 0
    failed_tool_calls = 0
    final_text = ""
    prepared_execution = None

    def public_summary_json(self):
        return "{}"


def _xyz(path: Path, symbol: str, x: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"1\nsource coordinate\n{symbol} {x:.3f} 0.0 0.0\n",
        encoding="utf-8",
    )
    return path


def _manifest(tmp_path: Path, workspace: Path) -> Path:
    neutral = _xyz(workspace / "inputs" / "neutral.xyz", "H", 0.0)
    cation = _xyz(workspace / "inputs" / "cation.xyz", "H", 0.1)
    records = []
    for input_id, name, role, charge, multiplicity, geometry in (
        (
            "hydrogen-neutral",
            "hydrogen atom",
            "neutral starting geometry",
            0,
            2,
            neutral,
        ),
        (
            "hydrogen-cation",
            "hydrogen cation",
            "cation starting geometry",
            1,
            1,
            cation,
        ),
    ):
        source_record = {
            "input_id": input_id,
            "identity": name,
            "geometry_role": role,
            "charge": charge,
            "multiplicity": multiplicity,
        }
        records.append(
            {
                "input_id": input_id,
                "identity_id": input_id,
                "approved_names": [name],
                "geometry_file": str(geometry.relative_to(workspace)),
                "geometry_sha256": file_sha256(geometry),
                "coordinate_units": "angstrom",
                "geometry_role": role,
                "charge": charge,
                "multiplicity": multiplicity,
                "source_locator": "https://example.test/source",
                "source_record_sha256": canonical_sha256(source_record),
                "state_source_locator": f"benchmark-case:test#{input_id}",
            }
        )
    path = tmp_path / "identity-manifest.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "schema_version": (
                    "chemsmart.approved-molecular-input-manifest.v1"
                ),
                "inputs": records,
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return path


def test_manifest_loads_multiple_path_free_identity_and_state_records(
    tmp_path,
):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    manifest = _manifest(tmp_path, workspace)

    records = load_approved_molecular_input_manifest(
        manifest, workspace=workspace
    )

    assert tuple(item.input_id for item in records) == (
        "hydrogen-cation",
        "hydrogen-neutral",
    )
    public = tuple(item.public_record() for item in records)
    assert public[0]["approved_names"] == ["hydrogen cation"]
    assert public[0]["geometry_role"] == "cation starting geometry"
    assert public[0]["charge"] == 1
    assert public[0]["multiplicity"] == 1
    assert public[0]["electronic_state_status"] == (
        "user_approved_for_exact_geometry"
    )
    assert "geometry_file" not in str(public)


def test_manifest_rejects_changed_geometry_bytes(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    manifest = _manifest(tmp_path, workspace)
    geometry = workspace / "inputs" / "neutral.xyz"
    geometry.write_text(
        "1\nchanged coordinate\nH 9.0 0.0 0.0\n", encoding="utf-8"
    )

    with pytest.raises(ContractError, match="geometry digest differs"):
        load_approved_molecular_input_manifest(manifest, workspace=workspace)


def test_manifest_states_are_prebound_to_exact_runtime_artifacts(tmp_path):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    manifest = _manifest(tmp_path, workspace)
    approved = load_approved_molecular_input_manifest(
        manifest, workspace=workspace
    )
    observations = _scan_xyz_artifacts(workspace)
    task_sha256 = _task_spec_sha256(
        "compare the two approved states",
        observations,
        tuple(item.molecular_identity for item in approved),
        approved_molecular_inputs=approved,
    )

    bindings = _approved_input_state_bindings(
        approved,
        observations=observations,
        task_spec_sha256=task_sha256,
    )

    assert {
        (item.charge, item.multiplicity) for item in bindings.values()
    } == {
        (0, 2),
        (1, 1),
    }
    assert all(
        item.task_spec_sha256 == task_sha256 for item in bindings.values()
    )


def test_agent_run_loads_identity_manifest_into_existing_session_plumbing(
    tmp_path, monkeypatch
):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    manifest = _manifest(tmp_path, workspace)
    secret = tmp_path / "secret.env"
    secret.write_text("KEY=value\n", encoding="utf-8")
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return _Result()

    import chemsmart.agent.live_session as live_session

    monkeypatch.setattr(live_session, "run_live_agent_session", fake_run)
    result = CliRunner().invoke(
        agent,
        [
            "plan",
            "--task",
            "compare the two approved states",
            "--workspace",
            str(workspace),
            "--identity-manifest",
            str(manifest),
        ],
    )

    assert result.exit_code == 0, result.output
    approved = captured["approved_molecular_inputs"]
    assert len(approved) == 2
    assert {item.charge for item in approved} == {0, 1}
    assert captured["execution_enabled"] is False
