from __future__ import annotations

import yaml

from chemsmart.agent.harness.command_semantics import evaluate_command_semantics
from chemsmart.agent.harness.generated_invariants import (
    check_generated_input_invariants,
)
from chemsmart.agent.harness.intent import IntentSpec, evaluate_intent


COMMAND = (
    "chemsmart run gaussian -p enzyme -f enzyme.xyz -c -1 -m 2 "
    "opt qmmm -ha 1-2 -ct -1 -mt 2 -ch 0 -mh 1"
)
THREE_LAYER_COMMAND = (
    "chemsmart run gaussian -p enzyme -f enzyme.xyz -c 0 -m 1 "
    "sp qmmm -ha 1 -ma 2 -la 3-4 -ct 0 -mt 1 "
    "-ci -1 -mi 2 -ch 0 -mh 1"
)


def _generated(**overrides):
    row = {
        "path": "enzyme.com",
        "route": "# opt oniom(b3lyp/6-31g(d):uff)",
        "content_tail": (
            "# opt oniom(b3lyp/6-31g(d):uff)\n\nenzyme\n\n"
            "-1 2 0 1 0 1\nFe 0 0 0 H\nN 0 0 1 H\n"
            "C 0 0 2 L\nH 0 1 2 L\n"
        ),
        "charge": -1,
        "multiplicity": 2,
        "element_counts": {"C": 1, "Fe": 1, "H": 1, "N": 1},
        "charge_multiplicity_pairs": [[-1, 2], [0, 1], [0, 1]],
        "atom_layers": ["H", "H", "L", "L"],
        "layer_atoms": {"H": [1, 2], "M": [], "L": [3, 4]},
    }
    row.update(overrides)
    return row


def _failed_ids(generated):
    return {
        issue.rule_id
        for issue in check_generated_input_invariants(COMMAND, [generated])
    }


def test_gaussian_qmmm_generated_partition_and_state_pass():
    assert _failed_ids(_generated()) == set()


def test_gaussian_qmmm_generated_high_layer_drift_is_rejected():
    generated = _generated(
        atom_layers=["H", "L", "H", "L"],
        layer_atoms={"H": [1, 3], "M": [], "L": [2, 4]},
    )

    assert "input.gaussian.qmmm.high_level_atoms" in _failed_ids(generated)
    assert "input.gaussian.qmmm.low_level_atoms" in _failed_ids(generated)


def test_gaussian_qmmm_generated_missing_layer_marker_is_rejected():
    generated = _generated(
        atom_layers=["H", None, "L", "L"],
        layer_atoms={"H": [1], "M": [], "L": [3, 4]},
    )

    assert "input.gaussian.qmmm.partition_coverage" in _failed_ids(generated)


def test_gaussian_qmmm_generated_total_state_drift_is_rejected():
    generated = _generated(
        charge=0,
        multiplicity=1,
        charge_multiplicity_pairs=[[0, 1], [0, 1], [0, 1]],
    )

    assert "input.gaussian.qmmm.total_state" in _failed_ids(generated)


def test_gaussian_qmmm_three_layer_states_pass():
    generated = _generated(
        charge=0,
        multiplicity=1,
        charge_multiplicity_pairs=[
            [0, 1],
            [-1, 2],
            [-1, 2],
            [0, 1],
            [0, 1],
            [0, 1],
        ],
        atom_layers=["H", "M", "L", "L"],
        layer_atoms={"H": [1], "M": [2], "L": [3, 4]},
    )

    assert (
        check_generated_input_invariants(THREE_LAYER_COMMAND, [generated]) == ()
    )


def test_gaussian_qmmm_three_layer_intermediate_state_drift_is_rejected():
    generated = _generated(
        charge=0,
        multiplicity=1,
        charge_multiplicity_pairs=[
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
            [0, 1],
        ],
        atom_layers=["H", "M", "L", "L"],
        layer_atoms={"H": [1], "M": [2], "L": [3, 4]},
    )
    failed = {
        issue.rule_id
        for issue in check_generated_input_invariants(
            THREE_LAYER_COMMAND,
            [generated],
        )
    }

    assert "input.gaussian.qmmm.intermediate_state" in failed


def test_gaussian_qmmm_intent_observes_parent_job_and_layer_state():
    result = evaluate_intent(
        COMMAND,
        IntentSpec.from_dict(
            {
                "action": "run",
                "program": "gaussian",
                "kind": "gaussian.qmmm",
                "project": "enzyme",
                "input_path": "enzyme.xyz",
                "charge": -1,
                "multiplicity": 2,
                "chemistry": {
                    "parent_job": "opt",
                    "high_level_atoms": "1-2",
                    "charge_total": -1,
                    "mult_total": 2,
                    "charge_high": 0,
                    "mult_high": 1,
                },
            }
        ),
    )

    assert result.verdict == "ok"


def test_gaussian_qmmm_fake_runtime_preserves_partition_and_state(
    tmp_path, monkeypatch
):
    home = tmp_path / "home"
    server_dir = home / ".chemsmart" / "server"
    server_dir.mkdir(parents=True)
    (server_dir / "local.yaml").write_text(
        yaml.safe_dump(
            {
                "SERVER": {"SCHEDULER": None, "NUM_CORES": 1},
                "GAUSSIAN": {
                    "EXEFOLDER": "/tmp",
                    "LOCAL_RUN": True,
                    "SCRATCH": False,
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("HOME", str(home))

    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "enzyme.yaml").write_text(
        yaml.safe_dump(
            {
                "gas": {"functional": "b3lyp", "basis": "6-31g(d)"},
                "solv": {"functional": "b3lyp", "basis": "6-31g(d)"},
                "qmmm": {
                    "high_level_functional": "b3lyp",
                    "high_level_basis": "6-31g(d)",
                    "low_level_force_field": "uff",
                    "freq": False,
                },
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "enzyme.xyz").write_text(
        "4\nenzyme fragment\nFe 0 0 0\nN 0 0 1.8\nC 0 0 3.0\nH 0 1 3.0\n",
        encoding="utf-8",
    )

    result = evaluate_command_semantics(COMMAND, cwd=tmp_path)

    assert result.verdict == "ok", result.to_dict()
    generated = result.generated_inputs[0]
    assert generated["charge_multiplicity_pairs"] == [
        [-1, 2],
        [0, 1],
        [0, 1],
    ]
    assert generated["layer_atoms"] == {"H": [1, 2], "M": [], "L": [3, 4]}
    assert "freq" not in str(generated["route"]).lower()


def test_gaussian_qmmm_rejects_explicit_frequency_policy_drift(tmp_path):
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    (project_dir / "enzyme.yaml").write_text(
        yaml.safe_dump({"qmmm": {"freq": False}}), encoding="utf-8"
    )
    generated = _generated()
    generated["route"] = "# opt freq oniom(b3lyp/6-31g(d):uff)"
    generated["content_tail"] = generated["route"]

    failed = {
        issue.rule_id
        for issue in check_generated_input_invariants(
            COMMAND,
            [generated],
            cwd=str(tmp_path),
        )
    }

    assert "input.gaussian.qmmm.frequency_policy" in failed
