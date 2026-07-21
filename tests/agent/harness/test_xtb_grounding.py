"""End-to-end harness coverage for the three xTB job kinds.

These tests exercise the real safe-execution gate (``run xtb ...`` with
``--fake --no-scratch``), the deterministic contract and generated-input
invariants, and the intent/parser layer, mirroring the Gaussian and ORCA
grounding suites.
"""

from __future__ import annotations

import yaml

from chemsmart.agent.harness.command_semantics import (
    evaluate_command_semantics,
)
from chemsmart.agent.harness.generated_invariants import (
    check_generated_input_invariants,
)
from chemsmart.agent.harness.intent import ObservedIntent
from chemsmart.agent.harness.sub_intent import build_sub_intent_assertions

WATER_XYZ = "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n"


def _write_xtb_workspace(tmp_path, project="test"):
    project_dir = tmp_path / ".chemsmart" / "xtb"
    project_dir.mkdir(parents=True)
    (project_dir / f"{project}.yaml").write_text(
        yaml.safe_dump(
            {
                "opt": {
                    "gfn_version": "gfn2",
                    "optimization_level": "tight",
                    "charge": 0,
                    "multiplicity": 1,
                },
                "sp": {"gfn_version": "gfn2", "charge": 0, "multiplicity": 1},
                "hess": {
                    "gfn_version": "gfn2",
                    "charge": 0,
                    "multiplicity": 1,
                },
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "water.xyz").write_text(WATER_XYZ, encoding="utf-8")
    return tmp_path


def _failed_ids(rows):
    return {row["id"] for row in rows if row["status"] == "fail"}


# --- Safe-execution gate: the three kinds pass and preserve method ---------


def test_xtb_opt_fake_runtime_passes_and_renders_method(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz opt"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "ok", result.to_dict()
    route = result.generated_inputs[0]["route"]
    assert "--gfn 2" in route
    assert "--opt tight" in route


def test_xtb_sp_fake_runtime_passes(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz sp"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "ok", result.to_dict()
    assert "--opt" not in result.generated_inputs[0]["route"]


def test_xtb_hess_fake_runtime_passes(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz hess"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "ok", result.to_dict()
    assert "--hess" in result.generated_inputs[0]["route"]


def test_xtb_solvation_is_preserved_in_generated_call(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz -sm alpb -si toluene sp"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "ok", result.to_dict()
    assert "--alpb toluene" in result.generated_inputs[0]["route"]


# --- Reject cases -----------------------------------------------------------


def test_xtb_half_specified_solvation_rejects(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz -sm alpb sp"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "reject", result.to_dict()
    assert "cmd.contract.xtb_solvent_pair" in result.failed_rule_ids


def test_xtb_unknown_leaf_rejects(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz md"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "reject", result.to_dict()


def test_xtb_missing_leaf_rejects(tmp_path):
    _write_xtb_workspace(tmp_path)
    command = "chemsmart run xtb -p test -f water.xyz"

    result = evaluate_command_semantics(command, cwd=tmp_path)

    assert result.verdict == "reject", result.to_dict()
    assert "cmd.contract.job_subcommand_required" in result.failed_rule_ids


# --- Generated-input invariants (rendered call preservation) ----------------


def _xtb_generated(route: str) -> dict:
    return {
        "path": "/tmp/water_sp/water_sp_fake.out",
        "route": route,
        "content_tail": route,
        "charge": 0,
        "multiplicity": 1,
    }


def test_xtb_invariant_rejects_dropped_solvation():
    command = "chemsmart run xtb -p test -f water.xyz -sm alpb -si toluene sp"
    generated = _xtb_generated("xtb water.xyz --gfn 2 --chrg 0 --uhf 0")

    failed = {
        issue.rule_id
        for issue in check_generated_input_invariants(command, [generated])
    }

    assert "input.xtb.solvent_preservation" in failed


def test_xtb_invariant_rejects_changed_gfn_version():
    command = "chemsmart run xtb -p test -f water.xyz -g gfn2 sp"
    generated = _xtb_generated("xtb water.xyz --gfnff --chrg 0 --uhf 0")

    failed = {
        issue.rule_id
        for issue in check_generated_input_invariants(command, [generated])
    }

    assert "input.xtb.gfn_preservation" in failed


def test_xtb_invariant_accepts_preserved_method_and_solvation():
    command = "chemsmart run xtb -p test -f water.xyz -g gfn2 -sm alpb -si toluene sp"
    generated = _xtb_generated(
        "xtb water.xyz --gfn 2 --chrg 0 --uhf 0 --alpb toluene"
    )

    issues = check_generated_input_invariants(command, [generated])

    assert issues == ()


# --- Intent / parser parity -------------------------------------------------


def test_xtb_observed_intent_reads_kind_and_chemistry():
    command = (
        "chemsmart run xtb -p test -f water.xyz -g gfn2 "
        "-sm alpb -si toluene opt --optimization-level tight"
    )

    observed = ObservedIntent.from_command(command)

    assert observed.program == "xtb"
    assert observed.kind == "xtb.opt"
    assert observed.chemistry["gfn_version"] == "gfn2"
    assert observed.chemistry["solvent_model"] == "alpb"
    assert observed.chemistry["solvent_id"] == "toluene"
    assert observed.chemistry["optimization_level"] == "tight"


def test_xtb_sub_intent_preserves_program_kind_and_file():
    command = (
        "chemsmart sub -s mock-pbs --fake xtb -p test "
        "-f examples/xtb/water.xyz -g gfn2 hess"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "xtb",
            "job": "hess",
            "server": "mock-pbs",
            "filename": "examples/xtb/water.xyz",
        },
    )

    assert _failed_ids(rows) == set()


def test_xtb_sub_intent_rejects_pubchem_substitution():
    command = (
        "chemsmart sub -s mock-pbs --fake xtb -p test -P water -g gfn2 sp"
    )
    rows = build_sub_intent_assertions(
        command,
        {
            "program": "xtb",
            "job": "sp",
            "server": "mock-pbs",
            "filename": "examples/xtb/water.xyz",
        },
    )

    assert "sub.intent_filename" in _failed_ids(rows)
