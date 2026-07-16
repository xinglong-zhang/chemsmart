from __future__ import annotations

import json

from chemsmart.agent.harness.evaluation import (
    load_case_matrix,
    reliability_metrics,
)
from chemsmart.agent.harness.failure_taxonomy import classify_runtime_failure
from chemsmart.agent.harness.generated_invariants import (
    check_generated_input_invariants,
)
from chemsmart.agent.harness.command_semantics import CommandSemanticResult
from chemsmart.agent.harness.intent import IntentSpec, evaluate_intent
from chemsmart.agent.harness.terminal_state import (
    assertion,
    build_terminal_state,
    terminal_state_is_positive,
    validate_terminal_state,
)
from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    reset_workflow_state,
    select_project_from_request,
    select_workspace_project,
    workflow_state_scope,
)
from chemsmart.agent.project_yaml import (
    read_project_yaml,
    update_project_yaml,
    validate_project_yaml,
)
from chemsmart.agent.synthesis import SynthesisSession


def _project_text() -> str:
    return (
        "gas:\n  functional: b3lyp\n  basis: def2svp\n"
        "solv:\n  functional: b3lyp\n  basis: def2svp\n"
    )


def test_explicit_project_selection_wins_with_multiple_workspace_yaml(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    for name in ("demo", "co2cat"):
        path = tmp_path / ".chemsmart" / "gaussian" / f"{name}.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_project_text(), encoding="utf-8")

    reset_workflow_state()
    selected = select_project_from_request(
        "Prepare a Gaussian optimization using the co2cat project."
    )

    assert selected and selected["selected"] is True
    state = current_workflow_state()
    assert state.project is not None
    assert state.project.name == "co2cat"
    assert len(state.project.sha256) == 64


def test_workspace_project_selection_is_isolated_per_agent_session(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    for name in ("alpha", "beta"):
        path = tmp_path / ".chemsmart" / "gaussian" / f"{name}.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_project_text(), encoding="utf-8")
    reset_workflow_state()

    with workflow_state_scope("session-a"):
        select_workspace_project("alpha", "gaussian")
        assert current_workflow_state().project.name == "alpha"
    with workflow_state_scope("session-b"):
        select_workspace_project("beta", "gaussian")
        assert current_workflow_state().project.name == "beta"
    with workflow_state_scope("session-a"):
        assert current_workflow_state().project.name == "alpha"


def test_project_read_uses_session_selection_with_multiple_workspace_yaml(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    for name in ("alpha", "beta"):
        path = tmp_path / ".chemsmart" / "gaussian" / f"{name}.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_project_text(), encoding="utf-8")
    reset_workflow_state()

    with workflow_state_scope("session-a"):
        select_workspace_project("beta", "gaussian")
        loaded = read_project_yaml()

    assert loaded["ok"] is True
    assert loaded["project_name"] == "beta"
    assert loaded["path"].endswith("/.chemsmart/gaussian/beta.yaml")


def test_project_read_and_stringified_update_refresh_state(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    path = tmp_path / ".chemsmart" / "gaussian" / "demo.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(_project_text(), encoding="utf-8")
    reset_workflow_state()

    read = read_project_yaml("demo", "gaussian")
    before_hash = read["state_delta"]["project"]["sha256"]
    updated = update_project_yaml(
        updates=json.dumps({"gas.freq": True}),
        project_name="demo",
        program="gaussian",
    )

    assert read["ok"] is True
    assert updated["ok"] is True
    assert updated["state_delta"]["project"]["sha256"] != before_hash


def test_project_yaml_rejects_basis_name_outside_program_catalog():
    result = validate_project_yaml(
        "gas:\n  functional: b3lyp\n  basis: def2-not-real\n"
        "solv:\n  functional: b3lyp\n  basis: def2-not-real\n",
        program="gaussian",
        project_name="bad_basis",
    )

    assert result["verdict"] == "reject"
    assert any(
        issue["rule_id"] == "yaml.basis.unrecognized"
        for issue in result["issues"]
    )


def test_synthesis_uses_explicit_project_when_workspace_has_multiple_yaml(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    for name in ("demo", "co2cat"):
        path = tmp_path / ".chemsmart" / "gaussian" / f"{name}.yaml"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(_project_text(), encoding="utf-8")
    (tmp_path / "water.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n", encoding="utf-8"
    )

    class Provider:
        name = "local"

        def chat(self, _messages):
            return {
                "choices": [
                    {
                        "message": {
                            "content": json.dumps(
                                {
                                    "status": "ready",
                                    "command": (
                                        "chemsmart run gaussian -f water.xyz "
                                        "-c 0 -m 1 opt"
                                    ),
                                    "confidence": "high",
                                }
                            )
                        }
                    }
                ]
            }

    reset_workflow_state()
    session = SynthesisSession(
        provider=Provider(),
        default_project="",
        semantic_gate=False,
        enable_intent_router=False,
    )
    result = session.prepare_command(
        "Run a Gaussian optimization for water.xyz using the co2cat project, "
        "charge 0 multiplicity 1."
    )

    assert result["status"] == "ready"
    assert "gaussian -p co2cat " in result["command"]


def test_synthesis_exposes_semantic_and_intent_evidence(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    path = tmp_path / ".chemsmart" / "gaussian" / "demo.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(_project_text(), encoding="utf-8")
    (tmp_path / "water.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
        encoding="utf-8",
    )

    class Provider:
        name = "local"

        def chat(self, _messages):
            return {
                "choices": [
                    {
                        "message": {
                            "content": json.dumps(
                                {
                                    "status": "ready",
                                    "command": (
                                        "chemsmart run gaussian -p demo "
                                        "-f water.xyz -c 0 -m 1 opt"
                                    ),
                                    "confidence": "high",
                                }
                            )
                        }
                    }
                ]
            }

    monkeypatch.setattr(
        "chemsmart.agent.synthesis.evaluate_command_semantics",
        lambda command, **_kwargs: CommandSemanticResult(
            verdict="ok",
            command=command,
        ),
    )
    reset_workflow_state()
    result = SynthesisSession(
        provider=Provider(),
        default_project="",
        enable_intent_router=False,
    ).prepare_command(
        "Run a Gaussian optimization for water.xyz using demo, neutral singlet."
    )

    assert result["semantic"]["verdict"] == "ok"
    assert result["intent_assertion"]["verdict"] == "ok"
    assert result["intent_assertion"]["failed_rule_ids"] == []


def test_intent_contract_detects_file_drift_and_preserves_scan_fields():
    expected = IntentSpec.from_dict(
        {
            "action": "sub",
            "program": "orca",
            "kind": "orca.scan",
            "project": "demo",
            "server": "hpc1",
            "input_path": "complex.xyz",
            "charge": 0,
            "multiplicity": 1,
            "chemistry": {
                "coordinates": [2, 7],
                "dist_start": 1.8,
                "dist_end": 3.0,
                "num_steps": 13,
            },
        }
    )
    command = (
        "chemsmart sub -s hpc1 orca -p demo -f complex.xyz -c 0 -m 1 "
        "scan -c '[2,7]' -x 1.8 -y 3.0 -n 13"
    )
    valid = evaluate_intent(command, expected)
    drifted = evaluate_intent(command.replace("complex.xyz", "water.xyz"), expected)

    assert valid.verdict == "ok"
    assert drifted.verdict == "reject"
    assert "intent.input_path" in drifted.failed_rule_ids


def test_generated_input_invariants_cover_td_neb_scan_and_sp():
    td = check_generated_input_invariants(
        "chemsmart run gaussian -p photo -f dye.xyz -c 0 -m 1 td "
        "--states singlets --nstates 8 --root 2 --eqsolv",
        [{"path": "dye.com", "route": "# TD(singlets,nstates=8,root=2,eqsolv) CAM-B3LYP", "content_tail": ""}],
    )
    neb = check_generated_input_invariants(
        "chemsmart run orca -p demo -f reactant.xyz -c 0 -m 1 neb "
        "-e product.xyz --nimages 9 --joboption NEB-TS",
        [{"path": "neb.inp", "route": "! B3LYP def2-SVP NEB-TS", "content_tail": "%NEB\nNEB_END_XYZFILE \"product.xyz\"\nNImages 9\nend"}],
    )
    scan = check_generated_input_invariants(
        "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 scan "
        "-c '[[1,2]]' -s '[0.05]' -n '[10]' -cc '[[1,3]]'",
        [{"path": "scan.com", "route": "# opt=modredundant b3lyp/6-31g(d)", "content_tail": "B 1 2 S 10 0.05\nB 1 3 F\n"}],
    )
    sp = check_generated_input_invariants(
        "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 sp",
        [{"path": "sp.com", "route": "# opt freq b3lyp/6-31g(d)", "content_tail": ""}],
    )

    assert td == ()
    assert neb == ()
    assert scan == ()
    assert [issue.rule_id for issue in sp] == ["input.sp.unrequested_route"]


def test_generated_input_invariants_detect_coordinate_and_qmmm_drift():
    gaussian_scan = check_generated_input_invariants(
        "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 scan "
        "-c '[[1,2]]' -s '[0.05]' -n '[10]' -cc '[[1,3]]'",
        [
            {
                "path": "scan.com",
                "route": "# opt=modredundant b3lyp/6-31g(d)",
                "content_tail": "B 1 3 S 10 0.05\nB 1 2 F\n",
            }
        ],
    )
    orca_scan = check_generated_input_invariants(
        "chemsmart run orca -p demo -f water.xyz -c 0 -m 1 scan "
        "-c '[1,2]' -x 0.9 -y 1.5 -n 12",
        [
            {
                "path": "scan.inp",
                "route": "! Scan B3LYP def2-SVP",
                "content_tail": "%geom\nScan\nB 0 2 = 0.9, 1.5, 12\nend\nend",
            }
        ],
    )
    orca_qmmm = check_generated_input_invariants(
        "chemsmart run orca -p demo -f water.xyz -c 0 -m 1 opt "
        "qmmm -ha 1-2 -lm MM.prms",
        [{"path": "qmmm.inp", "route": "! QMMM B3LYP", "content_tail": "%qmmm\nend"}],
    )

    assert "input.gaussian.scan.coordinate_atoms" in {
        issue.rule_id for issue in gaussian_scan
    }
    assert "input.gaussian.scan.constraint" in {
        issue.rule_id for issue in gaussian_scan
    }
    assert "input.orca.scan.coordinate_atoms" in {
        issue.rule_id for issue in orca_scan
    }
    assert "input.orca.qmmm.low_level_method" in {
        issue.rule_id for issue in orca_qmmm
    }


def test_orca_scan_exact_row_rejects_point_count_drift():
    command = (
        "chemsmart run orca -p demo -f complex.xyz -c 0 -m 1 scan "
        "--coordinates '[2,7]' --dist-start 1.8 --dist-end 3.0 "
        "--num-steps 13"
    )
    base = {
        "path": "complex_scan.inp",
        "route": "! B3LYP def2-SVP",
    }
    correct = check_generated_input_invariants(
        command,
        [
            {
                **base,
                "content_tail": (
                    "%geom\n  Scan\n  B 1 6 = 1.8, 3.0, 13\n"
                    "  end\nend\n* xyz 0 1\nH 0 0 0\n*\n"
                ),
            }
        ],
    )
    drifted = check_generated_input_invariants(
        command,
        [
            {
                **base,
                # The unrelated comment contains 13, proving that a global
                # number search cannot validate the scan definition.
                "content_tail": (
                    "%geom\n  Scan\n  B 1 6 = 1.8, 3.0, 12\n"
                    "  end\nend\n# requested 13 points\n"
                    "* xyz 0 1\nH 0 0 0\n*\n"
                ),
            }
        ],
    )

    assert correct == ()
    assert [issue.rule_id for issue in drifted] == [
        "input.orca.scan.num_steps"
    ]


def test_gaussian_scan_exact_row_rejects_step_drift():
    command = (
        "chemsmart run gaussian -p demo -f ethanol.xyz -c 0 -m 1 scan "
        "--coordinates '[[1,2]]' --step-size '[0.05]' --num-steps '[10]'"
    )
    correct = check_generated_input_invariants(
        command,
        [
            {
                "path": "ethanol_scan.com",
                "route": "# opt=modredundant b3lyp/6-31g(d)",
                "content_tail": "B 1 2 S 10 0.05\n",
            }
        ],
    )
    drifted = check_generated_input_invariants(
        command,
        [
            {
                "path": "ethanol_scan.com",
                "route": "# opt=modredundant b3lyp/6-31g(d)",
                "content_tail": "B 1 2 S 9 0.05\n# requested 10 steps\n",
            }
        ],
    )

    assert correct == ()
    assert [issue.rule_id for issue in drifted] == [
        "input.gaussian.scan.num_steps"
    ]


def test_generated_input_invariants_accept_orca_modred_constraint_indices():
    issues = check_generated_input_invariants(
        "chemsmart run orca -p demo -f probe.xyz -c 0 -m 1 modred "
        "--coordinates '[[2,5]]'",
        [
            {
                "path": "probe_modred.inp",
                "route": "! Opt B3LYP def2-SVP",
                "content_tail": (
                    "%geom\n  Constraints\n  {B 1 4 C}\n  end\nend\n"
                    "* xyz 0 1\nC 0.0 0.0 0.0\n*\n"
                ),
            }
        ],
    )

    assert issues == ()


def test_runtime_failure_taxonomy_decomposes_generic_failure():
    project = classify_runtime_failure(
        stderr="No project settings implemented. Currently available projects: ['demo']",
        returncode=1,
    )
    dependency = classify_runtime_failure(
        stderr="ModuleNotFoundError: No module named openbabel",
        returncode=1,
    )

    assert project.rule_id == "cmd.runtime.project_not_found"
    assert dependency.rule_id == "cmd.runtime.dependency_missing"


def test_terminal_profile_rejects_missing_required_assertion():
    state = build_terminal_state(
        action="submit_job",
        command="chemsmart sub -s hpc1 gaussian -p demo -f a.xyz -c 0 -m 1 opt",
        assertions=[assertion("command.returncode", expected=0, observed=0)],
        returncode=0,
        required_assertion_ids=["command.returncode", "sub.scheduler_marker"],
    )

    assert terminal_state_is_positive(state) is False
    assert "terminal_state.required_missing:sub.scheduler_marker" in validate_terminal_state(state)


def test_frozen_matrix_has_48_cases_and_reliability_metrics():
    cases = load_case_matrix(
        "tests/agent/harness/fixtures/high_risk_matrix.json"
    )
    rows = [
        {"case_id": "a", "passed": True},
        {"case_id": "a", "passed": True},
        {"case_id": "a", "passed": False},
        {"case_id": "b", "passed": False},
        {"case_id": "b", "passed": True},
        {"case_id": "b", "passed": True},
    ]
    metrics = reliability_metrics(rows, k=3)

    assert len(cases) == 48
    assert len({case.family for case in cases}) == 12
    assert metrics["pass_at_1"] == 4 / 6
    assert metrics["pass_at_3"] == 1.0
    assert metrics["pass_power_3"] == 0.0


def test_matrix_runner_executes_ready_command_in_safe_test_mode(
    monkeypatch, tmp_path
):
    import scripts.harness.run_agentic_matrix as matrix_runner

    matrix = tmp_path / "matrix.json"
    matrix.write_text(
        json.dumps(
            {
                "families": [
                    {
                        "id": "gaussian_opt",
                        "fixture": {
                            "input": "water.xyz",
                            "workspace_projects": ["demo"],
                        },
                        "intent": {
                            "action": "run",
                            "program": "gaussian",
                            "kind": "gaussian.opt",
                            "project": "demo",
                            "input_path": "water.xyz",
                            "charge": 0,
                            "multiplicity": 1,
                        },
                        "variants": [
                            {
                                "id": "complete",
                                "turns": ["Optimize water.xyz with demo."],
                                "expected_outcome": "direct_pass",
                            }
                        ],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    command = (
        "chemsmart run gaussian -p demo -f water.xyz -c 0 -m 1 opt"
    )

    class Provider:
        name = "test"
        model = "test-model"

    class Session:
        def __init__(self, **_kwargs):
            self._last_raw_response = "ready"

        def prepare_command(self, _turn):
            return {
                "status": "ready",
                "command": command,
                "semantic": {"verdict": "ok", "failed_rule_ids": []},
            }

    seen: list[bool] = []

    def fake_execute(_command, *, test, timeout_s):
        seen.append(test)
        return {
            "ok": True,
            "status": "ok",
            "test": test,
            "returncode": 0,
            "executed_argv": ["chemsmart", "run", "--fake"],
            "terminal_state": build_terminal_state(
                action="execute_command",
                command=command,
                returncode=0,
                assertions=[
                    assertion("command.returncode", expected=0, observed=0),
                ],
            ),
        }

    monkeypatch.setattr(matrix_runner, "get_provider", lambda: Provider())
    monkeypatch.setattr(matrix_runner, "SynthesisSession", Session)
    monkeypatch.setattr(
        matrix_runner,
        "execute_chemsmart_command",
        fake_execute,
    )

    rows = matrix_runner.run_batch(matrix, offset=0, limit=1, repeats=1)

    assert seen == [True]
    assert rows[0]["passed"] is True
    assert rows[0]["terminal_state"]["status"] == "passed"


def test_matrix_runner_does_not_emit_terminal_failure_before_execution(
    monkeypatch, tmp_path
):
    import scripts.harness.run_agentic_matrix as matrix_runner

    matrix = tmp_path / "matrix.json"
    matrix.write_text(
        json.dumps(
            {
                "families": [
                    {
                        "id": "gaussian_opt",
                        "fixture": {"workspace_projects": ["demo"]},
                        "intent": {
                            "program": "gaussian",
                            "kind": "gaussian.opt",
                        },
                        "variants": [
                            {
                                "id": "bad_command",
                                "turns": ["Optimize water.xyz."],
                                "expected_outcome": "direct_pass",
                            }
                        ],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    class Provider:
        name = "test"
        model = "test-model"

    class Session:
        def __init__(self, **_kwargs):
            self._last_raw_response = "invalid"

        def prepare_command(self, _turn):
            return {
                "status": "needs_clarification",
                "command": "",
                "semantic": {
                    "verdict": "reject",
                    "failed_rule_ids": [
                        "cmd.contract.job_subcommand_required"
                    ],
                },
            }

    monkeypatch.setattr(matrix_runner, "get_provider", lambda: Provider())
    monkeypatch.setattr(matrix_runner, "SynthesisSession", Session)

    rows = matrix_runner.run_batch(matrix, offset=0, limit=1, repeats=1)

    assert rows[0]["passed"] is False
    assert rows[0]["execution_stage"] == "not_reached"
    assert rows[0]["terminal_state"] is None
    assert rows[0]["terminal_rule_ids"] == []
