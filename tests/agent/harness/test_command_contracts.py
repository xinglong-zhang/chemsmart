from __future__ import annotations

import subprocess

import yaml

from chemsmart.agent.harness.command_contracts import check_command_contracts
from chemsmart.agent.harness.command_semantics import evaluate_command_semantics


def _issues(
    program: str,
    job: str,
    *,
    program_tokens: list[str] | None = None,
    job_tokens: list[str] | None = None,
    cwd=None,
):
    return check_command_contracts(
        program=program,
        job=job,
        program_tokens=program_tokens or [],
        job_tokens=job_tokens or [job],
        cwd=cwd,
    )


def test_orca_scan_requires_all_runtime_parameters() -> None:
    issues = _issues(
        "orca",
        "scan",
        job_tokens=["scan", "-c", "[1,2]", "-n", "10"],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.scan_required_parameters"
    ]
    assert "dist_start" in issues[0].evidence["missing"][0]
    assert "dist_end" in issues[0].evidence["missing"][1]


def test_scan_rejects_space_separated_coordinate_string() -> None:
    issues = _issues(
        "orca",
        "scan",
        job_tokens=[
            "scan",
            "-c",
            "1 2",
            "-x",
            "1.4",
            "-y",
            "2.2",
            "-n",
            "10",
        ],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.coordinate_literal"
    ]


def test_scan_accepts_single_and_multi_coordinate_literals() -> None:
    single = _issues(
        "orca",
        "scan",
        job_tokens=[
            "scan",
            "-c",
            "[1,2]",
            "-x",
            "1.4",
            "-y",
            "2.2",
            "-n",
            "10",
        ],
    )
    multiple = _issues(
        "gaussian",
        "scan",
        job_tokens=[
            "scan",
            "-c",
            "[[1,2],[1,2,3]]",
            "-s",
            "[0.05,2.0]",
            "-n",
            "[10,18]",
        ],
    )

    assert single == ()
    assert multiple == ()


def test_scan_rejects_parameter_cardinality_mismatch() -> None:
    issues = _issues(
        "gaussian",
        "scan",
        job_tokens=[
            "scan",
            "-c",
            "[[1,2],[1,2,3]]",
            "-s",
            "[0.05,1.0,2.0]",
            "-n",
            "10",
        ],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.scan_parameter_cardinality"
    ]


def test_gaussian_traj_requires_jobtype() -> None:
    issues = _issues(
        "gaussian",
        "traj",
        job_tokens=["traj", "-ns", "5"],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.traj_jobtype_required"
    ]
    assert _issues(
        "gaussian",
        "traj",
        job_tokens=["traj", "-j", "opt", "-ns", "5"],
    ) == ()


def test_gaussian_traj_rejects_conflicting_selection_modes() -> None:
    issues = _issues(
        "gaussian",
        "traj",
        job_tokens=[
            "traj",
            "-j",
            "opt",
            "-g",
            "rmsd",
            "-ns",
            "5",
        ],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.traj_selection_conflict"
    ]


def test_gaussian_traj_rejects_silent_clamp_values() -> None:
    count_issues = _issues(
        "gaussian",
        "traj",
        job_tokens=["traj", "-j", "opt", "-ns", "0"],
    )
    proportion_issues = _issues(
        "gaussian",
        "traj",
        job_tokens=["traj", "-j", "opt", "-x", "1.2"],
    )

    assert [issue.rule_id for issue in count_issues] == [
        "cmd.contract.traj_structure_count"
    ]
    assert [issue.rule_id for issue in proportion_issues] == [
        "cmd.contract.traj_proportion_range"
    ]


def test_gaussian_crest_rejects_conflicting_or_invalid_selection() -> None:
    conflict = _issues(
        "gaussian",
        "crest",
        job_tokens=["crest", "-j", "opt", "-g", "rmsd", "-nc", "5"],
    )
    invalid_count = _issues(
        "gaussian",
        "crest",
        job_tokens=["crest", "-j", "opt", "-nc", "0"],
    )

    assert [issue.rule_id for issue in conflict] == [
        "cmd.contract.crest_selection_conflict"
    ]
    assert [issue.rule_id for issue in invalid_count] == [
        "cmd.contract.crest_conformer_count"
    ]


def test_dias_requires_flat_fragment_one_range() -> None:
    nested = _issues(
        "gaussian",
        "dias",
        job_tokens=["dias", "-i", "[[1,2],[3,4]]"],
    )
    missing = _issues("gaussian", "dias")
    valid = _issues(
        "gaussian",
        "dias",
        job_tokens=["dias", "-i", "1-4,7"],
    )

    assert [issue.rule_id for issue in nested] == [
        "cmd.contract.dias_fragment_indices"
    ]
    assert [issue.rule_id for issue in missing] == [
        "cmd.contract.dias_fragment_indices_required"
    ]
    assert valid == ()


def test_coordinate_and_dias_indices_must_fit_xyz(tmp_path) -> None:
    (tmp_path / "water.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
        encoding="utf-8",
    )
    scan = _issues(
        "gaussian",
        "scan",
        program_tokens=["-f", "water.xyz"],
        job_tokens=["scan", "-c", "[1,4]", "-s", "0.05", "-n", "10"],
        cwd=tmp_path,
    )
    dias = _issues(
        "gaussian",
        "dias",
        program_tokens=["-f", "water.xyz"],
        job_tokens=["dias", "-i", "1,4"],
        cwd=tmp_path,
    )

    assert [issue.rule_id for issue in scan] == [
        "cmd.contract.coordinate_atom_bounds"
    ]
    assert [issue.rule_id for issue in dias] == [
        "cmd.contract.dias_fragment_atom_bounds"
    ]
    assert scan[0].evidence["atom_count"] == 3


def test_qmmm_requires_parent_job_and_explicit_region() -> None:
    direct = _issues("orca", "qmmm")
    nested = _issues(
        "orca",
        "opt",
        program_tokens=["-c", "0", "-m", "1"],
        job_tokens=["opt", "qmmm"],
    )
    valid = _issues(
        "orca",
        "opt",
        program_tokens=["-c", "0", "-m", "1"],
        job_tokens=["opt", "qmmm", "-ha", "1-6,9"],
    )

    assert [issue.rule_id for issue in direct] == [
        "cmd.contract.qmmm_parent_job"
    ]
    assert [issue.rule_id for issue in nested] == [
        "cmd.contract.qmmm_high_level_atoms_required",
        "cmd.contract.qmmm_low_level_method_required",
    ]
    assert [issue.rule_id for issue in valid] == [
        "cmd.contract.qmmm_low_level_method_required"
    ]


def test_orca_qmmm_accepts_low_level_method_from_command_or_project(tmp_path) -> None:
    command_issues = _issues(
        "orca",
        "opt",
        program_tokens=["-c", "0", "-m", "1"],
        job_tokens=[
            "opt",
            "qmmm",
            "-ha",
            "1-6",
            "-lm",
            "MM.prms",
        ],
        cwd=tmp_path,
    )

    project_dir = tmp_path / ".chemsmart" / "orca"
    project_dir.mkdir(parents=True)
    (project_dir / "enzyme.yaml").write_text(
        yaml.safe_dump(
            {
                "gas": {"functional": "b3lyp", "basis": "def2svp"},
                "qmmm": {
                    "high_level_functional": "b3lyp",
                    "high_level_basis": "def2svp",
                    "low_level_method": "MM.prms",
                },
            }
        ),
        encoding="utf-8",
    )
    project_issues = _issues(
        "orca",
        "opt",
        program_tokens=["-p", "enzyme", "-c", "0", "-m", "1"],
        job_tokens=["opt", "qmmm", "-ha", "1-6"],
        cwd=tmp_path,
    )

    assert command_issues == ()
    assert project_issues == ()

    (project_dir / "enzyme.yaml").write_text(
        yaml.safe_dump(
            {
                "gas": {"functional": "b3lyp", "basis": "def2svp"},
                "qmmm": {"mm_force_field": "AMBER=HardFirst"},
            }
        ),
        encoding="utf-8",
    )
    mm_force_field_issues = _issues(
        "orca",
        "opt",
        program_tokens=["-p", "enzyme", "-c", "0", "-m", "1"],
        job_tokens=["opt", "qmmm", "-ha", "1-6"],
        cwd=tmp_path,
    )
    assert mm_force_field_issues == ()


def test_orca_qmmm_project_mm_force_field_reaches_runtime_writer(
    tmp_path, monkeypatch
) -> None:
    home = tmp_path / "home"
    server_dir = home / ".chemsmart" / "server"
    server_dir.mkdir(parents=True)
    (server_dir / "local.yaml").write_text(
        yaml.safe_dump(
            {
                "SERVER": {"SCHEDULER": None, "NUM_CORES": 1},
                "ORCA": {"EXEFOLDER": "/tmp", "LOCAL_RUN": True, "SCRATCH": False},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("HOME", str(home))

    project_dir = tmp_path / ".chemsmart" / "orca"
    project_dir.mkdir(parents=True)
    (project_dir / "enzyme.yaml").write_text(
        yaml.safe_dump(
            {
                "gas": {"functional": "b3lyp", "basis": "def2svp"},
                "solv": {"functional": "b3lyp", "basis": "def2svp"},
                "qmmm": {"mm_force_field": "AMBER=HardFirst"},
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "enzyme.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
        encoding="utf-8",
    )

    result = evaluate_command_semantics(
        "chemsmart run orca -p enzyme -f enzyme.xyz -c 0 -m 1 "
        "opt qmmm -j QMMM -ha 1-2 -ct 0 -mt 1",
        cwd=tmp_path,
    )

    assert result.verdict == "ok"
    assert result.generated_inputs


def test_qmmm_rejects_space_separated_atom_indices() -> None:
    issues = _issues(
        "gaussian",
        "opt",
        program_tokens=["-c", "0", "-m", "1"],
        job_tokens=[
            "opt",
            "qmmm",
            "-ha",
            "1 2 3",
            "-ct",
            "0",
            "-mt",
            "1",
        ],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.qmmm_high_level_atoms"
    ]


def test_qmmm_regions_must_be_in_bounds_and_non_overlapping(tmp_path) -> None:
    (tmp_path / "water.xyz").write_text(
        "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
        encoding="utf-8",
    )
    issues = _issues(
        "gaussian",
        "sp",
        program_tokens=["-f", "water.xyz"],
        job_tokens=[
            "sp",
            "qmmm",
            "-ha",
            "1-2",
            "-la",
            "2,4",
            "-ct",
            "0",
            "-mt",
            "1",
        ],
        cwd=tmp_path,
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.qmmm_atom_bounds",
        "cmd.contract.qmmm_region_overlap",
    ]
    assert issues[0].evidence["invalid_indices"] == [4]
    assert issues[1].evidence["overlap"] == [2]


def test_gaussian_qmmm_requires_state_in_qmmm_scope() -> None:
    issues = _issues(
        "gaussian",
        "sp",
        program_tokens=["-c", "-1", "-m", "1"],
        job_tokens=["sp", "qmmm", "-ha", "2-9", "-la", "1,10-16"],
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.qmmm_total_state_required"
    ]
    assert issues[0].missing_info == (
        "total charge (-ct/--charge-total after qmmm)",
        "total multiplicity (-mt/--mult-total after qmmm)",
    )


def test_gaussian_td_requires_workspace_td_block(tmp_path) -> None:
    project_dir = tmp_path / ".chemsmart" / "gaussian"
    project_dir.mkdir(parents=True)
    path = project_dir / "paper.yaml"
    path.write_text(
        yaml.safe_dump({"gas": {"functional": "b3lyp", "basis": "def2svp"}}),
        encoding="utf-8",
    )

    issues = _issues(
        "gaussian",
        "td",
        program_tokens=["-p", "paper"],
        cwd=tmp_path,
    )

    assert [issue.rule_id for issue in issues] == [
        "cmd.contract.gaussian_td_project_settings"
    ]

    path.write_text(
        yaml.safe_dump(
            {
                "gas": {"functional": "b3lyp", "basis": "def2svp"},
                "td": {"functional": "cam-b3lyp", "basis": "def2svp"},
            }
        ),
        encoding="utf-8",
    )
    assert _issues(
        "gaussian",
        "td",
        program_tokens=["-p", "paper"],
        cwd=tmp_path,
    ) == ()


def test_direct_qmmm_is_rejected_before_click_or_subprocess(
    monkeypatch,
    tmp_path,
) -> None:
    def should_not_run(*_args, **_kwargs):  # pragma: no cover - defensive
        raise AssertionError("contract should reject before subprocess")

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        should_not_run,
    )
    result = evaluate_command_semantics(
        "chemsmart run orca -f water.xyz -c 0 -m 1 qmmm -ha 1-2",
        cwd=tmp_path,
    )

    assert result.verdict == "reject"
    assert result.failed_rule_ids == ["cmd.contract.qmmm_parent_job"]
    assert "parent computational job before qmmm" in result.missing_info


def test_valid_contract_reaches_safe_execution(monkeypatch, tmp_path) -> None:
    def fake_run(argv, **_kwargs):
        from pathlib import Path

        (Path(_kwargs["cwd"]) / "scan.inp").write_text(
            "! Scan B3LYP def2-SVP\n%geom\n"
            "  Scan B 0 1 = 0.9, 1.5, 12\nend\n"
            "* xyz 0 1\nH 0 0 0\nH 0 0 1\n*\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(
        "chemsmart.agent.harness.command_semantics.subprocess.run",
        fake_run,
    )
    result = evaluate_command_semantics(
        "chemsmart run orca -f water.xyz -c 0 -m 1 scan "
        "-c '[1,2]' -x 0.9 -y 1.5 -n 12",
        cwd=tmp_path,
    )

    assert result.verdict == "ok"
