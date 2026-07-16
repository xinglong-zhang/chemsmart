from __future__ import annotations

from chemsmart.agent.model_command_parser import (
    format_model_command_explanation,
    parse_model_command,
)


def test_parse_gaussian_project_method_and_position_dependent_p(tmp_path):
    parsed = parse_model_command(
        "chemsmart run gaussian -p gas_solv -f examples/h2o.xyz "
        "-c 0 -m 1 -l h2o_sp sp",
        cwd=tmp_path,
    )

    assert parsed.parse_error is None
    assert parsed.action == "run"
    assert parsed.program == "gaussian"
    assert parsed.job == "sp"
    assert parsed.project == "gas_solv"
    assert parsed.project_p_flag_meaning is not None
    assert parsed.top_level_program is None
    assert parsed.filename == "examples/h2o.xyz"
    assert parsed.label == "h2o_sp"
    assert parsed.charge == "0"
    assert parsed.multiplicity == "1"
    assert parsed.functional == "b3lyp empiricaldispersion=gd3bj"
    assert parsed.basis == "def2tzvp"
    assert parsed.solvent_model == "smd"
    assert parsed.solvent_id == "toluene"
    assert parsed.dry_run is False


def test_parse_top_level_program_p_is_not_project(tmp_path):
    parsed = parse_model_command(
        "chemsmart run -p gaussian gaussian -p gas_solv "
        "-f examples/h2o.xyz sp",
        cwd=tmp_path,
    )

    assert parsed.parse_error is None
    assert parsed.top_level_program == "gaussian"
    assert parsed.project == "gas_solv"
    text = format_model_command_explanation(parsed.command, cwd=tmp_path)
    assert "top-level `-p/--program`" in text
    assert "program-level -p/--project" in text
    assert text.splitlines()[-1].startswith("Summary: This command will")
    # A Markdown renderer must not lazily glue the Summary onto the last
    # bullet: there must be a blank line separating them.
    text_lines = text.splitlines()
    summary_index = next(
        i for i, line in enumerate(text_lines) if line.startswith("Summary:")
    )
    assert text_lines[summary_index - 1] == ""


def test_parse_submission_dry_run_and_structural_options(tmp_path):
    parsed = parse_model_command(
        "chemsmart sub -s gadi --test gaussian -p gas_solv "
        "-f ts.xyz -c 0 -m 1 ts --freeze-atoms 1,2",
        cwd=tmp_path,
    )

    assert parsed.action == "sub"
    assert parsed.server == "gadi"
    assert parsed.dry_run is True
    assert parsed.job == "ts"
    assert parsed.structural_options["freeze_atoms"] == "1,2"
    assert parsed.functional == "b3lyp empiricaldispersion=gd3bj"
    assert parsed.basis == "def2svp"


def test_parse_orca_sp_project_method_from_repo_fixture(tmp_path):
    parsed = parse_model_command(
        "chemsmart run orca -p orca -f examples/h2o.xyz "
        "-c 0 -m 1 -l h2o_0_1 sp",
        cwd=tmp_path,
    )

    assert parsed.parse_error is None
    assert parsed.program == "orca"
    assert parsed.job == "sp"
    assert parsed.project == "orca"
    assert parsed.ab_initio is None
    assert parsed.functional == "m062x"
    assert parsed.basis == "def2-tzvp"
    assert parsed.solvent_model == "smd"
    assert parsed.solvent_id == "cyclohexane"

    text = format_model_command_explanation(parsed.command, cwd=tmp_path)
    assert "functional `m062x`" in text
    assert "basis `def2-tzvp`" in text
    assert "resolved solvent: `smd` / `cyclohexane`" in text


def test_parse_database_selectors_and_orca_aux_basis(tmp_path):
    parsed = parse_model_command(
        "chemsmart run orca -p orca --record-id abc123 --structure-index 2 "
        "-f results.db -c 0 -m 1 -B def2/J sp",
        cwd=tmp_path,
    )

    assert parsed.parse_error is None
    assert parsed.record_id == "abc123"
    assert parsed.structure_index == "2"
    assert parsed.filename == "results.db"
    assert parsed.aux_basis == "def2/J"

    text = format_model_command_explanation(parsed.command, cwd=tmp_path)
    assert "database selection: `record_id=abc123, structure_index=2`" in text
    assert "auxiliary basis `def2/J`" in text


def test_parse_gaussian_userjob_route_after_subcommand(tmp_path):
    parsed = parse_model_command(
        "chemsmart run gaussian -f anisole.xyz -c 0 -m 1 userjob --route nmr=giao",
        cwd=tmp_path,
    )

    assert parsed.job == "userjob"
    assert parsed.route_parameters == "nmr=giao"
    assert "unrecognized option `--route`" not in parsed.warnings


def test_parse_gaussian_td_short_options_after_subcommand(tmp_path):
    parsed = parse_model_command(
        "chemsmart run gaussian -f dye.xyz -c 0 -m 1 "
        "td -s singlets -n 8 -r 3 -e eqsolv",
        cwd=tmp_path,
    )

    assert parsed.job == "td"
    assert parsed.route_parameters is None
    assert parsed.structural_options == {
        "states": "singlets",
        "nstates": "8",
        "root": "3",
        "eqsolv": "eqsolv",
    }
    assert not parsed.warnings


def test_parse_gaussian_traj_selection_options_after_subcommand(tmp_path):
    parsed = parse_model_command(
        "chemsmart run gaussian -f md_frames.xyz -c 0 -m 1 "
        "traj -j opt -ns 6 -x 0.2",
        cwd=tmp_path,
    )

    assert parsed.job == "traj"
    assert parsed.structural_options == {
        "jobtype": "opt",
        "num_structures_to_run": "6",
        "proportion_structures_to_use": "0.2",
    }
    assert not parsed.warnings
