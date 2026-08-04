"""The agent's environment must be the one ChemSmart declares.

The harness once hard-coded which programs existed.  That made its view of the
host *different from* ChemSmart's rather than a narrower slice of it: a program
the server YAML declared could be invisible, and a program it never declared
could appear.  These tests pin the direction of the dependency -- declaration
first, observation second.
"""

import os
import stat

import pytest

from chemsmart.agent import live_session


@pytest.fixture
def declared_server(tmp_path, monkeypatch):
    """Point ChemSmart at a config dir holding one fixture server YAML."""

    server_dir = tmp_path / "server"
    server_dir.mkdir(parents=True)
    present = tmp_path / "present"
    present.mkdir()
    (present / "orca").write_text("#!/bin/sh\nexit 0\n")
    (present / "orca").chmod(0o700)
    (server_dir / "local.yaml").write_text(
        "SERVER:\n"
        "  SCHEDULER: SLURM\n"
        "  NUM_CORES: 8\n"
        "ORCA:\n"
        f"  EXEFOLDER: {present}\n"
        "GAUSSIAN:\n"
        f"  EXEFOLDER: {tmp_path / 'absent'}\n"
    )
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(tmp_path))
    monkeypatch.setenv("CHEMSMART_AGENT_SERVER", "local")
    return tmp_path


def test_declared_programs_follow_the_server_yaml(declared_server):
    declared = dict(live_session._declared_server_programs())
    assert set(declared) == {"orca", "gaussian"}
    assert declared["orca"].endswith("present")


def test_scheduler_key_is_not_mistaken_for_a_program(declared_server):
    assert "server" not in dict(live_session._declared_server_programs())


def test_a_declared_program_with_no_folder_is_observed_as_missing(
    declared_server,
):
    """Absent is not the same fact as undeclared, so it must not be dropped."""

    records = [
        item
        for item in live_session._observe_environments()[2]
        if item.get("record_kind") == "program_environment"
    ]
    by_program = {item["program"]: item for item in records}
    assert by_program["gaussian"]["status"] == "missing"
    assert by_program["orca"]["status"] == "available"


def test_discovery_targets_cover_every_declared_program(declared_server):
    targets = live_session._observe_environments()[0]
    programs = {item.program for item in targets}
    assert {"orca", "gaussian"} <= programs


def test_executable_name_comes_from_chemsmart_not_the_program_name(tmp_path):
    """Gaussian's binary is ``g16``; guessing ``gaussian`` would never resolve."""

    resolved = live_session._declared_executable_path(
        "gaussian", str(tmp_path)
    )
    assert resolved == str(tmp_path / "g16")
    assert live_session._declared_executable_path(
        "orca", str(tmp_path)
    ) == str(tmp_path / "orca")


def test_a_discovery_stub_is_reported_as_a_stub(tmp_path):
    """A placeholder that cannot compute must not read as a usable program."""

    stub = tmp_path / "g16"
    stub.write_text(
        "#!/bin/sh\n# ChemSmart agent-harness DISCOVERY STUB -- not Gaussian.\n"
        "exit 127\n"
    )
    stub.chmod(stub.stat().st_mode | stat.S_IXUSR)
    assert live_session._executable_is_discovery_stub(str(stub)) is True

    real = tmp_path / "orca"
    real.write_text("#!/bin/sh\necho real\n")
    assert live_session._executable_is_discovery_stub(str(real)) is False
    assert live_session._executable_is_discovery_stub("") is False
    assert (
        live_session._executable_is_discovery_stub(str(tmp_path / "nope"))
        is False
    )


def test_conformance_matrix_is_derived_from_declarations(declared_server):
    """Adding a program to ChemSmart must not require editing the session."""

    from chemsmart.settings.capabilities import EXECUTABLE_PROGRAMS

    programs = set(live_session._conformance_programs())
    assert programs <= set(EXECUTABLE_PROGRAMS)
    # A program declaring no core stage contributes no conformance work.
    assert "nciplot" not in programs
    # Per-engine capability is honoured: GPU PySCF previews no excited state.
    assert "td" in live_session._conformance_jobtypes("pyscf", "cpu")
    assert "td" not in live_session._conformance_jobtypes("pyscf", "gpu")


def test_project_requiring_programs_get_a_fixture(declared_server):
    for program in ("orca", "gaussian", "pyscf"):
        assert live_session._conformance_project_sections(program) is not None
    assert live_session._conformance_project_sections("xtb") is None


def test_a_fixture_declares_every_stage_it_claims_to_cover(declared_server):
    """A stage with no section makes the loader return ``None`` and the CLI crash."""

    for program in ("gaussian", "pyscf"):
        sections = live_session._conformance_project_sections(program)
        for stage in live_session._conformance_jobtypes(program, "cpu"):
            if stage in live_session._ROUTE_PROGRAM_STAGE_SECTIONS or (
                program == "pyscf"
            ):
                assert stage in sections, f"{program} fixture omits {stage}"


def test_a_fixture_carries_only_project_owned_keys(declared_server):
    """An unknown project key raises, so one bad key breaks every stage."""

    from chemsmart.jobs.gaussian.settings import GaussianJobSettings

    allowed = set(GaussianJobSettings.default().__dict__)
    for settings in live_session._conformance_project_sections(
        "gaussian"
    ).values():
        assert set(settings) <= allowed


def test_environment_observation_does_not_depend_on_the_process_path(
    declared_server, monkeypatch
):
    """A program on PATH but undeclared is not part of ChemSmart's environment."""

    monkeypatch.setenv("PATH", os.devnull)
    records = [
        item
        for item in live_session._observe_environments()[2]
        if item.get("record_kind") == "program_environment"
    ]
    assert {item["program"] for item in records} >= {"orca", "gaussian"}
