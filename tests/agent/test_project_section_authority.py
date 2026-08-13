from __future__ import annotations

import pytest
import yaml

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.projects import (
    project_document,
    project_effective_section_settings,
    project_section_application_observation,
    render_project_yaml,
)
from chemsmart.jobs.settings import (
    molecular_project_section_names,
    molecular_project_section_sources,
    read_molecular_job_yaml,
)
from chemsmart.settings.capabilities import (
    PROGRAM_CAPABILITIES,
    loader_project_section_names,
)


def test_agent_projects_use_molecular_loader_section_order():
    sections = {
        "gas": {"basis": "def2-SVP", "functional": "PBE0"},
        "sp": {"basis": "def2-TZVP"},
    }
    document = project_document(program="orca", sections=sections)

    assert molecular_project_section_sources(
        sections, program="orca", jobtype="sp"
    ) == ("gas", "sp")
    assert dict(
        project_effective_section_settings(document, jobtype="sp")
    ) == {"basis": "def2-TZVP", "functional": "PBE0"}
    observation = project_section_application_observation(
        document, jobtype="sp"
    )
    assert observation["used_sections"] == ("gas", "sp")


def test_pyscf_migration_section_semantics_come_from_loader_constants():
    document = project_document(
        program="pyscf",
        sections={
            "gas": {
                "basis": "def2-SVP",
                "functional": "B3LYP",
            }
        },
    )

    assert dict(
        project_effective_section_settings(document, jobtype="opt")
    ) == {"basis": "def2-SVP", "functional": "B3LYP"}
    assert project_section_application_observation(document, jobtype="opt")[
        "used_sections"
    ] == ("gas",)


@pytest.mark.parametrize("program", ("gaussian", "orca"))
def test_every_declared_molecular_section_is_readable_by_loader(
    tmp_path, program
):
    declared = molecular_project_section_names(program)
    assert declared == loader_project_section_names(program)
    assert declared == PROGRAM_CAPABILITIES[program].project_section_names

    for section in declared:
        path = tmp_path / f"{program}-{section}.yaml"
        path.write_text(
            yaml.safe_dump({section: {}}),
            encoding="utf-8",
        )
        read_molecular_job_yaml(path, program=program)


@pytest.mark.parametrize("program", ("gaussian", "orca"))
def test_undeclared_molecular_section_is_rejected_by_loader_and_agent(
    tmp_path, program
):
    path = tmp_path / f"{program}-unknown.yaml"
    path.write_text("invented_stage: {}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unknown .* project section"):
        read_molecular_job_yaml(path, program=program)

    document = project_document(
        program=program,
        sections={"invented_stage": {}},
    )
    with pytest.raises(ContractError, match="has no section"):
        render_project_yaml(document)
