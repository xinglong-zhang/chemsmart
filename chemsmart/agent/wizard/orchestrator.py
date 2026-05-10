"""High-level orchestration for the wizard probe flow."""

from __future__ import annotations

from dataclasses import dataclass

from chemsmart.agent.wizard.project import discover_project
from chemsmart.agent.wizard.render import ServerYamlPlan, render_server_yaml
from chemsmart.agent.wizard.scratch import discover_scratch
from chemsmart.agent.wizard.software import run_software_survey
from chemsmart.agent.wizard.survey import run_schedule_survey
from chemsmart.agent.wizard.topology import detect_topology
from chemsmart.agent.wizard.validate import (
    ValidationResult,
    validate_server_yaml,
)
from chemsmart.agent.wizard.write import write_server_yaml


@dataclass(frozen=True)
class WizardOutcome:
    plan: ServerYamlPlan
    validation: ValidationResult
    surveys: dict
    written_path: str | None


def run_wizard(
    runner,
    server_name: str,
    ssh_host_hint: str | None = None,
    write: bool = False,
    overwrite: bool = False,
) -> WizardOutcome:
    """Probe the target environment and build a validated server YAML plan."""

    topology = detect_topology(runner, ssh_host_hint=ssh_host_hint)
    schedule_survey = run_schedule_survey(runner, topology)
    software_survey = run_software_survey(runner, topology)
    scratch_finding = discover_scratch(runner, topology)
    project_finding = discover_project(
        runner,
        topology,
        schedule_survey.scheduler,
    )
    plan = render_server_yaml(
        topology,
        schedule_survey,
        software_survey,
        scratch_finding,
        project_finding,
        runner=runner,
    )
    validation = validate_server_yaml(plan.text, server_name=server_name)

    written_path = None
    if write and validation.ok:
        written_path = write_server_yaml(
            name=server_name,
            yaml_text=plan.text,
            overwrite=overwrite,
        )

    return WizardOutcome(
        plan=plan,
        validation=validation,
        surveys={
            "topology": topology,
            "schedule": schedule_survey,
            "software": software_survey,
            "scratch": scratch_finding,
            "project": project_finding,
        },
        written_path=written_path,
    )
