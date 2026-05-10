from __future__ import annotations

from chemsmart.agent.wizard.orchestrator import run_wizard
from chemsmart.agent.wizard.project import ProjectFinding
from chemsmart.agent.wizard.render import ServerYamlPlan
from chemsmart.agent.wizard.scratch import ScratchFinding
from chemsmart.agent.wizard.software import (
    ModuleSystem,
    ProgramFinding,
    SoftwareSurvey,
)
from chemsmart.agent.wizard.survey import ScheduleSurvey
from chemsmart.agent.wizard.topology import Topology
from chemsmart.agent.wizard.validate import ValidationResult


class StubRunner:
    pass


YAML_TEXT = """SERVER:
  SCHEDULER: SLURM
GAUSSIAN: {}
"""


def test_run_wizard_happy_path_writes(monkeypatch, tmp_path):
    calls: list[str] = []
    topology = Topology(mode="A", host="localhost", evidence=["env:SLURM"])
    schedule = ScheduleSurvey(
        scheduler="SLURM",
        submit_command="sbatch",
        queues=[],
        chosen_queue="debug",
        evidence={"sinfo --json": "parsed"},
    )
    software = SoftwareSurvey(
        module_system=ModuleSystem(kind="lmod", version="8.7"),
        programs={
            "gaussian": ProgramFinding(
                program="gaussian",
                exefolder="/apps/gaussian",
                source="path",
                module_candidates=[],
                on_path=True,
            )
        },
        conda_base=None,
        conda_env=None,
    )
    scratch = ScratchFinding(
        path="/scratch/user",
        source="env:SCRATCH",
        writable=True,
        candidates=[("env:SCRATCH", "/scratch/user")],
    )
    project = ProjectFinding(
        project="chem-123",
        source="groups",
        candidates=["chem-123"],
    )
    plan = ServerYamlPlan(
        text=YAML_TEXT,
        server_block={"SCHEDULER": "SLURM"},
        program_blocks={"GAUSSIAN": {}},
        notes=[],
    )
    validation = ValidationResult(
        ok=True,
        errors=[],
        parsed={"SERVER": {"SCHEDULER": "SLURM"}},
    )

    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.detect_topology",
        lambda runner, ssh_host_hint=None: calls.append("detect_topology")
        or topology,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.run_schedule_survey",
        lambda runner, topo: calls.append("run_schedule_survey") or schedule,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.run_software_survey",
        lambda runner, topo: calls.append("run_software_survey") or software,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.discover_scratch",
        lambda runner, topo: calls.append("discover_scratch") or scratch,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.discover_project",
        lambda runner, topo, scheduler: calls.append("discover_project")
        or project,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.render_server_yaml",
        lambda *args, **kwargs: calls.append("render_server_yaml") or plan,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.validate_server_yaml",
        lambda yaml_text, server_name: calls.append("validate_server_yaml")
        or validation,
    )
    monkeypatch.setenv("HOME", str(tmp_path))

    outcome = run_wizard(StubRunner(), server_name="perlmutter", write=True)

    assert calls == [
        "detect_topology",
        "run_schedule_survey",
        "run_software_survey",
        "discover_scratch",
        "discover_project",
        "render_server_yaml",
        "validate_server_yaml",
    ]
    expected = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert expected.read_text(encoding="utf-8") == plan.text
    assert outcome.written_path == str(expected)
    assert outcome.plan == plan
    assert outcome.validation == validation
    assert outcome.surveys["topology"] == topology
    assert outcome.surveys["schedule"] == schedule
