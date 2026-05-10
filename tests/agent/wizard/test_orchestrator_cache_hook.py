from __future__ import annotations

import logging

from chemsmart.agent.wizard.cache import CacheEntry, load_cache, write_cache
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
  HOST: localhost
GAUSSIAN: {}
"""


def _stub_surveys():
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
    validation = ValidationResult(
        ok=True,
        errors=[],
        parsed={"SERVER": {"SCHEDULER": "SLURM"}},
    )
    plan = ServerYamlPlan(
        text=YAML_TEXT,
        server_block={"SCHEDULER": "SLURM"},
        program_blocks={"GAUSSIAN": {}},
        notes=[],
    )
    return topology, schedule, software, scratch, project, validation, plan


def _patch_probe_pipeline(monkeypatch):
    topology, schedule, software, scratch, project, validation, plan = (
        _stub_surveys()
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.detect_topology",
        lambda runner, ssh_host_hint=None: topology,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.run_schedule_survey",
        lambda runner, topo: schedule,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.run_software_survey",
        lambda runner, topo: software,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.discover_scratch",
        lambda runner, topo: scratch,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.discover_project",
        lambda runner, topo, scheduler: project,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.render_server_yaml",
        lambda *args, **kwargs: plan,
    )
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.validate_server_yaml",
        lambda yaml_text, server_name: validation,
    )


def test_run_wizard_write_updates_yaml_and_cache(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    _patch_probe_pipeline(monkeypatch)

    def fake_refresh_cache(runner, server_name, force=False):
        assert force is True
        entry = CacheEntry(
            server_name=server_name,
            host=None,
            mode="local",
            scheduler="SLURM",
            probed_at="2026-05-10T00:00:00Z",
            source_commands={"scheduler": "sinfo --json"},
            partitions=[],
            node_summary={"selected_queue": "debug"},
            program_candidates={},
            status="fresh",
            last_error=None,
        )
        write_cache(entry)
        return entry

    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.refresh_cache",
        fake_refresh_cache,
    )

    outcome = run_wizard(StubRunner(), server_name="perlmutter", write=True)

    yaml_path = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert outcome.written_path == str(yaml_path)
    assert yaml_path.read_text(encoding="utf-8") == YAML_TEXT
    assert load_cache("perlmutter") is not None


def test_cache_failure_does_not_fail_wizard(monkeypatch, tmp_path, caplog):
    monkeypatch.setenv("HOME", str(tmp_path))
    _patch_probe_pipeline(monkeypatch)
    monkeypatch.setattr(
        "chemsmart.agent.wizard.orchestrator.refresh_cache",
        lambda runner, server_name, force=False: (_ for _ in ()).throw(
            RuntimeError("cache down")
        ),
    )

    with caplog.at_level(logging.WARNING):
        outcome = run_wizard(
            StubRunner(), server_name="perlmutter", write=True
        )

    yaml_path = tmp_path / ".chemsmart" / "server" / "perlmutter.yaml"
    assert outcome.written_path == str(yaml_path)
    assert yaml_path.exists()
    assert "Wizard cache refresh failed" in caplog.text
