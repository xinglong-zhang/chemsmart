import yaml

from chemsmart.agent.wizard import (
    ModuleSystem,
    ProgramFinding,
    ProjectFinding,
    QueueFacts,
    ScheduleSurvey,
    ScratchFinding,
    SoftwareSurvey,
    Topology,
    render_server_yaml,
)
from chemsmart.agent.wizard.probe import ProbeResult


class StubRunner:
    def __init__(self, local_results=None):
        self.local_results = local_results or {}

    def run_local(self, command, timeout_s=15):
        return self.local_results.get(
            tuple(command),
            ProbeResult(
                command=" ".join(command),
                mode="local",
                host=None,
                returncode=1,
                stdout="",
                stderr="",
                duration_s=0.0,
                truncated=False,
            ),
        )


def _result(command, stdout="", returncode=0):
    return ProbeResult(
        command=command,
        mode="local",
        host=None,
        returncode=returncode,
        stdout=stdout,
        stderr="",
        duration_s=0.0,
        truncated=False,
    )


def _schedule_survey():
    return ScheduleSurvey(
        scheduler="SLURM",
        submit_command="sbatch",
        queues=[
            QueueFacts(
                name="debug",
                default=True,
                max_walltime_hours=48,
                default_walltime_hours=8,
                default_mem_gb=128,
                default_cores=32,
                gpus_per_node=0,
                enabled=True,
                started=True,
            )
        ],
        chosen_queue="debug",
        evidence={"scontrol show partition --oneliner": "parsed"},
    )


def _software_survey():
    return SoftwareSurvey(
        module_system=ModuleSystem(kind="lmod", version="Lmod 8.7.49"),
        programs={
            "gaussian": ProgramFinding(
                program="gaussian",
                exefolder="/apps/gaussian",
                source="path",
                module_candidates=[],
                on_path=True,
            ),
            "orca": ProgramFinding(
                program="orca",
                exefolder=None,
                source="module",
                module_candidates=["orca/6.0.1"],
                on_path=False,
            ),
            "nciplot": ProgramFinding(
                program="nciplot",
                exefolder="/apps/nciplot",
                source="path",
                module_candidates=[],
                on_path=True,
            ),
        },
        conda_base="/opt/conda",
        conda_env="/opt/conda/envs/chemsmart",
    )


def test_render_server_yaml_mode_a_applies_required_decisions():
    plan = render_server_yaml(
        Topology(mode="A", host="login.cluster.example", evidence=[]),
        _schedule_survey(),
        _software_survey(),
        ScratchFinding(
            path="/scratch/user",
            source="env:SCRATCH",
            writable=False,
            candidates=[("env:SCRATCH", "/scratch/user")],
        ),
        ProjectFinding(
            project="chem-123",
            source="sacctmgr",
            candidates=["chem-123"],
        ),
        runner=StubRunner(
            local_results={
                ("printenv", "CHEMSMART_BIN"): _result(
                    "printenv CHEMSMART_BIN",
                    "/opt/chemsmart/bin\n",
                )
            }
        ),
    )

    parsed = yaml.safe_load(plan.text)

    assert plan.server_block["HOST"] == "localhost"
    assert parsed["SERVER"]["HOST"] == "localhost"
    assert plan.server_block["PROJECT"] == "chem-123"
    assert parsed["SERVER"]["PROJECT"] == "chem-123"
    assert plan.program_blocks["GAUSSIAN"]["SCRATCH"] is True
    assert plan.program_blocks["GAUSSIAN"]["EXEFOLDER"] == "/apps/gaussian"
    assert plan.program_blocks["GAUSSIAN"]["MODULES"] == ""
    assert (
        plan.server_block["EXTRA_COMMANDS"]
        == "export PATH=/opt/chemsmart/bin:$PATH"
    )
    assert any("SCRATCH kept True" in note for note in plan.notes)


def test_render_server_yaml_comments_groups_project_in_text():
    plan = render_server_yaml(
        Topology(mode="B", host="cluster", evidence=[]),
        _schedule_survey(),
        _software_survey(),
        ScratchFinding(
            path="/scratch/user",
            source="env:SCRATCH",
            writable=True,
            candidates=[("env:SCRATCH", "/scratch/user")],
        ),
        ProjectFinding(
            project="chem-group",
            source="groups",
            candidates=["chem-group"],
        ),
    )

    parsed = yaml.safe_load(plan.text)

    assert "## PROJECT: chem-group" in plan.text
    assert "PROJECT" not in parsed["SERVER"]
