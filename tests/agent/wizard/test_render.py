from pathlib import Path

import yaml

from chemsmart.agent.wizard import (
    CondaEnvSurvey,
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
from chemsmart.agent.wizard.parsers import parse_sge_qconf_sq
from chemsmart.agent.wizard.probe import ProbeResult

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sge"


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
    queue = parse_sge_qconf_sq(
        (FIXTURE_DIR / "qconf_sq_20core.q.txt").read_text()
    )
    return ScheduleSurvey(
        scheduler="SGE",
        submit_command="qsub",
        queues=[
            QueueFacts(
                name=queue.name,
                default=queue.default,
                max_walltime_hours=queue.max_walltime_hours,
                default_walltime_hours=queue.default_walltime_hours,
                default_mem_gb=queue.default_mem_gb,
                default_cores=queue.default_cores,
                gpus_per_node=queue.gpus_per_node,
                enabled=queue.enabled,
                started=queue.started,
                slots_total=100,
            )
        ],
        chosen_queue="20core.q",
        evidence={"qconf -sq 20core.q": "parsed"},
    )


def _software_survey():
    conda_base = (FIXTURE_DIR / "conda_info_base.txt").read_text().strip()
    conda_env = (FIXTURE_DIR / "conda_prefix_env.txt").read_text().strip()
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
        conda=CondaEnvSurvey(
            base=conda_base,
            env_path=conda_env,
            env_name="chemsmart",
        ),
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
    assert parsed["SERVER"]["QUEUE_NAME"] == "20core.q"
    assert parsed["SERVER"]["SUBMIT_COMMAND"] == "qsub"
    assert plan.server_block["PROJECT"] == "chem-123"
    assert parsed["SERVER"]["PROJECT"] == "chem-123"
    assert plan.program_blocks["GAUSSIAN"]["SCRATCH"] is False
    assert parsed["GAUSSIAN"]["SCRATCH"] is False
    assert plan.program_blocks["GAUSSIAN"]["EXEFOLDER"] == "/apps/gaussian"
    assert plan.program_blocks["GAUSSIAN"]["MODULES"] == ""
    assert (
        plan.program_blocks["GAUSSIAN"]["CONDA_ENV"]
        == "source /opt/conda/etc/profile.d/conda.sh\n"
        "conda activate chemsmart"
    )
    assert (
        plan.server_block["EXTRA_COMMANDS"]
        == "export PATH=/opt/chemsmart/bin:$PATH"
    )
    assert any("SCRATCH set to False" in note for note in plan.notes)


def test_render_server_yaml_keeps_program_scratch_true_when_writable():
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
            project="chem-123",
            source="sacctmgr",
            candidates=["chem-123"],
        ),
    )

    parsed = yaml.safe_load(plan.text)

    assert all(
        block["SCRATCH"] is True for block in plan.program_blocks.values()
    )
    assert all(
        parsed[block_name]["SCRATCH"] is True
        for block_name in ["GAUSSIAN", "ORCA", "NCIPLOT"]
    )
    assert not any("SCRATCH set to False" in note for note in plan.notes)


def test_render_server_yaml_sets_all_program_scratch_false_when_unwritable():
    plan = render_server_yaml(
        Topology(mode="B", host="cluster", evidence=[]),
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
    )

    parsed = yaml.safe_load(plan.text)

    assert all(
        block["SCRATCH"] is False for block in plan.program_blocks.values()
    )
    assert all(
        parsed[block_name]["SCRATCH"] is False
        for block_name in ["GAUSSIAN", "ORCA", "NCIPLOT"]
    )
    assert any("SCRATCH set to False" in note for note in plan.notes)


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


def test_render_server_yaml_omits_conda_env_without_conda_base():
    survey = _software_survey()
    survey = SoftwareSurvey(
        module_system=survey.module_system,
        programs=survey.programs,
        conda=CondaEnvSurvey(base=None, env_path=None, env_name=None),
    )

    plan = render_server_yaml(
        Topology(mode="B", host="cluster", evidence=[]),
        _schedule_survey(),
        survey,
        ScratchFinding(
            path="/scratch/user",
            source="env:SCRATCH",
            writable=True,
            candidates=[("env:SCRATCH", "/scratch/user")],
        ),
        ProjectFinding(
            project="chem-123",
            source="sacctmgr",
            candidates=["chem-123"],
        ),
    )

    parsed = yaml.safe_load(plan.text)

    assert "CONDA_ENV" not in parsed["GAUSSIAN"]
