"""Render wizard findings into canonical server YAML.

Canonical wizard schema mirrors the normalized server templates under
``chemsmart/settings/templates/.chemsmart/server``:

- ``SERVER`` uses ``SCRATCH_DIR`` at the server level, never ``SCRATCH``.
- ``SERVER`` fields are ordered as ``SCHEDULER``, ``QUEUE_NAME``,
  ``NUM_HOURS``, ``MEM_GB``, ``NUM_CORES``, ``NUM_GPUS``,
  ``NUM_THREADS``, ``SUBMIT_COMMAND``, ``HOST``, optional ``PROJECT``,
  ``SCRATCH_DIR``, ``USE_HOSTS``, and ``EXTRA_COMMANDS``.
- Program blocks use ``EXEFOLDER``, ``LOCAL_RUN``, ``SCRATCH``,
  ``CONDA_ENV``, ``MODULES``, optional ``SCRIPTS``, and ``ENVARS``.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass

import yaml

from chemsmart.agent.wizard.normalize import (
    FALLBACK_NUM_HOURS,
    SCHEDULER_SUBMIT,
    normalize_resources,
    normalize_walltime,
)

_PROGRAM_NAME_MAP = {
    "gaussian": "GAUSSIAN",
    "orca": "ORCA",
    "nciplot": "NCIPLOT",
}
_LOCAL_RUN_DEFAULTS = {
    "GAUSSIAN": True,
    "ORCA": False,
    "NCIPLOT": False,
}
_PROJECT_COMMENT_KEY = "__WIZARD_PROJECT_COMMENT__"
_EXTRA_COMMANDS_PLACEHOLDER = (
    "#extra commands to activate chemsmart environment in submission script"
)


@dataclass(frozen=True)
class ServerYamlPlan:
    text: str
    server_block: dict
    program_blocks: dict[str, dict]
    notes: list[str]


def render_server_yaml(
    topology,
    schedule_survey,
    software_survey,
    scratch_finding,
    project_finding,
    runner=None,
) -> ServerYamlPlan:
    """Render a candidate server YAML from wizard findings."""

    notes: list[str] = []
    queue = _find_selected_queue(schedule_survey)
    resources = (
        normalize_resources(queue)
        if queue is not None
        else {
            "mem_gb": None,
            "cores": None,
            "gpus": 0,
            "threads": None,
        }
    )
    scheduler = schedule_survey.scheduler
    server_block: dict[str, object] = {
        "SCHEDULER": scheduler,
        "QUEUE_NAME": schedule_survey.chosen_queue,
        "NUM_HOURS": (
            normalize_walltime(queue)
            if queue is not None
            else FALLBACK_NUM_HOURS
        ),
        "MEM_GB": resources["mem_gb"],
        "NUM_CORES": resources["cores"],
        "NUM_GPUS": resources["gpus"],
        "NUM_THREADS": resources["threads"],
        "SUBMIT_COMMAND": schedule_survey.submit_command
        or SCHEDULER_SUBMIT.get(scheduler),
        "HOST": "localhost" if topology.mode == "A" else topology.host,
        "SCRATCH_DIR": scratch_finding.path,
        "USE_HOSTS": True,
        "EXTRA_COMMANDS": _render_extra_commands(
            topology=topology,
            runner=runner,
            notes=notes,
        ),
    }
    if (
        project_finding.project is not None
        and project_finding.source != "groups"
    ):
        server_block["PROJECT"] = project_finding.project

    if (
        project_finding.project is not None
        and project_finding.source == "groups"
    ):
        notes.append(
            "Rendered PROJECT as a comment because it came from groups."
        )

    if not scratch_finding.writable:
        notes.append(
            "Program SCRATCH kept True despite scratch probe not confirming "
            "writability."
        )

    program_blocks = _render_program_blocks(
        software_survey=software_survey,
        scratch_finding=scratch_finding,
        notes=notes,
    )

    yaml_doc = {"SERVER": dict(server_block), **program_blocks}
    if (
        project_finding.project is not None
        and project_finding.source == "groups"
    ):
        yaml_doc["SERVER"][
            _PROJECT_COMMENT_KEY
        ] = f"## PROJECT: {project_finding.project}"

    text = yaml.safe_dump(
        yaml_doc,
        sort_keys=False,
        default_flow_style=False,
    )
    text = _inject_project_comment(text)
    return ServerYamlPlan(
        text=text,
        server_block=server_block,
        program_blocks=program_blocks,
        notes=notes,
    )


def _find_selected_queue(schedule_survey):
    if schedule_survey.chosen_queue is None:
        return None
    for queue in schedule_survey.queues:
        if queue.name == schedule_survey.chosen_queue:
            return queue
    return None


def _render_program_blocks(
    software_survey,
    scratch_finding,
    notes: list[str],
) -> dict[str, dict]:
    program_blocks: dict[str, dict] = {}
    for program_key in ["gaussian", "orca", "nciplot"]:
        block_name = _PROGRAM_NAME_MAP[program_key]
        finding = software_survey.programs.get(program_key)
        if finding is None:
            notes.append(f"No {block_name} finding was available to render.")
            continue
        program_blocks[block_name] = _render_program_block(
            block_name=block_name,
            finding=finding,
            conda_env=software_survey.conda_env,
            scratch_dir=scratch_finding.path,
            notes=notes,
        )
    return program_blocks


def _render_program_block(
    block_name: str,
    finding,
    conda_env: str | None,
    scratch_dir: str | None,
    notes: list[str],
) -> dict[str, object]:
    block: dict[str, object] = {
        "EXEFOLDER": finding.exefolder,
        "LOCAL_RUN": _LOCAL_RUN_DEFAULTS[block_name],
        "SCRATCH": True,
        "CONDA_ENV": conda_env,
        "MODULES": "" if finding.on_path else _render_modules(finding),
        "ENVARS": _render_envars(block_name, finding.exefolder, scratch_dir),
    }
    if block_name == "GAUSSIAN":
        block["SCRIPTS"] = ""

    if finding.source == "module" and len(finding.module_candidates) > 1:
        notes.append(
            f"{block_name} selected module candidate "
            f"{finding.module_candidates[0]!r} from multiple options."
        )
    return block


def _render_modules(finding) -> str:
    if not finding.module_candidates:
        return ""
    lines = ["module purge", f"module load {finding.module_candidates[0]}"]
    return "\n".join(lines)


def _render_envars(
    block_name: str,
    exefolder: str | None,
    scratch_dir: str | None,
) -> str:
    lines: list[str] = []
    if scratch_dir:
        lines.append(f"export SCRATCH={scratch_dir}")
    if block_name == "GAUSSIAN" and exefolder:
        lines.append(f"export GAUSS_EXEDIR={exefolder}")
        lines.append(f"export g16root={exefolder}")
    if block_name == "NCIPLOT" and exefolder:
        lines.append(f"export NCIPLOT_HOME={exefolder}")
    return "\n".join(lines)


def _render_extra_commands(topology, runner, notes: list[str]) -> str:
    if (
        topology.mode != "A"
        or runner is None
        or not hasattr(runner, "run_local")
    ):
        notes.append(
            "Left EXTRA_COMMANDS as a placeholder because CHEMSMART_BIN "
            "was not probed locally."
        )
        return _EXTRA_COMMANDS_PLACEHOLDER

    result = runner.run_local(["printf", "%s\\n", "$CHEMSMART_BIN"])
    if result.returncode != 0:
        notes.append(
            "Left EXTRA_COMMANDS as a placeholder because CHEMSMART_BIN "
            "probe failed."
        )
        return _EXTRA_COMMANDS_PLACEHOLDER

    chemsmart_bin = _normalize_shell_value(result.stdout)
    if not chemsmart_bin:
        notes.append(
            "Left EXTRA_COMMANDS as a placeholder because CHEMSMART_BIN "
            "was empty."
        )
        return _EXTRA_COMMANDS_PLACEHOLDER

    return f"export PATH={shlex.quote(chemsmart_bin)}:$PATH"


def _normalize_shell_value(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped or stripped.startswith("$"):
        return None
    return stripped


def _inject_project_comment(text: str) -> str:
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(f"  {_PROJECT_COMMENT_KEY}:"):
            comment = line.split(": ", 1)[1].strip().strip("'\"")
            lines.append(f"  {comment}")
            continue
        lines.append(line)
    return "\n".join(lines) + "\n"
