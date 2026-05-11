"""Probe primitives for wizard transport-aware discovery."""

from __future__ import annotations

import re
import shlex
import string
import subprocess
import time
from dataclasses import dataclass
from typing import Callable

from chemsmart.agent.wizard.scheduler_env import build_scheduler_env

TRUNCATION_MARKER = "…[truncated]"
MAX_OUTPUT_BYTES = 1024 * 1024
_IDENTIFIER_PATTERN = re.compile(r"^[\w\-\.]+$")
_ENV_VAR_PATTERN = re.compile(r"^[A-Z_][A-Z0-9_]*$")
_PATH_PATTERN = re.compile(r"^[/~][\w\-./@:+]*$")
_FORBIDDEN_SLOT_CHARS = set(";&|`$(){}<>!\\'\"*\n\r\t")

SlotValidator = Callable[[str], None]


class ProbeError(Exception):
    """Raised when a probe request is invalid."""


@dataclass(frozen=True)
class ProbeSpec:
    """A closed, validated probe command template."""

    template_id: str
    argv_template: tuple[str, ...]
    slot_validators: dict[str, SlotValidator]


@dataclass(frozen=True)
class ProbeResult:
    command: str
    mode: str
    host: str | None
    returncode: int
    stdout: str
    stderr: str
    duration_s: float
    truncated: bool


ALL_PROBE_SPECS: dict[str, ProbeSpec] = {
    "common.printenv_all": ProbeSpec(
        template_id="common.printenv_all",
        argv_template=("printenv",),
        slot_validators={},
    ),
    "common.printenv_var": ProbeSpec(
        template_id="common.printenv_var",
        argv_template=("printenv", "{env_name}"),
        slot_validators={
            "env_name": lambda value: validate_env_var_name(value)
        },
    ),
    "topology.sinfo": ProbeSpec(
        template_id="topology.sinfo",
        argv_template=("sinfo",),
        slot_validators={},
    ),
    "topology.qstat": ProbeSpec(
        template_id="topology.qstat",
        argv_template=("qstat",),
        slot_validators={},
    ),
    "topology.bqueues": ProbeSpec(
        template_id="topology.bqueues",
        argv_template=("bqueues",),
        slot_validators={},
    ),
    "topology.qconf": ProbeSpec(
        template_id="topology.qconf",
        argv_template=("qconf",),
        slot_validators={},
    ),
    "survey.slurm.sinfo_json": ProbeSpec(
        template_id="survey.slurm.sinfo_json",
        argv_template=("sinfo", "--json"),
        slot_validators={},
    ),
    "survey.slurm.scontrol_partition": ProbeSpec(
        template_id="survey.slurm.scontrol_partition",
        argv_template=(
            "scontrol",
            "show",
            "partition",
            "--oneliner",
        ),
        slot_validators={},
    ),
    "survey.pbs.qstat_json": ProbeSpec(
        template_id="survey.pbs.qstat_json",
        argv_template=("qstat", "-Q", "-f", "-F", "json"),
        slot_validators={},
    ),
    "survey.pbs.qstat_text": ProbeSpec(
        template_id="survey.pbs.qstat_text",
        argv_template=("qstat", "-Q", "-f"),
        slot_validators={},
    ),
    "survey.lsf.bqueues_json": ProbeSpec(
        template_id="survey.lsf.bqueues_json",
        argv_template=(
            "bqueues",
            "-o",
            "queue_name status runlimit memlimit prolimit",
            "-json",
        ),
        slot_validators={},
    ),
    "survey.lsf.bqueues_text": ProbeSpec(
        template_id="survey.lsf.bqueues_text",
        argv_template=("bqueues", "-l"),
        slot_validators={},
    ),
    "survey.sge.qconf_sql": ProbeSpec(
        template_id="survey.sge.qconf_sql",
        argv_template=("qconf", "-sql"),
        slot_validators={},
    ),
    "survey.sge.qconf_sq": ProbeSpec(
        template_id="survey.sge.qconf_sq",
        argv_template=("qconf", "-sq", "{queue_name}"),
        slot_validators={
            "queue_name": lambda value: validate_identifier(value)
        },
    ),
    "survey.sge.qhost": ProbeSpec(
        template_id="survey.sge.qhost",
        argv_template=("qhost",),
        slot_validators={},
    ),
    "software.type_module": ProbeSpec(
        template_id="software.type_module",
        argv_template=("type", "module"),
        slot_validators={},
    ),
    "software.which_module": ProbeSpec(
        template_id="software.which_module",
        argv_template=("which", "module"),
        slot_validators={},
    ),
    "software.module_version": ProbeSpec(
        template_id="software.module_version",
        argv_template=("module", "--version"),
        slot_validators={},
    ),
    "software.command_v": ProbeSpec(
        template_id="software.command_v",
        argv_template=("command", "-v", "{exe_name}"),
        slot_validators={"exe_name": lambda value: validate_identifier(value)},
    ),
    "software.which_exe": ProbeSpec(
        template_id="software.which_exe",
        argv_template=("which", "{exe_name}"),
        slot_validators={"exe_name": lambda value: validate_identifier(value)},
    ),
    "software.readlink": ProbeSpec(
        template_id="software.readlink",
        argv_template=("readlink", "-f", "{path}"),
        slot_validators={"path": lambda value: validate_path_slot(value)},
    ),
    "software.module_avail": ProbeSpec(
        template_id="software.module_avail",
        argv_template=("module", "-t", "avail"),
        slot_validators={},
    ),
    "software.conda_base": ProbeSpec(
        template_id="software.conda_base",
        argv_template=("conda", "info", "--base"),
        slot_validators={},
    ),
    "scratch.test_dir_writable": ProbeSpec(
        template_id="scratch.test_dir_writable",
        argv_template=("test", "-d", "{path}", "-a", "-w", "{path}"),
        slot_validators={"path": lambda value: validate_path_slot(value)},
    ),
    "scratch.test_writable": ProbeSpec(
        template_id="scratch.test_writable",
        argv_template=("test", "-w", "{path}"),
        slot_validators={"path": lambda value: validate_path_slot(value)},
    ),
    "project.sacctmgr_show_user": ProbeSpec(
        template_id="project.sacctmgr_show_user",
        argv_template=(
            "sacctmgr",
            "-n",
            "-p",
            "show",
            "user",
            "{user}",
            "format=DefaultAccount,Account",
        ),
        slot_validators={"user": lambda value: validate_identifier(value)},
    ),
    "project.groups": ProbeSpec(
        template_id="project.groups",
        argv_template=("groups",),
        slot_validators={},
    ),
}


class ProbeRunner:
    """Run closed-catalog read-only probes locally or via SSH."""

    supports_probe_specs = True

    @classmethod
    def run_local(
        cls,
        spec: ProbeSpec | object,
        timeout_s: int = 15,
        **slots: str,
    ) -> ProbeResult:
        if not isinstance(spec, ProbeSpec):
            raise ProbeError("raw command string not allowed; use ProbeSpec")
        argv = resolve_probe_argv(spec, **slots)
        display_command = shlex.join(argv)
        return cls._run(
            args=["bash", "-lc", display_command],
            command_text=display_command,
            mode="local",
            host=None,
            timeout_s=timeout_s,
            env=build_scheduler_env(),
        )

    @classmethod
    def run_ssh(
        cls,
        host: str,
        spec: ProbeSpec | object,
        timeout_s: int = 15,
        **slots: str,
    ) -> ProbeResult:
        if not isinstance(spec, ProbeSpec):
            raise ProbeError("raw command string not allowed; use ProbeSpec")
        argv = resolve_probe_argv(spec, **slots)
        display_command = shlex.join(argv)
        ssh_command = [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ClearAllForwardings=yes",
            "-o",
            "ConnectTimeout=10",
            "-o",
            "StrictHostKeyChecking=yes",
            host,
            f"bash -lc {shlex.quote(display_command)}",
        ]
        return cls._run(
            args=ssh_command,
            command_text=display_command,
            mode="ssh",
            host=host,
            timeout_s=timeout_s,
        )

    @classmethod
    def _run(
        cls,
        args: list[str],
        command_text: str,
        mode: str,
        host: str | None,
        timeout_s: int,
        env: dict[str, str] | None = None,
    ) -> ProbeResult:
        start = time.monotonic()
        try:
            result = subprocess.run(
                args,
                capture_output=True,
                text=True,
                timeout=timeout_s,
                check=False,
                env=env,
            )
            stdout = cls._coerce_text(result.stdout)
            stderr = cls._coerce_text(result.stderr)
            returncode = result.returncode
        except subprocess.TimeoutExpired as exc:
            stdout = cls._coerce_text(getattr(exc, "stdout", None))
            stderr = cls._coerce_text(getattr(exc, "stderr", None))
            timeout_marker = "[probe-timeout]"
            stderr = (
                f"{stderr.rstrip()}\n{timeout_marker}"
                if stderr.strip()
                else timeout_marker
            )
            returncode = 124
        except FileNotFoundError as exc:
            stdout = ""
            stderr = str(exc)
            returncode = 127
        except OSError as exc:
            stdout = ""
            stderr = str(exc)
            returncode = 1

        stdout, stdout_truncated = cls._truncate_output(stdout)
        stderr, stderr_truncated = cls._truncate_output(stderr)
        duration_s = time.monotonic() - start
        return ProbeResult(
            command=command_text,
            mode=mode,
            host=host,
            returncode=returncode,
            stdout=stdout,
            stderr=stderr,
            duration_s=duration_s,
            truncated=stdout_truncated or stderr_truncated,
        )

    @staticmethod
    def _coerce_text(value: str | bytes | None) -> str:
        if value is None:
            return ""
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        return value

    @staticmethod
    def _truncate_output(text: str) -> tuple[str, bool]:
        if len(text) <= MAX_OUTPUT_BYTES:
            return text, False
        limit = MAX_OUTPUT_BYTES - len(TRUNCATION_MARKER)
        return f"{text[:limit]}{TRUNCATION_MARKER}", True


def resolve_probe_argv(spec: ProbeSpec, **slots: str) -> list[str]:
    """Resolve a probe spec into a validated argv list."""

    if not isinstance(spec, ProbeSpec):
        raise ProbeError("raw command string not allowed; use ProbeSpec")

    slot_names = _extract_slot_names(spec.argv_template)
    unexpected = sorted(set(slots) - slot_names)
    missing = sorted(slot_names - set(slots))
    if unexpected:
        raise ProbeError(
            "Unexpected probe slots for "
            f"{spec.template_id}: {', '.join(unexpected)}"
        )
    if missing:
        raise ProbeError(
            f"Missing probe slots for {spec.template_id}: {', '.join(missing)}"
        )

    resolved_slots: dict[str, str] = {}
    for slot_name in slot_names:
        validator = spec.slot_validators.get(slot_name)
        if validator is None:
            raise ProbeError(
                f"Probe slot {slot_name!r} has no validator in {spec.template_id}"
            )
        value = slots[slot_name]
        if not isinstance(value, str):
            value = str(value)
        validator(value)
        resolved_slots[slot_name] = value

    return [part.format_map(resolved_slots) for part in spec.argv_template]


def run_local_probe(runner, spec: ProbeSpec, **slots: str):
    """Dispatch a local probe through either spec-aware or legacy runners."""

    if getattr(runner, "supports_probe_specs", False):
        return runner.run_local(spec, **slots)
    return runner.run_local(resolve_probe_argv(spec, **slots))


def run_ssh_probe(runner, host: str, spec: ProbeSpec, **slots: str):
    """Dispatch an SSH probe through either spec-aware or legacy runners."""

    if getattr(runner, "supports_probe_specs", False):
        return runner.run_ssh(host, spec, **slots)
    return runner.run_ssh(host, shlex.join(resolve_probe_argv(spec, **slots)))


def validate_identifier(value: str) -> None:
    """Validate an identifier-style slot value."""

    _validate_slot_chars(value)
    if not _IDENTIFIER_PATTERN.fullmatch(value):
        raise ProbeError(f"Invalid identifier slot value: {value!r}")


def validate_env_var_name(value: str) -> None:
    """Validate an environment-variable name slot value."""

    _validate_slot_chars(value)
    if not _ENV_VAR_PATTERN.fullmatch(value):
        raise ProbeError(f"Invalid env var slot value: {value!r}")


def validate_path_slot(value: str) -> None:
    """Validate a path slot value."""

    _validate_slot_chars(value)
    if not _PATH_PATTERN.fullmatch(value):
        raise ProbeError(f"Invalid path slot value: {value!r}")


def _validate_slot_chars(value: str) -> None:
    if not value:
        raise ProbeError("Probe slot value cannot be empty.")
    if any(char in _FORBIDDEN_SLOT_CHARS for char in value):
        raise ProbeError(f"Forbidden probe slot characters in {value!r}")


def _extract_slot_names(argv_template: tuple[str, ...]) -> set[str]:
    formatter = string.Formatter()
    names: set[str] = set()
    for token in argv_template:
        for _, field_name, _, _ in formatter.parse(token):
            if field_name:
                names.add(field_name)
    return names
