from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

from chemsmart.agent.harness.command_contracts import check_command_contracts
from chemsmart.agent.harness.command_rules.db_selectors import (
    db_selector_issues,
)
from chemsmart.agent.harness.command_rules.tokens import flag_option_present
from chemsmart.agent.harness.failure_taxonomy import classify_runtime_failure
from chemsmart.agent.harness.generated_invariants import (
    check_generated_input_invariants,
)
from chemsmart.agent.harness.safe_runtime import (
    absolutize_file_args,
    generated_inputs as collect_generated_inputs,
    input_snapshot,
    prepare_safe_runtime_environment,
    safe_execution_argv,
)

CommandSemanticVerdict = Literal["ok", "warn", "reject"]
CommandSemanticSeverity = Literal["warn", "reject"]

_COMPUTATIONAL_TOP_LEVEL = {"run", "sub"}
_SOFTWARE_COMMANDS = {"gaussian", "orca"}


@dataclass(frozen=True)
class CommandSemanticIssue:
    rule_id: str
    severity: CommandSemanticSeverity
    message: str
    evidence: dict[str, Any] = field(default_factory=dict)
    missing_info: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class CommandSemanticResult:
    verdict: CommandSemanticVerdict
    command: str
    checked_argv: tuple[str, ...] = ()
    issues: tuple[CommandSemanticIssue, ...] = ()
    generated_inputs: tuple[dict[str, Any], ...] = ()
    stdout_tail: str = ""
    stderr_tail: str = ""

    @property
    def failed_rule_ids(self) -> list[str]:
        seen: set[str] = set()
        failed: list[str] = []
        for issue in self.issues:
            if issue.rule_id not in seen:
                seen.add(issue.rule_id)
                failed.append(issue.rule_id)
        return failed

    @property
    def missing_info(self) -> list[str]:
        seen: set[str] = set()
        missing: list[str] = []
        for issue in self.issues:
            for item in issue.missing_info:
                if item not in seen:
                    seen.add(item)
                    missing.append(item)
        return missing

    @property
    def notice(self) -> str:
        if self.verdict == "ok":
            return "runtime semantic validation passed"
        if self.verdict == "warn":
            return "runtime semantic validation produced warnings"
        return "runtime semantic validation rejected the command"

    def correction_prompt(self, original_request: str) -> str:
        issue_lines = "\n".join(
            f"- {issue.rule_id}: {issue.message}" for issue in self.issues
        )
        missing = ", ".join(self.missing_info) or "none"
        generated = "\n".join(
            f"- {item.get('path')}: {item.get('route') or 'route not found'}"
            for item in self.generated_inputs
        )
        if not generated:
            generated = "- none"
        return (
            f"Original request: {original_request}\n"
            f"Your command passed shallow schema validation but failed "
            f"runtime semantic validation.\n"
            f"Command: {self.command}\n"
            f"Safe validation argv: {shlex.join(self.checked_argv)}\n"
            f"Verdict: {self.verdict}\n"
            f"Missing/insufficient runtime information: {missing}\n"
            f"Issues:\n{issue_lines or '- none'}\n"
            f"Generated input evidence:\n{generated}\n"
            f"stdout tail:\n{self.stdout_tail or '(empty)'}\n"
            f"stderr tail:\n{self.stderr_tail or '(empty)'}\n"
            "Return ONLY corrected JSON matching the synthesis schema. "
            "If the command needs runtime information such as a server, "
            "include the appropriate legal chemsmart flags in the command. "
            "Do not invent options."
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "verdict": self.verdict,
            "command": self.command,
            "checked_argv": list(self.checked_argv),
            "issues": [issue.to_dict() for issue in self.issues],
            "failed_rule_ids": self.failed_rule_ids,
            "missing_info": self.missing_info,
            "generated_inputs": list(self.generated_inputs),
            "stdout_tail": self.stdout_tail,
            "stderr_tail": self.stderr_tail,
            "notice": self.notice,
        }


def evaluate_command_semantics(
    command: str,
    *,
    cwd: str | os.PathLike[str] | None = None,
    timeout_s: float = 30.0,
) -> CommandSemanticResult:
    """Validate a synthesized chemsmart command through a safe execution path.

    This gate is intentionally stricter than click/schema validation. For
    computational commands it executes only a safe variant: ``run`` receives
    ``--fake --no-scratch`` and ``sub`` receives ``--test --fake``. The command
    still has to resolve runtime configuration, parse real files, build the job,
    and write Gaussian/ORCA input evidence.
    """

    tokenized = _tokenize_chemsmart_command(command)
    if isinstance(tokenized, CommandSemanticResult):
        return tokenized
    tokens = tokenized
    base_cwd = Path(cwd or os.getcwd()).resolve()
    top_index, top_level = _top_level_command(tokens)
    preflight = _preflight_result(
        command,
        tokens,
        top_index=top_index,
        top_level=top_level,
        cwd=base_cwd,
    )
    if preflight is not None:
        return preflight
    return _execute_safe_semantic_runtime(
        command,
        tokens,
        top_index=top_index,
        top_level=top_level,
        base_cwd=base_cwd,
        timeout_s=timeout_s,
    )


def _tokenize_chemsmart_command(
    command: str,
) -> list[str] | CommandSemanticResult:
    try:
        tokens = shlex.split(command)
    except ValueError as exc:
        return _single_issue_result(
            command,
            "cmd.semantic.shlex",
            f"command could not be tokenized: {exc}",
        )
    if tokens and tokens[0] == "chemsmart":
        return tokens
    return _single_issue_result(
        command,
        "cmd.semantic.not_chemsmart",
        "command must start with 'chemsmart'",
    )


def _preflight_result(
    command: str,
    tokens: list[str],
    *,
    top_index: int,
    top_level: str | None,
    cwd: Path,
) -> CommandSemanticResult | None:
    if top_level in _COMPUTATIONAL_TOP_LEVEL:
        contract_issues = _command_contract_issues(tokens, top_index, cwd=cwd)
        rejected = _rejected_result(command, tokens, contract_issues)
        if rejected is not None:
            return rejected
    parse_issue = _strict_parser_issue(command, tokens)
    if parse_issue is not None:
        return _rejected_result(command, tokens, (parse_issue,))
    if top_level not in _COMPUTATIONAL_TOP_LEVEL:
        return _non_computational_result(command, tokens, top_level)
    for issue in (
        _option_order_issue(tokens, top_index),
        _safety_issue(tokens, top_index, top_level),
    ):
        if issue is not None:
            return _rejected_result(command, tokens, (issue,))
    return _rejected_result(
        command,
        tokens,
        _preflight_semantic_issues(tokens, top_index),
    )


def _rejected_result(
    command: str,
    argv: list[str],
    issues: tuple[CommandSemanticIssue, ...],
) -> CommandSemanticResult | None:
    if not any(issue.severity == "reject" for issue in issues):
        return None
    return CommandSemanticResult(
        verdict="reject",
        command=command,
        checked_argv=tuple(argv),
        issues=issues,
    )


def _non_computational_result(
    command: str,
    tokens: list[str],
    top_level: str | None,
) -> CommandSemanticResult:
    return CommandSemanticResult(
        verdict="warn",
        command=command,
        checked_argv=tuple(tokens),
        issues=(
            CommandSemanticIssue(
                rule_id="cmd.semantic.not_computational",
                severity="warn",
                message=(
                    "runtime execution gate only hard-checks chemsmart "
                    "run/sub computational commands"
                ),
                evidence={"top_level": top_level},
            ),
        ),
    )


def _execute_safe_semantic_runtime(
    command: str,
    tokens: list[str],
    *,
    top_index: int,
    top_level: str,
    base_cwd: Path,
    timeout_s: float,
) -> CommandSemanticResult:
    safe_argv = safe_execution_argv(tokens, top_index, top_level)
    # Run in an isolated temp dir so a stale <label>.out in the real cwd cannot
    # trigger chemsmart's skip-completed path (which writes no input and yields a
    # false generated_input_missing reject). The caller's cwd is still the
    # workspace of record: relative input files are resolved from it, and its
    # workspace-local .chemsmart settings are mirrored into the temp dir so
    # project/server YAML lookup sees the same configuration the user sees.
    workdir = Path(tempfile.mkdtemp(prefix="chemsmart-gate-"))
    safe_argv = absolutize_file_args(safe_argv, base_cwd)
    runtime_env = prepare_safe_runtime_environment(
        base_cwd=base_cwd,
        workdir=workdir,
        top_level=top_level,
    )
    try:
        before = input_snapshot(workdir)
        completed = _run_safe_command(
            command,
            safe_argv,
            workdir=workdir,
            runtime_env=runtime_env,
            timeout_s=timeout_s,
        )
        if isinstance(completed, CommandSemanticResult):
            return completed
        stdout_tail = _tail(completed.stdout)
        stderr_tail = _tail(completed.stderr)
        generated_inputs = collect_generated_inputs(workdir, before)
        issues = _runtime_semantic_issues(
            tokens=tokens,
            top_index=top_index,
            top_level=top_level,
            safe_argv=safe_argv,
            workdir=workdir,
            returncode=completed.returncode,
            stdout_tail=stdout_tail,
            stderr_tail=stderr_tail,
            generated_inputs=generated_inputs,
        )
        issues.extend(_generated_route_issues(generated_inputs))
        issues.extend(
            _generated_invariant_issues(command, generated_inputs, base_cwd)
        )
        return CommandSemanticResult(
            verdict=_semantic_verdict(issues),
            command=command,
            checked_argv=tuple(safe_argv),
            issues=tuple(issues),
            generated_inputs=tuple(generated_inputs),
            stdout_tail=stdout_tail,
            stderr_tail=stderr_tail,
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _run_safe_command(
    command: str,
    safe_argv: list[str],
    *,
    workdir: Path,
    runtime_env: dict[str, str],
    timeout_s: float,
) -> subprocess.CompletedProcess[str] | CommandSemanticResult:
    try:
        return subprocess.run(
            safe_argv,
            cwd=str(workdir),
            env=runtime_env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return CommandSemanticResult(
            verdict="reject",
            command=command,
            checked_argv=tuple(safe_argv),
            issues=(
                CommandSemanticIssue(
                    rule_id="cmd.semantic.timeout",
                    severity="reject",
                    message=(
                        "safe runtime validation timed out after "
                        f"{timeout_s:g}s"
                    ),
                    evidence={"timeout_s": timeout_s},
                ),
            ),
            stdout_tail=_tail(exc.stdout),
            stderr_tail=_tail(exc.stderr),
        )


def _runtime_semantic_issues(
    *,
    tokens: list[str],
    top_index: int,
    top_level: str,
    safe_argv: list[str],
    workdir: Path,
    returncode: int,
    stdout_tail: str,
    stderr_tail: str,
    generated_inputs: list[dict[str, Any]],
) -> list[CommandSemanticIssue]:
    if returncode != 0:
        failure = classify_runtime_failure(
            stdout=stdout_tail,
            stderr=stderr_tail,
            returncode=returncode,
        )
        return [
            CommandSemanticIssue(
                rule_id=failure.rule_id,
                severity="reject",
                message=failure.message,
                evidence={
                    "category": failure.category,
                    "returncode": returncode,
                    "argv": safe_argv,
                    "stdout_tail": stdout_tail,
                    "stderr_tail": stderr_tail,
                },
                missing_info=(
                    failure.missing_info
                    or tuple(
                        _missing_info_from_output(stderr_tail, stdout_tail)
                    )
                ),
            )
        ]
    if not _is_software_command(tokens, top_index) or generated_inputs:
        return []
    if top_level == "run":
        return [
            CommandSemanticIssue(
                rule_id="cmd.semantic.generated_input_missing",
                severity="reject",
                message=(
                    "safe runtime validation succeeded but no Gaussian/ORCA "
                    "input file was generated"
                ),
                evidence={"cwd": str(workdir), "argv": safe_argv},
            )
        ]
    return [
        CommandSemanticIssue(
            rule_id="cmd.semantic.submit_generated_input_not_observed",
            severity="warn",
            message=(
                "safe submit validation succeeded but no Gaussian/ORCA "
                "input file was observed in the working directory"
            ),
            evidence={"cwd": str(workdir), "argv": safe_argv},
        )
    ]


def _generated_route_issues(
    generated_inputs: list[dict[str, Any]],
) -> list[CommandSemanticIssue]:
    return [
        CommandSemanticIssue(
            rule_id="cmd.semantic.generated_route_missing",
            severity="reject",
            message="generated computational input is missing a route line",
            evidence={"path": generated.get("path")},
        )
        for generated in generated_inputs
        if not generated.get("route")
    ]


def _generated_invariant_issues(
    command: str,
    generated_inputs: list[dict[str, Any]],
    base_cwd: Path,
) -> list[CommandSemanticIssue]:
    return [
        CommandSemanticIssue(
            rule_id=issue.rule_id,
            severity=issue.severity,
            message=issue.message,
            evidence=issue.evidence,
        )
        for issue in check_generated_input_invariants(
            command,
            generated_inputs,
            cwd=str(base_cwd),
        )
    ]


def _semantic_verdict(
    issues: list[CommandSemanticIssue],
) -> CommandSemanticVerdict:
    if any(issue.severity == "reject" for issue in issues):
        return "reject"
    return "warn" if issues else "ok"


def _strict_parser_issue(
    command: str,
    tokens: list[str],
) -> CommandSemanticIssue | None:
    try:
        from chemsmart.cli.main import entry_point

        entry_point.make_context(
            "chemsmart",
            tokens[1:],
            resilient_parsing=False,
        )
    except Exception as exc:
        return CommandSemanticIssue(
            rule_id="cmd.semantic.strict_parser",
            severity="reject",
            message=f"real click parser rejected command: {exc}",
            evidence={"command": command},
        )
    return None


def _top_level_command(tokens: list[str]) -> tuple[int, str | None]:
    index = 1
    while index < len(tokens):
        token = tokens[index]
        if token in {"--verbose", "--no-verbose", "--version", "-h", "--help"}:
            index += 1
            continue
        return index, token
    return index, None


def _safety_issue(
    tokens: list[str],
    top_index: int,
    top_level: str | None,
) -> CommandSemanticIssue | None:
    scoped = tokens[top_index + 1 :]
    if "--no-fake" in scoped:
        return CommandSemanticIssue(
            rule_id="cmd.semantic.unsafe_no_fake",
            severity="reject",
            message=(
                "semantic validation refuses to execute a command that "
                "explicitly disables fake runners"
            ),
            evidence={"top_level": top_level},
        )
    return None


def _option_order_issue(
    tokens: list[str],
    top_index: int,
) -> CommandSemanticIssue | None:
    software_index = _software_index(tokens, top_index)
    if software_index is None or software_index + 1 >= len(tokens):
        return None
    job_index = _first_job_token_index(tokens, software_index + 1)
    if job_index is None:
        return None
    before_job = tokens[software_index + 1 : job_index]
    after_job = tokens[job_index + 1 :]
    has_filename_before_job = any(
        token in {"-f", "--filename"} for token in before_job
    )
    has_filename_after_job = any(
        token in {"-f", "--filename"} for token in after_job
    )
    if not has_filename_before_job and has_filename_after_job:
        return CommandSemanticIssue(
            rule_id="cmd.semantic.option_order",
            severity="reject",
            message=(
                "input filename option appears after the Gaussian/ORCA job "
                "subcommand; chemsmart expects program-level options such as "
                "-f/--filename before the job subcommand"
            ),
            evidence={
                "program": tokens[software_index],
                "job_subcommand": tokens[job_index],
                "expected_shape": (
                    "chemsmart run <run-opts> "
                    f"{tokens[software_index]} -f molecule.xyz "
                    f"<program-opts> {tokens[job_index]}"
                ),
            },
            missing_info=(
                "program-level input filename before job subcommand",
            ),
        )
    return None


def _preflight_semantic_issues(
    tokens: list[str],
    top_index: int,
) -> tuple[CommandSemanticIssue, ...]:
    software_index = _software_index(tokens, top_index)
    if software_index is None:
        return ()
    program = tokens[software_index]
    if program not in _SOFTWARE_COMMANDS:
        return ()
    job_index = _first_job_token_index(tokens, software_index + 1)
    program_tokens = tokens[software_index + 1 : job_index or len(tokens)]

    issues: list[CommandSemanticIssue] = []
    if program == "orca" and flag_option_present(program_tokens, ("-a",)):
        issues.append(
            CommandSemanticIssue(
                rule_id="cmd.semantic.orca_aux_basis_short_flag",
                severity="reject",
                message=(
                    "ORCA auxiliary basis uses -B/--aux-basis in chemsmart "
                    "2.0.1; -a is append-label and is ambiguous for agent "
                    "generated commands"
                ),
                evidence={"program": program, "flag": "-a"},
            )
        )

    issues.extend(
        CommandSemanticIssue(
            rule_id=issue.rule_id,
            severity=issue.severity,
            message=issue.message,
            evidence=issue.evidence,
            missing_info=issue.missing_info,
        )
        for issue in db_selector_issues(program, program_tokens)
    )
    return tuple(issues)


def _command_contract_issues(
    tokens: list[str],
    top_index: int,
    *,
    cwd: Path,
) -> tuple[CommandSemanticIssue, ...]:
    software_index = _software_index(tokens, top_index)
    if software_index is None:
        return ()
    program = tokens[software_index]
    if program not in _SOFTWARE_COMMANDS:
        return ()
    job_index = _first_job_token_index(tokens, software_index + 1)
    if job_index is None:
        return (
            CommandSemanticIssue(
                rule_id="cmd.contract.job_subcommand_required",
                severity="reject",
                message=(
                    f"{program} requires an explicit ChemSmart job "
                    "subcommand such as opt, sp, ts, scan, or modred"
                ),
                evidence={
                    "program": program,
                    "command_tokens": tokens,
                },
                missing_info=(
                    f"explicit {program} computational job subcommand",
                ),
            ),
        )
    contract_issues = check_command_contracts(
        program=program,
        job=tokens[job_index],
        program_tokens=tokens[software_index + 1 : job_index],
        job_tokens=tokens[job_index:],
        cwd=cwd,
    )
    return tuple(
        CommandSemanticIssue(
            rule_id=issue.rule_id,
            severity=issue.severity,
            message=issue.message,
            evidence=issue.evidence,
            missing_info=issue.missing_info,
        )
        for issue in contract_issues
    )


def _software_index(tokens: list[str], top_index: int) -> int | None:
    for index, token in enumerate(
        tokens[top_index + 1 :], start=top_index + 1
    ):
        if token in _SOFTWARE_COMMANDS:
            return index
    return None


def _first_job_token_index(tokens: list[str], start: int) -> int | None:
    index = start
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("-"):
            index += (
                2 if _option_consumes_value(token, tokens[index + 1 :]) else 1
            )
            continue
        return index
    return None


def _option_consumes_value(token: str, following: list[str]) -> bool:
    if "=" in token:
        return False
    if not following:
        return False
    # Boolean flag handling is intentionally conservative. Consuming one value
    # here only helps locate the first positional job token for order checks.
    return not following[0].startswith("-")


def _is_software_command(tokens: list[str], top_index: int) -> bool:
    return any(
        token in _SOFTWARE_COMMANDS for token in tokens[top_index + 1 :]
    )


def _missing_info_from_output(*chunks: str) -> list[str]:
    text = "\n".join(chunks)
    missing: list[str] = []
    has_project_error = "No project settings implemented" in text
    has_server_error = "No server implemented" in text
    if "No server implemented" in text:
        missing.append("valid chemsmart server configuration")
    if "Currently available servers:" in text:
        tail = text.split("Currently available servers:", 1)[1].splitlines()[0]
        missing.append(f"available servers: {tail.strip()}")
    if has_project_error:
        missing.append("valid chemsmart project configuration")
    if "Currently available projects:" in text:
        tail = text.split("Currently available projects:", 1)[1].splitlines()[
            0
        ]
        missing.append(f"available projects: {tail.strip()}")
    if "No such file or directory" in text or (
        "FileNotFoundError" in text
        and not has_project_error
        and not has_server_error
    ):
        missing.append("existing local input file path")
    return missing


def _tail(value: Any, limit: int = 2000) -> str:
    if value is None:
        return ""
    text = (
        value.decode("utf-8", errors="replace")
        if isinstance(value, bytes)
        else str(value)
    )
    return text[-limit:]


def _single_issue_result(
    command: str,
    rule_id: str,
    message: str,
) -> CommandSemanticResult:
    return CommandSemanticResult(
        verdict="reject",
        command=command,
        issues=(
            CommandSemanticIssue(
                rule_id=rule_id,
                severity="reject",
                message=message,
            ),
        ),
    )


__all__ = [
    "CommandSemanticIssue",
    "CommandSemanticResult",
    "evaluate_command_semantics",
]
