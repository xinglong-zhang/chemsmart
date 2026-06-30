from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Literal

from chemsmart.agent.harness.extractors import (
    extract_gaussian_route,
    extract_orca_route,
)

CommandSemanticVerdict = Literal["ok", "warn", "reject"]
CommandSemanticSeverity = Literal["warn", "reject"]

_COMPUTATIONAL_TOP_LEVEL = {"run", "sub"}
_SOFTWARE_COMMANDS = {"gaussian", "orca"}
_INPUT_SUFFIXES = (".com", ".gjf", ".inp")


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
            f"- {issue.rule_id}: {issue.message}"
            for issue in self.issues
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


def _absolutize_file_args(argv: list[str], base: Path) -> list[str]:
    """Rewrite relative input-file tokens to absolute paths against ``base``.

    Needed when the safe exec runs in an isolated temp cwd: a ``-f examples/x.xyz``
    style token would not resolve there. Only file-looking tokens (with a path
    separator or an extension) that actually exist are rewritten, so option values
    like a project name (``-p test``) or charge/multiplicity are never touched.
    """
    out: list[str] = []
    for tok in argv:
        if (
            tok
            and not tok.startswith("-")
            and not os.path.isabs(tok)
            and (os.sep in tok or os.path.splitext(tok)[1])
            and (base / tok).exists()
        ):
            out.append(str((base / tok).resolve()))
        else:
            out.append(tok)
    return out


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

    try:
        tokens = shlex.split(command)
    except ValueError as exc:
        return _single_issue_result(
            command,
            "cmd.semantic.shlex",
            f"command could not be tokenized: {exc}",
        )
    if not tokens or tokens[0] != "chemsmart":
        return _single_issue_result(
            command,
            "cmd.semantic.not_chemsmart",
            "command must start with 'chemsmart'",
        )

    parse_issue = _strict_parser_issue(command, tokens)
    if parse_issue is not None:
        return CommandSemanticResult(
            verdict="reject",
            command=command,
            checked_argv=tuple(tokens),
            issues=(parse_issue,),
        )

    top_index, top_level = _top_level_command(tokens)
    if top_level not in _COMPUTATIONAL_TOP_LEVEL:
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

    option_order_issue = _option_order_issue(tokens, top_index)
    if option_order_issue is not None:
        return CommandSemanticResult(
            verdict="reject",
            command=command,
            checked_argv=tuple(tokens),
            issues=(option_order_issue,),
        )

    safety_issue = _safety_issue(tokens, top_index, top_level)
    if safety_issue is not None:
        return CommandSemanticResult(
            verdict="reject",
            command=command,
            checked_argv=tuple(tokens),
            issues=(safety_issue,),
        )

    safe_argv = _safe_execution_argv(tokens, top_index, top_level)
    if cwd is None:
        # Run in an isolated temp dir so a stale <label>.out in the real cwd cannot
        # trigger chemsmart's skip-completed path (which writes no input and yields a
        # false generated_input_missing reject). Absolutize input files so relative
        # -f/-e tokens still resolve from the temp cwd.
        workdir = Path(tempfile.mkdtemp(prefix="chemsmart-gate-"))
        safe_argv = _absolutize_file_args(safe_argv, Path(os.getcwd()))
        cleanup_dir: Path | None = workdir
    else:
        workdir = Path(cwd)
        cleanup_dir = None
    try:
        before = _input_snapshot(workdir)
        try:
            completed = subprocess.run(
                safe_argv,
                cwd=str(workdir),
                env=_subprocess_env(),
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
                        message=f"safe runtime validation timed out after {timeout_s:g}s",
                        evidence={"timeout_s": timeout_s},
                    ),
                ),
                stdout_tail=_tail(exc.stdout),
                stderr_tail=_tail(exc.stderr),
            )

        stdout_tail = _tail(completed.stdout)
        stderr_tail = _tail(completed.stderr)
        generated_inputs = _generated_inputs(workdir, before)
        issues: list[CommandSemanticIssue] = []
        if completed.returncode != 0:
            issues.append(
                CommandSemanticIssue(
                    rule_id="cmd.semantic.safe_execution_failed",
                    severity="reject",
                    message=(
                        "safe chemsmart runtime validation failed with exit code "
                        f"{completed.returncode}"
                    ),
                    evidence={
                        "returncode": completed.returncode,
                        "argv": safe_argv,
                        "stdout_tail": stdout_tail,
                        "stderr_tail": stderr_tail,
                    },
                    missing_info=tuple(_missing_info_from_output(stderr_tail, stdout_tail)),
                )
            )
        elif (
            top_level == "run"
            and _is_software_command(tokens, top_index)
            and not generated_inputs
        ):
            issues.append(
                CommandSemanticIssue(
                    rule_id="cmd.semantic.generated_input_missing",
                    severity="reject",
                    message=(
                        "safe runtime validation succeeded but no Gaussian/ORCA "
                        "input file was generated"
                    ),
                    evidence={"cwd": str(workdir), "argv": safe_argv},
                )
            )
        elif _is_software_command(tokens, top_index) and not generated_inputs:
            issues.append(
                CommandSemanticIssue(
                    rule_id="cmd.semantic.submit_generated_input_not_observed",
                    severity="warn",
                    message=(
                        "safe submit validation succeeded but no Gaussian/ORCA "
                        "input file was observed in the working directory"
                    ),
                    evidence={"cwd": str(workdir), "argv": safe_argv},
                )
            )

        for generated in generated_inputs:
            if not generated.get("route"):
                issues.append(
                    CommandSemanticIssue(
                        rule_id="cmd.semantic.generated_route_missing",
                        severity="reject",
                        message="generated computational input is missing a route line",
                        evidence={"path": generated.get("path")},
                    )
                )

        verdict: CommandSemanticVerdict = "ok"
        if any(issue.severity == "reject" for issue in issues):
            verdict = "reject"
        elif issues:
            verdict = "warn"
        return CommandSemanticResult(
            verdict=verdict,
            command=command,
            checked_argv=tuple(safe_argv),
            issues=tuple(issues),
            generated_inputs=tuple(generated_inputs),
            stdout_tail=stdout_tail,
            stderr_tail=stderr_tail,
        )
    finally:
        if cleanup_dir is not None:
            shutil.rmtree(cleanup_dir, ignore_errors=True)


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
            missing_info=("program-level input filename before job subcommand",),
        )
    return None


def _software_index(tokens: list[str], top_index: int) -> int | None:
    for index, token in enumerate(tokens[top_index + 1 :], start=top_index + 1):
        if token in _SOFTWARE_COMMANDS:
            return index
    return None


def _first_job_token_index(tokens: list[str], start: int) -> int | None:
    index = start
    while index < len(tokens):
        token = tokens[index]
        if token.startswith("-"):
            index += 2 if _option_consumes_value(token, tokens[index + 1 :]) else 1
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


def _safe_execution_argv(
    tokens: list[str],
    top_index: int,
    top_level: str,
) -> list[str]:
    argv = _with_no_verbose(tokens)
    # Recompute because _with_no_verbose may insert a token before top_index.
    top_index, _top_level = _top_level_command(argv)
    insert_at = top_index + 1
    additions: list[str] = []
    if top_level == "run":
        if "--fake" not in argv[insert_at:]:
            additions.append("--fake")
        if (
            "--scratch" not in argv[insert_at:]
            and "--no-scratch" not in argv[insert_at:]
        ):
            additions.append("--no-scratch")
    elif top_level == "sub":
        if "--test" not in argv[insert_at:]:
            additions.append("--test")
        if "--fake" not in argv[insert_at:]:
            additions.append("--fake")
    cli_args = argv[1:insert_at] + additions + argv[insert_at:]
    return [sys.executable, "-m", "chemsmart.cli.main", *cli_args]


def _subprocess_env() -> dict[str, str]:
    env = dict(os.environ)
    source_root = str(Path(__file__).resolve().parents[3])
    existing = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        source_root if not existing else f"{source_root}{os.pathsep}{existing}"
    )
    return env


def _with_no_verbose(tokens: list[str]) -> list[str]:
    if "--verbose" in tokens[:3] or "--no-verbose" in tokens[:3]:
        return list(tokens)
    return [tokens[0], "--no-verbose", *tokens[1:]]


def _input_snapshot(workdir: Path) -> dict[Path, int]:
    snapshot: dict[Path, int] = {}
    if not workdir.exists():
        return snapshot
    for suffix in _INPUT_SUFFIXES:
        for path in workdir.glob(f"*{suffix}"):
            try:
                snapshot[path.resolve()] = path.stat().st_mtime_ns
            except FileNotFoundError:
                continue
    return snapshot


def _generated_inputs(
    workdir: Path,
    before: dict[Path, int],
) -> list[dict[str, Any]]:
    generated: list[dict[str, Any]] = []
    if not workdir.exists():
        return generated
    for suffix in _INPUT_SUFFIXES:
        for path in sorted(workdir.glob(f"*{suffix}")):
            try:
                resolved = path.resolve()
                mtime = path.stat().st_mtime_ns
            except FileNotFoundError:
                continue
            if before.get(resolved) == mtime:
                continue
            content = path.read_text(encoding="utf-8", errors="replace")
            route = (
                extract_gaussian_route(content)
                if suffix in {".com", ".gjf"}
                else extract_orca_route(content)
            )
            generated.append(
                {
                    "path": str(path),
                    "route": route,
                    "content_tail": _tail(content),
                }
            )
    return generated


def _is_software_command(tokens: list[str], top_index: int) -> bool:
    return any(token in _SOFTWARE_COMMANDS for token in tokens[top_index + 1 :])


def _missing_info_from_output(*chunks: str) -> list[str]:
    text = "\n".join(chunks)
    missing: list[str] = []
    if "No server implemented" in text:
        missing.append("valid chemsmart server configuration")
    if "Currently available servers:" in text:
        tail = text.split("Currently available servers:", 1)[1].splitlines()[0]
        missing.append(f"available servers: {tail.strip()}")
    if "No project settings implemented" in text:
        missing.append("valid chemsmart project configuration")
    if "Currently available projects:" in text:
        tail = text.split("Currently available projects:", 1)[1].splitlines()[0]
        missing.append(f"available projects: {tail.strip()}")
    if "No such file or directory" in text or "FileNotFoundError" in text:
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
