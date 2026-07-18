from __future__ import annotations

import os
import shutil
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import yaml

from chemsmart.agent.harness.extractors import (
    extract_cartesian_state,
    extract_gaussian_route,
    extract_orca_route,
)

_INPUT_SUFFIXES = (".com", ".gjf", ".inp")


def absolutize_file_args(argv: list[str], base: Path) -> list[str]:
    """Resolve existing file arguments before entering an isolated cwd."""
    resolved: list[str] = []
    for token in argv:
        if (
            token
            and not token.startswith("-")
            and not os.path.isabs(token)
            and (os.sep in token or os.path.splitext(token)[1])
            and (base / token).exists()
        ):
            resolved.append(str((base / token).resolve()))
        else:
            resolved.append(token)
    return resolved


def safe_execution_argv(
    tokens: list[str],
    top_index: int,
    top_level: str,
) -> list[str]:
    """Build a non-destructive CLI invocation for runtime validation."""
    argv = _with_no_verbose(tokens)
    # ``--no-verbose`` may have shifted the top-level command by one token.
    top_index = _top_level_index(argv, top_level)
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


def prepare_safe_runtime_environment(
    *,
    base_cwd: Path,
    workdir: Path,
    top_level: str,
) -> dict[str, str]:
    """Mirror workspace state and isolate local fake-run configuration.

    ``chemsmart run`` falls back to ``~/.chemsmart/server/local.yaml``. A
    semantic gate must not pass or fail according to an unrelated user HOME,
    so fake local runs receive a minimal run-local profile. ``sub`` keeps the
    caller HOME because its named server fixture is part of submission intent.
    """
    _mirror_workspace_config(base_cwd, workdir)
    env = _subprocess_env()
    if top_level != "run":
        return env

    gate_home = workdir / ".gate-home"
    server_dir = gate_home / ".chemsmart" / "server"
    server_dir.mkdir(parents=True, exist_ok=True)
    fake_bin = workdir / ".fake-executables"
    fake_bin.mkdir(exist_ok=True)
    local_server = {
        "SERVER": {
            "SCHEDULER": None,
            "QUEUE_NAME": None,
            "NUM_HOURS": None,
            "MEM_GB": 40,
            "NUM_CORES": 12,
            "NUM_GPUS": 0,
            "NUM_THREADS": 12,
            "SUBMIT_COMMAND": None,
            "SCRATCH_DIR": None,
            "USE_HOSTS": False,
        },
        "GAUSSIAN": {
            "EXEFOLDER": str(fake_bin),
            "LOCAL_RUN": True,
            "SCRATCH": False,
        },
        "ORCA": {
            "EXEFOLDER": str(fake_bin),
            "LOCAL_RUN": True,
            "SCRATCH": False,
        },
    }
    (server_dir / "local.yaml").write_text(
        yaml.safe_dump(local_server, sort_keys=False),
        encoding="utf-8",
    )
    env["HOME"] = str(gate_home)
    # Windows resolves ``~`` through USERPROFILE, not HOME; keep the gate
    # home authoritative on every platform.
    env["USERPROFILE"] = str(gate_home)
    return env


def input_snapshot(workdir: Path) -> dict[Path, int]:
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


def generated_inputs(
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
            software = "gaussian" if suffix in {".com", ".gjf"} else "orca"
            route = (
                extract_gaussian_route(content)
                if software == "gaussian"
                else extract_orca_route(content)
            )
            state = extract_cartesian_state(content, software=software)
            state_evidence: dict[str, Any] = {}
            if state:
                state_evidence = {
                    "charge": state["charge"],
                    "multiplicity": state["multiplicity"],
                    "element_counts": dict(
                        sorted(Counter(state["element_symbols"]).items())
                    ),
                }
                for key in (
                    "charge_multiplicity_pairs",
                    "atom_layers",
                    "layer_atoms",
                ):
                    if key in state:
                        state_evidence[key] = state[key]
            generated.append(
                {
                    "path": str(path),
                    "route": route,
                    "content_tail": input_excerpt(content),
                    **state_evidence,
                }
            )
    return generated


def _top_level_index(argv: list[str], top_level: str) -> int:
    try:
        return argv.index(top_level, 1)
    except ValueError as exc:  # pragma: no cover - caller already validated it
        raise ValueError(f"missing top-level command: {top_level}") from exc


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


def _mirror_workspace_config(base_cwd: Path, workdir: Path) -> None:
    source = base_cwd / ".chemsmart"
    if source.is_dir():
        shutil.copytree(source, workdir / ".chemsmart", dirs_exist_ok=True)


def input_excerpt(value: Any, limit: int = 4000) -> str:
    if value is None:
        return ""
    text = (
        value.decode("utf-8", errors="replace")
        if isinstance(value, bytes)
        else str(value)
    )
    if len(text) <= limit:
        return text
    half = limit // 2
    return f"{text[:half]}\n...<content omitted>...\n{text[-half:]}"


__all__ = [
    "absolutize_file_args",
    "generated_inputs",
    "input_excerpt",
    "input_snapshot",
    "prepare_safe_runtime_environment",
    "safe_execution_argv",
]
