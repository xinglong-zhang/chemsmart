#!/usr/bin/env python3
"""Measure ChemSmart agent structure without importing production modules."""

from __future__ import annotations

import argparse
import ast
import json
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class SymbolMetric:
    path: str
    qualified_name: str
    kind: str
    start_line: int
    end_line: int
    lines: int


def python_files(root: Path) -> list[Path]:
    return sorted(path for path in root.rglob("*.py") if "__pycache__" not in path.parts)


def module_name(path: Path, repo_root: Path) -> str:
    relative = path.relative_to(repo_root).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def source_lines(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def symbol_metrics(path: Path, repo_root: Path, tree: ast.AST) -> list[SymbolMetric]:
    relative = path.relative_to(repo_root).as_posix()
    metrics: list[SymbolMetric] = []

    def visit(body: Iterable[ast.stmt], prefix: tuple[str, ...] = ()) -> None:
        for node in body:
            if not isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            end = int(getattr(node, "end_lineno", node.lineno))
            kind = "class" if isinstance(node, ast.ClassDef) else "function"
            name = ".".join((*prefix, node.name))
            metrics.append(
                SymbolMetric(relative, name, kind, node.lineno, end, end - node.lineno + 1)
            )
            visit(node.body, (*prefix, node.name))

    visit(getattr(tree, "body", ()))
    return metrics


def _relative_import_name(
    source_module: str,
    *,
    source_is_package: bool,
    level: int,
    imported_module: str | None,
) -> str:
    package = source_module.split(".")
    if not source_is_package:
        package.pop()
    if level > 1:
        package = package[: -(level - 1)]
    if imported_module:
        package.extend(imported_module.split("."))
    return ".".join(package)


def imported_agent_modules(
    tree: ast.AST,
    source_module: str,
    *,
    source_is_package: bool,
) -> set[str]:
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names = (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported = node.module
            if node.level:
                imported = _relative_import_name(
                    source_module,
                    source_is_package=source_is_package,
                    level=node.level,
                    imported_module=imported,
                )
            names = (imported,) if imported else ()
        else:
            continue
        imports.update(name for name in names if name.startswith("chemsmart.agent"))
    return imports


def resolve_import(name: str, modules: set[str]) -> str | None:
    candidate = name
    while candidate.startswith("chemsmart.agent"):
        if candidate in modules:
            return candidate
        candidate = candidate.rpartition(".")[0]
    return None


def strongly_connected_components(graph: dict[str, set[str]]) -> list[list[str]]:
    index = 0
    indices: dict[str, int] = {}
    lowlinks: dict[str, int] = {}
    stack: list[str] = []
    on_stack: set[str] = set()
    components: list[list[str]] = []

    def connect(node: str) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)
        for target in graph.get(node, set()):
            if target not in indices:
                connect(target)
                lowlinks[node] = min(lowlinks[node], lowlinks[target])
            elif target in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[target])
        if lowlinks[node] != indices[node]:
            return
        component: list[str] = []
        while stack:
            member = stack.pop()
            on_stack.remove(member)
            component.append(member)
            if member == node:
                break
        if len(component) > 1 or node in graph.get(node, set()):
            components.append(sorted(component))

    for node in sorted(graph):
        if node not in indices:
            connect(node)
    return sorted(components)


def _run_json_command(command: list[str], cwd: Path) -> list[dict[str, object]] | None:
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            check=False,
            capture_output=True,
            text=True,
        )
        payload = json.loads(result.stdout or "[]")
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, list) else None


def _ruff_findings(
    repo_root: Path,
    agent_root: Path,
    *,
    select: str | None = None,
) -> list[dict[str, object]] | None:
    executable = shutil.which("ruff")
    if executable is None:
        return None
    command = [
        executable,
        "check",
        str(agent_root),
        "--output-format=json",
    ]
    if select:
        command.extend(("--select", select))
    findings = _run_json_command(command, repo_root)
    if findings is None:
        return None
    for finding in findings:
        filename = finding.get("filename")
        if not isinstance(filename, str):
            continue
        try:
            finding["filename"] = Path(filename).resolve().relative_to(
                repo_root
            ).as_posix()
        except ValueError:
            finding["filename"] = Path(filename).name
    return findings


def _pylint_duplicate_findings(
    repo_root: Path,
    agent_root: Path,
) -> list[dict[str, object]] | None:
    if shutil.which("pylint") is None:
        return None
    return _run_json_command(
        [
            sys.executable,
            "-m",
            "pylint",
            "--disable=all",
            "--enable=duplicate-code",
            "--min-similarity-lines=8",
            "--output-format=json",
            str(agent_root),
        ],
        repo_root,
    )


def _coverage_metrics(path: Path | None) -> dict[str, object] | None:
    if path is None:
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    totals = payload.get("totals")
    if not isinstance(totals, dict):
        return None
    return {
        key: totals[key]
        for key in (
            "covered_lines",
            "num_statements",
            "percent_covered",
            "covered_branches",
            "num_branches",
        )
        if key in totals
    }


def _git_value(repo_root: Path, *args: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    value = result.stdout.strip()
    return value or None


def build_report(
    repo_root: Path,
    agent_root: Path,
    *,
    coverage_json: Path | None = None,
) -> dict[str, object]:
    files = python_files(agent_root)
    parsed: dict[Path, ast.AST] = {}
    loc: dict[str, int] = {}
    symbols: list[SymbolMetric] = []
    for path in files:
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        parsed[path] = tree
        relative = path.relative_to(repo_root).as_posix()
        loc[relative] = source_lines(path)
        symbols.extend(symbol_metrics(path, repo_root, tree))

    module_by_path = {path: module_name(path, repo_root) for path in files}
    modules = set(module_by_path.values())
    graph: dict[str, set[str]] = {module: set() for module in modules}
    for path, tree in parsed.items():
        source = module_by_path[path]
        for imported in imported_agent_modules(
            tree,
            source,
            source_is_package=path.name == "__init__.py",
        ):
            target = resolve_import(imported, modules)
            if target is not None and target != source:
                graph[source].add(target)

    fan_in: dict[str, int] = defaultdict(int)
    for targets in graph.values():
        for target in targets:
            fan_in[target] += 1

    oversized_files = [
        {"path": path, "lines": lines}
        for path, lines in sorted(loc.items(), key=lambda item: (-item[1], item[0]))
        if lines > 800
    ]
    long_classes = [asdict(metric) for metric in symbols if metric.kind == "class" and metric.lines > 500]
    long_functions = [asdict(metric) for metric in symbols if metric.kind == "function" and metric.lines > 100]
    hotspots = sorted(
        (
            {
                "module": module,
                "fan_in": fan_in[module],
                "fan_out": len(graph[module]),
            }
            for module in modules
        ),
        key=lambda item: (-(item["fan_in"] + item["fan_out"]), item["module"]),
    )
    cycles = strongly_connected_components(graph)
    ruff_default = _ruff_findings(repo_root, agent_root)
    ruff_c901 = _ruff_findings(repo_root, agent_root, select="C901")
    pylint_duplicates = _pylint_duplicate_findings(repo_root, agent_root)
    return {
        "schema_version": 1,
        "git": {
            "branch": _git_value(repo_root, "branch", "--show-current"),
            "head": _git_value(repo_root, "rev-parse", "HEAD"),
        },
        "root": agent_root.relative_to(repo_root).as_posix(),
        "summary": {
            "python_files": len(files),
            "physical_lines": sum(loc.values()),
            "functions": sum(metric.kind == "function" for metric in symbols),
            "classes": sum(metric.kind == "class" for metric in symbols),
            "files_over_800_lines": len(oversized_files),
            "classes_over_500_lines": len(long_classes),
            "functions_over_100_lines": len(long_functions),
            "import_cycles": len(cycles),
            "ruff_default_violations": (
                len(ruff_default) if ruff_default is not None else None
            ),
            "ruff_c901_violations": (
                len(ruff_c901) if ruff_c901 is not None else None
            ),
            "pylint_duplicate_groups": (
                len(pylint_duplicates)
                if pylint_duplicates is not None
                else None
            ),
        },
        "files": dict(sorted(loc.items(), key=lambda item: (-item[1], item[0]))),
        "oversized_files": oversized_files,
        "long_classes": long_classes,
        "long_functions": long_functions,
        "dependency_hotspots": hotspots[:25],
        "import_cycles": cycles,
        "coverage": _coverage_metrics(coverage_json),
        "tooling": {
            "ruff_default": ruff_default,
            "ruff_c901": ruff_c901,
            "pylint_duplicates": pylint_duplicates,
        },
    }


def markdown(report: dict[str, object]) -> str:
    summary = report["summary"]
    assert isinstance(summary, dict)
    lines = [
        "# ChemSmart Agent Architecture Audit",
        "",
        (
            "This is the measured pre-refactor baseline. It is not evidence "
            "that the structural completion gates have been met."
        ),
        "",
    ]
    git = report.get("git")
    if isinstance(git, dict):
        lines.extend(
            [
                f"- Branch: `{git.get('branch') or 'unknown'}`",
                f"- HEAD: `{git.get('head') or 'unknown'}`",
                "",
            ]
        )
    lines.extend(
        [
        "| Metric | Value |",
        "|---|---:|",
        ]
    )
    labels = {
        "python_files": "Python files",
        "physical_lines": "Physical lines",
        "functions": "Functions and methods",
        "classes": "Classes",
        "files_over_800_lines": "Files over 800 lines",
        "classes_over_500_lines": "Classes over 500 lines",
        "functions_over_100_lines": "Functions over 100 lines",
        "import_cycles": "Import cycles",
        "ruff_default_violations": "Ruff default violations",
        "ruff_c901_violations": "Ruff C901 violations",
        "pylint_duplicate_groups": "Pylint duplicate groups",
    }
    lines.extend(
        f"| {labels[key]} | "
        f"{summary[key] if summary[key] is not None else 'not measured'} |"
        for key in labels
    )
    coverage = report.get("coverage")
    if isinstance(coverage, dict):
        percent = float(coverage.get("percent_covered", 0.0))
        lines.append(f"| Branch-aware test coverage | {percent:.2f}% |")
    lines.extend(["", "## Oversized Files", "", "| File | Lines |", "|---|---:|"])
    oversized = report.get("oversized_files") or []
    assert isinstance(oversized, list)
    lines.extend(f"| `{row['path']}` | {row['lines']} |" for row in oversized)
    lines.extend(
        [
            "",
            "## Long Classes",
            "",
            "| Class | File | Lines |",
            "|---|---|---:|",
        ]
    )
    long_classes = report.get("long_classes") or []
    assert isinstance(long_classes, list)
    lines.extend(
        f"| `{row['qualified_name']}` | `{row['path']}` | {row['lines']} |"
        for row in sorted(long_classes, key=lambda item: -item["lines"])
    )
    lines.extend(
        [
            "",
            "## Longest Functions",
            "",
            "| Function | File | Lines |",
            "|---|---|---:|",
        ]
    )
    long_functions = report.get("long_functions") or []
    assert isinstance(long_functions, list)
    lines.extend(
        f"| `{row['qualified_name']}` | `{row['path']}` | {row['lines']} |"
        for row in sorted(
            long_functions,
            key=lambda item: (-item["lines"], item["path"]),
        )[:20]
    )
    lines.extend(
        [
            "",
            "## Dependency Hotspots",
            "",
            "| Module | Fan-in | Fan-out |",
            "|---|---:|---:|",
        ]
    )
    hotspots = report.get("dependency_hotspots") or []
    assert isinstance(hotspots, list)
    lines.extend(
        f"| `{row['module']}` | {row['fan_in']} | {row['fan_out']} |"
        for row in hotspots[:15]
    )
    lines.extend(["", "## Import Cycles", ""])
    cycles = report.get("import_cycles") or []
    assert isinstance(cycles, list)
    if cycles:
        lines.extend(
            f"- `{' -> '.join(str(item) for item in cycle)}`"
            for cycle in cycles
        )
    else:
        lines.append("None.")
    lines.extend(["", "## Duplicate Groups", ""])
    tooling = report.get("tooling") or {}
    assert isinstance(tooling, dict)
    duplicates = tooling.get("pylint_duplicates") or []
    assert isinstance(duplicates, list)
    if duplicates:
        for index, item in enumerate(duplicates, start=1):
            message = str(item.get("message") or "")
            locations = [
                line.removeprefix("==")
                for line in message.splitlines()
                if line.startswith("==")
            ]
            lines.append(f"{index}. `{'` and `'.join(locations)}`")
    else:
        lines.append("None.")
    lines.extend(
        [
            "",
            "## Completion Targets",
            "",
            (
                "The final branch must reduce every production agent file to "
                "800 lines or fewer, every class to 500 lines or fewer, and "
                "every function to 100 lines or fewer. C901 findings, duplicate "
                "groups, and import cycles must all reach zero without "
                "increasing total agent production LOC or changing the frozen "
                "behavior contracts."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--json", type=Path)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--coverage-json", type=Path)
    args = parser.parse_args()
    repo_root = args.repo.resolve()
    report = build_report(
        repo_root,
        repo_root / "chemsmart/agent",
        coverage_json=args.coverage_json,
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown(report), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
