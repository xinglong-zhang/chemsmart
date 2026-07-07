"""Live evaluation of ``chemsmart agent ask`` command synthesis quality.

Runs ten curated natural-language prompts through ``SynthesisSession`` with
the configured provider (``~/.chemsmart/agent/agent.yaml`` preferred; legacy
``api.env`` fallback), then scores each synthesized command against
expected substring tokens.

PASS criterion per case:
  * ``result["status"] == "ready"``
  * ``SynthesisSession.validate_command(result["command"])`` succeeds
  * Every token in ``expected_substrings`` appears verbatim in the command

This script is intentionally standalone (not on the package import path).
Run from the chemsmart repo root:

    python scripts/eval_synthesis.py

The script writes a per-case report to ``/tmp/eval_synthesis_run.txt`` and
prints a Rich summary table to stdout. API keys are never echoed.
"""

from __future__ import annotations

import re
import sys
import traceback
from pathlib import Path

from rich.console import Console
from rich.table import Table

# Cases curated to span the chemsmart subcommand surface and to mix
# Korean/English natural phrasing. Substrings are the non-negotiable tokens.
CASES: list[tuple[str, str, list[str]]] = [
    (
        "local single-point",
        "single point of h2o on local with b3lyp/6-31g*",
        ["chemsmart", "gaussian", "-b", "6-31g*"],
    ),
    (
        "sub opt project+basis+cores",
        "opt of oxetane with b3lyp/def2-svp on chemnode1, 8 cores",
        ["chemsmart", "sub", "gaussian", "opt", "-p", "-b", "def2-svp"],
    ),
    (
        "transition state",
        "transition state search for cope rearrangement, def2-tzvp basis, "
        "24 cores on chemnode1",
        ["chemsmart", "sub", "gaussian", "ts", "-b", "def2-tzvp"],
    ),
    (
        "IRC follow-up",
        "IRC from my TS at /tmp/cope_ts.log on chemnode1, 16 cores",
        ["chemsmart", "sub", "gaussian", "irc"],
    ),
    (
        "NCI plot",
        "NCI analysis of dimer.wfn locally",
        # Raw .wfn input goes through the `nciplot` subcommand, NOT
        # `gaussian nci`. See system-prompt NCI family guidance.
        ["chemsmart", "run", "nciplot", "-f", "dimer.wfn"],
    ),
    (
        "TDDFT",
        "TDDFT 5 roots of pyridine in water solvent, def2-svp, "
        "on chemnode1, 16 cores",
        ["chemsmart", "sub", "gaussian", "td", "-b", "def2-svp"],
    ),
    (
        "ORCA opt",
        "ORCA opt of caffeine, b3lyp/def2-svp basis, run locally on 4 cores",
        ["chemsmart", "orca", "opt", "-b", "def2-svp"],
    ),
    (
        "PubChem input",
        "build benzene from PubChem and optimize with b3lyp/6-31g* locally",
        ["chemsmart", "gaussian", "opt", "--pubchem", "benzene"],
    ),
    (
        "Boltzmann thermochem",
        "Boltzmann weighting of conformers in conformers/ directory",
        ["chemsmart", "thermochemistry", "boltzmann"],
    ),
    (
        "modredundant scan",
        "modredundant scan of bond between atoms 1 and 2 from 1.0 to 2.0 "
        "in 0.1 steps for oxetane.xyz",
        ["chemsmart", "gaussian", "modred"],
    ),
]


_KEY_PATTERN = re.compile(r"(api[_-]?key[^=:]*[=:]\s*)([^\s\"']+)")


def _sanitize(text: str) -> str:
    return _KEY_PATTERN.sub(r"\1<REDACTED>", text)


def _evaluate_one(session, prompt: str, expected: list[str]) -> dict:
    """Run one synthesis case and return a result dict."""
    record = {
        "prompt": prompt,
        "expected": expected,
        "status": "error",
        "command": "",
        "valid": False,
        "validation_error": "",
        "missing": [],
        "ok": False,
        "exception": "",
    }
    try:
        result = session.synthesize(prompt)
    except Exception as exc:  # provider failure / parse failure
        record["exception"] = f"{type(exc).__name__}: {exc}"
        return record

    record["status"] = str(result.get("status", ""))
    record["command"] = str(result.get("command", ""))

    if record["status"] != "ready":
        return record

    valid, err = session.validate_command(record["command"])
    record["valid"] = bool(valid)
    record["validation_error"] = err

    missing = [tok for tok in expected if tok not in record["command"]]
    record["missing"] = missing
    record["ok"] = valid and not missing
    return record


def main() -> int:
    console = Console()

    from chemsmart.agent.providers import get_provider
    from chemsmart.agent.synthesis import SynthesisSession

    try:
        provider = get_provider()
    except Exception as exc:
        console.print(f"[red]Provider unavailable:[/red] {_sanitize(str(exc))}")
        return 2

    console.print(
        f"Provider type: [cyan]{getattr(provider, 'name', type(provider).__name__)}[/cyan]"
    )
    session = SynthesisSession(provider=provider)

    records: list[dict] = []
    table = Table(title="chemsmart agent ask synthesis evaluation")
    table.add_column("#", justify="right")
    table.add_column("category")
    table.add_column("synthesized command", overflow="fold")
    table.add_column("missing tokens", overflow="fold")
    table.add_column("verdict")

    for index, (category, prompt, expected) in enumerate(CASES, start=1):
        record = _evaluate_one(session, prompt, expected)
        record["category"] = category
        records.append(record)

        if record["exception"]:
            verdict = f"[red]ERR[/red] {record['exception']}"
        elif not record["status"] == "ready":
            verdict = f"[red]FAIL[/red] status={record['status']}"
        elif not record["valid"]:
            verdict = f"[red]FAIL[/red] invalid: {record['validation_error']}"
        elif record["missing"]:
            verdict = f"[red]FAIL[/red] missing={record['missing']}"
        else:
            verdict = "[green]PASS[/green]"

        table.add_row(
            str(index),
            category,
            record["command"] or "(no command)",
            ", ".join(record["missing"]) if record["missing"] else "—",
            verdict,
        )

    passed = sum(1 for r in records if r["ok"])
    console.print(table)
    console.print(
        f"\n[bold]{passed}/{len(records)} passed[/bold]"
    )

    out_path = Path("/tmp/eval_synthesis_run.txt")
    lines = ["chemsmart agent ask synthesis evaluation\n", "=" * 60 + "\n"]
    for index, record in enumerate(records, start=1):
        lines.append(f"\nCase {index}: {record['category']}\n")
        lines.append(f"  prompt: {record['prompt']}\n")
        lines.append(f"  expected_substrings: {record['expected']}\n")
        lines.append(f"  status: {record['status']}\n")
        lines.append(f"  command: {_sanitize(record['command'])}\n")
        if record["exception"]:
            lines.append(f"  exception: {record['exception']}\n")
        if record["validation_error"]:
            lines.append(f"  validation_error: {record['validation_error']}\n")
        if record["missing"]:
            lines.append(f"  missing: {record['missing']}\n")
        lines.append(f"  PASS: {record['ok']}\n")
    lines.append(f"\nSummary: {passed}/{len(records)} passed\n")
    out_path.write_text("".join(lines), encoding="utf-8")
    console.print(f"\nFull report: {out_path}")

    return 0 if passed == len(records) else 1


if __name__ == "__main__":
    sys.exit(main())
