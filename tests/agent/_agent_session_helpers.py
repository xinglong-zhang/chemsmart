from __future__ import annotations

import json
from typing import Any


class FakeProvider:
    def __init__(self, responses: list[dict[str, Any]]):
        self._responses = list(responses)
        self.calls: list[dict[str, Any]] = []

    def chat(self, messages, tools=None):
        self.calls.append({"messages": messages, "tools": tools})
        response = self._responses.pop(0)
        return {"content": json.dumps(response)}

    def ping(self):
        return {
            "ok": True,
            "resolved_model": "gpt-5.4-mock",
            "latency_ms": 0,
        }


def planner_plan(filepath: str, label: str, final_tool: str | None = None):
    steps = [
        {
            "tool": "build_molecule",
            "args": {"filepath": filepath},
            "rationale": "Load the structure.",
        },
        {
            "tool": "build_gaussian_settings",
            "args": {"functional": "b3lyp", "basis": "6-31g*"},
            "rationale": "Use a simple DFT level.",
        },
        {
            "tool": "build_job",
            "args": {
                "kind": "gaussian.opt",
                "molecule": "$step1",
                "settings": "$step2",
                "label": label,
            },
            "rationale": "Assemble the Gaussian optimization job.",
        },
        {
            "tool": "dry_run_input",
            "args": {"job": "$step3"},
            "rationale": "Write the input file for inspection.",
        },
        {
            "tool": "validate_runtime",
            "args": {"job": "$step3", "server": None},
            "rationale": "Check local/runtime prerequisites.",
        },
    ]
    if final_tool is not None:
        steps.append(
            {
                "tool": final_tool,
                "args": {"job": "$step3"},
                "rationale": f"Execute {final_tool} after approval.",
            }
        )
    return {
        "steps": steps,
        "rationale": "Prepare a single Gaussian optimization workflow.",
        "estimated_cost": "low",
    }


def critic_ok():
    return {
        "verdict": "ok",
        "issues": [],
        "rationale": "The generated input looks reasonable.",
    }


def critic_warn(issue: str):
    return {
        "verdict": "warn",
        "issues": [issue],
        "rationale": "Proceed only with explicit override.",
    }


def critic_reject(issue: str):
    return {
        "verdict": "reject",
        "issues": [issue],
        "rationale": "Do not execute this workflow.",
    }
