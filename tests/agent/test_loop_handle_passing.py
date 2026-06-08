from __future__ import annotations

import json

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop

from ._loop_helpers import (
    DummyMolecule,
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


class HandleAwareProvider:
    name = "openai"
    default_model = "gpt-5.4-mock"

    def __init__(self):
        self.calls = []

    def chat(self, messages, tools=None, timeout_s=30):
        self.calls.append(
            {"messages": messages, "tools": tools, "timeout_s": timeout_s}
        )
        if len(self.calls) == 1:
            return openai_tool_call_response(
                tool_call(
                    "call_1", "build_molecule", {"filepath": "water.xyz"}
                )
            )

        handle_payload = json.loads(messages[-1]["content"])
        handle_id = handle_payload["handle_id"]
        if len(self.calls) == 2:
            return openai_tool_call_response(
                tool_call(
                    "call_2",
                    "build_job",
                    {"kind": "gaussian.opt", "molecule": handle_id},
                )
            )

        return openai_final_response("Prepared the job.")


class DummyJob:
    def __init__(self, molecule):
        self.molecule = molecule

    def __repr__(self) -> str:
        return f"DummyJob({self.molecule!r})"


def test_tool_loop_resolves_handle_ids_into_original_objects(tmp_path):
    molecule = DummyMolecule("water")
    provider = HandleAwareProvider()
    registry = ScriptedRegistry(
        {
            "build_molecule": molecule,
            "build_job": lambda args: DummyJob(args["molecule"]),
        }
    )
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Prepare water."}],
        tool_defs=registry.openai_tool_defs(),
    )

    first_outcome, second_outcome = result["tool_outcomes"]
    assert first_outcome.handle_id is not None
    assert first_outcome.handle_id.startswith("mol_")
    assert registry.calls[1][0] == "build_job"
    assert registry.calls[1][1]["molecule"] is molecule
    assert second_outcome.handle_id is not None
    assert second_outcome.handle_id.startswith("job_")
