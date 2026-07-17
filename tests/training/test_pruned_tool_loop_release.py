from __future__ import annotations

import json

from chemsmart.agent.registry import ToolRegistry
from scripts.training.build_pruned_tool_loop_release import (
    project_tool_loop_record,
)


class WordTokenizer:
    def apply_chat_template(self, conversation, *, tools, **kwargs):
        del kwargs
        text = json.dumps(
            {"messages": conversation, "tools": tools}, ensure_ascii=False
        )
        return text.split()


class MappingTokenizer(WordTokenizer):
    def apply_chat_template(self, conversation, *, tools, **kwargs):
        return {
            "input_ids": super().apply_chat_template(
                conversation, tools=tools, **kwargs
            ),
            "attention_mask": [],
        }


def _record(tool_content: dict, *, system: str = "old prompt"):
    return {
        "messages": [
            {"role": "system", "content": system},
            {"role": "user", "content": "optimize h2o.xyz"},
            {
                "role": "assistant",
                "content": "",
                "tool_calls": [
                    {
                        "id": "call-1",
                        "type": "function",
                        "function": {
                            "name": "synthesize_command",
                            "arguments": '{"request":"optimize h2o.xyz"}',
                        },
                    }
                ],
            },
            {
                "role": "tool",
                "tool_call_id": "call-1",
                "content": json.dumps(tool_content),
            },
            {"role": "assistant", "content": "Prepared the command."},
        ],
        "meta": {"family": "tool_loop", "session_id": "s1"},
    }


def test_projection_removes_raw_reasoning_and_preserves_evidence():
    record = _record(
        {
            "ok": True,
            "status": "ready",
            "command": "chemsmart run gaussian -f h2o.xyz opt",
            "raw_response": "very long duplicate",
            "reasoning": "private provider trace",
            "workflow_state": {"cwd": "/Users/person/work"},
            "semantic": {
                "verdict": "ok",
                "failed_rule_ids": [],
                "generated_inputs": [{"path": "/Users/person/a.com"}],
            },
        }
    )
    rows, reason = project_tool_loop_record(
        record,
        registry=ToolRegistry.default(groups=["synthesis"]),
        tokenizer=WordTokenizer(),
    )
    assert reason == "ok"
    assert len(rows) == 1
    serialized = json.dumps(rows[0])
    assert "raw_response" not in serialized
    assert "private provider trace" not in serialized
    assert "workflow_state" not in serialized
    assert "<workspace-home>" in serialized
    assert rows[0]["meta"]["qwen3_14b_tokens"] <= 8192
    assert len(rows[0]["messages"][0]["content"]) <= 4096


def test_full_schema_and_unprovenanced_reasoning_are_excluded():
    rows, reason = project_tool_loop_record(
        _record({}, system="x" * 100_000),
        registry=ToolRegistry.default(groups=["synthesis"]),
        tokenizer=WordTokenizer(),
    )
    assert rows == []
    assert reason == "full_schema_system_prompt"

    rows, reason = project_tool_loop_record(
        _record({"reasoning": "x" * 8193}),
        registry=ToolRegistry.default(groups=["synthesis"]),
        tokenizer=WordTokenizer(),
    )
    assert rows == []
    assert reason == "unprovenanced_long_reasoning"


def test_projection_has_no_orphan_tool_pairs():
    record = _record({"ok": True, "status": "ready"})
    record["messages"] = [
        message
        for message in record["messages"]
        if message.get("role") != "tool"
    ]
    rows, reason = project_tool_loop_record(
        record,
        registry=ToolRegistry.default(groups=["synthesis"]),
        tokenizer=WordTokenizer(),
    )
    assert rows == []
    assert reason == "orphan_tool_call"


def test_mapping_tokenizer_counts_input_ids_not_mapping_keys():
    rows, reason = project_tool_loop_record(
        _record({"ok": True, "status": "ready"}),
        registry=ToolRegistry.default(groups=["synthesis"]),
        tokenizer=MappingTokenizer(),
    )
    assert reason == "ok"
    assert rows[0]["meta"]["qwen3_14b_tokens"] > 2
