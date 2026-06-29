from __future__ import annotations

import json

from chemsmart.agent.provider_adapter import (
    ToolOutcome,
    build_tool_result_messages,
    normalize_response,
)


def test_normalize_anthropic_text_only_response():
    response = {
        "content": [{"type": "text", "text": "Hello from Claude."}],
        "stop_reason": "end_turn",
    }

    assistant_text, requests, stop_reason = normalize_response(
        "anthropic", response
    )

    assert assistant_text == "Hello from Claude."
    assert requests == []
    assert stop_reason == "end_turn"


def test_normalize_anthropic_single_tool_call_response():
    response = {
        "content": [
            {"type": "text", "text": "I will inspect the structure."},
            {
                "type": "tool_use",
                "id": "toolu_01A",
                "name": "build_molecule",
                "input": {"filepath": "water.xyz", "index": "1"},
            },
        ],
        "stop_reason": "tool_use",
    }

    assistant_text, requests, stop_reason = normalize_response(
        "anthropic", response
    )

    assert assistant_text == "I will inspect the structure."
    assert stop_reason == "tool_use"
    assert len(requests) == 1
    request = requests[0]
    assert request.request_id == "anthropic:toolu_01A"
    assert request.provider == "anthropic"
    assert request.provider_call_id == "toolu_01A"
    assert request.name == "build_molecule"
    assert request.arguments == {"filepath": "water.xyz", "index": "1"}
    assert json.loads(request.arguments_json) == request.arguments


def test_normalize_anthropic_parallel_tool_calls_response():
    response = {
        "content": [
            {
                "type": "tool_use",
                "id": "toolu_01A",
                "name": "recommend_method",
                "input": {"task": "opt"},
            },
            {
                "type": "tool_use",
                "id": "toolu_01B",
                "name": "build_molecule",
                "input": {"filepath": "water.xyz"},
            },
        ],
        "stop_reason": "tool_use",
    }

    assistant_text, requests, stop_reason = normalize_response(
        "anthropic", response
    )

    assert assistant_text == ""
    assert stop_reason == "tool_use"
    assert [request.provider_call_id for request in requests] == [
        "toolu_01A",
        "toolu_01B",
    ]
    assert [request.name for request in requests] == [
        "recommend_method",
        "build_molecule",
    ]


def test_normalize_openai_text_only_response():
    response = {
        "choices": [
            {
                "message": {"role": "assistant", "content": "Done."},
                "finish_reason": "stop",
            }
        ]
    }

    assistant_text, requests, stop_reason = normalize_response(
        "openai", response
    )

    assert assistant_text == "Done."
    assert requests == []
    assert stop_reason == "stop"


def test_normalize_openai_single_tool_call_response():
    response = {
        "choices": [
            {
                "message": {
                    "role": "assistant",
                    "content": "I'll prepare the job.",
                    "tool_calls": [
                        {
                            "id": "call_01A",
                            "type": "function",
                            "function": {
                                "name": "build_job",
                                "arguments": (
                                    '{"kind": "gaussian.opt", '
                                    '"label": "water"}'
                                ),
                            },
                        }
                    ],
                },
                "finish_reason": "tool_calls",
            }
        ]
    }

    assistant_text, requests, stop_reason = normalize_response(
        "openai", response
    )

    assert assistant_text == "I'll prepare the job."
    assert stop_reason == "tool_calls"
    assert len(requests) == 1
    request = requests[0]
    assert request.request_id == "openai:call_01A"
    assert request.provider == "openai"
    assert request.provider_call_id == "call_01A"
    assert request.name == "build_job"
    assert request.arguments == {"kind": "gaussian.opt", "label": "water"}
    assert json.loads(request.arguments_json) == request.arguments


def test_normalize_openai_parallel_tool_calls_response():
    response = {
        "choices": [
            {
                "message": {
                    "role": "assistant",
                    "content": "",
                    "tool_calls": [
                        {
                            "id": "call_01A",
                            "type": "function",
                            "function": {
                                "name": "recommend_method",
                                "arguments": '{"task": "opt"}',
                            },
                        },
                        {
                            "id": "call_01B",
                            "type": "function",
                            "function": {
                                "name": "build_molecule",
                                "arguments": '{"filepath": "water.xyz"}',
                            },
                        },
                    ],
                },
                "finish_reason": "tool_calls",
            }
        ]
    }

    assistant_text, requests, stop_reason = normalize_response(
        "openai", response
    )

    assert assistant_text == ""
    assert stop_reason == "tool_calls"
    assert [request.provider_call_id for request in requests] == [
        "call_01A",
        "call_01B",
    ]
    assert [request.name for request in requests] == [
        "recommend_method",
        "build_molecule",
    ]


def test_build_tool_result_messages_for_anthropic():
    outcomes = [
        ToolOutcome(
            request_id="anthropic:toolu_01A",
            provider_call_id="toolu_01A",
            name="build_molecule",
            status="ok",
            result={"handle_id": "mol_abcd"},
        ),
        ToolOutcome(
            request_id="anthropic:toolu_01B",
            provider_call_id="toolu_01B",
            name="submit_hpc",
            status="denied",
            error_type="PermissionDenied",
            error_message="Denied by user",
        ),
    ]

    messages = build_tool_result_messages("anthropic", outcomes)

    assert len(messages) == 1
    assert messages[0]["role"] == "user"
    assert len(messages[0]["content"]) == 2
    ok_block, denied_block = messages[0]["content"]
    assert ok_block["tool_use_id"] == "toolu_01A"
    assert json.loads(ok_block["content"]) == {"handle_id": "mol_abcd"}
    assert "is_error" not in ok_block
    assert denied_block["tool_use_id"] == "toolu_01B"
    assert denied_block["is_error"] is True
    assert json.loads(denied_block["content"]) == {
        "ok": False,
        "error": {
            "type": "PermissionDenied",
            "message": "Denied by user",
            "tool": "submit_hpc",
        },
    }


def test_build_tool_result_messages_for_openai():
    outcomes = [
        ToolOutcome(
            request_id="openai:call_01A",
            provider_call_id="call_01A",
            name="build_job",
            status="ok",
            display_result={"handle_id": "job_abcd"},
        ),
        ToolOutcome(
            request_id="openai:call_01B",
            provider_call_id="call_01B",
            name="run_local",
            status="error",
            error_type="RuntimeError",
            error_message="Executable not found",
        ),
    ]

    messages = build_tool_result_messages("openai", outcomes)

    assert messages == [
        {
            "role": "tool",
            "tool_call_id": "call_01A",
            "content": '{"handle_id": "job_abcd"}',
        },
        {
            "role": "tool",
            "tool_call_id": "call_01B",
            "content": (
                '{"error": {"message": "Executable not found", '
                '"tool": "run_local", "type": "RuntimeError"}, '
                '"ok": false}'
            ),
        },
    ]
