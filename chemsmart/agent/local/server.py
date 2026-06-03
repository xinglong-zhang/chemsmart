"""OpenAI-compatible FastAPI server that fronts the V4 LoRA model.

Exposes ``POST /v1/chat/completions`` so that ``chemsmart agent ask`` can call
the local model with ``type: openai`` + a custom ``base_url`` in
``~/.chemsmart/agent/agent.yaml`` — no changes needed to provider plumbing.

Usage (Colab):

    from chemsmart.agent.local.loader import load_lora_model
    from chemsmart.agent.local.server import build_app
    import uvicorn
    bundle = load_lora_model()
    app = build_app(bundle)
    uvicorn.run(app, host="0.0.0.0", port=8000)

Then in ``~/.chemsmart/agent/agent.yaml``::

    active: local_chemsmart_v4
    providers:
      local_chemsmart_v4:
        type: openai
        model: chemsmart-qwen2.5-7b-lora
        base_url: https://<ngrok-tunnel>/v1
        api_key: local-not-used
"""

from __future__ import annotations

import json
import time
import uuid
from typing import Any

from chemsmart.agent.local.adapter import plan_to_synthesis_result
from chemsmart.agent.local.generator import generate_plan

_MODEL_ID = "chemsmart-qwen2.5-7b-lora"


def build_app(bundle: Any) -> Any:
    """Return a FastAPI app wrapping ``bundle`` with OpenAI-compatible endpoints.

    The ``/v1/chat/completions`` response wraps a SynthesisSession-shaped JSON
    object as the assistant message content, so that
    :func:`chemsmart.agent.synthesis._normalize_result` accepts it directly.
    """
    try:
        from fastapi import FastAPI, HTTPException, Request
    except ImportError as exc:
        raise RuntimeError(
            "fastapi is required to serve the local model. "
            "Install with: pip install fastapi uvicorn"
        ) from exc

    app = FastAPI(title="chemsmart-local-llm", version="1.0.0")

    @app.get("/v1/models")
    def list_models() -> dict[str, Any]:
        return {
            "object": "list",
            "data": [
                {
                    "id": _MODEL_ID,
                    "object": "model",
                    "created": int(time.time()),
                    "owned_by": "chemsmart",
                }
            ],
        }

    @app.post("/v1/chat/completions")
    async def chat_completions(request: Request) -> dict[str, Any]:
        payload = await request.json()
        messages = payload.get("messages")
        if not isinstance(messages, list) or not messages:
            raise HTTPException(status_code=400, detail="messages required")

        user_query, history = _extract_user_query(messages)
        if not user_query:
            raise HTTPException(
                status_code=400, detail="no user message found"
            )

        try:
            plan = generate_plan(bundle, user_query, history=history)
        except ValueError as exc:
            raise HTTPException(status_code=502, detail=str(exc)) from exc

        result = plan_to_synthesis_result(plan, user_query)
        content = json.dumps(result, ensure_ascii=False)
        return _wrap_openai_completion(content)

    return app


def _extract_user_query(
    messages: list[dict[str, Any]]
) -> tuple[str, list[dict[str, str]]]:
    """Return the last user message and the prior non-system turns as history.

    System prompts from the caller are discarded — the local server always
    uses the trained-on ``planner_compact.md`` prompt.
    """
    non_system: list[dict[str, str]] = []
    for msg in messages:
        if not isinstance(msg, dict):
            continue
        role = str(msg.get("role", ""))
        if role not in {"user", "assistant"}:
            continue
        content = msg.get("content", "")
        if not isinstance(content, str):
            content = json.dumps(content, ensure_ascii=False)
        non_system.append({"role": role, "content": content})

    last_user_index = -1
    for index in range(len(non_system) - 1, -1, -1):
        if non_system[index]["role"] == "user":
            last_user_index = index
            break
    if last_user_index < 0:
        return "", []

    user_query = non_system[last_user_index]["content"]
    history = non_system[:last_user_index]
    return user_query, history


def _wrap_openai_completion(content: str) -> dict[str, Any]:
    return {
        "id": f"chatcmpl-{uuid.uuid4().hex}",
        "object": "chat.completion",
        "created": int(time.time()),
        "model": _MODEL_ID,
        "choices": [
            {
                "index": 0,
                "message": {"role": "assistant", "content": content},
                "finish_reason": "stop",
            }
        ],
        "usage": {
            "prompt_tokens": 0,
            "completion_tokens": 0,
            "total_tokens": 0,
        },
    }
