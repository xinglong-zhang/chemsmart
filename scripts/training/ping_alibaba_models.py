"""Ping Alibaba/DashScope model ids without printing secrets."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from dotenv import dotenv_values

BASE_URL = "https://dashscope-intl.aliyuncs.com/compatible-mode/v1"
REPO = Path(__file__).resolve().parents[2]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("models", nargs="+")
    parser.add_argument("--timeout", type=float, default=20.0)
    args = parser.parse_args(argv)

    key = dotenv_values(REPO / "api.env").get("alibaba-cloud-api-key", "")
    if not key:
        raise SystemExit("api.env is missing alibaba-cloud-api-key")

    import openai

    client = openai.OpenAI(api_key=key, base_url=BASE_URL)
    for model in args.models:
        started = time.perf_counter()
        result = _ping_model(client, model, timeout=args.timeout)
        result["latency_s"] = round(time.perf_counter() - started, 2)
        print(
            json.dumps(result, ensure_ascii=False, sort_keys=True), flush=True
        )
    return 0


def _ping_model(client: Any, model: str, *, timeout: float) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "model": model,
        "messages": [{"role": "user", "content": "ping"}],
        "max_tokens": 4,
        "timeout": timeout,
    }
    try:
        response = client.chat.completions.create(**kwargs)
        return {
            "model": model,
            "ok": True,
            "resolved_model": getattr(response, "model", model),
        }
    except Exception as exc:  # pragma: no cover - live API helper
        message = str(exc)
        if "enable_think" in message or "enable_thinking" in message:
            try:
                retry_kwargs = dict(kwargs)
                retry_kwargs["extra_body"] = {"enable_thinking": False}
                response = client.chat.completions.create(**retry_kwargs)
                return {
                    "model": model,
                    "ok": True,
                    "resolved_model": getattr(response, "model", model),
                    "retry": "enable_thinking_false",
                }
            except Exception as retry_exc:
                message = str(retry_exc)
        return {
            "model": model,
            "ok": False,
            "error_class": type(exc).__name__,
            "error": message[:220],
            "quota": "free quota" in message.lower(),
        }


if __name__ == "__main__":
    raise SystemExit(main())
