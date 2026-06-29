"""Regression tests for the local v13.1 generator postprocessor wiring."""

from __future__ import annotations

import json
import sys
from types import SimpleNamespace
from typing import Any

from chemsmart.agent.local import generator, postprocessor
from chemsmart.agent.local import postprocess as exported_postprocess


def test_generator_uses_local_postprocessor() -> None:
    assert generator.postprocess is postprocessor.postprocess
    assert exported_postprocess is postprocessor.postprocess


def test_generate_plan_calls_two_argument_local_postprocessor(
    monkeypatch,
) -> None:
    class _NoGrad:
        def __enter__(self) -> None:
            return None

        def __exit__(self, *_exc: object) -> None:
            return None

    class _Inputs(dict):
        def to(self, _device: str) -> "_Inputs":
            return self

    class _InputIds:
        shape = (1, 1)

    class _Generated:
        def __getitem__(self, _item: object) -> list[int]:
            return [1, 2, 3]

    class _Tokenizer:
        def apply_chat_template(self, *_args: object, **_kwargs: object) -> str:
            return "prompt"

        def __call__(self, *_args: object, **_kwargs: object) -> _Inputs:
            return _Inputs({"input_ids": _InputIds()})

        def decode(self, *_args: object, **_kwargs: object) -> str:
            return json.dumps(
                {
                    "intent": "workflow",
                    "steps": [
                        {
                            "tool": "build_molecule",
                            "args": {"filepath": "examples/wrong.xyz"},
                        }
                    ],
                }
            )

        eos_token_id = 0

    class _Model:
        device = "cpu"

        def generate(self, **_kwargs: Any) -> list[_Generated]:
            return [_Generated()]

    monkeypatch.setitem(
        sys.modules,
        "torch",
        SimpleNamespace(no_grad=lambda: _NoGrad()),
    )

    plan = generator.generate_plan(
        SimpleNamespace(tokenizer=_Tokenizer(), model=_Model()),
        "run examples/water.xyz",
    )

    assert plan["steps"][0]["args"]["filepath"] == "examples/water.xyz"
