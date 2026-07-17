from __future__ import annotations

import sys
import types

import pytest


class _FakeModel:
    def __init__(self) -> None:
        self.moves: list[str] = []
        self.eval_called = False

    def to(self, device: str):
        self.moves.append(device)
        return self

    def eval(self) -> None:
        self.eval_called = True


@pytest.mark.parametrize(
    ("runtime", "dtype", "device_map", "moves"),
    [
        ("cuda", "bf16", "auto", []),
        ("mps", "fp16", None, ["mps"]),
        ("cpu", "fp32", None, []),
    ],
)
def test_merged_transformers_loader_preserves_platform_paths(
    monkeypatch,
    runtime,
    dtype,
    device_map,
    moves,
):
    model = _FakeModel()
    model_calls: list[dict] = []

    class FakeAutoModel:
        @classmethod
        def from_pretrained(cls, _model_id, **kwargs):
            model_calls.append(kwargs)
            return model

    class FakeAutoTokenizer:
        @classmethod
        def from_pretrained(cls, model_id, **kwargs):
            return {"model_id": model_id, **kwargs}

    fake_torch = types.SimpleNamespace(
        bfloat16="bf16",
        float16="fp16",
        float32="fp32",
    )
    fake_transformers = types.SimpleNamespace(
        AutoModelForCausalLM=FakeAutoModel,
        AutoTokenizer=FakeAutoTokenizer,
    )
    monkeypatch.setitem(sys.modules, "torch", fake_torch)
    monkeypatch.setitem(sys.modules, "transformers", fake_transformers)

    from chemsmart.agent.local import loader

    monkeypatch.setattr(loader, "_enforce_huggingface_hub_compat", lambda: None)
    bundle = loader.load_transformers_model(
        base_model_id="lab/merged-model",
        runtime=runtime,
    )

    assert model_calls == [
        {
            "torch_dtype": dtype,
            "device_map": device_map,
            "token": None,
        }
    ]
    assert model.moves == moves
    assert model.eval_called is True
    assert bundle.base_model_id == "lab/merged-model"
    assert bundle.model_repo_id == "lab/merged-model"


def test_legacy_loader_alias_loads_only_merged_models():
    from chemsmart.agent.local import loader

    assert loader.load_lora_model is loader.load_transformers_model
