from __future__ import annotations

import sys
import types


def test_mlx_loader_enables_mistral_regex_fix(monkeypatch):
    calls: list[tuple[str, dict]] = []

    def fake_load(model_id, **kwargs):
        calls.append((model_id, kwargs))
        return object(), object()

    fake_mlx = types.SimpleNamespace(load=fake_load)
    monkeypatch.setitem(sys.modules, "mlx_lm", fake_mlx)

    from chemsmart.agent.local.mlx_loader import load_mlx_model

    bundle = load_mlx_model("local-model")

    assert bundle.backend == "mlx"
    assert calls == [
        (
            "local-model",
            {"tokenizer_config": {"fix_mistral_regex": True}},
        )
    ]
