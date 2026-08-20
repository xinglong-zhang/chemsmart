"""Credentials resolve without a flag: environment first, then the store.

The one-use lease and its provider binding are unchanged; what changed is
where the secret may come from when no file is named. An exported variable
matching the profile's label wins; otherwise the managed key store that
`chemsmart config agent` writes; otherwise a refusal that names both
routes.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.api_access import (
    DEFAULT_KEY_LABELS,
    PROVIDER_KEY_LABEL_TOKENS,
    default_agent_keys_path,
    load_secret_lease,
)


def test_an_exported_variable_wins(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("CHEMSMART_AGENT_KEYS", raising=False)
    monkeypatch.setenv("OPENAI_API_KEY", "sk-env-wins")

    lease = load_secret_lease(
        provider="openai", label="OPENAI_API_KEY", ttl_seconds=30
    )

    assert lease.invoke(lambda secret: secret) == "sk-env-wins"


def test_the_managed_store_is_the_fallback(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("CHEMSMART_AGENT_KEYS", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    store = tmp_path / ".chemsmart" / "agent" / "keys.env"
    store.parent.mkdir(parents=True)
    store.write_text("OPENAI_API_KEY=sk-from-store\n", encoding="utf-8")

    lease = load_secret_lease(
        provider="openai", label="OPENAI_API_KEY", ttl_seconds=30
    )

    assert lease.invoke(lambda secret: secret) == "sk-from-store"
    assert default_agent_keys_path() == store


def test_nothing_found_names_both_routes(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.delenv("CHEMSMART_AGENT_KEYS", raising=False)
    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

    with pytest.raises(ContractError) as refusal:
        load_secret_lease(
            provider="openai", label="OPENAI_API_KEY", ttl_seconds=30
        )

    message = str(refusal.value)
    assert "export OPENAI_API_KEY" in message
    assert "chemsmart config agent" in message


def test_an_explicit_file_is_still_honoured(tmp_path):
    named = tmp_path / "explicit.env"
    named.write_text("OPENAI_API_KEY=sk-explicit\n", encoding="utf-8")

    lease = load_secret_lease(
        provider="openai",
        path=named,
        label="OPENAI_API_KEY",
        ttl_seconds=30,
    )

    assert lease.invoke(lambda secret: secret) == "sk-explicit"


def test_the_new_providers_are_scrub_protected():
    assert "OPENAI_API_KEY" in DEFAULT_KEY_LABELS["openai"]
    assert "ANTHROPIC_API_KEY" in DEFAULT_KEY_LABELS["anthropic"]
    assert PROVIDER_KEY_LABEL_TOKENS["openai"] == "OPENAI"
    assert PROVIDER_KEY_LABEL_TOKENS["anthropic"] == "ANTHROPIC"
