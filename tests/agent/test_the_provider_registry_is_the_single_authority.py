"""Every provider surface derives from the one registry declaration."""

from chemsmart.agent.api_access import (
    DEFAULT_KEY_LABELS,
    PROVIDER_KEY_LABEL_TOKENS,
)
from chemsmart.agent.provider_config import (
    ALIBABA_TOKEN_PLAN_ENDPOINT,
    ANTHROPIC_OFFICIAL_ENDPOINT,
    OPENAI_OFFICIAL_ENDPOINT,
)
from chemsmart.agent.providers import (
    PROVIDERS,
    declaration_for_endpoint,
    runnable_provider_names,
)


def test_the_four_declarations_carry_the_registered_identities():
    assert tuple(PROVIDERS) == (
        "openai",
        "anthropic",
        "alibaba-token-plan",
        "deepseek",
    )
    assert runnable_provider_names() == {
        "openai",
        "alibaba-token-plan",
        "deepseek",
    }
    refusing = PROVIDERS["anthropic"]
    assert not refusing.runnable
    assert "not registered in this release" in refusing.refusal


def test_key_tables_and_endpoints_derive_from_the_registry():
    for name, declaration in PROVIDERS.items():
        assert DEFAULT_KEY_LABELS[name] == declaration.key_labels
        assert PROVIDER_KEY_LABEL_TOKENS[name] == declaration.key_label_token
    assert OPENAI_OFFICIAL_ENDPOINT is PROVIDERS["openai"].endpoint
    assert ANTHROPIC_OFFICIAL_ENDPOINT is PROVIDERS["anthropic"].endpoint
    assert (
        ALIBABA_TOKEN_PLAN_ENDPOINT is PROVIDERS["alibaba-token-plan"].endpoint
    )
    assert declaration_for_endpoint("https://api.example.com") is None


def test_the_config_command_menu_matches_the_registry():
    from chemsmart.cli.config import _AGENT_PROVIDERS

    assert _AGENT_PROVIDERS == tuple(PROVIDERS)
    for name in _AGENT_PROVIDERS:
        assert PROVIDERS[name].cli_efforts, name
        assert PROVIDERS[name].model_hint, name
