"""A profile's declared ``api_key_env`` must be the key that gets billed.

Two defects, found while pointing a session at a second DeepSeek key:

* ``_build_profile`` enumerated three accepted labels, so a lab key was
  refused as "unapproved" purely for being new -- and the refusal did not say
  what would be approved.
* ``load_secret_lease`` then ignored ``api_key_env`` altogether. It selected
  the first entry of ``DEFAULT_KEY_LABELS`` present in the file, so an account
  holding both a personal and a lab key billed whichever label sorted first,
  silently, whatever the profile said.

The second is the dangerous one: the configuration looked correct, the session
ran, and the wrong account paid. Nothing in the run recorded which key was used
until ``SecretLease.label`` was inspected.

No test here reads a secret value; only labels are asserted.
"""

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.api_access import (
    load_secret_lease,
    require_provider_key_label,
)


@pytest.fixture
def two_key_file(tmp_path):
    """A file holding a default-priority key and a lab key for one provider."""

    path = tmp_path / "api.env"
    path.write_text(
        "DEEPSEEK-api-key=default-priority-value\n"
        "LAB-DEEPSEEK-API-KEY=lab-value\n",
        encoding="utf-8",
    )
    return path


def test_declared_label_wins_over_default_priority_order(two_key_file):
    """The regression that made a session bill the wrong account."""

    lease = load_secret_lease(
        provider="deepseek",
        path=two_key_file,
        label="LAB-DEEPSEEK-API-KEY",
        ttl_seconds=30,
    )
    assert lease.label == "LAB-DEEPSEEK-API-KEY"


def test_omitting_the_label_keeps_the_previous_behaviour(two_key_file):
    """Callers with no profile to consult must not change behaviour."""

    lease = load_secret_lease(
        provider="deepseek", path=two_key_file, ttl_seconds=30
    )
    assert lease.label == "DEEPSEEK-api-key"


def test_a_new_key_for_the_same_provider_is_accepted():
    """An enumeration cannot express a second key; a token rule can."""

    for label in (
        "DEEPSEEK-api-key",
        "DEEPSEEK_API_KEY",
        "LAB-DEEPSEEK-API-KEY",
        "deepseek_group_2_key",
    ):
        require_provider_key_label(label, provider="deepseek")


def test_another_providers_key_is_still_refused():
    """The hazard the enumeration existed for must survive generalising it."""

    for provider, label in (
        ("deepseek", "OPEN-AI-api-key"),
        ("deepseek", "ALIBABA_TOKEN_PLAN_KEY"),
        ("alibaba-token-plan", "DEEPSEEK-api-key"),
    ):
        with pytest.raises(ContractError) as excinfo:
            require_provider_key_label(label, provider=provider)
        # The refusal must say what would be accepted, not only that this
        # was not.
        assert provider in str(excinfo.value)


def test_a_missing_label_names_the_label_that_was_missing(tmp_path):
    """Regression: a shadowed loop variable named the file's last key."""

    path = tmp_path / "api.env"
    # The requested label is absent, and the last line is a *different*
    # provider -- which is what made the shadowing bug visible.
    path.write_text(
        "DEEPSEEK-api-key=value\nTAVILY_API_KEY=other\n", encoding="utf-8"
    )
    with pytest.raises(ContractError) as excinfo:
        load_secret_lease(
            provider="deepseek",
            path=path,
            label="LAB-DEEPSEEK-API-KEY",
            ttl_seconds=30,
        )
    assert "LAB-DEEPSEEK-API-KEY" in str(excinfo.value)
    assert "TAVILY" not in str(excinfo.value)


def test_the_live_session_passes_the_profile_label_through():
    """The fix is worthless if the one production caller keeps dropping it."""

    import inspect

    from chemsmart.agent.live_session import run_live_agent_session

    source = inspect.getsource(run_live_agent_session)
    assert "label=profile.api_key_env" in source, (
        "the profile declares which key it bills; load_secret_lease must be "
        "told, or it falls back to priority order"
    )
