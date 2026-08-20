"""Non-sourcing secret-file parser and one-use credential leases."""

from __future__ import annotations

import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, TypeVar

from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.providers import PROVIDERS

_LABEL = re.compile(r"^[A-Za-z0-9_-]+$")
_T = TypeVar("_T")

#: Instrument credentials have no runtime adapter; they live beside the
#: provider-derived entries so one mapping still answers every label.
_INSTRUMENT_KEY_LABELS = {
    "elsevier": ("ELSEVIER_API_KEY", "ELSIVIER_api_key"),
    "serpapi": ("SERPAPI_API_KEY", "SerpApi_api_key"),
    "tavily": ("TAVILY_API_KEY", "Tavily_api_key"),
}

_INSTRUMENT_KEY_LABEL_TOKENS = {
    "elsevier": "ELS",
    "serpapi": "SERPAPI",
    "tavily": "TAVILY",
}

DEFAULT_KEY_LABELS = {
    **{
        name: declaration.key_labels
        for name, declaration in PROVIDERS.items()
    },
    **_INSTRUMENT_KEY_LABELS,
}

#: The substring a provider's key label must carry.  Provider entries
#: derive from the registry so the engine-env scrub set cannot lag a
#: newly registered provider.
PROVIDER_KEY_LABEL_TOKENS = {
    **{
        name: declaration.key_label_token
        for name, declaration in PROVIDERS.items()
    },
    **_INSTRUMENT_KEY_LABEL_TOKENS,
}


def normalize_key_label(label: str) -> str:
    """Fold separators and case so ``A-b_C`` and ``a_B-c`` compare equal."""

    return re.sub(r"[-_]", "", str(label)).upper()


def require_provider_key_label(label: str, *, provider: str) -> None:
    """Require a key label to name the provider it will bill.

    An enumeration of accepted labels cannot express a second key for the same
    provider: a lab, team or project key was refused purely for being new, and
    the refusal did not say what would be accepted.  The hazard the
    enumeration defended against is narrower than it looks -- pointing one
    provider's profile at another provider's secret, which bills the wrong
    account and sends the key to the wrong endpoint -- and carrying the
    provider's name is enough to prevent exactly that.
    """

    token = PROVIDER_KEY_LABEL_TOKENS.get(str(provider).lower())
    if token is None:
        return
    if token not in normalize_key_label(label):
        raise ContractError(
            f"a {provider} profile must name a {provider} key: label "
            f"{str(label)!r} does not contain {token!r}. Any prefix or "
            "suffix is accepted, so a lab, team or project key is fine, but "
            "the label must identify the provider it bills."
        )


class SecretLease:
    """One-use in-memory secret whose repr never exposes its value."""

    def __init__(
        self,
        *,
        provider: str,
        label: str,
        secret: str,
        expires_at_monotonic: float,
        clock: Callable[[], float] = time.monotonic,
    ) -> None:
        self.provider = provider
        self.label = label
        self._secret = secret
        self._expires_at = expires_at_monotonic
        self._clock = clock
        self._consumed = False

    def __repr__(self) -> str:
        return f"SecretLease(provider={self.provider!r}, value=<redacted>)"

    def invoke(self, operation: Callable[[str], _T]) -> _T:
        if self._consumed:
            raise ContractError("credential lease has already been consumed")
        if self._clock() >= self._expires_at:
            self._secret = ""
            self._consumed = True
            raise ContractError("credential lease has expired")
        secret = self._secret
        try:
            return operation(secret)
        finally:
            self._secret = ""
            self._consumed = True


def parse_secret_file(path: str | Path) -> dict[str, str]:
    """Parse simple assignments as data; never source or expand the file."""

    result = {}
    for ordinal, line in enumerate(
        Path(path).read_text(encoding="utf-8").splitlines(), start=1
    ):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("export "):
            stripped = stripped[7:].lstrip()
        if "=" not in stripped:
            raise ContractError(
                f"secret file line {ordinal} is not an assignment"
            )
        label, value = stripped.split("=", 1)
        label = label.strip()
        value = value.strip()
        if _LABEL.fullmatch(label) is None:
            raise ContractError(
                f"secret file line {ordinal} has an invalid label"
            )
        if (
            len(value) >= 2
            and value[:1] == value[-1:]
            and value[0] in {'"', "'"}
        ):
            value = value[1:-1]
        if not value:
            raise ContractError(f"secret file label {label!r} has no value")
        result[label] = value
    return result


def default_agent_keys_path() -> Path:
    """The managed key store `chemsmart config agent` writes."""

    explicit = os.environ.get("CHEMSMART_AGENT_KEYS", "").strip()
    if explicit:
        return Path(explicit).expanduser()
    return Path.home() / ".chemsmart" / "agent" / "keys.env"


def _environment_values(labels: tuple[str, ...]) -> dict[str, str]:
    """Exported credentials matching the candidate labels, exact or folded."""

    wanted = {normalize_key_label(label): label for label in labels}
    found: dict[str, str] = {}
    for name, value in os.environ.items():
        label = wanted.get(normalize_key_label(name))
        if label is not None and value.strip():
            found[label] = value.strip()
    return found


def load_secret_lease(
    *,
    provider: str,
    path: str | Path | None = None,
    label: str | None = None,
    ttl_seconds: float = 60.0,
    clock: Callable[[], float] = time.monotonic,
) -> SecretLease:
    """Lease one secret, preferring the label the caller actually declared.

    ``label`` exists because a provider profile declares ``api_key_env`` and
    this loader used to ignore it, selecting instead the first entry of
    ``DEFAULT_KEY_LABELS`` present in the file.  An account with two keys --
    a personal key and a lab key, say -- therefore billed whichever label
    happened to sort first, silently, no matter which one the profile named.
    A declared label is now honoured exactly; omitting it keeps the previous
    priority-order behaviour for callers that have no profile to consult.
    """

    if ttl_seconds <= 0:
        raise ContractError("credential lease TTL must be positive")
    if label is not None and not str(label).strip():
        raise ContractError("requested secret label must not be blank")
    requested = None if label is None else str(label)
    if requested is not None:
        require_provider_key_label(requested, provider=provider)
    labels = (
        (requested,)
        if requested is not None
        else DEFAULT_KEY_LABELS.get(str(provider).lower(), ())
    )
    if not labels:
        raise ContractError("provider has no approved secret labels")
    if path is not None:
        values = parse_secret_file(path)
    else:
        # No file was named: an exported environment variable wins, then
        # the managed key store `chemsmart config agent` maintains.
        values = _environment_values(tuple(labels))
        if not values:
            store = default_agent_keys_path()
            if store.is_file():
                values = parse_secret_file(store)
        if not values:
            raise ContractError(
                f"no credential found for provider {provider!r}: export "
                f"{labels[0]} in your shell, or store it with "
                "'chemsmart config agent'"
            )
    normalized = {}
    for file_label, value in values.items():
        normalized_label = _normalize_label(file_label)
        if normalized_label in normalized:
            raise ContractError("secret file contains ambiguous key labels")
        normalized[normalized_label] = (file_label, value)
    selected = next(
        (
            normalized[_normalize_label(candidate)]
            for candidate in labels
            if _normalize_label(candidate) in normalized
        ),
        None,
    )
    if selected is None:
        if requested is not None:
            raise ContractError(
                f"secret file has no entry for the requested label "
                f"{requested!r}; the profile's api_key_env names the key it "
                "bills, so the label must exist in the secret file"
            )
        raise ContractError("approved provider key label was not found")
    label, secret = selected
    return SecretLease(
        provider=str(provider).lower(),
        label=label,
        secret=secret,
        expires_at_monotonic=clock() + ttl_seconds,
        clock=clock,
    )


def _normalize_label(label: str) -> str:
    # One implementation, so matching and provider-token checking can never
    # disagree about whether two spellings are the same label.
    return normalize_key_label(label).lower()


__all__ = [
    "DEFAULT_KEY_LABELS",
    "SecretLease",
    "load_secret_lease",
    "parse_secret_file",
]
