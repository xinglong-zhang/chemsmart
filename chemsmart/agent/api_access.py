"""Non-sourcing secret-file parser and one-use credential leases."""

from __future__ import annotations

import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, TypeVar

from chemsmart.agent._contracts import ContractError, canonical_sha256


_LABEL = re.compile(r"^[A-Za-z0-9_-]+$")
_T = TypeVar("_T")

DEFAULT_KEY_LABELS = {
    "deepseek": (
        "DEEPSEEK-api-key",
        "DEEPSEEK_API_KEY",
        "DEEPSEEK-API-KEY",
    ),
    "elsevier": ("ELSEVIER_API_KEY", "ELSIVIER_api_key"),
    "serpapi": ("SERPAPI_API_KEY", "SerpApi_api_key"),
    "tavily": ("TAVILY_API_KEY", "Tavily_api_key"),
}


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


@dataclass(frozen=True)
class ApiCredentialStatusReceiptV1:
    schema_version: str
    provider: str
    key_validation_status: str
    quota_status: str
    entitlement_status: str
    nonsecret_error_class: str
    receipt_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != "chemsmart.api-credential-status.v1":
            raise ContractError("unsupported API credential receipt schema")
        body = {
            "schema_version": self.schema_version,
            "provider": self.provider,
            "key_validation_status": self.key_validation_status,
            "quota_status": self.quota_status,
            "entitlement_status": self.entitlement_status,
            "nonsecret_error_class": self.nonsecret_error_class,
        }
        if self.receipt_sha256 != canonical_sha256(body):
            raise ContractError("API credential status digest mismatch")


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
            raise ContractError(f"secret file line {ordinal} is not an assignment")
        label, value = stripped.split("=", 1)
        label = label.strip()
        value = value.strip()
        if _LABEL.fullmatch(label) is None:
            raise ContractError(f"secret file line {ordinal} has an invalid label")
        if len(value) >= 2 and value[:1] == value[-1:] and value[0] in {'"', "'"}:
            value = value[1:-1]
        if not value:
            raise ContractError(f"secret file label {label!r} has no value")
        result[label] = value
    return result


def load_secret_lease(
    *,
    provider: str,
    path: str | Path,
    ttl_seconds: float = 60.0,
    clock: Callable[[], float] = time.monotonic,
) -> SecretLease:
    if ttl_seconds <= 0:
        raise ContractError("credential lease TTL must be positive")
    labels = DEFAULT_KEY_LABELS.get(str(provider).lower(), ())
    if not labels:
        raise ContractError("provider has no approved secret labels")
    values = parse_secret_file(path)
    normalized = {
        _normalize_label(label): (label, value)
        for label, value in values.items()
    }
    selected = next(
        (
            normalized[_normalize_label(label)]
            for label in labels
            if _normalize_label(label) in normalized
        ),
        None,
    )
    if selected is None:
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
    return re.sub(r"[-_]", "", str(label)).lower()


__all__ = [
    "ApiCredentialStatusReceiptV1",
    "DEFAULT_KEY_LABELS",
    "SecretLease",
    "load_secret_lease",
    "parse_secret_file",
]
