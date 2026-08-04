from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.adaptive_api_campaign import AdaptiveApiCampaignPolicyV1
from chemsmart.agent.api_access import load_secret_lease, parse_secret_file
from chemsmart.agent.permissions import (
    PermissionLedgerV1,
    build_approval_resolution,
    build_permission_request,
)


def test_api_env_hyphen_label_and_whitespace_are_parsed_without_sourcing(
    tmp_path,
):
    path = tmp_path / "api.env"
    label = "DEEPSEEK" + "-api-key"
    path.write_text(f'{label}   =   "private-value"\n')

    parsed = parse_secret_file(path)
    lease = load_secret_lease(provider="deepseek", path=path)

    assert parsed["DEEPSEEK-api-key"] == "private-value"
    assert "private-value" not in repr(lease)
    assert lease.invoke(lambda secret: len(secret)) == len("private-value")
    with pytest.raises(ContractError, match="consumed"):
        lease.invoke(lambda secret: secret)


def test_in_memory_permission_ledger_refuses_material_approval_consume():
    request = build_permission_request(
        request_id="execute-1",
        action="execute_local",
        scope_sha256="a" * 64,
        command_sha256="b" * 64,
        input_sha256s=("c" * 64,),
        project_sha256="d" * 64,
        environment_sha256="e" * 64,
    )
    approval = build_approval_resolution(
        approval_id="approval-1",
        request=request,
        allow=True,
        actor="user",
    )
    ledger = PermissionLedgerV1()

    with pytest.raises(ContractError, match="persistent"):
        ledger.resolve(request, approval=approval)


def test_alibaba_token_plan_lease_selects_only_its_exact_label(tmp_path):
    path = tmp_path / "api.env"
    path.write_text(
        "ALIBABA_TOKEN_PLAN_KEY=sk-sp-token-plan\n"
        "DEEPSEEK-api-key=deepseek-secret\n",
        encoding="utf-8",
    )

    lease = load_secret_lease(provider="alibaba-token-plan", path=path)

    assert lease.provider == "alibaba-token-plan"
    assert "sk-sp-token-plan" not in repr(lease)
    assert lease.invoke(lambda secret: secret.startswith("sk-sp-")) is True


def test_secret_lease_rejects_normalized_label_collision(tmp_path):
    path = tmp_path / "api.env"
    path.write_text(
        "DEEPSEEK-api-key=first\nDEEPSEEK_API_KEY=second\n",
        encoding="utf-8",
    )

    with pytest.raises(ContractError, match="ambiguous"):
        load_secret_lease(provider="deepseek", path=path)


def test_adaptive_campaign_has_no_transport_call_cap_or_topup():
    policy = AdaptiveApiCampaignPolicyV1()

    assert policy.transport_attempt_limit is None
    assert policy.top_up_allowed is False
    assert policy.unique_hypothesis_required is True
