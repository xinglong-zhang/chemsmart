from __future__ import annotations

from datetime import datetime, timezone
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

from chemsmart.agent._contracts import canonical_sha256, file_sha256
from chemsmart.agent.experiments.qwen_pyscf_campaign import (
    approved_xyz_source,
    bind_case_to_approved_xyz,
    build_frozen_transfer_manifest,
    build_qwen_campaign_window,
)
from chemsmart.agent.experiments.qwen_pyscf_dfc import (
    build_episode_plans,
    build_qwen_dfc_arm,
)
from chemsmart.agent.experiments.qwen_pyscf_fixtures import qwen_pyscf_cases_v1


_LAUNCHER = Path(
    "/Users/hongjiseung/.codex/runs/chemsmart-pyscf-agent-24h/"
    "campaign/run_primary_dfc.py"
)


def _load_launcher():
    spec = importlib.util.spec_from_file_location(
        "chemsmart_qwen_campaign_launcher_recovery_test",
        _LAUNCHER,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _inputs(tmp_path):
    xyz = tmp_path / "approved.xyz"
    xyz.write_text(
        "3\napproved water\n"
        "O 0.0 0.0 0.0\n"
        "H 0.0 0.757 0.586\n"
        "H 0.0 -0.757 0.586\n",
        encoding="utf-8",
    )
    source = approved_xyz_source(xyz, expected_sha256=file_sha256(xyz))
    registry = {item.case_id: item for item in qwen_pyscf_cases_v1()}
    cases = tuple(
        bind_case_to_approved_xyz(registry[case_id], source)
        for case_id in ("QP-DEV-005", "QP-TR-001")
    )
    arm = build_qwen_dfc_arm(
        decomposition=False,
        feedback_projection="full-v1",
        critic=False,
        max_concurrency=1,
    )
    plans = build_episode_plans(
        cases=cases,
        arms=(arm,),
        repeats=1,
        prompt_sha256=canonical_sha256("prompt"),
        tool_schema_sha256=canonical_sha256("tools"),
        task_order_sha256=canonical_sha256("order"),
        token_budget=1_000_000,
        tool_call_budget=256,
        wall_time_seconds=5400,
    )
    window = build_qwen_campaign_window(
        campaign_id="recovery-test",
        activation=datetime(2026, 8, 4, tzinfo=timezone.utc),
    )
    freeze = build_frozen_transfer_manifest(
        window=window,
        source=source,
        cases=cases,
        plans=plans,
    )
    implementation = SimpleNamespace(
        manifest_sha256=canonical_sha256("implementation-manifest"),
        tree_sha256=canonical_sha256("implementation-tree"),
    )
    return source, cases, plans, window, freeze, implementation


def _reset_capture(module, monkeypatch, tmp_path):
    monkeypatch.setattr(module, "DISPATCH_ROOT", tmp_path / "dispatch")
    monkeypatch.setattr(module, "PUBLIC_EVIDENCE_ROOT", tmp_path / "public")
    monkeypatch.setattr(module, "RECOVERY_LEDGER_ROOT", tmp_path / "recovery")
    monkeypatch.setattr(
        module,
        "_CAPTURED_RESULTS",
        {"development": {}, "transfer": {}},
    )
    monkeypatch.setattr(
        module,
        "_CAPTURED_GRADES",
        {"development": {}, "transfer": {}},
    )
    monkeypatch.setattr(
        module,
        "_CAPTURED_LEDGERS",
        {"development": {}, "transfer": {}},
    )


def test_dispatch_journal_prevents_a_second_provider_call(
    tmp_path, monkeypatch
):
    module = _load_launcher()
    _reset_capture(module, monkeypatch, tmp_path)
    source, cases, plans, window, freeze, implementation = _inputs(tmp_path)
    plan = next(item for item in plans if item.case_sha256 == cases[0].case_sha256)
    calls = []
    monkeypatch.setattr(
        module,
        "run_live_agent_session",
        lambda **kwargs: calls.append(kwargs) or "visible-result",
    )
    runner = module._JournaledLiveRunner(
        plans=plans,
        window=window,
        freeze=freeze,
        source=source,
        implementation_freeze=implementation,
    )
    kwargs = {
        "experiment_config": plan.experiment_config,
        "experiment_case": cases[0],
    }

    assert runner(**kwargs) == "visible-result"
    with pytest.raises(RuntimeError, match="refusing to overwrite"):
        runner(**kwargs)
    assert len(calls) == 1


def test_dispatch_without_envelope_recovers_once_as_inconclusive(
    tmp_path, monkeypatch
):
    module = _load_launcher()
    _reset_capture(module, monkeypatch, tmp_path)
    source, cases, plans, window, freeze, implementation = _inputs(tmp_path)
    plan = next(item for item in plans if item.case_sha256 == cases[0].case_sha256)
    case = cases[0]
    body = {
        "schema_version": "chemsmart.qwen-episode-dispatch-intent.v1",
        "episode_id": plan.episode_id,
        "split": case.split,
        "plan_sha256": plan.plan_sha256,
        "case_sha256": case.case_sha256,
        "hypothesis_sha256": plan.hypothesis.hypothesis_sha256,
        "experiment_config_sha256": plan.experiment_config.config_sha256,
        "campaign_window_sha256": window.window_sha256,
        "freeze_manifest_sha256": freeze.manifest_sha256,
        "source_binding_sha256": source.source_binding_sha256,
        "workspace_binding_sha256": module._workspace_binding_sha256(
            episode_id=plan.episode_id,
            source=source,
        ),
        "implementation_manifest_sha256": implementation.manifest_sha256,
        "implementation_tree_sha256": implementation.tree_sha256,
        "dispatched_at_utc": "2026-08-04T00:00:01Z",
    }
    module._write_once(
        module.DISPATCH_ROOT / f"{plan.episode_id}.json",
        {**body, "dispatch_sha256": canonical_sha256(body)},
    )

    dispatch = module._load_dispatch_records(
        plans=plans,
        window=window,
        freeze=freeze,
        source=source,
        implementation_freeze=implementation,
    )
    recovered = module._restore_episode_evidence(
        plans=plans,
        cases=cases,
        dispatch_records=dispatch,
        source=source,
    )
    ledger = module._build_recovery_split_ledger(
        split="development",
        already_covered=set(),
        dispatch_records=dispatch,
        window=window,
        freeze=freeze,
    )

    assert recovered == {plan.episode_id}
    assert ledger is not None
    episode = ledger.episode_ledgers[0]
    assert episode.failure_class == "interrupted_after_dispatch"
    assert episode.verdict == "inconclusive"
    assert episode.factor_realization_status == "not_observed"
    assert module._CAPTURED_RESULTS["development"] == {}
    assert module._require_captured_split(
        split="development",
        plans=(plan,),
        ledgers=(ledger,),
    ) == {}

    first_digest = episode.ledger_sha256
    _reset_capture(module, monkeypatch, tmp_path)
    dispatch = module._load_dispatch_records(
        plans=plans,
        window=window,
        freeze=freeze,
        source=source,
        implementation_freeze=implementation,
    )
    module._restore_episode_evidence(
        plans=plans,
        cases=cases,
        dispatch_records=dispatch,
        source=source,
    )
    assert (
        module._CAPTURED_LEDGERS["development"][plan.episode_id].ledger_sha256
        == first_digest
    )
