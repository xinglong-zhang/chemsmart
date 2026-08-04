from __future__ import annotations

from dataclasses import asdict
import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import canonical_sha256
from chemsmart.agent.experiments import implementation_freeze as freeze
from chemsmart.agent.experiments.implementation_freeze import (
    GitWorktreeObservationV1,
    ImplementationFreezeError,
    ImplementationFreezeIntegrityError,
    create_implementation_freeze,
    verify_implementation_freeze,
)


def _git_observation() -> GitWorktreeObservationV1:
    status = b" M chemsmart/core.py\x00?? tools/run_campaign.py\x00"
    diff = b"diff --git a/chemsmart/core.py b/chemsmart/core.py\n"
    body = {
        "schema_version": "chemsmart.git-worktree-observation.v1",
        "head": "1" * 40,
        "branch": "integration/example",
        "status_nonempty": True,
        "status_entry_count": 2,
        "status_sha256": hashlib.sha256(status).hexdigest(),
        "dirty_diff_sha256": hashlib.sha256(diff).hexdigest(),
    }
    return GitWorktreeObservationV1(
        **body, observation_sha256=canonical_sha256(body)
    )


def _workspace(root: Path) -> Path:
    workspace = root / "workspace"
    files = {
        "chemsmart/__init__.py": "__version__ = '3.1.4'\n",
        "chemsmart/core.py": "VALUE = 1\n",
        "tests/agent/test_core.py": "def test_core():\n    assert True\n",
        "tests/data/case.json": '{"case":"water"}\n',
        "docs/evaluation/protocol.md": "# Protocol\n",
        "pyproject.toml": '[project]\nname = "chemsmart"\n',
        "MANIFEST.in": "recursive-include chemsmart *\n",
    }
    for relative, content in files.items():
        path = workspace / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755 if relative == "tools/run_campaign.py" else 0o644)
    (workspace / "chemsmart" / "__pycache__").mkdir()
    (workspace / "chemsmart" / "__pycache__" / "core.pyc").write_bytes(
        b"not-public-source"
    )
    return workspace


def _campaign_root(root: Path) -> Path:
    campaign = root / "campaign"
    files = {
        "run_campaign.py": "print('public launcher')\n",
        "config/qwen-public.yaml": "model: qwen3.8-max\n",
    }
    for relative, content in files.items():
        path = campaign / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        path.chmod(0o755 if relative == "run_campaign.py" else 0o644)
    return campaign


def _create(monkeypatch, workspace: Path, evidence: Path):
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())
    campaign = _campaign_root(evidence.parent / "campaign-source")
    return create_implementation_freeze(
        workspace_root=workspace,
        evidence_directory=evidence,
        relevant_paths=(
            "tests/agent/test_core.py",
            "tests/data",
            "docs/evaluation/protocol.md",
        ),
        campaign_root=campaign,
        campaign_input_allowlist=(
            "run_campaign.py",
            "config/qwen-public.yaml",
        ),
    )


def test_freeze_copies_and_restores_exact_public_source(tmp_path, monkeypatch):
    workspace = _workspace(tmp_path / "source")
    evidence = tmp_path / "evidence"

    manifest = _create(monkeypatch, workspace, evidence)

    paths = tuple((item.namespace, item.relative_path) for item in manifest.files)
    assert paths == tuple(sorted(paths))
    assert ("worktree", "chemsmart/__pycache__/core.pyc") not in paths
    assert {
        ("worktree", "chemsmart/__init__.py"),
        ("worktree", "chemsmart/core.py"),
        ("worktree", "tests/agent/test_core.py"),
        ("worktree", "tests/data/case.json"),
        ("worktree", "docs/evaluation/protocol.md"),
        ("worktree", "pyproject.toml"),
        ("worktree", "MANIFEST.in"),
        ("campaign-inputs", "run_campaign.py"),
        ("campaign-inputs", "config/qwen-public.yaml"),
    } == set(paths)
    assert all(not Path(value).is_absolute() for _namespace, value in paths)
    launcher = next(
        item
        for item in manifest.files
        if item.namespace == "campaign-inputs"
        and item.relative_path == "run_campaign.py"
    )
    assert launcher.mode == "0755"
    assert manifest.git.status_nonempty is True
    assert manifest.restoration_verified is True
    manifest_path = evidence / "implementation-freeze-manifest.json"
    assert verify_implementation_freeze(manifest_path) == manifest
    assert (
        evidence
        / "implementation-snapshot"
        / "worktree"
        / "chemsmart"
        / "core.py"
    ).read_bytes() == (workspace / "chemsmart" / "core.py").read_bytes()
    assert (
        evidence
        / "implementation-snapshot"
        / "campaign-inputs"
        / "run_campaign.py"
    ).read_text(encoding="utf-8") == "print('public launcher')\n"


def test_freeze_rejects_source_mutation_during_copy(tmp_path, monkeypatch):
    workspace = _workspace(tmp_path / "source")
    evidence = tmp_path / "evidence"
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())
    original = freeze._copy_file_bytes
    mutated = False

    def mutate_after_copy(source_descriptor, destination):
        nonlocal mutated
        result = original(source_descriptor, destination)
        if destination.as_posix().endswith("chemsmart/core.py") and not mutated:
            mutated = True
            (workspace / "chemsmart" / "core.py").write_text(
                "VALUE = 2\n", encoding="utf-8"
            )
        return result

    monkeypatch.setattr(freeze, "_copy_file_bytes", mutate_after_copy)

    with pytest.raises(
        ImplementationFreezeIntegrityError, match="changed during copy"
    ):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=evidence,
        )
    assert not (evidence / "implementation-freeze-manifest.json").exists()


def test_verification_rejects_snapshot_hash_mismatch(tmp_path, monkeypatch):
    workspace = _workspace(tmp_path / "source")
    evidence = tmp_path / "evidence"
    _create(monkeypatch, workspace, evidence)
    copied = (
        evidence
        / "implementation-snapshot"
        / "worktree"
        / "chemsmart"
        / "core.py"
    )
    copied.write_text("VALUE = 999\n", encoding="utf-8")

    with pytest.raises(ImplementationFreezeIntegrityError):
        verify_implementation_freeze(
            evidence / "implementation-freeze-manifest.json"
        )


def test_freeze_rejects_secret_symlink_and_path_escape(tmp_path, monkeypatch):
    workspace = _workspace(tmp_path / "source")
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())
    campaign = _campaign_root(tmp_path / "campaign-source")
    secret = campaign / "config" / "api.env"
    secret.write_text("DO_NOT_COPY=value\n", encoding="utf-8")

    with pytest.raises(ImplementationFreezeError, match="secret"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "secret-evidence",
            campaign_root=campaign,
            campaign_input_allowlist=("config/api.env",),
        )
    with pytest.raises(ImplementationFreezeError, match="canonical|escapes"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "escape-evidence",
            relevant_paths=("../outside.py",),
        )

    external = tmp_path / "external.py"
    external.write_text("EXTERNAL = True\n", encoding="utf-8")
    (workspace / "chemsmart" / "linked.py").symlink_to(external)
    with pytest.raises(ImplementationFreezeError, match="symlink"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "symlink-evidence",
        )


def test_campaign_namespace_rejects_unapproved_or_symlink_sources(
    tmp_path, monkeypatch
):
    workspace = _workspace(tmp_path / "source")
    campaign = _campaign_root(tmp_path / "campaign-source")
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())

    with pytest.raises(ImplementationFreezeError, match="absolute"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "relative-root-evidence",
            campaign_root=Path("relative-campaign"),
            campaign_input_allowlist=("run_campaign.py",),
        )

    linked_root = tmp_path / "linked-campaign"
    linked_root.symlink_to(campaign, target_is_directory=True)
    with pytest.raises(ImplementationFreezeError, match="non-symlink"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "linked-root-evidence",
            campaign_root=linked_root,
            campaign_input_allowlist=("run_campaign.py",),
        )

    external = tmp_path / "external-config.yaml"
    external.write_text("model: forbidden-link\n", encoding="utf-8")
    (campaign / "config" / "linked.yaml").symlink_to(external)
    with pytest.raises(ImplementationFreezeError, match="symlink"):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=tmp_path / "linked-file-evidence",
            campaign_root=campaign,
            campaign_input_allowlist=("config/linked.yaml",),
        )


def test_freeze_rejects_campaign_input_mutation_during_copy(
    tmp_path, monkeypatch
):
    workspace = _workspace(tmp_path / "source")
    campaign = _campaign_root(tmp_path / "campaign-source")
    evidence = tmp_path / "campaign-mutation-evidence"
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())
    original = freeze._copy_file_bytes
    mutated = False

    def mutate_campaign_after_copy(source_descriptor, destination):
        nonlocal mutated
        result = original(source_descriptor, destination)
        if destination.as_posix().endswith(
            "campaign-inputs/run_campaign.py"
        ) and not mutated:
            mutated = True
            (campaign / "run_campaign.py").write_text(
                "print('changed')\n", encoding="utf-8"
            )
        return result

    monkeypatch.setattr(
        freeze, "_copy_file_bytes", mutate_campaign_after_copy
    )

    with pytest.raises(
        ImplementationFreezeIntegrityError, match="changed during copy"
    ):
        create_implementation_freeze(
            workspace_root=workspace,
            evidence_directory=evidence,
            campaign_root=campaign,
            campaign_input_allowlist=("run_campaign.py",),
        )


def test_manifest_identity_is_independent_of_absolute_paths(tmp_path, monkeypatch):
    first_workspace = _workspace(tmp_path / "first")
    second_workspace = _workspace(tmp_path / "elsewhere" / "second")
    monkeypatch.setattr(freeze, "_observe_git", lambda _root: _git_observation())
    first_campaign = _campaign_root(tmp_path / "first-campaign-source")
    second_campaign = _campaign_root(tmp_path / "second-campaign-source")

    first = create_implementation_freeze(
        workspace_root=first_workspace,
        evidence_directory=tmp_path / "first-evidence",
        relevant_paths=("tests", "docs"),
        campaign_root=first_campaign,
        campaign_input_allowlist=(
            "config/qwen-public.yaml",
            "run_campaign.py",
        ),
    )
    second = create_implementation_freeze(
        workspace_root=second_workspace,
        evidence_directory=tmp_path / "second-evidence",
        relevant_paths=("tests", "docs"),
        campaign_root=second_campaign,
        campaign_input_allowlist=(
            "config/qwen-public.yaml",
            "run_campaign.py",
        ),
    )

    assert first.tree_sha256 == second.tree_sha256
    assert first.manifest_sha256 == second.manifest_sha256
    assert asdict(first) == asdict(second)
