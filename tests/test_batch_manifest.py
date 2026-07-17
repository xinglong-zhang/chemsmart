"""Tests for batch manifest helpers used by array submissions."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.cli.pka import rewrite_pka_batch_cli_args
from chemsmart.jobs.batch_manifest import (
    batch_manifest_filename,
    build_manifest_children,
    get_job_batch_entry,
    load_batch_manifest_entry,
    resolve_array_cli_args,
    set_job_batch_entry,
    write_batch_manifest,
)


def test_resolve_array_cli_args_homogeneous_keeps_shared_list():
    jobs = [SimpleNamespace(label="a"), SimpleNamespace(label="b")]
    shared = ["gaussian", "-f", "mols.xyz", "-i", "1,2", "opt"]
    assert resolve_array_cli_args(jobs, shared) == shared


def test_resolve_array_cli_args_heterogeneous_requires_rewrite_cli():
    job = SimpleNamespace(label="acid1")
    set_job_batch_entry(
        job,
        {
            "filepath": "acid1.xyz",
            "proton_index": 2,
            "charge": 0,
            "multiplicity": 1,
            "scheme": "direct",
            "label": "acid1",
        },
    )
    with pytest.raises(ValueError, match="rewrite_cli"):
        resolve_array_cli_args([job], ["gaussian", "pka", "batch"])


def test_resolve_array_cli_args_heterogeneous_uses_rewrite_callback():
    job1 = SimpleNamespace(label="acid1")
    job2 = SimpleNamespace(label="acid2")
    set_job_batch_entry(
        job1,
        {
            "filepath": "acid1.xyz",
            "proton_index": 2,
            "charge": 0,
            "multiplicity": 1,
            "scheme": "direct",
            "label": "acid1",
        },
    )
    set_job_batch_entry(
        job2,
        {
            "filepath": "acid2.xyz",
            "proton_index": 2,
            "charge": 1,
            "multiplicity": 2,
            "scheme": "direct",
            "label": "acid2",
        },
    )
    shared = ["gaussian", "-f", "table.csv", "pka", "-s", "direct", "batch"]
    cli_lists = resolve_array_cli_args(
        [job1, job2], shared, rewrite_cli=rewrite_pka_batch_cli_args
    )
    assert len(cli_lists) == 2
    assert "acid1.xyz" in cli_lists[0]
    assert "acid2.xyz" in cli_lists[1]
    assert all("submit" in args for args in cli_lists)


def test_write_and_load_batch_manifest(tmp_path):
    job1 = SimpleNamespace(label="acid1")
    set_job_batch_entry(
        job1,
        {
            "filepath": "acid1.xyz",
            "proton_index": 2,
            "charge": 0,
            "multiplicity": 1,
            "scheme": "direct",
            "label": "acid1",
        },
    )
    shared = ["gaussian", "-f", "table.csv", "pka", "batch"]
    children = build_manifest_children(
        [job1], shared, rewrite_cli=rewrite_pka_batch_cli_args
    )
    path = write_batch_manifest(
        batch_label="pka_scale_pka_batch",
        program="gaussian",
        children=children,
        directory=tmp_path,
    )
    assert path.name == batch_manifest_filename("pka_scale_pka_batch")
    entry = load_batch_manifest_entry(path, 1)
    assert entry["label"] == "acid1"
    assert entry["task_id"] == 1
    assert get_job_batch_entry(job1)["filepath"] == "acid1.xyz"
    assert "submit" in entry["cli_args"]
