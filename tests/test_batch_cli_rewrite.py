"""Tests for batch submit-time CLI rewrite helpers used by array submissions."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.cli.pka import rewrite_pka_batch_cli_args
from chemsmart.jobs.batch import (
    get_job_batch_entry,
    prepare_batch_jobs,
    prepare_nestable_batch_jobs,
    resolve_array_cli_args,
    rewrite_batch_cli_args,
    rewrite_nestable_cli_args,
    set_job_batch_entry,
)


def test_resolve_array_cli_args_without_batch_entry_keeps_shared_list():
    jobs = [SimpleNamespace(label="a"), SimpleNamespace(label="b")]
    shared = ["gaussian", "-f", "mols.xyz", "-i", "1,2", "opt"]
    assert resolve_array_cli_args(jobs, shared) == shared


def test_resolve_array_cli_args_with_batch_entry_uses_shared_rewrite():
    job1 = SimpleNamespace(label="mol_idx1")
    job2 = SimpleNamespace(label="mol_idx2")
    set_job_batch_entry(job1, {"filepath": "mols.xyz", "molecule_index": 1})
    set_job_batch_entry(job2, {"filepath": "mols.xyz", "molecule_index": 2})
    shared = ["gaussian", "-f", "mols.xyz", "-i", "1,2", "opt"]
    cli_lists = resolve_array_cli_args(
        [job1, job2],
        shared,
        rewrite_cli=rewrite_batch_cli_args,
    )
    assert cli_lists[0][cli_lists[0].index("-i") + 1] == "1"
    assert cli_lists[1][cli_lists[1].index("-i") + 1] == "2"
    assert all("1,2" not in args for args in cli_lists)


def test_rewrite_batch_cli_args_narrows_shared_index():
    shared = ["gaussian", "-f", "mols.xyz", "-i", "1,2,3", "opt"]
    rewritten = rewrite_batch_cli_args(
        shared,
        {"filepath": "mols.xyz", "molecule_index": 2},
    )
    assert rewritten[rewritten.index("-i") + 1] == "2"
    assert "1,2,3" not in rewritten


def test_rewrite_pka_batch_cli_args_builds_on_shared_rewrite():
    entry = {
        "filepath": "acid.xyz",
        "proton_index": 2,
        "charge": 0,
        "multiplicity": 1,
        "scheme": "direct",
        "label": "acid",
    }
    shared = ["gaussian", "-f", "table.csv", "pka", "batch"]
    rewritten = rewrite_pka_batch_cli_args(shared, entry)
    assert "submit" in rewritten
    assert "acid.xyz" in rewritten
    assert rewritten[rewritten.index("--proton-index") + 1] == "2"


def test_prepare_batch_jobs_attaches_entries():
    jobs = [SimpleNamespace(label="a"), SimpleNamespace(label="b")]
    rewrite_cli = prepare_batch_jobs(jobs, [1, 2], filepath="mols.xyz")
    assert rewrite_cli is rewrite_batch_cli_args
    assert get_job_batch_entry(jobs[0]) == {
        "filepath": "mols.xyz",
        "molecule_index": 1,
    }
    assert get_job_batch_entry(jobs[1])["molecule_index"] == 2


def test_prepare_batch_jobs_returns_none_for_single_job():
    jobs = [SimpleNamespace(label="a")]
    assert prepare_batch_jobs(jobs, [1], filepath="mols.xyz") is None
    assert get_job_batch_entry(jobs[0]) is None


def test_prepare_batch_jobs_rejects_length_mismatch():
    jobs = [
        SimpleNamespace(label="a"),
        SimpleNamespace(label="b"),
        SimpleNamespace(label="c"),
    ]
    with pytest.raises(ValueError, match="3 job\\(s\\) vs 2 molecule index"):
        prepare_batch_jobs(jobs, [1, 2], filepath="mols.xyz")
    assert get_job_batch_entry(jobs[0]) is None
    assert get_job_batch_entry(jobs[1]) is None
    assert get_job_batch_entry(jobs[2]) is None


def test_resolve_array_cli_args_requires_rewrite_cli_when_entries_present():
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


def test_resolve_array_cli_args_pka_entries_use_rewrite_callback():
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


def test_prepare_nestable_batch_jobs_injects_child_index():
    jobs = [SimpleNamespace(label="qf"), SimpleNamespace(label="qr")]
    rewrite_cli = prepare_nestable_batch_jobs(jobs)
    assert rewrite_cli is rewrite_nestable_cli_args
    assert get_job_batch_entry(jobs[0])["child_index"] == 1
    assert get_job_batch_entry(jobs[1])["child_index"] == 2

    shared = ["gaussian", "-f", "ts.log", "qrc"]
    cli_lists = resolve_array_cli_args(
        jobs, shared, rewrite_cli=rewrite_nestable_cli_args
    )
    assert cli_lists[0][cli_lists[0].index("--child-index") + 1] == "1"
    assert cli_lists[1][cli_lists[1].index("--child-index") + 1] == "2"
    assert all("qrc" in args for args in cli_lists)


def test_rewrite_nestable_cli_args_without_entry_keeps_shared():
    shared = ["gaussian", "-f", "ts.log", "qrc"]
    assert rewrite_nestable_cli_args(shared, None) == shared
