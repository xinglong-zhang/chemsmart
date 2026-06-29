from __future__ import annotations

import json

import pytest

from chemsmart.agent.handles import HandleStore, is_handle_id


@pytest.mark.parametrize(
    ("kind", "prefix"),
    [
        ("mol", "mol_"),
        ("gset", "gset_"),
        ("oset", "oset_"),
        ("job", "job_"),
        ("dryrun", "dryrun_"),
        ("runtime", "runtime_"),
        ("runresult", "runresult_"),
        ("geom", "geom_"),
        ("submit", "submit_"),
        ("recmethod", "recmethod_"),
    ],
)
def test_handle_store_uses_expected_type_prefixes(tmp_path, kind, prefix):
    store = HandleStore(tmp_path)

    handle_id = store.put(kind, object(), {"kind": kind})

    assert handle_id.startswith(prefix)
    assert is_handle_id(handle_id) is True


def test_handle_store_mints_unique_ids_and_returns_same_object(tmp_path):
    store = HandleStore(tmp_path)

    handles = set()
    for index in range(25):
        obj = object()
        handle_id = store.put("job", obj, {"index": index})
        assert handle_id not in handles
        assert store.get(handle_id) is obj
        handles.add(handle_id)


def test_handle_store_round_trips_summaries_after_reload(tmp_path):
    store = HandleStore(tmp_path)
    summary = {
        "inputfile": str(tmp_path / "job.com"),
        "content": "#p opt\n",
        "meta": {"turn": 1},
    }

    handle_id = store.put("dryrun", object(), summary)
    reloaded = HandleStore(tmp_path)

    assert store.get_summary(handle_id) == summary
    assert reloaded.get_summary(handle_id) == summary


def test_handle_store_persists_jsonl_record_on_every_put(tmp_path):
    store = HandleStore(tmp_path)

    first = store.put("runtime", object(), {"ok": True})
    second = store.put("submit", object(), {"job_id": "12345"})

    lines = (
        (tmp_path / "handles.jsonl").read_text(encoding="utf-8").splitlines()
    )

    assert len(lines) == 2
    records = [json.loads(line) for line in lines]
    assert records[0]["handle_id"] == first
    assert records[0]["kind"] == "runtime"
    assert records[0]["summary"] == {"ok": True}
    assert records[1]["handle_id"] == second
    assert records[1]["kind"] == "submit"
    assert records[1]["summary"] == {"job_id": "12345"}
