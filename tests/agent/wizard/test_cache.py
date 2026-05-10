from __future__ import annotations

from datetime import UTC, datetime, timedelta

from chemsmart.agent.wizard.cache import (
    CacheEntry,
    cache_path,
    is_stale,
    load_cache,
    mark_status,
    write_cache,
)


def _entry(
    *, status: str = "fresh", probed_at: str | None = None
) -> CacheEntry:
    return CacheEntry(
        server_name="perlmutter",
        host="login.cluster",
        mode="ssh",
        scheduler="SLURM",
        probed_at=probed_at
        or datetime.now(UTC).isoformat().replace("+00:00", "Z"),
        source_commands={"scheduler": "sinfo --json"},
        partitions=[{"name": "debug", "default": True}],
        node_summary={"cpu": 16, "mem_gb": 64, "gpu": 0},
        program_candidates={"gaussian": {"source": "path"}},
        status=status,
        last_error=None,
    )


def test_write_and_load_cache_round_trip(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    entry = _entry()

    written = write_cache(entry)
    loaded = load_cache("perlmutter")

    assert written == str(cache_path("perlmutter"))
    assert loaded == entry


def test_load_cache_returns_none_for_missing_or_corrupt(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))

    assert load_cache("missing") is None

    path = cache_path("broken")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("{not-json", encoding="utf-8")

    assert load_cache("broken") is None


def test_is_stale_true_when_entry_exceeds_ttl(monkeypatch, tmp_path):
    monkeypatch.setenv("HOME", str(tmp_path))
    old = datetime.now(UTC) - timedelta(hours=25)
    entry = _entry(probed_at=old.isoformat().replace("+00:00", "Z"))

    assert is_stale(entry, ttl_hours=24) is True


def test_mark_status_returns_new_instance():
    entry = _entry()

    updated = mark_status(entry, "stale", last_error="boom")

    assert updated is not entry
    assert updated.status == "stale"
    assert updated.last_error == "boom"
    assert entry.status == "fresh"
    assert entry.last_error is None
