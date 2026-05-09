"""File indexing helpers for composer popups."""

from __future__ import annotations

from concurrent.futures import Future
from pathlib import Path
from threading import Lock, Thread

_BIASED_EXTENSIONS = [".xyz", ".log", ".com", ".inp", ".gjf", ".out"]
_FILE_INDEX_LOCK = Lock()
_FILE_INDEX_CACHE: dict[tuple[str, int], list[Path]] = {}
_FILE_INDEX_FUTURES: dict[tuple[str, int], Future] = {}


def iter_candidate_files(root: str | Path, *, depth: int = 3) -> list[Path]:
    root_path = Path(root)
    if not root_path.exists():
        return []

    request_candidate_file_refresh(root_path, depth=depth)
    key = _file_index_key(root_path, depth)
    with _FILE_INDEX_LOCK:
        return list(_FILE_INDEX_CACHE.get(key, []))


def request_candidate_file_refresh(
    root: str | Path,
    *,
    depth: int = 3,
) -> Future | None:
    root_path = Path(root)
    if not root_path.exists():
        return None

    key = _file_index_key(root_path, depth)
    with _FILE_INDEX_LOCK:
        future = _FILE_INDEX_FUTURES.get(key)
        if future is not None and not future.done():
            return future
        future = Future()
        _FILE_INDEX_FUTURES[key] = future
        future.add_done_callback(
            lambda done, cache_key=key: _clear_candidate_file_future(
                cache_key, done
            )
        )
        Thread(
            target=_run_candidate_file_refresh,
            args=(future, root_path, depth),
            daemon=True,
            name="chemsmart-file-index",
        ).start()
        return future


def refresh_candidate_file_cache(
    root: str | Path,
    *,
    depth: int = 3,
) -> list[Path]:
    root_path = Path(root)
    candidates = _scan_candidate_files(root_path, depth=depth)
    with _FILE_INDEX_LOCK:
        _FILE_INDEX_CACHE[_file_index_key(root_path, depth)] = list(candidates)
    return candidates


def _scan_candidate_files(root_path: Path, *, depth: int) -> list[Path]:
    candidates: list[Path] = []
    for path in root_path.rglob("*"):
        if not path.is_file():
            continue
        try:
            relative = path.relative_to(root_path)
        except ValueError:
            continue
        if len(relative.parts) > depth:
            continue
        candidates.append(path)

    def sort_key(path: Path):
        suffix = path.suffix.lower()
        bias = (
            _BIASED_EXTENSIONS.index(suffix)
            if suffix in _BIASED_EXTENSIONS
            else len(_BIASED_EXTENSIONS)
        )
        return (bias, len(path.parts), str(path).lower())

    return sorted(candidates, key=sort_key)


def _file_index_key(root: Path, depth: int) -> tuple[str, int]:
    return (str(root.resolve()), depth)


def _clear_candidate_file_future(
    cache_key: tuple[str, int],
    future: Future,
) -> None:
    with _FILE_INDEX_LOCK:
        existing = _FILE_INDEX_FUTURES.get(cache_key)
        if existing is future:
            _FILE_INDEX_FUTURES.pop(cache_key, None)


def _run_candidate_file_refresh(
    future: Future,
    root_path: Path,
    depth: int,
) -> None:
    if not future.set_running_or_notify_cancel():
        return
    try:
        result = refresh_candidate_file_cache(root_path, depth=depth)
    except Exception as exc:
        future.set_exception(exc)
    else:
        future.set_result(result)
