"""Immutable implementation snapshots for causal agent experiments.

An experiment configuration is not reproducible when it records only prompt
and tool-schema digests while the implementation under test remains mutable.
This module copies the exact public implementation bytes into a local evidence
directory *before* outcomes are observed.  It deliberately records Git state as
an observation, including a dirty worktree, rather than treating Git cleanliness
as an admission requirement.

The snapshot never executes project code. Source paths are relative to either
the worktree or an explicitly approved campaign-input root and are always
regular, non-symlink files. The copy operation checks each source before and
after streaming, rescans and rehashes the complete source set, and verifies the
copied tree by rereading it before publishing the manifest.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
from typing import Any, Iterable, Mapping, Sequence

from chemsmart.agent._contracts import (
    ContractError,
    canonical_json,
    canonical_sha256,
    require_sha256,
)


_SCHEMA_VERSION = "chemsmart.implementation-freeze-manifest.v1"
_FILE_SCHEMA_VERSION = "chemsmart.implementation-freeze-file.v1"
_GIT_SCHEMA_VERSION = "chemsmart.git-worktree-observation.v1"
_SNAPSHOT_DIRECTORY = "implementation-snapshot"
_MANIFEST_NAME = "implementation-freeze-manifest.json"
_WORKTREE_NAMESPACE = "worktree"
_CAMPAIGN_NAMESPACE = "campaign-inputs"
_DEFAULT_SOURCE_ROOTS = ("chemsmart",)
_PACKAGING_METADATA = (
    "pyproject.toml",
    "setup.py",
    "setup.cfg",
    "MANIFEST.in",
)
_EXCLUDED_DIRECTORY_NAMES = frozenset(
    {
        ".git",
        ".codex",
        ".chemsmart-agent",
        ".eggs",
        ".mypy_cache",
        ".nox",
        ".pytest_cache",
        ".ruff_cache",
        ".tox",
        ".venv",
        "__pycache__",
        "build",
        "dist",
        "env",
        "environments",
        "node_modules",
        "provider-private",
        "provider-raw",
        "provider-transcripts",
        "raw-provider",
        "raw_provider",
        "runs",
        "site-packages",
        "venv",
    }
)
_SECRET_FILE_NAMES = frozenset(
    {
        ".env",
        ".netrc",
        "api.env",
        "credentials.json",
        "secrets.json",
    }
)
_SECRET_SUFFIXES = frozenset({".key", ".p12", ".pem", ".pfx"})
_RAW_PROVIDER_FILE_RE = re.compile(
    r"(?:^|[-_.])(?:private[-_]?reasoning|provider[-_]?raw|raw[-_]?provider|"
    r"raw[-_]?(?:request|response)|reasoning[-_]?content|request[-_]?headers)"
    r"(?:[-_.]|$)",
    flags=re.IGNORECASE,
)
_GIT_OBJECT_RE = re.compile(r"^[0-9a-f]{40}(?:[0-9a-f]{24})?$")
_COPY_CHUNK_BYTES = 1024 * 1024


class ImplementationFreezeError(ContractError):
    """The implementation could not be frozen without ambiguity."""


class ImplementationFreezeIntegrityError(ImplementationFreezeError):
    """A source or copied snapshot changed during the freeze operation."""


def _require_relative_workspace_path(value: str | Path, field: str) -> str:
    raw = str(value)
    if not raw or "\x00" in raw or "\\" in raw:
        raise ImplementationFreezeError(
            f"{field} must be a canonical workspace-relative POSIX path"
        )
    candidate = PurePosixPath(raw)
    if candidate.is_absolute() or raw != candidate.as_posix():
        raise ImplementationFreezeError(
            f"{field} must be a canonical workspace-relative POSIX path"
        )
    if any(part in {"", ".", ".."} for part in candidate.parts):
        raise ImplementationFreezeError(f"{field} escapes the workspace")
    return candidate.as_posix()


def _is_secret_or_raw_provider_file(path: PurePosixPath) -> bool:
    name = path.name.lower()
    return (
        name in _SECRET_FILE_NAMES
        or name.startswith(".env.")
        or path.suffix.lower() in _SECRET_SUFFIXES
        or _RAW_PROVIDER_FILE_RE.search(name) is not None
    )


def _path_has_excluded_directory(path: PurePosixPath) -> bool:
    return any(part.lower() in _EXCLUDED_DIRECTORY_NAMES for part in path.parts)


def _mode_text(mode: int) -> str:
    return f"{stat.S_IMODE(mode):04o}"


def _parse_mode(value: str) -> int:
    text = str(value)
    if re.fullmatch(r"[0-7]{4}", text) is None:
        raise ImplementationFreezeError("snapshot mode must be four octal digits")
    return int(text, 8)


@dataclass(frozen=True)
class ImplementationFreezeFileV1:
    """One exact namespace-relative source file in the frozen snapshot."""

    schema_version: str
    namespace: str
    relative_path: str
    size_bytes: int
    mode: str
    sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _FILE_SCHEMA_VERSION:
            raise ImplementationFreezeError(
                "unsupported implementation-freeze file record"
            )
        if self.namespace not in {
            _WORKTREE_NAMESPACE,
            _CAMPAIGN_NAMESPACE,
        }:
            raise ImplementationFreezeError("unsupported snapshot namespace")
        _require_relative_workspace_path(self.relative_path, "relative_path")
        if self.size_bytes < 0:
            raise ImplementationFreezeError("snapshot size must be non-negative")
        _parse_mode(self.mode)
        require_sha256(self.sha256, "sha256")


@dataclass(frozen=True)
class GitWorktreeObservationV1:
    """Read-only Git observations; this contract never asserts cleanliness."""

    schema_version: str
    head: str
    branch: str
    status_nonempty: bool
    status_entry_count: int
    status_sha256: str
    dirty_diff_sha256: str
    observation_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _GIT_SCHEMA_VERSION:
            raise ImplementationFreezeError("unsupported Git observation")
        if _GIT_OBJECT_RE.fullmatch(str(self.head)) is None:
            raise ImplementationFreezeError("Git HEAD is not a supported object ID")
        branch = str(self.branch)
        if (
            not branch
            or len(branch) > 255
            or any(ord(char) < 32 for char in branch)
        ):
            raise ImplementationFreezeError("Git branch observation is invalid")
        if self.status_entry_count < 0:
            raise ImplementationFreezeError(
                "Git status entry count must be non-negative"
            )
        if self.status_nonempty != (self.status_entry_count > 0):
            raise ImplementationFreezeError("Git status observation is inconsistent")
        require_sha256(self.status_sha256, "status_sha256")
        require_sha256(self.dirty_diff_sha256, "dirty_diff_sha256")
        require_sha256(self.observation_sha256, "observation_sha256")
        body = {
            key: value
            for key, value in asdict(self).items()
            if key != "observation_sha256"
        }
        if self.observation_sha256 != canonical_sha256(body):
            raise ImplementationFreezeError("Git observation digest mismatch")


@dataclass(frozen=True)
class ImplementationFreezeManifestV1:
    """Path-independent manifest for one verified implementation snapshot."""

    schema_version: str
    snapshot_directory: str
    source_roots: tuple[str, ...]
    relevant_paths: tuple[str, ...]
    campaign_input_allowlist: tuple[str, ...]
    files: tuple[ImplementationFreezeFileV1, ...]
    file_count: int
    total_size_bytes: int
    tree_sha256: str
    restored_tree_sha256: str
    restoration_verified: bool
    git: GitWorktreeObservationV1
    manifest_sha256: str

    def __post_init__(self) -> None:
        if self.schema_version != _SCHEMA_VERSION:
            raise ImplementationFreezeError(
                "unsupported implementation-freeze manifest"
            )
        _require_relative_workspace_path(
            self.snapshot_directory, "snapshot_directory"
        )
        for field, values in (
            ("source_roots", self.source_roots),
            ("relevant_paths", self.relevant_paths),
            ("campaign_input_allowlist", self.campaign_input_allowlist),
        ):
            canonical = tuple(
                _require_relative_workspace_path(value, field) for value in values
            )
            if canonical != tuple(sorted(set(canonical))):
                raise ImplementationFreezeError(f"{field} must be sorted and unique")
        paths = tuple(
            (record.namespace, record.relative_path) for record in self.files
        )
        if paths != tuple(sorted(set(paths))):
            raise ImplementationFreezeError(
                "snapshot file records must be sorted and unique"
            )
        if self.file_count != len(self.files):
            raise ImplementationFreezeError("snapshot file count mismatch")
        if self.total_size_bytes != sum(item.size_bytes for item in self.files):
            raise ImplementationFreezeError("snapshot byte count mismatch")
        expected_tree = canonical_sha256(
            tuple(asdict(record) for record in self.files)
        )
        if self.tree_sha256 != expected_tree:
            raise ImplementationFreezeError("snapshot tree digest mismatch")
        if not self.restoration_verified:
            raise ImplementationFreezeError("snapshot restoration is not verified")
        if self.restored_tree_sha256 != self.tree_sha256:
            raise ImplementationFreezeError("restored snapshot tree digest mismatch")
        body = {
            key: value
            for key, value in asdict(self).items()
            if key != "manifest_sha256"
        }
        if self.manifest_sha256 != canonical_sha256(body):
            raise ImplementationFreezeError("implementation manifest digest mismatch")

    def public_record(self) -> dict[str, Any]:
        return asdict(self)


def _git_bytes(workspace: Path, arguments: Sequence[str]) -> bytes:
    try:
        completed = subprocess.run(
            ("git", "-C", os.fspath(workspace), *arguments),
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env={**os.environ, "LC_ALL": "C", "LANG": "C"},
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ImplementationFreezeError("Git observation failed") from exc
    return completed.stdout


def _observe_git(workspace: Path) -> GitWorktreeObservationV1:
    head = _git_bytes(workspace, ("rev-parse", "HEAD")).decode("ascii").strip()
    branch = _git_bytes(
        workspace, ("branch", "--show-current")
    ).decode("utf-8", errors="strict").strip()
    if not branch:
        branch = "(detached)"
    status = _git_bytes(
        workspace,
        (
            "status",
            "--porcelain=v1",
            "-z",
            "--untracked-files=all",
            "--ignored=no",
        ),
    )
    dirty_diff = _git_bytes(
        workspace,
        (
            "diff",
            "--binary",
            "--full-index",
            "--no-ext-diff",
            "--no-textconv",
            "--no-color",
            "HEAD",
            "--",
            ".",
        ),
    )
    status_count = len(tuple(item for item in status.split(b"\x00") if item))
    body = {
        "schema_version": _GIT_SCHEMA_VERSION,
        "head": head,
        "branch": branch,
        "status_nonempty": bool(status),
        "status_entry_count": status_count,
        "status_sha256": hashlib.sha256(status).hexdigest(),
        "dirty_diff_sha256": hashlib.sha256(dirty_diff).hexdigest(),
    }
    return GitWorktreeObservationV1(
        **body, observation_sha256=canonical_sha256(body)
    )


def _resolve_workspace_entry(workspace: Path, relative: str) -> Path:
    normalized = _require_relative_workspace_path(relative, "source path")
    lexical = workspace.joinpath(*PurePosixPath(normalized).parts)
    cursor = workspace
    for part in PurePosixPath(normalized).parts:
        cursor = cursor / part
        try:
            component_stat = os.lstat(cursor)
        except OSError as exc:
            raise ImplementationFreezeError(
                f"implementation source does not exist: {normalized}"
            ) from exc
        if stat.S_ISLNK(component_stat.st_mode):
            raise ImplementationFreezeError(
                f"implementation source uses a symlink: {normalized}"
            )
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as exc:
        raise ImplementationFreezeError(
            f"implementation source does not exist: {normalized}"
        ) from exc
    if not resolved.is_relative_to(workspace):
        raise ImplementationFreezeError(
            f"implementation source escapes the workspace: {normalized}"
        )
    if lexical.is_symlink() or resolved != lexical:
        raise ImplementationFreezeError(
            f"implementation source uses a symlink: {normalized}"
        )
    return lexical


def _collect_from_path(
    *, workspace: Path, relative: str, explicit: bool
) -> tuple[str, ...]:
    normalized = _require_relative_workspace_path(relative, "source path")
    pure = PurePosixPath(normalized)
    if _path_has_excluded_directory(pure):
        if explicit:
            raise ImplementationFreezeError(
                f"explicit source path is excluded: {normalized}"
            )
        return ()
    if _is_secret_or_raw_provider_file(pure):
        raise ImplementationFreezeError(
            f"secret or raw provider source is forbidden: {normalized}"
        )
    source = _resolve_workspace_entry(workspace, normalized)
    source_stat = os.lstat(source)
    if stat.S_ISLNK(source_stat.st_mode):
        raise ImplementationFreezeError(
            f"implementation source uses a symlink: {normalized}"
        )
    if stat.S_ISREG(source_stat.st_mode):
        return (normalized,)
    if not stat.S_ISDIR(source_stat.st_mode):
        raise ImplementationFreezeError(
            f"implementation source is not a regular file or directory: {normalized}"
        )

    collected: list[str] = []

    def walk(directory: Path, directory_relative: PurePosixPath) -> None:
        try:
            entries = sorted(os.scandir(directory), key=lambda item: item.name)
        except OSError as exc:
            raise ImplementationFreezeError(
                f"cannot inspect implementation source: {directory_relative}"
            ) from exc
        for entry in entries:
            child_relative = directory_relative / entry.name
            if entry.is_symlink():
                raise ImplementationFreezeError(
                    f"implementation source uses a symlink: {child_relative}"
                )
            if entry.is_dir(follow_symlinks=False):
                if entry.name.lower() in _EXCLUDED_DIRECTORY_NAMES:
                    continue
                walk(Path(entry.path), child_relative)
                continue
            if not entry.is_file(follow_symlinks=False):
                raise ImplementationFreezeError(
                    f"implementation source is not regular: {child_relative}"
                )
            if _is_secret_or_raw_provider_file(child_relative):
                raise ImplementationFreezeError(
                    "secret or raw provider source is forbidden: "
                    f"{child_relative}"
                )
            collected.append(child_relative.as_posix())

    walk(source, pure)
    return tuple(collected)


def _collect_source_paths(
    *,
    workspace: Path,
    relevant_paths: Sequence[str],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    roots: list[tuple[str, bool]] = [
        (value, False) for value in _DEFAULT_SOURCE_ROOTS
    ]
    roots.extend(
        (value, False)
        for value in _PACKAGING_METADATA
        if workspace.joinpath(value).exists()
    )
    roots.extend((value, True) for value in relevant_paths)

    collected: set[str] = set()
    for relative, explicit in roots:
        collected.update(
            _collect_from_path(
                workspace=workspace, relative=relative, explicit=explicit
            )
        )
    if not collected:
        raise ImplementationFreezeError("implementation snapshot is empty")
    return tuple(sorted(collected)), tuple(
        sorted(
            value
            for value in _PACKAGING_METADATA
            if workspace.joinpath(value).exists()
        )
    )


def _collect_campaign_inputs(
    *, campaign_root: Path, allowlist: Sequence[str]
) -> tuple[str, ...]:
    collected: list[str] = []
    for relative in allowlist:
        normalized = _require_relative_workspace_path(
            relative, "campaign_input_allowlist"
        )
        pure = PurePosixPath(normalized)
        if _path_has_excluded_directory(pure):
            raise ImplementationFreezeError(
                f"campaign input is excluded: {normalized}"
            )
        if _is_secret_or_raw_provider_file(pure):
            raise ImplementationFreezeError(
                f"secret or raw provider campaign input is forbidden: {normalized}"
            )
        source = _resolve_workspace_entry(campaign_root, normalized)
        source_stat = os.lstat(source)
        if not stat.S_ISREG(source_stat.st_mode) or stat.S_ISLNK(
            source_stat.st_mode
        ):
            raise ImplementationFreezeError(
                f"campaign input is not a regular file: {normalized}"
            )
        collected.append(normalized)
    return tuple(sorted(set(collected)))


def _stat_signature(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _write_all(descriptor: int, data: bytes) -> None:
    view = memoryview(data)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise OSError("snapshot write made no progress")
        view = view[written:]


def _copy_file_bytes(source_descriptor: int, destination: Path) -> tuple[int, str]:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    destination_descriptor = os.open(destination, flags, 0o600)
    digest = hashlib.sha256()
    size = 0
    try:
        while True:
            chunk = os.read(source_descriptor, _COPY_CHUNK_BYTES)
            if not chunk:
                break
            digest.update(chunk)
            size += len(chunk)
            _write_all(destination_descriptor, chunk)
        os.fsync(destination_descriptor)
    finally:
        os.close(destination_descriptor)
    return size, digest.hexdigest()


def _hash_stable_source(source: Path) -> tuple[int, str, os.stat_result]:
    before = os.lstat(source)
    if not stat.S_ISREG(before.st_mode) or stat.S_ISLNK(before.st_mode):
        raise ImplementationFreezeIntegrityError(
            "implementation source is no longer a regular file"
        )
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(source, flags)
    digest = hashlib.sha256()
    size = 0
    try:
        opened = os.fstat(descriptor)
        if _stat_signature(before) != _stat_signature(opened):
            raise ImplementationFreezeIntegrityError(
                "implementation source changed before hashing"
            )
        while True:
            chunk = os.read(descriptor, _COPY_CHUNK_BYTES)
            if not chunk:
                break
            digest.update(chunk)
            size += len(chunk)
        after_open = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = os.lstat(source)
    signature = _stat_signature(before)
    if (
        signature != _stat_signature(after_open)
        or signature != _stat_signature(after_path)
        or size != before.st_size
    ):
        raise ImplementationFreezeIntegrityError(
            "implementation source changed while hashing"
        )
    return size, digest.hexdigest(), before


def _copy_stable_source(
    *,
    source_root: Path,
    namespace: str,
    relative: str,
    snapshot_root: Path,
) -> ImplementationFreezeFileV1:
    source = _resolve_workspace_entry(source_root, relative)
    before = os.lstat(source)
    if not stat.S_ISREG(before.st_mode) or stat.S_ISLNK(before.st_mode):
        raise ImplementationFreezeIntegrityError(
            f"implementation source is not stable: {relative}"
        )
    destination = snapshot_root.joinpath(
        namespace, *PurePosixPath(relative).parts
    )
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    if not destination.parent.resolve(strict=True).is_relative_to(snapshot_root):
        raise ImplementationFreezeError("snapshot destination escapes evidence root")
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(source, flags)
    try:
        opened = os.fstat(descriptor)
        if _stat_signature(before) != _stat_signature(opened):
            raise ImplementationFreezeIntegrityError(
                f"implementation source changed before copy: {relative}"
            )
        size, digest = _copy_file_bytes(descriptor, destination)
        after_open = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    after_path = os.lstat(source)
    signature = _stat_signature(before)
    if (
        signature != _stat_signature(after_open)
        or signature != _stat_signature(after_path)
        or size != before.st_size
    ):
        destination.unlink(missing_ok=True)
        raise ImplementationFreezeIntegrityError(
            f"implementation source changed during copy: {relative}"
        )
    second_size, second_digest, second_stat = _hash_stable_source(source)
    if (
        second_size != size
        or second_digest != digest
        or _stat_signature(second_stat) != signature
    ):
        destination.unlink(missing_ok=True)
        raise ImplementationFreezeIntegrityError(
            f"implementation source changed during copy verification: {relative}"
        )
    mode = stat.S_IMODE(before.st_mode)
    os.chmod(destination, mode, follow_symlinks=False)
    return ImplementationFreezeFileV1(
        schema_version=_FILE_SCHEMA_VERSION,
        namespace=namespace,
        relative_path=relative,
        size_bytes=size,
        mode=_mode_text(before.st_mode),
        sha256=digest,
    )


def _verify_snapshot_records(
    *, snapshot_root: Path, records: Sequence[ImplementationFreezeFileV1]
) -> str:
    observed: list[ImplementationFreezeFileV1] = []
    for expected in records:
        target = snapshot_root.joinpath(
            expected.namespace, *PurePosixPath(expected.relative_path).parts
        )
        try:
            resolved = target.resolve(strict=True)
        except OSError as exc:
            raise ImplementationFreezeIntegrityError(
                "snapshot file is missing: "
                f"{expected.namespace}/{expected.relative_path}"
            ) from exc
        if (
            target.is_symlink()
            or resolved != target
            or not resolved.is_relative_to(snapshot_root)
        ):
            raise ImplementationFreezeIntegrityError(
                "snapshot path is unsafe: "
                f"{expected.namespace}/{expected.relative_path}"
            )
        size, digest, target_stat = _hash_stable_source(target)
        observed.append(
            ImplementationFreezeFileV1(
                schema_version=_FILE_SCHEMA_VERSION,
                namespace=expected.namespace,
                relative_path=expected.relative_path,
                size_bytes=size,
                mode=_mode_text(target_stat.st_mode),
                sha256=digest,
            )
        )
    if tuple(observed) != tuple(records):
        raise ImplementationFreezeIntegrityError(
            "restored snapshot bytes or modes do not match the manifest"
        )
    return canonical_sha256(tuple(asdict(record) for record in observed))


def _verify_sources_unchanged(
    *,
    workspace: Path,
    expected_paths: Sequence[str],
    records: Sequence[ImplementationFreezeFileV1],
    relevant_paths: Sequence[str],
) -> None:
    observed_paths, _metadata = _collect_source_paths(
        workspace=workspace,
        relevant_paths=relevant_paths,
    )
    if tuple(observed_paths) != tuple(expected_paths):
        raise ImplementationFreezeIntegrityError(
            "implementation source set changed during snapshot"
        )
    by_path = {
        record.relative_path: record
        for record in records
        if record.namespace == _WORKTREE_NAMESPACE
    }
    for relative in observed_paths:
        size, digest, source_stat = _hash_stable_source(
            _resolve_workspace_entry(workspace, relative)
        )
        expected = by_path[relative]
        if (
            size != expected.size_bytes
            or digest != expected.sha256
            or _mode_text(source_stat.st_mode) != expected.mode
        ):
            raise ImplementationFreezeIntegrityError(
                f"implementation source changed after copy: {relative}"
            )


def _verify_campaign_inputs_unchanged(
    *,
    campaign_root: Path,
    expected_paths: Sequence[str],
    records: Sequence[ImplementationFreezeFileV1],
) -> None:
    observed_paths = _collect_campaign_inputs(
        campaign_root=campaign_root, allowlist=expected_paths
    )
    if tuple(observed_paths) != tuple(expected_paths):
        raise ImplementationFreezeIntegrityError(
            "campaign input set changed during snapshot"
        )
    by_path = {
        record.relative_path: record
        for record in records
        if record.namespace == _CAMPAIGN_NAMESPACE
    }
    for relative in observed_paths:
        size, digest, source_stat = _hash_stable_source(
            _resolve_workspace_entry(campaign_root, relative)
        )
        expected = by_path[relative]
        if (
            size != expected.size_bytes
            or digest != expected.sha256
            or _mode_text(source_stat.st_mode) != expected.mode
        ):
            raise ImplementationFreezeIntegrityError(
                f"campaign input changed after copy: {relative}"
            )


def _manifest_from_mapping(value: Mapping[str, Any]) -> ImplementationFreezeManifestV1:
    expected_manifest_fields = {
        "schema_version",
        "snapshot_directory",
        "source_roots",
        "relevant_paths",
        "campaign_input_allowlist",
        "files",
        "file_count",
        "total_size_bytes",
        "tree_sha256",
        "restored_tree_sha256",
        "restoration_verified",
        "git",
        "manifest_sha256",
    }
    expected_file_fields = {
        "schema_version",
        "namespace",
        "relative_path",
        "size_bytes",
        "mode",
        "sha256",
    }
    expected_git_fields = {
        "schema_version",
        "head",
        "branch",
        "status_nonempty",
        "status_entry_count",
        "status_sha256",
        "dirty_diff_sha256",
        "observation_sha256",
    }
    if set(value) != expected_manifest_fields:
        raise ImplementationFreezeError(
            "implementation manifest fields do not match the schema"
        )
    try:
        if not isinstance(value["files"], list) or any(
            not isinstance(item, Mapping) or set(item) != expected_file_fields
            for item in value["files"]
        ):
            raise ImplementationFreezeError(
                "implementation file records do not match the schema"
            )
        if not isinstance(value["git"], Mapping) or set(value["git"]) != (
            expected_git_fields
        ):
            raise ImplementationFreezeError(
                "Git observation fields do not match the schema"
            )
        files = tuple(
            ImplementationFreezeFileV1(**item) for item in value["files"]
        )
        git = GitWorktreeObservationV1(**value["git"])
        return ImplementationFreezeManifestV1(
            schema_version=value["schema_version"],
            snapshot_directory=value["snapshot_directory"],
            source_roots=tuple(value["source_roots"]),
            relevant_paths=tuple(value["relevant_paths"]),
            campaign_input_allowlist=tuple(
                value["campaign_input_allowlist"]
            ),
            files=files,
            file_count=value["file_count"],
            total_size_bytes=value["total_size_bytes"],
            tree_sha256=value["tree_sha256"],
            restored_tree_sha256=value["restored_tree_sha256"],
            restoration_verified=value["restoration_verified"],
            git=git,
            manifest_sha256=value["manifest_sha256"],
        )
    except (KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, ImplementationFreezeError):
            raise
        raise ImplementationFreezeError(
            "implementation manifest has an invalid shape"
        ) from exc


def verify_implementation_freeze(
    manifest_path: str | Path,
) -> ImplementationFreezeManifestV1:
    """Reread and verify a copied implementation snapshot.

    ``manifest_path`` is a runtime locator only; no absolute path is persisted
    in the returned contract.
    """

    path = Path(manifest_path).expanduser()
    if not path.is_absolute():
        path = path.absolute()
    if path.is_symlink() or not path.is_file():
        raise ImplementationFreezeIntegrityError(
            "implementation manifest must be a regular non-symlink file"
        )
    try:
        decoded = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ImplementationFreezeIntegrityError(
            "implementation manifest cannot be read"
        ) from exc
    if not isinstance(decoded, Mapping):
        raise ImplementationFreezeIntegrityError(
            "implementation manifest root must be an object"
        )
    manifest = _manifest_from_mapping(decoded)
    evidence_root = path.parent.resolve(strict=True)
    snapshot_root = evidence_root.joinpath(
        *PurePosixPath(manifest.snapshot_directory).parts
    )
    if (
        snapshot_root.is_symlink()
        or not snapshot_root.is_dir()
        or snapshot_root.resolve(strict=True) != snapshot_root
        or not snapshot_root.is_relative_to(evidence_root)
    ):
        raise ImplementationFreezeIntegrityError(
            "implementation snapshot directory is unsafe"
        )
    observed_snapshot_paths: set[str] = set()
    for directory, directory_names, file_names in os.walk(
        snapshot_root, topdown=True, followlinks=False
    ):
        directory_path = Path(directory)
        for name in tuple(directory_names):
            candidate = directory_path / name
            if candidate.is_symlink():
                raise ImplementationFreezeIntegrityError(
                    "implementation snapshot contains a symlink"
                )
        for name in file_names:
            candidate = directory_path / name
            if candidate.is_symlink() or not candidate.is_file():
                raise ImplementationFreezeIntegrityError(
                    "implementation snapshot contains an unsafe file"
                )
            observed_snapshot_paths.add(
                candidate.relative_to(snapshot_root).as_posix()
            )
    expected_snapshot_paths = {
        f"{record.namespace}/{record.relative_path}"
        for record in manifest.files
    }
    if observed_snapshot_paths != expected_snapshot_paths:
        raise ImplementationFreezeIntegrityError(
            "implementation snapshot contains missing or unexpected files"
        )
    restored = _verify_snapshot_records(
        snapshot_root=snapshot_root, records=manifest.files
    )
    if restored != manifest.restored_tree_sha256:
        raise ImplementationFreezeIntegrityError(
            "implementation restoration digest mismatch"
        )
    return manifest


def create_implementation_freeze(
    *,
    workspace_root: str | Path,
    evidence_directory: str | Path,
    relevant_paths: Iterable[str | Path] = (),
    campaign_root: str | Path | None = None,
    campaign_input_allowlist: Iterable[str | Path] = (),
) -> ImplementationFreezeManifestV1:
    """Copy and verify the exact implementation used by an experiment.

    ``relevant_paths`` names the worktree tests, data, and documentation used by
    the study. ``campaign_root`` is an explicitly approved absolute local root;
    ``campaign_input_allowlist`` admits only named public launcher and
    configuration files beneath it.  The two source namespaces are copied to
    separate subtrees.  Absolute source locators, secrets, provider-private
    traffic, and symlinks never enter the manifest.
    """

    workspace_candidate = Path(workspace_root).expanduser()
    try:
        workspace = workspace_candidate.resolve(strict=True)
    except OSError as exc:
        raise ImplementationFreezeError("workspace root does not exist") from exc
    if not workspace.is_dir() or workspace_candidate.is_symlink():
        raise ImplementationFreezeError(
            "workspace root must be a regular non-symlink directory"
        )
    normalized_relevant = tuple(
        sorted(
            {
                _require_relative_workspace_path(value, "relevant_paths")
                for value in relevant_paths
            }
        )
    )
    normalized_campaign_inputs = tuple(
        sorted(
            {
                _require_relative_workspace_path(
                    value, "campaign_input_allowlist"
                )
                for value in campaign_input_allowlist
            }
        )
    )
    approved_campaign_root: Path | None = None
    if normalized_campaign_inputs and campaign_root is None:
        raise ImplementationFreezeError(
            "campaign root is required for campaign inputs"
        )
    if campaign_root is not None:
        campaign_candidate = Path(campaign_root).expanduser()
        if not campaign_candidate.is_absolute():
            raise ImplementationFreezeError(
                "approved campaign root must be absolute"
            )
        try:
            approved_campaign_root = campaign_candidate.resolve(strict=True)
        except OSError as exc:
            raise ImplementationFreezeError(
                "approved campaign root does not exist"
            ) from exc
        if (
            campaign_candidate.is_symlink()
            or approved_campaign_root != campaign_candidate
            or not approved_campaign_root.is_dir()
        ):
            raise ImplementationFreezeError(
                "approved campaign root must be a non-symlink directory"
            )
        if approved_campaign_root.is_relative_to(workspace) or workspace.is_relative_to(
            approved_campaign_root
        ):
            raise ImplementationFreezeError(
                "campaign and implementation roots must be disjoint"
            )

    evidence_candidate = Path(evidence_directory).expanduser()
    if evidence_candidate.exists() and evidence_candidate.is_symlink():
        raise ImplementationFreezeError(
            "evidence directory must not be a symlink"
        )
    evidence_candidate.mkdir(mode=0o700, parents=True, exist_ok=True)
    evidence = evidence_candidate.resolve(strict=True)
    if evidence.is_relative_to(workspace):
        raise ImplementationFreezeError(
            "evidence directory must remain outside the implementation workspace"
        )
    snapshot_root = evidence / _SNAPSHOT_DIRECTORY
    manifest_path = evidence / _MANIFEST_NAME
    if snapshot_root.exists() or manifest_path.exists():
        raise ImplementationFreezeError(
            "implementation freeze destination already exists"
        )

    git_before = _observe_git(workspace)
    source_paths, metadata_paths = _collect_source_paths(
        workspace=workspace,
        relevant_paths=normalized_relevant,
    )
    campaign_paths = (
        _collect_campaign_inputs(
            campaign_root=approved_campaign_root,
            allowlist=normalized_campaign_inputs,
        )
        if approved_campaign_root is not None
        else ()
    )
    snapshot_root.mkdir(mode=0o700)
    records: list[ImplementationFreezeFileV1] = []
    try:
        for relative in source_paths:
            records.append(
                _copy_stable_source(
                    source_root=workspace,
                    namespace=_WORKTREE_NAMESPACE,
                    relative=relative,
                    snapshot_root=snapshot_root,
                )
            )
        if approved_campaign_root is not None:
            for relative in campaign_paths:
                records.append(
                    _copy_stable_source(
                        source_root=approved_campaign_root,
                        namespace=_CAMPAIGN_NAMESPACE,
                        relative=relative,
                        snapshot_root=snapshot_root,
                    )
                )
        frozen_records = tuple(
            sorted(records, key=lambda item: (item.namespace, item.relative_path))
        )
        _verify_sources_unchanged(
            workspace=workspace,
            expected_paths=source_paths,
            records=frozen_records,
            relevant_paths=normalized_relevant,
        )
        if approved_campaign_root is not None:
            _verify_campaign_inputs_unchanged(
                campaign_root=approved_campaign_root,
                expected_paths=campaign_paths,
                records=frozen_records,
            )
        git_after = _observe_git(workspace)
        if git_after != git_before:
            raise ImplementationFreezeIntegrityError(
                "Git worktree observation changed during snapshot"
            )
        restored_tree = _verify_snapshot_records(
            snapshot_root=snapshot_root, records=frozen_records
        )
        source_roots = tuple(
            sorted(set(_DEFAULT_SOURCE_ROOTS).union(metadata_paths))
        )
        body = {
            "schema_version": _SCHEMA_VERSION,
            "snapshot_directory": _SNAPSHOT_DIRECTORY,
            "source_roots": source_roots,
            "relevant_paths": normalized_relevant,
            "campaign_input_allowlist": normalized_campaign_inputs,
            "files": frozen_records,
            "file_count": len(frozen_records),
            "total_size_bytes": sum(item.size_bytes for item in frozen_records),
            "tree_sha256": canonical_sha256(
                tuple(asdict(record) for record in frozen_records)
            ),
            "restored_tree_sha256": restored_tree,
            "restoration_verified": True,
            "git": git_before,
        }
        manifest = ImplementationFreezeManifestV1(
            **body, manifest_sha256=canonical_sha256(body)
        )
        temporary_manifest = evidence / f".{_MANIFEST_NAME}.tmp"
        with temporary_manifest.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(canonical_json(manifest.public_record()))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary_manifest.replace(manifest_path)
        verified = verify_implementation_freeze(manifest_path)
        if verified != manifest:
            raise ImplementationFreezeIntegrityError(
                "implementation manifest changed after publication"
            )
        return manifest
    except BaseException:
        # A failed directory is intentionally retained as non-authoritative
        # forensic evidence.  It has no published manifest and cannot be used
        # as an experiment freeze.
        raise


__all__ = [
    "GitWorktreeObservationV1",
    "ImplementationFreezeError",
    "ImplementationFreezeFileV1",
    "ImplementationFreezeIntegrityError",
    "ImplementationFreezeManifestV1",
    "create_implementation_freeze",
    "verify_implementation_freeze",
]
