#!/usr/bin/env python3
"""Audit agent file lineage and source similarity against pinned repositories."""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable, Sequence


TEXT_SUFFIXES = {
    ".css",
    ".js",
    ".jsx",
    ".md",
    ".py",
    ".rs",
    ".scss",
    ".sh",
    ".tcss",
    ".toml",
    ".ts",
    ".tsx",
    ".txt",
    ".yaml",
    ".yml",
}
SKIP_PARTS = {
    ".git",
    "__pycache__",
    "assets",
    "dist",
    "fixtures",
    "generated",
    "node_modules",
    "snapshots",
    "target",
    "vendor",
}
SKIP_NAMES = {
    "bun.lock",
    "package-lock.json",
    "pnpm-lock.yaml",
    "uv.lock",
    "yarn.lock",
}
TOKEN_RE = re.compile(
    r"[A-Za-z_][A-Za-z0-9_]*|\d+(?:\.\d+)?|"
    r"==|!=|<=|>=|:=|->|=>|::|&&|\|\|"
)
BLAME_HEADER_RE = re.compile(r"^([0-9a-f^]{40}) \d+ \d+(?: \d+)?$")


@dataclass(frozen=True)
class TextDocument:
    source: str
    path: str
    text: str

    @property
    def tokens(self) -> tuple[str, ...]:
        return tokenize(self.text)

    @property
    def nonblank_lines(self) -> tuple[tuple[int, str], ...]:
        rows = []
        for line_number, line in enumerate(self.text.splitlines(), start=1):
            normalized = normalize_line(line)
            if normalized:
                rows.append((line_number, normalized))
        return tuple(rows)


@dataclass(frozen=True)
class Match:
    kind: str
    candidate_path: str
    candidate_start: int
    source: str
    source_path: str
    source_start: int
    units: int
    tokens: int
    similarity: float


@dataclass(frozen=True)
class Unit:
    document_index: int
    start: int
    tokens: tuple[str, ...]


def run_git(repo: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize_line(line: str) -> str:
    return " ".join(line.strip().split())


def tokenize(text: str) -> tuple[str, ...]:
    return tuple(token.lower() for token in TOKEN_RE.findall(text))


def token_hash(tokens: Sequence[str]) -> str:
    payload = "\x1f".join(tokens).encode("utf-8")
    return hashlib.blake2b(payload, digest_size=12).hexdigest()


def tracked_paths(repo: Path, *prefixes: str) -> list[Path]:
    output = run_git(repo, "ls-files", "--", *prefixes)
    return [repo / line for line in output.splitlines() if line]


def is_text_candidate(path: Path) -> bool:
    return (
        path.suffix.lower() in TEXT_SUFFIXES
        and path.name not in SKIP_NAMES
        and not any(part in SKIP_PARTS for part in path.parts)
    )


def read_text_document(root: Path, path: Path, source: str) -> TextDocument | None:
    if not is_text_candidate(path) or path.stat().st_size > 1_000_000:
        return None
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None
    return TextDocument(source, path.relative_to(root).as_posix(), text)


def candidate_documents(repo: Path) -> list[TextDocument]:
    paths = tracked_paths(
        repo,
        "chemsmart/agent",
        "docs",
        "README.md",
    )
    documents: list[TextDocument] = []
    for path in paths:
        relative = path.relative_to(repo).as_posix()
        if (
            relative.startswith("docs/")
            and not relative.startswith("docs/research/")
            and not any(
                keyword in relative.lower()
                for keyword in ("agent", "harness", "tui", "wizard")
            )
        ):
            continue
        document = read_text_document(repo, path, "chemsmart")
        if document is not None:
            documents.append(document)
    return documents


def source_documents(root: Path, source: str) -> list[TextDocument]:
    documents: list[TextDocument] = []
    for path in tracked_paths(root, "."):
        document = read_text_document(root, path, source)
        if document is not None:
            documents.append(document)
    return documents


def _commit_rows(repo: Path, path: str) -> list[dict[str, str]]:
    output = run_git(
        repo,
        "log",
        "--follow",
        "--reverse",
        "--format=%H%x09%aI%x09%an%x09%s",
        "--",
        path,
    )
    rows: list[dict[str, str]] = []
    for line in output.splitlines():
        commit, date, author, subject = line.split("\t", 3)
        rows.append(
            {
                "commit": commit,
                "date": date,
                "author": author,
                "subject": subject,
            }
        )
    return rows


def _blame_counts(repo: Path, head: str, path: str) -> list[dict[str, object]]:
    output = run_git(repo, "blame", "--line-porcelain", head, "--", path)
    counts: Counter[str] = Counter()
    authors: dict[str, str] = {}
    current_commit: str | None = None
    for line in output.splitlines():
        header = BLAME_HEADER_RE.fullmatch(line)
        if header:
            current_commit = header.group(1).lstrip("^")
            counts[current_commit] += 1
        elif line.startswith("author ") and current_commit:
            authors.setdefault(current_commit, line.removeprefix("author "))
    return [
        {"commit": commit, "author": authors.get(commit), "lines": lines}
        for commit, lines in counts.most_common()
    ]


def production_lineage(repo: Path, head: str, boundary: str) -> list[dict[str, object]]:
    boundary_commit = run_git(repo, "rev-parse", boundary)
    rows: list[dict[str, object]] = []
    for path in tracked_paths(repo, "chemsmart/agent"):
        if path.suffix != ".py":
            continue
        relative = path.relative_to(repo).as_posix()
        commits = _commit_rows(repo, relative)
        rows.append(
            {
                "path": relative,
                "lines": len(path.read_text(encoding="utf-8").splitlines()),
                "introduction": commits[0] if commits else None,
                "latest_change": commits[-1] if commits else None,
                "commit_count": len(commits),
                "introduced_after_boundary": bool(
                    commits
                    and subprocess.run(
                        [
                            "git",
                            "merge-base",
                            "--is-ancestor",
                            boundary_commit,
                            commits[0]["commit"],
                        ],
                        cwd=repo,
                        check=False,
                    ).returncode
                    == 0
                ),
                "blame": _blame_counts(repo, head, relative),
            }
        )
    return rows


def _line_shingle_index(
    documents: Sequence[TextDocument],
    width: int,
) -> dict[str, list[tuple[int, int]]]:
    index: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for document_index, document in enumerate(documents):
        lines = document.nonblank_lines
        for position in range(len(lines) - width + 1):
            key = token_hash(tuple(value for _, value in lines[position : position + width]))
            index[key].append((document_index, position))
    return index


def exact_line_matches(
    candidates: Sequence[TextDocument],
    externals: Sequence[TextDocument],
    *,
    min_lines: int = 6,
    min_tokens: int = 50,
) -> list[Match]:
    candidate_index = _line_shingle_index(candidates, min_lines)
    seen: set[tuple[int, int, int, int, int]] = set()
    matches: list[Match] = []
    for external_index, external in enumerate(externals):
        external_lines = external.nonblank_lines
        for source_position in range(len(external_lines) - min_lines + 1):
            key = token_hash(
                tuple(
                    value
                    for _, value in external_lines[
                        source_position : source_position + min_lines
                    ]
                )
            )
            for candidate_index_value, candidate_position in candidate_index.get(key, ()):
                candidate_lines = candidates[candidate_index_value].nonblank_lines
                left = 0
                while (
                    candidate_position - left > 0
                    and source_position - left > 0
                    and candidate_lines[candidate_position - left - 1][1]
                    == external_lines[source_position - left - 1][1]
                ):
                    left += 1
                right = min_lines
                while (
                    candidate_position + right < len(candidate_lines)
                    and source_position + right < len(external_lines)
                    and candidate_lines[candidate_position + right][1]
                    == external_lines[source_position + right][1]
                ):
                    right += 1
                start_candidate = candidate_position - left
                start_source = source_position - left
                key_tuple = (
                    candidate_index_value,
                    start_candidate,
                    external_index,
                    start_source,
                    right + left,
                )
                if key_tuple in seen:
                    continue
                seen.add(key_tuple)
                block = "\n".join(
                    value
                    for _, value in candidate_lines[
                        start_candidate : candidate_position + right
                    ]
                )
                tokens = len(tokenize(block))
                units = right + left
                if units < min_lines or tokens < min_tokens:
                    continue
                candidate = candidates[candidate_index_value]
                matches.append(
                    Match(
                        "exact_lines",
                        candidate.path,
                        candidate_lines[start_candidate][0],
                        external.source,
                        external.path,
                        external_lines[start_source][0],
                        units,
                        tokens,
                        1.0,
                    )
                )
    return _drop_contained_matches(matches)


def _token_anchor_index(
    documents: Sequence[TextDocument],
    width: int,
) -> dict[str, list[tuple[int, int]]]:
    index: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for document_index, document in enumerate(documents):
        tokens = document.tokens
        for position in range(len(tokens) - width + 1):
            index[token_hash(tokens[position : position + width])].append(
                (document_index, position)
            )
    return index


def exact_token_matches(
    candidates: Sequence[TextDocument],
    externals: Sequence[TextDocument],
    *,
    min_tokens: int = 50,
    anchor_tokens: int = 12,
) -> list[Match]:
    candidate_index = _token_anchor_index(candidates, anchor_tokens)
    seen: set[tuple[int, int, int, int, int]] = set()
    matches: list[Match] = []
    for external_index, external in enumerate(externals):
        external_tokens = external.tokens
        for source_position in range(len(external_tokens) - anchor_tokens + 1):
            anchor = token_hash(
                external_tokens[source_position : source_position + anchor_tokens]
            )
            occurrences = candidate_index.get(anchor, ())
            if len(occurrences) > 200:
                continue
            for candidate_index_value, candidate_position in occurrences:
                candidate_tokens = candidates[candidate_index_value].tokens
                left = 0
                while (
                    candidate_position - left > 0
                    and source_position - left > 0
                    and candidate_tokens[candidate_position - left - 1]
                    == external_tokens[source_position - left - 1]
                ):
                    left += 1
                right = anchor_tokens
                while (
                    candidate_position + right < len(candidate_tokens)
                    and source_position + right < len(external_tokens)
                    and candidate_tokens[candidate_position + right]
                    == external_tokens[source_position + right]
                ):
                    right += 1
                units = left + right
                if units < min_tokens:
                    continue
                start_candidate = candidate_position - left
                start_source = source_position - left
                key = (
                    candidate_index_value,
                    start_candidate,
                    external_index,
                    start_source,
                    units,
                )
                if key in seen:
                    continue
                seen.add(key)
                matches.append(
                    Match(
                        "exact_tokens",
                        candidates[candidate_index_value].path,
                        start_candidate + 1,
                        external.source,
                        external.path,
                        start_source + 1,
                        units,
                        units,
                        1.0,
                    )
                )
    return _drop_contained_matches(matches)


def _python_units(document_index: int, text: str) -> list[Unit]:
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    lines = text.splitlines()
    units: list[Unit] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        end = int(getattr(node, "end_lineno", node.lineno))
        tokens = tokenize("\n".join(lines[node.lineno - 1 : end]))
        if len(tokens) >= 100:
            units.extend(_split_unit(document_index, node.lineno, tokens))
    return units


def _split_unit(document_index: int, start: int, tokens: tuple[str, ...]) -> list[Unit]:
    if len(tokens) <= 800:
        return [Unit(document_index, start, tokens)]
    return [
        Unit(document_index, start + position, tokens[position : position + 400])
        for position in range(0, len(tokens) - 99, 200)
    ]


def comparison_units(documents: Sequence[TextDocument]) -> list[Unit]:
    units: list[Unit] = []
    for document_index, document in enumerate(documents):
        if document.path.endswith(".py"):
            units.extend(_python_units(document_index, document.text))
            continue
        tokens = document.tokens
        for width in (120, 240, 480):
            if len(tokens) < width:
                continue
            for position in range(0, len(tokens) - width + 1, width // 2):
                units.append(Unit(document_index, position + 1, tokens[position : position + width]))
    return units


def near_token_matches(
    candidates: Sequence[TextDocument],
    externals: Sequence[TextDocument],
    *,
    minimum_similarity: float = 0.85,
    anchor_tokens: int = 8,
) -> list[Match]:
    candidate_units = comparison_units(candidates)
    external_units = comparison_units(externals)
    anchors: dict[str, set[int]] = defaultdict(set)
    for unit_index, unit in enumerate(candidate_units):
        for position in range(0, len(unit.tokens) - anchor_tokens + 1, 4):
            anchors[token_hash(unit.tokens[position : position + anchor_tokens])].add(
                unit_index
            )
    seen_pairs: set[tuple[int, int]] = set()
    matches: list[Match] = []
    for external_unit_index, external_unit in enumerate(external_units):
        candidates_for_unit: Counter[int] = Counter()
        for position in range(0, len(external_unit.tokens) - anchor_tokens + 1, 4):
            key = token_hash(
                external_unit.tokens[position : position + anchor_tokens]
            )
            unit_ids = anchors.get(key, ())
            if len(unit_ids) <= 40:
                candidates_for_unit.update(unit_ids)
        for candidate_unit_index, shared_anchors in candidates_for_unit.most_common(20):
            if shared_anchors < 2:
                continue
            pair = (candidate_unit_index, external_unit_index)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            candidate_unit = candidate_units[candidate_unit_index]
            length_ratio = min(len(candidate_unit.tokens), len(external_unit.tokens)) / max(
                len(candidate_unit.tokens), len(external_unit.tokens)
            )
            if length_ratio < minimum_similarity:
                continue
            matcher = SequenceMatcher(
                None,
                candidate_unit.tokens,
                external_unit.tokens,
                autojunk=False,
            )
            if matcher.quick_ratio() < minimum_similarity:
                continue
            ratio = matcher.ratio()
            if ratio < minimum_similarity:
                continue
            candidate = candidates[candidate_unit.document_index]
            external = externals[external_unit.document_index]
            matches.append(
                Match(
                    "near_tokens",
                    candidate.path,
                    candidate_unit.start,
                    external.source,
                    external.path,
                    external_unit.start,
                    min(len(candidate_unit.tokens), len(external_unit.tokens)),
                    min(len(candidate_unit.tokens), len(external_unit.tokens)),
                    ratio,
                )
            )
    return sorted(
        matches,
        key=lambda match: (-match.similarity, -match.tokens, match.candidate_path),
    )


def _drop_contained_matches(matches: Iterable[Match]) -> list[Match]:
    ordered = sorted(
        matches,
        key=lambda match: (-match.units, match.candidate_path, match.candidate_start),
    )
    kept: list[Match] = []
    for match in ordered:
        if any(
            current.kind == match.kind
            and current.candidate_path == match.candidate_path
            and current.source == match.source
            and current.source_path == match.source_path
            and current.candidate_start <= match.candidate_start
            and current.candidate_start + current.units
            >= match.candidate_start + match.units
            for current in kept
        ):
            continue
        kept.append(match)
    return sorted(kept, key=lambda item: (item.candidate_path, item.candidate_start))


def load_and_verify_sources(
    manifest_path: Path,
    source_root: Path,
) -> tuple[list[dict[str, str]], dict[str, Path]]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = payload["sources"]
    roots: dict[str, Path] = {}
    for row in rows:
        root = source_root / row["directory"]
        actual_commit = run_git(root, "rev-parse", "HEAD")
        if actual_commit != row["commit"]:
            raise ValueError(
                f"{row['name']} commit mismatch: {actual_commit} != {row['commit']}"
            )
        license_path = root / row["license_path"]
        actual_license_sha = sha256(license_path)
        if actual_license_sha != row["license_sha256"]:
            raise ValueError(f"{row['name']} license hash mismatch")
        roots[row["name"]] = root
    return rows, roots


def build_report(
    repo: Path,
    boundary: str,
    manifest_path: Path,
    source_root: Path,
) -> dict[str, object]:
    head = run_git(repo, "rev-parse", "HEAD")
    source_rows, roots = load_and_verify_sources(manifest_path, source_root)
    candidates = candidate_documents(repo)
    externals: list[TextDocument] = []
    source_counts: dict[str, int] = {}
    for name, root in roots.items():
        documents = source_documents(root, name)
        source_counts[name] = len(documents)
        externals.extend(documents)
    exact_lines = exact_line_matches(candidates, externals)
    exact_tokens = exact_token_matches(candidates, externals)
    near_tokens = near_token_matches(candidates, externals)
    lineage = production_lineage(repo, head, boundary)
    return {
        "schema_version": 1,
        "audit": {
            "head": head,
            "boundary": run_git(repo, "rev-parse", boundary),
            "candidate_documents": len(candidates),
            "production_python_files": len(lineage),
            "external_documents": source_counts,
            "thresholds": {
                "exact_nonblank_lines": 6,
                "exact_tokens": 50,
                "near_tokens": 100,
                "near_similarity": 0.85,
            },
        },
        "sources": source_rows,
        "lineage": lineage,
        "matches": {
            "exact_lines": [asdict(match) for match in exact_lines],
            "exact_tokens": [asdict(match) for match in exact_tokens],
            "near_tokens": [asdict(match) for match in near_tokens],
        },
    }


def markdown(report: dict[str, object]) -> str:
    audit = report["audit"]
    sources = report["sources"]
    matches = report["matches"]
    lineage = report["lineage"]
    assert isinstance(audit, dict)
    assert isinstance(sources, list)
    assert isinstance(matches, dict)
    assert isinstance(lineage, list)
    lines = [
        "# ChemSmart Agent Provenance Audit",
        "",
        (
            "This is an evidence-based source-similarity review, not a legal "
            "opinion. A match is a review lead; absence of a match is not proof "
            "of independent authorship."
        ),
        "",
        f"- Audit HEAD: `{audit['head']}`",
        f"- Agent-history boundary: `{audit['boundary']}`",
        f"- Production Python files traced: {audit['production_python_files']}",
        f"- Candidate code/document files compared: {audit['candidate_documents']}",
        "",
        "## Pinned Reference Snapshots",
        "",
        "| Source | Commit | License evidence | Review posture |",
        "|---|---|---|---|",
    ]
    for source in sources:
        lines.append(
            f"| {source['name']} | `{source['commit']}` | "
            f"`{source['license_path']}` / `{source['license_sha256'][:12]}...` | "
            f"{source['review_posture']} |"
        )
    lines.extend(
        [
            "",
            "## Automated Similarity Leads",
            "",
            "| Match class | Count | Threshold |",
            "|---|---:|---|",
            f"| Exact nonblank lines | {len(matches['exact_lines'])} | >=6 lines and >=50 tokens |",
            f"| Exact token sequence | {len(matches['exact_tokens'])} | >=50 normalized tokens |",
            f"| Near token unit | {len(matches['near_tokens'])} | >=100 tokens and ratio >=0.85 |",
            "",
        ]
    )
    all_matches = [
        row
        for category in ("exact_lines", "exact_tokens", "near_tokens")
        for row in matches[category]
    ]
    if all_matches:
        lines.extend(
            [
                "Every lead requires manual context review before attribution or "
                "clean-room action is decided.",
                "",
                "| Kind | Candidate | Reference | Tokens | Similarity |",
                "|---|---|---|---:|---:|",
            ]
        )
        for row in all_matches:
            lines.append(
                f"| {row['kind']} | `{row['candidate_path']}:{row['candidate_start']}` | "
                f"{row['source']} `{row['source_path']}:{row['source_start']}` | "
                f"{row['tokens']} | {row['similarity']:.3f} |"
            )
    else:
        lines.extend(
            [
                "No automated lead reached the configured thresholds. This does "
                "not replace manual review of history, comments, design citations, "
                "or lower-similarity adaptations.",
                "",
            ]
        )
    introductions = Counter(
        row["introduction"]["author"]
        for row in lineage
        if isinstance(row.get("introduction"), dict)
    )
    before_boundary = [
        row for row in lineage if not row["introduced_after_boundary"]
    ]
    lines.extend(
        [
            "## Lineage Summary",
            "",
            f"- Files whose first recorded commit predates the boundary: {len(before_boundary)}",
            "- Introduction authors: "
            + ", ".join(
                f"{author} ({count})" for author, count in introductions.most_common()
            ),
            "- Full per-file introduction, latest-change, and blame distributions "
            "are retained in `agent-provenance-audit.json`.",
            "",
            "## Interpretation Rules",
            "",
            "- Conceptual similarity is recorded as design influence, not copied code.",
            "- Compatible reuse still requires the applicable notices and attribution.",
            "- Restrictive, unknown, or unattributed matching blocks require behavior "
            "tests followed by clean-room replacement before release.",
            "- No legal conclusion is made by this report.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--boundary", default="7f1cb674")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    report = build_report(
        repo,
        args.boundary,
        args.manifest.resolve(),
        args.source_root.resolve(),
    )
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(payload, encoding="utf-8")
    args.markdown.write_text(markdown(report), encoding="utf-8")
    print(f"wrote {args.json} and {args.markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
