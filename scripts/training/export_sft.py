"""Export accumulated ChemSmart agent episodes to SFT JSONL files.

The raw cross-session logs are evidence ledgers. This exporter turns them into
separate, conservative training families:

* ``tool_loop_sft.jsonl`` preserves API/frontier tool-use conversations.
* ``compact_spec_sft.jsonl`` keeps local-model compact SPEC targets only after
  adapter round-trip validation.
* ``command_answer_sft.jsonl`` captures grounded command explanations.
* ``project_yaml_sft.jsonl`` captures project-YAML build/validation turns.
* ``reasoning_synthesis_review.jsonl`` keeps non-training reasoning candidates
  with skip reasons for later audit.
* ``wrong_route_contrast.jsonl`` pairs project-YAML wrong-route episodes with
  corrected command-synthesis episodes of the same workflow/kind.
* ``repair_pairs.jsonl`` stores reject->repair contrast examples.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shlex
import tempfile
from collections import Counter
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from chemsmart.agent.harness.terminal_state import (
    terminal_state_is_positive,
    validate_terminal_state,
)

JsonDict = dict[str, Any]
_PUBLIC_MESSAGE_KEYS = {
    "role",
    "content",
    "tool_calls",
    "tool_call_id",
    "name",
    "function_call",
}

_POSITIVE_GATES = {"ok", "warn", "none", None}
_COMMAND_POSITIVE_GATES = {"ok", "warn"}
# Canonical-rule violations that must never reach positive SFT records:
# freq is a project-YAML (runtime-owned) setting — a command that smuggles
# it through route parameters or a synthetic --freq flag teaches the model
# the opt+freq anti-pattern; explicit method flags are likewise yaml-owned.
_CANONICAL_VIOLATIONS = (
    (
        "freq_in_route",
        re.compile(
            r"(?:-r|--additional-route-parameters)\s+\S*freq|\bfreq=|--freq\b",
            re.I,
        ),
    ),
    ("method_flag", re.compile(r"--functional\b|--basis\b")),
    (
        "empty_qmmm_layer_override",
        re.compile(
            r"(?:-mx|--medium-level-functional|-mb|--medium-level-basis)\s+(?:\"\"|''|\"\s*\")",
            re.I,
        ),
    ),
)


# Runtime-owned (project-YAML) fields the model must never emit as CLI flags.
# The regex guard above only catches LONG forms (--functional/--basis); the
# short forms -x/-b/-sm/-si are real Click options (cli/gaussian/gaussian.py,
# cli/orca/orca.py) so they parse cleanly and would otherwise reach positive
# SFT. We reuse the real parser (which resolves -x->functional, -b->basis) so
# short forms are caught AND scan's numeric -x (dist_start) is NOT a false
# positive. Set matches the existing runtime-owned policy in
# `_compact_spec_from_command`, applied consistently to every positive family.
# Solvent flags (`-sm/--solvent-model`, `-si/--solvent-id`, `-so/--solvent-
# options`) are real ORCA subcommand options (cli/orca/orca.py:214-236) but the
# parser does NOT map them to solvent_model/solvent_id, so a parser-field check
# misses them — a regex is required. Solvent (theory-level) is project-YAML-
# owned, matching `_compact_spec_from_command`. `--states`/`--nstates` are
# deliberately NOT here: the number of excited states is a job-level structural
# parameter, not a runtime-owned theory field, so the model may emit it.
_SOLVENT_FLAG_RE = re.compile(
    r"(?<!\S)(?:-sm|-si|-so|--solvent-model|--solvent-id|--solvent-options)(?=[\s=]|$)"
)


def _runtime_owned_leak(command: str) -> str | None:
    from chemsmart.agent.model_command_parser import parse_model_command

    try:
        parsed = parse_model_command(command)
    except Exception:
        parsed = None
    if parsed is not None and not getattr(parsed, "parse_error", None):
        if parsed.functional or parsed.ab_initio:
            return "runtime_owned_method_flag"
        if parsed.basis or parsed.aux_basis or parsed.extrapolation_basis:
            return "runtime_owned_basis_flag"
        if parsed.solvent_model or parsed.solvent_id:
            return "runtime_owned_solvent_flag"
    if _SOLVENT_FLAG_RE.search(command):
        return "runtime_owned_solvent_flag"
    return None


def _canonical_violation(command: str) -> str | None:
    for reason, pattern in _CANONICAL_VIOLATIONS:
        if pattern.search(command):
            return reason
    return _runtime_owned_leak(command)


def _semantic_pair_violation(user: str, command: str) -> str | None:
    """Catch request/command mismatches that parser-only gates cannot see."""

    if re.search(r"\bmax\s*steps?\b|\bmaxstep\b", user, re.I) and re.search(
        r"\bmaxcycle\s*=", command, re.I
    ):
        return "maxstep_mapped_to_maxcycle"
    return None


_SECRET_RE = re.compile(
    r"\b(?:hf_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9_-]{20,}|sk-proj-[A-Za-z0-9_-]{20,})\b"
)
_TEMP_PATH_RE = re.compile(
    r"(?:/private)?/var/folders/[^\s`'\"),;]+|(?:/private)?/tmp/[^\s`'\"),;]+"
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "episodes",
        nargs="*",
        help=(
            "Episode JSONL files. Defaults to the configured ChemSmart "
            "agent-training store."
        ),
    )
    parser.add_argument(
        "--training-dir",
        default=None,
        help=(
            "Directory containing prompt hashes and default episode logs. "
            "Defaults to the repo-local var/agent-training store."
        ),
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help=(
            "Output directory. Defaults to "
            "<training-dir>/exports/<UTC timestamp>."
        ),
    )
    parser.add_argument(
        "--include-runs",
        action="store_true",
        help=(
            "Also read model-specific episode stores under "
            "<training-dir>/runs/*/episodes/*.jsonl."
        ),
    )
    parser.add_argument(
        "--include-rejected",
        action="store_true",
        help="Include rejected semantic-gate examples in tool_loop_sft.jsonl.",
    )
    parser.add_argument(
        "--allow-frontier-distill",
        action="store_true",
        help=(
            "Permit an optional configured frontier provider call to distill "
            "a compact SPEC when deterministic recovery fails."
        ),
    )
    args = parser.parse_args(argv)

    training_dir = _training_dir(args.training_dir)
    out_dir = (
        Path(args.out_dir).expanduser()
        if args.out_dir
        else training_dir
        / "exports"
        / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    )
    counts = export_sft(
        episode_paths=_episode_paths(
            args.episodes,
            training_dir,
            include_runs=args.include_runs,
        ),
        training_dir=training_dir,
        out_dir=out_dir,
        include_rejected=args.include_rejected,
        allow_frontier_distill=args.allow_frontier_distill,
    )
    print(json.dumps(counts, sort_keys=True))
    return 0


def export_sft(
    *,
    episode_paths: Iterable[Path],
    training_dir: Path,
    out_dir: Path,
    include_rejected: bool = False,
    allow_frontier_distill: bool = False,
) -> JsonDict:
    """Export episode logs and return a manifest dict."""

    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "tool_loop": out_dir / "tool_loop_sft.jsonl",
        "tool_loop_review": out_dir / "tool_loop_review.jsonl",
        "compact_spec": out_dir / "compact_spec_sft.jsonl",
        "reasoning_synthesis": out_dir / "reasoning_synthesis_sft.jsonl",
        "reasoning_synthesis_review": out_dir
        / "reasoning_synthesis_review.jsonl",
        "wrong_route_contrast": out_dir / "wrong_route_contrast.jsonl",
        "command_answer": out_dir / "command_answer_sft.jsonl",
        "project_yaml": out_dir / "project_yaml_sft.jsonl",
        "repair_pairs": out_dir / "repair_pairs.jsonl",
        "terminal_state": out_dir / "terminal_state_assertions.jsonl",
    }
    handles = {
        name: path.open("w", encoding="utf-8") for name, path in paths.items()
    }
    stats = _ExportStats()
    try:
        episodes = list(_dedup_episodes(_iter_episodes(episode_paths)))
        chains = list(_session_chains(episodes))
        repair_pair_keys: set[tuple[str, str]] = set()

        for chain in chains:
            stats.note_chain(chain)
            tool_record, reason = _tool_loop_chain_record(
                chain,
                training_dir=training_dir,
                include_rejected=include_rejected,
            )
            if tool_record is not None:
                _write_jsonl(handles["tool_loop"], tool_record)
                stats.written["tool_loop"] += 1
            else:
                stats.skipped[f"tool_loop:{reason}"] += 1
                review = _tool_loop_review_record(chain, reason=reason)
                if review is not None:
                    _write_jsonl(handles["tool_loop_review"], review)
                    stats.written["tool_loop_review"] += 1

        for episode in episodes:
            stats.episodes_seen += 1
            stats.note_episode(episode)

            terminal_record, terminal_reason = _terminal_state_record(episode)
            if terminal_record is not None:
                _write_jsonl(handles["terminal_state"], terminal_record)
                stats.written["terminal_state"] += 1
            elif episode.get("terminal_state") is not None:
                stats.skipped[f"terminal_state:{terminal_reason}"] += 1

            compact_record, reason = _compact_spec_record(
                episode,
                allow_frontier_distill=allow_frontier_distill,
            )
            if compact_record is not None:
                _write_jsonl(handles["compact_spec"], compact_record)
                stats.written["compact_spec"] += 1
                stats.note_compact(compact_record)
            else:
                stats.skipped[f"compact_spec:{reason}"] += 1

            reasoning_record, reason = _reasoning_synthesis_record(episode)
            if reasoning_record is not None:
                _write_jsonl(handles["reasoning_synthesis"], reasoning_record)
                stats.written["reasoning_synthesis"] += 1
            else:
                stats.skipped[f"reasoning_synthesis:{reason}"] += 1
                review_record = _reasoning_synthesis_review_record(
                    episode, skip_reason=reason
                )
                if review_record is not None:
                    _write_jsonl(
                        handles["reasoning_synthesis_review"],
                        review_record,
                    )
                    stats.written["reasoning_synthesis_review"] += 1

            command_record, reason = _command_answer_record(episode)
            if command_record is not None:
                _write_jsonl(handles["command_answer"], command_record)
                stats.written["command_answer"] += 1
            else:
                stats.skipped[f"command_answer:{reason}"] += 1

            yaml_record, reason = _project_yaml_record(episode)
            if yaml_record is not None:
                _write_jsonl(handles["project_yaml"], yaml_record)
                stats.written["project_yaml"] += 1
            else:
                stats.skipped[f"project_yaml:{reason}"] += 1

            repair_record, reason = _repair_pair_record(episode)
            if repair_record is not None:
                pair_key = _repair_pair_key(repair_record)
                if pair_key not in repair_pair_keys:
                    repair_pair_keys.add(pair_key)
                    _write_jsonl(handles["repair_pairs"], repair_record)
                    stats.written["repair_pairs"] += 1
            else:
                stats.skipped[f"repair_pairs:{reason}"] += 1

        for chain in chains:
            repair_record, reason = _cross_turn_repair_pair_record(chain)
            if repair_record is None:
                stats.skipped[f"repair_pairs_cross_turn:{reason}"] += 1
                continue
            pair_key = _repair_pair_key(repair_record)
            if pair_key in repair_pair_keys:
                continue
            repair_pair_keys.add(pair_key)
            _write_jsonl(handles["repair_pairs"], repair_record)
            stats.written["repair_pairs"] += 1

        for contrast_record in _wrong_route_contrast_records(episodes):
            _write_jsonl(handles["wrong_route_contrast"], contrast_record)
            stats.written["wrong_route_contrast"] += 1
    finally:
        for handle in handles.values():
            handle.close()

    manifest = stats.manifest()
    manifest["output_dir"] = _portable_path(out_dir, training_dir)
    manifest["files"] = {
        name: _portable_path(path, training_dir)
        for name, path in paths.items()
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return _legacy_counts(manifest)


class _ExportStats:
    def __init__(self) -> None:
        self.episodes_seen = 0
        self.session_chains_seen = 0
        self.multi_turn_chains = 0
        self.written: Counter[str] = Counter()
        self.skipped: Counter[str] = Counter()
        self.query_skeletons: Counter[str] = Counter()
        self.tool_trajectories: Counter[str] = Counter()
        self.job_kinds: Counter[str] = Counter()
        self.providers: Counter[str] = Counter()

    def note_chain(self, chain: list[JsonDict]) -> None:
        self.session_chains_seen += 1
        if _user_message_count(_merge_chain_messages(chain)) >= 2:
            self.multi_turn_chains += 1

    def note_episode(self, episode: JsonDict) -> None:
        user = _first_user_message(episode.get("messages"))
        if user:
            self.query_skeletons[_query_skeleton(user)] += 1
        tools = _invoked_tools(episode)
        if tools:
            self.tool_trajectories[" -> ".join(tools)] += 1
        provider = episode.get("provider")
        if isinstance(provider, dict):
            key = f"{provider.get('name') or 'unknown'}:{provider.get('model') or 'unknown'}"
            self.providers[key] += 1

    def note_compact(self, record: JsonDict) -> None:
        assistant = _last_assistant_content(record.get("messages"))
        spec = _load_json_object(assistant)
        if not isinstance(spec, dict):
            return
        for job in spec.get("jobs") or []:
            if isinstance(job, dict) and job.get("kind"):
                self.job_kinds[str(job["kind"])] += 1

    def manifest(self) -> JsonDict:
        distinct = len(self.query_skeletons)
        total = sum(self.query_skeletons.values())
        return {
            "episodes_seen": self.episodes_seen,
            "session_chains_seen": self.session_chains_seen,
            "multi_turn_chains": self.multi_turn_chains,
            "multi_turn_chain_share": (
                round(self.multi_turn_chains / self.session_chains_seen, 4)
                if self.session_chains_seen
                else 0.0
            ),
            "written": dict(sorted(self.written.items())),
            "skipped": dict(sorted(self.skipped.items())),
            "providers": dict(sorted(self.providers.items())),
            "diversity": {
                "query_count": total,
                "distinct_query_skeletons": distinct,
                "distinct_query_skeleton_ratio": (
                    round(distinct / total, 4) if total else 0.0
                ),
                "top_query_skeletons": self.query_skeletons.most_common(10),
                "tool_trajectories": self.tool_trajectories.most_common(20),
                "job_kinds": dict(sorted(self.job_kinds.items())),
            },
        }


def _legacy_counts(manifest: JsonDict) -> JsonDict:
    written = manifest.get("written") or {}
    skipped = manifest.get("skipped") or {}
    compact_skipped = sum(
        count
        for key, count in skipped.items()
        if str(key).startswith("compact_spec:")
    )
    return {
        "episodes_seen": manifest.get("episodes_seen", 0),
        "tool_loop_written": written.get("tool_loop", 0),
        "compact_spec_written": written.get("compact_spec", 0),
        "compact_spec_skipped": compact_skipped,
        "reasoning_synthesis_written": written.get("reasoning_synthesis", 0),
        "wrong_route_contrast_written": written.get(
            "wrong_route_contrast", 0
        ),
        "command_answer_written": written.get("command_answer", 0),
        "project_yaml_written": written.get("project_yaml", 0),
        "repair_pairs_written": written.get("repair_pairs", 0),
        "terminal_state_written": written.get("terminal_state", 0),
    }


def _training_dir(value: str | None) -> Path:
    if value:
        return Path(value).expanduser()
    from chemsmart.agent.training_log import _default_training_dir

    return _default_training_dir()


def _episode_paths(
    paths: list[str],
    training_dir: Path,
    *,
    include_runs: bool = False,
) -> list[Path]:
    if paths:
        return [Path(path).expanduser() for path in paths]
    pattern = str(training_dir / "episodes" / "*.jsonl")
    found = [Path(path) for path in sorted(glob.glob(pattern))]
    legacy = sorted(glob.glob(str(training_dir / "episodes-*.jsonl")))
    found.extend(Path(path) for path in legacy)
    if include_runs:
        for run_pattern in (
            "runs/*/episodes/*.jsonl",
            "runs/*/episodes-*.jsonl",
        ):
            found.extend(sorted(training_dir.glob(run_pattern)))
    return found


def _dedup_episodes(
    episodes: Iterable[JsonDict],
) -> Iterable[JsonDict]:
    """Keep only the last episode per ``(session_id, turn)``.

    A turn paused on ask_user is snapshotted (``paused: true``) and written
    again in full when it resumes and completes; exporting both would emit a
    near-duplicate prefix record. Later lines win (append order is
    chronological), so a resumed completion supersedes its pause snapshot
    while a never-resumed pause survives as the only record of that turn.
    """

    buffered: dict[tuple[Any, ...], JsonDict] = {}
    order: list[tuple[Any, ...]] = []
    for index, episode in enumerate(episodes):
        session_id = str(episode.get("session_id") or "")
        turn = episode.get("turn")
        key: tuple[Any, ...] = (session_id, turn)
        if not session_id or turn is None:
            key = ("", None, index)
        if key not in buffered:
            order.append(key)
        buffered[key] = episode
    for key in order:
        yield buffered[key]


def _iter_episodes(paths: Iterable[Path]) -> Iterable[JsonDict]:
    for path in paths:
        try:
            handle = path.open(encoding="utf-8")
        except OSError:
            continue
        with handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                try:
                    value = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if isinstance(value, dict):
                    yield value


def _session_chains(
    episodes: Iterable[JsonDict],
) -> Iterable[list[JsonDict]]:
    """Group deduplicated turn snapshots into deterministic session chains."""

    grouped: dict[str, list[tuple[int, JsonDict]]] = {}
    order: list[str] = []
    for index, episode in enumerate(episodes):
        session_id = str(episode.get("session_id") or "")
        key = session_id or f"__anonymous__:{index}"
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append((index, episode))

    for key in order:
        entries = grouped[key]
        entries.sort(key=lambda item: (_turn_sort_key(item[1]), item[0]))
        yield [episode for _, episode in entries]


def _turn_sort_key(episode: JsonDict) -> tuple[int, str]:
    turn = episode.get("turn")
    try:
        return int(turn), str(turn)
    except (TypeError, ValueError):
        return 10**9, str(turn or "")


def _merge_chain_messages(chain: list[JsonDict]) -> list[JsonDict]:
    """Merge per-turn traces while removing repeated history prefixes."""

    merged: list[JsonDict] = []
    system_message: JsonDict | None = None
    for episode in chain:
        episode_messages = [
            _public_message(message)
            for message in episode.get("messages") or []
            if isinstance(message, dict)
            and isinstance(message.get("role"), str)
        ]
        for message in episode_messages:
            if message.get("role") == "system" and system_message is None:
                system_message = message
        messages = [
            message
            for message in episode_messages
            if message.get("role") != "system"
        ]
        if not messages:
            continue
        overlap = 0
        maximum = min(len(merged), len(messages))
        for size in range(maximum, 0, -1):
            if merged[-size:] == messages[:size]:
                overlap = size
                break
        merged.extend(messages[overlap:])

    final_answer = str(chain[-1].get("final_answer") or "").strip()
    if final_answer and _last_assistant_content(merged).strip() != final_answer:
        merged.append({"role": "assistant", "content": final_answer})
    return ([system_message] if system_message is not None else []) + merged


def _user_message_count(messages: Any) -> int:
    if not isinstance(messages, list):
        return 0
    return sum(
        1
        for message in messages
        if isinstance(message, dict) and message.get("role") == "user"
    )


def _tool_loop_record(
    episode: JsonDict,
    *,
    training_dir: Path,
    include_rejected: bool,
) -> tuple[JsonDict | None, str]:
    return _tool_loop_chain_record(
        [episode],
        training_dir=training_dir,
        include_rejected=include_rejected,
    )


def _tool_loop_chain_record(
    chain: list[JsonDict],
    *,
    training_dir: Path,
    include_rejected: bool,
) -> tuple[JsonDict | None, str]:
    if not chain:
        return None, "empty_chain"
    reason = _chain_terminal_rejection_reason(
        chain, include_rejected=include_rejected
    )
    if reason:
        return None, reason

    messages = _merge_chain_messages(chain)
    if not messages:
        return None, "missing_messages"
    prompt_sha = next(
        (
            str(episode.get("system_prompt_sha") or "")
            for episode in reversed(chain)
            if episode.get("system_prompt_sha")
        ),
        "",
    )
    system_prompt = _load_prompt(training_dir, prompt_sha)
    exported_messages: list[JsonDict] = []
    if system_prompt:
        exported_messages.append({"role": "system", "content": system_prompt})
    exported_messages.extend(messages)
    if not exported_messages:
        return None, "empty_messages"
    system_messages = [
        message
        for message in exported_messages
        if message.get("role") == "system"
    ]
    if len(system_messages) == 0:
        return None, "missing_system_prompt"
    if len(system_messages) > 1:
        return None, "duplicate_system_prompt"

    return {
        "messages": exported_messages,
        "meta": _chain_meta(chain, family="tool_loop"),
    }, "ok"


def _public_message(message: JsonDict) -> JsonDict:
    """Strip provider transport and hidden-reasoning fields before export."""

    return {
        key: message[key]
        for key in _PUBLIC_MESSAGE_KEYS
        if key in message
    }


def _chain_terminal_rejection_reason(
    chain: list[JsonDict],
    *,
    include_rejected: bool,
) -> str:
    if any(_contains_secret(episode) for episode in chain):
        return "secret_detected"
    terminal = chain[-1]
    if terminal.get("paused"):
        return "unresolved_pause"
    outcome = terminal.get("outcome") or {}
    if isinstance(outcome, dict):
        if outcome.get("denied"):
            return "terminal_user_denied"
        rc = outcome.get("execute_rc")
        if isinstance(rc, int) and rc != 0:
            return "terminal_execute_failed"

    terminal_events = _tool_events(terminal)
    if terminal_events and _tool_event_failed(terminal_events[-1]):
        return "terminal_tool_error"

    synthesis = terminal.get("synthesis")
    if isinstance(synthesis, dict):
        status = str(synthesis.get("status") or "")
        command = str(synthesis.get("command") or "").strip()
        gate = _command_gate(terminal)
        if status == "needs_clarification":
            return "unresolved_clarification"
        if _is_submission_command(command):
            terminal_reason = _terminal_state_rejection_reason(terminal)
            if terminal_reason and not include_rejected:
                return terminal_reason
        if command:
            if not command.startswith("chemsmart "):
                return "malformed_command"
            violation = _canonical_violation(command)
            if violation:
                return f"canonical_{violation}"
            if gate not in _COMMAND_POSITIVE_GATES and not include_rejected:
                return "terminal_semantic_reject"
            if (
                gate in _COMMAND_POSITIVE_GATES
                and not _episode_has_generated_input_evidence(terminal)
            ):
                return "missing_generated_input_evidence"
        elif status == "ready":
            return "ready_without_command"

    messages = _merge_chain_messages(chain)
    if (
        not messages
        or messages[-1].get("role") != "assistant"
        or not str(messages[-1].get("content") or "").strip()
    ):
        return "terminal_no_response"
    final_response = str(terminal.get("final_answer") or "").strip()
    if not final_response:
        final_response = _last_assistant_content(messages).strip()
    if not final_response:
        return "terminal_no_response"
    return ""


def _tool_event_failed(event: JsonDict) -> bool:
    if str(event.get("status") or "").lower() == "error":
        return True
    result = event.get("result_summary")
    if not isinstance(result, dict):
        return False
    return result.get("ok") is False or bool(result.get("error"))


def _tool_loop_review_record(
    chain: list[JsonDict],
    *,
    reason: str,
) -> JsonDict | None:
    if not chain or any(_contains_secret(episode) for episode in chain):
        return None
    return {
        "messages": _merge_chain_messages(chain),
        "meta": {
            **_chain_meta(chain, family="tool_loop_review"),
            "skip_reason": reason,
        },
    }


def _chain_meta(chain: list[JsonDict], *, family: str) -> JsonDict:
    terminal = chain[-1]
    merged_messages = _merge_chain_messages(chain)
    tools = [
        tool
        for episode in chain
        for tool in _invoked_tools(episode)
    ]
    return {
        **_meta(terminal, family=family),
        "source_episode_count": len(chain),
        "turns": [episode.get("turn") for episode in chain],
        "user_turn_count": _user_message_count(merged_messages),
        "tool_trajectory": tools,
        "recovered_reject": any(
            _command_gate(episode) == "reject" for episode in chain[:-1]
        )
        and _command_gate(terminal) in _COMMAND_POSITIVE_GATES,
        "terminal_state": terminal.get("terminal_state"),
    }


def _compact_spec_record(
    episode: JsonDict,
    *,
    allow_frontier_distill: bool,
) -> tuple[JsonDict | None, str]:
    if _contains_secret(episode):
        return None, "secret_detected"
    if _command_gate(episode) not in _COMMAND_POSITIVE_GATES:
        return None, "semantic_gate_not_positive"

    synthesis = episode.get("synthesis")
    if not isinstance(synthesis, dict):
        return None, "missing_synthesis"

    raw = synthesis.get("raw_response")
    spec = _load_json_object(raw)
    default_project = _episode_project(episode) or _command_project(synthesis)
    if _looks_like_compact_spec(spec):
        record, reason = _compact_record_from_spec(
            episode, spec, default_project=default_project
        )
        if record is not None:
            return record, reason

    command = str(synthesis.get("command") or "").strip()
    spec, reason = _compact_spec_from_command(command)
    if spec is not None:
        record, record_reason = _compact_record_from_spec(
            episode, spec, default_project=default_project
        )
        if record is not None:
            return record, record_reason
        reason = record_reason

    if allow_frontier_distill:
        spec, reason = _frontier_distill_compact_spec(episode, command)
        if spec is not None:
            return _compact_record_from_spec(
                episode, spec, default_project=default_project
            )
    return None, reason


def _compact_record_from_spec(
    episode: JsonDict,
    spec: JsonDict | None,
    *,
    default_project: str | None,
) -> tuple[JsonDict | None, str]:
    if not _looks_like_compact_spec(spec):
        return None, "not_compact_spec"
    from chemsmart.agent.local.generator import load_system_prompt
    from chemsmart.agent.v8_adapter import adapt

    user = _first_user_message(episode.get("messages"))
    if not user:
        return None, "missing_user"

    adapted = adapt(spec, validate=True, default_project=default_project)
    if not adapted.get("valid") or not adapted.get("commands"):
        return None, "adapter_invalid"
    for rendered in adapted.get("commands") or []:
        violation = _canonical_violation(str(rendered))
        if violation:
            return None, f"canonical_{violation}"
        semantic_violation = _semantic_pair_violation(user, str(rendered))
        if semantic_violation:
            return None, f"semantic_{semantic_violation}"

    return {
        "messages": [
            {"role": "system", "content": load_system_prompt()},
            {"role": "user", "content": user},
            {
                "role": "assistant",
                "content": json.dumps(
                    adapted.get("spec"),
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
            },
        ],
        "meta": {
            **_meta(episode, family="compact_spec"),
            "commands": adapted.get("commands") or [],
            "source": (
                "raw_response"
                if _looks_like_compact_spec(
                    _load_json_object(
                        (episode.get("synthesis") or {}).get("raw_response")
                    )
                )
                else "deterministic_command_reverse"
            ),
        },
    }, "ok"


def _compact_spec_from_command(command: str) -> tuple[JsonDict | None, str]:
    if not command:
        return None, "missing_command"
    from chemsmart.agent.model_command_parser import parse_model_command

    parsed = parse_model_command(command)
    if parsed.parse_error:
        return None, "command_parse_error"
    if parsed.program not in {"gaussian", "orca"} or not parsed.job:
        return None, "unsupported_command"
    if parsed.functional or parsed.basis or parsed.ab_initio:
        return None, "runtime_owned_method_flags"
    if parsed.solvent_model or parsed.solvent_id:
        return None, "runtime_owned_solvent_flags"

    kind_job = (
        "tddft"
        if parsed.program == "gaussian" and parsed.job == "td"
        else parsed.job
    )
    job: JsonDict = {
        "id": 1,
        "kind": f"{parsed.program}.{kind_job}",
        "file": parsed.filename,
        "charge": _int_or_default(parsed.charge, 0),
        "mult": _int_or_default(parsed.multiplicity, 1),
    }
    if parsed.entrypoint == "sub":
        job["execution"] = "submit"
    if parsed.server:
        job["server"] = parsed.server
    if parsed.label:
        job["label"] = parsed.label

    settings: JsonDict = {}
    if parsed.route_parameters:
        settings["additional_route_parameters"] = parsed.route_parameters
    if parsed.opt_options:
        settings["additional_opt_options_in_route"] = parsed.opt_options
    for key, value in (parsed.structural_options or {}).items():
        if value not in (None, "", [], {}):
            settings[key] = value
    if settings:
        job["settings"] = settings
    if not job.get("file"):
        return None, "missing_file"
    return {"intent": "workflow", "jobs": [job]}, "ok"


# Compact synthesis contract for the reasoning-synthesis family. Kept short
# (no full CLI-schema dump) so training sequences stay small; the local model
# learns the CLI grammar from the reasoning->command pairs. Serving must send
# this same compact prompt to the local model (path R).
_REASONING_SYNTH_SYSTEM = (
    "You are the ChemSmart CLI command synthesizer. Write a short public, "
    "auditable decision trace inside a <think>...</think> block, then return "
    "ONE JSON object:\n"
    '{"status":"ready"|"needs_clarification"|"infeasible",'
    '"command":"chemsmart ...","explanation":"...",'
    '"confidence":"low"|"medium"|"high",'
    '"missing_info":[...],"alternatives":[...]}\n'
    "The command must start with `chemsmart run` or `chemsmart sub`, name the "
    "engine (gaussian|orca) and one job subcommand, and carry only structural "
    "options (-f file, -c charge, -m mult, -p project, plus job-specific "
    "structural flags). Method/basis/solvent/dispersion are runtime-owned "
    "(project YAML) and must NOT appear as flags. The trace may contain only "
    "observable request, command, project, and deterministic-gate facts; "
    "never include private chain-of-thought. Emit no text outside the <think> "
    "block and the JSON."
)


def _reasoning_source(episode: JsonDict) -> tuple[str, str]:
    """Return public decision trace or a quarantined legacy provenance."""

    synthesis = episode.get("synthesis")
    if isinstance(synthesis, dict):
        reasoning = synthesis.get("reasoning")
        if isinstance(reasoning, str) and reasoning.strip():
            provenance = str(
                synthesis.get("reasoning_provenance") or ""
            ).strip()
            if provenance == "public_decision_trace":
                return reasoning.strip(), provenance
            return reasoning.strip(), "untrusted_synthesis"
    messages = episode.get("messages")
    if isinstance(messages, list):
        for message in reversed(messages):
            if not isinstance(message, dict):
                continue
            if message.get("role") != "assistant":
                continue
            rc = message.get("reasoning_content")
            if isinstance(rc, str) and rc.strip():
                return rc.strip(), "assistant_message"
    return "", "none"


def _has_generated_input_evidence(synthesis: JsonDict) -> bool:
    evidence = synthesis.get("generated_input_evidence")
    return isinstance(evidence, list) and any(
        isinstance(item, dict) for item in evidence
    )


def _episode_has_generated_input_evidence(episode: JsonDict) -> bool:
    synthesis = episode.get("synthesis")
    if isinstance(synthesis, dict) and _has_generated_input_evidence(synthesis):
        return True
    for event in _tool_events(episode):
        evidence = event.get("generated_input_evidence")
        if isinstance(evidence, list) and any(
            isinstance(item, dict) for item in evidence
        ):
            return True
        result = event.get("result_summary")
        semantic = result.get("semantic") if isinstance(result, dict) else None
        generated = (
            semantic.get("generated_inputs")
            if isinstance(semantic, dict)
            else None
        )
        if isinstance(generated, list) and any(
            isinstance(item, dict) for item in generated
        ):
            return True
    return False


def _reasoning_synthesis_record(
    episode: JsonDict,
) -> tuple[JsonDict | None, str]:
    """Emit a public-trace then synthesize target for local training."""

    if _contains_secret(episode):
        return None, "secret_detected"
    if _command_gate(episode) not in _COMMAND_POSITIVE_GATES:
        return None, "semantic_gate_not_positive"
    synthesis = episode.get("synthesis")
    if not isinstance(synthesis, dict):
        return None, "missing_synthesis"
    command = str(synthesis.get("command") or "").strip()
    if not command:
        return None, "missing_command"
    if not command.startswith("chemsmart "):
        return None, "malformed_command"
    violation = _canonical_violation(command)
    if violation:
        return None, f"canonical_{violation}"
    user = _first_user_message(episode.get("messages"))
    if not user:
        return None, "missing_user"
    semantic_violation = _semantic_pair_violation(user, command)
    if semantic_violation:
        return None, f"semantic_{semantic_violation}"
    if not _has_generated_input_evidence(synthesis):
        return None, "missing_generated_input_evidence"
    reasoning, source = _reasoning_source(episode)
    if not reasoning:
        return None, "no_reasoning_trace"
    if source != "public_decision_trace":
        return None, f"untrusted_reasoning:{source}"

    target = {
        "status": synthesis.get("status") or "ready",
        "command": command,
        "explanation": str(synthesis.get("explanation") or "").strip()
        or "Synthesized ChemSmart command.",
        "confidence": "high",
        "missing_info": synthesis.get("missing_info") or [],
        "alternatives": synthesis.get("alternatives") or [],
    }
    assistant = f"<think>{reasoning}</think>\n" + json.dumps(
        target, ensure_ascii=False
    )
    return {
        "messages": [
            {"role": "system", "content": _REASONING_SYNTH_SYSTEM},
            {"role": "user", "content": user},
            {"role": "assistant", "content": assistant},
        ],
        "meta": {
            **_meta(episode, family="reasoning_synthesis"),
            "command": command,
            "reasoning_chars": len(reasoning),
            "reasoning_provenance": source,
        },
    }, "ok"


def _reasoning_synthesis_review_record(
    episode: JsonDict,
    *,
    skip_reason: str,
) -> JsonDict | None:
    """Non-SFT branch for skipped reasoning-synthesis candidates.

    This keeps the old accumulated data inspectable without letting ambiguous
    or post-tool reasoning leak into the positive local-model training target.
    """

    if _contains_secret(episode):
        return None
    synthesis = episode.get("synthesis")
    reasoning, source = _reasoning_source(episode)
    if not isinstance(synthesis, dict) and not reasoning:
        return None
    command = ""
    semantic_verdict = ""
    generated_input_count = 0
    if isinstance(synthesis, dict):
        command = str(synthesis.get("command") or "").strip()
        semantic_verdict = str(synthesis.get("semantic_verdict") or "")
        evidence = synthesis.get("generated_input_evidence")
        if isinstance(evidence, list):
            generated_input_count = len(
                [item for item in evidence if isinstance(item, dict)]
            )
    return {
        "meta": {
            **_meta(episode, family="reasoning_synthesis_review"),
            "skip_reason": skip_reason,
            "command": command,
            "semantic_verdict": semantic_verdict,
            "reasoning_source": source,
            "reasoning_chars": len(reasoning),
            "generated_input_count": generated_input_count,
        },
        "user": _first_user_message(episode.get("messages")) or "",
    }


def _frontier_distill_compact_spec(
    episode: JsonDict,
    command: str,
) -> tuple[JsonDict | None, str]:
    # Kept optional and off by default: export should be deterministic unless a
    # caller explicitly permits outbound provider calls.
    if not command:
        return None, "frontier_distill_missing_command"
    try:
        from chemsmart.agent.providers import get_provider

        provider = get_provider()
        raw = provider.chat(
            [
                {
                    "role": "system",
                    "content": (
                        "Return ONLY a ChemSmart compact SPEC JSON object. "
                        "Do not include project, functional, basis, solvent, "
                        "or other runtime-owned fields."
                    ),
                },
                {
                    "role": "user",
                    "content": json.dumps(
                        {
                            "user_query": _first_user_message(
                                episode.get("messages")
                            ),
                            "command": command,
                            "semantic": (episode.get("synthesis") or {}).get(
                                "semantic"
                            ),
                        },
                        ensure_ascii=False,
                    ),
                },
            ]
        )
    except Exception:
        return None, "frontier_distill_failed"
    spec = _load_json_object(_extract_text(raw))
    return (
        (spec, "ok")
        if _looks_like_compact_spec(spec)
        else (None, "frontier_distill_not_compact")
    )


def _command_answer_record(episode: JsonDict) -> tuple[JsonDict | None, str]:
    synthesis = episode.get("synthesis")
    if not isinstance(synthesis, dict):
        return None, "missing_synthesis"
    command = str(synthesis.get("command") or "").strip()
    if not command:
        return None, "missing_command"
    if _command_gate(episode) not in _COMMAND_POSITIVE_GATES:
        return None, "semantic_gate_not_positive"
    violation = _canonical_violation(command)
    if violation:
        return None, f"canonical_{violation}"
    user = _first_user_message(episode.get("messages"))
    if not user:
        return None, "missing_user"
    semantic_violation = _semantic_pair_violation(user, command)
    if semantic_violation:
        return None, f"semantic_{semantic_violation}"
    from chemsmart.agent.model_command_parser import (
        format_parsed_model_command,
    )

    parsed = _parse_command_for_episode(command, episode)
    workspace = _episode_workspace_label(episode)
    if workspace:
        parsed = replace(parsed, workspace=workspace)
    deterministic = _mask_current_paths(format_parsed_model_command(parsed))
    parsed_dict = _mask_current_paths(parsed.to_dict())
    final_answer = str(episode.get("final_answer") or "").strip()
    assistant = final_answer or deterministic
    return {
        "messages": [
            {
                "role": "system",
                "content": (
                    "Explain a grounded ChemSmart command using only the "
                    "deterministic parser facts."
                ),
            },
            {"role": "user", "content": user},
            {"role": "assistant", "content": assistant},
        ],
        "meta": {
            **_meta(episode, family="command_answer"),
            "command": command,
            "parsed_command": parsed_dict,
            "deterministic_answer": deterministic,
        },
    }, "ok"


def _parse_command_for_episode(command: str, episode: JsonDict) -> Any:
    from chemsmart.agent.model_command_parser import parse_model_command

    yaml_payload = _workspace_yaml_payload(episode)
    if yaml_payload is None:
        return parse_model_command(command)

    with tempfile.TemporaryDirectory(prefix="chemsmart-export-parse-") as tmp:
        root = Path(tmp)
        folder = root / ".chemsmart" / yaml_payload["program"]
        folder.mkdir(parents=True, exist_ok=True)
        (folder / f"{yaml_payload['project_name']}.yaml").write_text(
            yaml_payload["yaml_text"],
            encoding="utf-8",
        )
        original_cwd = os.getcwd()
        os.chdir(root)
        try:
            return parse_model_command(command, cwd=str(root))
        finally:
            os.chdir(original_cwd)


def _workspace_yaml_payload(episode: JsonDict) -> JsonDict | None:
    workspace = episode.get("workspace")
    workspace_project = ""
    workspace_program = ""
    if isinstance(workspace, dict):
        workspace_project = str(workspace.get("project") or "").strip()
        workspace_program = str(workspace.get("program") or "").strip()
    for event in reversed(_tool_events(episode)):
        result = event.get("result_summary")
        if not isinstance(result, dict):
            continue
        yaml_text = result.get("yaml_text")
        if not isinstance(yaml_text, str) or not yaml_text.strip():
            continue
        project_name = str(
            result.get("project_name") or workspace_project or "candidate"
        ).strip()
        program = str(
            result.get("program") or workspace_program or "gaussian"
        ).strip()
        if program not in {"gaussian", "orca"}:
            continue
        return {
            "program": program,
            "project_name": Path(project_name).stem or "candidate",
            "yaml_text": yaml_text,
        }
    return None


def _project_yaml_record(episode: JsonDict) -> tuple[JsonDict | None, str]:
    events = _tool_events(episode)
    if not _has_project_yaml_build_event(events):
        return None, "not_project_yaml_workflow"
    yaml_text = ""
    verdict = ""
    for event in events:
        result = event.get("result_summary")
        if not isinstance(result, dict):
            continue
        if isinstance(result.get("yaml_text"), str):
            yaml_text = result["yaml_text"]
        if result.get("verdict"):
            verdict = str(result["verdict"])
        if result.get("validation_verdict"):
            verdict = str(result["validation_verdict"])
    if not yaml_text:
        return None, "missing_yaml"
    if verdict and verdict not in {"ok", "warn"}:
        return None, "validation_not_positive"
    user = _first_user_message(episode.get("messages"))
    if not user:
        return None, "missing_user"
    final_answer = str(episode.get("final_answer") or "").strip()
    assistant = final_answer or yaml_text
    return {
        "messages": [
            {
                "role": "system",
                "content": (
                    "Build and validate a ChemSmart project YAML from a "
                    "reported computational protocol."
                ),
            },
            {"role": "user", "content": user},
            {"role": "assistant", "content": assistant},
        ],
        "meta": {
            **_meta(episode, family="project_yaml"),
            "yaml_text": yaml_text,
            "validation_verdict": verdict or "unknown",
            "tool_trajectory": _invoked_tools(episode),
        },
    }, "ok"


def _has_project_yaml_build_event(events: list[JsonDict]) -> bool:
    """Return whether this is a project-YAML authoring trajectory.

    ``read_project_yaml`` can appear in ordinary command-synthesis turns to
    inspect the active workspace. Those turns may carry ``yaml_text`` evidence,
    but they are not project-YAML SFT examples.
    """

    authoring_tools = {
        "extract_project_protocol",
        "render_project_yaml",
        "validate_project_yaml",
        "critic_project_yaml",
        "write_project_yaml",
        "update_project_yaml",
    }
    return any(
        str(event.get("tool") or "") in authoring_tools for event in events
    )


def _repair_pair_record(episode: JsonDict) -> tuple[JsonDict | None, str]:
    events = _tool_events(episode)
    rejected = ""
    repaired = ""
    reject_semantic: Any = None
    repair_semantic: Any = None
    for event in events:
        result = event.get("result_summary")
        if not isinstance(result, dict):
            continue
        if event.get("tool") == "synthesize_command":
            semantic = result.get("semantic")
            if (
                isinstance(semantic, dict)
                and semantic.get("verdict") == "reject"
            ):
                rejected = str(result.get("command") or "")
                reject_semantic = semantic
        if event.get("tool") == "repair_command":
            repaired = str(result.get("command") or "")
            repair_semantic = result.get("semantic")
            # repair_command results carry the broken input themselves:
            # original_command + repaired flag. That mints a contrast pair
            # even when no synthesize_command reject event exists in the
            # same episode (e.g. "fix this command: ..." requests).
            if not rejected:
                original = str(result.get("original_command") or "")
                if (
                    original
                    and result.get("repaired")
                    and original != (repaired)
                ):
                    rejected = original
                    reject_semantic = {"verdict": "reject", "issues": []}
    if not rejected:
        return None, "missing_reject"
    if not repaired:
        return None, "missing_repair"
    if not isinstance(repair_semantic, dict) or repair_semantic.get(
        "verdict"
    ) not in {"ok", "warn"}:
        return None, "repair_not_positive"
    violation = _canonical_violation(repaired)
    if violation:
        return None, f"canonical_{violation}"
    return {
        "prompt": _first_user_message(episode.get("messages")),
        "rejected": {
            "command": rejected,
            "semantic": reject_semantic,
        },
        "chosen": {
            "command": repaired,
            "semantic": repair_semantic,
        },
        "meta": _meta(episode, family="repair_pair"),
    }, "ok"


def _cross_turn_repair_pair_record(
    chain: list[JsonDict],
) -> tuple[JsonDict | None, str]:
    """Pair a rejected command with a later gated command in one session."""

    rejected: tuple[JsonDict, str, JsonDict] | None = None
    for episode in chain:
        synthesis = episode.get("synthesis")
        if not isinstance(synthesis, dict):
            continue
        command = str(synthesis.get("command") or "").strip()
        gate = _command_gate(episode)
        if command and gate == "reject":
            rejected = (
                episode,
                command,
                {
                    "verdict": "reject",
                    "failed_rule_ids": synthesis.get("failed_rule_ids") or [],
                },
            )
            continue
        if (
            rejected is None
            or not command
            or gate not in _COMMAND_POSITIVE_GATES
            or command == rejected[1]
        ):
            continue
        if _canonical_violation(command):
            continue
        if not _episode_has_generated_input_evidence(episode):
            continue
        rejected_episode, rejected_command, rejected_semantic = rejected
        prompts = [
            str(message.get("content") or "").strip()
            for message in _merge_chain_messages(chain)
            if message.get("role") == "user"
            and str(message.get("content") or "").strip()
        ]
        return {
            "prompt": "\n".join(
                f"Turn {index}: {prompt}"
                for index, prompt in enumerate(prompts, start=1)
            ),
            "rejected": {
                "command": rejected_command,
                "semantic": rejected_semantic,
            },
            "chosen": {
                "command": command,
                "semantic": {
                    "verdict": gate,
                    "failed_rule_ids": synthesis.get("failed_rule_ids") or [],
                    "generated_inputs": synthesis.get(
                        "generated_input_evidence"
                    )
                    or [],
                },
            },
            "meta": {
                **_chain_meta(chain, family="repair_pair"),
                "source": "cross_turn",
                "rejected_turn": rejected_episode.get("turn"),
                "chosen_turn": episode.get("turn"),
            },
        }, "ok"
    return None, "no_cross_turn_repair"


def _repair_pair_key(record: JsonDict) -> tuple[str, str]:
    rejected = record.get("rejected") or {}
    chosen = record.get("chosen") or {}
    return (
        str(rejected.get("command") or ""),
        str(chosen.get("command") or ""),
    )


_AUTHORING_TOOLS = {
    "extract_project_protocol",
    "render_project_yaml",
    "validate_project_yaml",
    "critic_project_yaml",
    "write_project_yaml",
    "update_project_yaml",
}


def _wrong_route_contrast_records(episodes: list[JsonDict]) -> list[JsonDict]:
    """Pair command requests routed to project-YAML tools with corrections.

    This is a contrast family, not positive SFT. It teaches the router that a
    request can mention loaded project settings without asking to author YAML.
    """

    wrong_by_key: dict[str, list[JsonDict]] = {}
    corrected_by_key: dict[str, list[JsonDict]] = {}
    for episode in episodes:
        if _contains_secret(episode):
            continue
        wrong_key = _wrong_route_key(episode)
        if wrong_key:
            wrong_by_key.setdefault(wrong_key, []).append(episode)
            continue
        corrected_key = _corrected_route_key(episode)
        if corrected_key:
            corrected_by_key.setdefault(corrected_key, []).append(episode)

    records: list[JsonDict] = []
    used_corrected: set[tuple[Any, Any]] = set()
    seen_pairs: set[tuple[str, str]] = set()
    for key, wrong_episodes in sorted(wrong_by_key.items()):
        corrected = corrected_by_key.get(key) or []
        if not corrected:
            continue
        for wrong in wrong_episodes:
            # Prefer the correction from the same session.  A global fallback
            # is still useful when the teacher ended the YAML route before a
            # later command-only example, but cross-session pairing can hide a
            # valid repair with a different molecule, state, or region.
            same_session = [
                candidate
                for candidate in corrected
                if candidate.get("session_id") == wrong.get("session_id")
            ]
            candidates = same_session or corrected
            chosen = None
            for candidate in candidates:
                candidate_id = (
                    candidate.get("session_id"),
                    candidate.get("turn"),
                )
                if candidate_id not in used_corrected:
                    chosen = candidate
                    used_corrected.add(candidate_id)
                    break
            if chosen is None:
                chosen = corrected[0]
            record = _wrong_route_contrast_record(key, wrong, chosen)
            if record is not None:
                pair_key = (
                    _query_skeleton(str(record.get("prompt") or "")),
                    str((record.get("chosen") or {}).get("command") or ""),
                )
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                records.append(record)
    return records


def _wrong_route_key(episode: JsonDict) -> str:
    tools = set(_invoked_tools(episode))
    if not tools.intersection(_AUTHORING_TOOLS):
        return ""
    synthesis = episode.get("synthesis")
    if isinstance(synthesis, dict):
        command = str(synthesis.get("command") or "").strip()
        if command and _command_gate(episode) in _COMMAND_POSITIVE_GATES:
            return ""
    user = _first_user_message(episode.get("messages"))
    if not _looks_like_command_request(user):
        return ""
    kind = _infer_kind_from_text(user)
    program = _infer_program_from_text(user)
    return f"{program}:{kind}" if program and kind else ""


def _corrected_route_key(episode: JsonDict) -> str:
    if _command_gate(episode) not in _COMMAND_POSITIVE_GATES:
        return ""
    synthesis = episode.get("synthesis")
    if not isinstance(synthesis, dict):
        return ""
    command = str(synthesis.get("command") or "").strip()
    if not command or _canonical_violation(command):
        return ""
    program = _infer_program_from_command(command)
    kind = _infer_kind_from_command(command)
    return f"{program}:{kind}" if program and kind else ""


def _wrong_route_contrast_record(
    key: str, wrong: JsonDict, corrected: JsonDict
) -> JsonDict | None:
    user = _first_user_message(wrong.get("messages"))
    corrected_synthesis = corrected.get("synthesis")
    if not user or not isinstance(corrected_synthesis, dict):
        return None
    command = str(corrected_synthesis.get("command") or "").strip()
    if not command:
        return None
    wrong_tools = _invoked_tools(wrong)
    corrected_tools = _invoked_tools(corrected)
    return {
        "prompt": user,
        "rejected": {
            "route": "project_yaml",
            "tools": wrong_tools,
            "reason": "job command request was routed to project-YAML authoring",
            "synthesis_status": (wrong.get("synthesis") or {}).get("status")
            if isinstance(wrong.get("synthesis"), dict)
            else None,
        },
        "chosen": {
            "route": "synthesize_command",
            "tools": corrected_tools,
            "command": command,
            "semantic_gate": _command_gate(corrected),
            "reasoning_source": _reasoning_source(corrected)[1],
        },
        "meta": {
            **_meta(corrected, family="wrong_route_contrast"),
            "contrast_key": key,
            "wrong_session_id": wrong.get("session_id"),
            "wrong_turn": wrong.get("turn"),
        },
    }


def _looks_like_command_request(text: str) -> bool:
    lowered = text.lower()
    if not lowered:
        return False
    command_markers = (
        "command",
        "cli",
        "run",
        "submit",
        "set up",
        "prepare",
        "make",
        "create",
        "using the loaded",
        "using the active",
        "current project",
        "active project",
    )
    yaml_markers = ("yaml", "project file", "project yaml")
    if not any(marker in lowered for marker in command_markers):
        return False
    if not any(marker in lowered for marker in yaml_markers):
        return True
    # A project-YAML authoring request is intentionally distinct from a
    # command request. Keep only the ambiguous multi-stage queries that also
    # ask for the actual job/CLI artifact; these are the source examples for
    # YAML-route -> command-route contrast training.
    command_after_yaml = (
        "actual job",
        "actual cli",
        "actual command",
        "tell me how to run",
        "then prepare",
        "then run",
        "prepare the calculation",
    )
    return any(marker in lowered for marker in command_after_yaml)


def _infer_program_from_text(text: str) -> str:
    lowered = text.lower()
    if "orca" in lowered:
        return "orca"
    if "gaussian" in lowered or "gauss" in lowered:
        return "gaussian"
    return ""


def _infer_kind_from_text(text: str) -> str:
    lowered = text.lower()
    patterns = (
        ("neb", ("neb", "neb-ts", "neb-ci")),
        # Gaussian literature often calls the same multiscale route ONIOM;
        # retain that synonym so YAML-authoring -> CLI corrections can form
        # the intended qmmm contrast family.
        ("qmmm", ("qm/mm", "qmmm", "oniom")),
        ("td", ("td-dft", "tddft", "td dft", "excited", "absorption")),
        ("dias", ("dias", "activation strain")),
        ("resp", ("resp", "charge fitting")),
        ("wbi", ("wiberg", "bond indices", "bond order")),
        ("nci", ("nci", "non-covalent", "noncovalent")),
        ("irc", ("irc", "connectivity")),
        ("scan", ("scan", "scants", "relaxed scan")),
        ("modred", ("freeze", "freezing", "frozen", "constraint")),
        ("ts", ("transition state", "ts optimization", "optts")),
        ("crest", ("crest", "conformer")),
        ("sp", ("single point", "energy")),
        ("opt", ("optimize", "optimization", "relax")),
    )
    for kind, markers in patterns:
        if any(marker in lowered for marker in markers):
            return kind
    return ""


def _infer_program_from_command(command: str) -> str:
    tokens = command.split()
    for token in tokens:
        if token in {"gaussian", "orca"}:
            return token
    return ""


def _infer_kind_from_command(command: str) -> str:
    tokens = command.split()
    for token in reversed(tokens):
        if token in {
            "crest",
            "dias",
            "irc",
            "modred",
            "neb",
            "nci",
            "opt",
            "qmmm",
            "resp",
            "scan",
            "sp",
            "td",
            "ts",
            "wbi",
        }:
            return token
    return ""


def _load_prompt(training_dir: Path, digest: str) -> str:
    if not digest:
        return ""
    candidates = [training_dir / "prompts" / f"{digest}.txt"]
    # Parallel harness runs write episodes (and their prompt hashes) into
    # per-model stores under runs/<model>/; resolve there too so exported
    # tool_loop records keep their system prompt.
    candidates.extend(
        sorted(training_dir.glob(f"runs/*/prompts/{digest}.txt"))
    )
    for path in candidates:
        try:
            return path.read_text(encoding="utf-8")
        except OSError:
            continue
    return ""


def _portable_path(path: Path, base: Path) -> str:
    try:
        return os.path.relpath(path.resolve(), base.resolve())
    except OSError:
        return str(path)


def _load_json_object(value: Any) -> JsonDict | None:
    if isinstance(value, dict):
        return value
    if not isinstance(value, str) or not value.strip():
        return None
    try:
        loaded = json.loads(value)
    except json.JSONDecodeError:
        return None
    return loaded if isinstance(loaded, dict) else None


def _looks_like_compact_spec(value: JsonDict | None) -> bool:
    if not isinstance(value, dict):
        return False
    intent = value.get("intent")
    return isinstance(intent, str) and intent in {
        "workflow",
        "advisory",
        "decline",
        "chitchat",
    }


def _first_user_message(messages: Any) -> str:
    if not isinstance(messages, list):
        return ""
    for message in messages:
        if not isinstance(message, dict):
            continue
        if message.get("role") != "user":
            continue
        content = message.get("content")
        if isinstance(content, str) and content.strip():
            return content.strip()
    return ""


def _last_assistant_content(messages: Any) -> str:
    if not isinstance(messages, list):
        return ""
    for message in reversed(messages):
        if not isinstance(message, dict):
            continue
        if message.get("role") == "assistant" and isinstance(
            message.get("content"), str
        ):
            return str(message["content"])
    return ""


def _meta(episode: JsonDict, *, family: str) -> JsonDict:
    outcome = episode.get("outcome") or {}
    synthesis = episode.get("synthesis") or {}
    meta: JsonDict = {
        "family": family,
        "episode_schema_version": episode.get("v"),
        "session_id": episode.get("session_id"),
        "turn": episode.get("turn"),
        "provider": episode.get("provider"),
        "schema_variant": (
            synthesis.get("schema_variant")
            if isinstance(synthesis, dict)
            else None
        ),
        "semantic_gate": _command_gate(episode),
        "workspace": episode.get("workspace") or {},
        "outcome": outcome if isinstance(outcome, dict) else {},
        "terminal_state": episode.get("terminal_state"),
    }
    # Carry sanitized provenance through to the exported record so the split
    # builder can group by scenario_family. Added only when present to keep
    # legacy (provenance-free) exports byte-identical.
    provenance = episode.get("dataset_provenance")
    if isinstance(provenance, dict) and provenance:
        meta["dataset_provenance"] = provenance
    return meta


def _is_submission_command(command: str) -> bool:
    try:
        tokens = shlex.split(command)
    except ValueError:
        return False
    return len(tokens) > 1 and tokens[0] == "chemsmart" and tokens[1] == "sub"


def _terminal_state_rejection_reason(episode: JsonDict) -> str | None:
    state = episode.get("terminal_state")
    if state is None:
        return "missing_terminal_state"
    if validate_terminal_state(state):
        return "terminal_state_assertion_failed"
    if not terminal_state_is_positive(state):
        return "terminal_state_expected_failure"
    return None


def _terminal_state_record(
    episode: JsonDict,
) -> tuple[JsonDict | None, str]:
    state = episode.get("terminal_state")
    if state is None:
        return None, "missing_terminal_state"
    issues = validate_terminal_state(state)
    if issues:
        return None, "invalid:" + ",".join(issues)
    return {
        "terminal_state": state,
        "meta": _meta(episode, family="terminal_state"),
    }, "ok"


def _episode_workspace_label(episode: JsonDict) -> str:
    cwd_masked = episode.get("cwd_masked")
    if isinstance(cwd_masked, str) and cwd_masked.strip():
        return cwd_masked.strip()
    workspace = episode.get("workspace")
    if isinstance(workspace, dict):
        path = workspace.get("path")
        if isinstance(path, str) and path.strip():
            return str(Path(path).parent.parent.parent or ".")
    return ""


def _write_jsonl(handle: Any, record: JsonDict) -> None:
    record = _mask_current_paths(record)
    handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def _contains_secret(value: Any) -> bool:
    try:
        text = json.dumps(value, ensure_ascii=False)
    except TypeError:
        text = repr(value)
    return bool(_SECRET_RE.search(text))


def _mask_current_paths(value: Any) -> Any:
    """Mask exporter-machine paths introduced by deterministic utilities."""

    cwd = str(Path.cwd())
    home = str(Path.home())

    def mask(item: Any) -> Any:
        if isinstance(item, str):
            if cwd and cwd != "/":
                item = item.replace(cwd, ".")
            if home and home != "/":
                item = item.replace(home, "~")
            item = _TEMP_PATH_RE.sub("<runtime-temp-path>", item)
            return item
        if isinstance(item, dict):
            return {key: mask(child) for key, child in item.items()}
        if isinstance(item, list):
            return [mask(child) for child in item]
        return item

    return mask(value)


def _command_gate(episode: JsonDict) -> str | None:
    synthesis = episode.get("synthesis")
    if isinstance(synthesis, dict) and synthesis.get("semantic_verdict"):
        return str(synthesis["semantic_verdict"])
    outcome = episode.get("outcome")
    if isinstance(outcome, dict):
        gate = outcome.get("gate")
        return str(gate) if gate is not None else None
    return None


def _episode_project(episode: JsonDict) -> str:
    workspace = episode.get("workspace")
    if isinstance(workspace, dict) and workspace.get("project"):
        return str(workspace["project"])
    return ""


def _command_project(synthesis: JsonDict) -> str:
    command = str(synthesis.get("command") or "")
    if not command:
        return ""
    try:
        from chemsmart.agent.model_command_parser import parse_model_command

        parsed = parse_model_command(command)
    except Exception:
        return ""
    return parsed.project or ""


def _tool_events(episode: JsonDict) -> list[JsonDict]:
    events = episode.get("tool_events")
    if isinstance(events, list):
        return [event for event in events if isinstance(event, dict)]
    records = []
    for index, message in enumerate(episode.get("messages") or [], start=1):
        if not isinstance(message, dict) or message.get("role") != "tool":
            continue
        payload = _load_json_object(message.get("content"))
        records.append(
            {
                "index": index,
                "tool": "",
                "status": "",
                "result_summary": payload,
            }
        )
    return records


def _invoked_tools(episode: JsonDict) -> list[str]:
    tools = episode.get("invoked_tools")
    if isinstance(tools, list):
        return [str(tool) for tool in tools if str(tool)]
    events = _tool_events(episode)
    names: list[str] = []
    seen: set[str] = set()
    for event in events:
        name = str(event.get("tool") or "")
        if name and name not in seen:
            seen.add(name)
            names.append(name)
    if names:
        return names
    legacy = episode.get("tools")
    if episode.get("v") == 2 and isinstance(legacy, list):
        return [str(tool) for tool in legacy if str(tool)]
    return []


def _query_skeleton(text: str) -> str:
    value = text.lower()
    value = re.sub(
        r"\b\S+\.(?:xyz|com|gjf|inp|out|log|yaml|yml|db)\b", "<file>", value
    )
    value = re.sub(r"\d+(?:\.\d+)?", "#", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _int_or_default(value: Any, default: int) -> int:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return default


def _extract_text(response: Any) -> str:
    if isinstance(response, str):
        return response
    if isinstance(response, dict):
        if isinstance(response.get("raw_plan"), str):
            return str(response["raw_plan"])
        if isinstance(response.get("content"), str):
            return str(response["content"])
        choices = response.get("choices")
        if isinstance(choices, list) and choices:
            message = (
                choices[0].get("message")
                if isinstance(choices[0], dict)
                else None
            )
            if isinstance(message, dict) and isinstance(
                message.get("content"), str
            ):
                return str(message["content"])
    return json.dumps(response, ensure_ascii=False, default=str)


if __name__ == "__main__":
    raise SystemExit(main())
