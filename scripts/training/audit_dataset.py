"""Audit the accumulated agent-training episode store before export/training.

Three question families, matched to what actually breaks small-model SFT:

1. **Pipeline counts** — how many talks (episodes) exist per tool
   trajectory, per export family, per job kind, per provider. Curriculum
   planning for the local model needs these counts, not vibes.
2. **Low-quality / over-frequent data** — flag episodes that must not be
   trained on (gate reject without repair, failed execution, denials,
   secrets, malformed commands) and query skeletons that dominate the
   corpus (the v1→v4 boilerplate trap: near-duplicate phrasings teach
   templates, not the task).
3. **Multi-turn stability** — for the multi-turn tool-using local agent:
   ask_user pause/resume resolution, tool-order sanity (execute without a
   prior synthesize), repeated identical tool calls (loop thrash), and the
   share of genuinely multi-turn conversations.

Usage::

    python scripts/training/audit_dataset.py [episodes.jsonl ...]
        [--training-dir DIR] [--out report.json]
        [--max-per-skeleton N] [--max-skeleton-share 0.05]
        [--write-clean clean.jsonl] [--strict]

With no positional args it audits the configured training store (same
resolution as export_sft.py). ``--write-clean`` emits a deduped, quality-
filtered, skeleton-capped per-turn episode file. It is useful for a separate
per-turn corpus, but is not the authoritative source for session-chain export.
``--strict`` exits 1 when a health threshold fails (CI gate).
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

from chemsmart.agent.harness.terminal_state import (
    terminal_state_is_positive,
    validate_terminal_state,
)

JsonDict = dict[str, Any]

PLANNER_CANONICAL_AGENT_KINDS = (
    "gaussian.sp",
    "gaussian.opt",
    "gaussian.ts",
    "gaussian.freq",
    "gaussian.irc",
    "gaussian.scan",
    "gaussian.modred",
    "gaussian.nci",
    "gaussian.resp",
    "gaussian.tddft",
    "gaussian.dias",
    "gaussian.crest",
    "gaussian.traj",
    "gaussian.wbi",
    "gaussian.qmmm",
    "gaussian.qrc",
    "orca.sp",
    "orca.opt",
    "orca.ts",
    "orca.freq",
    "orca.irc",
    "orca.scan",
    "orca.modred",
    "orca.neb",
    "orca.qmmm",
    "orca.qrc",
)

# Frequency is a project/YAML setting or adapter-routed frequency job. There
# is no standalone ``chemsmart ... freq`` subcommand, so command coverage
# must not reward an impossible model output.
PROJECT_YAML_ONLY_KINDS = ("gaussian.freq", "orca.freq")
COMMAND_CANONICAL_AGENT_KINDS = tuple(
    kind
    for kind in PLANNER_CANONICAL_AGENT_KINDS
    if kind not in PROJECT_YAML_ONLY_KINDS
)

_EXPORT = None


def _export_module():
    """Load sibling export_sft.py for its episode helpers (not a package)."""

    global _EXPORT
    if _EXPORT is None:
        path = Path(__file__).resolve().parent / "export_sft.py"
        spec = importlib.util.spec_from_file_location("export_sft", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _EXPORT = module
    return _EXPORT


# ---------------------------------------------------------------- quality


# A talk is UNTRAINABLE when any of these hold. Keep reasons machine-stable:
# the skill and CI assert on them.
def quality_flags(episode: JsonDict) -> list[str]:
    flags: list[str] = []
    export = _export_module()

    outcome = episode.get("outcome") or {}
    synthesis = episode.get("synthesis")
    messages = episode.get("messages")

    if export._contains_secret(episode):
        flags.append("secret_detected")
    if not isinstance(messages, list) or not export._first_user_message(
        messages
    ):
        flags.append("missing_user_message")

    gate = export._command_gate(episode)
    repaired_ok = _has_positive_repair(episode)
    if gate == "reject" and not repaired_ok:
        flags.append("gate_reject_unrepaired")
    if isinstance(outcome, dict):
        rc = outcome.get("execute_rc")
        if isinstance(rc, int) and rc != 0:
            flags.append("execute_failed")
        if outcome.get("denied"):
            flags.append("user_denied")

    if isinstance(synthesis, dict):
        command = str(synthesis.get("command") or "").strip()
        if command and not command.startswith("chemsmart"):
            flags.append("malformed_command")
        if export._is_submission_command(command):
            terminal_reason = export._terminal_state_rejection_reason(episode)
            if terminal_reason:
                flags.append(terminal_reason)

    final_answer = str(episode.get("final_answer") or "").strip()
    assistant_answer = export._last_assistant_content(messages).strip()
    if not episode.get("paused") and not final_answer and not assistant_answer:
        # A completed turn that never answered the user teaches silence.
        if _invoked(episode):
            flags.append("empty_final_answer")
        else:
            flags.append("terminal_nosynth")

    return flags


def _has_positive_repair(episode: JsonDict) -> bool:
    export = _export_module()
    for event in export._tool_events(episode):
        if event.get("tool") != "repair_command":
            continue
        result = event.get("result_summary")
        semantic = result.get("semantic") if isinstance(result, dict) else None
        if isinstance(semantic, dict) and semantic.get("verdict") in {
            "ok",
            "warn",
        }:
            return True
    return False


def _invoked(episode: JsonDict) -> list[str]:
    return _export_module()._invoked_tools(episode)


# ------------------------------------------------------------- stability


def stability_flags(episode: JsonDict) -> list[str]:
    """Multi-turn / tool-loop stability signals for one episode."""

    flags: list[str] = []
    export = _export_module()
    events = export._tool_events(episode)

    saw_synthesize = False
    signatures: Counter[str] = Counter()
    for event in events:
        tool = str(event.get("tool") or "")
        if tool == "synthesize_command":
            saw_synthesize = True
        if tool == "execute_chemsmart_command" and not saw_synthesize:
            flags.append("execute_without_synthesize")
        signature = json.dumps(
            [tool, event.get("normalized_args") or event.get("args") or {}],
            sort_keys=True,
            default=str,
        )
        signatures[signature] += 1
    if any(count >= 3 for count in signatures.values()):
        flags.append("repeated_identical_tool_call")
    if (
        any(str(event.get("status") or "") == "error" for event in events)
        and not str(episode.get("final_answer") or "").strip()
    ):
        flags.append("tool_error_unresolved")
    return sorted(set(flags))


def _user_turns(episode: JsonDict) -> int:
    messages = episode.get("messages")
    if not isinstance(messages, list):
        return 0
    return sum(
        1
        for message in messages
        if isinstance(message, dict) and message.get("role") == "user"
    )


# ------------------------------------------------------------- sota metrics


def trajectory_fingerprint(episode: JsonDict) -> str:
    """Deterministic, order-preserving fingerprint of the tool-call sequence.

    Semantic trajectory dedup (Action-Space Hash). The previous
    ``hash(frozenset(bigrams))`` had two defects proven on the real store:
    (1) Python's ``hash`` is salted per process (PYTHONHASHSEED) so the
    fingerprint was non-reproducible across runs; (2) every <2-tool episode
    has zero bigrams and collapsed into ONE bucket, inflating the duplicate
    count (821 of 1161 reported dups were this artifact). A stable hash of the
    ORDERED sequence fixes both: empty and single-tool trajectories are now
    distinct by their actual content.
    """

    export = _export_module()
    tools = [
        str(event.get("tool") or "") for event in export._tool_events(episode)
    ]
    key = "|".join(tools) if tools else "∅"
    return hashlib.sha1(key.encode("utf-8")).hexdigest()


def schema_efficiency_score(episode: JsonDict) -> float | None:
    """Fraction of tool calls that succeeded without an error status.

    ``1.0`` = perfect adherence (no schema/runtime errors forcing retries).
    Returns ``None`` for episodes that called no tools (chitchat/advisory) —
    efficiency is undefined there; forcing 0.0 wrongly deflated the average
    (0.77 vs the true tool-bearing 0.90 on the real store).
    """

    export = _export_module()
    events = export._tool_events(episode)
    total_actions = len(events)
    if total_actions == 0:
        return None
    schema_errors = sum(
        1 for event in events if str(event.get("status") or "") == "error"
    )
    return max(0.0, (total_actions - schema_errors) / total_actions)


# ------------------------------------------------------------------ audit


def audit(
    episode_paths: Iterable[Path],
    *,
    max_per_skeleton: int = 20,
    max_skeleton_share: float = 0.05,
) -> JsonDict:
    export = _export_module()
    raw_count = 0
    episodes: list[JsonDict] = []

    def _counting(source):
        nonlocal raw_count
        for episode in source:
            raw_count += 1
            yield episode

    episodes = list(
        export._dedup_episodes(_counting(export._iter_episodes(episode_paths)))
    )
    chains = list(export._session_chains(episodes))
    recovered_cross_turn_rejects: set[tuple[str, Any]] = set()
    for chain in chains:
        repair_record, _ = export._cross_turn_repair_pair_record(chain)
        if repair_record is None:
            continue
        meta = repair_record.get("meta") or {}
        recovered_cross_turn_rejects.add(
            (str(meta.get("session_id") or ""), meta.get("rejected_turn"))
        )

    trajectories: Counter[str] = Counter()
    kinds: Counter[str] = Counter()
    providers: Counter[str] = Counter()
    skeletons: Counter[str] = Counter()
    skeleton_members: dict[str, list[int]] = defaultdict(list)
    quality: Counter[str] = Counter()
    stability: Counter[str] = Counter()
    flagged: list[JsonDict] = []
    clean_indices: list[int] = []
    record_level_multi_turn = 0
    paused_orphans = 0
    exact_pairs: Counter[tuple[str, str]] = Counter()

    # SOTA Metrics tracking
    trajectory_transitions: list[tuple[str, str]] = []
    semantic_fingerprints: Counter[str] = Counter()
    multistep_fingerprints: Counter[str] = Counter()
    efficiency_scores: list[float] = []
    terminal_state_counts: Counter[str] = Counter()
    terminal_state_failures: Counter[str] = Counter()

    for index, episode in enumerate(episodes):
        trajectory = " -> ".join(_invoked(episode)) or "(no tools)"
        trajectories[trajectory] += 1

        provider = episode.get("provider") or {}
        providers[
            f"{provider.get('name') or '?'}:{provider.get('model') or '?'}"
        ] += 1

        synthesis = episode.get("synthesis")
        if isinstance(synthesis, dict):
            command = str(synthesis.get("command") or "")
            kind = _kind_of_command(command)
            if kind:
                kinds[kind] += 1
            user = export._first_user_message(episode.get("messages"))
            if user and command:
                exact_pairs[(user, command)] += 1
            if export._is_submission_command(command):
                if episode.get("terminal_state") is None:
                    terminal_state_counts["missing_for_submit"] += 1

        state = episode.get("terminal_state")
        if state is not None:
            terminal_state_counts["present"] += 1
            state_issues = validate_terminal_state(state)
            if state_issues:
                terminal_state_counts["invalid"] += 1
                terminal_state_failures.update(state_issues)
            else:
                terminal_state_counts["passed"] += 1
                if terminal_state_is_positive(state):
                    terminal_state_counts["positive"] += 1
                else:
                    terminal_state_counts["expected_failure"] += 1

        user = export._first_user_message(episode.get("messages"))
        if user:
            skeleton = export._query_skeleton(user)
            skeletons[skeleton] += 1
            skeleton_members[skeleton].append(index)

        if _user_turns(episode) >= 2:
            record_level_multi_turn += 1
        if episode.get("paused"):
            paused_orphans += 1  # dedup already removed superseded pauses

        # SOTA Metrics gathering
        tool_sequence = [
            str(event.get("tool") or "")
            for event in export._tool_events(episode)
        ]
        trajectory_transitions.extend(
            list(zip(tool_sequence, tool_sequence[1:]))
        )
        fp = trajectory_fingerprint(episode)
        semantic_fingerprints[fp] += 1
        # Single/zero-tool "duplicates" are the expected command-synthesis
        # norm, not dataset flooding; only multi-step trajectories signal it.
        if len(tool_sequence) >= 2:
            multistep_fingerprints[fp] += 1
        eff = schema_efficiency_score(episode)
        if eff is not None:
            efficiency_scores.append(eff)

        q_flags = quality_flags(episode)
        episode_key = (
            str(episode.get("session_id") or ""),
            episode.get("turn"),
        )
        if episode_key in recovered_cross_turn_rejects:
            q_flags = [
                flag for flag in q_flags if flag != "gate_reject_unrepaired"
            ]
        s_flags = stability_flags(episode)
        for flag in q_flags:
            quality[flag] += 1
        for flag in s_flags:
            stability[flag] += 1
        if q_flags or s_flags:
            flagged.append(
                {
                    "session_id": episode.get("session_id"),
                    "turn": episode.get("turn"),
                    "quality": q_flags,
                    "stability": s_flags,
                }
            )
        if not q_flags:
            clean_indices.append(index)

    total = len(episodes)
    multi_turn_chains = [
        chain
        for chain in chains
        if export._user_message_count(export._merge_chain_messages(chain)) >= 2
    ]
    episodes_in_multi_turn_chains = sum(
        len(chain) for chain in multi_turn_chains
    )
    chain_skip_reasons: Counter[str] = Counter()
    trainable_chains = 0
    for chain in chains:
        reason = export._chain_terminal_rejection_reason(
            chain,
            include_rejected=False,
        )
        if reason:
            chain_skip_reasons[reason] += 1
        else:
            trainable_chains += 1
    distinct = len(skeletons)
    over_cap = {
        skeleton: count
        for skeleton, count in skeletons.items()
        if count > max_per_skeleton
        or (total and count / total > max_skeleton_share)
    }
    duplicate_pairs = sum(
        count - 1 for count in exact_pairs.values() if count > 1
    )

    # Skeleton cap for the clean set: keep at most max_per_skeleton of each.
    capped: set[int] = set()
    for skeleton, members in skeleton_members.items():
        for member in members[max_per_skeleton:]:
            capped.add(member)
    clean_after_cap = [i for i in clean_indices if i not in capped]

    # Calculate Tool-Use Diversity Score (Shannon Entropy)
    transition_counts = Counter(trajectory_transitions)
    total_transitions = sum(transition_counts.values())
    trajectory_entropy = 0.0
    if total_transitions > 0:
        trajectory_entropy = -sum(
            (c / total_transitions) * math.log2(c / total_transitions)
            for c in transition_counts.values()
        )

    avg_efficiency = (
        sum(efficiency_scores) / len(efficiency_scores)
        if efficiency_scores
        else 0.0
    )
    # Duplicate multi-step approaches only (single-step is the expected norm).
    semantic_duplicates = sum(
        count - 1 for count in multistep_fingerprints.values() if count > 1
    )
    distinct_trajectories = len(semantic_fingerprints)

    report = {
        "episodes_raw": raw_count,
        "episodes_after_dedup": total,
        "pipelines": {
            "by_trajectory": dict(trajectories.most_common()),
            "by_job_kind": dict(sorted(kinds.items())),
            "canonical_kind_coverage": {
                "coverage_scope": "CLI-emittable command kinds",
                "command_emittable": list(COMMAND_CANONICAL_AGENT_KINDS),
                "project_yaml_only": list(PROJECT_YAML_ONLY_KINDS),
                "covered": sorted(
                    kind
                    for kind in kinds
                    if kind in COMMAND_CANONICAL_AGENT_KINDS
                ),
                "missing": sorted(
                    kind
                    for kind in COMMAND_CANONICAL_AGENT_KINDS
                    if kind not in kinds
                ),
                "noncanonical": sorted(
                    kind
                    for kind in kinds
                    if kind not in COMMAND_CANONICAL_AGENT_KINDS
                ),
            },
            "by_provider": dict(sorted(providers.items())),
        },
        "quality": {
            "flag_counts": dict(sorted(quality.items())),
            "flagged_episodes": flagged[:200],
            "clean_episodes": len(clean_indices),
            "exact_duplicate_pairs": duplicate_pairs,
            "recovered_cross_turn_rejects": len(recovered_cross_turn_rejects),
            "trainable_session_chains": trainable_chains,
            "review_session_chains": len(chains) - trainable_chains,
            "session_chain_skip_reasons": dict(
                sorted(chain_skip_reasons.items())
            ),
        },
        "diversity": {
            "distinct_query_skeletons": distinct,
            "distinct_query_skeleton_ratio": (
                round(distinct / total, 4) if total else 0.0
            ),
            "top_query_skeletons": skeletons.most_common(10),
            "over_frequent_skeletons": over_cap,
            "capped_out_episodes": len(capped),
            "trajectory_entropy": round(trajectory_entropy, 4),
            "distinct_tool_trajectories": distinct_trajectories,
            "semantic_trajectory_duplicates": semantic_duplicates,
        },
        "stability": {
            "flag_counts": dict(sorted(stability.items())),
            "session_count": len(chains),
            "multi_turn_sessions": len(multi_turn_chains),
            "multi_turn_session_share": (
                round(len(multi_turn_chains) / len(chains), 4)
                if chains
                else 0.0
            ),
            "episodes_in_multi_turn_sessions": episodes_in_multi_turn_chains,
            "episode_share_in_multi_turn_sessions": (
                round(episodes_in_multi_turn_chains / total, 4)
                if total
                else 0.0
            ),
            # Back-compatible aliases now use the session-level definition.
            "multi_turn_episodes": len(multi_turn_chains),
            "multi_turn_share": (
                round(len(multi_turn_chains) / len(chains), 4)
                if chains
                else 0.0
            ),
            "record_level_multi_turn_episodes": record_level_multi_turn,
            "unresolved_ask_user_pauses": paused_orphans,
            "average_schema_efficiency_score": round(avg_efficiency, 4),
            "schema_efficiency_sampled_episodes": len(efficiency_scores),
        },
        "terminal_state": {
            "counts": dict(sorted(terminal_state_counts.items())),
            "failure_rule_counts": dict(
                sorted(terminal_state_failures.items())
            ),
        },
        "clean_after_cap": len(clean_after_cap),
    }
    report["_clean_indices"] = clean_after_cap
    report["_episodes"] = episodes
    return report


def _kind_of_command(command: str) -> str | None:
    try:
        from chemsmart.agent.model_command_parser import parse_model_command

        parsed = parse_model_command(command)
    except Exception:
        return None
    if parsed.parse_error or parsed.program not in {"gaussian", "orca"}:
        return None
    job = str(parsed.job or "")
    if not job:
        return None
    # QM/MM is nested below a parent job, e.g. ``opt qmmm``. The parser keeps
    # the parent as ``job``; coverage needs the effective workflow instead.
    if "qmmm" in parsed.tokens:
        job = "qmmm"
    elif parsed.program == "gaussian" and job == "td":
        job = "tddft"
    if parsed.program == "gaussian" and job == "modred":
        scan_definition = str(
            (parsed.structural_options or {}).get("scan_definition") or ""
        )
        if any(
            " S " in f" {line.strip()} "
            for line in scan_definition.splitlines()
        ):
            job = "scan"
    return f"{parsed.program}.{job}"


HEALTH_THRESHOLDS = {
    # From the v1→v4 boilerplate-trap lesson: below this the corpus is
    # template-dominated and hyperparameter tuning will not save it.
    "min_distinct_query_skeleton_ratio": 0.35,
    "max_flagged_share": 0.30,
}


def health_failures(report: JsonDict) -> list[str]:
    failures: list[str] = []
    total = report["episodes_after_dedup"]
    if not total:
        return ["empty_dataset"]
    ratio = report["diversity"]["distinct_query_skeleton_ratio"]
    if ratio < HEALTH_THRESHOLDS["min_distinct_query_skeleton_ratio"]:
        failures.append(
            f"skeleton_ratio {ratio} < "
            f"{HEALTH_THRESHOLDS['min_distinct_query_skeleton_ratio']}"
        )
    flagged = total - report["quality"]["clean_episodes"]
    if flagged / total > HEALTH_THRESHOLDS["max_flagged_share"]:
        failures.append(
            f"flagged_share {round(flagged / total, 4)} > "
            f"{HEALTH_THRESHOLDS['max_flagged_share']}"
        )
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("episodes", nargs="*")
    parser.add_argument("--training-dir", default=None)
    parser.add_argument(
        "--include-runs",
        action="store_true",
        help=(
            "Also audit model-specific run stores under "
            "<training-dir>/runs/*/episodes/*.jsonl. Use this for "
            "multi-teacher API accumulation."
        ),
    )
    parser.add_argument("--out", default=None)
    parser.add_argument("--max-per-skeleton", type=int, default=20)
    parser.add_argument("--max-skeleton-share", type=float, default=0.05)
    parser.add_argument("--write-clean", default=None)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    export = _export_module()
    training_dir = export._training_dir(args.training_dir)
    paths = export._episode_paths(
        args.episodes,
        training_dir,
        include_runs=args.include_runs,
    )

    report = audit(
        paths,
        max_per_skeleton=args.max_per_skeleton,
        max_skeleton_share=args.max_skeleton_share,
    )
    episodes = report.pop("_episodes")
    clean_indices = report.pop("_clean_indices")

    if args.write_clean:
        clean_path = Path(args.write_clean).expanduser()
        clean_path.parent.mkdir(parents=True, exist_ok=True)
        with clean_path.open("w", encoding="utf-8") as handle:
            for index in clean_indices:
                handle.write(
                    json.dumps(episodes[index], sort_keys=True, default=str)
                    + "\n"
                )
        report["clean_output"] = str(clean_path)

    failures = health_failures(report)
    report["health_failures"] = failures

    payload = json.dumps(report, indent=2, sort_keys=True, default=str)
    if args.out:
        Path(args.out).expanduser().write_text(
            payload + "\n", encoding="utf-8"
        )
    print(payload)
    if args.strict and failures:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
