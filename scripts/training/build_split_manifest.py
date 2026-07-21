"""Build a deterministic, leakage-aware train/eval split for SFT exports.

The split unit is a connected component of source session IDs, scenario
families, and *all* normalized user-turn skeletons.  This is stricter than
random record splitting and stricter than the earlier first-turn-only
construction: a multi-turn session, any repeated query shape (first turn or
later), a shared scenario family, and its family derivatives cannot straddle
train and evaluation.

Component membership decides the split; nothing is force-routed to eval.
Preference families (repair/wrong-route) and evidence families
(terminal-state) follow their source component's split into clearly named
``preference``/``evidence`` sub-directories, so a repaired positive chain or a
terminal-success chain can legitimately land in train while the contrast/
receipt records for the same component stay out of the positive SFT stream.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Callable, Iterable

JsonDict = dict[str, Any]

# Positive SFT families are written at the split root; preference/evidence
# families are written into named sub-directories under the same split so they
# never mix into the positive stream but still honour the component split.
POSITIVE_SFT_FAMILIES = (
    "tool_loop_sft",
    "compact_spec_sft",
    "command_answer_sft",
    "project_yaml_sft",
    "reasoning_synthesis_sft",
)
PREFERENCE_FAMILIES = ("repair_pairs", "wrong_route_contrast")
EVIDENCE_FAMILIES = ("terminal_state_assertions",)
FAMILIES = POSITIVE_SFT_FAMILIES + PREFERENCE_FAMILIES + EVIDENCE_FAMILIES

_KIND_TOKENS = frozenset(
    {
        "crest",
        "dias",
        "irc",
        "modred",
        "neb",
        "nci",
        "opt",
        "qmmm",
        "qrc",
        "resp",
        "scan",
        "sp",
        "tddft",
        "td",
        "ts",
        "wbi",
        "freq",
    }
)


def _family_subdir(family: str) -> str:
    if family in PREFERENCE_FAMILIES:
        return "preference"
    if family in EVIDENCE_FAMILIES:
        return "evidence"
    return ""


class _UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, key: str) -> None:
        self.parent.setdefault(key, key)

    def find(self, key: str) -> str:
        self.add(key)
        parent = self.parent[key]
        if parent != key:
            self.parent[key] = self.find(parent)
        return self.parent[key]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def _skeleton(text: str) -> str:
    value = str(text or "").lower()
    value = re.sub(
        r"[A-Za-z0-9_./+-]+\.(?:xyz|log|out|inp|com|gjf|yaml|yml)",
        "<file>",
        value,
    )
    value = re.sub(r"\b\d+(?:\.\d+)?\b", "#", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _record_meta(record: JsonDict) -> JsonDict:
    value = record.get("meta")
    return value if isinstance(value, dict) else {}


def _all_user_texts(record: JsonDict) -> list[str]:
    """Collect every user-turn text a record carries, not just the first.

    Tool-loop rows merge a multi-turn conversation, cross-turn repair prompts
    pack ``Turn N: ...`` lines, and review rows keep a bare ``user`` field.
    Skeletons of *all* of these become linking keys so a repeated later-turn
    query shape cannot straddle the split.
    """

    texts: list[str] = []
    prompt = record.get("prompt")
    if isinstance(prompt, str) and prompt.strip():
        for line in prompt.splitlines():
            line = re.sub(r"^\s*turn\s+\d+:\s*", "", line.strip(), flags=re.I)
            if line:
                texts.append(line)
    user = record.get("user")
    if isinstance(user, str) and user.strip():
        texts.append(user.strip())
    for message in record.get("messages") or []:
        if not isinstance(message, dict) or message.get("role") != "user":
            continue
        content = message.get("content")
        if isinstance(content, str) and content.strip():
            texts.append(content.strip())
    return texts


def _sessions(record: JsonDict) -> set[str]:
    meta = _record_meta(record)
    values = {
        str(meta.get("session_id") or "").strip(),
        str(meta.get("wrong_session_id") or "").strip(),
    }
    return {value for value in values if value}


def _scenario_families(record: JsonDict) -> set[str]:
    provenance = _record_meta(record).get("dataset_provenance")
    if not isinstance(provenance, dict):
        return set()
    family = str(provenance.get("scenario_family") or "").strip()
    return {family} if family else set()


def _record_key(family: str, index: int) -> str:
    return f"record:{family}:{index}"


def _record_info(family: str, index: int, record: JsonDict) -> JsonDict:
    skeletons = sorted(
        {
            skeleton
            for text in _all_user_texts(record)
            if (skeleton := _skeleton(text))
        }
    )
    sessions = _sessions(record)
    families = _scenario_families(record)
    keys: set[str] = set()
    keys.update(f"skel:{skeleton}" for skeleton in skeletons)
    keys.update(f"session:{session}" for session in sessions)
    keys.update(f"family:{family_id}" for family_id in families)
    if not keys:
        keys.add(_record_key(family, index))
    return {
        "family": family,
        "index": index,
        "record": record,
        "session_ids": sorted(sessions),
        "scenario_families": sorted(families),
        "skeletons": skeletons,
        "keys": sorted(keys),
        "has_provenance": bool(_record_meta(record).get("dataset_provenance")),
    }


def _load_records(export_dir: Path) -> list[JsonDict]:
    records: list[JsonDict] = []
    for family in FAMILIES:
        path = export_dir / f"{family}.jsonl"
        if not path.exists():
            continue
        for index, line in enumerate(
            path.read_text(encoding="utf-8").splitlines()
        ):
            if line.strip():
                records.append(_record_info(family, index, json.loads(line)))
    return records


def _assign(records: list[JsonDict], eval_fraction: int) -> None:
    """Union records into components; split each component deterministically.

    No record is force-routed to eval.  The split of a component is a stable
    hash of its canonical (lexicographically smallest) member key, so it is
    independent of union order and reproducible across runs.
    """

    union = _UnionFind()
    for info in records:
        keys = info["keys"]
        for key in keys:
            union.add(key)
        for key in keys[1:]:
            union.union(keys[0], key)

    grouped: dict[str, list[JsonDict]] = defaultdict(list)
    for info in records:
        grouped[union.find(info["keys"][0])].append(info)

    for members in grouped.values():
        all_keys = sorted({key for info in members for key in info["keys"]})
        canonical = all_keys[0]
        digest = hashlib.sha256(canonical.encode("utf-8")).digest()[0] % 100
        split = "eval" if digest < eval_fraction else "train"
        for info in members:
            info["component"] = canonical
            info["split"] = split


def _write_split(export_dir: Path, records: list[JsonDict]) -> dict[str, Any]:
    split_root = export_dir / "splits"
    counts: dict[str, Counter[str]] = {"train": Counter(), "eval": Counter()}
    for info in records:
        split = info["split"]
        family = info["family"]
        subdir = _family_subdir(family)
        out_dir = split_root / split / subdir if subdir else split_root / split
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / f"{family}.jsonl"
        with path.open("a", encoding="utf-8") as handle:
            handle.write(
                json.dumps(info["record"], ensure_ascii=False, sort_keys=True)
                + "\n"
            )
        counts[split][family] += 1
    return {
        split: dict(sorted(counter.items()))
        for split, counter in counts.items()
    }


def _cross_split(
    records: list[JsonDict], selector: Callable[[JsonDict], Iterable[str]]
) -> list[str]:
    mapping: dict[str, set[str]] = defaultdict(set)
    for info in records:
        for value in selector(info):
            mapping[value].add(info["split"])
    return sorted(key for key, splits in mapping.items() if len(splits) > 1)


def _leakage(records: list[JsonDict]) -> dict[str, Any]:
    session_leaks = _cross_split(records, lambda info: info["session_ids"])
    skeleton_leaks = _cross_split(records, lambda info: info["skeletons"])
    family_leaks = _cross_split(
        records, lambda info: info["scenario_families"]
    )
    return {
        "session_leaks": session_leaks,
        "all_turn_skeleton_leaks": skeleton_leaks,
        "scenario_family_leaks": family_leaks,
        "session_leak_count": len(session_leaks),
        "all_turn_skeleton_leak_count": len(skeleton_leaks),
        "scenario_family_leak_count": len(family_leaks),
    }


def _command_template(command: str) -> str:
    value = command.lower()
    value = re.sub(
        r"[a-z0-9_./+-]+\.(?:xyz|log|out|inp|com|gjf|yaml|yml)",
        "<file>",
        value,
    )
    value = re.sub(r"\b\d+(?:\.\d+)?\b", "#", value)
    value = re.sub(r"\s+", " ", value).strip()
    return value


def _record_commands(record: JsonDict) -> list[str]:
    commands: list[str] = []
    meta = _record_meta(record)
    single = meta.get("command")
    if isinstance(single, str) and single.strip():
        commands.append(single.strip())
    for item in meta.get("commands") or []:
        if isinstance(item, str) and item.strip():
            commands.append(item.strip())
    chosen = record.get("chosen")
    if isinstance(chosen, dict):
        command = chosen.get("command")
        if isinstance(command, str) and command.strip():
            commands.append(command.strip())
    terminal = record.get("terminal_state")
    if isinstance(terminal, dict):
        command = terminal.get("command")
        if isinstance(command, str) and command.strip():
            commands.append(command.strip())
    return commands


def _command_template_overlap(records: list[JsonDict]) -> dict[str, Any]:
    """Diagnostic only: canonical commands legitimately recur across phrasings.

    A shared normalized command template across splits is reported, not gated;
    the same gated command can arise from genuinely different query shapes.
    """

    templates: dict[str, set[str]] = defaultdict(set)
    for info in records:
        for command in _record_commands(info["record"]):
            template = _command_template(command)
            if template:
                templates[template].add(info["split"])
    overlap = sorted(
        template for template, splits in templates.items() if len(splits) > 1
    )
    return {
        "distinct_command_templates": len(templates),
        "overlapping_command_templates": len(overlap),
        "examples": overlap[:20],
    }


def _per_family_eval_shares(counts: dict[str, Any]) -> dict[str, Any]:
    shares: dict[str, Any] = {}
    for family in FAMILIES:
        train = counts["train"].get(family, 0)
        evaluation = counts["eval"].get(family, 0)
        total = train + evaluation
        shares[family] = {
            "train": train,
            "eval": evaluation,
            "total": total,
            "eval_share": round(evaluation / total, 4) if total else 0.0,
        }
    return shares


def _provenance_coverage(records: list[JsonDict]) -> dict[str, Any]:
    total = len(records)
    with_provenance = sum(1 for info in records if info["has_provenance"])
    with_family = sum(1 for info in records if info["scenario_families"])
    by_family: dict[str, Any] = {}
    for family in FAMILIES:
        members = [info for info in records if info["family"] == family]
        covered = sum(1 for info in members if info["has_provenance"])
        by_family[family] = {
            "records": len(members),
            "with_dataset_provenance": covered,
        }
    return {
        "records": total,
        "with_dataset_provenance": with_provenance,
        "with_scenario_family": with_family,
        "dataset_provenance_ratio": (
            round(with_provenance / total, 4) if total else 0.0
        ),
        "scenario_family_ratio": (
            round(with_family / total, 4) if total else 0.0
        ),
        "by_family": by_family,
    }


def _recovered_reject_counts(records: list[JsonDict]) -> dict[str, int]:
    counts = {"train": 0, "eval": 0}
    for info in records:
        if info["family"] != "tool_loop_sft":
            continue
        meta = _record_meta(info["record"])
        if meta.get("recovered_reject"):
            counts[info["split"]] += 1
    return counts


def _terminal_state_counts(records: list[JsonDict]) -> dict[str, Any]:
    assertions = {"train": 0, "eval": 0}
    chains = {"train": 0, "eval": 0}
    for info in records:
        split = info["split"]
        if info["family"] == "terminal_state_assertions":
            assertions[split] += 1
        elif info["family"] == "tool_loop_sft":
            if _record_meta(info["record"]).get("terminal_state"):
                chains[split] += 1
    return {
        "terminal_state_assertions": assertions,
        "tool_loop_terminal_success_chains": chains,
    }


def build_split_manifest(
    export_dir: Path, eval_fraction: int
) -> dict[str, Any]:
    """Assign splits, (re)write split JSONL files, and return the manifest."""

    export_dir = Path(export_dir).expanduser().resolve()
    records = _load_records(export_dir)
    _assign(records, eval_fraction)

    # The output directory is regenerated from the same source export, so no
    # stale rows can survive a rerun.
    split_root = export_dir / "splits"
    if split_root.exists():
        for path in sorted(split_root.rglob("*.jsonl")):
            path.unlink()
    counts = _write_split(export_dir, records)

    leakage = _leakage(records)
    components = {info["component"] for info in records}
    manifest = {
        "schema_version": 2,
        "source_export": str(export_dir),
        "eval_fraction_percent": eval_fraction,
        "record_count": len(records),
        "component_count": len(components),
        "counts": counts,
        "per_family_eval_shares": _per_family_eval_shares(counts),
        "provenance_coverage": _provenance_coverage(records),
        "leakage": leakage,
        "command_template_overlap": _command_template_overlap(records),
        "recovered_rejects_by_split": _recovered_reject_counts(records),
        "terminal_state_by_split": _terminal_state_counts(records),
        "production_gate": {
            "no_session_cross_split": leakage["session_leak_count"] == 0,
            "no_all_turn_skeleton_cross_split": (
                leakage["all_turn_skeleton_leak_count"] == 0
            ),
            "no_scenario_family_cross_split": (
                leakage["scenario_family_leak_count"] == 0
            ),
            "deterministic": True,
        },
    }
    (export_dir / "split_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export-dir", required=True)
    parser.add_argument("--eval-fraction", type=int, default=20)
    args = parser.parse_args()
    if not 1 <= args.eval_fraction <= 99:
        raise SystemExit("--eval-fraction must be between 1 and 99")
    manifest = build_split_manifest(Path(args.export_dir), args.eval_fraction)
    print(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True))
    return 0 if all(manifest["production_gate"].values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
