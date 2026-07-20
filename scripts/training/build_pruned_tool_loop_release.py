#!/usr/bin/env python3
"""Build an immutable, tokenizer-bounded tool-loop release from frozen v1.

The source release and append-only ledger are never modified.  This builder
replaces only tool-loop rows with deterministic projections and decision
windows, then records every exclusion and content hash in the derived release.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from collections import Counter
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Protocol

from chemsmart.agent.loop import ASK_USER_TOOL_DEF
from chemsmart.agent.prompts import load_prompt
from chemsmart.agent.prompts.identity import build_system_prompt
from chemsmart.agent.registry import ToolRegistry
from scripts.training.export_sft import _contains_secret

JsonDict = dict[str, Any]
MAX_SYSTEM_CHARS = 4096
MAX_SYNTHESIS_REASONING_CHARS = 8192
DEFAULT_MAX_TOKENS = 8192
_HOME_RE = re.compile(r"/(?:Users|home)/[^/\s]+")
_HIDDEN_KEYS = {
    "chain_of_thought",
    "cot",
    "hidden_reasoning",
    "reasoning_content",
}


class ChatTokenizer(Protocol):
    def apply_chat_template(self, conversation: Any, **kwargs: Any) -> Any: ...


def project_tool_loop_record(
    record: JsonDict,
    *,
    registry: ToolRegistry,
    tokenizer: ChatTokenizer,
    max_tokens: int = DEFAULT_MAX_TOKENS,
) -> tuple[list[JsonDict], str]:
    """Project one v1 chain into one or more bounded v2 windows."""

    messages = record.get("messages")
    if not isinstance(messages, list) or not messages:
        return [], "missing_messages"
    system = next(
        (
            str(message.get("content") or "")
            for message in messages
            if isinstance(message, dict) and message.get("role") == "system"
        ),
        "",
    )
    if len(system) >= 100_000:
        return [], "full_schema_system_prompt"
    if _has_unprovenanced_long_reasoning(messages):
        return [], "unprovenanced_long_reasoning"

    projected = _project_messages(messages)
    if not projected:
        return [], "empty_after_projection"
    pair_error = _tool_pair_error(projected)
    if pair_error:
        return [], pair_error

    first_user = next(
        (
            str(message.get("content") or "")
            for message in projected
            if message.get("role") == "user"
        ),
        "",
    )
    workspace = record.get("meta", {}).get("workspace") or {}
    system_prompt = _canonical_system_prompt(
        registry, request=first_user, workspace=None
    )
    projected = [
        {"role": "system", "content": system_prompt},
        *[message for message in projected if message.get("role") != "system"],
    ]
    windows = _decision_windows(
        projected,
        registry=registry,
        tokenizer=tokenizer,
        max_tokens=max_tokens,
        workspace=workspace,
    )
    if not windows:
        return [], "oversize_atomic_window"

    meta = dict(record.get("meta") or {})
    output: list[JsonDict] = []
    for index, (window_messages, tools, token_count) in enumerate(windows, 1):
        if (
            len(str(window_messages[0].get("content") or ""))
            > MAX_SYSTEM_CHARS
        ):
            return [], "system_prompt_over_budget"
        if _tool_pair_error(window_messages):
            return [], "window_tool_pair_error"
        row_meta = {
            **meta,
            "family": "tool_loop_pruned_v2",
            "source_family": meta.get("family") or "tool_loop",
            "window_index": index,
            "window_count": len(windows),
            "qwen3_14b_tokens": token_count,
            "projection_schema": 2,
        }
        output.append(
            {"messages": window_messages, "tools": tools, "meta": row_meta}
        )
    return output, "ok"


def _project_messages(messages: list[Any]) -> list[JsonDict]:
    call_names: dict[str, str] = {}
    projected: list[JsonDict] = []
    for raw in messages:
        if not isinstance(raw, dict):
            continue
        role = raw.get("role")
        if role not in {"system", "user", "assistant", "tool"}:
            continue
        message: JsonDict = {"role": role}
        content = raw.get("content")
        if role == "tool":
            call_id = str(raw.get("tool_call_id") or "")
            name = str(raw.get("name") or call_names.get(call_id) or "")
            message["tool_call_id"] = call_id
            if name:
                message["name"] = name
            message["content"] = _project_tool_content(name, content)
        else:
            if isinstance(content, str):
                message["content"] = _mask_text(content)
            calls = raw.get("tool_calls")
            if role == "assistant" and isinstance(calls, list) and calls:
                clean_calls: list[JsonDict] = []
                for call in calls:
                    clean = _clean_tool_call(call)
                    if clean is None:
                        continue
                    clean_calls.append(clean)
                    call_names[clean["id"]] = clean["function"]["name"]
                if clean_calls:
                    message["tool_calls"] = clean_calls
        if (
            message.get("content")
            or message.get("tool_calls")
            or role == "tool"
        ):
            projected.append(message)
    return projected


def _clean_tool_call(value: Any) -> JsonDict | None:
    if not isinstance(value, dict):
        return None
    function = value.get("function")
    call_id = str(value.get("id") or "")
    if not call_id or not isinstance(function, dict):
        return None
    name = str(function.get("name") or "")
    if not name:
        return None
    arguments = function.get("arguments")
    if isinstance(arguments, str):
        arguments = _mask_text(arguments)
    elif isinstance(arguments, dict):
        arguments = json.dumps(
            _bounded_value(arguments), ensure_ascii=False, sort_keys=True
        )
    else:
        arguments = "{}"
    return {
        "id": call_id,
        "type": "function",
        "function": {"name": name, "arguments": arguments},
    }


def _project_tool_content(name: str, content: Any) -> str:
    if not isinstance(content, str):
        return json.dumps(
            _bounded_value(content), ensure_ascii=False, sort_keys=True
        )
    try:
        payload = json.loads(content)
    except (TypeError, ValueError, json.JSONDecodeError):
        return _mask_text(content)[:2000]
    if not isinstance(payload, dict):
        return json.dumps(
            _bounded_value(payload), ensure_ascii=False, sort_keys=True
        )

    if name == "synthesize_command":
        keys = (
            "ok",
            "status",
            "action",
            "command",
            "explanation",
            "confidence",
            "missing_info",
            "project",
            "schema_variant",
            "semantic",
            "intent",
        )
    elif name == "repair_command":
        keys = (
            "ok",
            "status",
            "repaired",
            "original_command",
            "command",
            "issues",
            "semantic",
        )
    elif "project_yaml" in name or name in {
        "extract_project_protocol",
        "render_project_yaml",
        "search_basis_sets",
    }:
        keys = (
            "ok",
            "verdict",
            "program",
            "project_name",
            "method",
            "td",
            "protocol_notes",
            "unsupported_yaml_features",
            "yaml_text",
            "issues",
            "runtime_summary",
            "written_path",
        )
    else:
        keys = (
            "ok",
            "status",
            "verdict",
            "command",
            "returncode",
            "issues",
            "failed_rule_ids",
            "generated_inputs",
            "terminal_state",
            "result",
            "summary",
            "question",
            "options",
            "error",
        )
    projected = {
        key: _project_result_field(key, payload[key])
        for key in keys
        if key in payload
    }
    return json.dumps(projected, ensure_ascii=False, sort_keys=True)


def _project_result_field(key: str, value: Any) -> Any:
    if key == "semantic" and isinstance(value, dict):
        generated = value.get("generated_inputs") or []
        return {
            field: _bounded_value(value[field])
            for field in (
                "verdict",
                "failed_rule_ids",
                "missing_info",
                "notice",
            )
            if field in value
        } | {
            "generated_inputs": [
                _compact_generated_input(item)
                for item in generated[:4]
                if isinstance(item, dict)
            ]
        }
    if key == "intent" and isinstance(value, dict):
        assertions = value.get("assertions") or []
        return {
            field: _bounded_value(value[field])
            for field in ("verdict", "failed_rule_ids")
            if field in value
        } | {
            "assertions": [
                {
                    field: _bounded_value(assertion[field])
                    for field in ("id", "status", "expected", "observed")
                    if field in assertion
                }
                for assertion in assertions[:16]
                if isinstance(assertion, dict)
            ]
        }
    if key == "issues" and isinstance(value, list):
        return [
            {
                field: _bounded_value(issue[field])
                for field in ("rule_id", "severity", "message", "missing_info")
                if field in issue
            }
            for issue in value[:12]
            if isinstance(issue, dict)
        ]
    if key == "generated_inputs" and isinstance(value, list):
        return [
            _compact_generated_input(item)
            for item in value[:4]
            if isinstance(item, dict)
        ]
    if key == "runtime_summary" and isinstance(value, dict):
        return {
            str(job): {
                field: _bounded_value(settings[field])
                for field in (
                    "functional",
                    "basis",
                    "solvent_model",
                    "solvent_id",
                    "freq",
                    "heavy_elements",
                    "heavy_elements_basis",
                    "light_elements_basis",
                )
                if isinstance(settings, dict) and field in settings
            }
            for job, settings in list(value.items())[:8]
        }
    return _bounded_value(value)


def _compact_generated_input(value: JsonDict) -> JsonDict:
    keep = (
        "path",
        "input_path",
        "program",
        "kind",
        "route",
        "route_line",
        "route_lines",
        "verdict",
        "failed_rule_ids",
    )
    projected = {
        key: _bounded_value(value[key]) for key in keep if key in value
    }
    content = value.get("content")
    if isinstance(content, str):
        # The full generated input is an artifact, not prompt context. Keep
        # enough chemistry syntax to audit route/directive preservation.
        projected["content_excerpt"] = _mask_text(content)[:1000]
    return projected


def _bounded_value(value: Any, *, depth: int = 0) -> Any:
    if depth >= 5:
        return "<nested-value>"
    if isinstance(value, str):
        return _mask_text(value)[:2000]
    if isinstance(value, list):
        return [_bounded_value(item, depth=depth + 1) for item in value[:16]]
    if isinstance(value, dict):
        drop = {
            "raw_response",
            "reasoning",
            "workflow_state",
            "stdout_tail",
            "stderr_tail",
            "checked_argv",
            "decision_trace",
        }
        return {
            str(key): _bounded_value(item, depth=depth + 1)
            for key, item in value.items()
            if str(key) not in drop
        }
    return value


def _has_unprovenanced_long_reasoning(messages: list[Any]) -> bool:
    for message in messages:
        if not isinstance(message, dict) or message.get("role") != "tool":
            continue
        content = message.get("content")
        if not isinstance(content, str):
            continue
        try:
            payload = json.loads(content)
        except (TypeError, ValueError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict):
            continue
        reasoning = payload.get("reasoning")
        if (
            isinstance(reasoning, str)
            and len(reasoning) > MAX_SYNTHESIS_REASONING_CHARS
            and payload.get("reasoning_provenance") != "public_decision_trace"
        ):
            return True
    return False


def _decision_windows(
    messages: list[JsonDict],
    *,
    registry: ToolRegistry,
    tokenizer: ChatTokenizer,
    max_tokens: int,
    workspace: Any,
) -> list[tuple[list[JsonDict], list[JsonDict], int]]:
    tools = _tool_defs(messages, registry)
    total = _token_count(tokenizer, messages, tools)
    if total <= max_tokens:
        return [(messages, tools, total)]

    system = messages[0]
    turns = _split_user_turns(messages[1:])
    windows: list[tuple[list[JsonDict], list[JsonDict], int]] = []
    for turn_index, turn in enumerate(turns):
        if not turn or turn[0].get("role") != "user":
            return []
        scoped_system = dict(system)
        if turn_index:
            scoped_system["content"] = _canonical_system_prompt(
                registry,
                request=str(turn[0].get("content") or ""),
                workspace=workspace,
            )
        blocks = _atomic_blocks(turn[1:])
        current = [scoped_system, turn[0]]
        for block in blocks:
            candidate = [*current, *block]
            candidate_tools = _tool_defs(candidate, registry)
            count = _token_count(tokenizer, candidate, candidate_tools)
            if count <= max_tokens:
                current = candidate
                continue
            if len(current) > 2:
                current_tools = _tool_defs(current, registry)
                current_count = _token_count(tokenizer, current, current_tools)
                windows.append((current, current_tools, current_count))
            current = [scoped_system, turn[0], *block]
            current_tools = _tool_defs(current, registry)
            if _token_count(tokenizer, current, current_tools) > max_tokens:
                return []
        if len(current) > 2:
            current_tools = _tool_defs(current, registry)
            windows.append(
                (
                    current,
                    current_tools,
                    _token_count(tokenizer, current, current_tools),
                )
            )
    return windows


def _split_user_turns(messages: list[JsonDict]) -> list[list[JsonDict]]:
    turns: list[list[JsonDict]] = []
    for message in messages:
        if message.get("role") == "user":
            turns.append([message])
        elif turns:
            turns[-1].append(message)
        else:
            return []
    return turns


def _atomic_blocks(messages: list[JsonDict]) -> list[list[JsonDict]]:
    blocks: list[list[JsonDict]] = []
    index = 0
    while index < len(messages):
        message = messages[index]
        block = [message]
        calls = (
            message.get("tool_calls")
            if message.get("role") == "assistant"
            else None
        )
        if isinstance(calls, list) and calls:
            expected = {str(call.get("id") or "") for call in calls}
            index += 1
            while (
                index < len(messages) and messages[index].get("role") == "tool"
            ):
                block.append(messages[index])
                expected.discard(
                    str(messages[index].get("tool_call_id") or "")
                )
                index += 1
            if expected:
                return []
            blocks.append(block)
            continue
        blocks.append(block)
        index += 1
    return blocks


def _tool_defs(
    messages: list[JsonDict], registry: ToolRegistry
) -> list[JsonDict]:
    names: list[str] = []
    for message in messages:
        for call in message.get("tool_calls") or []:
            name = str((call.get("function") or {}).get("name") or "")
            if name and name not in names:
                names.append(name)
    specs = [registry.get_tool(name) for name in names if name != "ask_user"]
    defs = registry.openai_tool_defs(
        [spec for spec in specs if spec is not None]
    )
    if "ask_user" in names:
        defs.append(ASK_USER_TOOL_DEF)
    return defs


def _tool_pair_error(messages: list[JsonDict]) -> str:
    calls = {
        str(call.get("id") or "")
        for message in messages
        for call in message.get("tool_calls") or []
    }
    results = {
        str(message.get("tool_call_id") or "")
        for message in messages
        if message.get("role") == "tool"
    }
    if "" in calls or "" in results:
        return "empty_tool_call_id"
    if results - calls:
        return "orphan_tool_result"
    if calls - results:
        return "orphan_tool_call"
    return ""


def _canonical_system_prompt(
    registry: ToolRegistry,
    *,
    request: str,
    workspace: Any,
) -> str:
    context = None
    if isinstance(workspace, dict) and workspace:
        context = {"entities": _bounded_value(workspace)}
    return build_system_prompt(
        registry=registry,
        stage_instructions=load_prompt("unified_agent.md"),
        session_meta={"stage": "tool_loop", "approval_mode": "driving"},
        conversation_context=context,
        request=request,
        max_chars=MAX_SYSTEM_CHARS,
    )


def _token_count(
    tokenizer: ChatTokenizer,
    messages: list[JsonDict],
    tools: list[JsonDict],
) -> int:
    # Count the exact JSONL representation used by training. Provider schemas
    # can contain tuples; JSON serialization converts them to lists and may
    # change the target chat template by a token.
    canonical_messages = json.loads(
        json.dumps(messages, ensure_ascii=False, sort_keys=True)
    )
    canonical_tools = json.loads(
        json.dumps(tools, ensure_ascii=False, sort_keys=True)
    )
    tokens = tokenizer.apply_chat_template(
        canonical_messages,
        tools=canonical_tools or None,
        tokenize=True,
        add_generation_prompt=False,
    )
    if isinstance(tokens, Mapping) and "input_ids" in tokens:
        tokens = tokens["input_ids"]
    if hasattr(tokens, "shape"):
        return int(tokens.shape[-1])
    if isinstance(tokens, list) and tokens and isinstance(tokens[0], list):
        return len(tokens[0])
    return len(tokens)


def _mask_text(value: str) -> str:
    return _HOME_RE.sub("<workspace-home>", value)


def _load_tokenizer(model: str) -> ChatTokenizer:
    from transformers import AutoTokenizer

    return AutoTokenizer.from_pretrained(model, trust_remote_code=True)


def _read_jsonl(path: Path) -> list[JsonDict]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _write_jsonl(path: Path, rows: list[JsonDict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(
            json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n"
            for row in rows
        ),
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_source(source: Path) -> None:
    if (source / "STATUS").read_text(encoding="utf-8").strip() != "FROZEN":
        raise ValueError("source release is not FROZEN")
    lock = json.loads(
        (source / "release.lock.json").read_text(encoding="utf-8")
    )
    failures = [
        relative
        for relative, receipt in lock.get("files", {}).items()
        if not (source / relative).is_file()
        or _sha256(source / relative) != receipt.get("sha256")
    ]
    if failures:
        raise ValueError(f"source release lock mismatch: {failures[:3]}")


def build_release(
    *,
    source: Path,
    output: Path,
    tokenizer: ChatTokenizer,
    tokenizer_name: str,
    max_tokens: int,
) -> JsonDict:
    _verify_source(source)
    if output.exists():
        raise FileExistsError(f"derived release already exists: {output}")
    shutil.copytree(source, output)
    for stale in (output / "STATUS", output / "release.lock.json"):
        stale.unlink(missing_ok=True)
    (output / "STATUS").write_text("BUILDING\n", encoding="utf-8")

    registry = ToolRegistry.default()
    report: JsonDict = {
        "schema_version": 2,
        "source_release": source.name,
        "release_id": output.name,
        "tokenizer": tokenizer_name,
        "max_tokens": max_tokens,
        "splits": {},
    }
    combined: list[JsonDict] = []
    review_rows: list[JsonDict] = []
    for split in ("train", "eval"):
        path = source / "dataset/splits" / split / "tool_loop_sft.jsonl"
        rows: list[JsonDict] = []
        reasons: Counter[str] = Counter()
        source_rows = _read_jsonl(path)
        for source_index, record in enumerate(source_rows):
            windows, reason = project_tool_loop_record(
                record,
                registry=registry,
                tokenizer=tokenizer,
                max_tokens=max_tokens,
            )
            reasons[reason] += 1
            rows.extend(windows)
            if reason != "ok":
                review_rows.append(
                    _review_receipt(
                        record,
                        split=split,
                        source_index=source_index,
                        reason=reason,
                    )
                )
        _write_jsonl(
            output / "dataset/splits" / split / "tool_loop_sft.jsonl",
            rows,
        )
        combined.extend(rows)
        report["splits"][split] = {
            "source_rows": len(source_rows),
            "output_windows": len(rows),
            "source_outcomes": dict(sorted(reasons.items())),
            "max_tokens": max(
                (row["meta"]["qwen3_14b_tokens"] for row in rows), default=0
            ),
            "max_system_chars": max(
                (len(row["messages"][0]["content"]) for row in rows), default=0
            ),
        }
    _write_jsonl(output / "dataset/tool_loop_sft.jsonl", combined)
    _write_jsonl(
        output / "dataset/tool_loop_pruning_review.jsonl", review_rows
    )
    report["review_receipts"] = len(review_rows)
    validation = output / "validation/tool_loop_pruning.json"
    validation.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    _refresh_manifests(output, report)
    _write_release_readme(output, report)
    postcheck = _validate_derived_dataset(
        output, tokenizer=tokenizer, max_tokens=max_tokens
    )
    (output / "validation/tool_loop_postcheck.json").write_text(
        json.dumps(postcheck, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    _freeze(output, report)
    return report


def _review_receipt(
    record: JsonDict,
    *,
    split: str,
    source_index: int,
    reason: str,
) -> JsonDict:
    meta = record.get("meta") or {}
    canonical = json.dumps(record, ensure_ascii=False, sort_keys=True)
    family = {
        "full_schema_system_prompt": "full_schema_fallback",
        "unprovenanced_long_reasoning": "long_reasoning_quarantine",
        "oversize_atomic_window": "oversize_atomic_window",
    }.get(reason, "projection_reject")
    return {
        "review_family": family,
        "reason": reason,
        "source_split": split,
        "source_row_index": source_index,
        "source_record_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "session_id": meta.get("session_id"),
        "schema_variant": meta.get("schema_variant"),
        "tool_trajectory": meta.get("tool_trajectory"),
    }


def _validate_derived_dataset(
    output: Path,
    *,
    tokenizer: ChatTokenizer,
    max_tokens: int,
) -> JsonDict:
    sessions: dict[str, set[str]] = {"train": set(), "eval": set()}
    counts: dict[str, int] = {}
    max_seen = 0
    max_system = 0
    forbidden_result_keys = {
        "raw_response",
        "reasoning",
        "workflow_state",
        "reasoning_content",
    }
    for split in ("train", "eval"):
        path = output / "dataset/splits" / split / "tool_loop_sft.jsonl"
        rows = _read_jsonl(path)
        counts[split] = len(rows)
        for index, row in enumerate(rows):
            messages = row.get("messages") or []
            tools = row.get("tools") or []
            if not messages or messages[0].get("role") != "system":
                raise ValueError(f"missing system prompt in {split}:{index}")
            system_chars = len(str(messages[0].get("content") or ""))
            max_system = max(max_system, system_chars)
            if system_chars > MAX_SYSTEM_CHARS:
                raise ValueError(
                    f"system prompt over budget in {split}:{index}"
                )
            actual_tokens = _token_count(tokenizer, messages, tools)
            max_seen = max(max_seen, actual_tokens)
            if actual_tokens > max_tokens:
                raise ValueError(f"token budget exceeded in {split}:{index}")
            if row.get("meta", {}).get("qwen3_14b_tokens") != actual_tokens:
                raise ValueError(f"token receipt mismatch in {split}:{index}")
            pair_error = _tool_pair_error(messages)
            if pair_error:
                raise ValueError(f"{pair_error} in {split}:{index}")
            for message in messages:
                if message.get("role") != "tool":
                    continue
                try:
                    payload = json.loads(str(message.get("content") or "{}"))
                except json.JSONDecodeError as exc:
                    raise ValueError(
                        f"invalid projected tool JSON in {split}:{index}"
                    ) from exc
                if isinstance(
                    payload, dict
                ) and forbidden_result_keys.intersection(payload):
                    raise ValueError(
                        f"unprojected provider payload in {split}:{index}"
                    )
            session = str(row.get("meta", {}).get("session_id") or "")
            if session:
                sessions[split].add(session)
    overlap = sorted(sessions["train"].intersection(sessions["eval"]))
    if overlap:
        raise ValueError(f"session leakage after windowing: {overlap[:3]}")
    return {
        "valid": True,
        "counts": counts,
        "max_qwen3_14b_tokens": max_seen,
        "max_system_chars": max_system,
        "tool_pair_orphans": 0,
        "session_leak_count": 0,
        "unprojected_provider_payloads": 0,
    }


def _refresh_manifests(output: Path, report: JsonDict) -> None:
    dataset = output / "dataset"
    manifest_path = dataset / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["release_id"] = output.name
    manifest["derived_from"] = report["source_release"]
    manifest["tool_loop_projection_schema"] = 2
    manifest["file_sha256"] = {
        path.name: _sha256(path) for path in sorted(dataset.glob("*.jsonl"))
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )

    split_path = dataset / "split_manifest.json"
    split = json.loads(split_path.read_text(encoding="utf-8"))
    for name in ("train", "eval"):
        split["counts"][name]["tool_loop_sft"] = report["splits"][name][
            "output_windows"
        ]
    split["record_count"] = sum(
        sum(values.values()) for values in split["counts"].values()
    )
    split["derived_from"] = report["source_release"]
    split["tool_loop_windowing"] = "source-split-preserving"
    split_path.write_text(
        json.dumps(split, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_release_readme(output: Path, report: JsonDict) -> None:
    prior = (output / "README.md").read_text(encoding="utf-8")
    header = f"""# ChemSmart agent SFT 2026-07-15 v2 pruned

Derived immutably from `{report['source_release']}`. Tool-loop rows use a
<=4,096-character outer system prompt, native per-window tool schemas,
deterministic tool-result projection, and <=8,192-token Qwen3-14B decision
windows. Full-schema prompts and unprovenanced long reasoning are excluded.

See `validation/tool_loop_pruning.json` for exact counts and exclusions.

---

"""
    (output / "README.md").write_text(header + prior, encoding="utf-8")


def _freeze(output: Path, report: JsonDict) -> None:
    for path in output.rglob("*.jsonl"):
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(), 1
        ):
            value = json.loads(line)
            if _contains_secret(value):
                raise ValueError(f"secret in {path}:{line_number}")
            serialized = json.dumps(value, ensure_ascii=False, sort_keys=True)
            if "/Users/" in serialized or "/home/" in serialized:
                raise ValueError(f"home path in {path}:{line_number}")
            if _find_hidden_key(value):
                raise ValueError(f"hidden reasoning in {path}:{line_number}")
    (output / "STATUS").write_text("FROZEN\n", encoding="utf-8")
    files = [
        path
        for path in sorted(output.rglob("*"))
        if path.is_file()
        and path.name != "release.lock.json"
        and "__pycache__" not in path.parts
    ]
    lock = {
        "schema_version": 2,
        "release_id": output.name,
        "status": "FROZEN",
        "derived_from": report["source_release"],
        "tokenizer": report["tokenizer"],
        "max_tokens": report["max_tokens"],
        "files": {
            str(path.relative_to(output)): {
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in files
        },
    }
    (output / "release.lock.json").write_text(
        json.dumps(lock, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _find_hidden_key(value: Any) -> str | None:
    if isinstance(value, dict):
        for key, item in value.items():
            if str(key).lower() in _HIDDEN_KEYS:
                return str(key)
            found = _find_hidden_key(item)
            if found:
                return found
    elif isinstance(value, list):
        for item in value:
            found = _find_hidden_key(item)
            if found:
                return found
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tokenizer", default="Qwen/Qwen3-14B")
    parser.add_argument("--max-tokens", type=int, default=DEFAULT_MAX_TOKENS)
    args = parser.parse_args()
    tokenizer = _load_tokenizer(args.tokenizer)
    report = build_release(
        source=args.source.resolve(),
        output=args.output.resolve(),
        tokenizer=tokenizer,
        tokenizer_name=args.tokenizer,
        max_tokens=args.max_tokens,
    )
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
