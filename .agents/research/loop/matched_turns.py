"""Matched turns: an archived prefix through a chosen tree, then real turns.

One sample is one live session on a fresh workspace: the archived
assistant turns 1..K are served verbatim to the host of the tree this
process imported (host-computed digests in their arguments translated to
this host's), so every reply the model reads is this tree's; then up to P
turns go to the real provider through the adapter the session built, with
the credential its own lease resolved; then the transport answers
"Ending." until the session ends. Two trees, or two provider profiles,
give two arms whose only difference is the host or the model.

The tree is whatever ``chemsmart`` PYTHONPATH resolves -- this file never
edits ``sys.path`` -- and every sample row names it (``chemsmart.__file__``
and a digest of the package's Python sources); ``--expect-tree`` refuses
to run when the import is not the tree named. HOME is fenced to
``--home`` before ``chemsmart`` is imported, so nothing a session writes
or reads by default lives in the developer's home; the credential file is
passed through by path (``--secret-file``) and never opened here.

    PYTHONPATH=<tree> python .agents/research/loop/matched_turns.py \\
        --home DIR --out DIR --arm NAME --task FILE --inputs FILE [...] \\
        --envelope FILE --provider-config FILE --profile NAME \\
        [--prefix TRANSCRIPT --cut K] [--probe-turns P] [--samples N] \\
        [--goal GOAL_ID --granted-by LABEL --max-revisions R] \\
        [--secret-file PATH | --stub] [--expect-tree SUBSTRING]

Setup: without ``--goal`` the session is ``run_live_agent_session`` over
the envelope (a planning session, as R10 Q26's counterfactual turn ran);
with ``--goal`` it is the goal driver's first cycle, whose planning
session is the sample and after which the driver is stopped -- nothing is
decided, approved or launched either way.

What a row records, all read from the host's own stream and transcript:
per real turn the requested and observed model, finish reason, tool calls
and token usage; the calls in view (``callable`` of every
``exposure_planned`` event, and of every notice that carries it); the
notices delivered; every act after the cut with its reply status; and
prefix faithfulness -- each replayed call's reply status beside the
archived one. A sample is INFRA (reported, never counted, never re-rolled)
when no real turn was observed, the transport failed, a turn deadline was
exceeded, or the replayed prefix did not reproduce the archived statuses.
``--stub`` answers every probe turn with a fixed text and contacts no
provider: a free check of the instrument, never a sample. The session's
stream and transcript are copied beside the row and the workspace is
deleted.

Adapted from R10 Q26 (``cf_probe.py``: digests paired by JSON path) and
Q32 (``matched_turns.py``: digests paired by position, the goal driver
stopped after the session).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any

HEX64 = re.compile(r"(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])")
STUB_EXECUTABLE = (
    "#!/bin/sh\n# ChemSmart agent-harness DISCOVERY STUB\nexit 1\n"
)
#: Executables a fenced server profile declares as discovery stubs, so a
#: program the envelope allows is discoverable on a machine without it.
STUB_PROGRAMS = {"orca": "orca", "gaussian": "g16", "xtb": "xtb"}
ENDING = "Ending."


def _arguments(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--home", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--arm", required=True)
    parser.add_argument("--task", required=True, type=Path)
    parser.add_argument("--inputs", nargs="+", default=[], type=Path)
    parser.add_argument("--envelope", required=True, type=Path)
    parser.add_argument("--provider-config", required=True, type=Path)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--prefix", type=Path)
    parser.add_argument("--cut", type=int, default=0)
    parser.add_argument("--probe-turns", type=int, default=3)
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--first-sample", type=int, default=0)
    parser.add_argument("--goal")
    parser.add_argument("--granted-by", default="")
    parser.add_argument("--max-revisions", type=int, default=5)
    parser.add_argument("--exposure-mode", default="")
    parser.add_argument("--secret-file", type=Path)
    parser.add_argument("--stub", action="store_true")
    parser.add_argument("--expect-tree", default="")
    args = parser.parse_args(argv)
    if args.cut and args.prefix is None:
        parser.error("--cut needs --prefix")
    if not args.stub and args.secret_file is None:
        parser.error("a real sample needs --secret-file (passed, never read)")
    if args.goal and not args.granted_by:
        parser.error("--goal needs --granted-by (a delegated label)")
    return args


def _fence(home: Path, envelope: Path) -> Path:
    """Fence HOME before chemsmart is imported; declare the envelope's
    programs through discovery stubs so they are discoverable here."""

    home = home.resolve()
    config = home / ".chemsmart"
    (config / "server").mkdir(parents=True, exist_ok=True)
    os.environ["HOME"] = str(home)
    os.environ["CHEMSMART_CONFIG_DIR"] = str(config)
    os.environ.pop("CHEMSMART_AGENT_SERVER", None)
    os.environ.pop("CHEMSMART_AGENT_CONFIG", None)
    os.environ.pop("CHEMSMART_AGENT_KEYS", None)
    text = envelope.read_text(encoding="utf-8")
    blocks = []
    for program, executable in sorted(STUB_PROGRAMS.items()):
        if not re.search(rf"^\s*{program}\s*:", text, re.MULTILINE):
            continue
        folder = home / "stubs" / program
        folder.mkdir(parents=True, exist_ok=True)
        stub = folder / executable
        stub.write_text(STUB_EXECUTABLE)
        stub.chmod(0o755)
        blocks.append(
            f"{program.upper()}:\n    EXEFOLDER: {folder}\n"
            "    LOCAL_RUN: true\n    SCRATCH: false\n"
        )
    (config / "server" / "local.yaml").write_text(
        "SERVER:\n    SCHEDULER: null\n    NUM_CORES: 4\n    MEM_GB: 8\n"
        "    NUM_HOURS: 1\n" + "".join(blocks),
        encoding="utf-8",
    )
    placeholder = home / "stub-keys.env"
    placeholder.write_text("ALIBABA_TOKEN_PLAN_KEY=sk-sp-stub-not-a-key\n")
    placeholder.chmod(0o600)
    return placeholder


def _tree_record(expect: str) -> dict[str, Any]:
    import chemsmart

    package = Path(chemsmart.__file__).resolve().parent
    digest = hashlib.sha256()
    for path in sorted(package.rglob("*.py")):
        digest.update(str(path.relative_to(package)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
    record = {
        "chemsmart_file": str(Path(chemsmart.__file__).resolve()),
        "package_py_sha256": digest.hexdigest(),
    }
    if expect and expect not in record["chemsmart_file"]:
        raise SystemExit(
            f"chemsmart imported from {record['chemsmart_file']}, not a "
            f"tree matching {expect!r}; set PYTHONPATH and run by path"
        )
    return record


def _messages(path: Path) -> list[dict[str, Any]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(document, dict):
        return list(document.get("transcript") or [])
    return list(document)


def _archived_turns(path: Path | None) -> tuple[list[dict], dict[str, str]]:
    if path is None:
        return [], {}
    messages = _messages(path)
    replies = {
        str(m.get("tool_call_id")): str(m.get("content") or "")
        for m in messages
        if m.get("role") == "tool"
    }
    turns = [m for m in messages if m.get("role") == "assistant"]
    return turns, replies


def _parsed(text: str) -> Any:
    try:
        return json.loads(text)
    except (TypeError, ValueError):
        return None


def _pair_by_path(archived: Any, replayed: Any, found: dict[str, str]) -> None:
    if isinstance(archived, dict) and isinstance(replayed, dict):
        for key in archived.keys() & replayed.keys():
            _pair_by_path(archived[key], replayed[key], found)
    elif isinstance(archived, list) and isinstance(replayed, list):
        for left, right in zip(archived, replayed):
            _pair_by_path(left, right, found)
    elif isinstance(archived, str) and isinstance(replayed, str):
        left, right = HEX64.findall(archived), HEX64.findall(replayed)
        if len(left) == len(right):
            for old, new in zip(left, right):
                if old != new:
                    found.setdefault(old, new)


def _learn(archived: str, replayed: str, found: dict[str, str]) -> None:
    """Pair host digests at the same place in two replies to one call."""

    before = len(found)
    _pair_by_path(_parsed(archived), _parsed(replayed), found)
    if len(found) == before:
        left, right = HEX64.findall(archived), HEX64.findall(replayed)
        if len(left) == len(right):
            for old, new in zip(left, right):
                if old != new:
                    found.setdefault(old, new)


def _statuses(text: str) -> tuple[str, str]:
    reply = _parsed(text)
    if not isinstance(reply, dict):
        return "", ""
    inner = reply.get("result")
    return (
        str(reply.get("status") or ""),
        str(inner.get("status") or "") if isinstance(inner, dict) else "",
    )


def _response(message: dict[str, Any], model: str, ordinal: int) -> dict:
    return {
        "id": f"matched-{ordinal}",
        "object": "chat.completion",
        "model": model,
        "choices": [
            {
                "index": 0,
                "finish_reason": (
                    "tool_calls" if message.get("tool_calls") else "stop"
                ),
                "message": message,
            }
        ],
        "usage": {
            "prompt_tokens": 1,
            "completion_tokens": 1,
            "total_tokens": 2,
            "completion_tokens_details": {"reasoning_tokens": 0},
        },
    }


def _hybrid(real_class, turns, replies, cut, probe_turns, stub, record):
    """The transport: archived turns, then real ones, then an ending."""

    state = {"served": 0, "pending": [], "map": {}}
    record["prefix_calls"] = []
    record["probe"] = []

    class Hybrid:
        def __init__(self, **kwargs):
            self.real = None if stub else real_class(**kwargs)

        def set_timeout_seconds(self, value):
            setter = getattr(self.real, "set_timeout_seconds", None)
            if callable(setter):
                setter(value)

        def close(self):
            if self.real is not None:
                self.real.close()

        def public_deadline_record(self):
            getter = getattr(self.real, "public_deadline_record", None)
            return getter() if callable(getter) else {}

        def __call__(self, payload):
            seen = {
                str(m.get("tool_call_id")): str(m.get("content") or "")
                for m in payload.get("messages") or ()
                if m.get("role") == "tool"
            }
            for call_id, name in state["pending"]:
                replayed = seen.get(call_id, "")
                archived = replies.get(call_id, "")
                _learn(archived, replayed, state["map"])
                record["prefix_calls"].append(
                    {
                        "call_id": call_id,
                        "tool": name,
                        "archived_status": _statuses(archived),
                        "replayed_status": _statuses(replayed),
                    }
                )
            state["pending"] = []
            index = state["served"]
            state["served"] += 1
            if index < cut:
                archived_turn = turns[index]
                message: dict[str, Any] = {
                    "role": "assistant",
                    "content": str(archived_turn.get("content") or ""),
                    "reasoning_content": "",
                }
                calls = []
                for call in archived_turn.get("tool_calls") or ():
                    call = json.loads(json.dumps(call))
                    function = call.setdefault("function", {})
                    text = str(function.get("arguments") or "")
                    for old, new in state["map"].items():
                        text = text.replace(old, new)
                    function["arguments"] = text
                    state["pending"].append(
                        (str(call.get("id")), str(function.get("name")))
                    )
                    calls.append(call)
                if calls:
                    message["tool_calls"] = calls
                elif not message["content"]:
                    message["content"] = "..."
                # The protocol checks that a response names the model it
                # was asked for, so a served turn echoes it; replayed and
                # real turns are told apart by position, never by name.
                return _response(message, payload.get("model"), index)
            if index < cut + probe_turns:
                if self.real is None:
                    response = _response(
                        {"role": "assistant", "content": "stub probe"},
                        payload.get("model"),
                        index,
                    )
                else:
                    response = self.real(payload)
                choice = (response.get("choices") or [{}])[0]
                message = choice.get("message") or {}
                record["probe"].append(
                    {
                        "turn": index - cut + 1,
                        "request_messages": len(payload.get("messages") or ()),
                        "finish_reason": choice.get("finish_reason"),
                        "response_model": response.get("model"),
                        "usage": response.get("usage"),
                    }
                )
                return response
            return _response(
                {"role": "assistant", "content": ENDING},
                payload.get("model"),
                index,
            )

    return Hybrid, state


def _read_session(
    workspace: Path, cut: int, served_real: int, sample_dir: Path
) -> dict:
    """What the host recorded, read back from the session's own files.

    Observed provider turns are the replayed prefix (the first ``cut``),
    then the real turns (one per response the provider served), then the
    transport's closing turns.
    """

    runs = sorted((workspace / ".chemsmart-agent" / "runs").glob("live-*"))
    if not runs:
        return {"session_found": False}
    session = runs[-1]
    events = [
        json.loads(line)
        for line in (session / "events.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    sample_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(session / "events.jsonl", sample_dir / "events.jsonl")
    transcripts = sorted(session.glob("public-transcript-*.json"))
    messages: list[dict] = []
    if transcripts:
        shutil.copyfile(transcripts[-1], sample_dir / "transcript.json")
        messages = _messages(transcripts[-1])
    turns = []
    views = []
    notices = []
    for event in events:
        kind = event.get("kind")
        payload = event.get("payload") or {}
        if kind == "provider_turn_observed":
            turns.append(
                {
                    "requested_model": payload.get("requested_model"),
                    "observed_model": payload.get("observed_model"),
                    "finish_reason": payload.get("finish_reason"),
                    "tool_calls_present": payload.get("tool_calls_present"),
                    "reasoning_tokens": payload.get("reasoning_tokens"),
                }
            )
        elif kind == "exposure_planned":
            views.append(
                {
                    "turn_index": len(turns),
                    "exposure_sha256": payload.get("exposure_sha256"),
                    "callable": payload.get("callable")
                    or payload.get("tools"),
                }
            )
        elif kind in (
            "execution_wave_decision_pending",
            "termination_notice_delivered",
            "host_context_reinjected",
        ):
            notices.append(
                {
                    "kind": kind,
                    "turn_index": len(turns),
                    "content_sha256": payload.get("content_sha256"),
                    "named_tools_in_view": payload.get("named_tools_in_view"),
                    "callable_recorded": "callable" in payload,
                }
            )
    for position, turn in enumerate(turns):
        turn["origin"] = (
            "replayed"
            if position < cut
            else ("real" if position < cut + served_real else "closing")
        )
    real = [t for t in turns if t["origin"] == "real"]
    replies = {
        str(m.get("tool_call_id")): str(m.get("content") or "")
        for m in messages
        if m.get("role") == "tool"
    }
    assistants = [m for m in messages if m.get("role") == "assistant"]
    acts = []
    texts = []
    for position, turn in enumerate(assistants[cut:], start=1):
        content = str(turn.get("content") or "")
        if content and content != ENDING:
            texts.append({"turn": position, "content": content})
        for call in turn.get("tool_calls") or ():
            function = call.get("function") or {}
            reply = replies.get(str(call.get("id")), "")
            acts.append(
                {
                    "turn": position,
                    "name": function.get("name"),
                    "arguments": str(function.get("arguments") or ""),
                    "status": _statuses(reply),
                }
            )
    blocked = [
        event.get("payload") or {}
        for event in events
        if event.get("kind") in ("runtime_terminated", "turn_blocked")
    ]
    return {
        "session_found": True,
        "provider_turns": turns,
        "real_turns": len(real),
        "models_observed": sorted({str(t["observed_model"]) for t in real}),
        "models_requested": sorted({str(t["requested_model"]) for t in real}),
        "views": views,
        "notices": notices,
        "acts_after_cut": acts,
        "texts_after_cut": texts,
        "deadline_exceeded": any(
            "turn_deadline_exceeded" in json.dumps(item) for item in blocked
        ),
        "terminated": [
            str(item.get("terminal_state") or item.get("reason") or "")
            for item in blocked
        ],
    }


def _local_envelope(source: Path, out: Path) -> Path:
    """The archived envelope with its engine scratch moved under ``out``.

    An archived envelope names the scratch of the machine it ran on; a
    planning session launches nothing, and the path is the one line that
    must exist here.
    """

    out.mkdir(parents=True, exist_ok=True)
    lines = [
        (
            f"scratch_root: {out.resolve() / 'engine-scratch'}"
            if line.startswith("scratch_root:")
            else line
        )
        for line in source.read_text(encoding="utf-8").splitlines()
    ]
    target = out / "envelope.yaml"
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def main(argv: list[str] | None = None) -> None:
    args = _arguments(argv)
    args.envelope = _local_envelope(args.envelope, args.out)
    placeholder = _fence(args.home, args.envelope)
    tree = _tree_record(args.expect_tree)
    print("chemsmart from", tree["chemsmart_file"], flush=True)

    import chemsmart.agent.runtime.alibaba as alibaba
    from chemsmart.agent.live_session import run_live_agent_session

    real_class = alibaba.AlibabaTokenPlanHttpsTransport
    turns, replies = _archived_turns(args.prefix)
    if args.cut > len(turns):
        raise SystemExit(f"--cut {args.cut} exceeds {len(turns)} turns")
    task = args.task.read_text(encoding="utf-8").strip()
    secret = placeholder if args.stub else args.secret_file
    args.out.mkdir(parents=True, exist_ok=True)

    class _StopAfterSession(Exception):
        pass

    for sample in range(args.first_sample, args.first_sample + args.samples):
        name = f"{args.arm}-{sample:02d}"
        root = args.out / "work" / name
        if root.exists():
            shutil.rmtree(root)
        workspace = root / "ws"
        workspace.mkdir(parents=True)
        for item in args.inputs:
            shutil.copyfile(item, workspace / item.name)
        record: dict[str, Any] = {
            "arm": args.arm,
            "sample": sample,
            "tree": tree,
            "profile": args.profile,
            "provider_config_sha256": hashlib.sha256(
                args.provider_config.read_bytes()
            ).hexdigest(),
            "envelope_sha256": hashlib.sha256(
                args.envelope.read_bytes()
            ).hexdigest(),
            "task_sha256": hashlib.sha256(args.task.read_bytes()).hexdigest(),
            "prefix_sha256": (
                hashlib.sha256(args.prefix.read_bytes()).hexdigest()
                if args.prefix
                else ""
            ),
            "cut": args.cut,
            "probe_turns": args.probe_turns,
            "goal": args.goal or "",
            "stub": bool(args.stub),
        }
        hybrid, state = _hybrid(
            real_class,
            turns,
            replies,
            args.cut,
            args.probe_turns,
            args.stub,
            record,
        )
        alibaba.AlibabaTokenPlanHttpsTransport = hybrid
        started = time.time()
        common = {
            "provider": args.profile,
            "provider_config_file": args.provider_config,
            "secret_file": secret,
            "workspace": workspace,
            "execution_enabled": False,
            "approval_file": None,
            "execution_envelope_file": args.envelope,
        }
        try:
            if args.goal:
                from chemsmart.agent.driver import GoalDriver

                def plan_session(**kwargs):
                    kwargs["secret_file"] = secret
                    if args.exposure_mode:
                        kwargs["exposure_mode"] = args.exposure_mode
                    result = run_live_agent_session(**kwargs)
                    record["terminal_state"] = result.terminal_state
                    raise _StopAfterSession

                driver = GoalDriver(
                    task=task,
                    workspace=workspace,
                    execution_envelope_file=args.envelope,
                    goal_id=args.goal,
                    granted_by=args.granted_by,
                    max_revisions=args.max_revisions,
                    provider=args.profile,
                    provider_config_file=args.provider_config,
                    plan_session=plan_session,
                )
                driver.run()
            else:
                result = run_live_agent_session(
                    task=task,
                    exposure_mode=args.exposure_mode,
                    **common,
                )
                record["terminal_state"] = result.terminal_state
        except _StopAfterSession:
            pass
        except Exception as exc:  # reported, never hidden
            record["error"] = f"{type(exc).__name__}: {exc}"[:2000]
        finally:
            alibaba.AlibabaTokenPlanHttpsTransport = real_class
        record["seconds"] = round(time.time() - started, 1)
        record["digests_translated"] = len(state["map"])
        record.update(
            _read_session(
                workspace,
                args.cut,
                0 if args.stub else len(record["probe"]),
                args.out / name,
            )
        )
        record["prefix_faithful"] = all(
            call["archived_status"] == call["replayed_status"]
            for call in record["prefix_calls"]
        ) and len(record["prefix_calls"]) == sum(
            len(t.get("tool_calls") or ()) for t in turns[: args.cut]
        )
        reasons = []
        if not args.stub and not record.get("real_turns"):
            reasons.append("no real provider turn")
        if record.get("deadline_exceeded"):
            reasons.append("turn deadline exceeded")
        if "error" in record:
            reasons.append("session raised")
        if not record["prefix_faithful"]:
            reasons.append("prefix did not reproduce the archived statuses")
        record["infra"] = reasons
        with (args.out / "samples.jsonl").open("a", encoding="utf-8") as out:
            out.write(json.dumps(record, default=str) + "\n")
        print(
            name,
            "real_turns",
            record.get("real_turns"),
            "models",
            record.get("models_observed"),
            "acts",
            [a["name"] for a in record.get("acts_after_cut", [])][:12],
            "infra",
            reasons,
            flush=True,
        )
        shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
