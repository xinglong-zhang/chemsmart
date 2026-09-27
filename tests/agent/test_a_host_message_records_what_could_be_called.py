"""A message the host places in a goal session records what could be called.

A session under a goal is told things mid-session -- a pending execution
decision, what remains undelivered before it ends, the goal's terms again
at cadence -- and each of those messages names acts. Whether the act a
message names was one the model could call when it read it is a fact about
the host, and it was recorded nowhere: R10 Q32 learned that 5 of 6 apparent
prose decisions were made with the named tool not in view only by replaying
streams and pairing each notice with the last exposure record before it.

Driven through the goal driver's first cycle and the live session with only
the provider transport replaced by a script, so the notice asserted on is
the one a model reads. The oracle is the stream's own exposure record, and
the catalogue names the delivered text mentions.
"""

from __future__ import annotations

import json
import re

import pytest

pytestmark = pytest.mark.capability("rule:wake.termination_notice")

WATER = (
    "3\nwater\nO 0.0 0.0 0.1173\nH 0.0 0.7572 -0.4692\n"
    "H 0.0 -0.7572 -0.4692\n"
)
STUB = "#!/bin/sh\n# ChemSmart agent-harness DISCOVERY STUB\nexit 1\n"
OBSERVABLE = {
    "observable_id": "water-total-energy",
    "unit": "hartree",
    "meaning": "Total electronic energy of the water molecule.",
    "role": "requested",
}


class _Stop(Exception):
    """Raised after the session, so the driver goes no further."""


def _fenced(monkeypatch, tmp_path):
    home = tmp_path / "home"
    config = home / ".chemsmart"
    (config / "server").mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config))
    monkeypatch.delenv("CHEMSMART_AGENT_SERVER", raising=False)
    folder = home / "stubs" / "orca"
    folder.mkdir(parents=True)
    (folder / "orca").write_text(STUB)
    (folder / "orca").chmod(0o755)
    (config / "server" / "local.yaml").write_text(
        "SERVER:\n    SCHEDULER: null\n    NUM_CORES: 4\n    MEM_GB: 8\n"
        f"    NUM_HOURS: 1\nORCA:\n    EXEFOLDER: {folder}\n"
        "    LOCAL_RUN: true\n    SCRATCH: false\n"
    )
    agent = home / "agent.yaml"
    agent.write_text(
        "active: alibaba-token-plan\nfallback: []\nproviders:\n"
        "  alibaba-token-plan:\n    type: openai\n"
        "    api_key_env: ALIBABA_TOKEN_PLAN_KEY\n"
        "    model: deepseek-v4-flash-0731\n    context_tokens: 1000000\n"
        "    max_output_tokens: 262144\n"
        "    base_url: https://token-plan.ap-southeast-1.maas.aliyuncs.com"
        "/compatible-mode/v1\n"
        "    reasoning_effort: max\n    preserve_thinking: true\n"
    )
    keys = home / "keys.env"
    keys.write_text("ALIBABA_TOKEN_PLAN_KEY=sk-sp-placeholder-not-a-key\n")
    keys.chmod(0o600)
    envelope = tmp_path / "envelope.yaml"
    envelope.write_text(
        "schema_version: chemsmart.bounded-execution-envelope.v1\n"
        "mode: bounded-local\nallowed_program_engines:\n  orca: [cpu]\n"
        "resources:\n  execution_target: run\n  cores: 4\n  memory_gb: 8\n"
        "  gpu_count: 0\n  scratch_policy: server\n"
        "  node_timeout_seconds: 3000\nepisode_wall_time_seconds: 5400\n"
        "postprocess_reserve_seconds: 600\nmax_engine_calls: 4\n"
        f"scratch_root: {tmp_path / 'scratch'}\n"
    )
    return agent, keys, envelope


def _scripted_transport(turns):
    """Plays the model: each turn is a list of calls or a closing text.

    A call the host answers by loading its schema is issued again, as a
    model under host search does.
    """

    state = {"pending": [], "turn": 0}

    class Transport:
        def __init__(self, **_kwargs):
            pass

        def set_timeout_seconds(self, _value):
            pass

        def close(self):
            pass

        def __call__(self, payload):
            replies = {
                m.get("tool_call_id"): m.get("content")
                for m in payload.get("messages") or ()
                if m.get("role") == "tool"
            }
            again = []
            for call_id, name, args in state["pending"]:
                try:
                    reply = json.loads(replies.get(call_id) or "{}")
                except ValueError:
                    reply = {}
                if reply.get("status") == "schema_loaded":
                    again.append((name, args))
            state["pending"] = []
            state["turn"] += 1
            step = again or (turns.pop(0) if turns else "Done.")
            if isinstance(step, str):
                message = {"role": "assistant", "content": step}
            else:
                message = {
                    "role": "assistant",
                    "content": "",
                    "reasoning_content": "",
                    "tool_calls": [],
                }
                for number, (name, args) in enumerate(step):
                    call_id = f"call-{state['turn']}-{number}"
                    state["pending"].append((call_id, name, args))
                    message["tool_calls"].append(
                        {
                            "id": call_id,
                            "type": "function",
                            "function": {
                                "name": name,
                                "arguments": json.dumps(args),
                            },
                        }
                    )
            return {
                "id": f"scripted-{state['turn']}",
                "model": payload.get("model"),
                "choices": [
                    {
                        "finish_reason": (
                            "tool_calls"
                            if message.get("tool_calls")
                            else "stop"
                        ),
                        "message": message,
                    }
                ],
                "usage": {"prompt_tokens": 1, "completion_tokens": 1},
            }

    return Transport


def _events(workspace):
    streams = sorted(
        (workspace / ".chemsmart-agent" / "runs").glob("live-*/events.jsonl")
    )
    assert streams, "the session left no stream"
    return [
        json.loads(line)
        for line in streams[-1].read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _delivered_host_messages(workspace):
    transcripts = sorted(
        (workspace / ".chemsmart-agent" / "runs").glob(
            "live-*/public-transcript-*.json"
        )
    )
    assert transcripts, "the session left no public transcript"
    document = json.loads(transcripts[-1].read_text(encoding="utf-8"))
    messages = (
        document.get("transcript") if isinstance(document, dict) else document
    )
    return [
        str(m.get("content") or "") for m in messages if m["role"] == "user"
    ]


def test_a_notice_under_a_goal_records_the_calls_in_view(
    monkeypatch, tmp_path
):
    import chemsmart.agent.runtime.alibaba as alibaba
    from chemsmart.agent.catalogue import build_tool_catalogue
    from chemsmart.agent.driver import GoalDriver
    from chemsmart.agent.live_session import run_live_agent_session

    agent, keys, envelope = _fenced(monkeypatch, tmp_path)
    workspace = tmp_path / "ws"
    workspace.mkdir()
    (workspace / "water.xyz").write_text(WATER)
    monkeypatch.setattr(
        alibaba,
        "AlibabaTokenPlanHttpsTransport",
        _scripted_transport(
            [
                [
                    (
                        "declare_requested_observable",
                        {"observables": [OBSERVABLE]},
                    )
                ],
                "I will stop here.",
                "Still done.",
            ]
        ),
    )

    def plan_session(**kwargs):
        run_live_agent_session(secret_file=keys, **kwargs)
        raise _Stop

    driver = GoalDriver(
        task="What is the total electronic energy of water?",
        workspace=workspace,
        execution_envelope_file=envelope,
        goal_id="goal-notice-view",
        granted_by="test-delegated",
        provider="alibaba-token-plan",
        provider_config_file=agent,
        plan_session=plan_session,
    )
    with pytest.raises(_Stop):
        driver.run()

    events = _events(workspace)
    notices = [
        (index, event)
        for index, event in enumerate(events)
        if event["kind"] == "termination_notice_delivered"
    ]
    assert len(notices) == 1, "the goal session was not told what remains"
    index, notice = notices[0]
    in_force = [
        event
        for event in events[:index]
        if event["kind"] == "exposure_planned"
    ][-1]["payload"]
    payload = notice["payload"]
    assert (
        payload.get("exposure_sha256") == in_force["exposure_sha256"]
    ), "the notice does not say which exposure was in force when it was shown"
    assert (
        payload.get("callable") == in_force["callable"]
    ), "the notice does not record the calls the model could make"
    text = next(
        message
        for message in _delivered_host_messages(workspace)
        if message.startswith("Informational, from the host, once")
    )
    mentioned = set(re.findall(r"[a-z][a-z0-9_]*", text))
    named = {
        name for name in build_tool_catalogue().names() if name in mentioned
    }
    assert payload.get("named_tools_in_view") == {
        name: name in in_force["callable"] for name in sorted(named)
    }, "the calls the notice names are not each marked in or out of view"
