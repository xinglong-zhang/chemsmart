"""Tests for the v8 spec-emission adapter and its wiring into SynthesisSession."""
import json

from chemsmart.agent import v8_adapter
from chemsmart.agent import synthesis as S


def test_postprocess_strips_invalid_and_repairs_ts_route():
    # `freq` is not a valid nci key -> stripped; a mangled TS route is repaired to the canonical list.
    spec = {
        "intent": "workflow",
        "jobs": [
            {"id": 1, "kind": "gaussian.nci", "file": "a.xyz", "charge": 0, "mult": 1,
             "settings": {"freq": True}},
            {"id": 2, "kind": "gaussian.ts", "file": "b.xyz", "charge": 0, "mult": 1,
             "settings": {"ts": "calcfc,noeigentest", "freq": True}},
        ],
    }
    out = v8_adapter.postprocess(spec)
    assert "settings" not in out["jobs"][0]  # freq stripped from nci
    assert out["jobs"][1]["settings"]["additional_opt_options_in_route"] == ["ts", "calcfc", "noeigentest"]


def test_adapt_renders_valid_chemsmart_command():
    spec = {"intent": "workflow", "jobs": [
        {"id": 1, "kind": "gaussian.opt", "file": "mol/water.xyz", "charge": 0, "mult": 1,
         "execution": "submit", "server": "cpuq"}]}
    out = v8_adapter.adapt(json.dumps(spec))
    assert out["intent"] == "workflow"
    assert out["commands"] == ["chemsmart sub -s cpuq gaussian opt -f mol/water.xyz -c 0 -m 1"]
    assert out["valid"] is True, out["errors"]


def test_adapt_chain_renders_ordered_commands():
    spec = {"intent": "workflow", "jobs": [
        {"id": 1, "kind": "gaussian.opt", "file": "m.xyz", "charge": 0, "mult": 1},
        {"id": 2, "kind": "gaussian.freq", "geom_from": 1, "charge": 0, "mult": 1}]}
    out = v8_adapter.adapt(json.dumps(spec))
    assert len(out["commands"]) == 2
    assert out["valid"] is True, out["errors"]


def test_adapt_non_workflow_has_message_no_command():
    out = v8_adapter.adapt(json.dumps({"intent": "decline", "message": "need atom indices"}))
    assert out["intent"] == "decline"
    assert out["commands"] == []
    assert out["message"] == "need atom indices"


def test_synthesis_routes_v8_spec():
    wf = {"intent": "workflow", "jobs": [
        {"id": 1, "kind": "gaussian.opt", "file": "a.xyz", "charge": 0, "mult": 1}]}
    assert S._is_v8_spec(wf) is True
    assert S._is_v8_spec({"status": "ready", "command": "chemsmart run gaussian opt -f a.xyz"}) is False
    res = S._normalize_v8_spec(wf)
    assert res["status"] == "ready"
    assert res["command"].startswith("chemsmart ")


def test_synthesis_v8_decline_is_infeasible():
    res = S._normalize_v8_spec({"intent": "decline", "message": "missing fragment indices"})
    assert res["status"] == "infeasible"
    assert res["command"] == ""
    assert "fragment" in res["explanation"]


def test_native_multiturn_history_is_replayed():
    """A follow-up turn must reach the model as native [system, user, assistant, user] chat
    history carrying the prior raw SPEC — so the model can edit its own previous spec."""
    from chemsmart.agent.synthesis import SynthesisSession
    from chemsmart.agent.services.conversation_memory import ConversationMemory

    spec1 = json.dumps({"intent": "workflow", "jobs": [
        {"id": 1, "kind": "gaussian.opt", "file": "m.xyz", "charge": 0, "mult": 1}]})
    spec2 = json.dumps({"intent": "workflow", "jobs": [
        {"id": 1, "kind": "gaussian.opt", "file": "m.xyz", "charge": -1, "mult": 1}]})
    seen = {}

    class MockProvider:
        def __init__(self):
            self.calls = 0

        def chat(self, messages):
            seen["last"] = messages
            self.calls += 1
            return spec1 if self.calls == 1 else spec2

    sess = SynthesisSession.__new__(SynthesisSession)
    sess.provider = MockProvider()
    sess.max_retries = 0
    sess.debug = False
    sess.schema = {"command": "chemsmart", "subcommands": {}}  # stub; only json-serialized into the prompt
    sess.memory = ConversationMemory()

    # turn 1
    r1 = sess.synthesize("optimize m.xyz")
    assert r1["status"] == "ready"
    assert seen["last"][-1]["role"] == "user"  # first call: just system + current user
    assert [m["role"] for m in seen["last"]] == ["system", "user"]
    sess._remember_turn("optimize m.xyz", r1["command"], assistant_message=r1["raw_response"])

    # turn 2: the follow-up must arrive as native history with the prior SPEC as an assistant turn
    r2 = sess.synthesize("change the charge to -1")
    assert [m["role"] for m in seen["last"]] == ["system", "user", "assistant", "user"]
    assert seen["last"][2]["content"] == spec1          # the prior raw SPEC, verbatim
    assert seen["last"][3]["content"] == "change the charge to -1"
    assert r2["status"] == "ready" and "-1" in r2["command"]
