"""The reached geometry is a route the menu named for a whole campaign.

``failed_nonconverged_geometry`` told sessions to "restart from the last
geometry the run reached"; ``timeout_terminated`` said the same. Nothing
could: the workspace scan bars the node directories, no tool read a
geometry out of a result, and a producer edge needs a validated node. A
live goal hit that wall twice, wrote "a failed-run end structure cannot
anchor a bound geometry_xyz input in this host" into its own decision
record, and asked for the reached geometry five times in three spellings
before giving up and building starts by hand.

What is pinned here is the invariant, not that case: the host reads the
structure through the shared selector plane, owns the bytes and the
lineage, names the ending its own streams recorded, binds no electronic
state, and refuses only for structural reasons.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.execution import build_reached_geometry

_OUT = Path(__file__).resolve().parents[1] / "data" / "ORCATests" / "outputs"


def _artifact(name: str, kind: str = "orca_output") -> TrustedArtifactRefV1:
    path = _OUT / name
    return TrustedArtifactRefV1(
        artifact_id=f"result.{name}",
        kind=kind,
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        size_bytes=path.stat().st_size,
        path=str(path.resolve()),
        cli_value=str(path.resolve()),
    )


@pytest.mark.capability("tool:bind_reached_geometry")
def test_the_structure_a_run_reached_becomes_a_starting_geometry(tmp_path):
    source = _artifact("sn2_ts.out")
    before = Path(source.path).read_bytes()
    artifact, receipt = build_reached_geometry(
        approved_workspace=tmp_path,
        reached_artifact_id="reached-a",
        result_artifact=source,
        program="orca",
    )

    assert artifact.kind == "geometry_xyz"
    assert artifact.artifact_id == "reached-a"
    assert receipt.source_result_sha256 == source.sha256
    assert receipt.atom_count > 0
    assert receipt.receipt_sha256

    lines = Path(artifact.path).read_text().splitlines()
    assert int(lines[0]) == receipt.atom_count
    assert len(lines) == receipt.atom_count + 2
    # The state is the next question's to answer, so the file says so
    # rather than carrying a charge and a multiplicity nobody chose.
    assert "electronic state deliberately unbound" in lines[1]
    assert "starting structure" in lines[1]

    # Reading a result never edits it: the source bytes are what the
    # engine wrote, before and after.
    assert Path(source.path).read_bytes() == before


@pytest.mark.capability("tool:bind_reached_geometry")
def test_the_host_names_the_ending_its_own_streams_recorded(tmp_path):
    """The ending comes from the workspace's record, never from the model.

    A session cannot talk a result into having succeeded: the lookup is
    keyed by the output bytes' own digest, and where the workspace holds
    no record the receipt says so instead of inventing one.
    """

    source = _artifact("sn2_ts.out")
    _artifact_out, receipt = build_reached_geometry(
        approved_workspace=tmp_path,
        reached_artifact_id="reached-b",
        result_artifact=source,
        program="orca",
        run_evidence_root=tmp_path,
    )
    assert receipt.recorded_terminal_state == ""
    assert receipt.recorded_node_id == ""
    header = Path(_artifact_out.path).read_text().splitlines()[1]
    assert "no ending recorded in this workspace" in header


@pytest.mark.capability("tool:bind_reached_geometry")
def test_the_refusals_are_structural_and_name_the_shape(tmp_path):
    with pytest.raises(ContractError, match="geometry_xyz"):
        build_reached_geometry(
            approved_workspace=tmp_path,
            reached_artifact_id="reached-c",
            result_artifact=_artifact("sn2_ts.out", kind="geometry_xyz"),
            program="orca",
        )
    with pytest.raises(ContractError, match="no result reader"):
        build_reached_geometry(
            approved_workspace=tmp_path,
            reached_artifact_id="reached-d",
            result_artifact=_artifact("sn2_ts.out"),
            program="nosuchprogram",
        )


@pytest.mark.capability("tool:bind_reached_geometry")
def test_carrying_a_geometry_forward_upgrades_nothing_about_the_source():
    """The producer-edge ban is what the charter actually forbids.

    This tool mints a new starting structure by an explicit session act.
    It must not become a side door by which a failed node satisfies a
    data edge, so the rule that keeps a failed node out of geometry
    handoff is checked here against the states the derivation can emit.
    """

    from chemsmart.agent.terminal_states import (
        NODE_TERMINAL_STATES,
        REPAIRABLE_NODE_STATES,
    )

    assert REPAIRABLE_NODE_STATES <= set(NODE_TERMINAL_STATES)
    assert "validated" not in REPAIRABLE_NODE_STATES, (
        "every ending this tool answers is a non-validated one, so none "
        "of them may satisfy a producer edge"
    )


@pytest.mark.capability("rule:recovery.restart_from_what_it_reached")
def test_the_menu_and_the_guide_name_a_route_that_exists():
    """A route named in host prose is a capability claim.

    Three repair-menu entries named restarts nothing could perform. The
    pin is that every tool the menu, the rules and the recovery guide
    name is a tool the model is actually handed.
    """

    import re

    from chemsmart.agent.driver import REPAIR_MENU
    from chemsmart.agent.guides import GUIDES
    from chemsmart.agent.rules import POLICY_RULES
    from chemsmart.agent.tool_specs import (
        build_command_compiled_tool_surface,
    )

    surface = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface(
            guides=tuple(guide.guide_id for guide in GUIDES)
        ).tool_definitions
    }
    assert "bind_reached_geometry" in surface

    prose = list(REPAIR_MENU.values())
    prose += [rule.text for rule in POLICY_RULES]
    prose += [guide.body for guide in GUIDES]

    # Tools are named by a verb; the rest of the host's snake_case
    # vocabulary (terminal states, selectors, artifact kinds) is not.
    # Matching on the verb rather than on membership is what lets this
    # fail: a sentence naming an action nobody implemented is caught,
    # which is precisely the defect that shipped.
    verbs = (
        "amend_",
        "append_",
        "bind_",
        "characterise_",
        "compile_",
        "compose_",
        "declare_",
        "derive_",
        "displace_",
        "edit_",
        "evaluate_",
        "extract_",
        "inspect_",
        "open_",
        "plan_",
        "project_",
        "record_",
    )
    named = {
        word
        for text in prose
        for word in re.findall(r"\b[a-z][a-z0-9_]*\b", text)
        if word.startswith(verbs) and "_" in word
    }
    assert len(named) > 10, "the guard must actually be reading the prose"
    missing = named - surface
    assert not missing, f"host prose names tools nobody is handed: {missing}"
