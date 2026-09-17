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

import numpy as np
import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.execution import build_reached_geometry

_OUT = Path(__file__).resolve().parents[1] / "data" / "ORCATests" / "outputs"
_XTB_OPT = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "XTBTests"
    / "outputs"
    / "co2_ohess"
    / "co2_ohess.out"
)


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


def _xtb_artifact() -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id="xtb-result-co2-opt",
        kind="xtb_output",
        sha256=hashlib.sha256(_XTB_OPT.read_bytes()).hexdigest(),
        size_bytes=_XTB_OPT.stat().st_size,
        path=str(_XTB_OPT.resolve()),
        cli_value=str(_XTB_OPT.resolve()),
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


@pytest.mark.capability("selector:xtb:opt:reached_positions")
@pytest.mark.capability("tool:bind_reached_geometry")
def test_a_validated_xtb_optimum_is_carried_from_its_exact_sidecar(tmp_path):
    """The later-cycle consumer must reach xTB's actual ``xtbopt.xyz``.

    The n-hexane goal proved that an xTB optimisation's sidecar is present
    and that an in-cycle producer edge can bind it.  At the wave barrier the
    next Agent session sees the archived ``xtb_output`` instead, so this
    public recovery route is the discriminating consumer: it must select the
    declared ``as_reached`` result role and write the exact final frame,
    rather than silently falling back to the input geometry or requiring the
    Agent to copy coordinates.
    """

    from chemsmart.analysis.result_readers import reader_for

    reader = reader_for("xtb")
    output = reader.open_output(_XTB_OPT)
    assert output.jobtype == "opt"
    assert "reached_positions" in reader.selectors_in_state_for_output(
        output, "as_reached"
    )

    artifact, receipt = build_reached_geometry(
        approved_workspace=tmp_path,
        reached_artifact_id="co2-optimum",
        result_artifact=_xtb_artifact(),
        program="xtb",
    )

    expected, _unit = reader.read(output, "reached_positions")
    observed = np.asarray(
        [
            line.split()[1:4]
            for line in Path(artifact.path).read_text().splitlines()[2:]
        ],
        dtype=float,
    )
    assert np.allclose(observed, np.asarray(expected), atol=1e-12)
    assert receipt.source_result_sha256 == _xtb_artifact().sha256
    assert dict(receipt.source_result_level) == {"method": "GFN2-xTB"}
    assert receipt.source_geometry_filename == "xtbopt.xyz"
    assert (
        receipt.source_geometry_sha256
        == hashlib.sha256(
            (_XTB_OPT.parent / "xtbopt.xyz").read_bytes()
        ).hexdigest()
    )
    assert receipt.normal_termination is True


@pytest.mark.capability("tool:bind_reached_geometry")
@pytest.mark.capability("selector:xtb:opt:reached_positions")
@pytest.mark.capability("program_jobtype:pyscf:cpu:sp")
def test_xtb_reached_geometry_can_anchor_a_new_pyscf_calculation(tmp_path):
    """One program-neutral handoff, then a distinct new electronic level.

    The test drives the public session tools from an archived xTB result into
    a PySCF calculation plan.  It is deliberately a *new* plan after an
    explicit charge/multiplicity binding: the geometry travels, but neither
    the xTB stationary-point verdict nor its electronic surface does.
    """

    from chemsmart.agent.runtime.event_store import RuntimeEventStore
    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    workspace = tmp_path / "workspace"
    workspace.mkdir()
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="cross-program"
        ),
        artifacts={"xtb-opt": _xtb_artifact()},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=workspace,
    )
    carried = host.dispatch(
        turn_id="t1",
        tool_name="bind_reached_geometry",
        arguments={
            "artifact_id": "xtb-opt",
            "reached_artifact_id": "co2-xtb-reached",
            "program": "xtb",
        },
    )["result"]
    assert carried["artifact"]["kind"] == "geometry_xyz"
    assert dict(carried["reached_geometry"]["source_result_level"]) == {
        "method": "GFN2-xTB"
    }
    host.dispatch(
        turn_id="t2",
        tool_name="bind_scientific_identity",
        arguments={
            "input_artifact_id": "co2-xtb-reached",
            "charge": 0,
            "multiplicity": 1,
        },
    )
    planned = host.dispatch(
        turn_id="t3",
        tool_name="plan_scientific_workflow",
        arguments={
            "plan_id": "xtb-geometry-pyscf-sp",
            "workflow_id": "xtb-geometry-pyscf-sp",
            "task_spec_id": "a" * 64,
            "required_output_ids": [],
            "analysis_nodes": [],
            "calculation_nodes": [
                {
                    "node_id": "pyscf-sp",
                    "program": "pyscf",
                    "jobtype": "sp",
                    "project_role": "project.pyscf",
                    "dependencies": [],
                    "inputs": [
                        {
                            "binding_id": "geometry.initial",
                            "artifact_id": "co2-xtb-reached",
                            "artifact_class": "xyz",
                            "producer_node_id": "",
                            "producer_output_id": "",
                        }
                    ],
                    "expected_outputs": [
                        {
                            "output_id": "structured-result",
                            "artifact_class": "pyscf_hdf5",
                        }
                    ],
                    "unresolved_fields": [],
                    "produces_observables": [],
                    "support_state": "planned",
                    "blocked_reason": "",
                }
            ],
        },
    )["result"]
    assert planned["workflow_frontier"]["actionable_node_ids"] == ["pyscf-sp"]
    draft = host._latest_program_workflows["xtb-geometry-pyscf-sp"].draft
    node = draft.nodes[0]
    assert node.program == "pyscf"
    assert node.inputs[0].artifact_id == carried["artifact"]["artifact_id"]
    assert any(
        binding.geometry_artifact_sha256 == carried["artifact"]["sha256"]
        and binding.charge == 0
        and binding.multiplicity == 1
        for binding in host.scientific_identities.values()
    )


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
