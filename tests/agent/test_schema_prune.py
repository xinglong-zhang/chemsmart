from __future__ import annotations

import copy
import json

import pytest

from chemsmart.agent.schema_prune import (
    prune_schema_for_request,
    schema_variant_id,
)

NODE_KEYS = {"name", "description", "options", "subcommands"}


def _node(name, **subcommands):
    return {
        "name": name,
        "description": f"help for {name}",
        "options": [],
        "subcommands": dict(subcommands),
    }


def _kind(name, with_qmmm=True):
    children = {"qmmm": _node("qmmm")} if with_qmmm else {}
    return _node(name, **children)


def make_schema():
    gaussian = _node(
        "gaussian",
        opt=_kind("opt"),
        sp=_kind("sp"),
        ts=_kind("ts"),
        td=_kind("td"),
        scan=_kind("scan"),
        wbi=_kind("wbi", with_qmmm=False),
    )
    orca = _node(
        "orca",
        opt=_kind("opt"),
        sp=_kind("sp"),
        ts=_kind("ts"),
        neb=_kind("neb", with_qmmm=False),
    )
    aux = {
        "mol": _node("mol"),
        "thermochemistry": _node("thermochemistry"),
        "nciplot": _node("nciplot"),
    }
    run = _node("run", gaussian=gaussian, orca=orca, **aux)
    sub = _node(
        "sub",
        gaussian=copy.deepcopy(gaussian),
        orca=copy.deepcopy(orca),
        **copy.deepcopy(aux),
    )
    return _node(
        "chemsmart",
        run=run,
        sub=sub,
        agent=_node("agent"),
        config=_node("config"),
        update=_node("update"),
    )


def _walk(node):
    yield node
    for child in (node.get("subcommands") or {}).values():
        yield from _walk(child)


def test_typical_opt_request_prunes_to_minimal_tree():
    schema = make_schema()
    pruned = prune_schema_for_request(
        schema, "geometry optimization of h2o.xyz, charge 0 multiplicity 1"
    )
    assert set(pruned["subcommands"]) == {"run"}
    programs = pruned["subcommands"]["run"]["subcommands"]
    # No engine named and no workspace: both engines stay.
    assert set(programs) == {"gaussian", "orca"}
    assert set(programs["gaussian"]["subcommands"]) == {"opt", "sp"}
    assert set(programs["orca"]["subcommands"]) == {"opt", "sp"}


def test_explicit_orca_drops_gaussian():
    pruned = prune_schema_for_request(
        make_schema(), "transition state search with orca for guess.xyz"
    )
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert set(programs) == {"orca"}
    assert set(programs["orca"]["subcommands"]) == {"opt", "sp", "ts"}


def test_explicit_gaussian_drops_orca():
    pruned = prune_schema_for_request(
        make_schema(), "gaussian tddft excited state calculation of dye.xyz"
    )
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert set(programs) == {"gaussian"}
    assert set(programs["gaussian"]["subcommands"]) == {"opt", "sp", "td"}


def test_workspace_program_used_when_request_is_silent():
    pruned = prune_schema_for_request(
        make_schema(),
        "optimize h2o.xyz charge 0 multiplicity 1",
        workspace_program="orca",
    )
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert set(programs) == {"orca"}


def test_explicit_program_beats_workspace_program():
    pruned = prune_schema_for_request(
        make_schema(),
        "optimize h2o.xyz with gaussian",
        workspace_program="orca",
    )
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert set(programs) == {"gaussian"}


def test_submit_keywords_select_sub_entry():
    pruned = prune_schema_for_request(
        make_schema(), "submit a wiberg bond index job to the cluster"
    )
    assert set(pruned["subcommands"]) == {"sub"}
    gaussian = pruned["subcommands"]["sub"]["subcommands"]["gaussian"]
    assert set(gaussian["subcommands"]) == {"opt", "sp", "wbi"}


def test_no_kind_cue_keeps_only_conservative_floor():
    pruned = prune_schema_for_request(
        make_schema(), "do something with my molecule please"
    )
    gaussian = pruned["subcommands"]["run"]["subcommands"]["gaussian"]
    assert set(gaussian["subcommands"]) == {"opt", "sp"}
    for kind_node in gaussian["subcommands"].values():
        assert "qmmm" not in kind_node["subcommands"]


def test_qmmm_request_keeps_qmmm_children():
    pruned = prune_schema_for_request(
        make_schema(), "oniom qm/mm optimization of protein.xyz"
    )
    opt = pruned["subcommands"]["run"]["subcommands"]["gaussian"][
        "subcommands"
    ]["opt"]
    assert "qmmm" in opt["subcommands"]


def test_qmmm_stripped_by_default():
    pruned = prune_schema_for_request(
        make_schema(), "optimize h2o.xyz with gaussian"
    )
    opt = pruned["subcommands"]["run"]["subcommands"]["gaussian"][
        "subcommands"
    ]["opt"]
    assert "qmmm" not in opt["subcommands"]


def test_aux_program_kept_only_when_named():
    silent = prune_schema_for_request(
        make_schema(), "optimize h2o.xyz with gaussian"
    )
    assert "thermochemistry" not in silent["subcommands"]["run"]["subcommands"]

    named = prune_schema_for_request(
        make_schema(), "thermochemistry analysis of my output files"
    )
    assert "thermochemistry" in named["subcommands"]["run"]["subcommands"]


def test_unknown_kind_for_program_falls_back_to_all_kinds():
    # ORCA has no td: with an orca-only request asking for tddft, the kind
    # filter would come up empty and must fall back to every orca jobkind.
    pruned = prune_schema_for_request(
        make_schema(), "orca excited state absorption spectrum"
    )
    orca = pruned["subcommands"]["run"]["subcommands"]["orca"]
    # td matched but orca lacks it; opt/sp floor still intersects, so the
    # selection is non-empty and stays minimal.
    assert set(orca["subcommands"]) == {"opt", "sp"}


def test_minimal_test_schema_passes_through_by_identity():
    minimal = {"subcommands": {}}
    assert prune_schema_for_request(minimal, "optimize h2o.xyz") is minimal
    stub = {"name": "chemsmart", "subcommands": {"foo": _node("foo")}}
    assert prune_schema_for_request(stub, "optimize h2o.xyz") is stub


def test_input_schema_is_not_mutated():
    schema = make_schema()
    snapshot = copy.deepcopy(schema)
    prune_schema_for_request(schema, "submit orca ts job to hpc queue")
    assert schema == snapshot


def test_node_shape_preserved_everywhere():
    pruned = prune_schema_for_request(
        make_schema(), "gaussian dihedral scan of butane.xyz"
    )
    for node in _walk(pruned):
        assert NODE_KEYS <= set(node)


def test_variant_id_formats():
    schema = make_schema()
    assert schema_variant_id(schema) == "full"
    assert schema_variant_id({"subcommands": {}}) == "empty"
    pruned = prune_schema_for_request(schema, "optimize h2o.xyz with gaussian")
    assert schema_variant_id(pruned) == "run/gaussian[opt,sp]"


@pytest.mark.slow
def test_real_schema_size_regression():
    from chemsmart.agent.cli_schema import build_chemsmart_cli_schema
    from chemsmart.agent.prompts.synthesis import (
        build_synthesis_system_prompt,
    )

    full = build_chemsmart_cli_schema()
    typical = prune_schema_for_request(
        full,
        "geometry optimization of h2o.xyz with gaussian, "
        "charge 0 multiplicity 1",
    )
    prompt = build_synthesis_system_prompt(typical)
    # Measured 2026-07: ~17k chars (~4k tokens). Guard against re-bloat.
    assert len(prompt) < 80_000

    # Worst case (no cues at all) must still fit 65k-token contexts.
    vague = prune_schema_for_request(full, "help me with my molecule")
    assert len(build_synthesis_system_prompt(vague)) < 150_000

    # The pruned tree must still be valid JSON-serializable schema.
    json.dumps(typical)


def test_filename_mol_xyz_does_not_pull_in_mol_program():
    pruned = prune_schema_for_request(
        make_schema(),
        "optimize mol.xyz with the bond between atoms 1 and 2 frozen, "
        "gaussian, charge 0 multiplicity 1",
    )
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert "mol" not in programs
    assert set(programs) == {"gaussian"}

    named = prune_schema_for_request(
        make_schema(), "visualize my molecule with mol"
    )
    assert "mol" in named["subcommands"]["run"]["subcommands"]


def _schema_with_qrc():
    schema = make_schema()
    for entry in ("run", "sub"):
        gaussian = schema["subcommands"][entry]["subcommands"]["gaussian"]
        gaussian["subcommands"]["qrc"] = _kind("qrc", with_qmmm=False)
    return schema


@pytest.mark.parametrize(
    "request_text",
    [
        "run a qrc job from ts_run.log",
        "quick reaction coordinate from the saddle point in ts.log",
        "quasi-reaction coordinate job from the saddle point in ts.log",
        "displace along the imaginary mode of barrier.log both directions",
    ],
)
def test_qrc_synonyms_keep_qrc_jobkind(request_text):
    # A request the pruner cannot map to qrc silently hides the subcommand,
    # so the model emits `ts` instead and the gate scores it WRONG.
    pruned = prune_schema_for_request(_schema_with_qrc(), request_text)
    gaussian = pruned["subcommands"]["run"]["subcommands"]["gaussian"]
    assert "qrc" in gaussian["subcommands"]


def _schema_with_xtb():
    schema = make_schema()
    for entry in ("run", "sub"):
        schema["subcommands"][entry]["subcommands"]["xtb"] = _node(
            "xtb",
            opt=_kind("opt", with_qmmm=False),
            sp=_kind("sp", with_qmmm=False),
            hess=_kind("hess", with_qmmm=False),
        )
    return schema


@pytest.mark.parametrize(
    "request_text",
    [
        "pre-optimize ligand.xyz cheaply before the DFT TS search",
        "preoptimization of the docked pose in complex.xyz",
        "semiempirical single point of aggregate.xyz",
        "tight-binding optimization of polymer.xyz",
        "water.xyz 사전최적화 해줘",
        "반경험 방법으로 에너지만 계산",
    ],
)
def test_preoptimization_phrasings_keep_xtb_program(request_text):
    # Missing a cue here prunes the xtb program out of the schema entirely,
    # leaving the model unable to emit the pre-optimization step at all.
    pruned = prune_schema_for_request(_schema_with_xtb(), request_text)
    programs = pruned["subcommands"]["run"]["subcommands"]
    assert "xtb" in programs


def test_frequency_wording_keeps_xtb_hess_leaf():
    pruned = prune_schema_for_request(
        _schema_with_xtb(), "xtb vibrational frequencies of water.xyz"
    )
    xtb = pruned["subcommands"]["run"]["subcommands"]["xtb"]
    assert "hess" in xtb["subcommands"]
