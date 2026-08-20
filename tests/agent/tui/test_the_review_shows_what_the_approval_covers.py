"""The one /approve covers the analysis chain and the composed lineage.

Charter step 3 requires the terminal to display the complete plan. These
render-capture tests pin that the analysis chain travels visibly with the
review, that a composed arrangement shows its parents and contact by NAME
(digests are provenance and never reach a human panel), and that a chainless
review says so explicitly instead of staying silent.
"""

from __future__ import annotations

from types import SimpleNamespace

import re

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.review import render_review_blocks  # noqa: E402


def _flatten(blocks) -> str:
    console = Console(record=True, width=220)
    for block in blocks:
        console.print(block)
    return " ".join(console.export_text().replace("│", " ").split())


def _node_review(**overrides):
    values = dict(
        node_id="sp-arrangement",
        program="xtb",
        engine="cpu",
        stage="sp",
        molecular_identity={
            "approved_names": (),
            "formula": "H5NO",
            "atom_order": ("O", "H", "H", "N", "H", "H", "H"),
            "charge": 0,
            "multiplicity": 1,
        },
        environment_summary={
            "status": "available",
            "target_kind": "executable",
            "observed_version": "6.7.1",
            "observation_method": "host probe",
        },
        project_settings_text="{}",
        real_execution_argv=("chemsmart", "run", "xtb", "sp"),
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def _review(**overrides):
    node = _node_review()
    values = dict(
        execution_resources=SimpleNamespace(
            cores=8, memory_gb=16.0, gpu_count=0, node_timeout_seconds=900
        ),
        execution_envelope={"max_engine_calls": 3},
        node_reviews=(node,),
        non_executable_node_ids=(),
        scientific_plan=SimpleNamespace(
            nodes=(node,), edges=(), plan_sha256="b" * 64
        ),
        scientific_toolchain_plan=None,
        review_sha256="a" * 64,
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def test_a_chainless_review_says_so_instead_of_staying_silent():
    text = _flatten(render_review_blocks(_review()))

    assert "No typed analysis chain is planned with this workflow" in text
    assert re.search(r"[0-9a-f]{64}", text) is None
    assert "sha256" not in text
    assert "The full review record is kept in the run evidence" in text


def test_the_analysis_chain_is_displayed_with_the_workflow():
    chain = SimpleNamespace(
        analysis_nodes=(
            SimpleNamespace(
                node_id="thermo",
                analysis_kind="thermochemistry",
                inputs=(
                    SimpleNamespace(
                        producer_node_id="freq",
                        producer_output_id="result-freq",
                    ),
                ),
                outputs=(
                    SimpleNamespace(output_id="gibbs", unit="hartree"),
                ),
                temperature_k=298.15,
                pressure_atm=1.0,
                support_state="planned",
                blocked_reason="",
            ),
            SimpleNamespace(
                node_id="verdict",
                analysis_kind="scientific_validation",
                inputs=(
                    SimpleNamespace(
                        producer_node_id="thermo",
                        producer_output_id="gibbs",
                    ),
                ),
                outputs=(
                    SimpleNamespace(output_id="ok", unit="1"),
                ),
                temperature_k=None,
                pressure_atm=None,
                support_state="blocked_unsupported",
                blocked_reason="rule family not in this release",
            ),
        )
    )
    text = _flatten(
        render_review_blocks(_review(scientific_toolchain_plan=chain))
    )

    assert "Typed analysis chain" in text
    assert "runs provider-free after every approved calculation node" in text
    assert "thermo" in text and "freq.result-freq" in text
    assert "gibbs (hartree)" in text
    assert "298.15 K" in text
    assert (
        "not executable in this release: rule family not in this release"
        in text
    )
    assert "The displayed analysis chain executes provider-free" in text


def test_a_composed_arrangement_shows_its_parents_and_contact():
    identity = {
        "identity_evidence_status": "composed-from-approved-parents",
        "approved_names": (),
        "formula": "H5NO",
        "atom_order": ("O", "H", "H", "N", "H", "H", "H"),
        "charge": 0,
        "multiplicity": 1,
        "composition": {
            "fragment_a_artifact_id": "water-monomer",
            "fragment_a_sha256": "1" * 64,
            "fragment_a_identity_sha256": "2" * 64,
            "fragment_b_artifact_id": "ammonia-monomer",
            "fragment_b_sha256": "3" * 64,
            "fragment_b_identity_sha256": "4" * 64,
            "placement": {
                "mode": "contact",
                "contact": {
                    "fragment_a_atom": 2,
                    "fragment_b_atom": 1,
                    "distance_angstrom": 1.94,
                },
            },
            "achieved_contact_distance_angstrom": 1.9400000001,
            "min_interfragment_distance_angstrom": 1.94,
            "atom_count": 7,
            "formula": "H5NO",
            "atom_order_note": "fragment A atoms first, then fragment B",
        },
    }
    node = _node_review(molecular_identity=identity)
    review = _review(
        node_reviews=(node,),
        scientific_plan=SimpleNamespace(
            nodes=(node,), edges=(), plan_sha256="b" * 64
        ),
    )
    text = _flatten(render_review_blocks(review))

    assert "composed arrangement lineage" in text
    assert "covered by this approval" in text
    assert "water-monomer" in text and "ammonia-monomer" in text
    assert "atom 2 of A to atom 1 of B" in text
    assert "requested distance: 1.94" in text
    assert "built from two approved parent structures" in text
    assert re.search(r"[0-9a-f]{64}", text) is None, "parent digests leaked"
