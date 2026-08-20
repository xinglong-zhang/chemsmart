"""The universal guard: no bare 64-hex digest ever reaches a human panel."""

from __future__ import annotations

import re
from types import SimpleNamespace

import pytest

pytest.importorskip("textual")

from rich.console import Console  # noqa: E402

from chemsmart.agent.tui.review import render_review_blocks  # noqa: E402

_HEX64 = re.compile(r"[0-9a-f]{64}")


def assert_no_bare_hex(text: str) -> None:
    match = _HEX64.search(text)
    assert match is None, f"digest leaked into a human panel: {match.group(0)}"


def test_a_digest_stuffed_review_renders_hash_free():
    node = SimpleNamespace(
        node_id="sp-arrangement",
        program="xtb",
        engine="cpu",
        stage="sp",
        project_artifact_sha256="c" * 64,
        molecular_identity={
            "approved_names": (),
            "formula": "H5NO",
            "atom_order": ("O", "H", "H", "N", "H", "H", "H"),
            "charge": 0,
            "multiplicity": 1,
            "identity_evidence_status": "composed-from-approved-parents",
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
                "achieved_contact_distance_angstrom": 1.94,
                "min_interfragment_distance_angstrom": 1.94,
                "atom_count": 7,
                "formula": "H5NO",
                "atom_order_note": "fragment A atoms first",
            },
        },
        environment_summary={"status": "available"},
        project_settings_text='{"a": 1}',
        real_execution_argv=(
            "chemsmart",
            "run",
            "xtb",
            "sp",
            "<geometry:sha256=" + "5" * 64 + ">",
        ),
    )
    chain = SimpleNamespace(
        analysis_nodes=(
            SimpleNamespace(
                node_id="extract",
                analysis_kind="result_extraction",
                inputs=(
                    SimpleNamespace(
                        producer_node_id="sp-arrangement",
                        producer_output_id="result-1",
                    ),
                ),
                outputs=(SimpleNamespace(output_id="e", unit="hartree"),),
                temperature_k=None,
                pressure_atm=None,
                support_state="planned",
                blocked_reason="",
            ),
        )
    )
    review = SimpleNamespace(
        execution_resources=SimpleNamespace(
            cores=8, memory_gb=16.0, gpu_count=0, node_timeout_seconds=900
        ),
        execution_envelope={"max_engine_calls": 3},
        node_reviews=(node,),
        non_executable_node_ids=(),
        scientific_plan=SimpleNamespace(
            nodes=(node,), edges=(), plan_sha256="b" * 64
        ),
        scientific_toolchain_plan=chain,
        review_sha256="a" * 64,
    )

    console = Console(record=True, width=220)
    for block in render_review_blocks(review):
        console.print(block)

    assert_no_bare_hex(console.export_text())
