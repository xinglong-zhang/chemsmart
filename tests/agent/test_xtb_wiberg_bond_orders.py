"""Tests for xTB Wiberg bond orders, dispersion energy, and evidence binding.

Pins the complete information path:
    real xTB output artifact
        ↓
    program result reader
        ↓
    public selector (wiberg_bond_orders, dispersion_energy)
        ↓
    typed extraction receipt
        ↓
    Agent-visible tool host & metadata inspection
"""

import shutil
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.postprocessing import (
    extract_trusted_result_quantities,
    typed_result_artifact_kind,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1
from chemsmart.analysis.result_quantities import (
    DIMENSIONLESS,
    QuantityExtractionError,
    QuantitySelectorV1,
)
from chemsmart.analysis.result_readers import (
    atom_resolved_selector_metadata,
    reader_for,
)

_WATER_OHESS = "tests/data/XTBTests/outputs/water_ohess/water_ohess.out"
_ACETALDEHYDE_HESS = (
    "tests/data/XTBTests/outputs/acetaldehyde_hess/acetaldehyde_hess.out"
)
_CO2_OHESS = "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
_BENZYNE_SP = "tests/data/XTBTests/outputs/p_benzyne_sp_alpb_toluene/p_benzyne_sp_alpb_toluene.out"
_BENZYNE_OPT = "tests/data/XTBTests/outputs/p_benzyne_opt_alpb_toluene/p_benzyne_opt_alpb_toluene.out"
_HE_HESS = "tests/data/XTBTests/outputs/he_hess/he_hess.out"


def _artifact(
    path: str | Path, program: str, artifact_id: str = "art-1"
) -> TrustedArtifactRefV1:
    resolved = Path(path).resolve()
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=typed_result_artifact_kind(program),
        sha256=file_sha256(resolved),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


@pytest.mark.capability("selector:xtb:hess:wiberg_bond_orders")
@pytest.mark.capability("selector:xtb:hess:dispersion_energy")
def test_water_wiberg_bond_orders_and_dispersion_extracted():
    """Water calculation yields sparse O-H WBO records and dispersion energy."""
    artifact = _artifact(_WATER_OHESS, "xtb", "water-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
            QuantitySelectorV1(
                quantity_id="disp",
                selector="dispersion_energy",
            ),
            QuantitySelectorV1(
                quantity_id="elements",
                selector="symbols",
            ),
        ),
    )
    assert receipt.status == "extracted"
    delivered = {item.quantity_id: item for item in receipt.quantities}

    # Symbols check: O, H, H
    assert delivered["elements"].value == ("O", "H", "H")

    # WBO check: sparse records of [atom_i, atom_j, bond_order]
    # In water: 0=O, 1=H, 2=H
    # Only pairs above xTB threshold (>0.10) are reported: O(0)-H(1) and O(0)-H(2)
    wbo_val = delivered["wbo"].value
    assert len(wbo_val) == 2
    assert delivered["wbo"].data_kind == "matrix"
    assert delivered["wbo"].unit == "1"
    assert delivered["wbo"].dimension == DIMENSIONLESS

    # O(0) - H(1)
    assert wbo_val[0][0] == 0
    assert wbo_val[0][1] == 1
    assert wbo_val[0][2] == pytest.approx(0.920214, abs=1e-5)

    # O(0) - H(2)
    assert wbo_val[1][0] == 0
    assert wbo_val[1][1] == 2
    assert wbo_val[1][2] == pytest.approx(0.920214, abs=1e-5)

    # Provenance binds to artifact
    assert (
        delivered["wbo"].evidence_ref
        == f"artifact:{artifact.artifact_id}#{artifact.sha256}"
    )
    assert receipt.level == {"method": "GFN2-xTB"}
    assert receipt.native_evidence == (
        (
            "wiberg_bond_orders",
            "wbo",
            file_sha256(Path(_WATER_OHESS).with_name("wbo")),
        ),
    )

    # Dispersion energy in Hartree
    assert delivered["disp"].unit == "hartree"
    assert delivered["disp"].source_unit == "Eh"
    assert delivered["disp"].value == pytest.approx(-0.000141082856, abs=1e-8)
    assert "dispersion_energy" not in {
        selector for selector, _filename, _sha256 in receipt.native_evidence
    }


@pytest.mark.capability("selector:xtb:hess:wiberg_bond_orders")
def test_acetaldehyde_non_symmetric_wbo_indices_and_order():
    """Non-symmetric acetaldehyde tests that bond orders map to exact atom pairs."""
    artifact = _artifact(_ACETALDEHYDE_HESS, "xtb", "acetaldehyde-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
            QuantitySelectorV1(
                quantity_id="disp",
                selector="dispersion_energy",
            ),
            QuantitySelectorV1(
                quantity_id="elements",
                selector="symbols",
            ),
        ),
    )
    delivered = {item.quantity_id: item.value for item in receipt.quantities}
    symbols = delivered["elements"]
    assert len(symbols) == 7
    wbo = delivered["wbo"]
    assert len(wbo) == 6

    # From wbo file:
    # 1 3 1.93034 -> atoms 0 and 2: carbonyl C=O
    assert wbo[0][0] == 0 and wbo[0][1] == 2
    assert wbo[0][2] == pytest.approx(1.930344, abs=1e-5)

    # 2 3 1.02744 -> atoms 1 and 2: C-C single bond
    assert wbo[1][0] == 1 and wbo[1][1] == 2
    assert wbo[1][2] == pytest.approx(1.027442, abs=1e-5)

    # 2 4 0.98250 -> atoms 1 and 3: methyl C-H
    assert wbo[2][0] == 1 and wbo[2][1] == 3
    assert wbo[2][2] == pytest.approx(0.982505, abs=1e-5)

    # 2 5 0.95552 -> atoms 1 and 4: methyl C-H
    assert wbo[3][0] == 1 and wbo[3][1] == 4
    assert wbo[3][2] == pytest.approx(0.955519, abs=1e-5)

    # 2 6 0.95499 -> atoms 1 and 5: methyl C-H
    assert wbo[4][0] == 1 and wbo[4][1] == 5
    assert wbo[4][2] == pytest.approx(0.954987, abs=1e-5)

    # 3 7 0.93161 -> atoms 2 and 6: formyl C-H
    assert wbo[5][0] == 2 and wbo[5][1] == 6
    assert wbo[5][2] == pytest.approx(0.931605, abs=1e-5)

    # Dispersion energy
    disp = delivered["disp"]
    assert disp == pytest.approx(-0.002041721444, abs=1e-8)


@pytest.mark.capability("selector:xtb:sp:wiberg_bond_orders")
@pytest.mark.capability("selector:xtb:sp:dispersion_energy")
def test_xtb_sp_jobtype_exposes_wbo_and_dispersion():
    """Single-point calculation legitimately exposes WBO and dispersion energy."""
    artifact = _artifact(_BENZYNE_SP, "xtb", "benzyne-sp")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
            QuantitySelectorV1(
                quantity_id="disp",
                selector="dispersion_energy",
            ),
        ),
    )
    delivered = {item.quantity_id: item for item in receipt.quantities}
    assert delivered["disp"].value == pytest.approx(-0.006740948916, abs=1e-8)
    assert len(delivered["wbo"].value) > 0


@pytest.mark.capability("selector:xtb:opt:wiberg_bond_orders")
@pytest.mark.capability("selector:xtb:opt:dispersion_energy")
def test_xtb_opt_jobtype_exposes_wbo_and_dispersion():
    """Optimization jobtype exposes reached WBO and dispersion energy."""
    artifact = _artifact(_BENZYNE_OPT, "xtb", "benzyne-opt")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
            QuantitySelectorV1(
                quantity_id="disp",
                selector="dispersion_energy",
            ),
        ),
    )
    delivered = {item.quantity_id: item for item in receipt.quantities}
    assert len(delivered["wbo"].value) > 0
    assert delivered["disp"].value == pytest.approx(-0.006740949744, abs=1e-8)


@pytest.mark.capability("selector:xtb:hess:wiberg_bond_orders")
def test_monoatomic_he_hess_wbo_is_absent_in_receipt():
    """Monoatomic system with empty wbo sidecar yields explicit absence."""
    artifact = _artifact(_HE_HESS, "xtb", "he-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
        ),
    )
    assert receipt.status == "partial"
    assert len(receipt.absent) == 1
    qid, sel, reason = receipt.absent[0]
    assert qid == "wbo"
    assert sel == "wiberg_bond_orders"
    assert "empty WBO sidecar" in reason


def test_missing_wbo_sidecar_yields_explicit_absence(tmp_path):
    """Result without a wbo sidecar yields explicit absence, not inferred zeros."""
    out_copy = tmp_path / "run.out"
    shutil.copy(Path(_WATER_OHESS), out_copy)
    # Note: do not copy wbo to tmp_path, so sidecar is truly missing
    artifact = _artifact(out_copy, "xtb", "missing-wbo")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
        ),
    )
    assert receipt.status == "partial"
    assert len(receipt.absent) == 1
    qid, sel, reason = receipt.absent[0]
    assert qid == "wbo"
    assert sel == "wiberg_bond_orders"
    assert "wrote no Wiberg bond order (wbo) sidecar" in reason


def test_wiberg_bond_orders_metadata_is_registered():
    """Metadata inspection reports Wiberg scheme, atom order, and sparsity."""
    meta = atom_resolved_selector_metadata("wiberg_bond_orders")
    assert meta == {
        "semantic_quantity": "bond_order",
        "population_scheme": "Wiberg",
        "atom_order": "zero-based molecular atom order",
        "data_shape": "rows of [atom_i, atom_j, wiberg_bond_order]",
        "sparsity": (
            "the native xTB sidecar is thresholded; omitted pairs have no "
            "reported Wiberg value and are not zero"
        ),
    }


def test_wbo_reaches_agent_inspect_run_and_durable_receipt(tmp_path):
    """Tool host inspect_run lists metadata and extraction mints valid receipt."""
    artifact = _artifact(_WATER_OHESS, "xtb", "water-hess")
    event_path = tmp_path / "events.jsonl"
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(event_path, session_id="test-wbo"),
        artifacts={artifact.artifact_id: artifact},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    inspected = host._inspect_run(
        "t1", {"program": "xtb", "artifact_id": artifact.artifact_id}
    )
    assert "wiberg_bond_orders" in inspected["requestable_selectors"]
    assert "dispersion_energy" in inspected["requestable_selectors"]
    assert inspected["atom_resolved_metadata"]["wiberg_bond_orders"] == {
        "semantic_quantity": "bond_order",
        "population_scheme": "Wiberg",
        "atom_order": "zero-based molecular atom order",
        "data_shape": "rows of [atom_i, atom_j, wiberg_bond_order]",
        "sparsity": (
            "the native xTB sidecar is thresholded; omitted pairs have no "
            "reported Wiberg value and are not zero"
        ),
    }

    receipt = host._extract_result_quantities(
        "t2",
        {
            "program": "xtb",
            "artifact_id": artifact.artifact_id,
            "selectors": [
                {
                    "quantity_id": "water_wbo",
                    "selector": "wiberg_bond_orders",
                },
                {
                    "quantity_id": "water_disp",
                    "selector": "dispersion_energy",
                },
            ],
        },
    )
    assert receipt.receipt_sha256
    rehydrated = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(event_path, session_id="test-wbo"),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    restored = rehydrated.quantity_extractions[receipt.receipt_sha256]
    assert restored.selector_bindings == receipt.selector_bindings
    assert restored.structural_states == (
        ("water_wbo", "as_reached"),
        ("water_disp", "as_reached"),
    )
    assert (
        restored.native_evidence
        == receipt.native_evidence
        == (
            (
                "wiberg_bond_orders",
                "wbo",
                file_sha256(Path(_WATER_OHESS).with_name("wbo")),
            ),
        )
    )

    public_host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "public-events.jsonl", session_id="public-wbo"
        ),
        artifacts={artifact.artifact_id: artifact},
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "public-workspace",
    )
    reply = public_host.dispatch(
        turn_id="t3",
        tool_name="extract_result_quantities",
        arguments={
            "program": "xtb",
            "artifact_id": artifact.artifact_id,
            "selectors": [
                {
                    "quantity_id": "water_wbo",
                    "selector": "wiberg_bond_orders",
                }
            ],
        },
    )
    assert reply["result"]["native_evidence"] == [
        [
            "wiberg_bond_orders",
            "wbo",
            file_sha256(Path(_WATER_OHESS).with_name("wbo")),
        ]
    ]


def test_xtb_scc_charges_bind_their_native_sidecar():
    """A second xTB sidecar uses the same receipt-level authority."""

    artifact = _artifact(_WATER_OHESS, "xtb", "water-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="scc_charges",
                selector="xtb_scc_atomic_charges",
            ),
        ),
    )
    assert receipt.status == "extracted"
    assert receipt.native_evidence == (
        (
            "xtb_scc_atomic_charges",
            "charges",
            file_sha256(Path(_WATER_OHESS).with_name("charges")),
        ),
    )


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("1 2\n", "malformed"),
        ("1 2 0.9\n2 1 0.9\n", "repeats unordered atom pair"),
        ("1 1 0.9\n", "self-pair"),
        ("1 2 nan\n", "non-finite or negative"),
        ("1 2 -0.1\n", "non-finite or negative"),
        ("1 4 0.9\n", "out of bounds"),
    ],
)
def test_malformed_or_ambiguous_wbo_sidecar_fails_closed(
    tmp_path, contents, message
):
    """A plausible partial pair table is not accepted as WBO evidence."""

    copied = tmp_path / "water"
    shutil.copytree(Path(_WATER_OHESS).parent, copied)
    (copied / "wbo").write_text(contents)
    artifact = _artifact(copied / Path(_WATER_OHESS).name, "xtb", "water")
    with pytest.raises(QuantityExtractionError, match=message):
        extract_trusted_result_quantities(
            artifact=artifact,
            program="xtb",
            selectors=(
                QuantitySelectorV1(
                    quantity_id="wbo", selector="wiberg_bond_orders"
                ),
            ),
        )


def test_changed_wbo_sidecar_during_extraction_is_refused(
    tmp_path, monkeypatch
):
    """The main log digest cannot stand in for mutable native WBO bytes."""

    copied = tmp_path / "water"
    shutil.copytree(Path(_WATER_OHESS).parent, copied)
    artifact = _artifact(copied / Path(_WATER_OHESS).name, "xtb", "water")
    reader = reader_for("xtb")
    assert reader is not None
    original = reader.accessors["wiberg_bond_orders"]

    def mutate_then_read(output):
        Path(output.wbo_file.filename).write_text("1 2 0.5\n1 3 0.5\n")
        return original(output)

    monkeypatch.setitem(
        reader.accessors, "wiberg_bond_orders", mutate_then_read
    )
    with pytest.raises(
        QuantityExtractionError, match="sidecar changed during extraction"
    ):
        extract_trusted_result_quantities(
            artifact=artifact,
            program="xtb",
            selectors=(
                QuantitySelectorV1(
                    quantity_id="wbo", selector="wiberg_bond_orders"
                ),
            ),
        )
