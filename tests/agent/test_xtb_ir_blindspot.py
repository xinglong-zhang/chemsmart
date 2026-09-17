"""Tests and witnesses for xTB IR absorption intensity (ir_intensities).

Covers:
1. Boundary A: Scientific vocabulary dimension and canonical unit km/mol.
2. Boundary B: Canonical ResultReader and Agent extraction for ir_intensities.
3. Level 2 Workflow: Index-aligned frequency and intensity table, dynamic maximum mode identification.
4. Source Precedence: vibspectrum_file priority over main_out (overcoming ****** fixed-width truncation).
5. Provenance: Sealing vibspectrum sidecar into native evidence paths.
6. Adversarial checks: Missing evidence, mismatched sequence lengths, non-finite values, and unsupported job types.
7. Existing selector regressions: Zero regressions across existing xTB selectors.
"""

import math
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.postprocessing import (
    extract_trusted_result_quantities,
    typed_result_artifact_kind,
)
from chemsmart.analysis.quantity_expressions import (
    canonical_unit_for_dimension,
    unit_dimension,
)
from chemsmart.analysis.result_quantities import (
    IR_INTENSITY,
    QuantityContractError,
    QuantitySelectorV1,
)
from chemsmart.analysis.result_readers import (
    RESULT_READERS,
    MissingQuantityError,
    _xtb_ir_intensities,
    _xtb_native_evidence_paths,
)
from chemsmart.io.xtb.output import XTBOutput

_ACETALDEHYDE_HESS = (
    "tests/data/XTBTests/outputs/acetaldehyde_hess/acetaldehyde_hess.out"
)
_CO2_HESS = "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"


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


def test_xtb_frequency_result_has_parsed_ir_intensities_in_output():
    """Verify that the xTB calculation folder genuinely contains parsed IR intensities."""
    artifact_path = Path(_ACETALDEHYDE_HESS)
    output = XTBOutput(artifact_path.parent)
    assert output.vibspectrum_file is not None
    intensities = output.vibspectrum_file.ir_intensities
    assert intensities is not None
    assert len(intensities) == 15
    assert all(isinstance(val, float) and val >= 0.0 for val in intensities)


def test_scientific_vocabulary_boundary_ir_intensity_unit():
    """Boundary A: prove that the typed scientific vocabulary accepts km/mol."""
    dim = unit_dimension("km/mol")
    assert len(dim) == 10
    assert dim == IR_INTENSITY
    assert canonical_unit_for_dimension(IR_INTENSITY) == "km/mol"

    # Compound and alias spellings
    assert unit_dimension("km mol^-1") == IR_INTENSITY
    assert unit_dimension("km·mol⁻¹") == IR_INTENSITY


@pytest.mark.capability("selector:xtb:hess:ir_intensities")
def test_canonical_agent_extraction_can_request_ir_intensities():
    """Boundary B: prove that canonical Agent extraction exposes ir_intensities."""
    artifact = _artifact(_ACETALDEHYDE_HESS, "xtb", "acetaldehyde-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="freqs",
                selector="vibrational_frequencies",
            ),
            QuantitySelectorV1(
                quantity_id="ir",
                selector="ir_intensities",
            ),
        ),
    )
    assert receipt.status == "extracted"
    delivered = {item.quantity_id: item for item in receipt.quantities}
    assert "ir" in delivered
    assert delivered["ir"].unit == "km/mol"
    assert len(delivered["ir"].value) == len(delivered["freqs"].value)


def test_co2_precedence_preserves_large_intensity_from_vibspectrum():
    """Section 5: verify vibspectrum_file priority overcomes main_out's ****** fixed-width truncation."""
    output = XTBOutput(Path(_CO2_HESS).parent)
    assert output.vibspectrum_file is not None
    # Mode 9 is asymmetric CO2 stretch: main_out printed '9:******', but vibspectrum printed 1046.66491
    ir_values = output.ir_intensities
    assert ir_values is not None
    assert len(ir_values) == 4
    # Mode 9 (4th non-zero mode) must be exactly preserved from vibspectrum
    assert math.isclose(ir_values[-1], 1046.66491, abs_tol=1e-3)


def test_truncated_main_output_is_absent_not_a_shifted_mode_table():
    """Asterisks in xTB's indexed table cannot become a fabricated zero."""
    output = XTBOutput(Path(_CO2_HESS).parent)
    assert output.main_out is not None
    assert output.main_out.ir_intensities is None


def test_main_output_preserves_indexed_ir_modes_when_complete():
    """The fallback parser accepts complete packed ``index:value`` rows."""
    output = XTBOutput(Path(_ACETALDEHYDE_HESS).parent)
    assert output.main_out is not None
    intensities = output.main_out.ir_intensities
    assert intensities is not None
    assert len(intensities) == 15
    assert math.isclose(intensities[10], 300.94, abs_tol=1e-6)


def test_level_2_workflow_frequency_and_ir_intensity_table_and_max_mode():
    """Level 2: Extract frequency and intensity together, build index-aligned table, find max mode."""
    artifact = _artifact(_ACETALDEHYDE_HESS, "xtb", "acetaldehyde-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="freqs",
                selector="vibrational_frequencies",
            ),
            QuantitySelectorV1(
                quantity_id="ir",
                selector="ir_intensities",
            ),
        ),
    )
    delivered = {item.quantity_id: item for item in receipt.quantities}
    frequencies = delivered["freqs"].value
    intensities = delivered["ir"].value
    assert len(frequencies) == len(intensities) == 15

    # Construct index-aligned table without hard-coding
    table = [
        (idx + 1, freq, intensity)
        for idx, (freq, intensity) in enumerate(zip(frequencies, intensities))
    ]
    # Identify mode with largest IR absorption intensity dynamically from extracted evidence
    max_mode_idx, max_freq, max_intensity = max(table, key=lambda row: row[2])
    # For acetaldehyde: mode 17 (11th vibrational mode) is C=O stretch: ~300.94 km/mol, ~1758 cm^-1
    assert max_intensity > 250.0
    assert 1700.0 <= max_freq <= 1800.0
    assert max_mode_idx == 11  # 11th vibrational mode out of 15


def test_provenance_sidecar_binding_for_ir_intensities():
    """The extraction receipt seals the paired sidecar bytes it consumed."""
    artifact = _artifact(_ACETALDEHYDE_HESS, "xtb", "acetaldehyde-hess")
    output = XTBOutput(Path(_ACETALDEHYDE_HESS).parent)
    paths = _xtb_native_evidence_paths(output, "ir_intensities")
    assert len(paths) == 1
    assert paths[0].name == "vibspectrum"
    assert paths[0].is_file()
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="freqs",
                selector="vibrational_frequencies",
            ),
            QuantitySelectorV1(
                quantity_id="ir",
                selector="ir_intensities",
            ),
        ),
    )
    expected_sha = file_sha256(paths[0])
    assert ("vibrational_frequencies", "vibspectrum", expected_sha) in (
        receipt.native_evidence
    )
    assert (
        "ir_intensities",
        "vibspectrum",
        expected_sha,
    ) in receipt.native_evidence


def test_unsupported_jobtypes_do_not_advertise_ir_intensities():
    """Section 6: Verify sp and opt do not advertise or admit ir_intensities."""
    reader = RESULT_READERS["xtb"]
    jobtypes = dict(reader.jobtype_selectors)
    assert "ir_intensities" in jobtypes["hess"]
    assert "ir_intensities" not in jobtypes["opt"]
    assert "ir_intensities" not in jobtypes["sp"]


def test_adversarial_missing_frequencies_fails_closed():
    """Section 10: Missing vibrational frequencies must fail closed."""
    fake_output = MagicMock()
    fake_output.vibrational_frequencies = None
    fake_output.ir_intensities = [10.0, 20.0]
    with pytest.raises(
        MissingQuantityError,
        match="require corresponding vibrational frequencies",
    ):
        _xtb_ir_intensities(fake_output)


def test_adversarial_mismatched_frequency_intensity_length_fails_closed():
    """Section 10: Length mismatch between frequencies and intensities must fail closed."""
    fake_output = MagicMock()
    fake_output.vibrational_frequencies = [100.0, 200.0, 300.0]
    fake_output.ir_intensities = [10.0, 20.0]  # length 2 != length 3
    with pytest.raises(QuantityContractError, match="count .* does not match"):
        _xtb_ir_intensities(fake_output)


def test_adversarial_non_finite_intensity_fails_closed():
    """Section 10: NaN or Inf intensity values must fail closed."""
    fake_output = MagicMock()
    fake_output.vibrational_frequencies = [100.0, 200.0]
    fake_output.ir_intensities = [10.0, float("nan")]
    with pytest.raises(QuantityContractError, match="non-finite IR intensity"):
        _xtb_ir_intensities(fake_output)


def test_existing_xtb_selectors_remain_unchanged():
    """Verify that existing xTB selectors continue to extract normally."""
    artifact = _artifact(_ACETALDEHYDE_HESS, "xtb", "acetaldehyde-hess")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xtb",
        selectors=(
            QuantitySelectorV1(
                quantity_id="freqs",
                selector="vibrational_frequencies",
            ),
            QuantitySelectorV1(
                quantity_id="wbo",
                selector="wiberg_bond_orders",
            ),
            QuantitySelectorV1(
                quantity_id="disp",
                selector="dispersion_energy",
            ),
            QuantitySelectorV1(
                quantity_id="charges",
                selector="xtb_scc_atomic_charges",
            ),
        ),
    )
    assert receipt.status == "extracted"
    delivered = {item.quantity_id: item for item in receipt.quantities}
    assert "freqs" in delivered and len(delivered["freqs"].value) == 15
    assert "wbo" in delivered and len(delivered["wbo"].value) > 0
    assert "disp" in delivered and delivered["disp"].value < 0.0
    assert "charges" in delivered and len(delivered["charges"].value) == 7
