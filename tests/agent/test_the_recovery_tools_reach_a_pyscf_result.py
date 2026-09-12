"""The repair menu's routes are reachable for the program they were refused on.

``displace_along_vibrational_mode`` and ``characterise_stationary_point``
required an artifact of kind ``f"{program}_output"`` and the first also
branched ``if program == "orca"``; a PySCF result is ``pyscf_hdf5``, so
both refused every PySCF Hessian while the repair menu for
``failed_wrong_stationary_point`` named both. Measured refusal, before:
"a stationary point characterisation on pyscf requires a pyscf_output
artifact, not 'pyscf_hdf5'". Driven here on the archived planar-ammonia
saddle, one imaginary mode at -830 cm-1.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from chemsmart.agent._contracts import (
    ContractError,
    TrustedArtifactRefV1,
    file_sha256,
)
from chemsmart.agent.execution import (
    build_stationary_point_characterisation,
    displace_trusted_geometry_along_mode,
)
from chemsmart.analysis.result_readers import reader_for

SADDLE = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "PySCFTests"
    / "outputs"
    / "nh3_planar_hess"
    / "nh3_planar_hess_gas_phase.h5"
)


def _artifact(path: Path = SADDLE) -> TrustedArtifactRefV1:
    return TrustedArtifactRefV1(
        artifact_id="pyscf-result-saddle",
        kind="pyscf_hdf5",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


@pytest.mark.capability("tool:characterise_stationary_point")
def test_a_pyscf_saddle_is_characterised_by_its_own_printed_modes():
    receipt = build_stationary_point_characterisation(
        result_artifact=_artifact(), program="pyscf", order_claimed=1
    )
    assert receipt.observed_imaginary_modes == 1
    assert receipt.order_claimed == 1
    assert receipt.lowest_imaginary_cm_1 < -800.0
    with pytest.raises(ContractError, match="carries 1 imaginary mode"):
        build_stationary_point_characterisation(
            result_artifact=_artifact(), program="pyscf", order_claimed=0
        )


@pytest.mark.capability("tool:displace_along_vibrational_mode")
def test_a_pyscf_saddle_is_stepped_along_its_umbrella_mode(tmp_path):
    reader = reader_for("pyscf")
    output = reader.open_output(SADDLE)
    shares = np.asarray(
        reader.read(output, "vibrational_mode_atom_participation")[0]
    )[0]
    carriers = {
        index + 1 for index, share in enumerate(shares) if share > 0.05
    }

    displaced, receipt = displace_trusted_geometry_along_mode(
        approved_workspace=tmp_path,
        displaced_artifact_id="nh3-stepped",
        result_artifact=_artifact(),
        program="pyscf",
        mode_index=1,
        amplitude_angstrom=0.1,
    )
    assert displaced.kind == "geometry_xyz"
    assert receipt.mode_is_imaginary is True
    assert receipt.mode_frequency_cm_1 < -800.0
    assert (
        set(receipt.moved_atoms) >= carriers
    ), "the atoms the participation table names carry the step"
    lines = Path(displaced.path).read_text().splitlines()
    stepped = np.asarray(
        [[float(v) for v in line.split()[1:4]] for line in lines[2:6]]
    )
    before = np.asarray(output.positions)
    # A planar start stepped along the umbrella mode leaves the plane:
    # the largest per-atom step is the requested amplitude and it is
    # out of plane.
    shift = stepped - before
    assert abs(np.linalg.norm(shift, axis=1).max() - 0.1) < 1e-6
    assert np.abs(shift[:, 2]).max() > 0.05

    # The other side of the saddle is the opposite sign.
    other, other_receipt = displace_trusted_geometry_along_mode(
        approved_workspace=tmp_path,
        displaced_artifact_id="nh3-stepped-back",
        result_artifact=_artifact(),
        program="pyscf",
        mode_index=1,
        amplitude_angstrom=-0.1,
    )
    lines = Path(other.path).read_text().splitlines()
    back = np.asarray(
        [[float(v) for v in line.split()[1:4]] for line in lines[2:6]]
    )
    assert np.allclose(back - before, -(stepped - before), atol=1e-5)


def test_the_kind_is_the_readers_word_not_the_programs_name():
    wrong = TrustedArtifactRefV1(
        artifact_id="pyscf-result-saddle",
        kind="pyscf_output",
        sha256=file_sha256(SADDLE),
        size_bytes=SADDLE.stat().st_size,
        path=str(SADDLE),
        cli_value=str(SADDLE),
    )
    with pytest.raises(ContractError, match="pyscf_hdf5"):
        build_stationary_point_characterisation(
            result_artifact=wrong, program="pyscf", order_claimed=1
        )
