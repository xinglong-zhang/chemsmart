"""Every declared PySCF selector is exercised on bytes PySCF wrote.

The archived results under ``tests/data/PySCFTests/outputs`` are real
PySCF 2.14 runs with their receipts and a ``reference.json`` PySCF itself
wrote (see the README there).  A synthetic fixture is the writer reading
itself; these tests are the referential half: what the reader declares
must hold on the program's own artifacts, and the structural facts --
supplied against reached, converged against stopped, open shell against
closed -- are checked on fixtures built to differ.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.execution import (
    TrustedArtifactRefV1,
    build_reached_geometry,
)
from chemsmart.analysis.result_quantities import (
    QuantityExtractionError,
    ThermochemistryRequestV1,
    derive_result_thermochemistry,
)
from chemsmart.analysis.result_readers import (
    MissingQuantityError,
    reader_for,
)

FIXTURES = Path(__file__).resolve().parent / "data" / "PySCFTests" / "outputs"

CASES = {
    "water_sp": ("sp", "water_sp_gas_phase.h5"),
    "water_opt": ("opt", "water_opt_gas_phase.h5"),
    "water_hess": ("hess", "water_hess_gas_phase.h5"),
    "nh3_planar_opt": ("opt", "nh3_planar_opt_gas_phase.h5"),
    "nh3_planar_hess": ("hess", "nh3_planar_hess_gas_phase.h5"),
    "hydroxyl_sp": ("sp", "hydroxyl_sp_gas_phase.h5"),
    "water_opt_maxsteps1": ("opt", "water_opt_maxsteps1_gas_phase.h5"),
    "water_stretched_hess": ("hess", "water_stretched_hess_gas_phase.h5"),
    "water_sp_scfmaxiter2": ("sp", "water_sp_scfmaxiter2_gas_phase.h5"),
}
GREEN = {
    "water_sp",
    "water_opt",
    "water_hess",
    "nh3_planar_opt",
    "nh3_planar_hess",
    "hydroxyl_sp",
    "water_stretched_hess",
}


def _path(case):
    jobtype, name = CASES[case]
    return FIXTURES / case / name


def _open(case):
    return reader_for("pyscf").open_output(_path(case))


def _artifact(case, artifact_id="pyscf-result-archived"):
    path = _path(case)
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind="pyscf_hdf5",
        sha256=file_sha256(path),
        size_bytes=path.stat().st_size,
        path=str(path),
        cli_value=str(path),
    )


def _reference(case):
    jobtype, name = CASES[case]
    return json.loads(
        (FIXTURES / case / name.replace(".h5", ".reference.json")).read_text()
    )


# ----------------------------------------------------------------------
# every declaration reads on real bytes, or refuses as an absence
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:*")
@pytest.mark.parametrize("case", sorted(GREEN))
def test_every_declared_selector_reads_or_is_honestly_absent(case):
    """A declared selector never fails on a real artifact for a host reason.

    ``MissingQuantityError`` is the honest answer for a quantity this run
    did not produce (spin populations on a closed shell, a Hartree-Fock
    method on a DFT run).  A unit disagreement or a shape error is a
    defect in the host, and it is what hid
    ``vibrational_mode_atom_participation`` on every real PySCF Hessian.
    """

    reader = reader_for("pyscf")
    output = _open(case)
    jobtype = CASES[case][0]
    assert output.jobtype == jobtype
    problems = {}
    for selector in reader.selectors_for_jobtype(jobtype):
        try:
            reader.read(output, selector)
        except MissingQuantityError:
            continue
        except Exception as exc:  # noqa: BLE001 - the defect class under test
            problems[selector] = f"{type(exc).__name__}: {exc}"
    assert not problems, problems


@pytest.mark.capability(
    "selector:pyscf:hess:vibrational_mode_atom_participation"
)
def test_participation_rows_sum_to_one_on_a_real_hessian():
    reader = reader_for("pyscf")
    for case in ("water_hess", "nh3_planar_hess"):
        output = _open(case)
        shares = np.asarray(
            reader.read(output, "vibrational_mode_atom_participation")[0]
        )
        assert shares.shape == (
            len(output.vibrational_frequencies),
            len(output.chemical_symbols),
        )
        assert np.allclose(shares.sum(axis=1), 1.0, atol=1e-9)
    # The umbrella mode of planar ammonia lives on the hydrogens, out of
    # plane, with the nitrogen carrying the balance: an observation, not
    # a verdict, but one a reader can check against the geometry.
    saddle = _open("nh3_planar_hess")
    umbrella = np.asarray(
        reader.read(saddle, "vibrational_mode_atom_participation")[0]
    )[0]
    assert saddle.vibrational_frequencies[0] < -800.0
    assert umbrella[1:].sum() > 0.5


# ----------------------------------------------------------------------
# supplied against reached
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:opt:supplied_positions")
@pytest.mark.capability("selector:pyscf:opt:reached_positions")
def test_an_optimisation_that_moved_serves_two_structures():
    reader = reader_for("pyscf")
    output = _open("water_opt")
    supplied = np.asarray(reader.read(output, "supplied_positions")[0])
    reached = np.asarray(reader.read(output, "reached_positions")[0])
    final = np.asarray(reader.read(output, "positions")[0])
    assert np.abs(reached - supplied).max() > 0.05, "the run moved"
    assert np.array_equal(reached, final), "one structure per result"
    assert reader.structural_state("supplied_positions") == "as_supplied"
    assert reader.structural_state("reached_positions") == "as_reached"
    assert reader.structural_state("positions") == "as_reached"
    # The supplied structure is the input file's bytes.
    xyz = (
        (FIXTURES / "inputs" / "water_distorted.xyz").read_text().splitlines()
    )
    from_input = np.asarray(
        [[float(v) for v in line.split()[1:4]] for line in xyz[2:5]]
    )
    assert np.allclose(supplied, from_input, atol=1e-8)


@pytest.mark.capability("selector:pyscf:sp:supplied_positions")
@pytest.mark.capability("selector:pyscf:hess:supplied_positions")
@pytest.mark.parametrize("case", ["water_sp", "water_hess", "hydroxyl_sp"])
def test_a_fixed_geometry_stage_reaches_what_it_was_handed(case):
    reader = reader_for("pyscf")
    output = _open(case)
    supplied = np.asarray(reader.read(output, "supplied_positions")[0])
    final = np.asarray(reader.read(output, "positions")[0])
    assert np.allclose(supplied, final, atol=1e-8)
    assert "reached_positions" not in reader.selectors_for_jobtype(
        CASES[case][0]
    ), "a fixed-geometry stage declares no reached structure"


@pytest.mark.capability("selector:pyscf:opt:converged")
def test_converged_is_read_from_the_drivers_own_stage_status():
    reader = reader_for("pyscf")
    assert reader.read(_open("water_opt"), "converged")[0] == 1
    assert reader.read(_open("water_opt_maxsteps1"), "converged")[0] == 0
    with pytest.raises(MissingQuantityError):
        reader.read(_open("water_sp"), "converged")
    assert _open("water_opt").molecule.is_optimized_structure is True
    assert (
        _open("water_opt_maxsteps1").molecule.is_optimized_structure is False
    ), "an optimiser stopped on its step limit reached a structure, not an optimised one"


# ----------------------------------------------------------------------
# a failed run is inspectable, its structure bindable, its numbers not
# ----------------------------------------------------------------------


@pytest.mark.capability("tool:bind_reached_geometry")
def test_a_failed_optimisation_opens_and_its_last_geometry_carries_forward(
    tmp_path,
):
    from chemsmart.analysis.result_quantities import (
        QuantitySelectorV1,
        ResultQuantityExtractionRequestV1,
    )
    from chemsmart.analysis.result_readers import extract_logged_quantities

    reader = reader_for("pyscf")
    output = _open("water_opt_maxsteps1")
    assert output.normal_termination is False
    assert output.converged is False
    artifact = _artifact("water_opt_maxsteps1")
    request = ResultQuantityExtractionRequestV1(
        schema_version="chemsmart.quantity-extraction-request.v1",
        program="pyscf",
        artifact_id=artifact.artifact_id,
        artifact_sha256=artifact.sha256,
        selectors=(QuantitySelectorV1(quantity_id="e", selector="energy"),),
    )
    with pytest.raises(QuantityExtractionError, match="normally terminated"):
        extract_logged_quantities(request=request, artifact_path=artifact.path)

    geometry, receipt = build_reached_geometry(
        approved_workspace=tmp_path,
        reached_artifact_id="reached-from-maxsteps1",
        result_artifact=artifact,
        program="pyscf",
    )
    assert geometry.kind == "geometry_xyz"
    assert receipt.normal_termination is False
    lines = Path(geometry.path).read_text().splitlines()
    carried = np.asarray(
        [[float(v) for v in line.split()[1:4]] for line in lines[2:5]]
    )
    assert np.allclose(carried, output.positions, atol=1e-6)
    supplied = np.asarray(reader.read(output, "supplied_positions")[0])
    assert np.abs(carried - supplied).max() > 0.01, (
        "the reached structure is the last geometry the optimiser evaluated, "
        "never the input"
    )


def test_an_unconverged_scf_is_inspectable_and_not_evidence():
    output = _open("water_sp_scfmaxiter2")
    assert output.normal_termination is False
    assert output.energies, "the printed number stays readable"
    from chemsmart.io.native_failure import summarize_pyscf_native_failure

    summary = summarize_pyscf_native_failure(output.status)
    assert summary is not None and summary.error_class == "scf_convergence"


# ----------------------------------------------------------------------
# identities and spin
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:sp:functional")
@pytest.mark.capability("selector:pyscf:sp:ab_initio")
def test_functional_is_the_requested_name_and_hf_is_refused():
    reader = reader_for("pyscf")
    output = _open("water_sp")
    assert reader.read(output, "functional")[0] == "b3lyp"
    assert output.spec["xc"] == "b3lypg", "the applied literal stays on spec"
    with pytest.raises(MissingQuantityError):
        reader.read(output, "ab_initio")


@pytest.mark.capability("selector:pyscf:sp:mulliken_atomic_spin_populations")
def test_spin_populations_sum_to_two_s_and_are_refused_on_a_closed_shell():
    reader = reader_for("pyscf")
    radical = _open("hydroxyl_sp")
    populations = reader.read(radical, "mulliken_atomic_spin_populations")[0]
    assert len(populations) == 2
    assert abs(sum(populations) - 1.0) < 0.05
    assert populations[0] > 0.9, "the spin lives on oxygen"
    assert abs(reader.read(radical, "spin_square")[0] - 0.75) < 0.01
    with pytest.raises(MissingQuantityError):
        reader.read(_open("water_sp"), "mulliken_atomic_spin_populations")


# ----------------------------------------------------------------------
# PySCF's own account of the same bytes
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:hess:vibrational_frequencies")
@pytest.mark.parametrize(
    "case", ["water_hess", "nh3_planar_hess", "water_stretched_hess"]
)
def test_stored_frequencies_are_what_pyscf_recomputes_from_the_stored_hessian(
    case,
):
    reference = _reference(case)
    output = _open(case)
    stored = np.sort(np.asarray(output.vibrational_frequencies))
    recomputed = np.sort(
        np.asarray(reference["harmonic_analysis_freq_wavenumber"])
    )
    assert np.allclose(stored, recomputed, atol=1e-6)
    assert output.point_group == reference["point_group_detect_symm"]


def test_the_saddle_and_the_non_stationary_hessian_are_told_apart_by_the_gradient():
    """Zero imaginary modes is not stationarity: the gradient says which."""

    saddle = _open("nh3_planar_hess")
    stretched = _open("water_stretched_hess")
    assert sum(f < -20.0 for f in saddle.vibrational_frequencies) == 1
    assert all(f > 0.0 for f in stretched.vibrational_frequencies)
    saddle_gradient = saddle.status["stages"]["hess"][
        "max_abs_gradient_eh_per_bohr"
    ]
    stretched_gradient = stretched.status["stages"]["hess"][
        "max_abs_gradient_eh_per_bohr"
    ]
    assert saddle_gradient < 4.5e-4
    assert stretched_gradient > 10 * 4.5e-4
    assert stretched.forces is not None and stretched.forces.shape == (3, 3)


@pytest.mark.capability("tool:derive_thermochemistry")
def test_the_host_rrho_and_pyscf_thermo_agree_where_they_must():
    """Must agree: ZPE and the vibrational terms, from identical frequencies.

    May differ, and is only recorded: the translational and rotational
    terms, where this engine applies most-abundant masses and PySCF
    isotope-averaged ones.  The symmetry number is common-mode (both
    read PySCF's detect_symm) so agreement there proves nothing.
    """

    artifact = _artifact("water_hess")
    receipt = derive_result_thermochemistry(
        request=ThermochemistryRequestV1(
            schema_version="chemsmart.thermochemistry-request.v1",
            program="pyscf",
            artifact_id=artifact.artifact_id,
            artifact_sha256=artifact.sha256,
            temperature_k=298.15,
            pressure_atm=1.0,
        ),
        artifact_path=artifact.path,
    )
    values = {item.quantity_id: item for item in receipt.quantities}
    reference = _reference("water_hess")["thermo_298K_1atm"]
    zpe = values["zero_point_energy"]
    assert zpe.unit == "hartree"
    assert abs(zpe.value - reference["ZPE"][0]) < 1e-8
    assert any(
        "standard state" in line for line in receipt.assumptions
    ), receipt.assumptions
    assert any(
        "isotope-averaged" in line for line in receipt.assumptions
    ), "a PySCF receipt names the mass table behind its frequencies"
    assert reference["sym_number"][0] == 2
