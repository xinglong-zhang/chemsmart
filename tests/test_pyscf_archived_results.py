"""Every declared PySCF selector is exercised on bytes PySCF wrote.

The archived results under ``tests/data/PySCFTests/outputs`` are real
PySCF 2.14 runs with their receipts and a ``reference.json`` PySCF itself
wrote (see the README there).  A synthetic fixture is the writer reading
itself; these tests are the referential half: what the reader declares
must hold on the program's own artifacts, and the structural facts --
supplied against reached, converged against stopped, open shell against
closed -- are checked on fixtures built to differ.
"""

import inspect
import json
from pathlib import Path

import numpy as np
import pytest

from chemsmart.agent._contracts import file_sha256
from chemsmart.agent.driver import REPAIR_MENU
from chemsmart.agent.execution import (
    TrustedArtifactRefV1,
    build_reached_geometry,
)
from chemsmart.agent.terminal_states import (
    REPAIRABLE_NODE_STATES,
    _classify_failure,
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
from chemsmart.io.native_failure import summarize_pyscf_native_failure
from chemsmart.jobs.pyscf.settings import PySCFJobSettings

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
    "water_hess_historical_unclassified": ("hess", "water_hess_gas_phase.h5"),
    # Contract v5 (expansion round): the response, excited-surface and
    # correlated stages on real bytes.
    "water_td_singlet": ("td", "water_td_singlet_gas_phase.h5"),
    "water_td_triplet": ("td", "water_td_triplet_gas_phase.h5"),
    "water_td_rpa": ("td", "water_td_rpa_gas_phase.h5"),
    "hydroxyl_td_unrestricted": (
        "td",
        "hydroxyl_td_unrestricted_gas_phase.h5",
    ),
    "water_td_cpcm_toluene": ("td", "water_td_cpcm_toluene_cpcm_toluene.h5"),
    "water_td_unconverged": ("td", "water_td_unconverged_gas_phase.h5"),
    "formaldehyde_s1_opt": ("opt", "formaldehyde_s1_opt_gas_phase.h5"),
    "formaldehyde_s1_opt_planar": (
        "opt",
        "formaldehyde_s1_opt_planar_gas_phase.h5",
    ),
    "water_s1_opt_degenerate": ("opt", "water_s1_opt_degenerate_gas_phase.h5"),
    "formaldehyde_s1_td": ("td", "formaldehyde_s1_td_gas_phase.h5"),
    "water_mp2_sp": ("sp", "water_mp2_sp_gas_phase.h5"),
    "water_mp2_sp_fc1": ("sp", "water_mp2_sp_fc1_gas_phase.h5"),
    "water_ccsd_sp": ("sp", "water_ccsd_sp_gas_phase.h5"),
    "water_ccsdt_sp": ("sp", "water_ccsdt_sp_gas_phase.h5"),
    "hydroxyl_ump2_sp": ("sp", "hydroxyl_ump2_sp_gas_phase.h5"),
    "water_mp2_opt": ("opt", "water_mp2_opt_gas_phase.h5"),
    "water_ccsd_opt": ("opt", "water_ccsd_opt_gas_phase.h5"),
    "water_ccsd_unconverged": ("sp", "water_ccsd_unconverged_gas_phase.h5"),
}
GREEN = {
    "water_sp",
    "water_opt",
    "water_hess",
    "nh3_planar_opt",
    "nh3_planar_hess",
    "hydroxyl_sp",
    "water_stretched_hess",
    "water_hess_historical_unclassified",
    "water_td_singlet",
    "water_td_triplet",
    "water_td_rpa",
    "hydroxyl_td_unrestricted",
    "water_td_cpcm_toluene",
    "formaldehyde_s1_opt",
    "formaldehyde_s1_opt_planar",
    "water_s1_opt_degenerate",
    "formaldehyde_s1_td",
    "water_mp2_sp",
    "water_mp2_sp_fc1",
    "water_ccsd_sp",
    "water_ccsdt_sp",
    "hydroxyl_ump2_sp",
    "water_mp2_opt",
    "water_ccsd_opt",
}
TD_CASES = (
    "water_td_singlet",
    "water_td_triplet",
    "water_td_rpa",
    "hydroxyl_td_unrestricted",
    "water_td_cpcm_toluene",
    "formaldehyde_s1_td",
)
EXCITED_OPT_CASES = (
    "formaldehyde_s1_opt",
    "formaldehyde_s1_opt_planar",
    "water_s1_opt_degenerate",
)
CORRELATED_CASES = (
    "water_mp2_sp",
    "water_mp2_sp_fc1",
    "water_ccsd_sp",
    "water_ccsdt_sp",
    "hydroxyl_ump2_sp",
    "water_mp2_opt",
    "water_ccsd_opt",
)
HARTREE_TO_EV = 27.211386245988


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


def test_a_historical_unclassified_receipt_still_admits_its_result():
    """The runner used to leave a Hessian ``unclassified`` when no per-plan
    policy told it what to expect; that policy was retired and the word is
    no longer produced, but every receipt that carries it was green on
    every invariant and stays analysis-ready."""

    from chemsmart.analysis.result_quantities import (
        result_file_sha256,
        validate_pyscf_analysis_artifact,
    )

    path = _path("water_hess_historical_unclassified")
    receipt = json.loads(path.with_suffix(".receipt.json").read_text())
    assert (receipt["state"], receipt["scientific_validation_state"]) == (
        "engine_complete",
        "unclassified",
    )
    output, admitted = validate_pyscf_analysis_artifact(
        path, expected_sha256=result_file_sha256(path)
    )
    assert admitted["result_sha256"] == result_file_sha256(path)
    current = json.loads(
        _path("water_hess").with_suffix(".receipt.json").read_text()
    )
    assert (current["state"], current["scientific_validation_state"]) == (
        "validated",
        "validated",
    ), "the shipped runner validates a green Hessian outright"


# ----------------------------------------------------------------------
# contract v5: the response stage on real bytes
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:td:excitation_energies")
@pytest.mark.capability("selector:pyscf:td:oscillator_strengths")
@pytest.mark.parametrize("case", TD_CASES)
def test_the_response_stage_serves_ascending_roots_pyscf_recomputes(case):
    """Roots are ascending indices within one manifold at this geometry,
    and PySCF's own recomputation from the applied spec returns them."""

    reader = reader_for("pyscf")
    output = _open(case)
    reference = _reference(case)
    excitations, unit = reader.read(output, "excitation_energies")
    assert unit == "Eh", "PySCF states its native unit; display is eV"
    assert excitations == sorted(excitations)
    tolerance = 1e-4 if case == "hydroxyl_td_unrestricted" else 1e-10
    assert np.allclose(
        excitations, reference["td_excitation_energies_eh"], atol=tolerance
    )
    strengths = reader.read(output, "oscillator_strengths")[0]
    # The UKS radical re-converges to an SCF a few 1e-5 Eh apart, which
    # moves its strengths by 1e-7; the closed shells recompute exactly.
    assert np.allclose(
        strengths,
        reference["td_oscillator_strengths"],
        atol=1e-5 if case == "hydroxyl_td_unrestricted" else 1e-8,
    )
    assert reader.read(output, "excited_state_indices")[0] == [1, 2, 3]
    assert reader.read(output, "excited_state_manifold_roots")[0] == [
        1,
        2,
        3,
    ]
    assert reader.read(output, "excited_state_converged")[0] == [1, 1, 1]
    dipoles = np.asarray(reader.read(output, "transition_dipole_moments")[0])
    assert dipoles.shape == (3, 3)
    # An eigenvector's sign is arbitrary, so a transition dipole is
    # compared per root by magnitude and by component up to sign; on the
    # linear radical the degenerate pi pair makes the direction within
    # the perpendicular plane arbitrary too, so only the magnitude holds.
    recomputed = np.asarray(reference["td_transition_dipole_moments_debye"])
    tolerance = 1e-3 if case == "hydroxyl_td_unrestricted" else 1e-6
    if case != "hydroxyl_td_unrestricted":
        assert np.allclose(np.abs(dipoles), np.abs(recomputed), atol=tolerance)
    assert np.allclose(
        np.linalg.norm(dipoles, axis=1),
        np.linalg.norm(recomputed, axis=1),
        atol=tolerance,
    )
    # A td result is one structure: the supplied geometry, held to 1e-8 A.
    supplied = np.asarray(reader.read(output, "supplied_positions")[0])
    final = np.asarray(reader.read(output, "positions")[0])
    assert np.allclose(supplied, final, atol=1e-8)
    assert "reached_positions" not in reader.selectors_for_jobtype("td")
    # ``energy`` on a spectrum is the reference, as ORCA's td serves it.
    assert (
        reader.read(output, "energy")[0]
        == reader.read(output, "scf_energy")[0]
    )


@pytest.mark.capability("selector:pyscf:td:singlet_excitation_energies")
@pytest.mark.capability("selector:pyscf:td:triplet_excitation_energies")
def test_the_manifold_selectors_answer_their_own_manifold_only():
    reader = reader_for("pyscf")
    singlet = _open("water_td_singlet")
    triplet = _open("water_td_triplet")
    radical = _open("hydroxyl_td_unrestricted")
    values, unit = reader.read(singlet, "singlet_excitation_energies")
    assert unit == "eV" and abs(values[0] - 7.552) < 1e-2
    assert reader.read(singlet, "excited_state_multiplicities")[0] == [1] * 3
    assert reader.read(triplet, "excited_state_multiplicities")[0] == [3] * 3
    assert reader.read(triplet, "triplet_oscillator_strengths")[0] == [0.0] * 3
    with pytest.raises(MissingQuantityError):
        reader.read(singlet, "triplet_excitation_energies")
    with pytest.raises(MissingQuantityError):
        reader.read(triplet, "singlet_excitation_energies")
    # An unrestricted reference has one spin-conserving manifold PySCF
    # labels neither: no multiplicities, no per-root <S^2>, both manifold
    # selectors refused, the aggregate ones served.
    assert radical.state_manifold == "unrestricted"
    for selector in (
        "excited_state_multiplicities",
        "singlet_excitation_energies",
        "triplet_excitation_energies",
    ):
        with pytest.raises(MissingQuantityError):
            reader.read(radical, selector)
    assert len(reader.read(radical, "oscillator_strengths")[0]) == 3
    assert "excited_state_spin_square" not in reader.selectors_for_jobtype(
        "td"
    ), "PySCF prints no per-root <S^2>; the selector is honestly undeclared"
    assert "absorption_wavelengths" not in reader.selectors_for_jobtype(
        "td"
    ), "one conversion authority: the photon_wavelength operation"


def test_the_solvated_response_records_the_dielectric_it_applied():
    """PySCF's non-equilibrium response uses eps 1.78 for every solvent;
    the artifact says so beside the static dielectric it applied."""

    output = _open("water_td_cpcm_toluene")
    solvent = output.td_stage["solvent"]
    assert abs(solvent["static_eps_applied"] - 2.3741) < 1e-3
    assert abs(solvent["response_eps_applied"] - 1.78) < 1e-9
    assert solvent["equilibrium_solvation"] is False
    assert (
        abs(_reference("water_td_cpcm_toluene")["td_response_eps"] - 1.78)
        < 1e-9
    )
    assert output.solvent_model == "cpcm" and output.solvent_id == "toluene"


@pytest.mark.capability("selector:pyscf:opt:excited_state_followed_root")
@pytest.mark.parametrize("case", EXCITED_OPT_CASES)
def test_an_excited_root_optimisation_names_its_root_and_its_gaps(case):
    """The surface an excited-root optimisation minimised is root k by
    index; the spectrum is re-evaluated at the reached geometry so every
    excited quantity belongs to the one structure, and the gaps the
    sensor step reports are recorded with their numbers."""

    reader = reader_for("pyscf")
    output = _open(case)
    reference = _reference(case)
    assert reader.read(output, "excited_state_followed_root")[0] == 1
    assert reader.read(output, "converged")[0] == 1
    energy = reader.read(output, "energy")[0]
    scf = reader.read(output, "scf_energy")[0]
    excitations = reader.read(output, "excitation_energies")[0]
    assert abs(energy - (scf + excitations[0])) < 1e-9
    record = output.excited_state_record
    assert abs(record["followed_root_total_energy_end"] - energy) < 1e-12
    assert abs(reference["followed_root_total_energy_eh"] - energy) < 1e-5
    assert record["final_gradient_max_eh_per_bohr"] < 4.5e-4
    assert (
        abs(
            record["root_gap_to_ground_end_ev"]
            - excitations[0] * HARTREE_TO_EV
        )
        < 1e-6
    )
    # Provenance: the total is the root's, the dipole is the reference's.
    assert reader.electronic_provenance_for_output(output, "energy") == (
        "excited_root"
    )
    assert reader.electronic_provenance_for_output(
        output, "dipole_moment"
    ) == ("reference")
    assert reader.electronic_provenance_for_output(output, "scf_energy") == (
        "reference"
    )


def test_the_three_excited_minima_are_three_different_observations():
    bent = _open("formaldehyde_s1_opt").excited_state_record
    planar = _open("formaldehyde_s1_opt_planar").excited_state_record
    water = _open("water_s1_opt_degenerate").excited_state_record
    # root == nstates: the case where PySCF 2.14's scanner property raises.
    assert bent["nstates"] == 1 and bent["root"] == 1
    assert bent["root_gap_to_neighbour_end_ev"] is None
    # The planar start stays planar: the symmetric stationary point no
    # Hessian here can characterise, with its neighbour far away.
    assert abs(planar["root_gap_to_ground_end_ev"] - 3.406) < 5e-3
    assert planar["root_gap_to_neighbour_end_ev"] > 4.0
    # The followed root ends degenerate with its neighbour: the index no
    # longer identifies a state, which the gap says without a verdict.
    assert water["root_gap_to_neighbour_end_ev"] < 1e-4
    assert water["roots_filtered_end"] == 0


# ----------------------------------------------------------------------
# contract v5: the correlated stage on real bytes
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:sp:correlation_energy")
@pytest.mark.capability("selector:pyscf:sp:reference_energy")
@pytest.mark.parametrize("case", CORRELATED_CASES)
def test_the_correlated_stage_serves_components_pyscf_recomputes(case):
    reader = reader_for("pyscf")
    output = _open(case)
    reference = _reference(case)
    energy = reader.read(output, "energy")[0]
    ref = reader.read(output, "reference_energy")[0]
    corr = reader.read(output, "correlation_energy")[0]
    assert abs(energy - (ref + corr)) < 1e-10
    assert abs(ref - reader.read(output, "scf_energy")[0]) < 1e-12
    assert abs(corr - reference["corr_correlation_energy_eh"]) < 1e-7
    assert output.frozen_core_applied == reference["corr_frozen_orbitals"]
    assert reader.electronic_provenance_for_output(output, "energy") == (
        "correlated"
    )
    assert reader.electronic_provenance_for_output(
        output, "dipole_moment"
    ) == ("reference")
    method = reader.read(output, "ab_initio")[0]
    if method == "mp2":
        for selector in ("ccsd_correlation_energy", "triples_correction"):
            with pytest.raises(MissingQuantityError):
                reader.read(output, selector)
    with pytest.raises(MissingQuantityError):
        reader.read(output, "functional")


@pytest.mark.capability("selector:pyscf:sp:triples_correction")
@pytest.mark.capability("selector:pyscf:sp:ccsd_correlation_energy")
def test_ccsd_t_correlation_is_ccsd_plus_triples_as_orca_means_it():
    reader = reader_for("pyscf")
    output = _open("water_ccsdt_sp")
    ccsd = reader.read(output, "ccsd_correlation_energy")[0]
    triples = reader.read(output, "triples_correction")[0]
    total = reader.read(output, "correlation_energy")[0]
    assert abs(total - (ccsd + triples)) < 1e-12
    assert triples < 0 and abs(triples) < 0.01
    assert output.frozen_core_applied == 1, "auto froze the oxygen 1s"
    assert output.status["stages"]["corr"]["frozen_core_requested"] == "auto"
    plain = _open("water_ccsd_sp")
    with pytest.raises(MissingQuantityError):
        reader.read(plain, "triples_correction")
    assert (
        abs(
            reader.read(plain, "ccsd_correlation_energy")[0]
            - reader.read(plain, "correlation_energy")[0]
        )
        < 1e-12
    )


def _orca_mp2_correlation(name):
    text = (FIXTURES / "orca_differential" / f"{name}.out").read_text()
    for line in text.splitlines():
        if "MP2 CORRELATION ENERGY" in line:
            return float(line.split(":")[1].split()[0])
    raise AssertionError(f"no MP2 correlation energy in {name}.out")


def test_frozen_core_is_a_convention_both_programs_agree_on_once_named():
    """ORCA freezes core by default and PySCF correlates every electron;
    matched conventions agree to 5e-8 Eh on the same relaxed water."""

    reader = reader_for("pyscf")
    frozen = reader.read(_open("water_mp2_sp_fc1"), "correlation_energy")[0]
    all_electron = reader.read(_open("water_mp2_sp"), "correlation_energy")[0]
    assert abs(frozen - _orca_mp2_correlation("water_mp2_fc")) < 1e-6
    assert abs(all_electron - _orca_mp2_correlation("water_mp2_nofc")) < 1e-6
    assert abs(frozen - all_electron) > 1e-3, "the convention moves 2.4 mEh"
    assert _open("water_mp2_sp").frozen_core_applied == 0
    assert _open("water_mp2_sp_fc1").frozen_core_applied == 1


def _orca_tda_spectrum(name):
    import re

    text = (FIXTURES / "orca_differential" / f"{name}.out").read_text()
    pattern = re.compile(
        r"^\s*0-1A\s+->\s+\d+-1A\s+([\d.]+)\s+[\d.]+\s+[\d.]+\s+([\d.]+)"
    )
    rows = []
    for line in text.splitlines():
        match = pattern.match(line)
        if match:
            rows.append((float(match.group(1)), float(match.group(2))))
    assert rows, f"no TDA spectrum in {name}.out"
    return rows[:3]


def test_the_orca_tda_differential_agrees_within_the_measured_band():
    """B3LYP/G in ORCA is the VWN3 functional PySCF's b3lypg names; on one
    geometry the two TDA spectra agree to 2 meV and 1e-3 in strength.
    The band is what was measured, recorded here so drift is visible."""

    reader = reader_for("pyscf")
    output = _open("water_td_singlet")
    pyscf_ev = reader.read(output, "singlet_excitation_energies")[0]
    pyscf_f = reader.read(output, "singlet_oscillator_strengths")[0]
    orca = _orca_tda_spectrum("water_td_b3lypg")
    for (orca_ev, orca_f), ev, f in zip(orca, pyscf_ev, pyscf_f):
        assert abs(orca_ev - ev) < 2e-3, (orca_ev, ev)
        assert abs(orca_f - f) < 1e-3, (orca_f, f)


# ----------------------------------------------------------------------
# an excited minimum is a structure producer: its reached geometry feeds
# a response consumer through the handoff every optimisation uses
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:opt:reached_positions")
@pytest.mark.capability("selector:pyscf:td:excitation_energies")
def test_an_excited_minimum_hands_its_geometry_to_a_response_consumer(
    tmp_path,
):
    """``formaldehyde_s1_td`` is a real run on the XYZ the host's own
    reached-geometry route wrote from ``formaldehyde_s1_opt``.  The
    validated handoff materialises that geometry from the producer's
    bytes under the identity gate; the consumer computed at it and holds
    to it; and its first root is the gap the producer recorded at its
    end -- the emission energy, read where the excited surface reached,
    on a reference energy the two runs agree on to the last digit."""

    from chemsmart.agent.execution import (
        build_program_execution_invocation,
        build_program_execution_receipt,
        handoff_optimized_pyscf_geometry,
    )
    from tests.agent.test_program_execution import _artifact as _bound_artifact
    from tests.agent.test_program_execution import (
        _test_approval,
        _test_resources,
    )

    approval = _test_approval(tmp_path)
    opt_node = approval.node("opt-initial")
    invocation = build_program_execution_invocation(
        node_id=opt_node.node_id,
        approval=approval,
        project_artifact=_bound_artifact(
            tmp_path / "water-pyscf.yaml",
            artifact_id="project.water.pyscf",
            kind="project_yaml",
        ),
        input_artifact=_bound_artifact(
            tmp_path / "water.xyz",
            artifact_id="geometry.water.initial",
            kind="geometry_xyz",
        ),
        scientific_identity_sha256=opt_node.scientific_identity_sha256,
        environment_receipt_sha256="b" * 64,
        resources=_test_resources(),
        argv=("chemsmart", "run", "pyscf", "opt"),
    )
    producer = _artifact("formaldehyde_s1_opt", artifact_id="result.s1.hdf5")
    receipt = build_program_execution_receipt(
        invocation,
        execution_state="validated",
        exit_status=0,
        engine_complete=True,
        validated=True,
        output_artifacts=(producer,),
        validator_receipt_sha256s=("e" * 64,),
        result_validation_receipt_sha256="e" * 64,
        started_at="2026-09-13T00:00:00+00:00",
        finished_at="2026-09-13T00:00:01+00:00",
    )
    start = FIXTURES / "inputs" / "formaldehyde_bent_start.xyz"
    supplied = TrustedArtifactRefV1(
        artifact_id="geometry.h2co.start",
        kind="geometry_xyz",
        sha256=file_sha256(start),
        size_bytes=start.stat().st_size,
        path=str(start),
        cli_value=str(start),
    )
    geometry, handoff = handoff_optimized_pyscf_geometry(
        producer_receipt=receipt,
        result_artifact=producer,
        producer_edge=approval.producer_edges[0],
        approved_workspace=tmp_path,
        geometry_artifact_id="geometry.h2co.s1",
        expected_charge=0,
        expected_multiplicity=1,
        input_artifact=supplied,
    )
    assert handoff.status == "validated_handoff"
    assert handoff.symbols == ("C", "O", "H", "H")
    text = Path(geometry.path).read_text(encoding="utf-8")
    assert "source_sha256=" + producer.sha256 in text
    carried = np.asarray(
        [
            [float(v) for v in line.split()[1:4]]
            for line in text.splitlines()[2:6]
        ]
    )

    reader = reader_for("pyscf")
    opt = _open("formaldehyde_s1_opt")
    td = _open("formaldehyde_s1_td")
    assert np.allclose(carried, np.asarray(opt.positions), atol=1e-9)
    assert np.allclose(np.asarray(td.supplied_positions), carried, atol=1e-9)
    assert np.allclose(
        np.asarray(td.positions), np.asarray(td.supplied_positions), atol=1e-12
    ), "a response stage holds to the geometry it was handed"
    emission_ev = td.excitation_energies[0] * HARTREE_TO_EV
    end_gap_ev = opt.excited_state_record["root_gap_to_ground_end_ev"]
    # Two independent response solves at one geometry -- the producer's
    # from the optimiser's own density, the consumer's from scratch --
    # agree to 7e-6 eV (2.5e-7 Eh), measured; the band is ten times that.
    assert abs(emission_ev - end_gap_ev) < 1e-4, (emission_ev, end_gap_ev)
    assert abs(td.scf_energy - opt.scf_energy) < 1e-9
    root_total = td.scf_energy + td.excitation_energies[0]
    assert abs(root_total - opt.total_energy) < 1e-4 / HARTREE_TO_EV, (
        "the producer's total is the root's total; the consumer's is the "
        "reference's, and the root's total is rebuilt from it"
    )
    assert reader.electronic_provenance_for_output(opt, "energy") == (
        "excited_root"
    )
    assert reader.electronic_provenance_for_output(td, "energy") == "reference"


@pytest.mark.capability("tool:bind_reached_geometry")
def test_the_reached_geometry_of_an_excited_minimum_is_what_its_consumer_ran_on(
    tmp_path,
):
    artifact = _artifact("formaldehyde_s1_opt")
    geometry, receipt = build_reached_geometry(
        approved_workspace=tmp_path,
        reached_artifact_id="reached-s1",
        result_artifact=artifact,
        program="pyscf",
    )
    assert receipt.normal_termination is True
    archived = FIXTURES / "inputs" / "formaldehyde_s1_reached.xyz"
    assert (
        Path(geometry.path).read_text().splitlines()[2:]
        == archived.read_text().splitlines()[2:]
    ), "the archived consumer input is this route's own output"


# ----------------------------------------------------------------------
# the level a result computed at is read from its own record
# ----------------------------------------------------------------------


@pytest.mark.capability("tool:inspect_run")
def test_the_level_names_the_convention_and_the_root():
    """Method and basis alone do not name a level: a correlated result's
    frozen-core count and an excited-surface result's response, manifold
    and followed root are what make two results at one functional and one
    basis different calculations, and the inspection reply says so from
    the artifact rather than from a project the session may not hold."""

    reader = reader_for("pyscf")
    assert reader.level_for_output(_open("water_sp")) == {
        "functional": "b3lyp",
        "basis": "def2-svp",
    }
    assert reader.level_for_output(_open("water_mp2_sp")) == {
        "ab_initio": "mp2",
        "basis": "def2-svp",
        "frozen_core": 0,
    }, "PySCF's all-electron default is a level, never an absence"
    assert reader.level_for_output(_open("water_ccsdt_sp")) == {
        "ab_initio": "ccsd(t)",
        "basis": "def2-svp",
        "frozen_core": 1,
    }, "'auto' is displayed as the count it applied"
    assert reader.level_for_output(_open("water_td_cpcm_toluene")) == {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "solvent_model": "cpcm",
        "solvent": "toluene",
        "response_method": "tda",
        "state_manifold": "singlet",
        "nstates": 3,
    }
    assert reader.level_for_output(_open("formaldehyde_s1_opt")) == {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "response_method": "tda",
        "state_manifold": "singlet",
        "nstates": 1,
        "excited_state_root": 1,
    }
    assert reader_for("orca").level_for_output(object()) == {}


# ----------------------------------------------------------------------
# an unconverged root or amplitude set is typed and is not evidence
# ----------------------------------------------------------------------


@pytest.mark.capability("setting:pyscf:td_max_cycle")
@pytest.mark.capability("setting:pyscf:cc_max_cycle")
@pytest.mark.parametrize(
    ("case", "stage"),
    [("water_td_unconverged", "td"), ("water_ccsd_unconverged", "corr")],
)
def test_an_unconverged_stage_is_inspectable_and_not_evidence(case, stage):
    from chemsmart.analysis.result_quantities import (
        result_file_sha256,
        validate_pyscf_analysis_artifact,
    )

    output = _open(case)
    assert output.normal_termination is False
    assert output.status["stages"][stage]["converged"] is False
    assert output.excitation_energies or output.correlation_energy is not None
    path = _path(case)
    with pytest.raises(QuantityExtractionError):
        validate_pyscf_analysis_artifact(
            path, expected_sha256=result_file_sha256(path)
        )
    if stage == "td":
        assert output.excited_state_converged == [False, False, False]
        assert output.td_stage["unconverged_roots"] == [1, 2, 3]
        assert output.td_stage["max_cycle_applied"] == 1
        native_class, word, control = (
            "excited_state_convergence",
            "failed_nonconverged_excited_state",
            "td_max_cycle",
        )
    else:
        assert output.status["stages"]["corr"]["max_cycle_applied"] == 1
        native_class, word, control = (
            "correlation_convergence",
            "failed_nonconverged_correlation",
            "cc_max_cycle",
        )
    # The word is derived from the artifact's own typed status -- the
    # summariser reads the stage flag, the classifier maps the class, and
    # the menu names the public control that answers it, which is a
    # control the settings object carries and the artifact recorded.
    summary = summarize_pyscf_native_failure(output.status)
    assert summary is not None and summary.error_class == native_class
    assert "max_cycle 1" in " ".join(summary.diagnostic_lines)
    assert (
        _classify_failure(
            jobtype=output.jobtype,
            findings=(f"pyscf.native_failure.{native_class}",),
            native_class=native_class,
            converged=None,
            reached=None,
            planned=None,
        )
        == word
    )
    assert word in REPAIRABLE_NODE_STATES
    assert control in REPAIR_MENU[word]
    assert control in inspect.signature(PySCFJobSettings.__init__).parameters
    assert output.spec[control] == 1


@pytest.mark.capability("gate:terminal_state_vocabulary")
def test_a_followed_root_that_vanishes_is_an_excited_state_convergence():
    """The driver raises ``FollowedRootFiltered`` inside the ``opt`` stage
    when the root it follows falls below PySCF's positive-eigenvalue
    filter, rather than switching roots; the summariser maps the exception
    by name, so the word is the response solver's and not the optimiser's.
    Every exception the summariser maps by name is one the generated
    driver defines."""

    from chemsmart.io.native_failure import _PYSCF_EXCEPTION_CLASSES
    from chemsmart.jobs.pyscf.writer import _SKELETON

    for name in _PYSCF_EXCEPTION_CLASSES:
        assert f"class {name}(" in _SKELETON, name
    summary = summarize_pyscf_native_failure(
        {
            "normal_termination": False,
            "failure": {
                "stage": "opt",
                "type": "FollowedRootFiltered",
                "message": "followed root 1 fell below the filter",
            },
            "stages": {"scf": {"converged": True}},
        }
    )
    assert summary.error_class == "excited_state_convergence"
    assert (
        _classify_failure(
            jobtype="opt",
            findings=("pyscf.native_failure.excited_state_convergence",),
            native_class=summary.error_class,
            converged=False,
            reached=None,
            planned=None,
        )
        == "failed_nonconverged_excited_state"
    )


# ----------------------------------------------------------------------
# every declared selector is requestable, dimensioned and provenanced
# ----------------------------------------------------------------------


@pytest.mark.capability("selector:pyscf:*")
def test_every_declared_pyscf_selector_is_requestable_and_provenanced():
    """A declaration is only reachable through the request gate, and a
    quantity that belongs to a density or a method says which."""

    from chemsmart.analysis import result_readers as readers_module
    from chemsmart.analysis.result_quantities import (
        SUPPORTED_SELECTORS,
        QuantitySelectorV1,
    )
    from chemsmart.analysis.result_readers import (
        ELECTRONIC_PROVENANCES,
        SELECTOR_UNITS,
    )

    reader = reader_for("pyscf")
    identity = {
        "ab_initio",
        "basis",
        "charge",
        "connectivity",
        "converged",
        "functional",
        "method",
        "multiplicity",
        "positions",
        "reached_positions",
        "supplied_positions",
        "symbols",
    }
    for jobtype, _selectors in reader.jobtype_selectors:
        for selector in reader.selectors_for_jobtype(jobtype):
            assert selector in SUPPORTED_SELECTORS, selector
            QuantitySelectorV1(quantity_id="q", selector=selector)
            assert selector in SELECTOR_UNITS, selector
            assert selector in readers_module._SELECTOR_DIMENSIONS, selector
            word = reader.electronic_provenance(selector)
            assert word in ELECTRONIC_PROVENANCES
            if selector not in identity:
                assert word != "stateless", (
                    f"{jobtype}:{selector} belongs to a density or a "
                    "method and declares neither"
                )
    assert reader.electronic_provenance("energy") == "computed_surface"


@pytest.mark.capability("tool:extract_result_quantities")
def test_the_extraction_receipt_carries_the_provenance_of_each_value():
    from chemsmart.analysis.result_quantities import (
        QuantityExtractionReceiptV1,
        QuantitySelectorV1,
        ResultQuantityExtractionRequestV1,
    )
    from chemsmart.analysis.result_readers import extract_logged_quantities

    for case, expected in (
        ("formaldehyde_s1_opt", "excited_root"),
        ("water_ccsdt_sp", "correlated"),
        ("water_sp", "reference"),
    ):
        artifact = _artifact(case)
        request = ResultQuantityExtractionRequestV1(
            schema_version="chemsmart.quantity-extraction-request.v1",
            program="pyscf",
            artifact_id=artifact.artifact_id,
            artifact_sha256=artifact.sha256,
            selectors=(
                QuantitySelectorV1(quantity_id="e", selector="energy"),
                QuantitySelectorV1(quantity_id="mu", selector="dipole_moment"),
                QuantitySelectorV1(quantity_id="sym", selector="symbols"),
            ),
        )
        receipt = extract_logged_quantities(
            request=request, artifact_path=artifact.path
        )
        assert dict(receipt.electronic_provenance) == {
            "e": expected,
            "mu": "reference",
        }, case
        # The provenance is inside the digest the receipt verifies.
        QuantityExtractionReceiptV1(
            **{
                key: getattr(receipt, key)
                for key in (
                    "schema_version",
                    "artifact_id",
                    "artifact_sha256",
                    "program",
                    "parser_id",
                    "quantities",
                    "status",
                    "receipt_sha256",
                    "absent",
                    "derived_adjacency",
                    "electronic_provenance",
                )
            }
        )
        with pytest.raises(Exception):
            QuantityExtractionReceiptV1(
                **{
                    key: getattr(receipt, key)
                    for key in (
                        "schema_version",
                        "artifact_id",
                        "artifact_sha256",
                        "program",
                        "parser_id",
                        "quantities",
                        "status",
                        "receipt_sha256",
                        "absent",
                        "derived_adjacency",
                    )
                },
                electronic_provenance=(),
            )
