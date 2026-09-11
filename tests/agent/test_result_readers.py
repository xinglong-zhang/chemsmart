"""Every program's results must enter the typed quantity chain.

Typed extraction was written against PySCF's structured HDF5.  While ORCA,
Gaussian and xTB had no reader, the only way to get a number out of them was
for a model to read it off a log -- the exact hallucination channel the
project-YAML hub exists to close.  A number a model typed is not a number
ChemSmart measured.

These tests run against real output fixtures, not synthetic strings, because
the thing being verified is that the shipped parsers answer the shared
selector vocabulary.
"""

import dataclasses
import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError, TrustedArtifactRefV1
from chemsmart.agent.capabilities import CapabilityQueryV1, query_capability
from chemsmart.agent.postprocessing import extract_trusted_result_quantities
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface
from chemsmart.analysis.result_quantities import (
    QuantitySelectorV1,
)
from chemsmart.analysis.result_readers import (
    RESULT_READERS,
    MissingQuantityError,
    reader_for,
    registered_reader_programs,
    registered_reader_selectors,
)

_GAUSSIAN_LOG = "tests/data/GaussianTests/boltzmann/udc3_mCF3_monomer_c1.log"
_GAUSSIAN_LOG_C4 = (
    "tests/data/GaussianTests/boltzmann/udc3_mCF3_monomer_c4.log"
)
_GAUSSIAN_TD_LOG = (
    "tests/data/GaussianTests/tddft/tddft_r1s50_gas_radical_anion.log"
)
_ORCA_DLPNO_LOG = "tests/data/ORCATests/outputs/water_dlpno_ccsdt_sp.out"
_ORCA_ERROR_LOG = "tests/data/ORCATests/error_files/GTOInt_error.out"

pytestmark = pytest.mark.capability("selector:*")


def _artifact(path, program, artifact_id="result"):
    resolved = Path(path).resolve()
    return TrustedArtifactRefV1(
        artifact_id=artifact_id,
        kind=reader_for(program).artifact_kind,
        sha256=hashlib.sha256(resolved.read_bytes()).hexdigest(),
        size_bytes=resolved.stat().st_size,
        path=str(resolved),
        cli_value=str(resolved),
    )


def _extract(path, program, selector, quantity_id="q"):
    return extract_trusted_result_quantities(
        artifact=_artifact(path, program, quantity_id),
        program=program,
        selectors=(
            QuantitySelectorV1(quantity_id=quantity_id, selector=selector),
        ),
    )


def test_every_result_program_is_registered():
    # PySCF is here because a structured HDF5 result answers the same
    # selector vocabulary through the same registry, job-type declaration
    # gate and discovery path as a parsed log.  It used to be a second
    # extraction plane reached by name, which is how a capability query
    # could report that this program answers no selector at all.
    assert registered_reader_programs() == (
        "gaussian",
        "orca",
        "pyscf",
        "xtb",
        "xyz",
    )
    assert reader_for("pyscf").artifact_kind == "pyscf_hdf5"
    assert registered_reader_selectors()["xyz"] == (
        "connectivity",
        "energy",
        "positions",
        "symbols",
        "trajectory_connectivity_changed",
        "trajectory_end_connectivity",
        "trajectory_end_positions",
        "trajectory_frame_count",
        "trajectory_start_connectivity",
        "trajectory_start_positions",
    )


def test_model_tool_surface_exposes_the_registered_result_plane():
    surface = build_command_compiled_tool_surface()
    tool = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "inspect_run"
    )
    # The per-program selector union is stated once, on the tool whose
    # job is listing selectors; extract_result_quantities points here.
    properties = tool["function"]["parameters"]["properties"]
    assert properties["program"]["enum"] == [
        "gaussian",
        "orca",
        "pyscf",
        "xtb",
        "xyz",
    ]
    assert (
        "xyz: connectivity, energy, positions, symbols"
        in properties["program"]["description"]
    )
    extract = next(
        item
        for item in surface.tool_definitions
        if item["function"]["name"] == "extract_result_quantities"
    )
    properties = extract["function"]["parameters"]["properties"]
    selectors = properties["selectors"]["items"]["properties"]["selector"][
        "enum"
    ]
    assert "excitation_energies" in selectors
    assert "oscillator_strengths" in selectors
    assert "connectivity" in selectors
    assert (
        "not an electronic bond-order"
        in properties["selectors"]["items"]["properties"]["selector"][
            "description"
        ]
    )


def test_xtb_capability_receipt_binds_exact_job_result_coverage():
    hess = query_capability(CapabilityQueryV1("xtb", "hess", "cpu"))
    opt = query_capability(CapabilityQueryV1("xtb", "opt", "cpu"))

    hess_coverage = hess.job_result_selector_coverage
    opt_coverage = opt.job_result_selector_coverage
    assert hess_coverage is not None
    assert opt_coverage is not None
    assert hess_coverage.artifact_kind == "xtb_output"
    assert hess_coverage.parser_id == "chemsmart.io.xtb.output.XTBOutput"
    assert {
        "dipole_moment",
        "dipole_moment_magnitude",
        "energy",
        "gap",
        "homo",
        "lumo",
        "vibrational_frequencies",
    } <= set(hess_coverage.selectors)
    assert "vibrational_frequencies" not in opt_coverage.selectors
    tampered_coverage = dataclasses.replace(
        hess_coverage,
        selectors=tuple(
            selector
            for selector in hess_coverage.selectors
            if selector != "gap"
        ),
    )
    with pytest.raises(ContractError, match="digest mismatch"):
        dataclasses.replace(
            hess, job_result_selector_coverage=tampered_coverage
        )


def test_registered_xyz_geometry_enters_the_typed_quantity_plane(tmp_path):
    endpoint = tmp_path / "endpoint.xyz"
    endpoint.write_text(
        "3\nCoordinates from ORCA-job endpoint E -1.632059341860\n"
        "H -1.8556865849 0.0 0.0\n"
        "H  0.5609660462 0.0 0.0\n"
        "H  1.2947205387 0.0 0.0\n",
        encoding="utf-8",
    )

    artifact = _artifact(endpoint, "xyz", "endpoint")
    receipt = extract_trusted_result_quantities(
        artifact=artifact,
        program="xyz",
        selectors=(
            QuantitySelectorV1(quantity_id="e", selector="energy"),
            QuantitySelectorV1(quantity_id="c", selector="connectivity"),
            QuantitySelectorV1(quantity_id="r", selector="positions"),
            QuantitySelectorV1(quantity_id="z", selector="symbols"),
        ),
    )

    values = {item.quantity_id: item for item in receipt.quantities}
    assert values["e"].value == pytest.approx(-1.632059341860)
    assert values["c"].value == (
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
    )
    assert values["r"].value[1][0] == pytest.approx(0.5609660462)
    assert values["z"].value == ("H", "H", "H")


def test_xyz_energy_refuses_an_unlabeled_negative_comment_number(tmp_path):
    """The reason still has to be stated; it now travels on the receipt.

    A comment line's bare negative number is a charge, not an energy, and
    saying so is the point.  Since a refusal settles per output rather than
    ending the extraction, the sentence reaches a reader as a named absence
    instead of as an exception -- which is the same claim, recorded where
    the evidence is rather than only where the traceback was.
    """

    endpoint = tmp_path / "charged.xyz"
    endpoint.write_text(
        "1\ncharge -1.0 multiplicity 2\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )

    receipt = _extract(endpoint, "xyz", "energy")
    assert receipt.status == "partial"
    assert not receipt.quantities
    ((quantity_id, selector, reason),) = receipt.absent
    assert selector == "energy"
    assert "explicit Energy" in reason

    positions = _extract(endpoint, "xyz", "positions").quantities[0]
    assert positions.value == ((0.0, 0.0, 0.0),)


@pytest.mark.parametrize("program", ("gaussian", "orca", "xtb"))
def test_the_analysis_registry_exposes_every_log_reader(program):
    reader = reader_for(program)
    assert reader.parser_id
    assert reader.artifact_kind
    assert "energy" in reader.selectors


def test_a_gaussian_energy_becomes_a_hash_bound_quantity():
    receipt = _extract(_GAUSSIAN_LOG, "gaussian", "energy")
    quantity = receipt.quantities[0]
    assert quantity.value == pytest.approx(-2189.63187379)
    assert quantity.unit == "hartree"
    assert receipt.parser_id.endswith("Gaussian16Output")
    # The value carries its own digest, so a later claim cannot quietly
    # substitute a different number for the one that was measured.
    assert len(quantity.value_sha256) == 64


def test_a_gaussian_free_energy_is_read_from_the_thermochemistry_block():
    receipt = _extract(_GAUSSIAN_LOG, "gaussian", "gibbs_free_energy")
    assert receipt.quantities[0].value == pytest.approx(-2189.409887)


def test_orca_post_hf_energy_preserves_total_and_reference_components():
    """A correlated ORCA result must not collapse to its SCF reference."""

    total = _extract(_ORCA_DLPNO_LOG, "orca", "energy", "total").quantities[0]
    scf = _extract(_ORCA_DLPNO_LOG, "orca", "scf_energy", "scf").quantities[0]
    reference = _extract(
        _ORCA_DLPNO_LOG,
        "orca",
        "reference_energy",
        "reference",
    ).quantities[0]
    correlation = _extract(
        _ORCA_DLPNO_LOG,
        "orca",
        "correlation_energy",
        "correlation",
    ).quantities[0]

    assert total.value == pytest.approx(-76.377481488944)
    assert scf.value == pytest.approx(-76.066479259)
    assert reference.value == pytest.approx(scf.value)
    # The tracked ORCA CBS fixture prints this native value explicitly.  The
    # selector must not reconstruct it from the finite-basis SCF record.
    assert correlation.value == pytest.approx(-0.311002230)
    assert correlation.value == pytest.approx(total.value - scf.value)
    assert total.value != pytest.approx(scf.value)


def test_orca_dft_d_total_scf_and_dispersion_are_distinct_components():
    path = "tests/data/ORCATests/outputs/KOH.out"
    total = _extract(path, "orca", "energy", "total").quantities[0]
    scf = _extract(path, "orca", "scf_energy", "scf").quantities[0]
    dispersion = _extract(
        path, "orca", "dispersion_energy", "dispersion"
    ).quantities[0]

    assert total.value == pytest.approx(-675.522805891018)
    assert scf.value == pytest.approx(-675.52250804211144)
    assert dispersion.value == pytest.approx(-0.000297849)
    assert total.value == pytest.approx(scf.value + dispersion.value)

    # A dispersion-corrected DFT total is not a correlation energy, and the
    # reason is now recorded as an absence rather than raised.
    refused = _extract(path, "orca", "correlation_energy", "correlation")
    assert refused.status == "partial"
    assert "post-SCF correlation" in refused.absent[0][2]


def test_orca_correlation_requires_a_native_printed_result_not_method_name():
    reader = reader_for("orca")

    class _MethodOnly:
        ab_initio = "DLPNO-CCSD(T)"
        contents = (
            "Total Energy       : -10.000000 Eh -272.1 eV",
            "FINAL SINGLE POINT ENERGY -10.200000",
        )
        final_dispersion_energy = None

    with pytest.raises(MissingQuantityError, match="explicit final"):
        reader.read(_MethodOnly(), "correlation_energy")


def test_orca_native_correlation_record_is_evidence_without_method_heuristic():
    reader = reader_for("orca")

    class _PrintedResult:
        ab_initio = None
        contents = (
            "Total Energy       : -10.000000 Eh -272.1 eV",
            "Final correlation energy ... -0.200000",
            "FINAL SINGLE POINT ENERGY -10.200000",
        )
        final_dispersion_energy = None

    value, unit = reader.read(_PrintedResult(), "correlation_energy")
    assert value == pytest.approx(-0.2)
    assert unit == "Eh"


@pytest.mark.parametrize(
    ("lines", "expected"),
    (
        (("RI-MP2 CORRELATION ENERGY: -0.123456 Eh",), -0.123456),
        (("E(CORR)(corrected) ... -0.211643112",), -0.211643112),
        (
            (
                "E(CORR)(corrected) ... -0.211643112",
                "Triples Correction (T) ... -0.002954307",
                "Final correlation energy ... -0.214597419",
            ),
            -0.214597419,
        ),
        (
            (
                "Final correlation energy ... -0.275446122",
                "Extrapolated CBS correlation energy (2/3) : "
                "-0.311002230 (-0.035556108)",
            ),
            -0.311002230,
        ),
    ),
)
def test_orca_native_correlation_record_precedence(lines, expected):
    class _PrintedResult:
        contents = lines

    value, unit = reader_for("orca").read(
        _PrintedResult(), "correlation_energy"
    )
    assert value == pytest.approx(expected)
    assert unit == "Eh"


def test_orca_compound_output_uses_last_chronological_correlation_record():
    """A prior job's CBS record cannot outrank a later native RI-MP2 result."""

    class _CompoundResult:
        contents = (
            "Final correlation energy ... -0.275446122",
            "Extrapolated CBS correlation energy (2/3) : "
            "-0.311002230 (-0.035556108)",
            "FINAL SINGLE POINT ENERGY -10.311002230",
            "---------------- later calculation ----------------",
            "RI-MP2 CORRELATION ENERGY: -0.123456 Eh",
            "FINAL SINGLE POINT ENERGY -20.123456000",
        )

    value, unit = reader_for("orca").read(
        _CompoundResult(), "correlation_energy"
    )

    assert value == pytest.approx(-0.123456)
    assert unit == "Eh"


def test_gaussian_missing_post_annihilation_spin_is_a_typed_absence():
    reader = reader_for("gaussian")

    class _NoPostAnnihilationSpin:
        spin_square_history = [
            {"before_annihilation": 2.1, "after_annihilation": None}
        ]

    with pytest.raises(
        MissingQuantityError,
        match=r"does not print <S\^2> after annihilation",
    ):
        reader.read(
            _NoPostAnnihilationSpin(), "spin_square_after_annihilation"
        )


def test_error_terminated_orca_output_cannot_supply_scientific_quantities():
    """A wrapper-created log is not a result when ORCA itself aborted."""

    with pytest.raises(Exception, match="normally terminated"):
        _extract(_ORCA_ERROR_LOG, "orca", "energy")


def test_gaussian_excited_state_results_enter_the_shared_quantity_plane():
    energies = _extract(
        _GAUSSIAN_TD_LOG,
        "gaussian",
        "excitation_energies",
        "excitation_energies",
    ).quantities[0]
    strengths = _extract(
        _GAUSSIAN_TD_LOG,
        "gaussian",
        "oscillator_strengths",
        "oscillator_strengths",
    ).quantities[0]
    wavelengths = _extract(
        _GAUSSIAN_TD_LOG,
        "gaussian",
        "absorption_wavelengths",
        "absorption_wavelengths",
    ).quantities[0]
    assert energies.source_value[0] == pytest.approx(0.7744)
    assert energies.source_unit == "eV"
    assert energies.unit == "hartree"
    assert strengths.value[0] == pytest.approx(0.0084)
    assert strengths.unit == "1"
    assert wavelengths.source_value[0] == pytest.approx(1601.13)
    assert wavelengths.unit == "angstrom"


def test_excitation_energy_converts_to_wavelength_without_model_math():
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionNodeV1,
        QuantityExpressionRequestV1,
        evaluate_quantity_expression,
    )

    energies = _extract(
        _GAUSSIAN_TD_LOG,
        "gaussian",
        "excitation_energies",
        "excitation_energies",
    ).quantities[0]
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="excitation_wavelengths",
            inputs=(energies,),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="wavelengths",
                    operation="photon_wavelength",
                    input_ids=("excitation_energies",),
                ),
            ),
            output_node_ids=("wavelengths",),
        )
    )
    assert receipt.outputs[0].unit == "angstrom"
    assert receipt.outputs[0].value[0] == pytest.approx(16011.3, rel=2e-4)


def test_extraction_refuses_an_artifact_of_the_wrong_kind():
    artifact = _artifact(_GAUSSIAN_LOG, "gaussian")
    wrong = TrustedArtifactRefV1(
        artifact_id=artifact.artifact_id,
        kind="pyscf_hdf5",
        sha256=artifact.sha256,
        size_bytes=artifact.size_bytes,
        path=artifact.path,
        cli_value=artifact.cli_value,
    )
    with pytest.raises(Exception, match="gaussian_output"):
        extract_trusted_result_quantities(
            artifact=wrong,
            program="gaussian",
            selectors=(
                QuantitySelectorV1(quantity_id="q", selector="energy"),
            ),
        )


def test_an_unregistered_program_is_refused_by_name():
    artifact = _artifact(_GAUSSIAN_LOG, "gaussian")
    with pytest.raises(
        Exception, match="no typed result reader is registered"
    ):
        extract_trusted_result_quantities(
            artifact=artifact,
            program="nciplot",
            selectors=(
                QuantitySelectorV1(quantity_id="q", selector="energy"),
            ),
        )


def test_a_quantity_the_run_never_produced_is_absent_not_a_crash():
    """A single point has no Gibbs energy; that is a fact, not a parse bug."""

    reader = reader_for("orca")

    class _NoThermo:
        gibbs_free_energy = None

    with pytest.raises(MissingQuantityError, match="contains no"):
        reader.read(_NoThermo(), "gibbs_free_energy")


def test_a_parser_exception_becomes_an_absent_quantity():
    """The parsers raise IndexError for a block the run never wrote."""

    reader = reader_for("orca")

    class _Empty:
        @property
        def thermochemistry_molecule(self):
            raise IndexError("list index out of range")

    with pytest.raises(MissingQuantityError, match="IndexError"):
        reader.read(_Empty(), "positions")


def test_a_selector_a_reader_does_not_implement_is_refused_by_name():
    """Distinct from absent: xTB implements no scan surface at all.

    The previous example selector was ``charge`` with a docstring claiming
    xTB exposes none -- the parser resolved it all along, and an expert
    review flagged the false comment. The subject of this test is the
    refusal message for a genuinely unimplemented selector.
    """

    with pytest.raises(ValueError, match="does not provide"):
        reader_for("xtb").read(object(), "scan_energies")


@pytest.mark.parametrize("program", sorted(RESULT_READERS))
def test_every_reader_declares_a_unit_for_each_selector_it_provides(program):
    from chemsmart.analysis.result_readers import SELECTOR_UNITS

    assert reader_for(program).selectors <= set(SELECTOR_UNITS)


def test_two_extracted_energies_drive_a_real_expression():
    """The point of registering readers: reaching a number ChemSmart measured."""

    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionNodeV1,
        QuantityExpressionRequestV1,
        evaluate_quantity_expression,
    )

    first = _extract(
        _GAUSSIAN_LOG, "gaussian", "gibbs_free_energy", "conf1"
    ).quantities[0]
    second = _extract(
        _GAUSSIAN_LOG_C4, "gaussian", "gibbs_free_energy", "conf4"
    ).quantities[0]
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="dG.conf4_minus_conf1",
            inputs=(first, second),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="dG",
                    operation="subtract",
                    input_ids=("conf4", "conf1"),
                ),
            ),
            output_node_ids=("dG",),
        )
    )
    # Arithmetic stays in canonical hartree by design; the display unit is the
    # claim contract's business, not the evaluator's.
    assert receipt.outputs[0].unit == "hartree"
    assert receipt.outputs[0].value == pytest.approx(
        second.value - first.value, abs=1e-12
    )


@pytest.mark.parametrize(
    ("selector", "expected"),
    (("homo", -14.5428), ("lumo", -6.0942), ("gap", 8.448655866329)),
)
def test_xtb_frontier_orbitals_enter_the_shared_quantity_plane(
    selector, expected
):
    path = "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
    quantity = _extract(path, "xtb", selector, selector).quantities[0]
    assert quantity.source_value == pytest.approx(expected)
    assert quantity.source_unit == "eV"
    assert quantity.unit == "hartree"


@pytest.mark.parametrize(
    ("selector", "expected"),
    (("homo", -12.077619), ("lumo", 2.453142), ("gap", 14.530761)),
)
def test_orca_frontier_orbitals_enter_the_shared_quantity_plane(
    selector, expected
):
    # ORCAOutput parsed these eigenvalues all along and no selector reached
    # them, so a session could run ORCA for a HOMO-LUMO gap, complete the
    # calculation, and bind nothing from it.
    path = "tests/data/ORCATests/outputs/CO2.out"
    quantity = _extract(path, "orca", selector, selector).quantities[0]
    assert quantity.source_value == pytest.approx(expected, abs=1e-4)
    assert quantity.source_unit == "eV"
    assert quantity.unit == "hartree"


def test_orca_frontier_orbitals_span_both_spin_channels():
    """An unrestricted frontier pair need not live in the alpha channel.

    This Fe(II) triplet keeps both frontier orbitals in beta.  Reading the
    alpha channel alone would report a 11.30 eV gap for a system whose real
    frontier separation is 9.25 eV, so the cross-channel extremum is the
    definition, and the spin-resolved selectors stay available for a question
    that is genuinely about one channel.
    """

    reader = reader_for("orca")
    output = reader.open_output(
        Path("tests/data/ORCATests/outputs/fe2_triplet.out")
    )
    assert output.is_unrestricted is True
    frontier = {
        selector: reader.read(output, selector)[0]
        for selector in (
            "homo",
            "lumo",
            "gap",
            "alpha_homo",
            "alpha_lumo",
            "beta_homo",
            "beta_lumo",
        )
    }
    assert frontier["homo"] == pytest.approx(frontier["beta_homo"])
    assert frontier["lumo"] == pytest.approx(frontier["beta_lumo"])
    assert frontier["alpha_homo"] < frontier["homo"]
    assert frontier["alpha_lumo"] > frontier["lumo"]
    assert frontier["gap"] == pytest.approx(
        frontier["lumo"] - frontier["homo"]
    )


def test_orca_declares_frontier_coverage_only_where_a_reference_converges():
    declared = dict(reader_for("orca").jobtype_selectors)
    family = {
        "alpha_homo",
        "alpha_lumo",
        "beta_homo",
        "beta_lumo",
        "gap",
        "homo",
        "lumo",
    }
    for jobtype in ("sp", "opt", "ts"):
        assert family <= set(declared[jobtype]), jobtype
    # A TD job's answer is its excited states; its ground-state orbitals are
    # not what it was run to establish, so the family stays undeclared there.
    assert not (family & set(declared["td"]))


def test_orca_irc_declares_its_endpoint_but_not_a_log_trajectory():
    """The path lives in a sidecar, so the log must not promise it.

    ORCA writes the reaction path to `_IRC_Full_trj.xyz` and leaves a single
    structure in the log, so the ORCA reader declares the endpoint it
    converged to plus the explicitly declared direction, and the trajectory
    family stays with the `xyz` reader that can actually answer it.
    """

    declared = dict(reader_for("orca").jobtype_selectors)
    assert "irc_direction" in declared["irc"]
    assert not [
        selector
        for selector in declared["irc"]
        if selector.startswith("trajectory_")
    ]
    assert {"trajectory_connectivity_changed", "trajectory_frame_count"} <= (
        reader_for("xyz").selectors
    )


def test_orca_scan_declares_the_surface_it_was_run_to_produce():
    """A relaxed scan could execute and promise nothing about its own profile.

    ORCA's reader has carried the scan accessors since the scan family was
    added, but no ``scan`` job type was declared, so coverage read *unknown*
    for the one job type whose entire purpose is producing a surface.
    """

    reader = reader_for("orca")
    declared = dict(reader.jobtype_selectors)
    surface = {"scan_coordinate_values", "scan_energies", "scan_point_indices"}
    assert surface <= set(declared["scan"])
    # A scan runs no frequency step, so nothing derived from a Hessian is
    # promised for it.
    assert not (
        {"vibrational_frequencies", "gibbs_free_energy"}
        & set(declared["scan"])
    )

    output = reader.open_output(
        Path("tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out")
    )
    coordinates, _ = reader.read(output, "scan_coordinate_values")
    energies, _ = reader.read(output, "scan_energies")
    assert list(coordinates) == [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]
    assert len(energies) == len(coordinates)
    # The barrier is the spread of the surface, which is why the profile is
    # exposed as two parallel vectors rather than one opaque record list.
    spread_kcal = (max(energies) - min(energies)) * 627.5094740631
    assert spread_kcal == pytest.approx(8.257, abs=0.01)


@pytest.mark.parametrize("jobtype", ("hess", "opt", "sp"))
def test_xtb_declares_the_geometry_its_reader_already_establishes(jobtype):
    # XTBOutput has always exposed the molecule, so an xTB optimisation could
    # run and still not promise that it yields a geometry.
    declared = dict(reader_for("xtb").jobtype_selectors)
    assert {"connectivity", "positions", "symbols"} <= set(declared[jobtype])


def test_xtb_optimised_geometry_enters_the_shared_quantity_plane():
    path = "tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out"
    receipt = extract_trusted_result_quantities(
        artifact=_artifact(path, "xtb", "geometry"),
        program="xtb",
        selectors=(
            QuantitySelectorV1(quantity_id="r", selector="positions"),
            QuantitySelectorV1(quantity_id="z", selector="symbols"),
            QuantitySelectorV1(quantity_id="c", selector="connectivity"),
        ),
    )
    values = {item.quantity_id: item.value for item in receipt.quantities}
    assert values["z"] == ("O", "O", "C")
    # Linear CO2: both oxygens bond the carbon and neither bonds the other.
    assert values["c"] == ((0, 0, 1), (0, 0, 1), (1, 1, 0))
    assert values["r"][0][0] == pytest.approx(-values["r"][1][0])


def test_gaussian_does_not_infer_multiplicity_from_open_shell_td_labels():
    output = reader_for("gaussian").open_output(Path(_GAUSSIAN_TD_LOG))
    labels, _ = reader_for("gaussian").read(output, "excited_state_labels")
    assert labels[0] == "2.316-A"
    with pytest.raises(MissingQuantityError, match="source-labelled"):
        reader_for("gaussian").read(output, "excited_state_multiplicities")


def test_orca_tda_roots_are_filtered_only_by_printed_multiplicity(tmp_path):
    output_path = tmp_path / "singlets-and-triplets.out"
    output_path.write_text(
        "STATE  1:  E= 0.123456 au  3.3594 eV  27096.0 cm**-1 "
        "<S**2> = 0.000000 Mult 1\n"
        "STATE  2:  E= 0.100000 au  2.7211 eV  21947.0 cm**-1 "
        "<S**2> = 2.000000 Sym: A' Mult 3\n"
        "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n"
        "  0-1A -> 1-1A 3.3594 27096.0 369.1 0.010000000 0.0 0.0\n"
        "  0-1A -> 1-3A 2.7211 21947.0 455.6 0.000000000 0.0 0.0\n"
        "------------------------------------------------------------\n",
        encoding="utf-8",
    )
    output = reader_for("orca").open_output(output_path)
    reader = reader_for("orca")

    assert reader.read(output, "singlet_excitation_energies") == (
        [3.3594],
        "eV",
    )
    assert reader.read(output, "triplet_excitation_energies") == (
        [2.7211],
        "eV",
    )
    assert reader.read(output, "excited_state_multiplicities") == (
        [1, 3],
        "",
    )
    assert reader.read(output, "triplet_oscillator_strengths") == (
        [0.0],
        "",
    )


def test_stability_and_spin_diagnostics_remain_distinct_observations():
    gaussian = reader_for("gaussian").open_output(
        Path("tests/data/GaussianTests/outputs/link/dna_link_sp.log")
    )
    verdict, _ = reader_for("gaussian").read(
        gaussian, "wavefunction_stability_verdict"
    )
    history, _ = reader_for("gaussian").read(
        gaussian, "wavefunction_stability_history"
    )
    assert history == [
        "internal_instability",
        "stable_under_considered_perturbations",
    ]
    assert verdict == "stable_under_considered_perturbations"

    orca = reader_for("orca").open_output(
        Path("tests/data/ORCATests/outputs/fe3_sextet.out")
    )
    observed, _ = reader_for("orca").read(orca, "spin_square")
    target, _ = reader_for("orca").read(orca, "spin_square_target")
    deviation, _ = reader_for("orca").read(orca, "spin_square_deviation")
    assert observed == pytest.approx(8.759007)
    assert target == pytest.approx(8.75)
    assert deviation == pytest.approx(observed - target)
    with pytest.raises(ValueError, match="does not provide"):
        reader_for("orca").read(orca, "wavefunction_stability_verdict")


def test_pyscf_spin_diagnostics_preserve_observation_target_and_deviation():
    from types import SimpleNamespace

    import numpy as np

    output = SimpleNamespace(
        multiplicity=2,
        result_units={
            "results/spin_square": "dimensionless",
            "results/spin_square_effective_multiplicity": ("dimensionless"),
        },
        results={
            "spin_square": np.asarray(0.80),
            "spin_square_effective_multiplicity": np.asarray(
                (1.0 + 4.0 * 0.80) ** 0.5
            ),
        },
    )

    def value(selector):
        read, _unit = reader_for("pyscf").read(output, selector)
        return read

    assert value("spin_square") == pytest.approx(0.80)
    assert value("spin_square_target") == pytest.approx(0.75)
    assert value("spin_square_deviation") == pytest.approx(0.05)
    assert value("effective_multiplicity") == pytest.approx(
        (1.0 + 4.0 * 0.80) ** 0.5
    )


def test_pyscf_connectivity_uses_the_structured_final_geometry():
    from types import SimpleNamespace

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule(
        symbols=["H", "H", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
    )
    output = SimpleNamespace(
        positions=molecule.positions,
        symbols=molecule.chemical_symbols,
        result_units={"results/positions": "Angstrom"},
        get_molecule=lambda: molecule,
    )
    value, _unit = reader_for("pyscf").read(output, "connectivity")
    assert value == [[0, 1, 0], [1, 0, 0], [0, 0, 0]]

    # The stored unit is part of the machine contract, so a geometry written
    # under a unit this reader does not read it as is a divergence to state
    # rather than a number to convert.
    drifted = SimpleNamespace(
        positions=molecule.positions,
        symbols=molecule.chemical_symbols,
        result_units={"results/positions": "Bohr"},
        get_molecule=lambda: molecule,
    )
    with pytest.raises(Exception, match="units are absent or incompatible"):
        reader_for("pyscf").read(drifted, "connectivity")


def test_geometry_connectivity_recognizes_a_peroxide_bond():
    from types import SimpleNamespace

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule(
        symbols=["H", "O", "O", "H"],
        positions=np.array(
            [
                [-0.96, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.45, 0.0, 0.0],
                [2.41, 0.0, 0.0],
            ]
        ),
    )

    connectivity, unit = reader_for("xyz").read(
        SimpleNamespace(molecule=molecule), "connectivity"
    )

    assert unit == ""
    assert connectivity == [
        [0, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [0, 0, 1, 0],
    ]


def test_geometry_connectivity_keeps_separated_atoms_disconnected():
    """Genuinely separated atoms stay disconnected.

    The fixture is re-derived from chemistry, because the previous one was
    geometrically impossible and its assertion was true only by accident.
    It placed O-O at 2.000 A -- far longer than a peroxide bond (1.45) and
    far shorter than any real non-bonded O...O contact (2.7+) -- and O-H
    at **1.200 A**, which is a proton in transit between two oxygens, not
    a separated atom. The all-zeros assertion held only because the
    perceiver shrank every X-H tolerance to 0.05 A, so the suite had
    encoded that defect as a requirement: the same 1.200 A O-H that a
    proton-transfer saddle carries was pinned as *not a bond*.

    These distances are ones no convention should call bonded: an O...O
    contact at 3.000 A and an O...H at 2.000 A, both well beyond a
    forming bond.
    """

    from types import SimpleNamespace

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule(
        symbols=["O", "O", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [5.0, 0.0, 0.0]]
        ),
    )

    connectivity, _ = reader_for("xyz").read(
        SimpleNamespace(molecule=molecule), "connectivity"
    )

    assert connectivity == [[0, 0, 0], [0, 0, 0], [0, 0, 0]]


def test_geometry_connectivity_sees_a_proton_in_transit():
    """The case the old fixture asserted away.

    An O-H of 1.200 A is what a proton-transfer transition state carries,
    and this laboratory computes aqueous pKa cycles. Perceiving it is the
    point: a convention that drops it reports a transferring hydrogen as
    bonded to nothing.
    """

    from types import SimpleNamespace

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule

    # O-H...O, the hydrogen 1.200 A from one oxygen and 1.300 A from the
    # other: both partial bonds, as a transfer saddle has.
    molecule = Molecule(
        symbols=["O", "H", "O"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [2.5, 0.0, 0.0]]
        ),
    )

    connectivity, _ = reader_for("xyz").read(
        SimpleNamespace(molecule=molecule), "connectivity"
    )

    assert connectivity[0][1] == 1, (
        "the transferring hydrogen is not perceived as bonded to the "
        "oxygen it is leaving, so a transfer saddle would be delivered "
        "with a hydrogen bonded to nothing"
    )
    # The two oxygens are 2.500 A apart, which is no O-O bond.
    assert connectivity[0][2] == 0

    # And the other half of the transfer, stated as it is rather than as
    # one would wish. The O-H cutoff is 1.261 A, so the approaching side
    # at 1.300 A falls outside by 3.1%. No distance convention resolves a
    # proton in transit -- that is the irreducible case -- so what the
    # host owes the reader is the number, not a verdict dressed as one.
    from chemsmart.io.molecules.perception import perceive_pairs

    pairs = {
        (pair.first_index, pair.second_index): pair
        for pair in perceive_pairs(
            molecule.chemical_symbols,
            molecule.positions,
            include_rejected=True,
        )
    }
    approaching = pairs[(1, 2)]
    assert not approaching.adjacent
    assert approaching.margin_angstrom < 0.0
    assert abs(approaching.relative_margin) < 0.05, (
        "the approaching O-H is rejected by more than 5% of the cutoff, "
        "so the host is no longer reporting a near-threshold decision as "
        f"near: margin {approaching.margin_angstrom:+.4f} A of cutoff "
        f"{approaching.cutoff_angstrom:.4f} A"
    )


def test_irc_trajectory_selectors_report_observed_endpoint_topology():
    from types import SimpleNamespace

    import numpy as np

    from chemsmart.io.molecules.structure import Molecule

    start = Molecule(
        symbols=["H", "H", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [0.7, 0.0, 0.0], [3.0, 0.0, 0.0]]
        ),
    )
    end = Molecule(
        symbols=["H", "H", "H"],
        positions=np.array(
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.7, 0.0, 0.0]]
        ),
    )
    output = SimpleNamespace(jobtype="ircf", all_structures=[start, end])
    reader = reader_for("gaussian")

    assert reader.read(output, "trajectory_frame_count") == (2, "")
    assert reader.read(output, "irc_direction") == ("forward", "")
    assert reader.read(output, "trajectory_connectivity_changed") == (1, "")
    start_graph, _ = reader.read(output, "trajectory_start_connectivity")
    end_graph, _ = reader.read(output, "trajectory_end_connectivity")
    assert start_graph != end_graph

    orca = SimpleNamespace(
        jobtype="irc",
        contents=["|  7>   Direction both"],
        all_structures=[start, end],
    )
    assert reader_for("orca").read(orca, "irc_direction") == ("both", "")
    orca.contents = []
    with pytest.raises(MissingQuantityError, match="explicitly establish"):
        reader_for("orca").read(orca, "irc_direction")


def test_multiframe_xyz_sidecar_exposes_trajectory_without_inferred_direction(
    tmp_path,
):
    trajectory = tmp_path / "reaction_IRC_Full_trj.xyz"
    trajectory.write_text(
        "3\nframe 1\n"
        "H 0.0 0.0 0.0\nH 0.7 0.0 0.0\nH 3.0 0.0 0.0\n"
        "3\nframe 2\n"
        "H 0.0 0.0 0.0\nH 3.0 0.0 0.0\nH 3.7 0.0 0.0\n",
        encoding="utf-8",
    )
    frame_count = _extract(
        trajectory, "xyz", "trajectory_frame_count", "frames"
    ).quantities[0]
    changed = _extract(
        trajectory,
        "xyz",
        "trajectory_connectivity_changed",
        "changed",
    ).quantities[0]
    assert frame_count.value == 2
    assert changed.value == 1
    with pytest.raises(ValueError, match="does not provide"):
        reader_for("xyz").read(
            reader_for("xyz").open_output(trajectory), "irc_direction"
        )


def test_connectivity_difference_count_matches_endpoints_by_atom_order():
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionError,
        QuantityExpressionNodeV1,
        QuantityExpressionRequestV1,
        evaluate_quantity_expression,
    )
    from chemsmart.analysis.result_quantities import (
        DIMENSIONLESS,
        make_quantity_value,
    )

    def numeric(quantity_id, value):
        return make_quantity_value(
            quantity_id=quantity_id,
            source_value=value,
            source_unit="1",
            value=value,
            unit="1",
            dimension=DIMENSIONLESS,
            evidence_ref=f"receipt:{'a' * 64};quantity:{quantity_id}",
        )

    def symbols(quantity_id, value):
        return make_quantity_value(
            quantity_id=quantity_id,
            source_value=value,
            source_unit="",
            value=value,
            unit="",
            dimension=DIMENSIONLESS,
            evidence_ref=f"receipt:{'b' * 64};quantity:{quantity_id}",
            data_kind="text_vector",
        )

    first = numeric("first-graph", ((0, 1, 0), (1, 0, 0), (0, 0, 0)))
    second = numeric("second-graph", ((0, 0, 0), (0, 0, 1), (0, 1, 0)))
    first_symbols = symbols("first-symbols", ("C", "O", "H"))
    second_symbols = symbols("second-symbols", ("C", "O", "H"))

    def evaluate(inputs):
        return evaluate_quantity_expression(
            QuantityExpressionRequestV1(
                schema_version="chemsmart.quantity-expression-request.v1",
                expression_id="compare-endpoints",
                inputs=inputs,
                nodes=(
                    QuantityExpressionNodeV1(
                        node_id="changed-edges",
                        operation="connectivity_difference_count",
                        input_ids=tuple(item.quantity_id for item in inputs),
                    ),
                ),
                output_node_ids=("changed-edges",),
            )
        )

    receipt = evaluate((first, first_symbols, second, second_symbols))
    assert receipt.outputs[0].value == 2
    assert receipt.outputs[0].unit == "1"

    reordered = symbols("second-symbols", ("C", "H", "O"))
    with pytest.raises(
        QuantityExpressionError, match="identity or atom order"
    ):
        evaluate((first, first_symbols, second, reordered))

    asymmetric = numeric("first-graph", ((0, 1, 0), (0, 0, 0), (0, 0, 0)))
    with pytest.raises(QuantityExpressionError, match="must be symmetric"):
        evaluate((asymmetric, first_symbols, second, second_symbols))

    smaller = numeric("second-graph", ((0, 1), (1, 0)))
    smaller_symbols = symbols("second-symbols", ("C", "O"))
    with pytest.raises(QuantityExpressionError, match="different atom counts"):
        evaluate((first, first_symbols, smaller, smaller_symbols))
