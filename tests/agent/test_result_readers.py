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

import hashlib
from pathlib import Path

import pytest

from chemsmart.agent._contracts import TrustedArtifactRefV1
from chemsmart.agent.analysis_nodes import (
    ResultQuantitySelectorV1,
    build_default_result_parser_registry,
)
from chemsmart.agent.postprocessing import extract_trusted_result_quantities
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface
from chemsmart.analysis.result_quantities import (
    QuantityExtractionError,
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


def test_the_log_parsing_programs_are_registered():
    assert registered_reader_programs() == ("gaussian", "orca", "xtb", "xyz")
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
        if item["function"]["name"] == "extract_result_quantities"
    )
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
    endpoint = tmp_path / "charged.xyz"
    endpoint.write_text(
        "1\ncharge -1.0 multiplicity 2\nH 0.0 0.0 0.0\n",
        encoding="utf-8",
    )

    with pytest.raises(QuantityExtractionError, match="explicit Energy"):
        _extract(endpoint, "xyz", "energy")

    positions = _extract(endpoint, "xyz", "positions").quantities[0]
    assert positions.value == ((0.0, 0.0, 0.0),)


@pytest.mark.parametrize("program", ("gaussian", "orca", "xtb"))
def test_the_analysis_registry_exposes_every_log_reader(program):
    reader = reader_for(program)
    adapter = build_default_result_parser_registry().resolve(
        program=program,
        artifact_kind=reader.artifact_kind,
        selectors=(
            ResultQuantitySelectorV1(
                schema_version="chemsmart.result-quantity-selector.v1",
                quantity_id="energy",
                selector="energy",
            ),
        ),
    )
    assert adapter.parser_id == reader.parser_id


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
    assert scf.value == pytest.approx(-76.05666270)
    assert reference.value == pytest.approx(scf.value)
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
    with pytest.raises(Exception, match="post-SCF correlation"):
        _extract(path, "orca", "correlation_energy", "correlation")


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
    with pytest.raises(Exception, match="no quantity reader is registered"):
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
    """Distinct from absent: xTB exposes no charge at all."""

    with pytest.raises(ValueError, match="does not provide"):
        reader_for("xtb").read(object(), "charge")


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

    from chemsmart.analysis.result_quantities import _extract_selector

    output = SimpleNamespace(
        multiplicity=2,
        results={
            "spin_square": np.asarray(0.80),
            "spin_square_effective_multiplicity": np.asarray(
                (1.0 + 4.0 * 0.80) ** 0.5
            ),
        },
    )

    def value(selector):
        return _extract_selector(
            output,
            QuantitySelectorV1(quantity_id=selector, selector=selector),
            "artifact:spin#" + "a" * 64,
        ).value

    assert value("spin_square") == pytest.approx(0.80)
    assert value("spin_square_target") == pytest.approx(0.75)
    assert value("spin_square_deviation") == pytest.approx(0.05)
    assert value("effective_multiplicity") == pytest.approx(
        (1.0 + 4.0 * 0.80) ** 0.5
    )


def test_pyscf_connectivity_uses_the_structured_final_geometry():
    from types import SimpleNamespace

    import numpy as np

    from chemsmart.analysis.result_quantities import _extract_selector
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
        get_molecule=lambda: molecule,
    )
    quantity = _extract_selector(
        output,
        QuantitySelectorV1(quantity_id="graph", selector="connectivity"),
        "artifact:geometry#" + "a" * 64,
    )
    assert quantity.value == (
        (0.0, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
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
