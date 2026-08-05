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
from chemsmart.agent.postprocessing import extract_trusted_result_quantities
from chemsmart.analysis.result_quantities import QuantitySelectorV1
from chemsmart.analysis.result_readers import (
    RESULT_READERS,
    MissingQuantityError,
    reader_for,
    registered_reader_programs,
)

_GAUSSIAN_LOG = "tests/data/GaussianTests/boltzmann/udc3_mCF3_monomer_c1.log"
_GAUSSIAN_LOG_C4 = (
    "tests/data/GaussianTests/boltzmann/udc3_mCF3_monomer_c4.log"
)


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
    assert registered_reader_programs() == ("gaussian", "orca", "xtb")


def test_a_gaussian_energy_becomes_a_hash_bound_quantity():
    receipt = _extract(_GAUSSIAN_LOG, "gaussian", "energy")
    quantity = receipt.quantities[0]
    assert quantity.value == pytest.approx(-2189.63187379)
    assert quantity.unit == "Eh"
    assert receipt.parser_id.endswith("Gaussian16Output")
    # The value carries its own digest, so a later claim cannot quietly
    # substitute a different number for the one that was measured.
    assert len(quantity.value_sha256) == 64


def test_a_gaussian_free_energy_is_read_from_the_thermochemistry_block():
    receipt = _extract(_GAUSSIAN_LOG, "gaussian", "gibbs_free_energy")
    assert receipt.quantities[0].value == pytest.approx(-2189.409887)


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
        def molecule(self):
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
