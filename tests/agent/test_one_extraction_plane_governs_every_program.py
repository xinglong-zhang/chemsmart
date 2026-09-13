"""A structured result answers the same vocabulary as a parsed log.

PySCF used to be a second extraction plane, reached by name.  It had its own
selector vocabulary written as an ``if``/``elif`` chain, its own unit table,
and -- the part that mattered -- no job-type declaration gate, because
``reader_for("pyscf")`` returned ``None`` and every consulting site was
written to skip it.  Three things followed.  A single point could be asked
for frequencies and only a runtime absence message stood between the plan and
the answer.  A capability query could not report that this program answers
any selector at all, so the twelve that did work were undiscoverable.  And a
plan naming a selector this program never implemented passed planning: two
such plans sat in this repository's own fixtures, one asking a Hessian for
``frequencies`` -- a name in no program's vocabulary -- and one asking a
PySCF single point for ORCA's ``scf_energy``.  Both would have been refused
only after every engine had finished.
"""

from __future__ import annotations

import pathlib

import pytest

from chemsmart.analysis.result_quantities import (
    SUPPORTED_SELECTORS,
    QuantityContractError,
    QuantitySelectorV1,
    ResultQuantityExtractionRequestV1,
)
from chemsmart.analysis.result_readers import (
    RESULT_READERS,
    SELECTOR_UNITS,
    reader_for,
    registered_reader_programs,
)

_ROOT = pathlib.Path(__file__).resolve().parents[2]


def test_pyscf_is_registered_like_every_other_program():
    assert "pyscf" in registered_reader_programs()
    reader = reader_for("pyscf")
    assert reader.artifact_kind == "pyscf_hdf5"
    assert reader.parser_id == "chemsmart.io.pyscf.output.PySCFOutput"


def test_the_extraction_dispatch_no_longer_branches_on_the_program():
    """The seam was four ``if program == "pyscf"`` sites, not a rewrite.

    What must not come back is a reader-registration check that exempts one
    program by name: those two lines are why a request naming an unusable
    selector could be constructed at all.  Branches that are genuinely about
    the HDF5 format -- the admission guard, a legacy receipt digest, an
    argument check on a named entry point -- are not dispatch and stay.
    """

    dispatch = (_ROOT / "chemsmart/agent/postprocessing.py").read_text(
        encoding="utf-8"
    )
    assert "pyscf" not in dispatch

    quantities = (_ROOT / "chemsmart/analysis/result_quantities.py").read_text(
        encoding="utf-8"
    )
    for exemption in (
        'if self.program != "pyscf":',
        'if normalized_program != "pyscf":',
    ):
        assert exemption not in quantities, exemption


def test_a_single_point_cannot_be_asked_for_a_hessian_quantity():
    """The declaration gate, which this plane never had."""

    reader = reader_for("pyscf")
    assert "vibrational_frequencies" in reader.selectors_for_jobtype("hess")
    assert "vibrational_frequencies" not in reader.selectors_for_jobtype("sp")


def test_a_session_can_discover_what_pyscf_answers():
    """Coverage is what the capability query projects; it was empty before."""

    reader = reader_for("pyscf")
    for jobtype in ("sp", "opt", "hess", "td"):
        assert reader.selectors_for_jobtype(jobtype)
    # ``td`` is executable (contract v5): an approved workflow emits an
    # excited state, so the response set is declared for it -- and for
    # ``opt``, whose excited-root form re-evaluates the spectrum at the
    # reached geometry -- while ``sp`` and ``hess`` still refuse it.
    td = reader.selectors_for_jobtype("td")
    assert "excitation_energies" in td and "oscillator_strengths" in td
    assert "excitation_energies" in reader.selectors_for_jobtype("opt")
    for jobtype in ("sp", "hess"):
        assert "excitation_energies" not in reader.selectors_for_jobtype(
            jobtype
        )


def test_an_unregistered_program_is_still_refused_by_name():
    with pytest.raises(QuantityContractError, match="no result reader"):
        ResultQuantityExtractionRequestV1(
            schema_version="chemsmart.quantity-extraction-request.v1",
            artifact_id="a1",
            artifact_sha256="a" * 64,
            program="nwchem",
            selectors=(
                QuantitySelectorV1(quantity_id="e", selector="energy"),
            ),
        )


def test_the_hartree_versus_electronvolt_disagreement_is_one_entry():
    """A recorded cross-program disagreement, closed where units belong.

    PySCF stores excitation energies in hartree; ORCA and Gaussian print
    electronvolts.  That is a unit, not a different quantity, so it belongs
    in the reader's native-unit override and never in the arithmetic.
    """

    assert SELECTOR_UNITS["excitation_energies"] == "eV"
    assert reader_for("pyscf").source_units["excitation_energies"] == "Eh"
    for program in ("orca", "gaussian"):
        overrides = RESULT_READERS[program].source_units
        assert "excitation_energies" not in overrides


def test_every_pyscf_accessor_is_in_the_shared_vocabulary():
    """The invariant the separate plane was outside of."""

    accessors = set(reader_for("pyscf").accessors)
    assert accessors <= set(SUPPORTED_SELECTORS)
    assert accessors <= set(SELECTOR_UNITS)
