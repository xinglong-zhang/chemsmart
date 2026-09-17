"""Per-atom charges, ordered by the molecule rather than by a label.

Every supported program parses atomic charges and the typed layer exposed
none of them, so the Fukui composition -- which the expression layer already
evaluates unmodified -- had no input, and no session could ask where a charge
sits.

Two facts decide the shape. Atom-label schemes disagree between programs:
ORCA and Gaussian number atoms globally, so carbon dioxide's carbon is "C3",
while xTB counts within each element and calls the same atom "C1". And the
typed layer cannot repair the order afterwards, because freezing a mapping
sorts it by label, turning O1,C2,H3 into C2,H3,O1. So the label is resolved
against the molecule's own symbols at the accessor, and what leaves is a
positional vector.

The scheme is in the selector name because Mulliken, Loewdin, Hirshfeld, CM5
and a tight-binding population are different quantities, not different
spellings of one. The last test here shows how different: on one anion the two
schemes disagree about the charge on oxygen by a factor of three.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.analysis.result_quantities import CHARGE, SUPPORTED_SELECTORS
from chemsmart.analysis.result_readers import (
    _SELECTOR_DIMENSIONS,
    RESULT_READERS,
    SELECTOR_UNITS,
    MissingQuantityError,
    _per_atom_vector,
)

_CO2 = Path("tests/data/ORCATests/outputs/CO2.out")
_PHENOXIDE = Path(
    "/home/chemsmart/agent-campaigns/ax41-refine-100/qualification"
    "/pcet-analysis-path/workspace/phenoxide-optfreq.out"
)
#: The two ORCA prints without being asked.
_DEFAULT_SCHEMES = ("mulliken_atomic_charges", "loewdin_atomic_charges")
#: Every declared scheme, including the one a route directive must request.
_SCHEMES = _DEFAULT_SCHEMES + ("hirshfeld_atomic_charges",)
_XTB_SCC = "xtb_scc_atomic_charges"
_HIRSHFELD = Path("tests/data/ORCATests/outputs/udc3_ts1_c15_sp_hirshfeld.out")


def test_a_global_index_scheme_orders_by_the_molecule():
    ordered = _per_atom_vector(
        {"O1": -0.42, "C2": 0.13, "H3": 0.11},
        ["O", "C", "H"],
        quantity="probe",
    )
    assert ordered == [-0.42, 0.13, 0.11]


def test_a_per_element_scheme_is_refused_rather_than_reordered():
    """xTB's own test asserts this exact mapping for CO2, order O,O,C.

    Its carbon is atom 3 and is keyed "C1". Read as a global index that names
    atom 1, which is an oxygen -- so the mismatch is caught by the molecule's
    own symbols rather than passing quietly and returning another atom's
    charge.
    """

    with pytest.raises(MissingQuantityError, match="atom-label scheme"):
        _per_atom_vector(
            {"O1": -0.232, "O2": -0.232, "C1": 0.464},
            ["O", "O", "C"],
            quantity="probe",
        )


@pytest.mark.parametrize(
    "labelled,symbols,reason",
    [
        ({"O1": -0.4}, ["O", "H"], "value count"),
        ({"O1": -0.4, "O9": 0.4}, ["O", "H"], "out-of-range index"),
        ({"O1": -0.4, "O1x": 0.4}, ["O", "H"], "unreadable label"),
    ],
)
def test_a_mapping_that_does_not_describe_this_molecule_is_refused(
    labelled, symbols, reason
):
    with pytest.raises(MissingQuantityError):
        _per_atom_vector(labelled, symbols, quantity="probe")


@pytest.mark.parametrize("selector", _SCHEMES)
def test_each_scheme_is_typed_as_a_charge(selector):
    assert selector in SUPPORTED_SELECTORS
    assert SELECTOR_UNITS[selector] == "e"
    import chemsmart.analysis.result_quantities as rq

    assert getattr(rq, _SELECTOR_DIMENSIONS[selector]) == CHARGE


def test_xtb_scc_population_is_typed_and_stays_out_of_mulliken_names():
    """The xTB sidecar is positional, but its SCC scheme is not Mulliken."""

    import chemsmart.analysis.result_quantities as rq

    assert _XTB_SCC in SUPPORTED_SELECTORS
    assert SELECTOR_UNITS[_XTB_SCC] == "e"
    assert getattr(rq, _SELECTOR_DIMENSIONS[_XTB_SCC]) == CHARGE
    reader = RESULT_READERS["xtb"]
    output = reader.open_output(
        Path("tests/data/XTBTests/outputs/co2_ohess/co2_ohess.out")
    )
    values, unit = reader.read(output, _XTB_SCC)
    symbols, _ = reader.read(output, "symbols")
    total, _ = reader.read(output, "charge")
    assert unit == "e"
    assert len(values) == len(symbols) == 3
    assert sum(values) == pytest.approx(float(total), abs=1e-3)
    assert "mulliken_atomic_charges" not in reader.accessors


@pytest.mark.parametrize("selector", _DEFAULT_SCHEMES)
def test_the_charges_sum_to_the_molecular_charge(selector):
    """The physical checksum that also proves the ordering is complete.

    A population analysis partitions all the electrons, so the per-atom
    charges must sum to the molecule's total charge. A vector that dropped,
    duplicated or misplaced an atom would not.
    """

    reader = RESULT_READERS["orca"]
    output = reader.open_output(_CO2)
    charges, unit = reader.read(output, selector)
    symbols, _ = reader.read(output, "symbols")
    total, _ = reader.read(output, "charge")
    assert unit == "e"
    assert len(charges) == len(symbols)
    assert sum(charges) == pytest.approx(float(total), abs=1e-3)


def test_declared_only_where_the_meaning_was_audited():
    declared = dict(RESULT_READERS["orca"].jobtype_selectors)
    for jobtype in ("opt", "sp", "ts"):
        for selector in _SCHEMES:
            assert selector in declared[jobtype], (jobtype, selector)
    for jobtype in ("scan", "td", "irc"):
        for selector in _SCHEMES:
            assert selector not in declared[jobtype], (jobtype, selector)


def test_each_program_declares_only_the_schemes_its_output_carries():
    """A scheme is declared where that program computes it, and nowhere else.

    Mulliken and Loewdin are ORCA defaults printed without being asked.
    PySCF's driver computes a Mulliken population and stores it under its own
    declared unit, so the same name is honest there -- for the scheme.  It is
    not an invitation to subtract one from the other: a Mulliken charge in a
    triple-zeta basis and one from a different program at a different basis
    are the same partition of different densities.  xTB's population comes
    from a minimal tight-binding density and is not Mulliken at all, which is
    why no xTB accessor answers to this name.
    """

    expected = {
        "orca": set(_SCHEMES),
        "pyscf": {"mulliken_atomic_charges"},
        "gaussian": set(),
        "xtb": {_XTB_SCC},
        "xyz": set(),
    }
    for program, reader in RESULT_READERS.items():
        carried = {
            name for name in (*_SCHEMES, _XTB_SCC) if name in reader.accessors
        }
        assert carried == expected[program], program

    # Parsed and deliberately undeclared: reaching Hirshfeld or CM5 needs a
    # print directive through the project route hatch, and the scheme's own
    # accessor is a separate question from whether that channel may carry it.
    orca = RESULT_READERS["orca"]
    for withheld in (
        "hirshfeld_cm5_charges",
        "loewdin_spin_densities",
    ):
        assert withheld not in orca.accessors


def test_the_two_schemes_disagree_and_that_is_why_the_name_carries_one():
    """Scheme-qualification is load-bearing, not decorative.

    On the phenoxide anion Mulliken places more than a whole electron of
    excess charge on the hydroxyl oxygen while Loewdin places about a third
    of one. Both are correct readings of their own definitions; neither is
    "the charge on the oxygen". A selector called `atomic_charges` would have
    invited exactly that reading.
    """

    if not _PHENOXIDE.exists():
        pytest.skip("campaign fixture not present on this host")
    reader = RESULT_READERS["orca"]
    output = reader.open_output(_PHENOXIDE)
    symbols, _ = reader.read(output, "symbols")
    oxygen = symbols.index("O")
    mulliken, _ = reader.read(output, "mulliken_atomic_charges")
    loewdin, _ = reader.read(output, "loewdin_atomic_charges")
    assert mulliken[oxygen] < -0.9
    assert loewdin[oxygen] > -0.5
    assert abs(mulliken[oxygen] - loewdin[oxygen]) > 0.5


_OPEN_SHELL = [
    ("fe3_doublet", 3, 2),
    ("fe3_quartet", 3, 4),
    ("fe3_sextet", 3, 6),
    ("fe2_triplet", 2, 3),
]


@pytest.mark.parametrize("stem,charge,multiplicity", _OPEN_SHELL)
@pytest.mark.parametrize("selector", _DEFAULT_SCHEMES)
def test_an_open_shell_result_reports_charges_not_spin(
    stem, charge, multiplicity, selector
):
    """The defect a live session found by adding up what it was given.

    ORCA prints one column for a closed shell and two for an open shell --
    "MULLIKEN ATOMIC CHARGES" becomes "... AND SPIN POPULATIONS", charge
    first and spin second -- and the header substring matches both. Reading
    the last token therefore returned the charge of a restricted result and
    the spin population of an unrestricted one, under one name.

    It went unseen because nothing in the typed layer had ever read these
    values; declaring the selector is what made it reachable, and the first
    session to use it noticed that a neutral phenoxyl radical's charges summed
    to +1.00 e. A doublet cation had looked correct only because its total
    spin and its formal charge are both +1 -- which is why the sextet below is
    the decisive case: its charge is +3 and its spin is +5.
    """

    reader = RESULT_READERS["orca"]
    output = reader.open_output(
        Path(f"tests/data/ORCATests/outputs/{stem}.out")
    )
    charges, unit = reader.read(output, selector)
    assert unit == "e"
    total = sum(charges)
    assert total == pytest.approx(float(charge), abs=1e-3)
    if multiplicity - 1 != charge:
        # the old reading; it must not come back
        assert total != pytest.approx(float(multiplicity - 1), abs=1e-3)


def test_the_charge_column_is_read_by_position_not_from_the_end():
    """Pinned directly against both printed layouts.

    A closed-shell row is "0 C : 1.034841" and an open-shell row is
    "0 C : 0.850767 0.111279". The charge is the first number after the colon
    in both; only the closed-shell case makes it the last one too.
    """

    reader = RESULT_READERS["orca"]
    open_shell = reader.open_output(
        Path("tests/data/ORCATests/outputs/fe3_sextet.out")
    )
    charges, _ = reader.read(open_shell, "mulliken_atomic_charges")
    text = Path("tests/data/ORCATests/outputs/fe3_sextet.out").read_text()
    # The colon may be attached to the element ("Fe:") or stand alone
    # ("O :"), which is itself why the parser locates it rather than
    # assuming a column count.
    start = max(
        index
        for index, line in enumerate(text.splitlines())
        if "MULLIKEN ATOMIC CHARGES" in line
    )
    first_row = text.splitlines()[start + 2].split()
    assert first_row[0] == "0", "expected the first atom's row"
    colon = next(
        index for index, token in enumerate(first_row) if token.endswith(":")
    )
    assert len(first_row) - colon - 1 == 2, "row must have two numbers"
    assert charges[0] == pytest.approx(float(first_row[colon + 1]))
    assert charges[0] != pytest.approx(float(first_row[-1]))


def test_hirshfeld_reads_every_atom_of_a_run_that_asked_for_it():
    """Reachability is three conditions, and this fixture carries all three.

    The parser produces it, the typed layer can be asked for it, and the
    Agent can make ORCA emit it -- ``Hirshfeld`` is a route token appended to
    the route line, which is why this archived job has the block at all.
    Hydrogens matter here in a way they do not for a heavy-atom summary: a
    condensed Fukui function needs every atom or it does not sum to one.
    """

    reader = RESULT_READERS["orca"]
    output = reader.open_output(_HIRSHFELD)
    charges, unit = reader.read(output, "hirshfeld_atomic_charges")
    symbols, _ = reader.read(output, "symbols")
    assert unit == "e"
    assert len(charges) == len(symbols) == 78
    assert "H" in symbols
    assert charges[0] == pytest.approx(-0.177356)


def test_a_basin_partition_closes_on_a_grid_and_a_basis_partition_exactly():
    """The checksum holds for all three, but not to the same tolerance.

    Mulliken and Loewdin divide a sum over basis functions, so their only
    residual is the six printed decimals -- 2e-6 over 78 atoms.  Hirshfeld
    divides real space against a promolecular reference on a numerical grid,
    and closes two orders of magnitude looser at 3e-4.  Reading that residual
    as a defect would be wrong, and reading it as licence to skip the
    checksum would be worse: the checksum is what proves the vector is
    complete and in molecular order, and 3e-4 is still nowhere near the
    0.5 e a dropped or duplicated atom would cost.
    """

    reader = RESULT_READERS["orca"]
    output = reader.open_output(_HIRSHFELD)
    total, _ = reader.read(output, "charge")
    sums = {}
    for selector in _SCHEMES:
        charges, _ = reader.read(output, selector)
        sums[selector] = sum(charges)
    for selector in _DEFAULT_SCHEMES:
        assert sums[selector] == pytest.approx(float(total), abs=1e-5)
    basin = abs(sums["hirshfeld_atomic_charges"] - float(total))
    basis = abs(sums["mulliken_atomic_charges"] - float(total))
    assert basin < 1e-3
    assert basin > 100 * basis


def test_three_schemes_three_answers_on_one_oxygen():
    """Why the scheme is in the name, extended to the recommended one."""

    reader = RESULT_READERS["orca"]
    output = reader.open_output(_HIRSHFELD)
    symbols, _ = reader.read(output, "symbols")
    assert symbols[0] == "O"
    values = {
        selector: reader.read(output, selector)[0][0] for selector in _SCHEMES
    }
    assert values["mulliken_atomic_charges"] == pytest.approx(
        -0.5444, abs=1e-3
    )
    assert values["loewdin_atomic_charges"] == pytest.approx(-0.1536, abs=1e-3)
    assert values["hirshfeld_atomic_charges"] == pytest.approx(
        -0.1774, abs=1e-3
    )
    spread = max(values.values()) - min(values.values())
    assert spread > 0.35


def test_a_run_that_never_asked_has_no_hirshfeld_analysis():
    """Absence is meaning, and it must not arrive as a parser traceback.

    ORCA prints the block only on request, so the accessor returns nothing
    rather than indexing an empty list -- which is what it used to do.
    """

    reader = RESULT_READERS["orca"]
    output = reader.open_output(_CO2)
    with pytest.raises(MissingQuantityError, match="hirshfeld"):
        reader.read(output, "hirshfeld_atomic_charges")
