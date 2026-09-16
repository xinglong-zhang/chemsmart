"""An exactly idealised coordinate is a saddle however it was built.

`idealised_internal_coordinate_count` exists for one hazard and names it
in its own docstring: *"A builder that types 60/180/300 for three methyl
hydrogens has placed a threefold rotor exactly on its saddle."* It took
append receipts and nothing else, and the geometry-chain walk stepped
*past* every edit receipt on its way to the appends.

Measured on the live CUHK run (goal `butane-wave-2`): the Agent built its
conformer set by **editing** the C1-C2-C3-C4 torsion to exactly 0.00° and
exactly 180.00°. Both are on the 60° lattice the function already tests
for, and the 0.00° one is butane's syn-periplanar saddle -- the top of
the rotational barrier, not a conformer. It started at 0.00°, did not
move, and finished 5.673 kcal/mol above anti. Nothing observed it,
because the receipt that would have said so was the wrong kind.

The charter's own sentence is narrow in the same way -- it says "the
count of *appended atoms* placed on the exact 60° torsion lattice" --
while the evidence that paragraph cites includes an idealised D4h start,
which is no append at all. A declaration that under-describes its own
motivating cases is the same defect one layer up.
"""

from __future__ import annotations

from types import SimpleNamespace


from chemsmart.agent.symmetry import (
    idealised_coordinate_observation,
    idealised_internal_coordinate_count,
)


def _edit(operation, requested, unit="degree"):
    return SimpleNamespace(
        operation=operation,
        value_requested=requested,
        value_unit=unit,
    )


def _append(dihedral, angle=109.4712206):
    return SimpleNamespace(dihedral_degrees=dihedral, angle_degrees=angle)


def test_an_edit_to_an_exact_lattice_torsion_is_counted():
    counts = idealised_internal_coordinate_count(
        (), edit_receipts=(_edit("dihedral", 0.0),)
    )
    assert counts["idealised_torsions"] == 1
    assert counts["built_coordinates"] == 1


def test_the_live_run_s_two_edits_are_both_counted():
    """0.00 and 180.00 degrees, exactly as the Agent typed them."""

    counts = idealised_internal_coordinate_count(
        (),
        edit_receipts=(
            _edit("dihedral", 0.0),
            _edit("dihedral", 180.0),
            _edit("dihedral", -60.0),
        ),
    )
    assert counts["idealised_torsions"] == 3
    observation = idealised_coordinate_observation(counts)
    assert "3 of 3" in observation
    assert "saddle" in observation


def test_a_torsion_off_the_lattice_is_not_counted():
    counts = idealised_internal_coordinate_count(
        (), edit_receipts=(_edit("dihedral", -65.4),)
    )
    assert counts["idealised_torsions"] == 0
    assert counts["built_coordinates"] == 1


def test_an_edited_angle_at_an_idealised_value_is_counted():
    counts = idealised_internal_coordinate_count(
        (), edit_receipts=(_edit("angle", 109.4712206),)
    )
    assert counts["idealised_angles"] == 1


def test_a_bond_length_edit_is_neither():
    """Only the two coordinates that have idealised values."""

    counts = idealised_internal_coordinate_count(
        (), edit_receipts=(_edit("bond", 1.09, unit="angstrom"),)
    )
    assert counts["idealised_torsions"] == 0
    assert counts["idealised_angles"] == 0
    assert counts["built_coordinates"] == 1


def test_appends_still_count_exactly_as_they_did():
    counts = idealised_internal_coordinate_count((_append(60.0),))
    assert counts["idealised_torsions"] == 1
    assert counts["idealised_angles"] == 1
    assert counts["built_coordinates"] == 1
    assert counts["appended_atoms"] == 1


def test_the_symmetry_walk_offers_its_edits_to_the_sensor():
    """The walk stepped past every edit on its way to the appends."""

    import inspect

    from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

    source = inspect.getsource(
        CommandCompiledToolHostV1._symmetry_observations
    )
    assert "edit_receipts=" in source, (
        "the chain walk collects appends and steps past edits, so a "
        "coordinate edited onto its own saddle is never observed"
    )
    assert "edits.append(edit_receipt)" in source


def test_the_observation_reports_the_lattice_and_not_the_chemistry():
    """"starts on its own saddle" is periodicity-dependent, and the host
    does not know the periodicity.

    For the threefold methyl rotor this function was written from,
    60/180/300 are all saddles and the sentence holds. For butane's
    C-C-C-C backbone torsion it is false in the most visible way
    possible: 180 degrees is the global minimum and 0 degrees is the
    saddle. The host cannot tell the two apart -- it sees a number on a
    lattice, not a rotor's order.

    On the live run this mattered immediately: `opt-anti-r2` was built by
    an edit to exactly 180.00 degrees, so it carried a host observation
    saying it starts on its own saddle while its Hessian was computing
    the frequencies that say otherwise. The charter's own line is that
    the host reports what its convention said and how narrowly, and the
    scientist draws the chemical conclusion.
    """

    counts = idealised_internal_coordinate_count(
        (), edit_receipts=(_edit("dihedral", 180.0),)
    )
    observation = idealised_coordinate_observation(counts)

    # The measurement stays.
    assert "1 of 1" in observation
    assert "60°" in observation

    # The conclusion goes.
    assert "starts on its own saddle" not in observation
    # What is actually known: it is a stationary point of the symmetry it
    # was built with, and which kind is a question for the frequencies.
    assert "stationary point" in observation
    assert "frequencies" in observation or "minimum or a saddle" in (
        observation
    )
