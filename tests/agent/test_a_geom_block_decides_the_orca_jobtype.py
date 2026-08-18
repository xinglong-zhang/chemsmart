"""ORCA writes a scan as an optimisation, so the route line cannot decide.

A relaxed scan and a constrained optimisation are both `! Opt` in ORCA's
simple input; what makes them what they are is the `%geom` block underneath.
The reader classified from the route line alone, so a correctly compiled scan
came back as `opt`.

The cost was measured, not imagined. In the scan qualification session the
Agent compiled `run orca scan --coordinates [[3, 1, 2, 4]] --dist-start 0.0
--dist-end 360.0 --num-steps 25`, ChemSmart emitted exactly the right native
input, and safe preview then refused it as `preview.semantic.mismatch`, jobtype
expected 'scan', observed 'opt'. The session tried seven project variants
looking for a setting that would change the emitted input, correctly reported
that the input was byte-identical every time, and ended with the stage blocked.
Nothing it could have done would have worked, because nothing was wrong with
what it had built.

The precedent for the fix is one property up in the same file: a fixed-geometry
TD calculation has no `TD` keyword either, and is recognised from its %tddft
block rather than its route.

The fixtures are inputs ChemSmart actually generated through the live CLI for
the three cases, kept as written.
"""

from __future__ import annotations

import pytest

from chemsmart.io.orca.input import ORCAInput
from chemsmart.io.orca.output import ORCAOutput

_WRITTEN = "tests/data/ORCATests/written_files"
_SCAN = f"{_WRITTEN}/hooh_relaxed_scan_geom_block.inp"
_MODRED = f"{_WRITTEN}/hooh_constrained_modred_geom_block.inp"
_PLAIN_OPT = f"{_WRITTEN}/hooh_plain_opt.inp"
_COMPLETED_SCAN = "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out"


def test_the_route_line_alone_still_says_opt():
    """The premise. If ORCA ever stops writing `! Opt` this test says so."""

    assert ORCAInput(_SCAN).route_object.jobtype == "opt"
    assert ORCAInput(_MODRED).route_object.jobtype == "opt"


def test_a_driven_coordinate_makes_it_a_scan():
    parsed = ORCAInput(_SCAN)

    assert parsed.jobtype == "scan"
    assert parsed.scan_coordinates == ("D 2 0 1 3 = 0.0, 360.0, 25",)


def test_a_held_coordinate_makes_it_a_constrained_optimisation():
    parsed = ORCAInput(_MODRED)

    assert parsed.jobtype == "modred"
    assert parsed.constrained_coordinates == ("{B 0 1 C}",)


def test_an_optimisation_with_no_geom_block_is_still_an_optimisation():
    """The control: this fix must not turn every opt into a scan."""

    parsed = ORCAInput(_PLAIN_OPT)

    assert parsed.jobtype == "opt"
    assert parsed.scan_coordinates == ()
    assert parsed.constrained_coordinates == ()


def test_the_same_reading_works_on_an_echoed_input_in_a_result():
    """ORCA echoes the submitted input into its output with `| n>` prefixes.

    One parser has to serve both, or a completed scan would be classified
    differently from the input that produced it.
    """

    completed = ORCAOutput(_COMPLETED_SCAN)

    assert completed.jobtype == "scan"
    assert completed.scan_coordinates


@pytest.mark.parametrize("path", (_SCAN, _MODRED, _PLAIN_OPT))
def test_the_chemistry_is_unchanged_by_the_reclassification(path):
    """Only the label moves; the method must read back as it was written."""

    parsed = ORCAInput(path)

    assert parsed.functional == "b3lyp"
    assert parsed.basis in {"def2-tzvp", "def2-svp"}
