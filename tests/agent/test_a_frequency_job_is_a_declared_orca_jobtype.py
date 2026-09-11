"""ORCA `! Freq` without `Opt` is a declared reader jobtype.

A budget-bound session folded its frequencies into a single-point stage
and every extraction from the four resulting outputs was refused: the
reader detected jobtype `freq` and nothing had declared what its printed
values mean. A frequency job at a fixed geometry means everything an
optimisation's output means minus the optimisation claim.
"""

from pathlib import Path

import pytest

from chemsmart.analysis.result_readers import RESULT_READERS

_FREQ = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "ORCATests"
    / "outputs"
    / "He_freq.out"
)


@pytest.mark.capability("selector:orca:freq:vibrational_frequencies")
def test_freq_declares_what_opt_declares_minus_the_optimisation_claim():
    """What an optimisation has that a fixed-geometry frequency job has not.

    Two things, and both are about having moved: the claim that the
    walk converged, and the structure the walk arrived at. A ``freq``
    job's printed structure is the one it was handed, so there is no
    second structure for ``reached_positions`` to mean, and declaring
    it would make the selector's own state declaration false here.
    """

    reader = RESULT_READERS["orca"]
    freq = set(reader.selectors_for_jobtype("freq"))
    opt = set(reader.selectors_for_jobtype("opt"))
    assert freq == opt - {"converged", "reached_positions"}
    assert {"energy", "positions", "vibrational_frequencies"} <= freq


@pytest.mark.capability("selector:orca:freq:vibrational_frequencies")
def test_an_archived_frequency_output_reads_as_freq():
    reader = RESULT_READERS["orca"]
    output = reader.open_output(str(_FREQ))
    assert output.jobtype == "freq"
    assert reader.selectors_for_jobtype(output.jobtype)
