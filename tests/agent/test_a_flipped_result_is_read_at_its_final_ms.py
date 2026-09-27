"""A broken-symmetry result is read at the Ms it converged to.

ORCA's FlipSpin and BrokenSym start from the high-spin determinant the
coordinate line names and converge to the Ms the request asks for. The
atomic spin populations of that determinant sum to 2 Ms, not to 2S of the
coordinate-line multiplicity, and the population reader, which checks the
sum so that a dropped or duplicated atom cannot pass as a population,
refused every such result: on ino2's Ni(II)2 both routes to the
antiferromagnetic state were read as "sums to 0.000 where 2S for
multiplicity 5 is 4.0" (R11 truth-2, CUHK 2157086). The Ms is the
program's own record -- ORCA prints the state it converges to -- and the
coordinate-line multiplicity stays what it is, the state the flip started
from.

The fixtures are those two outputs, ORCA 6.1.1, B3LYP/G def2-SVP, charge 2,
multiplicity 5: FlipSpin 1 with FinalMs 0, and BrokenSym 2,2.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chemsmart.analysis.result_readers import reader_for

_DATA = Path(__file__).resolve().parents[1] / "data/ORCATests/broken_symmetry"

#: The two nickel atoms' final spin populations, as each output prints them.
_FINAL_NICKEL_SPIN = {
    "ni2_fs1_gas_phase.out": {
        "mulliken_atomic_spin_populations": (1.726172, -1.722662),
        "loewdin_atomic_spin_populations": (1.710982, -1.707443),
    },
    "ni2_bs22_gas_phase.out": {
        "mulliken_atomic_spin_populations": (1.726159, -1.722658),
        "loewdin_atomic_spin_populations": (1.710973, -1.707440),
    },
}


@pytest.mark.capability("selector:orca:sp:mulliken_atomic_spin_populations")
@pytest.mark.capability("selector:orca:sp:loewdin_atomic_spin_populations")
@pytest.mark.parametrize("name", sorted(_FINAL_NICKEL_SPIN))
@pytest.mark.parametrize(
    "selector",
    ["mulliken_atomic_spin_populations", "loewdin_atomic_spin_populations"],
)
def test_a_flipped_result_is_read_at_the_ms_it_converged_to(name, selector):
    reader = reader_for("orca")
    output = reader.open_output(str(_DATA / name))

    values, _unit = reader.read(output, selector)
    values = [float(item) for item in values]

    assert len(values) == 29
    assert (values[0], values[1]) == pytest.approx(
        _FINAL_NICKEL_SPIN[name][selector], abs=1e-6
    )
    assert sum(values) == pytest.approx(0.0, abs=0.05)
    # The coordinate line still names the high-spin state the flip began
    # from, which is what the request's own record compares against.
    assert output.multiplicity == 5
